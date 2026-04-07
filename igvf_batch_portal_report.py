#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from downloading_from_samplesheet import load_rows, sanitize_label, stable_sort_rows
from igvf_portal import IGVFCredentials, generate_samplesheet_rows, load_credentials, portal_json
from seqspec_parser import canonicalize_samplesheet_rows


DEFAULT_REPORT_TSV = (
    "igvf_criprs_portal_all_deposited_mar30_2026_samples/"
    "igvf_data_analysis_set_report_2026_3_30_15h_13m.tsv"
)
DEFAULT_SORT_COLUMNS = ["file_modality", "measurement_sets", "sequencing_run", "lane"]
GOOD_FLAGS = {"perfect_match", "close_enough"}


@dataclass
class SelectedSubset:
    mode: str
    rows: List[Dict[str, str]]
    note: str
    samplesheet_path: Optional[Path]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the IGVF one-lane validator in batch against a portal multireport TSV, "
            "falling back to one representative pair per modality when no shared scRNA+gRNA lane exists."
        )
    )
    parser.add_argument("--input-report", default=DEFAULT_REPORT_TSV, help="Portal multireport TSV with analysis sets")
    parser.add_argument(
        "--analysis-root",
        default="igvf_mar30_2026_one_lane_report",
        help="Root directory for per-accession outputs and the final report",
    )
    parser.add_argument("--igvf-keypair", required=True, help="IGVF credential JSON with key/secret")
    parser.add_argument(
        "--accession",
        action="append",
        default=[],
        help="Optional analysis-set accession filter. Repeat to run a subset.",
    )
    parser.add_argument("--limit", type=int, help="Optional limit after filtering")
    parser.add_argument(
        "--barcode-sample-reads",
        type=int,
        default=10000,
        help="Barcode read sample depth passed through to seqspec_parser.py",
    )
    parser.add_argument(
        "--feature-sample-reads",
        type=int,
        default=20000,
        help="Feature read sample depth passed through to seqspec_parser.py",
    )
    parser.add_argument(
        "--ignore-metadata-md5",
        action="store_true",
        help="Pass through to seqspec_parser.py to skip MD5 verification for barcode/hash/guide metadata and seqspec files",
    )
    parser.add_argument(
        "--fresh-run",
        action="store_true",
        help="Delete existing contents under --analysis-root before starting the batch",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="If used with --fresh-run, clean --analysis-root and exit without running the batch",
    )
    parser.add_argument("--force", action="store_true", help="Re-run accessions even if an outcome JSON already exists")
    return parser.parse_args()


def extract_accession(path_or_accession: str) -> str:
    token = str(path_or_accession or "").strip()
    if token.startswith("/analysis-sets/"):
        return token.strip("/").split("/")[-1]
    return token


def load_multireport_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        metadata_line = handle.readline()
        if not metadata_line:
            raise ValueError(f"{path} is empty")
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Could not read the multireport header from {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    return rows


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def write_tsv(rows: Sequence[Dict[str, str]], path: Path) -> None:
    if not rows:
        raise ValueError("No rows available to write.")
    ensure_parent(path)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def field_or_empty(row: Dict[str, str], key: str) -> str:
    return str(row.get(key, "") or "").strip()


def try_fetch_lab(accession: str, credentials: IGVFCredentials) -> str:
    try:
        payload = portal_json(f"/analysis-sets/{accession}/@@embedded?format=json", credentials)
    except Exception:
        return ""
    lab = payload.get("lab")
    if isinstance(lab, dict):
        for key in ("title", "display_title", "name"):
            value = str(lab.get(key, "") or "").strip()
            if value:
                return value
    if isinstance(lab, str):
        return lab.strip()
    return ""


def first_matching_fallback(accession: str, modality_suffixes: Sequence[str]) -> str:
    base = Path.cwd()
    for suffix in modality_suffixes:
        for pattern in (f"*_{accession}_{suffix}.yaml", f"*_{accession}_{suffix}.yml"):
            for candidate in sorted(base.glob(pattern)):
                return str(candidate)
    return ""


def fallback_seqspec_paths(accession: str) -> Dict[str, str]:
    return {
        "rna_seqspec": first_matching_fallback(accession, ("rna", "scrna")),
        "sgrna_seqspec": first_matching_fallback(accession, ("guide", "grna", "sgrna")),
        "hash_seqspec": first_matching_fallback(accession, ("hash", "hto", "tag")),
    }


def fallback_seqspec_args(accession: str) -> List[str]:
    args: List[str] = []
    fallback_paths = fallback_seqspec_paths(accession)
    rna = fallback_paths["rna_seqspec"]
    guide = fallback_paths["sgrna_seqspec"]
    hash_seqspec = fallback_paths["hash_seqspec"]
    if rna:
        args.extend(["--rna-seqspec", rna])
    if guide:
        args.extend(["--sgrna-seqspec", guide])
    if hash_seqspec:
        args.extend(["--hash-seqspec", hash_seqspec])
    return args


def run_cli(
    cmd: Sequence[str],
    workdir: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> subprocess.CompletedProcess[str]:
    ensure_parent(stdout_path)
    ensure_parent(stderr_path)
    result = subprocess.run(
        list(cmd),
        cwd=str(workdir),
        text=True,
        capture_output=True,
        env=dict(os.environ),
    )
    stdout_path.write_text(result.stdout, encoding="utf-8")
    stderr_path.write_text(result.stderr, encoding="utf-8")
    return result


def sort_rows(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    return stable_sort_rows(rows, DEFAULT_SORT_COLUMNS)


def lane_key(row: Dict[str, str]) -> tuple[str, str]:
    return (field_or_empty(row, "sequencing_run"), field_or_empty(row, "lane"))


def sortable_value(value: str) -> tuple[int, object]:
    token = str(value or "").strip()
    if token.isdigit():
        return (0, int(token))
    return (1, token)


def select_shared_lane(rows: Sequence[Dict[str, str]]) -> Optional[SelectedSubset]:
    by_lane: Dict[tuple[str, str], Dict[str, object]] = {}
    for row in rows:
        key = lane_key(row)
        entry = by_lane.setdefault(key, {"rows": [], "modalities": set()})
        entry["rows"].append(row)
        entry["modalities"].add(field_or_empty(row, "file_modality"))
    complete_keys = [
        key for key, value in by_lane.items() if {"scRNA", "gRNA"}.issubset(value["modalities"])
    ]
    if not complete_keys:
        return None
    selected_key = min(complete_keys, key=lambda item: (sortable_value(item[0]), sortable_value(item[1])))
    selected_rows: List[Dict[str, str]] = []
    seen_modalities = set()
    for row in sort_rows(by_lane[selected_key]["rows"]):
        modality = field_or_empty(row, "file_modality")
        if modality in seen_modalities:
            continue
        seen_modalities.add(modality)
        selected_rows.append(dict(row))
    sequencing_run, lane = selected_key
    return SelectedSubset(
        mode="shared_lane",
        rows=selected_rows,
        note=f"Selected shared lane sequencing_run={sequencing_run}, lane={lane}.",
        samplesheet_path=None,
    )


def build_representative_subset(rows: Sequence[Dict[str, str]], output_path: Path) -> Optional[SelectedSubset]:
    selected: Dict[str, Dict[str, str]] = {}
    for row in sort_rows(rows):
        modality = field_or_empty(row, "file_modality")
        if modality not in {"scRNA", "gRNA", "hash"}:
            continue
        if modality not in selected:
            selected[modality] = dict(row)
    if not {"scRNA", "gRNA"}.issubset(selected):
        return None
    warning = (
        "Mixed lanes used intentionally because no shared scRNA+gRNA lane exists "
        "in this accession-derived samplesheet."
    )
    reduced_rows: List[Dict[str, str]] = []
    for modality in ("scRNA", "gRNA", "hash"):
        row = selected.get(modality)
        if not row:
            continue
        row["selection_group"] = "representative_mixed_modalities"
        row["selection_warning"] = warning
        reduced_rows.append(row)
    write_tsv(reduced_rows, output_path)
    pieces = []
    for row in reduced_rows:
        pieces.append(
            f"{field_or_empty(row, 'file_modality')}: measurement_sets={field_or_empty(row, 'measurement_sets')} "
            f"sequencing_run={field_or_empty(row, 'sequencing_run')} lane={field_or_empty(row, 'lane')}"
        )
    return SelectedSubset(
        mode="representative_mixed_modalities",
        rows=reduced_rows,
        note="Representative mixed-modality subset: " + "; ".join(pieces),
        samplesheet_path=output_path,
    )


def materialize_selected_subset(accession: str, accession_root: Path, subset: SelectedSubset) -> SelectedSubset:
    if subset.samplesheet_path is not None:
        return subset
    output_path = accession_root / f"{accession}_{subset.mode}_samplesheet.tsv"
    rows: List[Dict[str, str]] = []
    for row in subset.rows:
        item = dict(row)
        item["selection_group"] = subset.mode
        item["selection_warning"] = subset.note
        rows.append(item)
    write_tsv(rows, output_path)
    return SelectedSubset(
        mode=subset.mode,
        rows=rows,
        note=subset.note,
        samplesheet_path=output_path,
    )


def predictor_readiness(rows: Sequence[Dict[str, str]]) -> tuple[List[str], List[str]]:
    rows_by_modality = {field_or_empty(row, "file_modality"): row for row in rows}
    fatal: List[str] = []
    warnings: List[str] = []
    for modality in ("scRNA", "gRNA"):
        row = rows_by_modality.get(modality)
        if row is None:
            fatal.append(f"Missing required modality row: {modality}")
            continue
        for column in ("R1_path", "R2_path"):
            if not field_or_empty(row, column):
                fatal.append(f"{modality} is missing {column}")
    if not any(field_or_empty(row, "barcode_onlist") for row in rows):
        fatal.append("No barcode_onlist was present in the generated samplesheet subset")
    if not any(field_or_empty(row, "guide_design") for row in rows if field_or_empty(row, "file_modality") == "gRNA"):
        fatal.append("No guide_design was present for the gRNA subset")
    hash_row = rows_by_modality.get("hash")
    if hash_row is not None:
        for column in ("R1_path", "R2_path", "barcode_hashtag_map"):
            if not field_or_empty(hash_row, column):
                fatal.append(f"hash is missing {column}")
    for modality, row in rows_by_modality.items():
        if not field_or_empty(row, "seqspec"):
            warnings.append(f"{modality} is missing seqspec; prediction can run but comparison will be incomplete")
    return fatal, warnings


def locate_first(path_root: Path, pattern: str) -> Optional[Path]:
    matches = sorted(path_root.glob(pattern))
    return matches[0] if matches else None


def summarize_flags(comparison_rows: Sequence[Dict[str, object]]) -> str:
    if not comparison_rows:
        return ""
    counts = Counter(str(row.get("flag", "")) for row in comparison_rows)
    return "; ".join(f"{flag}={counts[flag]}" for flag in sorted(counts))


def extract_missing_seqspec_modalities(rows: Sequence[Dict[str, str]]) -> List[str]:
    missing = [field_or_empty(row, "file_modality") for row in rows if not field_or_empty(row, "seqspec")]
    return sorted(item for item in missing if item)


def classify_completed_run(
    selected_subset: SelectedSubset,
    analysis_summary: Dict[str, object],
) -> tuple[str, str, str]:
    comparison_rows = analysis_summary.get("comparison_rows", [])
    if not isinstance(comparison_rows, list):
        comparison_rows = []
    flags = [str(row.get("flag", "")) for row in comparison_rows if isinstance(row, dict)]
    missing_seqspec = extract_missing_seqspec_modalities(selected_subset.rows)
    if missing_seqspec:
        return (
            "samplesheet_incomplete",
            f"Prediction completed, but seqspec was missing for: {', '.join(missing_seqspec)}.",
            summarize_flags([row for row in comparison_rows if isinstance(row, dict)]),
        )
    if not comparison_rows:
        return (
            "samplesheet_incomplete",
            "Prediction completed, but no comparison rows were produced.",
            "",
        )
    red_flags = [flag for flag in flags if flag not in GOOD_FLAGS]
    if red_flags:
        return (
            "mismatch_needs_review",
            "Comparison contains one or more red flags.",
            summarize_flags([row for row in comparison_rows if isinstance(row, dict)]),
        )
    return (
        "good_no_red_flags",
        "All comparison rows were perfect_match or close_enough.",
        summarize_flags([row for row in comparison_rows if isinstance(row, dict)]),
    )


def tail_text(text: str, lines: int = 12) -> str:
    chunks = [line.rstrip() for line in text.splitlines() if line.strip()]
    if not chunks:
        return ""
    return " | ".join(chunks[-lines:])


def summarize_error_output(stdout: str, stderr: str) -> str:
    preferred_prefixes = (
        "RuntimeError:",
        "ValueError:",
        "requests.exceptions.",
        "urllib3.exceptions.",
        "ConnectionError:",
        "[ERROR]",
    )
    combined_lines = [line.strip() for line in (stderr + "\n" + stdout).splitlines() if line.strip()]
    for line in reversed(combined_lines):
        if any(line.startswith(prefix) for prefix in preferred_prefixes):
            return line
    return tail_text(stderr or stdout)


def write_outcome_json(path: Path, payload: Dict[str, object]) -> None:
    ensure_parent(path)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def build_report_row(
    source_row: Dict[str, str],
    lab: str,
    status: str,
    note: str,
    *,
    selection_mode: str = "",
    selection_note: str = "",
    flags_summary: str = "",
    samplesheet_path: str = "",
    analysis_summary_path: str = "",
    index_path: str = "",
    stdout_log: str = "",
    stderr_log: str = "",
) -> Dict[str, str]:
    analysis_set_path = field_or_empty(source_row, "ID")
    accession = extract_accession(analysis_set_path)
    return {
        "analysis_set_path": analysis_set_path,
        "analysis_set_accession": accession,
        "lab": lab,
        "status": status,
        "status_note": note,
        "selection_mode": selection_mode,
        "selection_note": selection_note,
        "comparison_flags": flags_summary,
        "summary": field_or_empty(source_row, "Summary"),
        "simplified_sample_summary": field_or_empty(source_row, "Simplified Sample Summary"),
        "sample_summary": field_or_empty(source_row, "Sample Summary"),
        "file_set_type": field_or_empty(source_row, "File Set Type"),
        "file_content_type": field_or_empty(source_row, "File Content Type"),
        "samplesheet_path": samplesheet_path,
        "analysis_summary_path": analysis_summary_path,
        "index_path": index_path,
        "stdout_log": stdout_log,
        "stderr_log": stderr_log,
    }


def print_progress(message: str) -> None:
    print(message, flush=True)


def prepare_analysis_root(analysis_root: Path) -> None:
    if not analysis_root.exists():
        analysis_root.mkdir(parents=True, exist_ok=True)
        return
    for child in analysis_root.iterdir():
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()


def make_accession_progress(idx: int, total: int, accession: str):
    prefix = f"[{idx}/{total}] {accession}: "

    def emit(message: str) -> None:
        print_progress(prefix + message)

    return emit


def summarize_modalities(rows: Sequence[Dict[str, str]]) -> str:
    counts = Counter(field_or_empty(row, "file_modality") for row in rows)
    return ", ".join(f"{modality}={counts[modality]}" for modality in sorted(counts) if modality)


def summarize_subset_rows(rows: Sequence[Dict[str, str]]) -> str:
    parts: List[str] = []
    for row in rows:
        modality = field_or_empty(row, "file_modality")
        measurement_sets = field_or_empty(row, "measurement_sets")
        sequencing_run = field_or_empty(row, "sequencing_run") or "NA"
        lane = field_or_empty(row, "lane") or "NA"
        piece = f"{modality} run={sequencing_run} lane={lane}"
        if measurement_sets:
            piece += f" measurement_sets={measurement_sets}"
        parts.append(piece)
    return "; ".join(parts)


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent
    input_report = (repo_root / args.input_report).resolve()
    analysis_root = (repo_root / args.analysis_root).resolve()
    if args.prepare_only and not args.fresh_run:
        raise ValueError("--prepare-only requires --fresh-run.")
    if args.fresh_run:
        print_progress(f"Preparing fresh analysis root: {analysis_root}")
        prepare_analysis_root(analysis_root)
    analysis_root.mkdir(parents=True, exist_ok=True)
    if args.prepare_only:
        print_progress(f"Prepared {analysis_root} for a new run.")
        return

    credentials = load_credentials(keypair_path=args.igvf_keypair, required=True)
    report_rows = load_multireport_rows(input_report)
    if args.accession:
        wanted = {extract_accession(item) for item in args.accession}
        report_rows = [row for row in report_rows if extract_accession(field_or_empty(row, "ID")) in wanted]
    if args.limit is not None:
        report_rows = report_rows[: args.limit]
    if not report_rows:
        raise ValueError("No analysis sets matched the requested filters.")

    aggregate_rows: List[Dict[str, str]] = []
    total = len(report_rows)
    aggregate_json_path = analysis_root / "report.json"
    aggregate_tsv_path = analysis_root / "report.tsv"

    for idx, source_row in enumerate(report_rows, start=1):
        analysis_set_path = field_or_empty(source_row, "ID")
        accession = extract_accession(analysis_set_path)
        emit = make_accession_progress(idx, total, accession)
        accession_root = analysis_root / accession
        accession_root.mkdir(parents=True, exist_ok=True)
        outcome_json = accession_root / "batch_outcome.json"

        if outcome_json.exists() and not args.force:
            payload = json.loads(outcome_json.read_text(encoding="utf-8"))
            aggregate_rows.append(payload["report_row"])
            emit(f"reused existing outcome ({payload['report_row']['status']})")
            continue

        lab = try_fetch_lab(accession, credentials)
        logs_dir = accession_root / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        cli_base = [
            sys.executable,
            str(repo_root / "seqspec_parser.py"),
        ]
        fallback_paths = fallback_seqspec_paths(accession)
        shared_samplesheet_path = accession_root / f"{sanitize_label(accession)}_samplesheet.tsv"

        emit("generating portal samplesheet (fetching portal metadata, file sets, barcodes, guide metadata, and seqspec links)")
        samplesheet_cmd: List[str] = ["local_generate_samplesheet_rows", accession]
        samplesheet_stdout_path = logs_dir / "samplesheet_only.stdout.log"
        samplesheet_stderr_path = logs_dir / "samplesheet_only.stderr.log"
        try:
            generated_rows = generate_samplesheet_rows(
                accession,
                credentials,
                hash_seqspec=fallback_paths["hash_seqspec"] or None,
                rna_seqspec=fallback_paths["rna_seqspec"] or None,
                sgrna_seqspec=fallback_paths["sgrna_seqspec"] or None,
                stop_after_first_complete_measurement_set=True,
                progress=lambda message: emit(f"portal samplesheet: {message}"),
            )
            write_tsv(generated_rows, shared_samplesheet_path)
            samplesheet_stdout_path.write_text(
                f"Generated {len(generated_rows)} rows via local portal generator.\n",
                encoding="utf-8",
            )
            samplesheet_stderr_path.write_text("", encoding="utf-8")
            samplesheet_error = ""
            emit(f"portal samplesheet ready: {len(generated_rows)} rows ({summarize_modalities(generated_rows)})")
        except Exception as exc:
            samplesheet_stdout_path.write_text("", encoding="utf-8")
            samplesheet_stderr_path.write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")
            samplesheet_error = f"{type(exc).__name__}: {exc}"
        if samplesheet_error or not shared_samplesheet_path.exists():
            note = samplesheet_error or "Samplesheet generation failed."
            report_row = build_report_row(
                source_row,
                lab,
                "samplesheet_extraction_failed",
                note,
                samplesheet_path=str(shared_samplesheet_path) if shared_samplesheet_path.exists() else "",
                stdout_log=str(samplesheet_stdout_path.resolve()),
                stderr_log=str(samplesheet_stderr_path.resolve()),
            )
            payload = {
                "report_row": report_row,
                "samplesheet_cmd": samplesheet_cmd,
                "analysis_cmd": None,
            }
            write_outcome_json(outcome_json, payload)
            aggregate_rows.append(report_row)
            write_tsv(aggregate_rows, aggregate_tsv_path)
            write_outcome_json(aggregate_json_path, {"rows": aggregate_rows})
            emit("samplesheet extraction failed")
            continue

        generated_rows = canonicalize_samplesheet_rows(load_rows(str(shared_samplesheet_path)))
        emit("selecting one analyzable subset (shared scRNA+gRNA lane first, otherwise one representative pair per modality)")
        selected_subset = select_shared_lane(generated_rows)
        if selected_subset is None:
            mixed_path = accession_root / f"{accession}_representative_mixed_modalities_samplesheet.tsv"
            selected_subset = build_representative_subset(generated_rows, mixed_path)
        if selected_subset is None:
            note = "No shared scRNA+gRNA lane and no representative scRNA+gRNA subset could be built."
            report_row = build_report_row(
                source_row,
                lab,
                "samplesheet_incomplete",
                note,
                samplesheet_path=str(shared_samplesheet_path),
                stdout_log=str(samplesheet_stdout_path.resolve()),
                stderr_log=str(samplesheet_stderr_path.resolve()),
            )
            payload = {
                "report_row": report_row,
                "samplesheet_cmd": samplesheet_cmd,
                "analysis_cmd": None,
            }
            write_outcome_json(outcome_json, payload)
            aggregate_rows.append(report_row)
            write_tsv(aggregate_rows, aggregate_tsv_path)
            write_outcome_json(aggregate_json_path, {"rows": aggregate_rows})
            emit("no analyzable scRNA+gRNA subset")
            continue

        selected_subset = materialize_selected_subset(accession, accession_root, selected_subset)
        emit(f"selected {selected_subset.mode}: {summarize_subset_rows(selected_subset.rows)}")
        fatal_issues, warnings = predictor_readiness(selected_subset.rows)
        if fatal_issues:
            note = "; ".join(fatal_issues + warnings)
            report_row = build_report_row(
                source_row,
                lab,
                "samplesheet_incomplete",
                note,
                selection_mode=selected_subset.mode,
                selection_note=selected_subset.note,
                samplesheet_path=str(selected_subset.samplesheet_path or shared_samplesheet_path),
                stdout_log=str(samplesheet_stdout_path.resolve()),
                stderr_log=str(samplesheet_stderr_path.resolve()),
            )
            payload = {
                "report_row": report_row,
                "samplesheet_cmd": samplesheet_cmd,
                "analysis_cmd": None,
            }
            write_outcome_json(outcome_json, payload)
            aggregate_rows.append(report_row)
            write_tsv(aggregate_rows, aggregate_tsv_path)
            write_outcome_json(aggregate_json_path, {"rows": aggregate_rows})
            emit("samplesheet incomplete for prediction")
            continue

        analysis_cmd = [
            *cli_base,
            "samplesheet",
            "--samplesheet",
            str(selected_subset.samplesheet_path),
            "--analysis-root",
            str(accession_root),
            "--igvf-keypair",
            args.igvf_keypair,
            "--group-by",
            "selection_group",
            "--barcode-sample-reads",
            str(args.barcode_sample_reads),
            "--feature-sample-reads",
            str(args.feature_sample_reads),
        ]
        if args.ignore_metadata_md5:
            analysis_cmd.append("--ignore-metadata-md5")

        emit(
            f"running {selected_subset.mode} analysis "
            f"(staging first-pair FASTQ chunks and metadata, then running predictor/comparison)"
        )
        analysis_result = run_cli(
            analysis_cmd,
            repo_root,
            logs_dir / "analysis.stdout.log",
            logs_dir / "analysis.stderr.log",
        )
        analysis_summary_path = locate_first(accession_root, "**/analysis_summary.json")
        index_path = accession_root / "index.html"

        if analysis_result.returncode != 0 or analysis_summary_path is None:
            note_parts = [selected_subset.note]
            if warnings:
                note_parts.extend(warnings)
            error_tail = summarize_error_output(analysis_result.stdout, analysis_result.stderr)
            if error_tail:
                note_parts.append(error_tail)
            report_row = build_report_row(
                source_row,
                lab,
                "analysis_failed",
                "; ".join(part for part in note_parts if part),
                selection_mode=selected_subset.mode,
                selection_note=selected_subset.note,
                samplesheet_path=str(selected_subset.samplesheet_path or shared_samplesheet_path),
                index_path=str(index_path) if index_path.exists() else "",
                stdout_log=str((logs_dir / "analysis.stdout.log").resolve()),
                stderr_log=str((logs_dir / "analysis.stderr.log").resolve()),
            )
            payload = {
                "report_row": report_row,
                "samplesheet_cmd": samplesheet_cmd,
                "analysis_cmd": analysis_cmd,
            }
            write_outcome_json(outcome_json, payload)
            aggregate_rows.append(report_row)
            write_tsv(aggregate_rows, aggregate_tsv_path)
            write_outcome_json(aggregate_json_path, {"rows": aggregate_rows})
            emit("analysis failed")
            continue

        analysis_summary = json.loads(analysis_summary_path.read_text(encoding="utf-8"))
        status, note, flags_summary = classify_completed_run(selected_subset, analysis_summary)
        note_parts = [selected_subset.note]
        if warnings:
            note_parts.extend(warnings)
        if note:
            note_parts.append(note)
        report_row = build_report_row(
            source_row,
            lab,
            status,
            "; ".join(part for part in note_parts if part),
            selection_mode=selected_subset.mode,
            selection_note=selected_subset.note,
            flags_summary=flags_summary,
            samplesheet_path=str(selected_subset.samplesheet_path or shared_samplesheet_path),
            analysis_summary_path=str(analysis_summary_path.resolve()),
            index_path=str(index_path.resolve()) if index_path.exists() else "",
            stdout_log=str((logs_dir / "analysis.stdout.log").resolve()),
            stderr_log=str((logs_dir / "analysis.stderr.log").resolve()),
        )
        payload = {
            "report_row": report_row,
            "samplesheet_cmd": samplesheet_cmd,
            "analysis_cmd": analysis_cmd,
        }
        write_outcome_json(outcome_json, payload)
        aggregate_rows.append(report_row)
        write_tsv(aggregate_rows, aggregate_tsv_path)
        write_outcome_json(aggregate_json_path, {"rows": aggregate_rows})
        emit(status)

    counts = Counter(row["status"] for row in aggregate_rows)
    summary_path = analysis_root / "summary.json"
    write_outcome_json(summary_path, {"status_counts": dict(counts), "rows": aggregate_rows})
    print_progress(
        "Completed batch run. "
        + ", ".join(f"{status}={counts[status]}" for status in sorted(counts))
        + f". Report: {aggregate_tsv_path}"
    )


if __name__ == "__main__":
    main()
