#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path
from typing import Dict, List, Sequence


DEFAULT_CHUNK_BYTES = 10 * 1024 * 1024
DEFAULT_METADATA_TAGS = {
    "barcode_onlist": "barcodeWhitelist",
    "barcode_hashtag_map": "hashMetadata",
    "guide_design": "guideMetadata",
    "seqspec": "seqspec",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a reusable shell script that downloads representative FASTQ chunks and metadata "
            "from an IGVF-style sample sheet."
        )
    )
    parser.add_argument("--samplesheet", required=True, help="CSV or TSV sample sheet")
    parser.add_argument("--output-script", default="download_test_data.sh", help="Shell script to write")
    parser.add_argument("--output-dir", default="test_data", help="Download destination directory")
    parser.add_argument(
        "--modalities",
        help="Comma-separated modality filter, for example scRNA,gRNA,hash",
    )
    parser.add_argument(
        "--filter",
        action="append",
        default=[],
        help="Repeatable COLUMN=VALUE row filter applied before grouping",
    )
    parser.add_argument(
        "--group-by",
        default="file_modality",
        help="Comma-separated columns used to collapse the sheet into one download group per key",
    )
    parser.add_argument(
        "--sort-by",
        default="file_modality,measurement_sets,sequencing_run,lane",
        help="Comma-separated columns used to choose the first row within each group",
    )
    parser.add_argument(
        "--chunk-bytes",
        type=int,
        default=DEFAULT_CHUNK_BYTES,
        help="If > 0, fetch only this many bytes from each FASTQ via gsutil cat -r; use 0 for full gsutil cp",
    )
    parser.add_argument(
        "--credentials",
        help="Optional GOOGLE_APPLICATION_CREDENTIALS path to export in the generated script",
    )
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="Print the generated shell script to stdout in addition to writing it",
    )
    return parser.parse_args()


def sniff_delimiter(path: str) -> str:
    return "\t" if path.endswith(".tsv") else ","


def load_rows(path: str) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=sniff_delimiter(path))
        if reader.fieldnames is None:
            raise ValueError(f"Could not read headers from {path}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def parse_csv_list(value: str | None) -> List[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_filters(raw_filters: Sequence[str]) -> List[tuple[str, str]]:
    filters: List[tuple[str, str]] = []
    for item in raw_filters:
        if "=" not in item:
            raise ValueError(f"Invalid filter '{item}'. Expected COLUMN=VALUE.")
        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            raise ValueError(f"Invalid filter '{item}'. Missing column name.")
        filters.append((key, value))
    return filters


def apply_filters(
    rows: Sequence[Dict[str, str]],
    modalities: Sequence[str],
    filters: Sequence[tuple[str, str]],
) -> List[Dict[str, str]]:
    selected = list(rows)
    if modalities:
        allowed = set(modalities)
        selected = [row for row in selected if row.get("file_modality", "") in allowed]
    for column, value in filters:
        selected = [row for row in selected if row.get(column, "") == value]
    return selected


def stable_sort_rows(rows: Sequence[Dict[str, str]], sort_columns: Sequence[str]) -> List[Dict[str, str]]:
    return sorted(rows, key=lambda row: tuple(row.get(column, "") for column in sort_columns))


def sanitize_label(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("_") or "sample"


def build_group_label(row: Dict[str, str], group_columns: Sequence[str]) -> str:
    parts = [sanitize_label(row.get(column, "")) for column in group_columns]
    return "_".join(part for part in parts if part)


def first_nonempty_extension(path: str) -> str:
    suffixes = Path(path).suffixes
    if not suffixes:
        return ""
    return "".join(suffixes)


def is_remote_path(value: str) -> bool:
    return bool(value) and value.startswith("gs://")


def chunk_suffix(chunk_bytes: int) -> str:
    if chunk_bytes <= 0:
        return ""
    if chunk_bytes % (1024 * 1024) == 0:
        return f"_{chunk_bytes // (1024 * 1024)}mb"
    if chunk_bytes % 1024 == 0:
        return f"_{chunk_bytes // 1024}kb"
    return f"_{chunk_bytes}b"


def build_fastq_command(remote_path: str, local_path: str, chunk_bytes: int) -> str:
    if chunk_bytes > 0:
        return (
            'gsutil -o "GSUtil:parallel_process_count=1" '
            f'cat -r 0-{chunk_bytes} "{remote_path}" > "{local_path}"'
        )
    return f'gsutil -o "GSUtil:parallel_process_count=1" cp "{remote_path}" "{local_path}"'


def build_metadata_commands(row: Dict[str, str], label: str, output_dir: str) -> List[str]:
    commands: List[str] = []
    for column, tag in DEFAULT_METADATA_TAGS.items():
        remote_path = row.get(column, "")
        if not is_remote_path(remote_path):
            continue
        extension = first_nonempty_extension(remote_path)
        local_path = os.path.join(output_dir, f"{label}_{tag}{extension}")
        commands.append(f"echo 'Downloading {tag} for {label}'")
        commands.append(
            'gsutil -o "GSUtil:parallel_process_count=1" '
            f'cp "{remote_path}" "{local_path}"'
        )
    return commands


def collapse_rows(rows: Sequence[Dict[str, str]], group_columns: Sequence[str]) -> List[tuple[str, Dict[str, str]]]:
    grouped: Dict[tuple[str, ...], Dict[str, str]] = {}
    for row in rows:
        key = tuple(row.get(column, "") for column in group_columns)
        if key not in grouped:
            grouped[key] = row
    collapsed: List[tuple[str, Dict[str, str]]] = []
    for key, row in grouped.items():
        label = build_group_label(row, group_columns)
        collapsed.append((label, row))
    return collapsed


def render_script(
    selected_rows: Sequence[tuple[str, Dict[str, str]]],
    output_dir: str,
    chunk_bytes: int,
    credentials: str | None,
) -> str:
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "# Auto-generated by downloading_from_samplesheet.py",
    ]

    if credentials:
        lines.append(f'export GOOGLE_APPLICATION_CREDENTIALS="{credentials}"')

    lines.append(f'mkdir -p "{output_dir}"')
    lines.append("")

    if chunk_bytes > 0:
        lines.append(f"# FASTQ files are fetched as range chunks ({chunk_bytes} bytes).")
        lines.append("# For .gz inputs this creates truncated gzip members by design.")
    else:
        lines.append("# FASTQ files are copied in full with gsutil cp.")
    lines.append("")

    suffix = chunk_suffix(chunk_bytes)

    for label, row in selected_rows:
        modality = row.get("file_modality", label)
        read1 = row.get("R1_path", "")
        read2 = row.get("R2_path", "")
        if not (is_remote_path(read1) and is_remote_path(read2)):
            raise ValueError(f"Row for {label} is missing a valid gs:// R1_path or R2_path")

        read1_ext = first_nonempty_extension(read1)
        read2_ext = first_nonempty_extension(read2)
        local_r1 = os.path.join(output_dir, f"{label}_R1{suffix}{read1_ext}")
        local_r2 = os.path.join(output_dir, f"{label}_R2{suffix}{read2_ext}")

        lines.append(f"### {label} ({modality}) ###")
        measurement_set = row.get("measurement_sets", "")
        sequencing_run = row.get("sequencing_run", "")
        lane = row.get("lane", "")
        summary = ", ".join(part for part in [
            f"measurement_sets={measurement_set}" if measurement_set else "",
            f"sequencing_run={sequencing_run}" if sequencing_run else "",
            f"lane={lane}" if lane else "",
        ] if part)
        if summary:
            lines.append(f"echo 'Preparing {label}: {summary}'")
        else:
            lines.append(f"echo 'Preparing {label}'")
        lines.append(build_fastq_command(read1, local_r1, chunk_bytes))
        lines.append(build_fastq_command(read2, local_r2, chunk_bytes))
        lines.extend(build_metadata_commands(row, label, output_dir))
        lines.append("")

    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    modalities = parse_csv_list(args.modalities)
    group_columns = parse_csv_list(args.group_by)
    sort_columns = parse_csv_list(args.sort_by)
    if not group_columns:
        raise ValueError("--group-by must contain at least one column")

    rows = load_rows(args.samplesheet)
    filtered = apply_filters(rows, modalities, parse_filters(args.filter))
    if not filtered:
        raise ValueError("No rows remained after applying modality/filter selection.")

    sorted_rows = stable_sort_rows(filtered, sort_columns)
    selected_rows = collapse_rows(sorted_rows, group_columns)
    script_text = render_script(
        selected_rows=selected_rows,
        output_dir=args.output_dir,
        chunk_bytes=args.chunk_bytes,
        credentials=args.credentials,
    )

    with open(args.output_script, "w", encoding="utf-8") as handle:
        handle.write(script_text)
    os.chmod(args.output_script, 0o755)

    if args.stdout:
        print(script_text, end="")

    print(f"Wrote {args.output_script} with {len(selected_rows)} download group(s).")
    for label, row in selected_rows:
        modality = row.get("file_modality", "")
        print(f"  - {label} ({modality})")


if __name__ == "__main__":
    main()
