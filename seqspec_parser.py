#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import yaml

from downloading_from_samplesheet import (
    DEFAULT_CHUNK_BYTES,
    apply_filters,
    build_group_label,
    first_nonempty_extension,
    is_remote_path as is_gcs_path,
    load_rows,
    parse_csv_list,
    parse_filters,
    sanitize_label,
    stable_sort_rows,
)
from igvf_portal import (
    IGVFCredentials,
    download_file as download_igvf_file,
    generate_samplesheet_rows,
    load_credentials as load_igvf_credentials,
    write_samplesheet as write_igvf_samplesheet,
)


class TagIgnoringLoader(yaml.SafeLoader):
    pass


def _ignore_yaml_tag(loader, tag_suffix, node):
    if isinstance(node, yaml.MappingNode):
        return loader.construct_mapping(node)
    if isinstance(node, yaml.SequenceNode):
        return loader.construct_sequence(node)
    return loader.construct_scalar(node)


TagIgnoringLoader.add_multi_constructor("!", _ignore_yaml_tag)


MODALITY_TO_FEATURE = {
    "scrna": "rna",
    "rna": "rna",
    "grna": "guide",
    "guide": "guide",
    "crispr": "guide",
    "hash": "hash",
    "tag": "hash",
    "hto": "hash",
    "hashing": "hash",
}

MODALITY_TO_SEQSPEC_INDEX = {
    "scrna": "rna",
    "rna": "rna",
    "grna": "crispr",
    "guide": "crispr",
    "crispr": "crispr",
    "hash": "tag",
    "tag": "tag",
    "hto": "tag",
    "hashing": "tag",
}

METADATA_COLUMNS = {
    "barcode_onlist": "barcodeWhitelist",
    "guide_design": "guideMetadata",
    "barcode_hashtag_map": "hashMetadata",
}

PORTAL_MODALITY_CANONICAL = {
    "scrna sequencing": "scRNA",
    "grna sequencing": "gRNA",
    "cell hashing barcode sequencing": "hash",
}

IGVF_FILE_TYPES = {
    "sequence": ("sequence-files", ".fastq.gz"),
    "configuration": ("configuration-files", ".yaml.gz"),
    "tabular": ("tabular-files", ".tsv.gz"),
}

FLAG_ORDER = {
    "perfect_match": 0,
    "close_enough": 1,
    "very_distant": 2,
    "missing_prediction": 3,
    "missing_seqspec": 4,
}

DNA_ALPHABET = set("ACGTN")


@dataclass
class ReadSpec:
    file_id: str
    filename: str
    read_id: str
    index_id: str
    read_label: str
    strand: str
    min_len: int
    max_len: int
    name: str


@dataclass
class IndexedRegion:
    file_id: str
    read_label: str
    strand: str
    region_id: str
    region_type: str
    canonical_name: Optional[str]
    raw_start: int
    raw_end: int
    length: int
    sequence_type: str
    template_sequence: str
    declared_min_len: int
    declared_max_len: int


@dataclass
class SeqspecSummary:
    raw_modality: str
    seqspec_modality: str
    seqspec_path: str
    kb_string: str
    read_specs: List[ReadSpec]
    indexed_regions: List[IndexedRegion]


@dataclass
class FeatureMetadataSummary:
    source_path: str
    sequence_column: str
    total_sequences: int
    length_counts: List[tuple[int, int]]


@dataclass
class ComparisonRow:
    modality: str
    region: str
    seqspec_interval: Optional[str]
    prediction_interval: Optional[str]
    seqspec_read: Optional[str]
    prediction_read: Optional[str]
    seqspec_strand: Optional[str]
    prediction_strand: Optional[str]
    strand_flag: str
    max_distance_bp: Optional[int]
    flag: str
    note: str


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--analysis-root", default="seqspec_analysis", help="Output directory for analysis runs")
    parser.add_argument(
        "--downloads-dir",
        help=(
            "Optional existing download/cache directory containing test_data-style modality files. "
            "Use this to avoid re-downloading data that is already present. In single-lane mode, "
            "explicit files are staged there if needed."
        ),
    )
    parser.add_argument(
        "--chunk-bytes",
        type=int,
        default=DEFAULT_CHUNK_BYTES,
        help="Chunk size for FASTQ downloads; use 0 for full remote file downloads",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Stage FASTQ chunks, seqspecs, and metadata, then stop before seqspec/prediction/report steps",
    )
    parser.add_argument("--force-download", action="store_true", help="Re-download assets even if local copies exist")
    parser.add_argument("--seqspec-bin", default="seqspec", help="Path to the seqspec executable")
    parser.add_argument("--barcode-sample-reads", type=int, default=10_000)
    parser.add_argument("--feature-sample-reads", type=int, default=100_000)
    parser.add_argument("--guide-upstream-bases", type=int, default=12)
    parser.add_argument("--igvf-keypair", help="JSON file containing IGVF API credentials for portal-backed downloads")
    parser.add_argument("--igvf-api-key", help="IGVF API key used for portal-backed downloads")
    parser.add_argument("--igvf-api-secret", help="IGVF API secret used for portal-backed downloads")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Unified seqspec validation and predictor comparison workflow. "
            "Downloads only missing assets, runs seqspec index, runs the predictor, "
            "and writes seqspec/prediction/comparison HTML outputs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Modes:\n"
            "  samplesheet   Run grouped analysis from a metadata table.\n"
            "  igvf          Generate a validator-ready samplesheet from an IGVF analysis set, then optionally validate it.\n"
            "  single-lane   Run one analysis directly from explicit FASTQs, seqspecs, and metadata files.\n\n"
            "Examples:\n"
            "  seqspec_parser.py samplesheet --samplesheet sample_metadata.csv --analysis-root out --filter sequencing_run=2 --filter lane=1\n"
            "  seqspec_parser.py igvf --accession IGVFDS9445RJOU --analysis-root out --samplesheet-only\n"
            "  seqspec_parser.py igvf --accession IGVFDS9445RJOU --analysis-root out --igvf-keypair igvf_key.json --one-lane\n"
            "  seqspec_parser.py single-lane --analysis-root out --group-label lane1 "
            "--scrna-r1 scRNA_R1.fastq.gz --scrna-r2 scRNA_R2.fastq.gz --scrna-seqspec scRNA.yaml "
            "--guide-r1 gRNA_R1.fastq.gz --guide-r2 gRNA_R2.fastq.gz --guide-seqspec gRNA.yaml "
            "--barcode-whitelist whitelist.tsv --guide-metadata guides.tsv"
        ),
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    samplesheet_parser = subparsers.add_parser(
        "samplesheet",
        help="Run analysis from a multi-modality samplesheet",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    samplesheet_parser.add_argument("--samplesheet", required=True, help="CSV or TSV sample sheet containing seqspec paths")
    samplesheet_parser.add_argument(
        "--group-by",
        default="sequencing_run,lane",
        help="Columns used to define one multi-modality analysis group",
    )
    samplesheet_parser.add_argument(
        "--sort-by",
        default="file_modality,measurement_sets,sequencing_run,lane",
        help="Columns used to deterministically pick one row per modality inside each group",
    )
    samplesheet_parser.add_argument("--modalities", default="scRNA,gRNA,hash", help="Comma-separated modality filter")
    samplesheet_parser.add_argument("--filter", action="append", default=[], help="Repeatable COLUMN=VALUE filter")
    add_common_args(samplesheet_parser)

    igvf_parser = subparsers.add_parser(
        "igvf",
        help="Generate a validator samplesheet from an IGVF analysis set and optionally validate it",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    igvf_parser.add_argument("--accession", required=True, help="IGVF analysis-set accession, for example IGVFDS9445RJOU")
    igvf_parser.add_argument(
        "--samplesheet-out",
        help="Path for the generated validator-compatible samplesheet TSV (defaults to <analysis-root>/<accession>_samplesheet.tsv)",
    )
    igvf_parser.add_argument("--samplesheet-only", action="store_true", help="Write the IGVF-derived samplesheet and exit")
    igvf_parser.add_argument(
        "--one-lane",
        action="store_true",
        help="After generating the IGVF samplesheet, process only the first complete sequencing_run/lane group",
    )
    igvf_parser.add_argument(
        "--group-by",
        default="sequencing_run,lane",
        help="Columns used to define one multi-modality analysis group after the samplesheet is generated",
    )
    igvf_parser.add_argument(
        "--sort-by",
        default="file_modality,measurement_sets,sequencing_run,lane",
        help="Columns used to deterministically pick one row per modality inside each group",
    )
    igvf_parser.add_argument("--modalities", default="scRNA,gRNA,hash", help="Comma-separated modality filter")
    igvf_parser.add_argument("--filter", action="append", default=[], help="Repeatable COLUMN=VALUE filter")
    igvf_parser.add_argument("--rna-seqspec", help="Fallback scRNA seqspec path if the portal record lacks one")
    igvf_parser.add_argument("--sgrna-seqspec", help="Fallback gRNA seqspec path if the portal record lacks one")
    igvf_parser.add_argument("--hash-seqspec", help="Fallback hash seqspec path if the portal record lacks one")
    add_common_args(igvf_parser)

    single_lane_parser = subparsers.add_parser(
        "single-lane",
        help="Run one analysis from explicit modality files without a samplesheet",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    single_lane_parser.add_argument("--group-label", default="single_lane", help="Label used for the output group directory")
    single_lane_parser.add_argument("--scrna-r1", required=True, help="scRNA read 1 FASTQ path")
    single_lane_parser.add_argument("--scrna-r2", required=True, help="scRNA read 2 FASTQ path")
    single_lane_parser.add_argument("--scrna-seqspec", required=True, help="scRNA seqspec YAML path")
    single_lane_parser.add_argument("--guide-r1", required=True, help="Guide read 1 FASTQ path")
    single_lane_parser.add_argument("--guide-r2", required=True, help="Guide read 2 FASTQ path")
    single_lane_parser.add_argument("--guide-seqspec", required=True, help="Guide seqspec YAML path")
    single_lane_parser.add_argument(
        "--barcode-whitelist",
        "--barcodes",
        dest="barcode_whitelist",
        required=True,
        help="Cell barcode whitelist TSV/CSV path",
    )
    single_lane_parser.add_argument("--guide-metadata", required=True, help="Guide metadata TSV/CSV path")
    single_lane_parser.add_argument("--hash-r1", help="Optional hash read 1 FASTQ path")
    single_lane_parser.add_argument("--hash-r2", help="Optional hash read 2 FASTQ path")
    single_lane_parser.add_argument("--hash-seqspec", help="Optional hash seqspec YAML path")
    single_lane_parser.add_argument("--hash-metadata", help="Optional hash metadata TSV/CSV path")
    add_common_args(single_lane_parser)

    return parser


def normalize_cli_argv(argv: Optional[Sequence[str]] = None) -> List[str]:
    tokens = list(sys.argv[1:] if argv is None else argv)
    if not tokens:
        return tokens
    if tokens[0] in {"samplesheet", "igvf", "single-lane", "-h", "--help"}:
        return tokens
    if "--samplesheet" in tokens:
        return ["samplesheet", *tokens]
    single_lane_markers = {"--scrna-r1", "--guide-r1", "--barcode-whitelist", "--barcodes"}
    if any(marker in tokens for marker in single_lane_markers):
        return ["single-lane", *tokens]
    return tokens


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    return build_parser().parse_args(normalize_cli_argv(argv))


def normalize_modality(raw_modality: str) -> str:
    key = raw_modality.strip().lower()
    return MODALITY_TO_FEATURE.get(key, key)


def canonicalize_samplesheet_modality(raw_modality: str) -> str:
    token = str(raw_modality or "").strip()
    lowered = token.lower()
    if lowered in PORTAL_MODALITY_CANONICAL:
        return PORTAL_MODALITY_CANONICAL[lowered]
    return token


def canonicalize_samplesheet_rows(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    normalized: List[Dict[str, str]] = []
    for row in rows:
        item = dict(row)
        item["file_modality"] = canonicalize_samplesheet_modality(item.get("file_modality", ""))
        normalized.append(item)
    return normalized


def normalize_seqspec_index_modality(raw_modality: str) -> str:
    key = raw_modality.strip().lower()
    return MODALITY_TO_SEQSPEC_INDEX.get(key, key)


def normalize_strand(value: str) -> str:
    token = str(value or "").strip().lower()
    if token in {"neg", "-", "reverse", "rev"}:
        return "reverse"
    return "forward"


def group_rows_by_analysis(
    rows: Sequence[Dict[str, str]],
    group_columns: Sequence[str],
) -> List[Dict[str, object]]:
    groups: Dict[tuple[str, ...], Dict[str, object]] = {}
    for row in rows:
        key = tuple(row.get(column, "") for column in group_columns)
        if key not in groups:
            label_seed = {column: row.get(column, "") for column in group_columns}
            groups[key] = {
                "key": key,
                "label": build_group_label(label_seed, group_columns),
                "rows_by_modality": {},
            }
        modality = row.get("file_modality", "")
        if modality not in groups[key]["rows_by_modality"]:
            groups[key]["rows_by_modality"][modality] = row
    return list(groups.values())


def validate_single_lane_args(args: argparse.Namespace) -> None:
    hash_values = [args.hash_r1, args.hash_r2, args.hash_seqspec, args.hash_metadata]
    if any(hash_values) and not all(hash_values):
        raise ValueError(
            "Hash mode is optional, but if used you must provide --hash-r1, --hash-r2, --hash-seqspec, and --hash-metadata together."
        )


def build_single_lane_group(args: argparse.Namespace) -> Dict[str, object]:
    validate_single_lane_args(args)
    rows_by_modality: Dict[str, Dict[str, str]] = {
        "scRNA": {
            "file_modality": "scRNA",
            "R1_path": args.scrna_r1,
            "R2_path": args.scrna_r2,
            "seqspec": args.scrna_seqspec,
            "barcode_onlist": args.barcode_whitelist,
            "guide_design": args.guide_metadata,
            "barcode_hashtag_map": "",
        },
        "gRNA": {
            "file_modality": "gRNA",
            "R1_path": args.guide_r1,
            "R2_path": args.guide_r2,
            "seqspec": args.guide_seqspec,
            "barcode_onlist": args.barcode_whitelist,
            "guide_design": args.guide_metadata,
            "barcode_hashtag_map": "",
        },
    }
    if args.hash_r1:
        rows_by_modality["hash"] = {
            "file_modality": "hash",
            "R1_path": args.hash_r1,
            "R2_path": args.hash_r2,
            "seqspec": args.hash_seqspec,
            "barcode_onlist": args.barcode_whitelist,
            "guide_design": "",
            "barcode_hashtag_map": args.hash_metadata,
        }
    return {
        "key": (args.group_label,),
        "label": args.group_label,
        "rows_by_modality": rows_by_modality,
    }


def ensure_seqspec_cli(seqspec_bin: str, env: Dict[str, str]) -> None:
    result = subprocess.run(
        [seqspec_bin, "index", "--help"],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Could not run '{seqspec_bin} index --help'. stderr:\n{result.stderr}")


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def is_http_path(value: str) -> bool:
    token = str(value or "").strip().lower()
    return token.startswith("http://") or token.startswith("https://")


def is_igvf_accession(value: str) -> bool:
    token = str(value or "").strip().upper()
    return token.startswith("IGVFFI") and token[6:].isalnum()


def is_remote_asset(value: str) -> bool:
    return is_gcs_path(value) or is_http_path(value) or is_igvf_accession(value)


def igvf_download_url(accession: str, file_type: str) -> str:
    route, extension = IGVF_FILE_TYPES[file_type]
    clean_accession = accession.strip()
    return f"https://api.data.igvf.org/{route}/{clean_accession}/@@download/{clean_accession}{extension}"


def asset_extension(remote_or_local: str, igvf_file_type: Optional[str]) -> str:
    extension = first_nonempty_extension(remote_or_local)
    if extension:
        return extension
    if is_igvf_accession(remote_or_local):
        if not igvf_file_type:
            raise ValueError(f"IGVF accession {remote_or_local} is missing a file type context.")
        return IGVF_FILE_TYPES[igvf_file_type][1]
    return ""


def igvf_credentials_from_args(args: argparse.Namespace, *, required: bool = False) -> Optional[IGVFCredentials]:
    return load_igvf_credentials(
        keypair_path=getattr(args, "igvf_keypair", None),
        api_key=getattr(args, "igvf_api_key", None),
        api_secret=getattr(args, "igvf_api_secret", None),
        required=required,
    )


def run_gsutil_copy(remote_path: str, local_path: Path, chunk_bytes: int, is_fastq: bool, env: Dict[str, str]) -> None:
    ensure_parent(local_path)
    if is_fastq and chunk_bytes > 0:
        with open(local_path, "wb") as handle:
            subprocess.run(
                ["gsutil", "-o", "GSUtil:parallel_process_count=1", "cat", "-r", f"0-{chunk_bytes}", remote_path],
                check=True,
                stdout=handle,
                env=env,
            )
    else:
        subprocess.run(
            ["gsutil", "-o", "GSUtil:parallel_process_count=1", "cp", remote_path, str(local_path)],
            check=True,
            env=env,
        )


def run_http_copy(
    remote_path: str,
    local_path: Path,
    chunk_bytes: int,
    is_fastq: bool,
    credentials: Optional[IGVFCredentials],
    expected_md5: Optional[str],
) -> None:
    range_end = chunk_bytes if is_fastq and chunk_bytes > 0 else None
    checksum = None if range_end is not None else expected_md5
    download_igvf_file(
        remote_path,
        local_path,
        credentials,
        byte_range_end=range_end,
        expected_md5=checksum,
    )


def run_igvf_accession_copy(
    accession: str,
    local_path: Path,
    chunk_bytes: int,
    is_fastq: bool,
    credentials: Optional[IGVFCredentials],
    expected_md5: Optional[str],
    igvf_file_type: str,
) -> None:
    run_http_copy(
        igvf_download_url(accession, igvf_file_type),
        local_path,
        chunk_bytes,
        is_fastq,
        credentials,
        expected_md5,
    )


def ensure_local_asset(
    remote_or_local: str,
    local_path: Path,
    chunk_bytes: int,
    is_fastq: bool,
    force_download: bool,
    env: Dict[str, str],
    credentials: Optional[IGVFCredentials],
    expected_md5: Optional[str] = None,
    igvf_file_type: Optional[str] = None,
) -> str:
    if local_path.exists() and local_path.stat().st_size > 0 and not force_download:
        return "reused"
    ensure_parent(local_path)
    if is_gcs_path(remote_or_local):
        run_gsutil_copy(remote_or_local, local_path, chunk_bytes, is_fastq, env)
        return "downloaded"
    if is_http_path(remote_or_local):
        run_http_copy(remote_or_local, local_path, chunk_bytes, is_fastq, credentials, expected_md5)
        return "downloaded"
    if is_igvf_accession(remote_or_local):
        if not igvf_file_type:
            raise ValueError(f"IGVF accession {remote_or_local} requires a file type context for download.")
        run_igvf_accession_copy(
            remote_or_local,
            local_path,
            chunk_bytes,
            is_fastq,
            credentials,
            expected_md5,
            igvf_file_type,
        )
        return "downloaded"
    source = Path(remote_or_local)
    if not source.exists():
        raise FileNotFoundError(f"Missing local asset: {remote_or_local}")
    if source.resolve() != local_path.resolve():
        shutil.copyfile(source, local_path)
    return "copied"


def ensure_group_assets(
    group_dir: Path,
    downloads_dir: Path,
    rows_by_modality: Dict[str, Dict[str, str]],
    chunk_bytes: int,
    force_download: bool,
    env: Dict[str, str],
    credentials: Optional[IGVFCredentials],
) -> Dict[str, Dict[str, str]]:
    resolved: Dict[str, Dict[str, str]] = {}
    for raw_modality, row in rows_by_modality.items():
        read1_ext = asset_extension(row["R1_path"], "sequence")
        read2_ext = asset_extension(row["R2_path"], "sequence")
        suffix = ""
        if chunk_bytes > 0:
            if chunk_bytes % (1024 * 1024) == 0:
                suffix = f"_{chunk_bytes // (1024 * 1024)}mb"
            else:
                suffix = f"_{chunk_bytes}b"
        local_r1 = downloads_dir / f"{raw_modality}_R1{suffix}{read1_ext}"
        local_r2 = downloads_dir / f"{raw_modality}_R2{suffix}{read2_ext}"
        ensure_local_asset(
            row["R1_path"],
            local_r1,
            chunk_bytes,
            True,
            force_download,
            env,
            credentials,
            expected_md5=row.get("R1_md5sum", ""),
            igvf_file_type="sequence",
        )
        ensure_local_asset(
            row["R2_path"],
            local_r2,
            chunk_bytes,
            True,
            force_download,
            env,
            credentials,
            expected_md5=row.get("R2_md5sum", ""),
            igvf_file_type="sequence",
        )

        modality_files = {
            "R1": str(local_r1),
            "R2": str(local_r2),
        }

        for column, tag in METADATA_COLUMNS.items():
            remote_path = row.get(column, "")
            if not remote_path:
                continue
            extension = asset_extension(remote_path, "tabular")
            local_path = downloads_dir / f"{raw_modality}_{tag}{extension}"
            ensure_local_asset(
                remote_path,
                local_path,
                chunk_bytes,
                False,
                force_download,
                env,
                credentials,
                expected_md5=row.get(f"{column}_md5sum", ""),
                igvf_file_type="tabular",
            )
            modality_files[column] = str(local_path)

        seqspec_remote = row.get("seqspec", "")
        if seqspec_remote:
            extension = asset_extension(seqspec_remote, "configuration") or ".yaml"
            local_seqspec = group_dir / "seqspec" / f"{raw_modality}_seqspec{extension}"
            ensure_local_asset(
                seqspec_remote,
                local_seqspec,
                chunk_bytes,
                False,
                force_download,
                env,
                credentials,
                expected_md5=row.get("seqspec_md5sum", ""),
                igvf_file_type="configuration",
            )
            modality_files["seqspec"] = str(local_seqspec)

        resolved[raw_modality] = modality_files
    return resolved


def load_seqspec_yaml(path: Path) -> Dict[str, object]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return yaml.load(handle, Loader=TagIgnoringLoader)


def derive_read_label(index: int, token: str) -> str:
    upper = token.upper()
    if "R1" in upper:
        return "R1"
    if "R2" in upper:
        return "R2"
    return f"R{index + 1}"


def extract_read_specs(seqspec_data: Dict[str, object]) -> List[ReadSpec]:
    read_specs: List[ReadSpec] = []
    for idx, read in enumerate(seqspec_data.get("sequence_spec", [])):
        files = read.get("files") or []
        file_info = files[0] if files else {}
        read_id = str(read.get("read_id") or read.get("name") or f"read_{idx + 1}")
        file_id = str(file_info.get("file_id") or file_info.get("filename") or read_id)
        filename = str(file_info.get("filename") or file_id)
        read_name = str(read.get("name") or filename or read_id)
        index_id = filename or file_id or read_id
        read_label = derive_read_label(idx, filename or read_id or read_name)
        read_specs.append(
            ReadSpec(
                file_id=file_id,
                filename=filename,
                read_id=read_id,
                index_id=index_id,
                read_label=read_label,
                strand=normalize_strand(read.get("strand", "forward")),
                min_len=int(read.get("min_len") or 0),
                max_len=int(read.get("max_len") or 0),
                name=read_name,
            )
        )
    return read_specs


def flatten_region_definitions(
    regions: Sequence[Dict[str, object]],
    definitions: Dict[str, Dict[str, object]],
) -> None:
    for region in regions or []:
        region_id = str(region.get("region_id") or "")
        if region_id:
            definitions[region_id] = region
        flatten_region_definitions(region.get("regions") or [], definitions)


def extract_region_definitions(seqspec_data: Dict[str, object]) -> Dict[str, Dict[str, object]]:
    definitions: Dict[str, Dict[str, object]] = {}
    flatten_region_definitions(seqspec_data.get("library_spec") or [], definitions)
    return definitions


def summarize_feature_metadata(path: Path) -> FeatureMetadataSummary:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    if path.suffix.lower() == ".gz":
        nested_suffixes = Path(path.stem).suffixes
        if nested_suffixes and nested_suffixes[-1].lower() == ".tsv":
            delimiter = "\t"
    with opener(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        rows = list(reader)
    if not rows:
        raise ValueError(f"Could not read headers from {path}")

    header = [field.strip() for field in rows[0]]
    sequence_index = next((idx for idx, column in enumerate(header) if column in {"spacer", "sequence", "guide", "seq"}), None)
    sequence_column = header[sequence_index] if sequence_index is not None else ""
    data_rows = rows[1:]

    if sequence_index is None:
        # Some hash metadata files are headerless two-column tables: id<TAB>sequence.
        for idx in (1, 0):
            if len(rows[0]) <= idx:
                continue
            value = str(rows[0][idx]).strip().upper()
            if value and set(value) <= DNA_ALPHABET:
                sequence_index = idx
                sequence_column = "sequence"
                data_rows = rows
                break
    if sequence_index is None:
        raise ValueError(f"No supported sequence column found in {path}")

    lengths: Counter[int] = Counter()
    total = 0
    for row in data_rows:
        if len(row) <= sequence_index:
            continue
        value = str(row[sequence_index]).strip().upper()
        if not value:
            continue
        lengths[len(value)] += 1
        total += 1
    return FeatureMetadataSummary(
        source_path=str(path),
        sequence_column=sequence_column,
        total_sequences=total,
        length_counts=lengths.most_common(),
    )


def canonical_region_name(region_id: str, region_type: str, seqspec_modality: str) -> Optional[str]:
    token_text = f"{region_id} {region_type}".lower()
    tokens = set(re.findall(r"[a-z0-9]+", token_text))
    region_type_token = str(region_type).strip().lower()
    if "barcode" in tokens:
        return "barcode"
    if "umi" in tokens:
        return "umi"
    if {"hash", "hashing", "hto", "tag"} & tokens:
        return "hash"
    if "guide" in tokens or "grna" in tokens or "spacer" in tokens or "protospacer" in tokens:
        return "guide"
    if "cdna" in tokens or token_text.strip() == "rna rna" or "mrna" in tokens:
        return "rna"
    ignored_types = {
        "custom_primer",
        "index5",
        "index7",
        "illumina_p5",
        "illumina_p7",
        "linker",
        "poly_t",
        "truseq_read1",
        "truseq_read2",
        "nextera_read1",
    }
    if seqspec_modality == "rna" and region_type_token not in {"barcode", "umi", *ignored_types}:
        return "rna"
    if seqspec_modality == "guide" and region_type_token not in {"barcode", "umi", *ignored_types}:
        return "guide"
    if seqspec_modality == "hash" and region_type_token not in {"barcode", "umi", *ignored_types}:
        return "hash"
    return None


def seqspec_region_priority(region: IndexedRegion) -> tuple[int, int, int]:
    token = f"{region.region_id} {region.region_type}".lower()
    score = 0
    if region.canonical_name in {"barcode", "umi"}:
        score += 100
    elif region.canonical_name == "rna":
        if region.region_type == "cdna" or region.region_id.lower() == "rna":
            score += 100
    elif region.canonical_name == "guide":
        if region.region_type in {"sgrna_target", "guide"} or region.region_id.lower() == "guide":
            score += 120
        if "scaffold" in token or region.region_type == "linker":
            score -= 120
    elif region.canonical_name == "hash":
        if any(name in token for name in ("hash", "hashing", "hto", "tag")):
            score += 100
    return score, -region.raw_start, -region.length


def run_seqspec_index(
    seqspec_bin: str,
    seqspec_path: Path,
    seqspec_modality: str,
    file_ids: Sequence[str],
    env: Dict[str, str],
) -> tuple[str, str]:
    joined_ids = ",".join(file_ids)
    base_cmd = [seqspec_bin, "index", "-m", seqspec_modality, "-s", "file", "-i", joined_ids, str(seqspec_path)]
    tab_result = subprocess.run(base_cmd, check=True, capture_output=True, text=True, env=env)
    kb_result = subprocess.run(base_cmd[:-1] + ["-t", "kb", str(seqspec_path)], check=True, capture_output=True, text=True, env=env)
    return tab_result.stdout.strip(), kb_result.stdout.strip()


def summarize_seqspec(raw_modality: str, seqspec_path: Path, seqspec_bin: str, env: Dict[str, str]) -> SeqspecSummary:
    data = load_seqspec_yaml(seqspec_path)
    seqspec_modalities = data.get("modalities") or []
    seqspec_mode_raw = str(seqspec_modalities[0] if seqspec_modalities else raw_modality)
    seqspec_modality = normalize_modality(seqspec_mode_raw)
    seqspec_index_modality = normalize_seqspec_index_modality(seqspec_mode_raw)
    read_specs = extract_read_specs(data)
    region_definitions = extract_region_definitions(data)
    file_ids = [read.index_id for read in read_specs]
    tab_output, kb_output = run_seqspec_index(seqspec_bin, seqspec_path, seqspec_index_modality, file_ids, env)
    read_map: Dict[str, ReadSpec] = {}
    for read in read_specs:
        for key in {read.file_id, read.filename, read.read_id, read.index_id, read.name}:
            if key:
                read_map[key] = read
    indexed_regions: List[IndexedRegion] = []
    for line in tab_output.splitlines():
        if not line.strip():
            continue
        file_id, region_id, region_type, start, stop = line.split("\t")
        read = read_map[file_id]
        region_definition = region_definitions.get(region_id, {})
        start_int = int(start)
        stop_int = int(stop)
        indexed_regions.append(
            IndexedRegion(
                file_id=file_id,
                read_label=read.read_label,
                strand=read.strand,
                region_id=region_id,
                region_type=region_type,
                canonical_name=canonical_region_name(region_id, region_type, seqspec_modality),
                raw_start=start_int,
                raw_end=stop_int - 1,
                length=stop_int - start_int,
                sequence_type=str(region_definition.get("sequence_type") or ""),
                template_sequence=str(region_definition.get("sequence") or ""),
                declared_min_len=int(region_definition.get("min_len") or 0),
                declared_max_len=int(region_definition.get("max_len") or 0),
            )
        )
    return SeqspecSummary(
        raw_modality=raw_modality,
        seqspec_modality=seqspec_modality,
        seqspec_path=str(seqspec_path),
        kb_string=kb_output,
        read_specs=read_specs,
        indexed_regions=indexed_regions,
    )


def seqspec_region_map(summary: SeqspecSummary) -> Dict[str, IndexedRegion]:
    region_map: Dict[str, IndexedRegion] = {}
    for row in summary.indexed_regions:
        if not row.canonical_name:
            continue
        existing = region_map.get(row.canonical_name)
        if existing is None or seqspec_region_priority(row) > seqspec_region_priority(existing):
            region_map[row.canonical_name] = row
    return region_map


def prediction_region_map(prediction_data: Dict[str, object]) -> Dict[str, Dict[str, object]]:
    region_map: Dict[str, Dict[str, object]] = {}
    if prediction_data.get("shared_barcode"):
        region_map["barcode"] = prediction_data["shared_barcode"]
    if prediction_data.get("shared_umi"):
        region_map["umi"] = prediction_data["shared_umi"]
    if prediction_data.get("rna_interval"):
        region_map["rna"] = prediction_data["rna_interval"]
    if prediction_data.get("guide_interval"):
        region_map["guide"] = prediction_data["guide_interval"]
    for item in prediction_data.get("hash_intervals", []):
        if item.get("region") == "hash":
            region_map["hash"] = item
    return region_map


def interval_text(read_label: str, strand: str, start: int, end: int) -> str:
    return f"{read_label} {strand} {start}-{end}"


def compare_region(
    modality: str,
    region: str,
    seqspec_region: Optional[IndexedRegion],
    prediction_region: Optional[Dict[str, object]],
) -> ComparisonRow:
    if seqspec_region is None:
        return ComparisonRow(modality, region, None, None, None, None, None, None, "na", None, "missing_seqspec", "Region not present in seqspec index output")
    seqspec_text = interval_text(seqspec_region.read_label, seqspec_region.strand, seqspec_region.raw_start, seqspec_region.raw_end)
    if prediction_region is None:
        return ComparisonRow(
            modality,
            region,
            seqspec_text,
            None,
            seqspec_region.read_label,
            None,
            seqspec_region.strand,
            None,
            "na",
            None,
            "missing_prediction",
            "Predictor did not emit this region",
        )

    prediction_text = interval_text(
        prediction_region["read_label"],
        prediction_region["strand"],
        int(prediction_region["raw_start"]),
        int(prediction_region["raw_end"]),
    )

    strand_flag = "match" if seqspec_region.strand == prediction_region["strand"] else "differs"

    if seqspec_region.read_label == prediction_region["read_label"]:
        start_distance = abs(seqspec_region.raw_start - int(prediction_region["raw_start"]))
        end_distance = abs(seqspec_region.raw_end - int(prediction_region["raw_end"]))
        max_distance = max(start_distance, end_distance)
        if max_distance == 0:
            flag = "perfect_match"
            note = "Exact raw interval match"
        elif max_distance <= 2:
            flag = "close_enough"
            note = "Within 2 bp of the seqspec interval"
        else:
            flag = "very_distant"
            note = "Same read, but interval differs by more than 2 bp"
        return ComparisonRow(
            modality,
            region,
            seqspec_text,
            prediction_text,
            seqspec_region.read_label,
            prediction_region["read_label"],
            seqspec_region.strand,
            prediction_region["strand"],
            strand_flag,
            max_distance,
            flag,
            note,
        )

    return ComparisonRow(
        modality,
        region,
        seqspec_text,
        prediction_text,
        seqspec_region.read_label,
        prediction_region["read_label"],
        seqspec_region.strand,
        prediction_region["strand"],
        strand_flag,
        None,
        "very_distant",
        "Read assignment differs",
    )


def build_comparison_rows(
    seqspec_summaries: Dict[str, SeqspecSummary],
    prediction_data: Dict[str, object],
) -> List[ComparisonRow]:
    predicted = prediction_region_map(prediction_data)
    rows: List[ComparisonRow] = []
    for raw_modality, summary in seqspec_summaries.items():
        feature_name = normalize_modality(raw_modality)
        region_order = ["barcode", "umi", feature_name]
        seqspec_regions = seqspec_region_map(summary)
        for region in region_order:
            rows.append(compare_region(raw_modality, region, seqspec_regions.get(region), predicted.get(region)))
    rows.sort(key=lambda item: (item.modality, FLAG_ORDER[item.flag], item.region))
    return rows


def write_seqspec_report(
    group_dir: Path,
    group_label: str,
    seqspec_summaries: Dict[str, SeqspecSummary],
    feature_metadata: Dict[str, FeatureMetadataSummary],
) -> Path:
    cards = []
    for raw_modality, summary in sorted(seqspec_summaries.items()):
        selected_regions = seqspec_region_map(summary)
        duplicated = [
            region
            for region in summary.indexed_regions
            if region.canonical_name and sum(item.canonical_name == region.canonical_name for item in summary.indexed_regions) > 1
        ]
        duplicate_note = ""
        if duplicated:
            rows = []
            seen = set()
            for region in duplicated:
                key = (region.canonical_name, region.region_id)
                if key in seen:
                    continue
                seen.add(key)
                rows.append(
                    f"<li><strong>{html.escape(region.canonical_name or '')}</strong>: "
                    f"{html.escape(region.region_id)} ({html.escape(region.region_type)}) "
                    f"{region.raw_start}-{region.raw_end}</li>"
                )
            duplicate_note = (
                "<div class='note'>"
                "<strong>Multiple seqspec regions map to the same comparison label.</strong> "
                "The comparison prefers the most specific target-like region, for example `guide` / `sgrna_target` over scaffold/linker regions."
                f"<ul>{''.join(rows)}</ul>"
                "</div>"
            )

        metadata_note = ""
        metadata_summary = feature_metadata.get(raw_modality)
        if metadata_summary is not None:
            lengths_text = ", ".join(f"{length} bp x {count}" for length, count in metadata_summary.length_counts)
            chosen_feature = selected_regions.get(normalize_modality(raw_modality))
            mismatch = ""
            if chosen_feature and metadata_summary.length_counts:
                valid_lengths = {length for length, _ in metadata_summary.length_counts}
                if chosen_feature.length not in valid_lengths:
                    mismatch = (
                        f" <strong>Length mismatch:</strong> seqspec selects {chosen_feature.length} bp "
                        f"but the metadata file contains {lengths_text}."
                    )
            metadata_note = (
                "<div class='note'>"
                f"<strong>Reference metadata:</strong> {html.escape(metadata_summary.sequence_column)} from "
                f"{html.escape(Path(metadata_summary.source_path).name)}. Lengths: {html.escape(lengths_text)}."
                f"{mismatch}"
                "</div>"
            )

        read_rows = "".join(
            "<tr>"
            f"<td>{html.escape(read.read_label)}</td>"
            f"<td>{html.escape(read.file_id)}</td>"
            f"<td>{html.escape(read.strand)}</td>"
            f"<td>{read.min_len}</td>"
            f"<td>{read.max_len}</td>"
            "</tr>"
            for read in summary.read_specs
        )
        region_rows = "".join(
            "<tr>"
            f"<td>{html.escape(region.read_label)}</td>"
            f"<td>{html.escape(region.region_id)}</td>"
            f"<td>{html.escape(region.region_type)}</td>"
            f"<td>{html.escape(region.canonical_name or '')}</td>"
            f"<td>{region.raw_start}</td>"
            f"<td>{region.raw_end}</td>"
            f"<td>{html.escape(region.sequence_type)}</td>"
            f"<td><code>{html.escape(region.template_sequence or '')}</code></td>"
            f"<td>{'yes' if region.canonical_name and selected_regions.get(region.canonical_name) == region else ''}</td>"
            "</tr>"
            for region in summary.indexed_regions
        )
        cards.append(
            "<section class='card'>"
            f"<h2>{html.escape(raw_modality)}</h2>"
            f"<p><strong>Seqspec file:</strong> {html.escape(summary.seqspec_path)}</p>"
            f"<p><strong>seqspec index -t kb:</strong> <code>{html.escape(summary.kb_string)}</code></p>"
            f"{duplicate_note}"
            f"{metadata_note}"
            "<h3>Read Specs</h3>"
            "<table><thead><tr><th>Read</th><th>File ID</th><th>Strand</th><th>Min Len</th><th>Max Len</th></tr></thead>"
            f"<tbody>{read_rows}</tbody></table>"
            "<h3>Indexed Regions</h3>"
            "<table><thead><tr><th>Read</th><th>Region ID</th><th>Region Type</th><th>Canonical</th><th>Start</th><th>End</th><th>Seq Type</th><th>Template</th><th>Used</th></tr></thead>"
            f"<tbody>{region_rows}</tbody></table>"
            "</section>"
        )

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Seqspec Report - {html.escape(group_label)}</title>
  <style>
    body {{ font-family: "Avenir Next", "Segoe UI", sans-serif; margin: 0; background: #f7f3eb; color: #2c221a; }}
    main {{ width: min(1180px, calc(100vw - 32px)); margin: 24px auto 40px; }}
    .hero, .card {{ background: #fffaf2; border: 1px solid #deceb7; border-radius: 22px; box-shadow: 0 14px 34px rgba(85, 64, 39, 0.08); }}
    .hero {{ padding: 24px 26px; }}
    .card {{ padding: 20px 22px; margin-top: 18px; }}
    .note {{ margin-top: 12px; padding: 12px 14px; background: #f5efe2; border: 1px solid #e5d8c1; border-radius: 14px; }}
    h1, h2, h3 {{ margin-top: 0; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
    th, td {{ text-align: left; padding: 9px 10px; border-bottom: 1px solid #eadfce; font-size: 14px; }}
    th {{ font-size: 12px; text-transform: uppercase; letter-spacing: 0.06em; color: #6f6355; }}
    code {{ font-family: "SFMono-Regular", "Menlo", "Consolas", monospace; background: #eef3ff; color: #264c98; padding: 2px 6px; border-radius: 6px; }}
  </style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>Seqspec Report</h1>
      <p>Analysis group: {html.escape(group_label)}. These intervals come directly from <code>seqspec index</code>.</p>
    </section>
    {''.join(cards)}
  </main>
</body>
</html>
"""
    out_path = group_dir / "seqspec_report.html"
    out_path.write_text(html_text, encoding="utf-8")
    return out_path


def write_comparison_report(
    group_dir: Path,
    group_label: str,
    comparison_rows: Sequence[ComparisonRow],
    seqspec_report_path: Path,
    prediction_report_path: Path,
) -> Path:
    counts = Counter(row.flag for row in comparison_rows)
    table_rows = "".join(
        "<tr>"
        f"<td>{html.escape(row.modality)}</td>"
        f"<td>{html.escape(row.region)}</td>"
        f"<td>{html.escape(row.seqspec_interval or 'NA')}</td>"
        f"<td>{html.escape(row.prediction_interval or 'NA')}</td>"
        f"<td><span class=\"strand-pill {html.escape(row.strand_flag)}\">{html.escape(row.strand_flag)}</span></td>"
        f"<td><span class=\"pill {html.escape(row.flag)}\">{html.escape(row.flag)}</span></td>"
        f"<td>{row.max_distance_bp if row.max_distance_bp is not None else 'NA'}</td>"
        f"<td>{html.escape(row.note)}</td>"
        "</tr>"
        for row in comparison_rows
    )
    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Comparison Report - {html.escape(group_label)}</title>
  <style>
    body {{ font-family: "Avenir Next", "Segoe UI", sans-serif; margin: 0; background: #f5f6f1; color: #1f2720; }}
    main {{ width: min(1180px, calc(100vw - 32px)); margin: 24px auto 40px; }}
    .hero, .card {{ background: #fcfdf8; border: 1px solid #d5dcc9; border-radius: 22px; box-shadow: 0 14px 34px rgba(56, 73, 47, 0.08); }}
    .hero {{ padding: 24px 26px; }}
    .card {{ padding: 20px 22px; margin-top: 18px; }}
    .summary {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; margin-top: 16px; }}
    .pill {{ display: inline-flex; padding: 5px 10px; border-radius: 999px; font-size: 12px; font-weight: 700; text-transform: uppercase; }}
    .strand-pill {{ display: inline-flex; padding: 5px 10px; border-radius: 999px; font-size: 12px; font-weight: 700; text-transform: uppercase; }}
    .perfect_match {{ background: #dff5e7; color: #146a34; }}
    .close_enough {{ background: #eef4d6; color: #556b12; }}
    .very_distant {{ background: #fbe4dd; color: #9a2f1f; }}
    .missing_prediction, .missing_seqspec {{ background: #eceaf3; color: #574c79; }}
    .match {{ background: #e5f2f1; color: #24645c; }}
    .differs {{ background: #f6eadb; color: #8d5b18; }}
    .na {{ background: #eceaf3; color: #574c79; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
    th, td {{ text-align: left; padding: 9px 10px; border-bottom: 1px solid #e4ead9; font-size: 14px; vertical-align: top; }}
    th {{ font-size: 12px; text-transform: uppercase; letter-spacing: 0.06em; color: #61705b; }}
    a {{ color: #2b5ab3; text-decoration: none; }}
  </style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>Seqspec vs Predictor Comparison</h1>
      <p>Analysis group: {html.escape(group_label)}. A region is <strong>perfect_match</strong> when both raw intervals are identical on the same read, <strong>close_enough</strong> when the max boundary difference is at most 2 bp on the same read, and <strong>very_distant</strong> otherwise. Strand differences are shown separately and do not count as wrong by themselves.</p>
      <p><a href="{html.escape(seqspec_report_path.name)}">Seqspec report</a> | <a href="{html.escape(prediction_report_path.name)}">Prediction report</a></p>
      <div class="summary">
        <div class="card"><strong>Perfect</strong><br />{counts['perfect_match']}</div>
        <div class="card"><strong>Close</strong><br />{counts['close_enough']}</div>
        <div class="card"><strong>Distant</strong><br />{counts['very_distant']}</div>
        <div class="card"><strong>Missing</strong><br />{counts['missing_prediction'] + counts['missing_seqspec']}</div>
      </div>
    </section>
    <section class="card">
      <table>
        <thead><tr><th>Modality</th><th>Region</th><th>Seqspec</th><th>Prediction</th><th>Strand</th><th>Flag</th><th>Max Distance</th><th>Note</th></tr></thead>
        <tbody>{table_rows}</tbody>
      </table>
    </section>
  </main>
</body>
</html>
"""
    out_path = group_dir / "comparison.html"
    out_path.write_text(html_text, encoding="utf-8")
    return out_path


def write_analysis_index(analysis_root: Path, group_summaries: Sequence[Dict[str, str]]) -> Path:
    rows = "".join(
        "<tr>"
        f"<td>{html.escape(item['label'])}</td>"
        f"<td><a href=\"{html.escape(item['relative_dir'])}/comparison.html\">comparison</a></td>"
        f"<td><a href=\"{html.escape(item['relative_dir'])}/seqspec_report.html\">seqspec</a></td>"
        f"<td><a href=\"{html.escape(item['relative_dir'])}/prediction_report.html\">prediction</a></td>"
        "</tr>"
        for item in group_summaries
    )
    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Seqspec Analysis Index</title>
  <style>
    body {{ font-family: "Avenir Next", "Segoe UI", sans-serif; margin: 0; background: #f4f1ec; color: #2b241d; }}
    main {{ width: min(980px, calc(100vw - 32px)); margin: 24px auto 40px; }}
    .card {{ background: #fffaf2; border: 1px solid #dccdb9; border-radius: 22px; box-shadow: 0 14px 34px rgba(76, 57, 35, 0.08); padding: 22px; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
    th, td {{ text-align: left; padding: 10px; border-bottom: 1px solid #eadfcc; }}
    th {{ font-size: 12px; text-transform: uppercase; letter-spacing: 0.06em; color: #6f6456; }}
    a {{ color: #2d57b0; text-decoration: none; }}
  </style>
</head>
<body>
  <main>
    <section class="card">
      <h1>Seqspec Analysis Index</h1>
      <table>
        <thead><tr><th>Group</th><th>Comparison</th><th>Seqspec</th><th>Prediction</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
  </main>
</body>
</html>
"""
    out_path = analysis_root / "index.html"
    out_path.write_text(html_text, encoding="utf-8")
    return out_path


def write_download_index(analysis_root: Path, group_summaries: Sequence[Dict[str, str]]) -> Path:
    rows = "".join(
        "<tr>"
        f"<td>{html.escape(item['label'])}</td>"
        f"<td>{html.escape(item['relative_dir'])}/downloads</td>"
        f"<td><a href=\"{html.escape(item['relative_dir'])}/download_summary.json\">summary</a></td>"
        "</tr>"
        for item in group_summaries
    )
    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Seqspec Download Index</title>
  <style>
    body {{ font-family: "Avenir Next", "Segoe UI", sans-serif; margin: 0; background: #f4f1ec; color: #2b241d; }}
    main {{ width: min(980px, calc(100vw - 32px)); margin: 24px auto 40px; }}
    .card {{ background: #fffaf2; border: 1px solid #dccdb9; border-radius: 22px; box-shadow: 0 14px 34px rgba(76, 57, 35, 0.08); padding: 22px; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
    th, td {{ text-align: left; padding: 10px; border-bottom: 1px solid #eadfcc; }}
    th {{ font-size: 12px; text-transform: uppercase; letter-spacing: 0.06em; color: #6f6456; }}
    a {{ color: #2d57b0; text-decoration: none; }}
  </style>
</head>
<body>
  <main>
    <section class="card">
      <h1>Seqspec Download Index</h1>
      <table>
        <thead><tr><th>Group</th><th>Downloads</th><th>Summary</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
  </main>
</body>
</html>
"""
    out_path = analysis_root / "download_index.html"
    out_path.write_text(html_text, encoding="utf-8")
    return out_path


def run_prediction_report(
    downloads_dir: Path,
    group_dir: Path,
    group_label: str,
    args: argparse.Namespace,
    env: Dict[str, str],
) -> Dict[str, object]:
    predictor_script = Path(__file__).with_name("seqspec_check.py")
    html_path = group_dir / "prediction_report.html"
    json_path = group_dir / "prediction_report.json"
    cmd = [
        sys.executable,
        str(predictor_script),
        "--input-dir",
        str(downloads_dir),
        "--barcode-sample-reads",
        str(args.barcode_sample_reads),
        "--feature-sample-reads",
        str(args.feature_sample_reads),
        "--guide-upstream-bases",
        str(args.guide_upstream_bases),
        "--title",
        f"Prediction Report - {group_label}",
        "--output-html",
        str(html_path),
        "--output-json",
        str(json_path),
    ]
    subprocess.run(cmd, check=True, env=env)
    return json.loads(json_path.read_text(encoding="utf-8"))


def process_group(
    group: Dict[str, object],
    analysis_root: Path,
    args: argparse.Namespace,
    env: Dict[str, str],
) -> Dict[str, str]:
    credentials = igvf_credentials_from_args(args, required=False)
    group_label = group["label"] or "analysis_group"
    group_dir = analysis_root / sanitize_label(group_label)
    if args.downloads_dir:
        downloads_dir = Path(args.downloads_dir)
    else:
        downloads_dir = group_dir / "downloads"
    downloads_dir.mkdir(parents=True, exist_ok=True)
    (group_dir / "seqspec").mkdir(parents=True, exist_ok=True)

    assets = ensure_group_assets(
        group_dir=group_dir,
        downloads_dir=downloads_dir,
        rows_by_modality=group["rows_by_modality"],
        chunk_bytes=args.chunk_bytes,
        force_download=args.force_download,
        env=env,
        credentials=credentials,
    )

    if args.download_only:
        summary_path = group_dir / "download_summary.json"
        summary_payload = {
            "group_label": group_label,
            "downloads_dir": str(downloads_dir),
            "staged_assets": assets,
            "chunk_bytes": args.chunk_bytes,
        }
        summary_path.write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")
        return {
            "label": group_label,
            "relative_dir": group_dir.name,
            "download_summary": str(summary_path),
            "downloads_dir": str(downloads_dir),
        }

    seqspec_summaries: Dict[str, SeqspecSummary] = {}
    feature_metadata: Dict[str, FeatureMetadataSummary] = {}
    for raw_modality, paths in assets.items():
        seqspec_path = paths.get("seqspec")
        if not seqspec_path:
            continue
        seqspec_summaries[raw_modality] = summarize_seqspec(
            raw_modality=raw_modality,
            seqspec_path=Path(seqspec_path),
            seqspec_bin=args.seqspec_bin,
            env=env,
        )
        feature_name = normalize_modality(raw_modality)
        metadata_path = paths.get("barcode_hashtag_map") if feature_name == "hash" else paths.get("guide_design")
        if metadata_path:
            feature_metadata[raw_modality] = summarize_feature_metadata(Path(metadata_path))

    prediction_data = run_prediction_report(downloads_dir, group_dir, group_label, args, env)
    comparison_rows = build_comparison_rows(seqspec_summaries, prediction_data)
    seqspec_report = write_seqspec_report(group_dir, group_label, seqspec_summaries, feature_metadata)
    comparison_report = write_comparison_report(
        group_dir,
        group_label,
        comparison_rows,
        seqspec_report,
        group_dir / "prediction_report.html",
    )

    summary_payload = {
        "group_label": group_label,
        "downloads_dir": str(downloads_dir),
        "seqspec_summaries": {key: asdict(value) for key, value in seqspec_summaries.items()},
        "feature_metadata": {key: asdict(value) for key, value in feature_metadata.items()},
        "prediction_summary": prediction_data,
        "comparison_rows": [asdict(row) for row in comparison_rows],
    }
    (group_dir / "analysis_summary.json").write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")

    return {
        "label": group_label,
        "relative_dir": group_dir.name,
        "comparison_report": str(comparison_report),
        "seqspec_report": str(seqspec_report),
        "prediction_report": str(group_dir / "prediction_report.html"),
    }


def build_runtime_env(analysis_root: Path) -> Dict[str, str]:
    analysis_root.mkdir(parents=True, exist_ok=True)
    mplconfig = analysis_root / ".mplconfig"
    mplconfig.mkdir(parents=True, exist_ok=True)
    gsutil_state = analysis_root / ".gsutil"
    gsutil_state.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["MPLCONFIGDIR"] = str(mplconfig)
    env["GSUTIL_STATE_DIR"] = str(gsutil_state)
    return env


def upstream_generate_per_sample_script() -> Optional[Path]:
    candidate = Path(__file__).resolve().parent / "vendor" / "CRISPR_Pipeline" / "download_development" / "generate_per_sample.py"
    if candidate.exists():
        return candidate
    return None


def write_igvf_samplesheet_via_upstream(
    args: argparse.Namespace,
    env: Dict[str, str],
    samplesheet_path: Path,
) -> bool:
    script_path = upstream_generate_per_sample_script()
    if script_path is None:
        return False

    cmd = [
        sys.executable,
        str(script_path),
        "--accession",
        args.accession,
        "--output",
        str(samplesheet_path),
    ]
    if args.igvf_keypair:
        cmd.extend(["--keypair", args.igvf_keypair])
    if args.hash_seqspec:
        cmd.extend(["--hash_seqspec", args.hash_seqspec])
    if args.rna_seqspec:
        cmd.extend(["--rna_seqspec", args.rna_seqspec])
    if args.sgrna_seqspec:
        cmd.extend(["--sgrna_seqspec", args.sgrna_seqspec])

    upstream_env = dict(env)
    if args.igvf_api_key:
        upstream_env["IGVF_API_KEY"] = args.igvf_api_key
    if args.igvf_api_secret:
        upstream_env["IGVF_SECRET_KEY"] = args.igvf_api_secret
    try:
        subprocess.run(cmd, check=True, env=upstream_env)
    except subprocess.CalledProcessError as exc:
        print(
            f"Upstream IGVF samplesheet generation failed for {args.accession} "
            f"(exit {exc.returncode}); falling back to the local generator."
        )
        return False
    return True


def sortable_group_value(value: str) -> tuple[int, object]:
    token = str(value or "").strip()
    if token.isdigit():
        return (0, int(token))
    return (1, token)


def select_first_complete_lane_filters(args: argparse.Namespace, samplesheet_path: Path) -> List[str]:
    rows = canonicalize_samplesheet_rows(load_rows(str(samplesheet_path)))
    filtered_rows = apply_filters(rows, parse_csv_list(args.modalities), parse_filters(args.filter))
    if not filtered_rows:
        raise ValueError("No IGVF rows remained after applying filters before one-lane selection.")

    by_lane: Dict[tuple[str, str], Dict[str, object]] = {}
    for row in filtered_rows:
        key = (row.get("sequencing_run", ""), row.get("lane", ""))
        entry = by_lane.setdefault(key, {"rows": [], "modalities": set()})
        entry["rows"].append(row)
        entry["modalities"].add(row.get("file_modality", ""))

    complete_keys = [
        key
        for key, entry in by_lane.items()
        if {"scRNA", "gRNA"}.issubset(entry["modalities"])
    ]
    if not complete_keys:
        raise ValueError("Could not find a complete scRNA+gRNA lane in the generated IGVF samplesheet.")

    selected_run, selected_lane = min(
        complete_keys,
        key=lambda item: (sortable_group_value(item[0]), sortable_group_value(item[1])),
    )
    return [f"sequencing_run={selected_run}", f"lane={selected_lane}"]


def run_samplesheet_path(
    args: argparse.Namespace,
    analysis_root: Path,
    env: Dict[str, str],
    samplesheet_path: Path,
) -> None:
    if not args.download_only:
        ensure_seqspec_cli(args.seqspec_bin, env)

    group_columns = parse_csv_list(args.group_by)
    sort_columns = parse_csv_list(args.sort_by)
    rows = canonicalize_samplesheet_rows(load_rows(str(samplesheet_path)))
    filtered_rows = apply_filters(rows, parse_csv_list(args.modalities), parse_filters(args.filter))
    if not filtered_rows:
        raise ValueError("No rows remained after applying filters.")

    analysis_groups = group_rows_by_analysis(stable_sort_rows(filtered_rows, sort_columns), group_columns)
    if args.downloads_dir and len(analysis_groups) > 1:
        raise ValueError("--downloads-dir can only be used when filters reduce the run to a single analysis group.")

    group_summaries: List[Dict[str, str]] = []
    for group in analysis_groups:
        row_modalities = set(group["rows_by_modality"])
        if not {"scRNA", "gRNA"}.issubset(row_modalities):
            continue
        group_summaries.append(process_group(group, analysis_root, args, env))

    if not group_summaries:
        raise ValueError("No complete analysis groups with at least scRNA and gRNA were found.")

    if args.download_only:
        index_path = write_download_index(analysis_root, group_summaries)
        print(f"Wrote download index to {index_path}")
        for item in group_summaries:
            print(f"  - {item['label']}: {item['downloads_dir']}")
        return

    index_path = write_analysis_index(analysis_root, group_summaries)
    print(f"Wrote analysis index to {index_path}")
    for item in group_summaries:
        print(f"  - {item['label']}: {item['relative_dir']}")


def run_samplesheet_mode(args: argparse.Namespace, analysis_root: Path, env: Dict[str, str]) -> None:
    run_samplesheet_path(args, analysis_root, env, Path(args.samplesheet))


def run_igvf_mode(args: argparse.Namespace, analysis_root: Path, env: Dict[str, str]) -> None:
    if args.samplesheet_out:
        samplesheet_path = Path(args.samplesheet_out)
    else:
        samplesheet_path = analysis_root / f"{sanitize_label(args.accession)}_samplesheet.tsv"
    if not write_igvf_samplesheet_via_upstream(args, env, samplesheet_path):
        credentials = igvf_credentials_from_args(args, required=True)
        rows = generate_samplesheet_rows(
            args.accession,
            credentials,
            hash_seqspec=args.hash_seqspec,
            rna_seqspec=args.rna_seqspec,
            sgrna_seqspec=args.sgrna_seqspec,
        )
        write_igvf_samplesheet(rows, samplesheet_path)
    print(f"Wrote IGVF samplesheet to {samplesheet_path}")
    if args.samplesheet_only:
        return
    effective_args = argparse.Namespace(**vars(args))
    if args.one_lane:
        lane_filters = select_first_complete_lane_filters(args, samplesheet_path)
        effective_args.filter = [*list(args.filter), *lane_filters]
        print(f"Selected one lane for processing: {lane_filters[0]}, {lane_filters[1]}")
    run_samplesheet_path(effective_args, analysis_root, env, samplesheet_path)


def run_single_lane_mode(args: argparse.Namespace, analysis_root: Path, env: Dict[str, str]) -> None:
    if not args.download_only:
        ensure_seqspec_cli(args.seqspec_bin, env)
    group = build_single_lane_group(args)
    group_summary = process_group(group, analysis_root, args, env)
    if args.download_only:
        index_path = write_download_index(analysis_root, [group_summary])
        print(f"Wrote download index to {index_path}")
        print(f"  - {group_summary['label']}: {group_summary['downloads_dir']}")
        return
    index_path = write_analysis_index(analysis_root, [group_summary])
    print(f"Wrote analysis index to {index_path}")
    print(f"  - {group_summary['label']}: {group_summary['relative_dir']}")


def main() -> None:
    args = parse_args()
    analysis_root = Path(args.analysis_root)
    analysis_root.mkdir(parents=True, exist_ok=True)
    env = build_runtime_env(analysis_root)

    if args.mode == "samplesheet":
        run_samplesheet_mode(args, analysis_root, env)
        return
    if args.mode == "igvf":
        run_igvf_mode(args, analysis_root, env)
        return
    if args.mode == "single-lane":
        run_single_lane_mode(args, analysis_root, env)
        return
    raise ValueError(f"Unsupported mode: {args.mode}")


if __name__ == "__main__":
    main()
