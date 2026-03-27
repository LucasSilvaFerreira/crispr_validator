#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
import os
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime
from glob import glob
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


DNA_ALPHABET = set("ACGTN")
BARCODE_SAMPLE_READS_DEFAULT = 10_000
FEATURE_SAMPLE_READS_DEFAULT = 100_000
UMI_LENGTH_CANDIDATES = (10, 12)
UMI_SHORT_LENGTH = 10
UMI_LONG_LENGTH = 12
UMI_MIN_LONG_SCORE_DELTA = -0.02
UMI_MIN_LONG_UNIQUE_GAIN = 0.015
UMI_MIN_LONG_COLLISION_REDUCTION = 0.10
UMI_MIN_LONG_TAIL_ENTROPY = 1.65
READ_LABELS = ("R1", "R2")
STRANDS = ("forward", "reverse")


@dataclass
class IntervalCall:
    region: str
    source: str
    read_label: str
    strand: str
    oriented_start: int
    oriented_end: int
    raw_start: int
    raw_end: int
    length: int
    confidence: str
    note: str


@dataclass
class BarcodeConfigResult:
    modality: str
    read_label: str
    strand: str
    read_length: int
    total_reads: int
    total_hits: int
    reads_with_hit: int
    dominant_position: Optional[int]
    dominant_count: int
    hit_ratio: float
    position_purity: float
    dominant_read_fraction: float
    score: float
    raw_start: Optional[int]
    raw_end: Optional[int]
    top_positions: List[Tuple[int, int]]
    position_counts: Dict[int, int]


@dataclass
class SequenceConfigResult:
    modality: str
    feature_name: str
    read_label: str
    strand: str
    read_length: int
    feature_length: int
    total_reads: int
    total_hits: int
    dominant_position: Optional[int]
    dominant_count: int
    hit_ratio: float
    position_purity: float
    flank_purity: float
    gini: float
    score: float
    raw_start: Optional[int]
    raw_end: Optional[int]
    top_positions: List[Tuple[int, int]]
    top_feature_hits: List[Tuple[str, int]]
    top_flanks: List[Tuple[str, int]]
    top_interval_sequences: List[Tuple[str, int]]
    position_counts: Dict[int, int]


@dataclass
class UmiCandidate:
    side: str
    gap: int
    start: int
    length: int
    raw_start: int
    raw_end: int
    unique_fraction: float
    unique_count: int
    top_sequence_fraction: float
    entropy_bits: float
    per_base_entropy: List[float]
    tail_entropy_bits: float
    n_fraction: float
    score: float
    top_sequences: List[Tuple[str, int]]
    selected: bool = False


def is_gzipped(filename: str) -> bool:
    with open(filename, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_file(filename: str, mode: str = "rt"):
    if is_gzipped(filename):
        return gzip.open(filename, mode)
    return open(filename, mode)


def get_reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(table)[::-1]


def clean_dna(token: str) -> Optional[str]:
    seq = token.strip().upper()
    if not seq:
        return None
    if set(seq) - DNA_ALPHABET:
        return None
    return seq


def read_fastq_sequences(filename: str, max_reads: int) -> List[str]:
    sequences: List[str] = []
    with open_file(filename, "rt") as handle:
        while len(sequences) < max_reads:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip().upper()
            handle.readline()
            handle.readline()
            if seq:
                sequences.append(seq)
    return sequences


def mode_read_length(reads: Sequence[str]) -> int:
    if not reads:
        return 0
    return Counter(len(read) for read in reads).most_common(1)[0][0]


def load_barcode_whitelist(path: str) -> Tuple[set[str], int]:
    barcodes: List[str] = []
    with open_file(path, "rt") as handle:
        for raw_line in handle:
            token = raw_line.strip().split("\t")[0].split(",")[0]
            seq = clean_dna(token)
            if seq:
                barcodes.append(seq)
    if not barcodes:
        raise ValueError(f"No barcode sequences found in {path}")
    length = Counter(len(seq) for seq in barcodes).most_common(1)[0][0]
    whitelist = {seq for seq in barcodes if len(seq) == length}
    return whitelist, length


def load_feature_sequences(path: str) -> List[str]:
    delimiter = "\t" if path.endswith(".tsv") or path.endswith(".tsv.gz") else ","
    with open_file(path, "rt") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        rows = list(reader)
    if not rows:
        raise ValueError(f"Could not read headers from {path}")

    header = [field.strip() for field in rows[0]]
    target_index = next((idx for idx, col in enumerate(header) if col in {"spacer", "sequence", "guide", "seq"}), None)
    data_rows = rows[1:]

    if target_index is None:
        # Some hash metadata files are headerless two-column tables: id<TAB>sequence.
        candidate_indices = [idx for idx in (1, 0) if len(rows[0]) > idx and clean_dna(rows[0][idx])]
        if not candidate_indices:
            raise ValueError(f"No sequence column found in {path}; expected one of spacer/sequence/guide/seq")
        target_index = candidate_indices[0]
        data_rows = rows

    sequences: List[str] = []
    for row in data_rows:
        token = row[target_index] if len(row) > target_index else ""
        seq = clean_dna(str(token))
        if seq:
            sequences.append(seq)
    if not sequences:
        raise ValueError(f"No valid feature sequences found in {path}")
    return list(dict.fromkeys(sequences))


def calculate_gini(values: Sequence[int]) -> float:
    if not values:
        return 0.0
    array = sorted(values)
    total = sum(array)
    if total == 0:
        return 0.0
    n_items = len(array)
    weighted = sum((2 * idx - n_items - 1) * value for idx, value in enumerate(array, start=1))
    return weighted / (n_items * total)


def shannon_entropy(chars: Sequence[str]) -> float:
    if not chars:
        return 0.0
    counts = Counter(chars)
    total = len(chars)
    entropy = 0.0
    for count in counts.values():
        frequency = count / total
        entropy -= frequency * math.log2(frequency)
    return entropy


def collect_window_values(reads: Sequence[str], start: int, length: int) -> List[str]:
    return [read[start : start + length] for read in reads if 0 <= start and len(read) >= start + length]


def oriented_to_raw(oriented_start: int, length: int, read_length: int, strand: str) -> Tuple[int, int]:
    if strand == "forward":
        raw_start = oriented_start
    else:
        raw_start = read_length - oriented_start - length
    raw_end = raw_start + length - 1
    return raw_start, raw_end


def raw_to_oriented(raw_start: int, length: int, read_length: int, strand: str) -> int:
    if strand == "forward":
        return raw_start
    return read_length - raw_start - length


def get_oriented_reads(reads: Sequence[str], strand: str) -> List[str]:
    if strand == "forward":
        return list(reads)
    return [get_reverse_complement(read) for read in reads]


def scan_barcode_config(
    modality: str,
    read_label: str,
    strand: str,
    reads: Sequence[str],
    whitelist: set[str],
    reverse_whitelist: set[str],
    barcode_length: int,
) -> BarcodeConfigResult:
    read_length = mode_read_length(reads)
    active_whitelist = whitelist if strand == "forward" else reverse_whitelist
    position_counts: Counter[int] = Counter()
    total_hits = 0
    reads_with_hit = 0

    for read in reads:
        found_any = False
        max_start = len(read) - barcode_length + 1
        if max_start < 1:
            continue
        for idx in range(max_start):
            if read[idx : idx + barcode_length] in active_whitelist:
                position_counts[idx] += 1
                total_hits += 1
                found_any = True
        if found_any:
            reads_with_hit += 1

    total_reads = len(reads)
    dominant_position: Optional[int] = None
    dominant_count = 0
    raw_start = None
    raw_end = None

    if position_counts:
        dominant_position, dominant_count = position_counts.most_common(1)[0]
        raw_start = dominant_position
        raw_end = dominant_position + barcode_length - 1

    hit_ratio = reads_with_hit / total_reads if total_reads else 0.0
    position_purity = dominant_count / total_hits if total_hits else 0.0
    dominant_read_fraction = dominant_count / total_reads if total_reads else 0.0
    score = (4.0 * dominant_read_fraction) + (2.0 * position_purity) + hit_ratio

    return BarcodeConfigResult(
        modality=modality,
        read_label=read_label,
        strand=strand,
        read_length=read_length,
        total_reads=total_reads,
        total_hits=total_hits,
        reads_with_hit=reads_with_hit,
        dominant_position=dominant_position,
        dominant_count=dominant_count,
        hit_ratio=hit_ratio,
        position_purity=position_purity,
        dominant_read_fraction=dominant_read_fraction,
        score=score,
        raw_start=raw_start,
        raw_end=raw_end,
        top_positions=position_counts.most_common(10),
        position_counts=dict(sorted(position_counts.items())),
    )


def pick_barcode_source_modalities(
    available_modalities: Sequence[str],
    requested_source: str,
) -> List[str]:
    if requested_source == "auto":
        ordered = [name for name in ("rna", "guide", "hash") if name in available_modalities]
        return ordered
    if requested_source in available_modalities:
        return [requested_source]
    return []


def confidence_from_barcode(result: BarcodeConfigResult) -> str:
    if result.dominant_read_fraction >= 0.60 and result.position_purity >= 0.90:
        return "high"
    if result.dominant_read_fraction >= 0.25 and result.position_purity >= 0.70:
        return "medium"
    return "low"


def analyze_feature_config(
    modality: str,
    feature_name: str,
    read_label: str,
    strand: str,
    reads: Sequence[str],
    features: Sequence[str],
) -> SequenceConfigResult:
    read_length = mode_read_length(reads)
    flank_window = 12
    feature_hits: Counter[str] = Counter({feature: 0 for feature in features})
    position_counts: Counter[int] = Counter()
    upstream_map: Dict[int, List[str]] = defaultdict(list)

    for read in reads:
        for feature in features:
            idx = read.find(feature)
            if idx != -1:
                position_counts[idx] += 1
                feature_hits[feature] += 1
                fragment = read[max(0, idx - flank_window) : idx]
                if len(fragment) < flank_window:
                    fragment = ("-" * (flank_window - len(fragment))) + fragment
                upstream_map[idx].append(fragment)
                break

    total_reads = len(reads)
    total_hits = sum(position_counts.values())
    dominant_position: Optional[int] = None
    dominant_count = 0
    raw_start = None
    raw_end = None
    feature_length = Counter(len(feature) for feature in features).most_common(1)[0][0]

    if position_counts:
        dominant_position, dominant_count = position_counts.most_common(1)[0]
        raw_start = dominant_position
        raw_end = dominant_position + feature_length - 1

    hit_ratio = total_hits / total_reads if total_reads else 0.0
    position_purity = dominant_count / total_hits if total_hits else 0.0
    if dominant_position is None or not upstream_map[dominant_position]:
        flank_purity = 0.0
        top_flanks: List[Tuple[str, int]] = []
    else:
        flank_purity = Counter(upstream_map[dominant_position]).most_common(1)[0][1] / len(upstream_map[dominant_position])
        top_flanks = Counter(upstream_map[dominant_position]).most_common(10)
    top_interval_sequences: List[Tuple[str, int]] = []
    if dominant_position is not None:
        oriented_reads = get_oriented_reads(reads, strand)
        oriented_start = raw_to_oriented(dominant_position, feature_length, read_length, strand)
        top_interval_sequences = Counter(
            collect_window_values(oriented_reads, oriented_start, feature_length)
        ).most_common(10)
    gini = calculate_gini(list(feature_hits.values()))
    score = (3.0 * hit_ratio) + (2.0 * position_purity) + flank_purity + (1.0 - gini)

    return SequenceConfigResult(
        modality=modality,
        feature_name=feature_name,
        read_label=read_label,
        strand=strand,
        read_length=read_length,
        feature_length=feature_length,
        total_reads=total_reads,
        total_hits=total_hits,
        dominant_position=dominant_position,
        dominant_count=dominant_count,
        hit_ratio=hit_ratio,
        position_purity=position_purity,
        flank_purity=flank_purity,
        gini=gini,
        score=score,
        raw_start=raw_start,
        raw_end=raw_end,
        top_positions=position_counts.most_common(10),
        top_feature_hits=feature_hits.most_common(10),
        top_flanks=top_flanks,
        top_interval_sequences=top_interval_sequences,
        position_counts=dict(sorted(position_counts.items())),
    )


def confidence_from_sequence_result(result: SequenceConfigResult) -> str:
    if result.hit_ratio >= 0.60 and result.position_purity >= 0.90:
        return "high"
    if result.hit_ratio >= 0.15 and result.position_purity >= 0.70:
        return "medium"
    return "low"


def score_umi_window(reads: Sequence[str], start: int, length: int) -> Optional[Tuple[int, float, float, float, float, List[float]]]:
    values = collect_window_values(reads, start, length)
    if not values:
        return None
    unique_count = len(set(values))
    unique_fraction = unique_count / len(values)
    top_sequence_fraction = Counter(values).most_common(1)[0][1] / len(values)
    per_base_entropy = [shannon_entropy([value[idx] for value in values]) for idx in range(length)]
    entropy_bits = sum(per_base_entropy) / length
    n_fraction = sum("N" in value for value in values) / len(values)
    return unique_count, unique_fraction, top_sequence_fraction, entropy_bits, n_fraction, per_base_entropy


def umi_anchor(candidate: UmiCandidate) -> Tuple[str, int, int]:
    return candidate.side, candidate.gap, candidate.start


def umi_collision_reduction(short_candidate: UmiCandidate, long_candidate: UmiCandidate) -> float:
    short_collision = max(1.0 - short_candidate.unique_fraction, 1e-9)
    long_collision = max(1.0 - long_candidate.unique_fraction, 0.0)
    return 1.0 - (long_collision / short_collision)


def choose_umi_winner(candidates: Sequence[UmiCandidate]) -> Tuple[UmiCandidate, str]:
    ranked = sorted(
        candidates,
        key=lambda item: (item.score, item.unique_fraction, item.entropy_bits, -item.length),
        reverse=True,
    )
    provisional = ranked[0]
    comparable = [item for item in ranked if umi_anchor(item) == umi_anchor(provisional)]
    by_length = {item.length: item for item in comparable}
    shorter = min(comparable, key=lambda item: item.length)
    short_candidate = by_length.get(UMI_SHORT_LENGTH)
    long_candidate = by_length.get(UMI_LONG_LENGTH)
    if short_candidate is not None and long_candidate is not None:
        score_gain = long_candidate.score - short_candidate.score
        unique_gain = long_candidate.unique_fraction - short_candidate.unique_fraction
        collision_reduction = umi_collision_reduction(short_candidate, long_candidate)
        if provisional == long_candidate and (
            score_gain < UMI_MIN_LONG_SCORE_DELTA
            or (
                unique_gain < UMI_MIN_LONG_UNIQUE_GAIN
                and collision_reduction < UMI_MIN_LONG_COLLISION_REDUCTION
            )
            or long_candidate.tail_entropy_bits < UMI_MIN_LONG_TAIL_ENTROPY
        ):
            return (
                short_candidate,
                "Preferred 10 bp because the added 2 bp did not add strong extra UMI evidence "
                f"(score delta {score_gain:.3f}, unique gain {unique_gain:.3f}, collision reduction {collision_reduction:.3f}, "
                f"tail entropy {long_candidate.tail_entropy_bits:.3f}).",
            )
        if provisional == short_candidate and (
            score_gain >= UMI_MIN_LONG_SCORE_DELTA
            and (
                unique_gain >= UMI_MIN_LONG_UNIQUE_GAIN
                or collision_reduction >= UMI_MIN_LONG_COLLISION_REDUCTION
            )
            and long_candidate.tail_entropy_bits >= UMI_MIN_LONG_TAIL_ENTROPY
        ):
            return (
                long_candidate,
                "Extended to 12 bp because the added 2 bp remained highly variable and materially improved "
                f"collision resolution over the 10 bp anchor (score delta {score_gain:.3f}, unique gain {unique_gain:.3f}, "
                f"collision reduction {collision_reduction:.3f}).",
            )
    if (
        shorter.length < provisional.length
        and shorter.score >= provisional.score - 0.05
        and shorter.unique_fraction >= provisional.unique_fraction - 0.02
        and shorter.entropy_bits >= provisional.entropy_bits - 0.03
    ):
        return shorter, "Preferred the shortest adjacent window because the longer candidate was only marginally better."
    return provisional, ""


def predict_umi(
    source_modality: str,
    source_reads: Dict[str, List[str]],
    barcode_call: BarcodeConfigResult,
    whitelist: set[str],
    barcode_length: int,
) -> Tuple[Optional[IntervalCall], List[UmiCandidate], int]:
    if barcode_call.dominant_position is None:
        return None, [], 0

    supporting_reads = []
    oriented_start = raw_to_oriented(
        barcode_call.dominant_position,
        barcode_length,
        barcode_call.read_length,
        barcode_call.strand,
    )
    oriented_reads = get_oriented_reads(source_reads[barcode_call.read_label], barcode_call.strand)
    for read in oriented_reads:
        start = oriented_start
        end = start + barcode_length
        if len(read) >= end and read[start:end] in whitelist:
            supporting_reads.append(read)

    candidates: List[UmiCandidate] = []
    start = oriented_start + barcode_length
    for length in UMI_LENGTH_CANDIDATES:
        values = collect_window_values(supporting_reads, start, length)
        if not values:
            continue
        metrics = score_umi_window(supporting_reads, start, length)
        if metrics is None:
            continue
        unique_count, unique_fraction, top_sequence_fraction, entropy_bits, n_fraction, per_base_entropy = metrics
        score = (
            (1.5 * unique_fraction)
            + (1.1 * (entropy_bits / 2.0))
            - (0.8 * top_sequence_fraction)
            - (0.5 * n_fraction)
        )
        raw_start, raw_end = oriented_to_raw(start, length, barcode_call.read_length, barcode_call.strand)
        tail_entropy_bits = 0.0
        if length > UMI_SHORT_LENGTH:
            tail_entropy_bits = sum(per_base_entropy[UMI_SHORT_LENGTH:length]) / (length - UMI_SHORT_LENGTH)
        candidates.append(
            UmiCandidate(
                side="downstream",
                gap=0,
                start=start,
                length=length,
                raw_start=raw_start,
                raw_end=raw_end,
                unique_fraction=unique_fraction,
                unique_count=unique_count,
                top_sequence_fraction=top_sequence_fraction,
                entropy_bits=entropy_bits,
                per_base_entropy=per_base_entropy,
                tail_entropy_bits=tail_entropy_bits,
                n_fraction=n_fraction,
                score=score,
                top_sequences=Counter(values).most_common(10),
            )
        )

    if not candidates:
        return None, [], len(supporting_reads)

    winner, selection_note = choose_umi_winner(candidates)
    candidates.sort(key=lambda item: (0 if item == winner else 1, -item.score))
    for candidate in candidates:
        candidate.selected = candidate == winner
    if winner.score >= 1.7 and winner.unique_fraction >= 0.70:
        confidence = "high"
    elif winner.score >= 1.25 and winner.unique_fraction >= 0.45:
        confidence = "medium"
    else:
        confidence = "low"

    umi_call = IntervalCall(
        region="umi",
        source=f"{source_modality} variability heuristic",
        read_label=barcode_call.read_label,
        strand=barcode_call.strand,
        oriented_start=winner.start,
        oriented_end=winner.start + winner.length - 1,
        raw_start=winner.raw_start,
        raw_end=winner.raw_end,
        length=winner.length,
        confidence=confidence,
        note=(
            f"Predicted from {len(supporting_reads)} barcode-supporting reads; "
            "compared 10 bp and 12 bp windows immediately after the barcode"
            + (f". {selection_note}" if selection_note else "")
        ),
    )
    return umi_call, candidates[:12], len(supporting_reads)


def build_barcode_interval(result: BarcodeConfigResult, barcode_length: int, source: str) -> Optional[IntervalCall]:
    if result.dominant_position is None or result.raw_start is None or result.raw_end is None:
        return None
    oriented_start = raw_to_oriented(result.dominant_position, barcode_length, result.read_length, result.strand)
    return IntervalCall(
        region="barcode",
        source=source,
        read_label=result.read_label,
        strand=result.strand,
        oriented_start=oriented_start,
        oriented_end=oriented_start + barcode_length - 1,
        raw_start=result.raw_start,
        raw_end=result.raw_end,
        length=barcode_length,
        confidence=confidence_from_barcode(result),
        note=f"Whitelist-dominant interval supported by {result.dominant_count} / {result.total_reads} sampled reads",
    )


def build_sequence_interval(result: SequenceConfigResult) -> Optional[IntervalCall]:
    if result.dominant_position is None or result.raw_start is None or result.raw_end is None:
        return None
    oriented_start = raw_to_oriented(result.dominant_position, result.feature_length, result.read_length, result.strand)
    return IntervalCall(
        region=result.feature_name,
        source=f"{result.modality} exact-match scan",
        read_label=result.read_label,
        strand=result.strand,
        oriented_start=oriented_start,
        oriented_end=oriented_start + result.feature_length - 1,
        raw_start=result.raw_start,
        raw_end=result.raw_end,
        length=result.feature_length,
        confidence=confidence_from_sequence_result(result),
        note=f"Dominant interval supported by {result.dominant_count} / {result.total_reads} sampled reads",
    )


def build_guide_prefix_interval(result: SequenceConfigResult, flank_length: int) -> Optional[IntervalCall]:
    if flank_length <= 0 or result.raw_start is None:
        return None
    raw_end = result.raw_start - 1
    if raw_end < 0:
        return None
    raw_start = max(0, raw_end - flank_length + 1)
    interval_length = raw_end - raw_start + 1
    top_flank = result.top_flanks[0][0] if result.top_flanks else ""
    note = f"Top {interval_length} bp prefix directly before the winning guide interval"
    if top_flank:
        note += f"; most common observed prefix: {top_flank}"
    oriented_start = raw_to_oriented(raw_start, interval_length, result.read_length, result.strand)
    return IntervalCall(
        region="guide_prefix",
        source="guide dominant-position flank",
        read_label=result.read_label,
        strand=result.strand,
        oriented_start=oriented_start,
        oriented_end=oriented_start + interval_length - 1,
        raw_start=raw_start,
        raw_end=raw_end,
        length=interval_length,
        confidence=confidence_from_sequence_result(result),
        note=note,
    )


def interval_to_kallisto_triplet(interval: IntervalCall, use_full_read: bool = False) -> str:
    file_index = 0 if interval.read_label == "R1" else 1
    start = 0 if use_full_read else interval.raw_start
    stop = 0 if use_full_read else interval.raw_end + 1
    return f"{file_index},{start},{stop}"


def build_kallisto_like_string(
    barcode: IntervalCall,
    umi: Optional[IntervalCall],
    feature: IntervalCall,
    use_full_feature_read: bool = False,
) -> Optional[Dict[str, str]]:
    if umi is None:
        return None
    technology = (
        f"{interval_to_kallisto_triplet(barcode)}:"
        f"{interval_to_kallisto_triplet(umi)}:"
        f"{interval_to_kallisto_triplet(feature, use_full_read=use_full_feature_read)}"
    )
    note = f"Feature read {feature.read_label} is interpreted on the {feature.strand} strand."
    if use_full_feature_read:
        note = f"Third triplet uses the full {feature.read_label} read; stop=0 follows kallisto convention."
    return {
        "label": feature.region,
        "technology": technology,
        "note": note,
    }


def find_first_existing(directory: str, patterns: Sequence[str]) -> Optional[str]:
    for pattern in patterns:
        matches = sorted(glob(os.path.join(directory, pattern)))
        if matches:
            return matches[0]
    return None


def resolve_inputs(args: argparse.Namespace) -> Dict[str, object]:
    input_dir = args.input_dir
    resolved: Dict[str, object] = {
        "barcode_whitelist": args.barcode_whitelist,
        "guide_metadata": args.guide_metadata,
        "hash_metadata": args.hash_metadata,
        "modalities": {},
    }

    if input_dir:
        if not resolved["barcode_whitelist"]:
            resolved["barcode_whitelist"] = find_first_existing(
                input_dir,
                [
                    "scRNA_barcodeWhitelist.tsv*",
                    "gRNA_barcodeWhitelist.tsv*",
                    "hash_barcodeWhitelist.tsv*",
                    "*barcodeWhitelist.tsv*",
                ],
            )
        if not resolved["guide_metadata"]:
            resolved["guide_metadata"] = find_first_existing(
                input_dir,
                [
                    "gRNA_guideMetadata.tsv*",
                    "scRNA_guideMetadata.tsv*",
                    "hash_guideMetadata.tsv*",
                    "*guideMetadata.tsv*",
                ],
            )
        if not resolved["hash_metadata"]:
            resolved["hash_metadata"] = find_first_existing(
                input_dir,
                [
                    "hash_hashMetadata.tsv*",
                    "*hashMetadata.tsv*",
                ],
            )

    modality_specs = {
        "rna": {
            "display_name": "RNA",
            "r1": args.rna_r1,
            "r2": args.rna_r2,
            "fallback_r1": ["scRNA_R1*.fastq*", "*rna*R1*.fastq*"],
            "fallback_r2": ["scRNA_R2*.fastq*", "*rna*R2*.fastq*"],
        },
        "guide": {
            "display_name": "Guide",
            "r1": args.guide_r1,
            "r2": args.guide_r2,
            "fallback_r1": ["gRNA_R1*.fastq*", "*guide*R1*.fastq*", "*gRNA*R1*.fastq*"],
            "fallback_r2": ["gRNA_R2*.fastq*", "*guide*R2*.fastq*", "*gRNA*R2*.fastq*"],
        },
        "hash": {
            "display_name": "Hash",
            "r1": args.hash_r1,
            "r2": args.hash_r2,
            "fallback_r1": ["hash_R1*.fastq*", "*hash*R1*.fastq*"],
            "fallback_r2": ["hash_R2*.fastq*", "*hash*R2*.fastq*"],
        },
    }

    for modality, spec in modality_specs.items():
        read1 = spec["r1"]
        read2 = spec["r2"]
        if input_dir:
            if not read1:
                read1 = find_first_existing(input_dir, spec["fallback_r1"])
            if not read2:
                read2 = find_first_existing(input_dir, spec["fallback_r2"])
        if bool(read1) != bool(read2):
            raise ValueError(f"{modality} requires both read1 and read2 paths")
        if read1 and read2:
            resolved["modalities"][modality] = {
                "display_name": spec["display_name"],
                "R1": read1,
                "R2": read2,
            }

    if not resolved["modalities"]:
        raise ValueError("No modalities were provided. Supply explicit paths or use --input-dir.")
    if not resolved["barcode_whitelist"]:
        raise ValueError("A barcode whitelist is required. Supply --barcode-whitelist or use --input-dir.")
    if "guide" in resolved["modalities"] and not resolved["guide_metadata"]:
        raise ValueError("Guide reads were provided but guide metadata was not found.")
    if "hash" in resolved["modalities"] and not resolved["hash_metadata"]:
        raise ValueError("Hash reads were provided but hash metadata was not found.")

    return resolved


def load_sampled_reads(
    modality_paths: Dict[str, Dict[str, str]],
    max_reads: int,
) -> Dict[str, Dict[str, object]]:
    loaded: Dict[str, Dict[str, object]] = {}
    for modality, paths in modality_paths.items():
        read1 = read_fastq_sequences(paths["R1"], max_reads)
        read2 = read_fastq_sequences(paths["R2"], max_reads)
        loaded[modality] = {
            "display_name": paths["display_name"],
            "R1_path": paths["R1"],
            "R2_path": paths["R2"],
            "R1": read1,
            "R2": read2,
            "R1_length": mode_read_length(read1),
            "R2_length": mode_read_length(read2),
        }
    return loaded


def choose_shared_barcode_call(
    sampled_reads: Dict[str, Dict[str, object]],
    whitelist: set[str],
    reverse_whitelist: set[str],
    barcode_length: int,
    requested_source: str,
) -> Tuple[str, BarcodeConfigResult, Dict[str, List[BarcodeConfigResult]], str]:
    available_modalities = list(sampled_reads)
    source_candidates = pick_barcode_source_modalities(available_modalities, requested_source)
    if not source_candidates:
        available = ", ".join(available_modalities)
        raise ValueError(f"Requested barcode source '{requested_source}' is unavailable. Available: {available}")

    all_results: Dict[str, List[BarcodeConfigResult]] = {}
    chosen_modality = source_candidates[0]
    chosen_best: Optional[BarcodeConfigResult] = None

    for modality in source_candidates:
        results: List[BarcodeConfigResult] = []
        for read_label in READ_LABELS:
            reads = sampled_reads[modality][read_label]
            for strand in STRANDS:
                results.append(
                    scan_barcode_config(
                        modality=modality,
                        read_label=read_label,
                        strand=strand,
                        reads=reads,
                        whitelist=whitelist,
                        reverse_whitelist=reverse_whitelist,
                        barcode_length=barcode_length,
                    )
                )
        results.sort(key=lambda item: item.score, reverse=True)
        all_results[modality] = results
        if chosen_best is None:
            chosen_modality = modality
            chosen_best = results[0]
        if results[0].dominant_read_fraction >= 0.20 and results[0].position_purity >= 0.70:
            chosen_modality = modality
            chosen_best = results[0]
            break

    assert chosen_best is not None
    explanation = (
        f"Shared barcode/UMI inferred from {chosen_modality} because it is the preferred available modality "
        f"and its best whitelist interval is concentrated at one region."
    )
    return chosen_modality, chosen_best, all_results, explanation


def call_modality_features(
    modality: str,
    reads: Dict[str, object],
    feature_name: str,
    feature_sequences: Sequence[str],
) -> Tuple[SequenceConfigResult, List[SequenceConfigResult]]:
    results: List[SequenceConfigResult] = []
    forward = list(feature_sequences)
    reverse = [get_reverse_complement(seq) for seq in feature_sequences]

    for read_label in READ_LABELS:
        read_values = reads[read_label]
        results.append(
            analyze_feature_config(modality, feature_name, read_label, "forward", read_values, forward)
        )
        results.append(
            analyze_feature_config(modality, feature_name, read_label, "reverse", read_values, reverse)
        )

    results.sort(key=lambda item: item.score, reverse=True)
    return results[0], results


def make_rna_interval(
    shared_barcode: IntervalCall,
    rna_lengths: Dict[str, int],
) -> IntervalCall:
    other_read = "R2" if shared_barcode.read_label == "R1" else "R1"
    other_length = rna_lengths[other_read]
    opposite_strand = "reverse" if shared_barcode.strand == "forward" else "forward"
    return IntervalCall(
        region="rna",
        source="RNA assumption from shared barcode orientation",
        read_label=other_read,
        strand=opposite_strand,
        oriented_start=0,
        oriented_end=max(other_length - 1, 0),
        raw_start=0,
        raw_end=max(other_length - 1, 0),
        length=other_length,
        confidence="assumed",
        note="Placed on the full opposite read/strand relative to the inferred barcode read",
    )


def clamp_interval(interval: IntervalCall, read_length: int) -> IntervalCall:
    raw_start = max(0, min(interval.raw_start, max(read_length - 1, 0)))
    raw_end = max(raw_start, min(interval.raw_end, max(read_length - 1, 0)))
    length = raw_end - raw_start + 1 if read_length else 0
    return IntervalCall(
        region=interval.region,
        source=interval.source,
        read_label=interval.read_label,
        strand=interval.strand,
        oriented_start=interval.oriented_start,
        oriented_end=interval.oriented_end,
        raw_start=raw_start,
        raw_end=raw_end,
        length=length,
        confidence=interval.confidence,
        note=interval.note,
    )


def render_svg_histogram(position_counts: Dict[int, int], read_length: int, color: str) -> str:
    width = 560
    height = 170
    padding_left = 38
    padding_right = 14
    padding_top = 14
    padding_bottom = 28
    plot_width = width - padding_left - padding_right
    plot_height = height - padding_top - padding_bottom

    if not position_counts:
        return (
            f"<svg viewBox='0 0 {width} {height}' class='histogram'>"
            f"<rect x='0' y='0' width='{width}' height='{height}' rx='16' fill='#fbf7ef' />"
            f"<text x='{width / 2}' y='{height / 2}' text-anchor='middle' fill='#7b6f63' font-size='15'>No dominant positions detected</text>"
            f"</svg>"
        )

    max_count = max(position_counts.values())
    max_index = max(read_length - 1, max(position_counts))
    bars = []
    for position, count in position_counts.items():
        x = padding_left + (position / max(1, max_index)) * plot_width
        bar_width = max(plot_width / max(1, max_index + 1), 2.0)
        bar_height = (count / max_count) * plot_height
        y = padding_top + (plot_height - bar_height)
        bars.append(
            f"<rect x='{x:.2f}' y='{y:.2f}' width='{bar_width:.2f}' height='{bar_height:.2f}' "
            f"rx='1.4' fill='{color}' opacity='0.85' />"
        )

    return (
        f"<svg viewBox='0 0 {width} {height}' class='histogram'>"
        f"<rect x='0' y='0' width='{width}' height='{height}' rx='16' fill='#fbf7ef' />"
        f"<line x1='{padding_left}' y1='{padding_top + plot_height}' x2='{width - padding_right}' y2='{padding_top + plot_height}' stroke='#bcae9a' stroke-width='1' />"
        f"<line x1='{padding_left}' y1='{padding_top}' x2='{padding_left}' y2='{padding_top + plot_height}' stroke='#bcae9a' stroke-width='1' />"
        f"{''.join(bars)}"
        f"<text x='{padding_left}' y='{height - 9}' font-size='11' fill='#7b6f63'>0</text>"
        f"<text x='{width - padding_right}' y='{height - 9}' text-anchor='end' font-size='11' fill='#7b6f63'>{max_index}</text>"
        f"<text x='12' y='{padding_top + 8}' font-size='11' fill='#7b6f63'>{max_count}</text>"
        f"</svg>"
    )


def render_read_tracks(
    read_lengths: Dict[str, int],
    intervals: Sequence[Tuple[IntervalCall, str]],
) -> str:
    width = 860
    row_height = 70
    header_height = 14
    total_height = 2 * row_height + header_height
    left_pad = 98
    right_pad = 28
    track_width = width - left_pad - right_pad

    svg_parts = [
        f"<svg viewBox='0 0 {width} {total_height}' class='read-map'>",
        f"<rect x='0' y='0' width='{width}' height='{total_height}' rx='20' fill='#fffaf2' />",
    ]

    for row_index, read_label in enumerate(READ_LABELS):
        y_base = header_height + row_index * row_height
        read_length = max(read_lengths.get(read_label, 0), 1)
        svg_parts.append(
            f"<text x='26' y='{y_base + 28}' fill='#40362b' font-size='16' font-weight='700'>{read_label}</text>"
        )
        svg_parts.append(
            f"<rect x='{left_pad}' y='{y_base + 14}' width='{track_width}' height='18' rx='9' fill='#ebe0cf' />"
        )
        svg_parts.append(
            f"<text x='{left_pad}' y='{y_base + 52}' fill='#7b6f63' font-size='11'>0</text>"
        )
        svg_parts.append(
            f"<text x='{left_pad + track_width}' y='{y_base + 52}' text-anchor='end' fill='#7b6f63' font-size='11'>{read_length - 1}</text>"
        )

        row_intervals = [(interval, color) for interval, color in intervals if interval.read_label == read_label]
        for item_index, (interval, color) in enumerate(row_intervals):
            start = max(0, min(interval.raw_start, read_length - 1))
            end = max(start, min(interval.raw_end, read_length - 1))
            width_ratio = (end - start + 1) / max(read_length, 1)
            x = left_pad + (start / max(read_length, 1)) * track_width
            rect_width = max(track_width * width_ratio, 4.0)
            y = y_base + 14 + (item_index % 2) * 7
            svg_parts.append(
                f"<rect x='{x:.2f}' y='{y:.2f}' width='{rect_width:.2f}' height='18' rx='6' fill='{color}' opacity='0.88' />"
            )
            label_x = min(x + rect_width + 8, width - 12)
            svg_parts.append(
                f"<text x='{label_x:.2f}' y='{y + 13:.2f}' fill='{color}' font-size='12' font-weight='700'>{html.escape(interval.region)}</text>"
            )

    svg_parts.append("</svg>")
    return "".join(svg_parts)


def format_interval(interval: IntervalCall) -> str:
    return f"{interval.raw_start}-{interval.raw_end}"


def render_interval_table(intervals: Sequence[IntervalCall]) -> str:
    rows = []
    for interval in intervals:
        rows.append(
            "<tr>"
            f"<td>{html.escape(interval.region)}</td>"
            f"<td>{html.escape(interval.read_label)}</td>"
            f"<td>{html.escape(interval.strand)}</td>"
            f"<td>{format_interval(interval)}</td>"
            f"<td>{interval.length}</td>"
            f"<td><span class='pill pill-{html.escape(interval.confidence)}'>{html.escape(interval.confidence)}</span></td>"
            f"<td>{html.escape(interval.note)}</td>"
            "</tr>"
        )
    return (
        "<table class='metrics-table'>"
        "<thead><tr><th>Region</th><th>Read</th><th>Strand</th><th>Raw Interval</th><th>Length</th><th>Confidence</th><th>Evidence</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def render_barcode_config_table(results: Sequence[BarcodeConfigResult]) -> str:
    rows = []
    for result in results:
        interval_text = "NA"
        if result.raw_start is not None and result.raw_end is not None:
            interval_text = f"{result.raw_start}-{result.raw_end}"
        rows.append(
            "<tr>"
            f"<td>{result.read_label} / {result.strand}</td>"
            f"<td>{result.score:.3f}</td>"
            f"<td>{result.hit_ratio:.3f}</td>"
            f"<td>{result.position_purity:.3f}</td>"
            f"<td>{result.dominant_read_fraction:.3f}</td>"
            f"<td>{interval_text}</td>"
            f"<td>{result.dominant_count}</td>"
            "</tr>"
        )
    return (
        "<table class='metrics-table'>"
        "<thead><tr><th>Config</th><th>Score</th><th>Hit Ratio</th><th>Pos Purity</th><th>Dominant / Reads</th><th>Raw Interval</th><th>Dominant Reads</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def render_sequence_config_table(results: Sequence[SequenceConfigResult]) -> str:
    rows = []
    for result in results:
        interval_text = "NA"
        if result.raw_start is not None and result.raw_end is not None:
            interval_text = f"{result.raw_start}-{result.raw_end}"
        rows.append(
            "<tr>"
            f"<td>{result.read_label} / {result.strand}</td>"
            f"<td>{result.score:.3f}</td>"
            f"<td>{result.hit_ratio:.3f}</td>"
            f"<td>{result.position_purity:.3f}</td>"
            f"<td>{result.flank_purity:.3f}</td>"
            f"<td>{1.0 - result.gini:.3f}</td>"
            f"<td>{interval_text}</td>"
            "</tr>"
        )
    return (
        "<table class='metrics-table'>"
        "<thead><tr><th>Config</th><th>Score</th><th>Hit Ratio</th><th>Pos Purity</th><th>Flank Purity</th><th>Diversity</th><th>Raw Interval</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def render_umi_candidates(candidates: Sequence[UmiCandidate]) -> str:
    if not candidates:
        return "<p class='empty-block'>No UMI candidate window cleared the adjacency search.</p>"
    rows = []
    for candidate in candidates:
        winner_badge = "<span class='pill pill-high'>selected</span>" if candidate.selected else ""
        tail_entropy = f"{candidate.tail_entropy_bits:.3f}" if candidate.length > UMI_SHORT_LENGTH else "-"
        rows.append(
            "<tr>"
            f"<td>{winner_badge}</td>"
            f"<td>{candidate.raw_start}-{candidate.raw_end}</td>"
            f"<td>{candidate.length}</td>"
            f"<td>{candidate.unique_fraction:.3f}</td>"
            f"<td>{candidate.entropy_bits:.3f}</td>"
            f"<td>{tail_entropy}</td>"
            f"<td>{candidate.top_sequence_fraction:.4f}</td>"
            f"<td>{candidate.score:.3f}</td>"
            "</tr>"
        )
    return (
        "<table class='metrics-table'>"
        "<thead><tr><th>Winner</th><th>Raw Interval</th><th>Length</th><th>Unique Fraction</th><th>Entropy</th><th>12 bp Tail Entropy</th><th>Top Seq Fraction</th><th>Score</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def render_sequence_spectrum(
    title: str,
    subtitle: str,
    items: Sequence[Tuple[str, int]],
    color: str,
) -> str:
    if not items:
        return ""
    max_count = max(count for _, count in items)
    rows = []
    for sequence, count in items:
        width_pct = (count / max_count) * 100 if max_count else 0
        rows.append(
            "<div class='sequence-row'>"
            f"<code>{html.escape(sequence)}</code>"
            "<div class='sequence-bar-wrap'>"
            f"<div class='sequence-bar' style='width: {width_pct:.2f}%; background: {color};'></div>"
            f"<span>{count}</span>"
            "</div>"
            "</div>"
        )
    return (
        "<div class='diagnostics-block spectrum-block'>"
        f"<div class='section-header'><h3>{html.escape(title)}</h3>"
        f"<p>{html.escape(subtitle)}</p></div>"
        f"{''.join(rows)}"
        "</div>"
    )


def render_umi_spectra(candidates: Sequence[UmiCandidate]) -> str:
    panels = []
    for candidate in candidates[:4]:
        label = f"{candidate.raw_start}-{candidate.raw_end} ({candidate.length} bp)"
        subtitle = "Selected UMI window." if candidate.selected else "Alternative adjacent UMI candidate."
        panels.append(render_sequence_spectrum(label, subtitle, candidate.top_sequences[:8], "#14866d"))
    return "".join(panels)


def render_feature_string_panels(
    feature_name: str,
    results: Sequence[SequenceConfigResult],
    color: str,
) -> str:
    panels = []
    for result in results:
        if result.dominant_position is None or result.total_hits == 0:
            continue
        interval_text = f"{result.raw_start}-{result.raw_end}"
        title = f"{result.read_label} / {result.strand} at {interval_text}"
        subtitle = (
            f"{feature_name.capitalize()} exact matches: {result.total_hits} total hits, "
            f"{result.dominant_count} at the dominant interval."
        )
        exact_hits = [(sequence, count) for sequence, count in result.top_feature_hits if count > 0][:8]
        if exact_hits:
            panels.append(
                render_sequence_spectrum(
                    f"{title} matched {feature_name} strings",
                    subtitle,
                    exact_hits,
                    color,
                )
            )
        observed = [(sequence, count) for sequence, count in result.top_interval_sequences if count > 0][:8]
        if observed:
            panels.append(
                render_sequence_spectrum(
                    f"{title} observed interval strings",
                    "Most common strings observed directly at the dominant interval, oriented to the called strand.",
                    observed,
                    color,
                )
            )
    return "".join(panels)


def render_explanation_section(guide_upstream_bases: int) -> str:
    items = [
        "Barcode: only the first 10k reads from the selected source modality are scanned against the whitelist across both reads and both strands; the winning call is the dominant on-list region, not just any match.",
        "UMI: only the 10 bp and 12 bp windows immediately after the barcode are scored. A 12 bp call is only accepted when the added 2 bp also carry strong tail entropy and materially reduce collisions relative to the anchored 10 bp window.",
        "RNA: the cDNA region is placed on the full opposite read/strand relative to the shared barcode call.",
        "Guide and hash: exact matches are scored across read and strand configurations using hit ratio, position purity, flank purity, and diversity, following the original guide-first script logic.",
        f"Guide prefix: when requested, the report also returns the {guide_upstream_bases} bp immediately before the winning guide interval in raw FASTQ coordinates.",
        "Kallisto-like strings: technology strings are emitted as bc:umi:feature triplets using 0-based start and exclusive stop coordinates, matching the kallisto-style syntax.",
    ]
    rows = "".join(f"<li>{html.escape(item)}</li>" for item in items)
    return (
        "<section class='card'>"
        "<div class='section-header'><h2>How The Regions Are Computed</h2>"
        "<p>Short summary of the inference rules used in this report.</p></div>"
        f"<ul class='explanation-list'>{rows}</ul>"
        "</section>"
    )


def render_kallisto_table(entries: Sequence[Dict[str, str]]) -> str:
    if not entries:
        return (
            "<section class='card'>"
            "<div class='section-header'><h2>Kallisto-Like Technology Strings</h2>"
            "<p>No technology string was emitted because the shared UMI could not be inferred.</p></div>"
            "</section>"
        )
    rows = []
    for entry in entries:
        rows.append(
            "<tr>"
            f"<td>{html.escape(entry['label'])}</td>"
            f"<td><code>{html.escape(entry['technology'])}</code></td>"
            f"<td>{html.escape(entry['note'])}</td>"
            "</tr>"
        )
    return (
        "<section class='card'>"
        "<div class='section-header'><h2>Kallisto-Like Technology Strings</h2>"
        "<p>These use the `bc:umi:feature` triplet syntax analogous to kallisto bus technology strings.</p></div>"
        "<table class='metrics-table'>"
        "<thead><tr><th>Region</th><th>Technology String</th><th>Note</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
        "</section>"
    )


def render_flank_table(title: str, flanks: Sequence[Tuple[str, int]]) -> str:
    if not flanks:
        return ""
    rows = []
    for sequence, count in flanks:
        rows.append(
            "<tr>"
            f"<td><code>{html.escape(sequence)}</code></td>"
            f"<td>{count}</td>"
            "</tr>"
        )
    return (
        "<div class='diagnostics-block'>"
        f"<div class='section-header'><h3>{html.escape(title)}</h3>"
        "<p>Most common raw-read prefixes directly before the winning guide interval.</p></div>"
        "<table class='metrics-table'>"
        "<thead><tr><th>Prefix Sequence</th><th>Count</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
        "</div>"
    )


def render_modality_section(
    title: str,
    subtitle: str,
    read_lengths: Dict[str, int],
    intervals: Sequence[IntervalCall],
    colors: Dict[str, str],
    diagnostics_html: str,
) -> str:
    track_items = [(interval, colors.get(interval.region, "#6b7280")) for interval in intervals]
    return (
        "<section class='card modality-card'>"
        f"<div class='section-header'><h2>{html.escape(title)}</h2><p>{html.escape(subtitle)}</p></div>"
        f"{render_read_tracks(read_lengths, track_items)}"
        f"{render_interval_table(intervals)}"
        f"{diagnostics_html}"
        "</section>"
    )


def build_html_report(report: Dict[str, object]) -> str:
    colors = {
        "barcode": "#c26a11",
        "umi": "#14866d",
        "rna": "#2262c6",
        "guide": "#a8224f",
        "hash": "#6d38bb",
    }

    summary_cards = [
        ("Barcode Source", report["shared_barcode_source"]),
        ("Shared Barcode", report["shared_barcode_interval"]),
        ("Predicted UMI", report["shared_umi_interval"]),
        ("RNA Assumption", report["rna_interval"]),
    ]

    summary_html = "".join(
        "<div class='summary-card'>"
        f"<div class='summary-label'>{html.escape(label)}</div>"
        f"<div class='summary-value'>{html.escape(value)}</div>"
        "</div>"
        for label, value in summary_cards
    )

    barcode_diag = (
        "<section class='card diagnostics-card'>"
        "<div class='section-header'><h2>Shared Barcode Diagnostics</h2>"
        f"<p>{html.escape(report['barcode_explanation'])}</p></div>"
        f"<p class='meta-line'>Barcode scan used the first {report['barcode_sample_reads']} reads from {report['shared_barcode_source']}.</p>"
        f"{render_barcode_config_table(report['barcode_config_results'])}"
        f"{render_svg_histogram(report['barcode_config_results'][0].position_counts, report['barcode_config_results'][0].read_length, colors['barcode'])}"
        "<div class='section-header umi-header'><h3>UMI Candidates</h3>"
        f"<p>{html.escape(report['umi_note'])}</p></div>"
        f"{render_umi_candidates(report['umi_candidates'])}"
        f"{render_umi_spectra(report['umi_candidates'])}"
        "</section>"
    )
    explanation_section = render_explanation_section(report["guide_upstream_bases"])
    kallisto_section = render_kallisto_table(report["kallisto_like_strings"])

    modality_sections = []
    modality_sections.append(
        render_modality_section(
            title="RNA Modality",
            subtitle="Barcode and UMI are shared design intervals; RNA is placed on the opposite full strand.",
            read_lengths=report["rna_read_lengths"],
            intervals=report["rna_intervals"],
            colors=colors,
            diagnostics_html="",
        )
    )

    guide_diag = (
        "<div class='diagnostics-block'>"
        "<div class='section-header'><h3>Guide Diagnostics</h3>"
        f"<p>Guide interval inferred with the same exact-match strategy as the original script. "
        f"Guide reference lengths observed in metadata: {html.escape(report['guide_reference_lengths'])}.</p></div>"
        f"{render_sequence_config_table(report['guide_config_results'])}"
        f"{render_svg_histogram(report['guide_config_results'][0].position_counts, report['guide_config_results'][0].read_length, colors['guide'])}"
        f"{render_flank_table('Guide Prefixes', report['guide_top_flanks'])}"
        f"{render_feature_string_panels('guide', report['guide_config_results'], colors['guide'])}"
        "</div>"
    )
    modality_sections.append(
        render_modality_section(
            title="Guide Modality",
            subtitle="Shared barcode/UMI plus guide interval inferred from dominant exact-match placement.",
            read_lengths=report["guide_read_lengths"],
            intervals=report["guide_intervals"],
            colors=colors,
            diagnostics_html=guide_diag,
        )
    )

    if report.get("hash_intervals"):
        hash_diag = (
            "<div class='diagnostics-block'>"
            "<div class='section-header'><h3>Hash Diagnostics</h3>"
            "<p>Hash interval inferred with the same exact-match strategy used for guide calls.</p></div>"
            f"{render_sequence_config_table(report['hash_config_results'])}"
            f"{render_svg_histogram(report['hash_config_results'][0].position_counts, report['hash_config_results'][0].read_length, colors['hash'])}"
            f"{render_feature_string_panels('hash', report['hash_config_results'], colors['hash'])}"
            "</div>"
        )
        modality_sections.append(
            render_modality_section(
                title="Hash Modality",
                subtitle="Optional hash interval, shown only when hash reads and metadata are provided.",
                read_lengths=report["hash_read_lengths"],
                intervals=report["hash_intervals"],
                colors=colors,
                diagnostics_html=hash_diag,
            )
        )

    css = """
    :root {
      --bg-top: #f2ead7;
      --bg-bottom: #f8f5ef;
      --ink: #2d241a;
      --muted: #6f6355;
      --card: rgba(255, 252, 246, 0.88);
      --line: #d7cab8;
      --accent: #0f6b62;
      --shadow: 0 18px 45px rgba(79, 59, 35, 0.09);
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: "Avenir Next", "Segoe UI", "Helvetica Neue", sans-serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(199, 129, 38, 0.12), transparent 36%),
        radial-gradient(circle at top right, rgba(17, 107, 98, 0.12), transparent 34%),
        linear-gradient(180deg, var(--bg-top), var(--bg-bottom));
    }
    .page {
      width: min(1240px, calc(100vw - 40px));
      margin: 28px auto 40px;
    }
    .hero {
      padding: 28px 30px 24px;
      border-radius: 28px;
      background: linear-gradient(135deg, rgba(255, 249, 240, 0.92), rgba(252, 244, 230, 0.82));
      box-shadow: var(--shadow);
      border: 1px solid rgba(191, 174, 149, 0.45);
    }
    .hero h1 {
      margin: 0 0 10px;
      font-size: 34px;
      line-height: 1.1;
      letter-spacing: -0.03em;
    }
    .hero p {
      margin: 0;
      font-size: 16px;
      color: var(--muted);
      max-width: 880px;
    }
    .summary-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 14px;
      margin-top: 20px;
    }
    .summary-card, .card {
      border-radius: 24px;
      background: var(--card);
      box-shadow: var(--shadow);
      border: 1px solid rgba(191, 174, 149, 0.42);
      backdrop-filter: blur(10px);
    }
    .summary-card {
      padding: 16px 18px;
    }
    .summary-label {
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.08em;
      color: var(--muted);
      margin-bottom: 8px;
    }
    .summary-value {
      font-size: 18px;
      font-weight: 700;
      line-height: 1.3;
    }
    .card {
      padding: 22px 22px 24px;
      margin-top: 18px;
    }
    .section-header {
      display: flex;
      justify-content: space-between;
      align-items: baseline;
      gap: 16px;
      margin-bottom: 14px;
      flex-wrap: wrap;
    }
    .section-header h2, .section-header h3 {
      margin: 0;
      letter-spacing: -0.02em;
    }
    .section-header p {
      margin: 0;
      color: var(--muted);
      max-width: 780px;
      font-size: 14px;
    }
    .umi-header {
      margin-top: 18px;
    }
    .spectrum-block {
      margin-top: 16px;
    }
    .sequence-row {
      display: grid;
      grid-template-columns: minmax(180px, 2fr) minmax(220px, 3fr);
      gap: 12px;
      align-items: center;
      margin-top: 8px;
    }
    .sequence-row code {
      overflow-wrap: anywhere;
      font-size: 12px;
      background: #f4efe6;
      padding: 6px 8px;
      border-radius: 8px;
      color: #43392d;
    }
    .sequence-bar-wrap {
      position: relative;
      background: #eee4d4;
      border-radius: 999px;
      min-height: 28px;
      overflow: hidden;
    }
    .sequence-bar {
      position: absolute;
      top: 0;
      bottom: 0;
      left: 0;
      border-radius: 999px;
      opacity: 0.82;
    }
    .sequence-bar-wrap span {
      position: relative;
      z-index: 1;
      display: inline-flex;
      align-items: center;
      min-height: 28px;
      padding: 0 10px;
      font-size: 12px;
      font-weight: 700;
      color: #2f281f;
    }
    @media (max-width: 720px) {
      .sequence-row {
        grid-template-columns: 1fr;
      }
    }
    .meta-line {
      margin: 8px 0 14px;
      color: var(--muted);
      font-size: 14px;
    }
    .metrics-table {
      width: 100%;
      border-collapse: collapse;
      margin-top: 12px;
      font-size: 14px;
      overflow: hidden;
    }
    .metrics-table th,
    .metrics-table td {
      text-align: left;
      padding: 10px 12px;
      border-bottom: 1px solid rgba(208, 195, 176, 0.65);
      vertical-align: top;
    }
    .metrics-table thead th {
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.07em;
      color: var(--muted);
    }
    .pill {
      display: inline-flex;
      align-items: center;
      border-radius: 999px;
      padding: 4px 10px;
      font-size: 12px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }
    .pill-high { background: rgba(20, 134, 109, 0.16); color: #0f6b62; }
    .pill-medium { background: rgba(194, 106, 17, 0.14); color: #9b540d; }
    .pill-low { background: rgba(163, 35, 67, 0.13); color: #8c1d43; }
    .pill-assumed { background: rgba(34, 98, 198, 0.14); color: #1f57b2; }
    .histogram, .read-map {
      width: 100%;
      height: auto;
      margin-top: 14px;
    }
    .empty-block {
      margin: 10px 0 0;
      color: var(--muted);
    }
    .explanation-list {
      margin: 0;
      padding-left: 20px;
      color: var(--ink);
      line-height: 1.55;
    }
    .explanation-list li + li {
      margin-top: 8px;
    }
    code {
      font-family: "SFMono-Regular", "Menlo", "Consolas", monospace;
      font-size: 12px;
      background: rgba(34, 98, 198, 0.08);
      padding: 2px 6px;
      border-radius: 6px;
      color: #1f57b2;
    }
    .diagnostics-block {
      margin-top: 18px;
    }
    .footer {
      margin-top: 18px;
      padding: 18px 20px;
      color: var(--muted);
      font-size: 13px;
      text-align: center;
    }
    @media (max-width: 900px) {
      .page { width: min(100vw - 18px, 1000px); margin: 14px auto 28px; }
      .hero { padding: 22px 18px 18px; }
      .card { padding: 18px 16px 20px; }
      .hero h1 { font-size: 28px; }
      .summary-value { font-size: 16px; }
    }
    """

    generated_at = html.escape(report["generated_at"])
    input_rows = "".join(
        f"<li><strong>{html.escape(label)}:</strong> {html.escape(path)}</li>"
        for label, path in report["input_files"]
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(report['title'])}</title>
  <style>{css}</style>
</head>
<body>
  <main class="page">
    <section class="hero">
      <h1>{html.escape(report['title'])}</h1>
      <p>One-pass seqspec-style inference for shared barcode/UMI design plus modality-specific RNA, guide, and optional hash intervals.</p>
      <div class="summary-grid">{summary_html}</div>
    </section>

    {explanation_section}
    {kallisto_section}
    {barcode_diag}
    {''.join(modality_sections)}

    <section class="card">
      <div class="section-header">
        <h2>Inputs</h2>
        <p>Generated at {generated_at}</p>
      </div>
      <ul>{input_rows}</ul>
    </section>

    <div class="footer">Report generated by seqspec_check.py</div>
  </main>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Infer shared barcode/UMI intervals once per well, then emit an HTML report with "
            "putative RNA, guide, and optional hash regions for each modality."
        )
    )
    parser.add_argument("--input-dir", help="Optional directory with test_data-like filenames for auto-discovery")
    parser.add_argument("--barcode-whitelist", help="Cell barcode whitelist used for shared barcode inference")
    parser.add_argument("--guide-metadata", help="Guide metadata TSV/CSV with a spacer/sequence column")
    parser.add_argument("--hash-metadata", help="Hash metadata TSV/CSV with a sequence column")
    parser.add_argument("--rna-r1", help="RNA read 1 FASTQ")
    parser.add_argument("--rna-r2", help="RNA read 2 FASTQ")
    parser.add_argument("--guide-r1", help="Guide read 1 FASTQ")
    parser.add_argument("--guide-r2", help="Guide read 2 FASTQ")
    parser.add_argument("--hash-r1", help="Optional hash read 1 FASTQ")
    parser.add_argument("--hash-r2", help="Optional hash read 2 FASTQ")
    parser.add_argument(
        "--barcode-source",
        choices=["auto", "rna", "guide", "hash"],
        default="auto",
        help="Which modality should drive shared barcode/UMI inference",
    )
    parser.add_argument(
        "--barcode-sample-reads",
        type=int,
        default=BARCODE_SAMPLE_READS_DEFAULT,
        help="Number of reads to subsample for barcode inference; defaults to the first 10k reads",
    )
    parser.add_argument(
        "--feature-sample-reads",
        type=int,
        default=FEATURE_SAMPLE_READS_DEFAULT,
        help="Number of reads to sample for guide/hash exact-match scans",
    )
    parser.add_argument(
        "--guide-upstream-bases",
        type=int,
        default=12,
        help="Return this many bases immediately before the winning guide interval; use 0 to disable",
    )
    parser.add_argument("--title", default="Seqspec Inference Report")
    parser.add_argument("--output-html", default="seqspec_inference_report.html")
    parser.add_argument("--output-json", help="Optional machine-readable summary JSON")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    resolved = resolve_inputs(args)
    whitelist, barcode_length = load_barcode_whitelist(resolved["barcode_whitelist"])
    reverse_whitelist = {get_reverse_complement(barcode) for barcode in whitelist}
    sample_reads = max(args.barcode_sample_reads, args.feature_sample_reads)
    sampled_reads = load_sampled_reads(resolved["modalities"], sample_reads)

    barcode_source, barcode_result, barcode_results_by_modality, barcode_explanation = choose_shared_barcode_call(
        sampled_reads=sampled_reads,
        whitelist=whitelist,
        reverse_whitelist=reverse_whitelist,
        barcode_length=barcode_length,
        requested_source=args.barcode_source,
    )
    shared_barcode = build_barcode_interval(
        barcode_result,
        barcode_length=barcode_length,
        source=f"{barcode_source} whitelist scan",
    )
    if shared_barcode is None:
        raise ValueError("Could not infer a shared barcode interval from the sampled reads.")

    barcode_source_reads = {
        "R1": sampled_reads[barcode_source]["R1"][: args.barcode_sample_reads],
        "R2": sampled_reads[barcode_source]["R2"][: args.barcode_sample_reads],
    }
    umi_source_reads = {
        "R1": sampled_reads[barcode_source]["R1"],
        "R2": sampled_reads[barcode_source]["R2"],
    }
    barcode_result_for_umi = scan_barcode_config(
        modality=barcode_source,
        read_label=barcode_result.read_label,
        strand=barcode_result.strand,
        reads=barcode_source_reads[barcode_result.read_label],
        whitelist=whitelist,
        reverse_whitelist=reverse_whitelist,
        barcode_length=barcode_length,
    )
    umi_call, umi_candidates, umi_support_reads = predict_umi(
        source_modality=barcode_source,
        source_reads=umi_source_reads,
        barcode_call=barcode_result_for_umi,
        whitelist=whitelist,
        barcode_length=barcode_length,
    )

    if "rna" not in sampled_reads:
        raise ValueError("RNA modality is required for the RNA interval section.")
    rna_lengths = {"R1": sampled_reads["rna"]["R1_length"], "R2": sampled_reads["rna"]["R2_length"]}
    rna_interval = make_rna_interval(shared_barcode, rna_lengths)
    rna_intervals = [clamp_interval(shared_barcode, rna_lengths[shared_barcode.read_label])]
    if umi_call is not None:
        rna_intervals.append(clamp_interval(umi_call, rna_lengths[umi_call.read_label]))
    rna_intervals.append(rna_interval)

    if "guide" not in sampled_reads:
        raise ValueError("Guide modality is required for the guide interval section.")
    guide_sequences = load_feature_sequences(resolved["guide_metadata"])
    guide_reference_lengths = ", ".join(
        f"{length} bp x {count}" for length, count in Counter(len(seq) for seq in guide_sequences).most_common()
    )
    guide_winner, guide_all_results = call_modality_features(
        modality="guide",
        reads={
            "R1": sampled_reads["guide"]["R1"][: args.feature_sample_reads],
            "R2": sampled_reads["guide"]["R2"][: args.feature_sample_reads],
        },
        feature_name="guide",
        feature_sequences=guide_sequences,
    )
    guide_call = build_sequence_interval(guide_winner)
    if guide_call is None:
        raise ValueError("Could not infer a dominant guide interval from the sampled reads.")
    guide_prefix_call = build_guide_prefix_interval(guide_winner, args.guide_upstream_bases)
    guide_lengths = {"R1": sampled_reads["guide"]["R1_length"], "R2": sampled_reads["guide"]["R2_length"]}
    guide_intervals = [clamp_interval(shared_barcode, guide_lengths[shared_barcode.read_label])]
    if umi_call is not None:
        guide_intervals.append(clamp_interval(umi_call, guide_lengths[umi_call.read_label]))
    guide_intervals.append(clamp_interval(guide_call, guide_lengths[guide_call.read_label]))
    if guide_prefix_call is not None:
        guide_intervals.append(clamp_interval(guide_prefix_call, guide_lengths[guide_prefix_call.read_label]))

    hash_intervals: List[IntervalCall] = []
    hash_all_results: List[SequenceConfigResult] = []
    hash_lengths: Dict[str, int] = {}
    hash_call: Optional[IntervalCall] = None
    if "hash" in sampled_reads:
        hash_sequences = load_feature_sequences(resolved["hash_metadata"])
        hash_winner, hash_all_results = call_modality_features(
            modality="hash",
            reads={
                "R1": sampled_reads["hash"]["R1"][: args.feature_sample_reads],
                "R2": sampled_reads["hash"]["R2"][: args.feature_sample_reads],
            },
            feature_name="hash",
            feature_sequences=hash_sequences,
        )
        hash_call = build_sequence_interval(hash_winner)
        hash_lengths = {"R1": sampled_reads["hash"]["R1_length"], "R2": sampled_reads["hash"]["R2_length"]}
        hash_intervals = [clamp_interval(shared_barcode, hash_lengths[shared_barcode.read_label])]
        if umi_call is not None:
            hash_intervals.append(clamp_interval(umi_call, hash_lengths[umi_call.read_label]))
        if hash_call is not None:
            hash_intervals.append(clamp_interval(hash_call, hash_lengths[hash_call.read_label]))

    input_files = [("Barcode whitelist", resolved["barcode_whitelist"])]
    if resolved.get("guide_metadata"):
        input_files.append(("Guide metadata", resolved["guide_metadata"]))
    if resolved.get("hash_metadata") and "hash" in sampled_reads:
        input_files.append(("Hash metadata", resolved["hash_metadata"]))
    for modality in ("rna", "guide", "hash"):
        if modality in sampled_reads:
            display_name = sampled_reads[modality]["display_name"]
            input_files.append((f"{display_name} R1", sampled_reads[modality]["R1_path"]))
            input_files.append((f"{display_name} R2", sampled_reads[modality]["R2_path"]))

    umi_text = "not found"
    umi_note = "No adjacent 10/12 nt UMI-like window survived the variability heuristic."
    if umi_call is not None:
        umi_text = f"{umi_call.read_label} {umi_call.strand} {format_interval(umi_call)}"
        umi_note = (
            f"UMI predicted from {umi_support_reads} barcode-supporting reads. "
            "Candidates were restricted to 10/12 nt windows adjacent to the inferred barcode. "
            f"{umi_call.note}"
        )

    report = {
        "title": args.title,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "input_files": input_files,
        "barcode_sample_reads": args.barcode_sample_reads,
        "guide_upstream_bases": args.guide_upstream_bases,
        "shared_barcode_source": f"{barcode_source} ({shared_barcode.read_label} {shared_barcode.strand})",
        "shared_barcode_interval": f"{shared_barcode.read_label} {shared_barcode.strand} {format_interval(shared_barcode)}",
        "shared_umi_interval": umi_text,
        "rna_interval": f"{rna_interval.read_label} {rna_interval.strand} {format_interval(rna_interval)}",
        "barcode_explanation": barcode_explanation,
        "barcode_config_results": barcode_results_by_modality[barcode_source],
        "umi_note": umi_note,
        "umi_candidates": umi_candidates,
        "kallisto_like_strings": [],
        "rna_read_lengths": rna_lengths,
        "rna_intervals": rna_intervals,
        "guide_read_lengths": guide_lengths,
        "guide_intervals": guide_intervals,
        "guide_config_results": guide_all_results,
        "guide_reference_lengths": guide_reference_lengths,
        "guide_top_flanks": guide_winner.top_flanks,
        "hash_read_lengths": hash_lengths,
        "hash_intervals": hash_intervals,
        "hash_config_results": hash_all_results,
    }

    rna_kallisto = build_kallisto_like_string(shared_barcode, umi_call, rna_interval, use_full_feature_read=True)
    if rna_kallisto is not None:
        rna_kallisto["label"] = "rna"
        report["kallisto_like_strings"].append(rna_kallisto)
    guide_kallisto = build_kallisto_like_string(shared_barcode, umi_call, guide_call)
    if guide_kallisto is not None:
        guide_kallisto["label"] = "guide"
        report["kallisto_like_strings"].append(guide_kallisto)
    if hash_call is not None:
        hash_kallisto = build_kallisto_like_string(shared_barcode, umi_call, hash_call)
        if hash_kallisto is not None:
            hash_kallisto["label"] = "hash"
            report["kallisto_like_strings"].append(hash_kallisto)

    html_report = build_html_report(report)
    with open(args.output_html, "w", encoding="utf-8") as handle:
        handle.write(html_report)

    if args.output_json:
        json_payload = {
            "shared_barcode": asdict(shared_barcode),
            "shared_umi": asdict(umi_call) if umi_call is not None else None,
            "rna_interval": asdict(rna_interval),
            "guide_interval": asdict(guide_call),
            "guide_prefix_interval": asdict(guide_prefix_call) if guide_prefix_call is not None else None,
            "hash_intervals": [asdict(interval) for interval in hash_intervals],
            "barcode_config_results": [asdict(item) for item in barcode_results_by_modality[barcode_source]],
            "guide_config_results": [asdict(item) for item in guide_all_results],
            "guide_reference_lengths": guide_reference_lengths,
            "guide_top_flanks": guide_winner.top_flanks,
            "hash_config_results": [asdict(item) for item in hash_all_results],
            "umi_candidates": [asdict(item) for item in umi_candidates],
            "kallisto_like_strings": report["kallisto_like_strings"],
        }
        with open(args.output_json, "w", encoding="utf-8") as handle:
            json.dump(json_payload, handle, indent=2)

    print(f"Saved HTML report to {args.output_html}")
    if args.output_json:
        print(f"Saved JSON summary to {args.output_json}")


if __name__ == "__main__":
    main()
