#!/usr/bin/env python3
"""
Triage low-information reactions that produce `none -> none`-style outputs.

The script focuses on rows where deterministic featurization cannot extract
reacted/formed motifs and returns no reliable bond evidence. It classifies such
rows into actionable data-quality buckets.

Example:
  python scripts/triage_no_motif_evidence.py \
    --input results/reaction_dataset_10x100_sample.csv \
    --reaction-column reaction_smiles \
    --output-json results/no_motif_evidence_triage.10x100.json \
    --output-csv results/no_motif_evidence_triage_rows.10x100.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util import rdkit_helpers


REACTION_COL_CANDIDATES: Tuple[str, ...] = (
    "reaction_smiles",
    "reaction",
    "rxn_smiles",
    "rxn",
    "smiles",
)
SOURCE_LABEL_CANDIDATES: Tuple[str, ...] = (
    "reaction_type_standardized",
    "reaction_type",
    "source_reaction_label",
    "source_reaction_family",
)

METAL_TOKENS = {
    "Li", "Na", "K", "Mg", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Ga", "In", "Sn", "Sb", "Bi", "Al", "B",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
}

DEFAULT_OUTPUT_JSON = "results/no_motif_evidence_triage.json"
DEFAULT_OUTPUT_CSV = "results/no_motif_evidence_triage_rows.csv"


def _find_column(fieldnames: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    normalized = {str(name).strip().lower(): str(name).strip() for name in fieldnames}
    for candidate in candidates:
        match = normalized.get(str(candidate).strip().lower())
        if match:
            return match
    return None


def _split_reaction_sides(reaction_smiles: str) -> Tuple[List[str], List[str]]:
    text = str(reaction_smiles or "").strip()
    if not text:
        return [], []
    if ">>" in text:
        left, right = text.split(">>", 1)
        return [x for x in left.split(".") if x], [x for x in right.split(".") if x]
    parts = text.split(">")
    if len(parts) == 3:
        return [x for x in parts[0].split(".") if x], [x for x in parts[2].split(".") if x]
    return [x for x in text.split(".") if x], []


def _iter_rows(
    input_path: Path,
    *,
    reaction_column: Optional[str],
    limit: Optional[int] = None,
) -> Iterator[Dict[str, Any]]:
    if input_path.is_dir():
        files = sorted(path for path in input_path.rglob("*.csv") if path.is_file())
    else:
        files = [input_path]

    emitted = 0
    for path in files:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = list(reader.fieldnames or [])
            if not fieldnames:
                continue
            rxn_col = reaction_column or _find_column(fieldnames, REACTION_COL_CANDIDATES)
            if not rxn_col:
                continue
            label_col = _find_column(fieldnames, SOURCE_LABEL_CANDIDATES)
            for row_num, row in enumerate(reader, start=2):
                if limit is not None and emitted >= max(0, int(limit)):
                    return
                if not isinstance(row, dict):
                    continue
                rxn = str(row.get(rxn_col) or "").strip()
                if not rxn:
                    continue
                emitted += 1
                yield {
                    "reaction_smiles": rxn,
                    "source_file": path.name,
                    "source_path": str(path),
                    "row_number": row_num,
                    "source_reaction_label": str(row.get(label_col) or "").strip() if label_col else "",
                }


def _extract_quality(result: Dict[str, Any]) -> Dict[str, Any]:
    events = result.get("reaction_events")
    if isinstance(events, dict):
        quality = events.get("reaction_key_quality")
        if isinstance(quality, dict) and quality:
            return quality
    detection = result.get("detection")
    if isinstance(detection, dict):
        quality = detection.get("reaction_key_quality")
        if isinstance(quality, dict) and quality:
            return quality
    return {}


def _reaction_type(result: Dict[str, Any]) -> str:
    value = result.get("reaction_type")
    if isinstance(value, dict):
        text = str(value.get("reaction_type") or value.get("name") or "").strip()
        return text or "Unknown"
    text = str(value or "").strip()
    return text or "Unknown"


def _aggregates(result: Dict[str, Any]) -> Dict[str, Any]:
    payload = result.get("aggregates")
    return payload if isinstance(payload, dict) else {}


def _contains_metal_token(text: str) -> Tuple[bool, List[str]]:
    tokens = re.findall(r"\[([A-Z][a-z]?)", str(text or ""))
    present = sorted({tok for tok in tokens if tok in METAL_TOKENS})
    return bool(present), present


def _count_unparseable(smiles_list: Iterable[str]) -> int:
    if not rdkit_helpers.rdkit_available():
        return 0
    failures = 0
    for smi in smiles_list:
        try:
            mol = rdkit_helpers.parse_smiles(str(smi))
        except Exception:
            mol = None
        if mol is None:
            failures += 1
    return failures


def _is_no_motif_evidence_case(result: Dict[str, Any]) -> bool:
    aggr = _aggregates(result)
    reacted = [str(m) for m in (aggr.get("reacted_motifs") or []) if str(m).strip()]
    formed = [
        str(m)
        for m in (aggr.get("formed_motifs_all") or aggr.get("formed_motifs") or [])
        if str(m).strip()
    ]
    if reacted or formed:
        return False
    quality = _extract_quality(result)
    reasons = {str(r).strip() for r in (quality.get("reasons") or []) if str(r).strip()}
    if {
        "missing_formed_bond_and_product_motif_evidence",
        "missing_bond_key",
    }.issubset(reasons):
        return True
    rt = _reaction_type(result).strip().lower()
    return rt in {"unknown", "unresolved:nomotifevidence"}


def classify_no_motif_evidence_case(
    reaction_smiles: str,
    result: Dict[str, Any],
) -> Tuple[str, Dict[str, Any]]:
    reactants, products = _split_reaction_sides(reaction_smiles)
    reaction_key = str(result.get("reaction_key") or "").strip()
    quality = _extract_quality(result)
    reasons = [str(r).strip() for r in (quality.get("reasons") or []) if str(r).strip()]
    quality_level = str(quality.get("level") or "").strip().lower()
    quality_score = float(quality.get("score_0_1") or 0.0)
    rt = _reaction_type(result)

    detail: Dict[str, Any] = {
        "reaction_type": rt,
        "reaction_key_empty": not bool(reaction_key),
        "quality_level": quality_level or "unknown",
        "quality_score": quality_score,
        "quality_reasons": reasons,
        "reactant_count": len(reactants),
        "product_count": len(products),
    }

    if not reactants or not products:
        detail["bucket_reason"] = "missing_reaction_side"
        return "missing_side_or_malformed_reaction_smiles", detail

    has_metal, metal_tokens = _contains_metal_token(reaction_smiles)
    detail["metal_tokens"] = metal_tokens
    if has_metal:
        detail["bucket_reason"] = "metal_or_coordination_smiles"
        return "organometallic_or_coordination_complex", detail

    unparseable = _count_unparseable(reactants + products)
    detail["unparseable_components"] = unparseable
    if unparseable > 0:
        detail["bucket_reason"] = "rdkit_parse_failures"
        return "unparseable_components", detail

    max_len = max((len(x) for x in (reactants + products)), default=0)
    detail["max_component_length"] = max_len
    if max_len >= 160:
        detail["bucket_reason"] = "very_long_component_without_motif_hits"
        return "very_large_or_complex_structure", detail

    if set(reactants) == set(products):
        detail["bucket_reason"] = "identical_sides_text"
        return "no_structural_change_or_unmapped_transform", detail

    detail["bucket_reason"] = "no_motif_no_bond_evidence"
    return "unsupported_or_missing_motif_coverage", detail


def build_triage_report(
    rows: Iterable[Dict[str, Any]],
    *,
    progress_every: int = 0,
) -> Dict[str, Any]:
    run_options = {
        "detailed": False,
        "include_reaction_type": True,
        "confirm_coupling_products": True,
    }

    processed = 0
    selected = 0
    errors = 0

    bucket_counts: Counter[str] = Counter()
    source_file_counts: Counter[str] = Counter()
    source_label_counts: Counter[str] = Counter()
    reaction_type_counts: Counter[str] = Counter()
    quality_reason_counts: Counter[str] = Counter()
    metal_token_counts: Counter[str] = Counter()

    selected_rows: List[Dict[str, Any]] = []

    for row in rows:
        processed += 1
        if progress_every and processed % max(1, int(progress_every)) == 0:
            print(f"[progress] processed={processed} selected={selected} errors={errors}")

        rxn = str(row.get("reaction_smiles") or "").strip()
        if not rxn:
            continue
        try:
            result = featurize_reaction(rxn, options=run_options)
        except Exception as exc:
            errors += 1
            selected_rows.append(
                {
                    **row,
                    "triage_bucket": "featurization_error",
                    "reaction_type": "error",
                    "quality_level": "unknown",
                    "quality_score": 0.0,
                    "quality_reasons": "",
                    "error": str(exc),
                }
            )
            bucket_counts["featurization_error"] += 1
            continue

        if not _is_no_motif_evidence_case(result):
            continue

        selected += 1
        bucket, detail = classify_no_motif_evidence_case(rxn, result)
        bucket_counts[bucket] += 1
        source_file = str(row.get("source_file") or "").strip()
        source_label = str(row.get("source_reaction_label") or "").strip()
        if source_file:
            source_file_counts[source_file] += 1
        if source_label:
            source_label_counts[source_label] += 1
        reaction_type_counts[str(detail.get("reaction_type") or "Unknown")] += 1
        for reason in detail.get("quality_reasons") or []:
            quality_reason_counts[str(reason)] += 1
        for token in detail.get("metal_tokens") or []:
            metal_token_counts[str(token)] += 1

        selected_rows.append(
            {
                **row,
                "triage_bucket": bucket,
                "reaction_type": detail.get("reaction_type"),
                "reaction_key_empty": detail.get("reaction_key_empty"),
                "quality_level": detail.get("quality_level"),
                "quality_score": detail.get("quality_score"),
                "quality_reasons": "|".join(detail.get("quality_reasons") or []),
                "reactant_count": detail.get("reactant_count"),
                "product_count": detail.get("product_count"),
                "unparseable_components": detail.get("unparseable_components", 0),
                "metal_tokens": "|".join(detail.get("metal_tokens") or []),
                "max_component_length": detail.get("max_component_length", 0),
                "bucket_reason": detail.get("bucket_reason", ""),
            }
        )

    selected_rows.sort(
        key=lambda r: (
            str(r.get("triage_bucket") or ""),
            str(r.get("source_file") or ""),
            int(r.get("row_number") or 0),
        )
    )
    denominator = max(1, processed)
    return {
        "summary": {
            "processed_reactions": processed,
            "selected_no_motif_evidence_rows": selected,
            "selection_rate": round(selected / denominator, 4),
            "featurization_errors": errors,
        },
        "bucket_counts": bucket_counts.most_common(),
        "source_file_counts": source_file_counts.most_common(),
        "source_reaction_label_counts": source_label_counts.most_common(),
        "reaction_type_counts": reaction_type_counts.most_common(),
        "quality_reason_counts": quality_reason_counts.most_common(),
        "metal_token_counts": metal_token_counts.most_common(),
        "rows": selected_rows,
    }


def write_rows_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    fields = [
        "source_file",
        "source_path",
        "row_number",
        "source_reaction_label",
        "triage_bucket",
        "bucket_reason",
        "reaction_type",
        "reaction_key_empty",
        "quality_level",
        "quality_score",
        "quality_reasons",
        "reactant_count",
        "product_count",
        "unparseable_components",
        "metal_tokens",
        "max_component_length",
        "reaction_smiles",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Triage none->none / no-motif-evidence reactions")
    parser.add_argument("--input", required=True, help="Input CSV file or directory")
    parser.add_argument(
        "--reaction-column",
        default="",
        help="Optional explicit reaction SMILES column",
    )
    parser.add_argument("--limit", type=int, default=None, help="Optional global row limit")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=200,
        help="Print progress every N processed rows (0 to disable)",
    )
    parser.add_argument("--output-json", default=DEFAULT_OUTPUT_JSON, help="Output JSON path")
    parser.add_argument("--output-csv", default=DEFAULT_OUTPUT_CSV, help="Output rows CSV path")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input path not found: {input_path}")

    rows = _iter_rows(
        input_path,
        reaction_column=str(args.reaction_column).strip() or None,
        limit=args.limit,
    )
    report = build_triage_report(rows, progress_every=max(0, int(args.progress_every)))

    out_json = Path(args.output_json).expanduser().resolve()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    out_csv = Path(args.output_csv).expanduser().resolve()
    write_rows_csv(out_csv, report.get("rows") or [])

    summary = report.get("summary") or {}
    print(f"Triage JSON: {out_json}")
    print(f"Triage CSV: {out_csv}")
    print(f"Processed: {summary.get('processed_reactions', 0)}")
    print(f"Selected no-motif-evidence rows: {summary.get('selected_no_motif_evidence_rows', 0)}")
    print(f"Selection rate: {summary.get('selection_rate', 0)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
