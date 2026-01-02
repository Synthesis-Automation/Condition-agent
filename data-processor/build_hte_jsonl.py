#!/usr/bin/env python3
"""
Convert the HTE CSV into a compact JSONL format using taxonomy normalization.
The HTE_0 CSV stores reactant labels (not SMILES), so this script canonicalizes
those labels against the taxonomy rather than running RDKit featurization.

Example:
  python data-processor/build_hte_jsonl.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from chemtools.reagent import (
    build_reactant_lookup,
    get_reactant_type_definitions,
    normalize_reaction_type,
)
from chemtools.featurizers.analysis.reactants import CSV_REACTANT_OVERRIDES
from chemtools.recommend.utils import canonical_family


_DEFAULT_INPUT = "data/HTE_db/HTE_0.csv"
_DEFAULT_OUTPUT = "data/HTE_db/HTE_0.jsonl"
_PROGRESS_EVERY = 5000


def _split_values(value: Optional[str]) -> List[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text:
        return []
    parts = [part.strip() for part in text.split(",")]
    return [part for part in parts if part]


def _unique_list(values: Iterable[str]) -> List[str]:
    seen: Set[str] = set()
    out: List[str] = []
    for value in values:
        if not value:
            continue
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _normalize_reaction(raw: str) -> Tuple[Optional[str], Optional[str]]:
    raw = (raw or "").strip()
    if not raw:
        return None, None
    normalized = normalize_reaction_type(raw)
    family = canonical_family(normalized or raw)
    return normalized or None, family or None


def _candidate_variants(token: str) -> Sequence[str]:
    token = token.strip()
    if not token:
        return []
    variants = [token]
    compact = token.replace(" ", "")
    if compact != token:
        variants.append(compact)
    if token.startswith("Ar-"):
        variants.append("Ar" + token[3:])
    if token.startswith("HetAr-"):
        variants.append("HetAr" + token[6:])
    return variants


def _canonical_member(
    token: str,
    *,
    alias_map: Dict[str, str],
    id_to_category: Dict[str, str],
    override_map: Dict[str, str],
    category_ids: Set[str],
) -> Optional[str]:
    if not token:
        return None

    for candidate in _candidate_variants(token):
        if candidate in id_to_category:
            return candidate
        if candidate in category_ids:
            return candidate
        mapped = alias_map.get(candidate.lower())
        if mapped:
            return mapped
        override = override_map.get(candidate.lower())
        if override:
            if override in id_to_category:
                return override
            if override in category_ids:
                return override
            mapped = alias_map.get(override.lower())
            if mapped:
                return mapped
    return None


def _categories_from_types(
    tokens: Iterable[str],
    *,
    id_to_category: Dict[str, str],
) -> List[str]:
    categories: List[str] = []
    for token in tokens:
        token = token.strip()
        if not token:
            continue
        categories.append(id_to_category.get(token, token))
    return _unique_list(categories)


def build_hte_jsonl(
    input_csv: Path,
    output_jsonl: Path,
    *,
    max_rows: Optional[int] = None,
) -> None:
    alias_map, id_to_category = build_reactant_lookup()
    override_map = {k.lower(): v for k, v in CSV_REACTANT_OVERRIDES.items()}
    category_ids = set(get_reactant_type_definitions().keys())

    total = 0
    written = 0
    unmapped_types: Set[str] = set()
    unmapped_reactions: Set[str] = set()

    with input_csv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        output_jsonl.parent.mkdir(parents=True, exist_ok=True)
        with output_jsonl.open("w", encoding="utf-8") as out_handle:
            for row_idx, row in enumerate(reader):
                total += 1
                if max_rows and total > max_rows:
                    break

                reaction_raw = (row.get("Reaction_Type_Standardized") or "").strip()
                reaction_type, _ = _normalize_reaction(reaction_raw)
                if not reaction_type and reaction_raw:
                    unmapped_reactions.add(reaction_raw)
                    reaction_type = reaction_raw

                reactant_a_tokens = _split_values(row.get("Reactant_A"))
                if not reactant_a_tokens:
                    reactant_a_tokens = _split_values(row.get("Reactant_A_Type"))
                reactant_b_tokens = _split_values(row.get("Reactant_B"))
                if not reactant_b_tokens:
                    reactant_b_tokens = _split_values(row.get("Reactant_B_Type"))
                reactant_tokens = _unique_list(reactant_a_tokens + reactant_b_tokens)

                reactant_types: List[str] = []
                for token in reactant_tokens:
                    canonical = _canonical_member(
                        token,
                        alias_map=alias_map,
                        id_to_category=id_to_category,
                        override_map=override_map,
                        category_ids=category_ids,
                    )
                    if canonical:
                        reactant_types.append(canonical)
                    else:
                        unmapped_types.add(token)

                reactant_types = _unique_list(reactant_types)

                reactant_categories = _categories_from_types(
                    reactant_types,
                    id_to_category=id_to_category,
                )

                conditions = {
                    "catalyst": _unique_list(_split_values(row.get("Catalyst"))),
                    "ligand": _unique_list(_split_values(row.get("Ligand"))),
                    "base": _unique_list(_split_values(row.get("Base"))),
                    "solvent": _unique_list(_split_values(row.get("Solvent"))),
                    "secondary_solvent": _unique_list(_split_values(row.get("Secondary Solvent"))),
                    "additive": _unique_list(_split_values(row.get("Additive"))),
                    "coupling_reagent": _unique_list(_split_values(row.get("Coupling Reagent"))),
                }
                conditions = {k: v for k, v in conditions.items() if v}

                metrics = {
                    "area_total_reduced": _to_float(row.get("AREA_TOTAL_REDUCED")),
                    "z_score": _to_float(row.get("z-Score")),
                }
                metrics = {k: v for k, v in metrics.items() if v is not None}

                entry = {
                    "reaction_type": reaction_type,
                    "reactant_types": reactant_types,
                    "reactant_categories": sorted(reactant_categories),
                    "conditions": conditions,
                    "metrics": metrics,
                }

                out_handle.write(json.dumps(entry, ensure_ascii=True) + "\n")
                written += 1

                if written % _PROGRESS_EVERY == 0:
                    print(f"Processed {written} rows...")

    print(f"Done. Wrote {written} JSONL records to {output_jsonl}")
    if unmapped_types:
        preview = ", ".join(sorted(unmapped_types)[:10])
        print(f"Unmapped reactant tokens: {len(unmapped_types)} (sample: {preview})")
    if unmapped_reactions:
        preview = ", ".join(sorted(unmapped_reactions)[:10])
        print(f"Unmapped reaction types: {len(unmapped_reactions)} (sample: {preview})")


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert HTE CSV to JSONL")
    parser.add_argument("--input", default=_DEFAULT_INPUT)
    parser.add_argument("--output", default=_DEFAULT_OUTPUT)
    parser.add_argument("--max-rows", type=int, default=None)
    args = parser.parse_args()

    build_hte_jsonl(
        Path(args.input),
        Path(args.output),
        max_rows=args.max_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
