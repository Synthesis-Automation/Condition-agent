import argparse
import json
import sys
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.featurizers.structural import featurize_molecule


@lru_cache(maxsize=20000)
def _cached_featurize(smiles: str) -> Dict[str, Any]:
    return featurize_molecule(smiles)


def _normalize_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        items = [str(v).strip() for v in value if str(v).strip()]
        return " / ".join(items)
    text = str(value).strip()
    if not text:
        return ""
    text = text.replace(" OR ", " / ").replace(" or ", " / ")
    text = text.replace("; ", " / ")
    return " ".join(text.split())


def _get_condition_value(conditions: Dict[str, Any], keys: Iterable[str]) -> str:
    for key in keys:
        value = conditions.get(key)
        if value:
            return _normalize_value(value)
    return ""


def _primary_motif(smiles: str) -> str:
    try:
        analysis = _cached_featurize(smiles)
    except Exception:
        return ""
    motifs = [m.get("compound_id", "") for m in analysis.get("motifs", []) if m.get("compound_id")]
    return motifs[0] if motifs else ""


def _derive_primary_reactants(reference_reactions: List[str]) -> Tuple[str, str]:
    counts: Dict[int, Counter] = defaultdict(Counter)
    for rxn in reference_reactions:
        if ">>" not in rxn:
            continue
        reactants_part = rxn.split(">>", 1)[0]
        reactants = [r for r in reactants_part.split(".") if r]
        for idx, r_smiles in enumerate(reactants[:2]):
            motif = _primary_motif(r_smiles)
            if motif:
                counts[idx][motif] += 1

    reactant_1 = counts[0].most_common(1)[0][0] if counts.get(0) else ""
    reactant_2 = counts[1].most_common(1)[0][0] if counts.get(1) else ""
    return reactant_1, reactant_2


def _extract_expr(data: Any) -> str:
    if isinstance(data, str):
        return data.strip()
    if isinstance(data, dict):
        expr = data.get("expr")
        if isinstance(expr, str):
            return expr.strip()
    return ""


def convert_rule_db_to_hte_csv(
    input_path: str,
    output_path: str,
    *,
    yield_value: float = 90.0,
    z_score_value: float = 2.0,
) -> None:
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    reaction = data.get("reaction", {})
    metadata = data.get("metadata", {})
    reaction_family = reaction.get("family") or metadata.get("id") or input_path.stem
    reference_reactions = reaction.get("reference_reactions") or []
    reaction_notes = reaction.get("notes", "")
    applies_if_expr = _extract_expr(data.get("applies_if"))

    reactant_1, reactant_2 = _derive_primary_reactants(reference_reactions)

    rows: List[Dict[str, Any]] = []
    rule_blocks = []

    default_rule = data.get("default_rule")
    if isinstance(default_rule, dict):
        rule_blocks.append(("default", default_rule))

    for rule in data.get("base_rules", []) or []:
        if isinstance(rule, dict):
            rule_blocks.append(("base", rule))

    for rule_type, rule in rule_blocks:
        conditions = rule.get("conditions", {}) if isinstance(rule.get("conditions"), dict) else {}

        catalyst = _get_condition_value(
            conditions,
            [
                "pd_precatalyst",
                "pd_source",
                "precatalyst",
                "catalyst",
                "cu_source",
                "ru_catalyst",
                "metal_catalyst",
                "ni_source",
                "catalyst_source",
            ],
        )
        ligand = _get_condition_value(conditions, ["ligand", "ligands"])
        base = _get_condition_value(conditions, ["base"])
        solvent = _get_condition_value(conditions, ["solvent", "solvents"])
        additive = _get_condition_value(
            conditions,
            ["additive", "additives", "oxidant", "cocatalyst", "co_catalyst"],
        )

        reactant_features_expr = _extract_expr(rule.get("reactant_features"))
        if not reactant_features_expr and isinstance(rule.get("reactant_features"), dict):
            reactant_features_expr = json.dumps(rule.get("reactant_features"), ensure_ascii=True)

        row = {
            # HTE canonical columns
            "reaction_type": reaction_family,
            "yield": float(yield_value),
            "z_score": float(z_score_value),
            "reactant_1": reactant_1,
            "reactant_2": reactant_2,
            "catalyst": catalyst,
            "ligand": ligand,
            "base": base,
            "solvent": solvent,
            "additive": additive,

            # Rule metadata
            "source": "rule",
            "rule_id": rule.get("id", ""),
            "rule_name": rule.get("name", ""),
            "rule_type": rule_type,
            "rule_description": rule.get("description", ""),
            "metadata_id": metadata.get("id", ""),
            "metadata_name": metadata.get("name", ""),
            "metadata_version": metadata.get("version", ""),
            "applies_if": applies_if_expr,
            "reactant_features": reactant_features_expr,
            "reference_reactions": " / ".join(reference_reactions),
            "reaction_notes": reaction_notes,
            "rule_source_file": input_path.name,
            "conditions_json": json.dumps(conditions, ensure_ascii=True),

            # Optional flattened fields
            "temperature_C": _normalize_value(conditions.get("temperature_C")),
            "time_h": _normalize_value(conditions.get("time_h")),
            "atmosphere": _normalize_value(conditions.get("atmosphere")),
            "catalyst_loading_molpct": _normalize_value(
                conditions.get("catalyst_loading_molpct") or conditions.get("pd_loading_molpct")
            ),
            "base_equiv": _normalize_value(conditions.get("base_equiv")),
            "solvent_volume": _normalize_value(
                conditions.get("solvent_volume_mL_per_mmol")
                or conditions.get("solvent_volume_vol_per_g")
            ),
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    canonical_cols = [
        "reaction_type",
        "yield",
        "z_score",
        "reactant_1",
        "reactant_2",
        "catalyst",
        "ligand",
        "base",
        "solvent",
        "additive",
    ]
    extra_cols = [col for col in df.columns if col not in canonical_cols]
    df = df[canonical_cols + extra_cols]
    df.to_csv(output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a rule DB JSON into an HTE-style CSV (single-row rules)."
    )
    parser.add_argument("--input", "-i", required=True, help="Path to rule_db_v2 JSON file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV path")
    parser.add_argument(
        "--yield",
        dest="yield_value",
        type=float,
        default=90.0,
        help="Yield value to assign to rules (default: 90.0)",
    )
    parser.add_argument(
        "--z-score",
        dest="z_score_value",
        type=float,
        default=2.0,
        help="Z-score value to assign to rules (default: 2.0)",
    )

    args = parser.parse_args()
    convert_rule_db_to_hte_csv(
        args.input,
        args.output,
        yield_value=args.yield_value,
        z_score_value=args.z_score_value,
    )
