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
import re

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from chemtools.reagent import normalize_reaction_type
from chemtools.recommend.utils import canonical_family


_DEFAULT_INPUT = "data/HTE_db/HTE_0.csv"
_DEFAULT_OUTPUT = "data/HTE_db/HTE_0.jsonl"
_PROGRESS_EVERY = 5000

_ORGANIC_COMPOUNDS_FILE = (
    ROOT / "chemtools" / "taxonomy" / "data" / "organic_compounds.v1.3.json"
)

_HTE_REACTION_OVERRIDES = {
    "c_n_coupling": "c_n_cross_coupling",
    "cn_coupling": "c_n_cross_coupling",
    "c_o_coupling": "c_o_coupling",
    "co_coupling": "c_o_coupling",
    "suzuki": "suzuki_miyaura",
    "suzuki_in_situ": "suzuki_miyaura",
    "sonogashira": "sonogashira",
    "heck": "heck",
    "amide_formation": "amide_coupling",
    "arylation_acidic_c_h": "c_h_activation",
    "ch_activation": "c_h_activation",
    "condensation": "condensation",
    "cyclization": "cyclization",
    "borylation_miyaura": "miyaura_borylation",
    "negishi_in_situ": "negishi_coupling",
    "negishi": "negishi_coupling",
}

_HTE_ORGANIC_OVERRIDES = {
    "rnh2-a-branch": "Any-NH2",
    "rnh2 a-branch": "Any-NH2",
    "r2nh-a-branch": "Any-NHR",
    "r2nh a-branch": "Any-NHR",
    "rnh2": "Any-NH2",
    "r2nh": "Any-NHR",
    "rconh2": "Any-CONHR",
    "rconhr": "Any-CONHR",
    "rco2h or m": "Any-CO2H",
    "rco2h": "Any-CO2H",
    "rco2r": "Any-CO2R",
    "alkyl-cl": "R-Cl",
    "alkyl-br": "R-Br",
    "alkyl-i": "R-I",
    "alkyl-oh a-branch": "R2CH-OH",
    "roh-a-branch": "R2CH-OH",
    "alkyl-oh primary": "RCH2-OH",
    "alkene": "Any-Alkene",
    "arh": "Ar-H",
    "arnh2": "Ar-NH2",
    "arnhr": "Ar-NHR",
    "arb(or)2": "Ar-B(OR)2",
    "arbf3k": "Ar-BF3K",
    "aroh": "Ar-OH",
    "aroso2r": "Ar-OTs",
    "alkyl-b(oh)2": "R-B(OH)2",
    "alkyl-b(or)2": "Any-B(OR)2",
    "urea": "Any-Urea",
    "lactam": "Any-Lactam",
    "rcn": "Any-CN",
    "r2co": "Any-CO-R",
    "r2cnr": "Any-Imine",
    "rch2pph3x": "Any-Phosphonium",
}

_CATALYST_TYPE_RULES = [
    ("Pd", re.compile(r"Pd(?![a-z])"), ["palladium"]),
    ("Cu", re.compile(r"Cu(?![a-z])"), ["copper"]),
    ("Ni", re.compile(r"Ni(?![a-z])"), ["nickel"]),
    ("Fe", re.compile(r"Fe(?![a-z])"), ["iron"]),
    ("Co", re.compile(r"Co(?![a-z])"), ["cobalt"]),
    ("Ru", re.compile(r"Ru(?![a-z])"), ["ruthenium"]),
    ("Rh", re.compile(r"Rh(?![a-z])"), ["rhodium"]),
    ("Ir", re.compile(r"Ir(?![a-z])"), ["iridium"]),
    ("Pt", re.compile(r"Pt(?![a-z])"), ["platinum"]),
    ("Au", re.compile(r"Au(?![a-z])"), ["gold"]),
    ("Ag", re.compile(r"Ag(?![a-z])"), ["silver"]),
    ("Zn", re.compile(r"Zn(?![a-z])"), ["zinc"]),
    ("Mn", re.compile(r"Mn(?![a-z])"), ["manganese"]),
    ("Cr", re.compile(r"Cr(?![a-z])"), ["chromium"]),
    ("Ti", re.compile(r"Ti(?![a-z])"), ["titanium"]),
    ("Mo", re.compile(r"Mo(?![a-z])"), ["molybdenum"]),
    ("W", re.compile(r"W(?![a-z])"), ["tungsten"]),
    ("Al", re.compile(r"Al(?![a-z])"), ["aluminum", "aluminium"]),
    ("Sn", re.compile(r"Sn(?![a-z])"), ["tin", "stannous", "stannic"]),
    ("Mg", re.compile(r"Mg(?![a-z])"), ["magnesium"]),
]


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
    raw_key = raw.lower()
    raw_slug = re.sub(r"[^0-9a-z]+", "_", raw_key).strip("_")
    override = _HTE_REACTION_OVERRIDES.get(raw_key) or _HTE_REACTION_OVERRIDES.get(raw_slug)
    if override:
        return override, canonical_family(override)
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
    override_map: Dict[str, str],
) -> Optional[str]:
    if not token:
        return None

    for candidate in _candidate_variants(token):
        key = candidate.lower()
        override = override_map.get(key)
        if override:
            return override
        mapped = alias_map.get(key)
        if mapped:
            return mapped
    return None


def _detect_catalyst_types(catalysts: Iterable[str]) -> List[str]:
    found: Set[str] = set()
    has_text = False
    for catalyst in catalysts:
        text = str(catalyst).strip()
        if not text:
            continue
        has_text = True
        text_lower = text.lower()
        if "organocatalyst" in text_lower or "organic catalyst" in text_lower:
            found.add("organocatalyst")
        for symbol, pattern, names in _CATALYST_TYPE_RULES:
            if pattern.search(text):
                found.add(symbol)
                continue
            if any(name in text_lower for name in names):
                found.add(symbol)

    if has_text and not found:
        found.add("organocatalyst")
    return sorted(found)


def _load_organic_compounds_aliases() -> Dict[str, str]:
    if not _ORGANIC_COMPOUNDS_FILE.exists():
        return {}
    payload = json.loads(_ORGANIC_COMPOUNDS_FILE.read_text(encoding="utf-8"))
    aliases: Dict[str, str] = {}
    for entry in payload.get("compounds", []):
        compound_id = entry.get("id")
        if not compound_id:
            continue
        aliases[compound_id.lower()] = compound_id
        name = entry.get("name")
        if name:
            aliases[str(name).lower()] = compound_id
        for alias in entry.get("aliases", []):
            if alias:
                aliases[str(alias).lower()] = compound_id
    return aliases


def build_hte_jsonl(
    input_csv: Path,
    output_jsonl: Path,
    *,
    max_rows: Optional[int] = None,
) -> None:
    alias_map = _load_organic_compounds_aliases()
    override_map = dict(_HTE_ORGANIC_OVERRIDES)

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
                if not reaction_type:
                    if reaction_raw:
                        unmapped_reactions.add(reaction_raw)
                    reaction_type = "Unknown"

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
                        override_map=override_map,
                    )
                    if canonical:
                        reactant_types.append(canonical)
                    else:
                        unmapped_types.add(token)

                reactant_types = _unique_list(reactant_types)

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
                catalyst_types = _detect_catalyst_types(conditions.get("catalyst", []))

                metrics = {
                    "area_total_reduced": _to_float(row.get("AREA_TOTAL_REDUCED")),
                    "z_score": _to_float(row.get("z-Score")),
                }
                metrics = {k: v for k, v in metrics.items() if v is not None}

                entry = {
                    "reaction_type": reaction_type,
                    "reactant_types": reactant_types,
                    "catalyst_type": catalyst_types,
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
