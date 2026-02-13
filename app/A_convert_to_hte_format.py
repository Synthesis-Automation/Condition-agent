import csv
import json
import re
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any, Tuple, Iterable, Optional
from functools import lru_cache
import argparse
import sys

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Import chemtools components
from chemtools.featurizers.unified import featurize_molecule, featurize_reaction
from chemtools.featurizers.formatters.reaction import get_crk_options
from chemtools.featurizers.spectator_rank import rank_spectator_groups
from chemtools.smiles import normalize
from chemtools.reagent.lookup import find_reagent

@lru_cache(maxsize=10000)
def cached_featurize(smiles: str):
    """Cache featurization results to speed up processing of repeated reactants."""
    # Use substituent-level filtering to avoid redundant labels (e.g., Ar-NR2 and RCH2-NR2)
    # also ensures that Aromatic scaffolds win over Aliphatic ones due to updated priorities.
    options = {"motif_site_filter": "substituent", "detailed": True}
    return featurize_molecule(smiles, options=options)


@lru_cache(maxsize=20000)
def cached_featurize_reaction(
    smiles: str,
    llm_assist_signature: str = "",
) -> Dict[str, Any]:
    """Cache reaction featurization to keep Reaction_Key generation consistent."""
    if not smiles:
        return {}
    try:
        return featurize_reaction(
            smiles,
            options=_build_reaction_options(llm_assist_signature),
        )
    except Exception:
        return {}


def _build_reaction_options(llm_assist_signature: str = "") -> Dict[str, Any]:
    options = get_crk_options()
    if not llm_assist_signature:
        return options
    try:
        llm_assist_cfg = json.loads(llm_assist_signature)
    except Exception:
        return options
    if isinstance(llm_assist_cfg, dict):
        options["llm_assist"] = llm_assist_cfg
    return options


def _llm_assist_signature(llm_assist_options: Optional[Dict[str, Any]]) -> str:
    if not isinstance(llm_assist_options, dict):
        return ""
    payload = dict(llm_assist_options)
    if not payload.get("enabled", True):
        return ""
    payload["enabled"] = True
    try:
        return json.dumps(payload, sort_keys=True, separators=(",", ":"))
    except Exception:
        return ""


@lru_cache(maxsize=20000)
def _detect_reaction_type(
    reaction_smiles: str,
    llm_assist_signature: str = "",
) -> str:
    """
    Detect reaction type using full featurization pipeline with taxonomy-driven validation.
    
    This now uses featurize_reaction() which includes:
    - Slot-based detection
    - Product-based validation (taxonomy-driven)
    - Reacted motif pattern matching
    
    This ensures accurate detection (e.g., Suzuki_miyaura instead of Arylation_Ar_H
    when organoboron + aryl halide → biaryl pattern is present).
    """
    if not reaction_smiles:
        return ""
    try:
        result = featurize_reaction(
            reaction_smiles,
            options=_build_reaction_options(llm_assist_signature),
        )
        reaction_type = result.get("reaction_type", {})
        
        # Extract reaction type from the result (handles both dict and string formats)
        if isinstance(reaction_type, dict):
            return reaction_type.get("reaction_type", "")
        return str(reaction_type) if reaction_type else ""
        
    except Exception:
        return ""


def _extract_reaction_type_from_bundle(bundle: Dict[str, Any]) -> str:
    """Extract normalized reaction type string from a featurize_reaction bundle."""
    if not isinstance(bundle, dict):
        return ""
    reaction_type = bundle.get("reaction_type")
    if isinstance(reaction_type, dict):
        value = reaction_type.get("reaction_type") or reaction_type.get("name") or ""
        return str(value).strip()
    if reaction_type:
        return str(reaction_type).strip()
    return ""


def _build_taxonomy_gap_entry(
    rxn_bundle: Dict[str, Any],
    *,
    source_dataset: str,
    source_format: str,
    source_row_index: int,
    reaction_smiles: str,
    reaction_key: str,
    detected_reaction_type: str,
    reference: str,
) -> Optional[Dict[str, Any]]:
    """Build an export row when LLM assist flags taxonomy coverage gaps."""
    if not isinstance(rxn_bundle, dict):
        return None
    meta = rxn_bundle.get("meta")
    if not isinstance(meta, dict):
        return None
    llm_meta = meta.get("llm_assist")
    if not isinstance(llm_meta, dict):
        return None
    if not bool(llm_meta.get("taxonomy_gap_suspected", False)):
        return None

    return {
        "source_dataset": source_dataset,
        "source_format": source_format,
        "source_row_index": int(source_row_index),
        "reaction_smiles": reaction_smiles,
        "reaction_key": reaction_key,
        "detected_reaction_type": detected_reaction_type or "Unknown",
        "suggested_reaction_type": str(llm_meta.get("suggested_reaction_type") or "Unknown"),
        "suggested_confidence": llm_meta.get("suggested_confidence"),
        "decision": str(llm_meta.get("decision") or ""),
        "status": str(llm_meta.get("status") or ""),
        "validation": str(llm_meta.get("validation") or ""),
        "provider": str(llm_meta.get("provider") or ""),
        "model": str(llm_meta.get("model") or ""),
        "requires_human_review": bool(llm_meta.get("requires_human_review", False)),
        "taxonomy_gap_suspected": True,
        "taxonomy_gap_note": str(llm_meta.get("taxonomy_gap_note") or ""),
        "mechanistic_family": str(llm_meta.get("mechanistic_family") or ""),
        "mechanistic_rationale": str(llm_meta.get("mechanistic_rationale") or ""),
        "tautomer_or_representation_issue": bool(
            llm_meta.get("tautomer_or_representation_issue", False)
        ),
        "uncertainty_flags": list(llm_meta.get("uncertainty_flags") or []),
        "uncertainty_reasons": list(llm_meta.get("uncertainty_reasons") or []),
        "deterministic_checks_used": list(llm_meta.get("deterministic_checks_used") or []),
        "reference": reference or "",
        "llm_assist_meta": llm_meta,
    }


def _write_taxonomy_gap_exports(
    output_path: Path,
    gap_entries: List[Dict[str, Any]],
) -> None:
    """Write taxonomy-gap sidecar reports in JSON and CSV formats."""
    if not gap_entries:
        return

    json_path = output_path.with_name(f"{output_path.stem}.taxonomy_gaps.json")
    csv_path = output_path.with_name(f"{output_path.stem}.taxonomy_gaps.csv")

    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(gap_entries, handle, ensure_ascii=False, indent=2)

    csv_rows: List[Dict[str, Any]] = []
    for entry in gap_entries:
        row = dict(entry)
        for key in (
            "uncertainty_flags",
            "uncertainty_reasons",
            "deterministic_checks_used",
            "llm_assist_meta",
        ):
            value = row.get(key)
            if isinstance(value, (list, dict)):
                row[key] = json.dumps(value, ensure_ascii=False, sort_keys=True)
        csv_rows.append(row)
    pd.DataFrame(csv_rows).to_csv(csv_path, index=False)

    print(f"Saved {len(gap_entries)} taxonomy-gap entries to {json_path}")
    print(f"Saved {len(gap_entries)} taxonomy-gap entries to {csv_path}")


CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")
CAS_INLINE_PATTERN = re.compile(r"\b\d{2,7}-\d{2}-\d\b")

def _normalize_smiles(smiles: str) -> Optional[str]:
    if not smiles:
        return None
    info = normalize(smiles)
    return info.get("smiles_norm") or info.get("largest_smiles") or None

def _split_reagent_names(value: str) -> List[str]:
    if not value:
        return []
    parts = []
    for token in value.replace(";", "/").split("/"):
        for item in token.split(","):
            name = item.strip()
            if name:
                parts.append(name)
    return parts

def _split_csv_list(value: Any) -> List[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text:
        return []
    return [part.strip() for part in text.split(",") if part.strip()]

def _normalize_reagent_text(value: str) -> str:
    return " ".join(str(value).strip().lower().split())

def _collect_reagent_smiles(
    record: Dict[str, Any],
    reagent_db: Optional[Dict[str, Any]] = None,
    unknown_cas: Optional[set[str]] = None,
) -> set[str]:
    reagent_names: List[str] = []
    reagent_cas: List[str] = []

    def normalize_entries(values: Any) -> List[Dict[str, Any]]:
        if values is None:
            return []
        if isinstance(values, dict):
            values = [values]
        entries = []
        for entry in values:
            if isinstance(entry, dict):
                entries.append(entry)
            elif isinstance(entry, str):
                entries.append({"name": entry})
        return entries

    for entry in normalize_entries(record.get("reagents")):
        name = entry.get("name")
        if name:
            reagent_names.extend(_split_reagent_names(str(name)))
        cas = entry.get("cas")
        if cas:
            reagent_cas.append(str(cas).strip())
    
    for entry in normalize_entries(record.get("solvents")):
        name = entry.get("name")
        if name:
            reagent_names.extend(_split_reagent_names(str(name)))
        cas = entry.get("cas")
        if cas:
            reagent_cas.append(str(cas).strip())

    catalytic_system = record.get("catalytic_system") or ""
    condition_core = record.get("condition_core") or ""
    if catalytic_system:
        reagent_names.extend(_split_reagent_names(str(catalytic_system)))
    if condition_core:
        reagent_names.extend(_split_reagent_names(str(condition_core)))

    reagent_smiles = set()
    unknown = unknown_cas if unknown_cas is not None else set()

    # Use the new reagent system for lookups
    for cas in _dedupe_list(reagent_cas):
        if not cas or cas.startswith("$") or not CAS_PATTERN.match(cas):
            continue
        
        # Try finding by CAS in any role with exact CAS match verification
        hit = None
        for r_type in ['metal_catalyst', 'ligand', 'base', 'solvent', 'additive', 'acid', 'oxidant', 'reductant']:
            hit = find_reagent(cas, r_type)
            if hit and hit.get('cas') == cas:
                # Exact CAS match - use this
                break
            elif hit:
                # Partial match (possibly due to normalization collision) - keep searching
                hit = None
        
        if hit:
            smiles = hit.get('smiles')
            if smiles:
                normalized = _normalize_smiles(smiles)
                if normalized:
                    reagent_smiles.add(normalized)
        else:
            unknown.add(cas)

    for name in _dedupe_list(reagent_names):
        # Check for inline CAS
        key = name.lower()
        found_inline = False
        for match in CAS_INLINE_PATTERN.findall(key):
            found_inline = True
            # Try finding by CAS with exact match verification
            hit = None
            for r_type in ['metal_catalyst', 'ligand', 'base', 'solvent', 'additive', 'acid', 'oxidant', 'reductant']:
                hit = find_reagent(match, r_type)
                if hit and hit.get('cas') == match:
                    # Exact CAS match
                    break
                elif hit:
                    # Partial match - keep searching
                    hit = None
            if hit and hit.get('smiles'):
                normalized = _normalize_smiles(hit['smiles'])
                if normalized:
                    reagent_smiles.add(normalized)
        
        if not found_inline:
            # Try finding by name
            hit = None
            for r_type in ['metal_catalyst', 'ligand', 'base', 'solvent', 'additive', 'acid', 'oxidant', 'reductant']:
                hit = find_reagent(name, r_type)
                if hit:
                    break
            if hit and hit.get('smiles'):
                normalized = _normalize_smiles(hit['smiles'])
                if normalized:
                    reagent_smiles.add(normalized)

    return reagent_smiles

def _write_new_reagents(path: Path, cas_values: set[str]) -> None:
    if not cas_values:
        return
    existing: set[str] = set()
    if path.exists():
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                cas = (row.get("cas") or "").strip()
                if cas:
                    existing.add(cas)
    new_only = {cas.strip() for cas in cas_values if cas.strip()} - existing
    if not new_only:
        return
    merged = sorted(existing | new_only)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cas"])
        writer.writeheader()
        for cas in merged:
            writer.writerow({"cas": cas})

@lru_cache(maxsize=1)
def _load_scaffold_motif_ids() -> set[str]:
    path = PROJECT_ROOT / "chemtools" / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"
    if not path.exists():
        return set()
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)
    except Exception:
        return set()
    motifs = set()
    for entry in payload.get("compounds", []) or []:
        if not isinstance(entry, dict):
            continue
        motif_id = str(entry.get("id") or "").strip()
        if motif_id:
            motifs.add(motif_id)
    return motifs


@lru_cache(maxsize=1)
def _load_nh_heterocycle_scaffolds() -> set[str]:
    path = PROJECT_ROOT / "chemtools" / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"
    if not path.exists():
        return set()
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)
    except Exception:
        return set()
    motifs: set[str] = set()
    for entry in payload.get("compounds", []) or []:
        if not isinstance(entry, dict):
            continue
        motif_id = str(entry.get("id") or "").strip()
        description = str(entry.get("description") or "").lower()
        if motif_id and "nh-heteroaromatic" in description:
            motifs.add(motif_id)
    return motifs


def _apply_aromn_scaffold_override(
    motifs: List[str],
    context_scaffolds: List[str],
) -> List[str]:
    if "AromN-H" not in motifs or not context_scaffolds:
        return motifs
    nh_scaffolds = _load_nh_heterocycle_scaffolds()
    scaffold = next((cid for cid in context_scaffolds if cid in nh_scaffolds), "")
    if not scaffold:
        return motifs
    idx = motifs.index("AromN-H")
    filtered = [m for m in motifs if m != "AromN-H"]
    if scaffold not in filtered:
        insert_at = min(idx, len(filtered))
        filtered.insert(insert_at, scaffold)
    return filtered


def _extract_motif_ids(entries: Any) -> List[str]:
    ids: List[str] = []
    for entry in entries or []:
        if isinstance(entry, dict):
            cid = str(entry.get("id") or entry.get("compound_id") or "").strip()
        else:
            cid = str(entry).strip()
        if cid:
            ids.append(cid)
    return ids


def _extract_context_scaffolds(analysis: Dict[str, Any]) -> List[str]:
    context_ids: List[str] = []
    context_ids.extend(_extract_motif_ids(analysis.get("context_motifs", [])))
    extended = analysis.get("extended")
    if isinstance(extended, dict):
        context_ids.extend(_extract_motif_ids(extended.get("context_motifs", [])))
    return _dedupe_list(context_ids)

def _dedupe_list(values: Iterable[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for value in values:
        if not value:
            continue
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out

def _reactant_key(values: Iterable[Optional[str]]) -> str:
    """Generate a standardized key for reactant combinations."""
    items = _dedupe_list([str(v).strip() for v in values if v])
    return "|".join(sorted(items))

def _select_primary_reactant_motifs(
    reactant_data: List[Dict[str, Any]],
    reacted_set: set[str],
) -> List[str]:
    primary: List[str] = []
    for r_info in reactant_data:
        r_motifs = r_info.get("motifs") or []
        reacted_here = _dedupe_list([m for m in r_motifs if m in reacted_set])
        if reacted_here:
            primary.append("|".join(reacted_here))
        elif r_motifs:
            primary.append(r_motifs[0])
        else:
            primary.append("")
    return primary

def _motif_group_ids(values: Iterable[str]) -> set[str]:
    group_ids: set[str] = set()
    for value in values:
        token = str(value).strip()
        if not token:
            continue
        if "-" in token:
            group_id = token.split("-")[-1]
        else:
            group_id = token
        group_id = group_id.strip()
        if group_id:
            group_ids.add(group_id)
    return group_ids

_SPECTATOR_GROUP_STOPLIST = {
    "Ar",
    "R",
    "Any",
    "Alkyl",
    "Alkenyl",
    "Alkynyl",
    "H",
}
_HETERO_C_H_SCAFFOLD_GROUPS = {"HeteroAr", "AromN"}

def _group_id_from_motif_id(motif_id: str) -> str:
    text = str(motif_id).strip()
    if not text:
        return ""
    if "-" in text:
        scaffold, group_id = text.rsplit("-", 1)
        scaffold = scaffold.strip()
        group_id = group_id.strip()
        if group_id == "H" and scaffold in _HETERO_C_H_SCAFFOLD_GROUPS:
            return scaffold
        return group_id
    return text

def _collect_spectator_groups(
    reactant_data: List[Dict[str, Any]],
    spectators_set: set[str],
) -> List[str]:
    scaffold_ids = _load_scaffold_motif_ids()
    seen = set()
    groups: List[str] = []

    for r_info in reactant_data:
        motifs = r_info.get("motifs") or []
        for motif_id in motifs:
            if motif_id not in spectators_set:
                continue
            group_id = _group_id_from_motif_id(motif_id)
            if not group_id or group_id in _SPECTATOR_GROUP_STOPLIST:
                continue
            if group_id in seen:
                continue
            seen.add(group_id)
            groups.append(group_id)

    for motif_id in spectators_set:
        if motif_id in scaffold_ids and motif_id not in seen:
            seen.add(motif_id)
            groups.append(motif_id)

    return groups

def extract_reagents(record: Dict[str, Any], csv_row: Optional[Dict[str, Any]] = None) -> Dict[str, str]:
    """Extract and categorize reagents by role using the reagent system."""
    # Collect all CAS numbers from the record
    all_cas = []
    
    # Use raw CSV row if provided, otherwise try from record
    source = csv_row if csv_row else record
    
    # From CSV columns
    for cas in _split_csv_list(source.get("reagent_cas", "")):
        all_cas.append(cas)
    for cas in _split_csv_list(source.get("catalyst_cas", "")):
        all_cas.append(cas)
    for cas in _split_csv_list(source.get("solvent_cas", "")):
        all_cas.append(cas)
    
    # Classify each CAS by role
    role_items = {
        "metal_catalyst": [],
        "ligand": [],
        "base": [],
        "acid": [],
        "oxidant": [],
        "reductant": [],
        "additive": [],
        "condensation_agent": [],
        "other_reagent": [],
        "solvent": [],
    }
    
    for cas in _dedupe_list(all_cas):
        if not cas or not CAS_PATTERN.match(cas):
            continue
            
        # Try to find this CAS in the reagent system
        # Check by exact CAS match to avoid normalization issues
        hit = None
        for r_type in ['metal_catalyst', 'ligand', 'base', 'acid', 'oxidant', 'reductant', 'additive', 'condensation_agent', 'other_reagent', 'solvent']:
            hit = find_reagent(cas, r_type)
            if hit and hit.get('cas') == cas:
                # Exact CAS match - use this
                break
            elif hit:
                # Partial match (possibly due to normalization collision) - keep searching
                hit = None
        
        if hit:
            role = hit.get('role', 'other_reagent')
            # Use abbreviation if available, otherwise use name, then CAS
            abbr_list = hit.get('abbreviation', [])
            if abbr_list and abbr_list[0]:
                display_name = abbr_list[0]
            elif hit.get('name'):
                display_name = hit.get('name')
            else:
                display_name = cas
            
            if role in role_items:
                role_items[role].append(display_name)
            else:
                role_items['other_reagent'].append(display_name)
        else:
            # Not found - use CAS as-is
            role_items['other_reagent'].append(cas)
    
    return {
        "catalyst": "/".join(role_items["metal_catalyst"]),
        "ligand": "/".join(role_items["ligand"]),
        "base": "/".join(role_items["base"]),
        "acid": "/".join(role_items["acid"]),
        "oxidant": "/".join(role_items["oxidant"]),
        "reductant": "/".join(role_items["reductant"]),
        "additive": "/".join(role_items["additive"]),
        "condensation_agent": "/".join(role_items["condensation_agent"]),
        "other_reagent": "/".join(role_items["other_reagent"]),
        "solvent": "/".join(role_items["solvent"]),
    }

def enrich_reaction_dataset_cas(input_path: str | Path) -> None:
    """
    Enrich the source reaction dataset CSV with CAS numbers from names
    using the new reagent system. Updates the input CSV file in place.
    """
    input_path = Path(input_path)
    if not input_path.exists():
        print(f"Error: {input_path} not found.")
        return

    print(f"Enriching CAS numbers in {input_path}...")
    try:
        with open(input_path, "r", encoding="utf-8", errors="replace", newline="") as handle:
            df = pd.read_csv(handle)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # Mappings from name column to CAS column and reagent type
    mappings = {
        'reagent_amd': ('reagent_cas', None),
        'catalyst_amd': ('catalyst_cas', 'metal_catalyst'),
        'solvent_amd': ('solvent_cas', 'solvent'),
    }

    modified = False

    lookup_cache: Dict[Tuple[str, Optional[str]], List[str]] = {}

    def lookup_cas(name_str: Any, preferred_type: Optional[str] = None) -> List[str]:
        if not name_str or pd.isna(name_str):
            return []
        cache_key = (str(name_str).strip(), preferred_type)
        cached = lookup_cache.get(cache_key)
        if cached is not None:
            return list(cached)

        # Split by typical separators
        names = _split_reagent_names(str(name_str))
        found_cas = []
        
        for name in names:
            # 1. Try preferred type if provided
            hit = None
            if preferred_type:
                hit = find_reagent(name, preferred_type)
            
            # 2. Try common roles if not found or no preferred type
            if not hit:
                for r_type in ['metal_catalyst', 'ligand', 'base', 'solvent', 'additive', 'acid', 'oxidant', 'reductant']:
                    hit = find_reagent(name, r_type)
                    if hit:
                        break
            
            if hit and hit.get('cas'):
                found_cas.append(hit['cas'])

        lookup_cache[cache_key] = found_cas
        return list(found_cas)

    for name_col, (cas_col, r_type) in mappings.items():
        if name_col not in df.columns:
            continue
        
        # If CAS column doesn't exist, create it
        if cas_col not in df.columns:
            df[cas_col] = ""

        def process_row(row):
            nonlocal modified
            names = row[name_col]
            existing_cas_str = str(row[cas_col]) if not pd.isna(row[cas_col]) else ""
            
            # Extract existing CAS
            existing_cas = [c.strip() for c in existing_cas_str.split(',') if c.strip()]
            
            found_cas = lookup_cas(names, r_type)
            if not found_cas:
                return existing_cas_str

            # Merge and deduplicate
            all_cas = _dedupe_list(existing_cas + found_cas)
            new_cas_str = ", ".join(all_cas)
            
            if new_cas_str != existing_cas_str:
                modified = True
                return new_cas_str
            return existing_cas_str

        df[cas_col] = df.apply(process_row, axis=1)

    if modified:
        df.to_csv(input_path, index=False)
        print(f"Successfully updated CAS numbers in {input_path}")
    else:
        print(f"No changes needed for {input_path}")

def process_reaction_dataset(
    input_path: str,
    output_path: str,
    drop_no_catalyst: bool = True,
    reagent_csv_path: Optional[str | Path] = None,
    new_reagents_path: Optional[str | Path] = None,
    llm_assist_options: Optional[Dict[str, Any]] = None,
):
    """Convert reaction dataset to a minimal HTE recommender CSV."""
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    source_label = input_path.stem

    reagent_csv = Path(reagent_csv_path) if reagent_csv_path else (
        PROJECT_ROOT / "data" / "reagent_db" / "reagents.csv"
    )
    # The new system uses a centralized registry, so we don't need to load it here
    # but we keep the variable for compatibility with the record collection
    unknown_cas: set[str] = set()
    new_reagents_csv = Path(new_reagents_path) if new_reagents_path else (
        PROJECT_ROOT / "data" / "reagent_db" / "new_reagents.csv"
    )
    llm_signature = _llm_assist_signature(llm_assist_options)
    if llm_signature:
        print("LLM assist enabled for reaction featurization.")
    
    rows = []
    taxonomy_gap_entries: List[Dict[str, Any]] = []
    print(f"Reading {input_path}...")
    if not input_path.exists():
        print(f"Error: Input file {input_path} not found.")
        return

    def _csv_row_to_record(row: Dict[str, Any]) -> Dict[str, Any]:
        def _entries(cas_value: Any, amd_value: Any) -> List[Dict[str, str]]:
            entries: List[Dict[str, str]] = []
            for cas in _split_csv_list(cas_value):
                entries.append({"name": cas, "cas": cas})
            for name in _split_csv_list(amd_value):
                entries.append({"name": name})
            return entries

        catalyst_amd = _split_csv_list(row.get("catalyst_amd", ""))
        catalyst_cas = _split_csv_list(row.get("catalyst_cas", ""))
        catalyst_parts = catalyst_amd or catalyst_cas
        catalytic_system = ", ".join(catalyst_parts)

        return {
            "smiles": row.get("reaction_smiles") or row.get("smiles") or "",
            "yield": row.get("yield_pct") or row.get("yield") or "",
            "reagents": _entries(row.get("reagent_cas", ""), row.get("reagent_amd", "")),
            "solvents": _entries(row.get("solvent_cas", ""), row.get("solvent_amd", "")),
            "catalytic_system": catalytic_system,
            "condition_core": "",
        }

    input_suffix = input_path.suffix.lower()
    if input_suffix == ".csv":
        with open(input_path, "r", encoding="utf-8", errors="replace", newline="") as f:
            total = max(sum(1 for _ in f) - 1, 0)
        print(f"Processing {total} reactions...")
        with open(input_path, "r", encoding="utf-8", errors="replace", newline="") as f:
            reader = csv.DictReader(f)
            processed_count = 0
            skipped_no_smiles = 0
            skipped_no_arrow = 0
            skipped_no_catalyst = 0
            skipped_other = 0
            
            for i, row in enumerate(reader):
                if total and i % 500 == 0:
                    print(f"Progress: {i}/{total} ({(i/total)*100:.1f}%)")
                record = _csv_row_to_record(row)
                if not record.get("smiles"):
                    skipped_no_smiles += 1
                    continue
                smiles = record.get("smiles", "")
                if ">>" not in smiles:
                    skipped_no_arrow += 1
                    continue
                reactants_part, _ = smiles.split(">>")
                reactants = reactants_part.split(".")
                reagents = extract_reagents(record, csv_row=row)
                if drop_no_catalyst and not reagents.get("catalyst"):
                    skipped_no_catalyst += 1
                    continue
                _collect_reagent_smiles(record, None, unknown_cas)

                reactant_data = []

                for r_smiles in reactants:
                    try:
                        analysis = cached_featurize(r_smiles)
                        motifs = analysis.get("motifs", [])
                        current_r_motifs = _extract_motif_ids(motifs)
                        context_scaffolds = _extract_context_scaffolds(analysis)

                        current_r_motifs = _apply_aromn_scaffold_override(
                            current_r_motifs,
                            context_scaffolds,
                        )
                        reactant_data.append({
                            "motifs": _dedupe_list(current_r_motifs),
                        })
                    except Exception:
                        continue

                if not reactant_data:
                    continue

                rxn_bundle = cached_featurize_reaction(smiles, llm_signature)
                aggregates = rxn_bundle.get("aggregates") or {}
                reacted_set = set(aggregates.get("reacted_motifs") or [])
                formed_set = set(aggregates.get("formed_motifs") or [])
                spectators_set = set(aggregates.get("spectator_motifs") or [])
                reaction_key = rxn_bundle.get("reaction_key") or ""

                formed_motifs_str = _reactant_key(list(formed_set)) or "None"

                primary_reactant_motifs = _select_primary_reactant_motifs(
                    reactant_data,
                    reacted_set,
                )

                if len(reactants) == 1:
                    primary_reactant_motifs = primary_reactant_motifs[:1]

                type_a = primary_reactant_motifs[0] if len(primary_reactant_motifs) > 0 else ""
                type_b = primary_reactant_motifs[1] if len(primary_reactant_motifs) > 1 else ""
                type_c = primary_reactant_motifs[2] if len(primary_reactant_motifs) > 2 else ""

                spectator_groups = rank_spectator_groups(
                    _collect_spectator_groups(reactant_data, spectators_set)
                )
                detected_reaction_type = _extract_reaction_type_from_bundle(rxn_bundle)
                if not detected_reaction_type:
                    # Fallback keeps compatibility for any unexpected bundle shape.
                    detected_reaction_type = _detect_reaction_type(smiles, llm_signature)
                gap_entry = _build_taxonomy_gap_entry(
                    rxn_bundle,
                    source_dataset=source_label,
                    source_format="csv",
                    source_row_index=i + 2,
                    reaction_smiles=smiles,
                    reaction_key=reaction_key,
                    detected_reaction_type=detected_reaction_type,
                    reference=str(row.get("reference", "") or ""),
                )
                if gap_entry:
                    taxonomy_gap_entries.append(gap_entry)

                row_out = {
                    "reaction_id": source_label,
                    "detected_reaction_type": detected_reaction_type,
                    "yield": record.get("yield", 0.0),
                    "z_score": 0.0,
                    "reactant_1": type_a,
                    "reactant_2": type_b,
                    "reactant_3": type_c,
                    "catalyst": reagents.get("catalyst", ""),
                    "ligand": reagents.get("ligand", ""),
                    "base": reagents.get("base", ""),
                    "acid": reagents.get("acid", ""),
                    "oxidant": reagents.get("oxidant", ""),
                    "reductant": reagents.get("reductant", ""),
                    "additive": reagents.get("additive", ""),
                    "condensation_agent": reagents.get("condensation_agent", ""),
                    "other_reagent": reagents.get("other_reagent", ""),
                    "solvent": reagents.get("solvent", ""),
                    "reaction_smiles": smiles,
                    "formed_motifs": formed_motifs_str,
                    "spectator_groups": " / ".join(spectator_groups),
                    "reference": row.get("reference", ""),
                    "Reaction_Key": reaction_key,
                    "_reaction_key": reaction_key,
                }
                rows.append(row_out)
                processed_count += 1
            
            print(f"Processing summary:")
            print(f"  Processed: {processed_count}")
            print(f"  Skipped - no SMILES: {skipped_no_smiles}")
            print(f"  Skipped - no >>: {skipped_no_arrow}")
            print(f"  Skipped - no catalyst: {skipped_no_catalyst}")
    else:
        with open(input_path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        total = len(lines)
        print(f"Processing {total} reactions...")
        
        processed_count = 0
        skipped_no_smiles = 0
        skipped_no_arrow = 0
        skipped_no_catalyst = 0

        for i, line in enumerate(lines):
            if i % 500 == 0:
                print(f"Progress: {i}/{total} ({(i/total)*100:.1f}%)")

            try:
                record = json.loads(line)
            except:
                continue
            
            smiles = record.get("smiles", "")
            if ">>" not in smiles:
                skipped_no_arrow += 1
                continue

            reactants_part, _ = smiles.split(">>")
            reactants = reactants_part.split(".")
            reagents = extract_reagents(record)
            if drop_no_catalyst and not reagents.get("catalyst"):
                skipped_no_catalyst += 1
                continue
            _collect_reagent_smiles(record, None, unknown_cas)

            reactant_data = []

            for r_smiles in reactants:
                try:
                    analysis = cached_featurize(r_smiles)
                    motifs = analysis.get("motifs", [])
                    current_r_motifs = _extract_motif_ids(motifs)
                    context_scaffolds = _extract_context_scaffolds(analysis)

                    current_r_motifs = _apply_aromn_scaffold_override(
                        current_r_motifs,
                        context_scaffolds,
                    )
                    reactant_data.append({
                        "motifs": _dedupe_list(current_r_motifs),
                    })
                except Exception:
                    continue

            if not reactant_data:
                continue

            rxn_bundle = cached_featurize_reaction(smiles, llm_signature)
            aggregates = rxn_bundle.get("aggregates") or {}
            reacted_set = set(aggregates.get("reacted_motifs") or [])
            formed_set = set(aggregates.get("formed_motifs") or [])
            spectators_set = set(aggregates.get("spectator_motifs") or [])
            reaction_key = rxn_bundle.get("reaction_key") or ""

            formed_motifs_str = _reactant_key(list(formed_set)) or "None"

            primary_reactant_motifs = _select_primary_reactant_motifs(
                reactant_data,
                reacted_set,
            )

            if len(reactants) == 1:
                primary_reactant_motifs = primary_reactant_motifs[:1]

            type_a = primary_reactant_motifs[0] if len(primary_reactant_motifs) > 0 else ""
            type_b = primary_reactant_motifs[1] if len(primary_reactant_motifs) > 1 else ""
            type_c = primary_reactant_motifs[2] if len(primary_reactant_motifs) > 2 else ""

            spectator_groups = rank_spectator_groups(
                _collect_spectator_groups(reactant_data, spectators_set)
            )
            detected_reaction_type = _extract_reaction_type_from_bundle(rxn_bundle)
            if not detected_reaction_type:
                # Fallback keeps compatibility for any unexpected bundle shape.
                detected_reaction_type = _detect_reaction_type(smiles, llm_signature)
            gap_entry = _build_taxonomy_gap_entry(
                rxn_bundle,
                source_dataset=source_label,
                source_format="jsonl",
                source_row_index=i + 1,
                reaction_smiles=smiles,
                reaction_key=reaction_key,
                detected_reaction_type=detected_reaction_type,
                reference=str(record.get("reference", "") or ""),
            )
            if gap_entry:
                taxonomy_gap_entries.append(gap_entry)

            row = {
                "reaction_id": source_label,
                "detected_reaction_type": detected_reaction_type,
                "yield": record.get("yield", 0.0),
                "z_score": 0.0,
                "reactant_1": type_a,
                "reactant_2": type_b,
                "reactant_3": type_c,
                "catalyst": reagents.get("catalyst", ""),
                "ligand": reagents.get("ligand", ""),
                "base": reagents.get("base", ""),
                "acid": reagents.get("acid", ""),
                "oxidant": reagents.get("oxidant", ""),
                "reductant": reagents.get("reductant", ""),
                "additive": reagents.get("additive", ""),
                "condensation_agent": reagents.get("condensation_agent", ""),
                "other_reagent": reagents.get("other_reagent", ""),
                "solvent": reagents.get("solvent", ""),
                "reaction_smiles": smiles,
                "formed_motifs": formed_motifs_str,
                "spectator_groups": " / ".join(spectator_groups),
                "reference": record.get("reference", ""),
                "Reaction_Key": reaction_key,
                "_reaction_key": reaction_key,
            }
            rows.append(row)
            processed_count += 1
        
        print(f"Processing summary:")
        print(f"  Processed: {processed_count}")
        print(f"  Skipped - no SMILES: {skipped_no_smiles}")
        print(f"  Skipped - no >>: {skipped_no_arrow}")
        print(f"  Skipped - no catalyst: {skipped_no_catalyst}")
        
    if not rows:
        print("Warning: No valid reactions were processed.")
        if skipped_no_catalyst > 0:
            print(f"  {skipped_no_catalyst} reactions were skipped due to missing catalyst.")
            print(f"  Consider running without --drop-no-catalyst or add catalyst information to the source data.")
        return
    
    df = pd.DataFrame(rows)
    
    # Ensure yield is numeric
    df["yield"] = pd.to_numeric(df["yield"], errors='coerce').fillna(0.0)
    
    # Calculate Pseudo Z-Score based on motif groups
    print("Calculating Pseudo Z-Scores...")
    df["z_score"] = 0.0
    df["z_score"] = df["z_score"].astype(float)
    
    if "_reaction_key" in df.columns:
        group_key = "_reaction_key"
    else:
        group_key = "reaction_id"

    for key, group in df.groupby(group_key):
        if len(group) > 1:
            yields = group["yield"]
            mean = yields.mean()
            std = yields.std()
            if std > 0:
                df.loc[group.index, "z_score"] = (yields - mean) / std
            else:
                # If all yields are same, Z-score is 0
                df.loc[group.index, "z_score"] = 0.0
        else:
            # Single experiment in group
            df.loc[group.index, "z_score"] = 0.0

    df["z_score"] = df["z_score"].round(2)

    canonical_cols = [
        "reaction_id", "detected_reaction_type", "reaction_smiles", "Reaction_Key",
        "reactant_1", "reactant_2", "reactant_3", "formed_motifs",
        "catalyst", "ligand", "base", "acid", "oxidant", "reductant",
        "additive", "condensation_agent", "other_reagent", "solvent",
        "spectator_groups", "reference", "yield", "z_score",
    ]
    df = df[canonical_cols]

    df.to_csv(output_path, index=False)
    print(f"Successfully saved {len(df)} reactions to {output_path}")
    if llm_signature:
        if taxonomy_gap_entries:
            _write_taxonomy_gap_exports(output_path, taxonomy_gap_entries)
        else:
            print("No taxonomy-gap entries were flagged by LLM assist.")
    if unknown_cas:
        _write_new_reagents(new_reagents_csv, unknown_cas)
        print(f"Saved {len(unknown_cas)} new reagent CAS entries to {new_reagents_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert reaction dataset to HTE-canonical CSV format.")
    parser.add_argument("--dataset", "-d", help="Name of the dataset (e.g., C_N_Coupling). If provided, processes 'data/reaction_dataset/{dataset}.jsonl' to 'data/HTE_db/{dataset}_canonical.csv'.")
    parser.add_argument("--input", "-i", help="Direct path to input JSONL file.")
    parser.add_argument("--output", "-o", help="Direct path to output CSV file.")
    parser.add_argument(
        "--reagent-csv",
        help="Path to reagent registry CSV (default: data/reagent_db/reagents.csv).",
    )
    parser.add_argument(
        "--new-reagents",
        help="Path to write missing reagent CAS list (default: data/reagent_db/new_reagents.csv).",
    )
    parser.add_argument(
        "--keep-no-catalyst",
        action="store_true",
        help="Keep reactions without a catalyst (by default, they are dropped).",
    )
    
    args = parser.parse_args()
    
    if args.dataset:
        input_file = f"data/reaction_dataset/{args.dataset}.jsonl"
        output_file = f"data/HTE_db/{args.dataset}_canonical.csv"
        process_reaction_dataset(
            input_file,
            output_file,
            drop_no_catalyst=not args.keep_no_catalyst,
            reagent_csv_path=args.reagent_csv,
            new_reagents_path=args.new_reagents,
        )
    elif args.input and args.output:
        process_reaction_dataset(
            args.input,
            args.output,
            drop_no_catalyst=not args.keep_no_catalyst,
            reagent_csv_path=args.reagent_csv,
            new_reagents_path=args.new_reagents,
        )
    else:
        # Default: Process known core datasets
        datasets = ["C_N_Coupling", "C_O_Coupling", "C_S_Coupling"]
        for ds in datasets:
            input_file = f"data/reaction_dataset/{ds}.jsonl"
            output_file = f"data/HTE_db/{ds}_canonical.csv"
            process_reaction_dataset(
                input_file,
                output_file,
                reagent_csv_path=args.reagent_csv,
                new_reagents_path=args.new_reagents,
            )
