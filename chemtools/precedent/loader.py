"""Dataset loading and transformation for precedent search.

This module handles loading reaction datasets from JSONL files, normalizing
family names, and transforming raw dataset records into the precedent format.
"""
from typing import Dict, Any, List, Optional, Tuple
import os
import json
from functools import lru_cache

# Helper to pick electrophile vs nucleophile from reactants list
def _pick_electrophile_nucleophile(reactants: List[str]) -> Tuple[str, str]:
    """Heuristically identify electrophile and nucleophile from reactant list."""
    def is_electrophile(s: str) -> bool:
        t = (s or "").lower()
        return (
            ("br" in t) or ("cl" in t) or (" i" in t)
            or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
        )
    if not reactants:
        return "", ""
    if len(reactants) == 1:
        return reactants[0], ""
    r0, r1 = reactants[0], reactants[1]
    if is_electrophile(r0):
        return r0, r1
    if is_electrophile(r1):
        return r1, r0
    return r0, r1


# Path to small demo dataset used for precedents retrieval
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
DATA_PATH = os.path.join(BASE_DIR, "data", "reactions_sample.jsonl")
# Allow overriding the dataset directory so data-processor can save directly there
_ENV_DIR = os.environ.get("CHEMTOOLS_DATASET_DIR", "").strip()
DATASET_DIR = (
    os.path.abspath(_ENV_DIR) if _ENV_DIR else os.path.join(BASE_DIR, "data", "reaction_dataset")
)


def _iter_dataset_files() -> List[str]:
    """List all JSONL files in the dataset directory."""
    files: List[str] = []
    if os.path.isdir(DATASET_DIR):
        for name in os.listdir(DATASET_DIR):
            if name.lower().endswith(".jsonl"):
                files.append(os.path.join(DATASET_DIR, name))
    return sorted(files)


def _dataset_family_map(raw: Optional[str], fallback: Optional[str] = None) -> str:
    """Normalize dataset reaction_type to API family text.
    
    Maps both legacy and new naming conventions to canonical family names.
    """
    t = (raw or "").strip()
    if not t and fallback:
        t = fallback.strip()
    if not t:
        return ""
    tl = t.lower()
    
    # New systematic naming (preferred)
    if tl in {
        "c_n_cross_coupling",
        "c_n_coupling",
        "c_n_coupling_cu",
        "c_n_coupling_cu_ullmann",
        "c_n_coupling_pd",
        "c_n_coupling_pd_buchwald",
        "c_n_coupling_ni",
    }:
        return "C_N_Coupling"
    if tl in {"snar_cn", "snar_cn_coupling", "snar c-n"}:
        return "SNAr-CN"
    if tl in {"snar_co", "snar_co_coupling", "snar c-o"}:
        return "SNAr-CO"
    if tl in {"snar_cs", "snar_cs_coupling", "snar c-s"}:
        return "SNAr-CS"
    
    # Legacy naming (supported for backward compatibility)
    if tl in {
        "ullman",
        "ullmann",
        "ullman-c-n",
        "ullmann-c-n",
        "ullmann c-n",
        "ullmann_cn",
    }:
        return "C_N_Coupling"
    if tl in {"buchwald", "buchwald-c-n", "buchwald c-n", "buchwald_cn"}:
        return "C_N_Coupling"
    
    # Other reaction types - use exact dataset file names
    if tl in {
        "suzuki",
        "suzuki-miyaura",
        "suzuki_miyaura",
        "suzuki cc",
        "suzuki_cc",
        "suzuki_coupling",
    }:
        return "Suzuki"
    if tl in {"sonogashira", "sonogashira_coupling", "sonogashira cc", "sonogashira_cc"}:
        return "Sonogashira_coupling"
    if tl in {"heck", "heck_mizoroki", "heckmizoroki", "heckmizoroki_coupling"}:
        return "HeckMizoroki_coupling"
    if tl in {"amide-formation", "amide formation", "amideformation", "amide", "amide coupling", "amide_coupling"}:
        return "Amide_formation"
    if tl in {"amide_formation"}:
        return "Amide_formation"
    
    return t


def _make_row_from_dataset(
    rec: Dict[str, Any],
    *,
    file_family: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """Transform a raw dataset record into the standardized precedent format.
    
    Args:
        rec: Raw record from JSONL dataset
        file_family: Optional filename-derived family name for fallback mapping
        
    Returns:
        Standardized precedent row dict, or None if transformation fails
    """
    try:
        # Import featurizers here to avoid circular dependency
        from ..featurizers import reaction_pair as feat_pair
        
        payload = rec.get("recommendation_payload") or {}
        rxn_id = rec.get("reaction_id")
        rt = rec.get("reaction_type") or payload.get("family") or rec.get("reaction_family") or rec.get("family")
        fam_txt = _dataset_family_map(rt, fallback=file_family)
        cond = rec.get("conditions") or {}
        y = cond.get("yield_pct")
        if y is None:
            y = cond.get("yield")
        if y is None:
            y = cond.get("yield_percent")
        T_C = cond.get("temperature_c")
        if T_C is None:
            T_C = cond.get("temp_c")
        if T_C is None:
            T_C = cond.get("T_C")
        time_h = cond.get("time_h")
        if time_h is None:
            time_h = cond.get("time_hours")
        core = rec.get("condition_core") or rec.get("core")
        
        def _coerce_list(value: Any) -> List[str]:
            if not value:
                return []
            if isinstance(value, list):
                out = []
                for item in value:
                    if isinstance(item, (str, int, float)):
                        text = str(item).strip()
                        if text:
                            out.append(text)
                    elif isinstance(item, dict):
                        text = item.get("name") or item.get("cas") or item.get("uid")
                        if text:
                            out.append(str(text).strip())
                return out
            if isinstance(value, str):
                text = value.strip()
                return [text] if text else []
            return []
        
        def _first_value(value: Any) -> Optional[str]:
            items = _coerce_list(value)
            return items[0] if items else None
        
        def _extract_uid(entry: Any) -> Optional[str]:
            if isinstance(entry, dict):
                return entry.get("cas") or entry.get("uid") or entry.get("name")
            if isinstance(entry, str):
                return entry.strip() or None
            return None
        
        # Base/Solvent: take first entries' CAS where present
        base_uid = None
        reagents = rec.get("reagents")
        if not isinstance(reagents, list):
            reagents = []
        elif reagents and isinstance(reagents[0], str):
            reagents = [{"name": name, "role": ""} for name in reagents if str(name).strip()]
        
        if not reagents:
            for key, role in {
                "base": "BASE",
                "additive": "ADDITIVE",
                "acid": "ACID",
                "oxidant": "OXIDANT",
                "reductant": "REDUCTANT",
            }.items():
                for name in _coerce_list(cond.get(key)):
                    reagents.append({"name": name, "role": role})
        
        for rg in reagents or []:
            if (rg.get("role") or "").upper() == "BASE":
                base_uid = _extract_uid(rg)
                if base_uid:
                    break
        if not base_uid:
            base_uid = _first_value(cond.get("base"))
        solvent_uid = None
        solvents = rec.get("solvents")
        if not isinstance(solvents, list):
            solvents = []
        elif solvents and isinstance(solvents[0], str):
            solvents = [{"name": name} for name in solvents if str(name).strip()]
        
        if not solvents:
            for name in _coerce_list(cond.get("solvent") or cond.get("solvents")):
                solvents.append({"name": name})
        
        if solvents:
            solvent_uid = _extract_uid(solvents[0])
        if not solvent_uid:
            solvent_uid = _first_value(cond.get("solvent") or cond.get("solvents"))
        
        # Reaction SMILES: try precomputed normalized first, fallback to building from raw
        precomputed = rec.get("precomputed") or {}
        rxn_smiles = rec.get("reaction_smiles") or precomputed.get("reaction_smiles")
        
        if not rxn_smiles:
            # Fallback: build from raw reactants>>products
            smiles_block = rec.get("smiles") or {}
            rcts = (smiles_block.get("reactants") or "").strip()
            prods = (smiles_block.get("products") or "").strip()
            rxn_smiles = f"{rcts}>>{prods}"
        
        # Features: use precomputed if available, otherwise compute
        features = precomputed.get("features") or payload.get("features") or rec.get("features")
        
        if not features:
            # Fallback: compute features on-the-fly (legacy datasets)
            smiles_block = rec.get("smiles") or {}
            rcts = (smiles_block.get("reactants") or "").strip()
            if rcts:
                reactants_list = [p for p in (rcts.split('.') if rcts else []) if p]
                elec, nuc = _pick_electrophile_nucleophile(reactants_list)
                feat_result = feat_pair.featurize_pair(elec, nuc)
                features = feat_result.get("flat", {})
            else:
                features = {}
        
        # Build uniform row
        catalyst_obj = rec.get("catalyst") or {}
        full_system = catalyst_obj.get("full_system") if isinstance(catalyst_obj, dict) else None
        
        # Load catalytic_system from dataset (new field added to datasets)
        catalytic_system = rec.get("catalytic_system")
        # Fallback to full_system if catalytic_system not present (backward compatibility)
        if not catalytic_system and full_system:
            catalytic_system = full_system
        
        return {
            "reaction_id": rxn_id,
            "rxn_type": fam_txt,
            "yield_value": y,
            "T_C": T_C,
            "time_h": time_h,
            "condition_core": core,
            "base_uid": base_uid,
            "solvent_uid": solvent_uid,
            "reagents": reagents,
            "solvents": solvents,
            "reference": rec.get("reference") or {},
            "conditions": cond,
            "catalyst": catalyst_obj,
            "full_system": full_system,
            "catalytic_system": catalytic_system,
            "features": features,
            "reaction_smiles": rxn_smiles,
            "precomputed": precomputed,  # Include precomputed field for DRFP fingerprints
        }
    except Exception:
        return None


def _load_selective(families: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Load dataset rows, optionally filtered by reaction family.
    
    Args:
        families: Optional list of family names to load (e.g., ["C_N_Coupling_Pd"]).
                 If None, loads all datasets. Supports both canonical names and file names.
    
    Returns:
        List of reaction rows matching the specified families.
    """
    from .core_utils import _family_text
    
    rows: List[Dict[str, Any]] = []
    
    # Normalize family filter
    family_filter = None
    if families:
        family_filter = set()
        for f in families:
            # Add both the input and canonical mapping
            family_filter.add(f)
            family_filter.add(_family_text(f))
            # Also add lowercase version
            family_filter.add(f.lower())
    
    # 1) Small demo file, if present
    if os.path.exists(DATA_PATH):
        try:
            with open(DATA_PATH, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rows.append(json.loads(line))
                    except Exception:
                        pass
        except Exception:
            pass
    
    # 2) Auto-load local dataset directory (transformed) when available, unless disabled.
    #    - Explicit override via CHEMTOOLS_LOAD_DATASET (1/true/on to enable; 0/false/off to disable)
    #    - During pytest (PYTEST_CURRENT_TEST set), default to disabled to keep unit tests deterministic
    _flag = str(os.environ.get("CHEMTOOLS_LOAD_DATASET", "")).strip().lower()
    if _flag in {"0", "false", "no", "off"}:
        use_dataset = False
    elif _flag in {"1", "true", "yes", "on"}:
        use_dataset = True
    else:
        use_dataset = ("PYTEST_CURRENT_TEST" not in os.environ)
    
    if not rows and use_dataset and os.path.isdir(DATASET_DIR):
        for path in _iter_dataset_files():
            # If family filter is active, check if this file should be loaded
            if family_filter:
                # Extract family name from filename (e.g., "C_N_Coupling_Pd.jsonl" -> "C_N_Coupling_Pd")
                file_name = os.path.basename(path)
                file_family = file_name.replace('.jsonl', '')
                
                # Check if this family matches the filter
                if not any(ff in family_filter for ff in [file_family, file_family.lower(), _family_text(file_family)]):
                    continue  # Skip this file
            
            try:
                with open(path, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        try:
                            rec = json.loads(line)
                        except Exception:
                            continue
                        row = _make_row_from_dataset(rec, file_family=file_family)
                        if row is not None:
                            # Double-check family filter on row data
                            if family_filter:
                                row_family = row.get("rxn_type") or ""
                                if not any(ff in family_filter for ff in [row_family, row_family.lower()]):
                                    continue
                            rows.append(row)
            except Exception:
                continue
    return rows


@lru_cache(maxsize=1)
def _load() -> List[Dict[str, Any]]:
    """Load all datasets. Cached for performance."""
    return _load_selective(families=None)
