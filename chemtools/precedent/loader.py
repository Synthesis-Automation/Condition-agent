"""Dataset loading and transformation for precedent search.

This module loads HTE source CSVs and transforms records into the
standardized precedent format.
"""
from typing import Dict, Any, List, Optional, Tuple
import csv
import os
from functools import lru_cache

from ..synthon import select_electrophile_nucleophile


def _pick_electrophile_nucleophile(reactants: List[str]) -> Tuple[str, str]:
    """Identify electrophile and nucleophile from reactant list."""
    return select_electrophile_nucleophile(reactants)


BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
_ENV_HTE_DIR = os.environ.get("CHEMTOOLS_HTE_DB_DIR", "").strip()
HTE_DB_DIR = (
    os.path.abspath(_ENV_HTE_DIR)
    if _ENV_HTE_DIR
    else os.path.join(BASE_DIR, "data", "HTE_db")
)
# Backward compatibility
_ENV_LIT_DIR = os.environ.get("CHEMTOOLS_LITERATURE_DIR", "").strip()
if _ENV_LIT_DIR:
    HTE_DB_DIR = os.path.dirname(os.path.abspath(_ENV_LIT_DIR))


def _iter_literature_files() -> List[str]:
    """List HTE CSV files across source folders (protocols are treated as literature)."""
    files: List[str] = []
    subdirs = ["literature", "protocols", "rules", "motif", "experiments"]
    
    for subdir in subdirs:
        subdir_path = os.path.join(HTE_DB_DIR, subdir)
        if os.path.isdir(subdir_path):
            for name in os.listdir(subdir_path):
                if name.lower().endswith(".csv"):
                    files.append(os.path.join(subdir_path, name))
    return sorted(files)


def _infer_source_group_from_path(path: str) -> str:
    text = os.path.normpath(path).lower()
    parts = text.split(os.sep)
    if "literature" in parts or "datasets" in parts or "dataset" in parts:
        return "literature"
    if "protocols" in parts or "protocol" in parts:
        return "literature"
    if "rules" in parts or "rule" in parts:
        return "rules"
    if "motif" in parts or "motifs" in parts or "experiments" in parts or "experiment" in parts or "experiements" in parts:
        return "experiments"
    return "unknown"


def _file_family_from_name(filename: str) -> str:
    stem = os.path.splitext(os.path.basename(filename))[0]
    if stem.lower().endswith("_canonical"):
        stem = stem[: -len("_canonical")]
    return stem


def _clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _parse_float(value: Any) -> Optional[float]:
    text = _clean_text(value)
    if not text:
        return None
    try:
        return float(text)
    except Exception:
        return None


def _split_items(text: str) -> List[str]:
    if not text:
        return []
    parts = [p.strip() for p in text.replace(";", "/").split("/") if p.strip()]
    return parts if parts else [text]


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
    if tl in {"c-n coupling", "c_n coupling", "c n coupling", "c-n-coupling"}:
        return "C_N_Coupling"
    if tl in {"c-o coupling", "c_o coupling", "c o coupling", "c-o-coupling"}:
        return "C_O_Coupling"
    if tl in {"c-s coupling", "c_s coupling", "c s coupling", "c-s-coupling"}:
        return "C_S_Coupling"
    if tl in {"snar_cn", "snar_cn_coupling", "snar c-n"}:
        return "C_N_Coupling"
    if tl in {"snar_co", "snar_co_coupling", "snar c-o"}:
        return "C_O_Coupling"
    if tl in {"snar_cs", "snar_cs_coupling", "snar c-s"}:
        return "C_S_Coupling"
    
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


def _make_row_from_csv(
    rec: Dict[str, Any],
    *,
    row_index: int,
    file_family: Optional[str] = None,
    source_group: Optional[str] = None,
    fast: bool = False,
) -> Optional[Dict[str, Any]]:
    """Transform a CSV record into the standardized precedent format.

    Args:
        fast: When True, skip ``normalize_reaction`` and ``featurize_pair``.
              Produces the same dict but with ``features={}`` and uses the raw
              reaction SMILES without normalization.  ~10-20x faster per row;
              sufficient for index-building pipelines that only need
              ``reaction_smiles``, ``reaction_id``, ``condition_core``, and
              ``rxn_type``.
    """
    try:
        if not fast:
            from ..featurizers import reaction_pair as feat_pair
            from ..smiles import normalize_reaction

        raw_reaction_id = _clean_text(rec.get("reaction_id"))
        fam_txt = _dataset_family_map(raw_reaction_id, fallback=file_family)
        row_id = _clean_text(file_family) or fam_txt or "reaction"
        rxn_id = f"{row_id}:{row_index}"

        yield_val = _parse_float(
            rec.get("yield") or rec.get("yield_pct") or rec.get("yield_percent")
        )
        z_score = _parse_float(rec.get("z_score"))

        catalyst = _clean_text(rec.get("catalyst"))
        ligand = _clean_text(rec.get("ligand"))
        base = _clean_text(rec.get("base"))
        acid = _clean_text(rec.get("acid"))
        oxidant = _clean_text(rec.get("oxidant"))
        reductant = _clean_text(rec.get("reductant"))
        additive = _clean_text(rec.get("additive"))
        condensation_agent = _clean_text(rec.get("condensation_agent"))
        other_reagent = _clean_text(rec.get("other_reagent"))
        solvent = _clean_text(rec.get("solvent"))

        core_parts = [p for p in (catalyst, ligand) if p]
        core = "/".join(core_parts)

        reagents = []
        for name, role in (
            (base, "BASE"),
            (acid, "ACID"),
            (oxidant, "OXIDANT"),
            (reductant, "REDUCTANT"),
            (additive, "ADDITIVE"),
            (condensation_agent, "COUPLING_REAGENT"),
            (other_reagent, "OTHER"),
        ):
            if name:
                reagents.append({"name": name, "role": role})

        solvents = []
        for name in _split_items(solvent):
            if name:
                solvents.append({"name": name})

        catalytic_system = []
        if catalyst:
            catalytic_system.append({"name": catalyst, "role": "CATALYST"})
        if ligand:
            catalytic_system.append({"name": ligand, "role": "LIGAND"})

        rxn_smiles = _clean_text(rec.get("reaction_smiles"))
        if not fast:
            normalized = normalize_reaction(rxn_smiles) if rxn_smiles else None
            if normalized:
                rxn_smiles = normalized.get("normalized") or rxn_smiles
        else:
            normalized = None

        features: Dict[str, Any] = {}
        reactant_smiles: List[str] = []
        if normalized:
            for entry in normalized.get("reactants", []) or []:
                if not isinstance(entry, dict):
                    continue
                smi = entry.get("smiles_norm") or entry.get("largest_smiles") or entry.get("input")
                if smi:
                    reactant_smiles.append(smi)
        if not reactant_smiles:
            reactant_smiles = [
                _clean_text(rec.get("reactant_1")),
                _clean_text(rec.get("reactant_2")),
                _clean_text(rec.get("reactant_3")),
            ]
            reactant_smiles = [s for s in reactant_smiles if s]

        if reactant_smiles and not fast:
            elec, nuc = _pick_electrophile_nucleophile(reactant_smiles)
            feat_result = feat_pair.featurize_pair(elec, nuc)
            features = feat_result.get("flat", {}) if isinstance(feat_result, dict) else {}

        conditions = {
            "catalyst": catalyst,
            "ligand": ligand,
            "base": base,
            "acid": acid,
            "oxidant": oxidant,
            "reductant": reductant,
            "additive": additive,
            "condensation_agent": condensation_agent,
            "other_reagent": other_reagent,
            "solvent": solvent,
            "yield_pct": yield_val,
            "z_score": z_score,
        }

        catalyst_obj = {"name": catalyst} if catalyst else {}

        return {
            "reaction_id": rxn_id,
            "dataset_reaction_id": raw_reaction_id,
            "source_file": _clean_text(file_family),
            "source_group": _clean_text(source_group),
            "rxn_type": fam_txt,
            "yield_value": yield_val,
            "T_C": None,
            "time_h": None,
            "condition_core": core,
            "base_uid": base or None,
            "solvent_uid": solvents[0]["name"] if solvents else (solvent or None),
            "reagents": reagents,
            "solvents": solvents,
            "reference": _clean_text(rec.get("reference")),
            "conditions": conditions,
            "catalyst": catalyst_obj,
            "full_system": None,
            "catalytic_system": catalytic_system,
            "features": features,
            "reaction_smiles": rxn_smiles,
            "precomputed": {},
        }
    except Exception:
        return None


def _read_csv_records(path: str) -> List[Dict[str, Any]]:
    encodings = ("utf-8", "utf-8-sig", "cp1252", "latin-1")
    last_exc: Optional[Exception] = None
    for encoding in encodings:
        try:
            with open(path, "r", encoding=encoding, newline="") as handle:
                reader = csv.DictReader(handle)
                return [dict(row) for row in reader]
        except Exception as exc:
            last_exc = exc
    if last_exc:
        raise last_exc
    return []


def _family_key(family_filter: Optional[set]) -> Tuple[str, ...]:
    if not family_filter:
        return tuple()
    return tuple(sorted(str(item) for item in family_filter))


@lru_cache(maxsize=8)
def _load_literature_cached(family_key: Tuple[str, ...]) -> List[Dict[str, Any]]:
    family_filter = set(family_key) if family_key else None
    family_lower = {f.lower() for f in family_filter} if family_filter else set()
    rows: List[Dict[str, Any]] = []

    for path in _iter_literature_files():
        file_family = _file_family_from_name(path)
        source_group = _infer_source_group_from_path(path)
        mapped_family = _dataset_family_map(file_family, fallback=file_family)
        if family_filter:
            candidates = {
                file_family,
                mapped_family,
                file_family.lower(),
                mapped_family.lower(),
            }
            if not any(c in family_filter or c in family_lower for c in candidates):
                continue

        try:
            records = _read_csv_records(path)
        except Exception:
            continue

        for row_index, rec in enumerate(records):
            row = _make_row_from_csv(
                rec,
                row_index=row_index,
                file_family=file_family,
                source_group=source_group,
            )
            if row is None:
                continue
            if family_filter:
                row_family = (row.get("rxn_type") or "")
                if row_family and row_family not in family_filter and row_family.lower() not in family_lower:
                    continue
            rows.append(row)
    return rows

def _load_selective(families: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Load HTE database CSV rows from all sources (protocols folded into literature),
    optionally filtered by reaction family.
    
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
    
    csv_files = _iter_literature_files()
    if csv_files:
        rows = list(_load_literature_cached(_family_key(family_filter)))
    return rows


@lru_cache(maxsize=1)
def _load() -> List[Dict[str, Any]]:
    """Load all datasets. Cached for performance."""
    return _load_selective(families=None)


def _load_all_with_progress() -> List[Dict[str, Any]]:
    """Like ``_load()`` but prints per-file tqdm progress.

    Used by build scripts where a long silent pause is unacceptable.
    Result is NOT cached (use ``_load()`` for cached access after build).
    """
    from .core_utils import _family_text  # noqa: F401 (kept for consistency)

    files = _iter_literature_files()
    rows: List[Dict[str, Any]] = []

    try:
        from tqdm import tqdm
        file_iter = tqdm(files, desc="  Loading CSVs", unit="file", dynamic_ncols=True)
    except ImportError:
        import logging as _logging
        _log = _logging.getLogger(__name__)
        _total = len(files)
        _report = max(1, _total // 10)

        def _plain_iter(fs):
            for _i, _f in enumerate(fs):
                if _i % _report == 0:
                    _log.info("  Loading file %d / %d …", _i + 1, _total)
                yield _f

        file_iter = _plain_iter(files)

    for path in file_iter:
        file_family = _file_family_from_name(path)
        source_group = _infer_source_group_from_path(path)
        try:
            records = _read_csv_records(path)
        except Exception:
            continue
        for row_index, rec in enumerate(records):
            row = _make_row_from_csv(
                rec,
                row_index=row_index,
                file_family=file_family,
                source_group=source_group,
                fast=True,  # skip normalize_reaction + featurize_pair
            )
            if row is not None:
                rows.append(row)
    return rows
