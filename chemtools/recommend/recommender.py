"""
HTE-based Condition Recommendation System

This module provides condition recommendations based on High-Throughput Experimentation (HTE) data.
Recommendations are primarily based on reactant types since no reaction SMILES is provided.

Key Features:
- Reactant type-based matching using existing chemtools detection
- Z-score based ranking (primary metric for condition success)
- Success-weighted condition selection (yield > 50% threshold)
- Statistical ranking with confidence scores
- Multi-reaction type support with automatic detection
- Condition component recommendations (catalyst, ligand, base, solvent)
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any, Iterable, Set, Iterator, Mapping
from collections import defaultdict, Counter, OrderedDict
from functools import lru_cache
import hashlib
import pickle
import itertools
import re
import logging
import time
import pandas as pd
import numpy as np
from pathlib import Path
import json

from chemtools.featurizers.reaction_path import (
    build_reaction_featurization_options,
    cleanup_reaction_smiles_for_featurization,
)
from chemtools.featurizers.unified import (
    featurize_molecule,
    featurize_reaction as _unified_featurize_reaction,
)
from chemtools.featurizers.analysis.reaction_record import ReactionRecord
from chemtools.featurizers.formatters.utils import normalize_motif_id
from chemtools.featurizers.spectator_rank import (
    spectator_group_weight,
    weighted_spectator_similarity,
)
from chemtools.smiles import normalize_reaction
from chemtools.recommend.reaction_key_utils import (
    build_reaction_events_payload,
    canonicalize_reaction_key_minimal,
    deserialize_reaction_events_text,
    normalize_reaction_events_text,
    serialize_reaction_events_payload,
)
from chemtools.taxonomy import loader as taxonomy_loader
from chemtools.taxonomy.substituent_composer import load_organic_groups_with_compositions

try:
    from chemtools.taxonomy import reaction_catalog as _reaction_catalog
except Exception:
    _reaction_catalog = None

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOGGER = logging.getLogger(__name__)
_EVENT_MATCH_PREFIX = "RXNEVT|"
_EVENT_SIGNATURE_PRIORITY = [
    "benzyl_o_alkylation_like",
    "ester_hydrolysis_like",
    "amidation_like",
    "ring_closure_or_annulation",
    "leaving_group_displacement",
    "c_n_bond_formation",
    "c_o_bond_formation",
    "c_s_bond_formation",
    "c_c_bond_formation",
]
_EVENT_SIGNATURE_CODE = {
    "benzyl_o_alkylation_like": "BzOAlk",
    "ester_hydrolysis_like": "EsterHyd",
    "amidation_like": "Amid",
    "ring_closure_or_annulation": "Ann",
    "leaving_group_displacement": "LGDisp",
    "c_n_bond_formation": "C-N",
    "c_o_bond_formation": "C-O",
    "c_s_bond_formation": "C-S",
    "c_c_bond_formation": "C-C",
}
_PRECEDENT_REACTION_EVENT_BLEND_WEIGHT = 0.7


def featurize_reaction(
    reaction_smiles: str,
    options: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    cleaned_smiles, _ = cleanup_reaction_smiles_for_featurization(reaction_smiles)
    return _unified_featurize_reaction(
        cleaned_smiles,
        options=build_reaction_featurization_options(base_options=options),
    )
_REACTION_EVENT_COMPONENT_WEIGHTS = {
    "formed": 0.35,
    "broken": 0.20,
    "signature": 0.10,
    "event_kinds": 0.05,
    "redox": 0.05,
    "reacted_context": 0.20,
    "formed_context": 0.05,
}


def _infer_source_group(source_path: Optional[Path]) -> str:
    if not source_path:
        return "unknown"
    parts = [part.lower() for part in source_path.parts]
    for part in parts:
        if part in ("literature", "datasets", "dataset"):
            return "literature"
        if "protocol" in part:  # Protocol datasets are folded into literature mode.
            return "literature"
        if "rule" in part:      # Matches rules, rule_db
            return "rules"
        if part in ("motif", "motifs", "experiments", "experiment", "experiements"):
            return "motif"
    return "other"


def _format_source_path(source_path: Optional[Path]) -> str:
    if not source_path:
        return ""
    try:
        return source_path.resolve().relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(source_path)


def _collect_hte_files(db_path: Path) -> List[Path]:
    if db_path.is_file():
        return [db_path]
    if not db_path.exists():
        return []

    candidates: List[Path] = []
    candidates.extend(db_path.glob("*.csv"))

    for subdir in ("literature", "datasets", "protocols", "rules", "motif", "experiments", "experiment", "experiements"):
        sub_path = db_path / subdir
        if not sub_path.exists():
            continue
        candidates.extend(sub_path.glob("*.csv"))

    seen = set()
    ordered: List[Path] = []
    for path in sorted(candidates, key=lambda p: str(p)):
        key = str(path.resolve())
        if key in seen:
            continue
        seen.add(key)
        ordered.append(path)
    return ordered


def _cache_key_for_path(path: Path) -> str:
    resolved = str(path.resolve())
    return hashlib.sha256(resolved.encode("utf-8")).hexdigest()[:16]


def _hte_cache_dir(db_path: Path) -> Path:
    cache_root = PROJECT_ROOT / "results" / "hte_cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    return cache_root / _cache_key_for_path(db_path)


# Tracks the most recent non-memory source used by `_load_hte_database_cached`
# for a given target path ("disk" vs "rebuilt"). Memory hits bypass the
# function body due to `lru_cache`, so warm_hte_cache detects those separately.
_HTE_CACHE_LAST_LOAD_SOURCE: Dict[str, str] = {}


def _compute_hte_manifest(file_paths: List[Path]) -> Dict[str, Any]:
    entries = []
    for path in file_paths:
        try:
            stat = path.stat()
        except OSError:
            continue
        entries.append(
            {
                "path": str(path.resolve()),
                "mtime_ns": stat.st_mtime_ns,
                "size": stat.st_size,
            }
        )
    entries.sort(key=lambda item: item["path"])
    return {"version": 6, "files": entries}


def _to_row_index_array(indexer: Any) -> np.ndarray:
    if isinstance(indexer, np.ndarray):
        return indexer.astype(np.int64, copy=False)
    if isinstance(indexer, pd.Index):
        return indexer.to_numpy(dtype=np.int64, copy=True)
    if isinstance(indexer, range):
        return np.fromiter(indexer, dtype=np.int64)
    if isinstance(indexer, (list, tuple, set)):
        return np.asarray(list(indexer), dtype=np.int64)
    raise TypeError(f"Unsupported row index payload: {type(indexer)!r}")


def _materialize_group_frame(df: pd.DataFrame, indexer: np.ndarray) -> pd.DataFrame:
    if len(indexer) == 0:
        return df.iloc[0:0]
    return df.loc[indexer]


class _IndexedFrameLookup(Mapping[str, pd.DataFrame]):
    """Lazy DataFrame materialization backed by cached row-index arrays."""

    def __init__(
        self,
        df: pd.DataFrame,
        row_map: Dict[str, Any],
        *,
        cache_size: int = 256,
    ) -> None:
        self._df = df
        self._row_map: Dict[str, np.ndarray] = {
            str(key): _to_row_index_array(value)
            for key, value in dict(row_map).items()
        }
        self._cache_size = max(1, int(cache_size))
        self._frame_cache: OrderedDict[str, pd.DataFrame] = OrderedDict()

    def __getitem__(self, key: str) -> pd.DataFrame:
        normalized_key = str(key)
        cached = self._frame_cache.get(normalized_key)
        if cached is not None:
            self._frame_cache.move_to_end(normalized_key)
            return cached
        indexer = self._row_map[normalized_key]
        frame = _materialize_group_frame(self._df, indexer)
        self._frame_cache[normalized_key] = frame
        if len(self._frame_cache) > self._cache_size:
            self._frame_cache.popitem(last=False)
        return frame

    def __iter__(self) -> Iterator[str]:
        return iter(self._row_map)

    def __len__(self) -> int:
        return len(self._row_map)

    def items(self) -> Iterator[Tuple[str, pd.DataFrame]]:
        for key in self._row_map:
            yield key, self[key]

    def row_map(self) -> Dict[str, np.ndarray]:
        return dict(self._row_map)


def _is_frame_mapping(mapping: Any) -> bool:
    if not isinstance(mapping, dict) or not mapping:
        return False
    sample = next(iter(mapping.values()))
    return isinstance(sample, pd.DataFrame)


def _coerce_frame_lookup(
    df: pd.DataFrame,
    mapping: Any,
    *,
    cache_size: int,
) -> Any:
    if isinstance(mapping, _IndexedFrameLookup):
        return mapping
    if _is_frame_mapping(mapping):
        return dict(mapping)
    return _IndexedFrameLookup(df, dict(mapping or {}), cache_size=cache_size)


def _load_hte_cache(
    cache_dir: Path,
    manifest: Dict[str, Any],
) -> Optional[
    Tuple[pd.DataFrame, Dict[str, np.ndarray], Dict[str, Counter], Dict[str, np.ndarray]]
]:
    manifest_path = cache_dir / "manifest.json"
    payload_path = cache_dir / "hte_cache.pkl"
    if not manifest_path.exists() or not payload_path.exists():
        return None
    try:
        stored_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if stored_manifest != manifest:
        return None
    try:
        with payload_path.open("rb") as handle:
            payload = pickle.load(handle)
    except Exception:
        return None
    if not isinstance(payload, dict):
        return None
    df = payload.get("df")
    indexed_data = payload.get("indexed_row_map")
    reaction_type_patterns = payload.get("reaction_type_patterns")
    transformation_indices = payload.get("transformation_row_map")
    if df is None or indexed_data is None or reaction_type_patterns is None or transformation_indices is None:
        return None
    try:
        normalized_indexed = {
            str(key): _to_row_index_array(value)
            for key, value in dict(indexed_data).items()
        }
        normalized_transformations = {
            str(key): _to_row_index_array(value)
            for key, value in dict(transformation_indices).items()
        }
    except Exception:
        return None
    return df, normalized_indexed, reaction_type_patterns, normalized_transformations


def _save_hte_cache(
    cache_dir: Path,
    manifest: Dict[str, Any],
    df: pd.DataFrame,
    indexed_data: Dict[str, np.ndarray],
    reaction_type_patterns: Dict[str, Counter],
    transformation_indices: Dict[str, np.ndarray],
) -> None:
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        payload = {
            "df": df,
            "indexed_row_map": indexed_data,
            "reaction_type_patterns": reaction_type_patterns,
            "transformation_row_map": transformation_indices,
        }
        payload_path = cache_dir / "hte_cache.pkl"
        with payload_path.open("wb") as handle:
            pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
        manifest_path = cache_dir / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    except Exception:
        return


def _read_hte_csv(path: Path) -> pd.DataFrame:
    encodings = ("utf-8", "utf-8-sig", "cp1252", "latin-1")
    last_exc: Optional[Exception] = None
    for encoding in encodings:
        try:
            return pd.read_csv(path, encoding=encoding)
        except Exception as exc:
            last_exc = exc
    if last_exc:
        raise last_exc
    return pd.read_csv(path)


def _looks_like_precanonical_hte_frame(
    df: pd.DataFrame,
    source_path: Optional[Path],
) -> bool:
    """
    Heuristic fast-path detector for converter-produced canonical HTE CSVs.

    We only skip expensive key/event canonicalization when the frame already has:
    - canonical-ish file naming (contains `canonical`), and
    - populated Reaction_Key + Reaction_Events for rows with reaction_smiles, and
    - reactant columns already left-packed (no A-empty/B-filled style rows).
    """
    if df is None or df.empty:
        return False
    file_name = str(getattr(source_path, "name", "") or "").lower()
    if "canonical" not in file_name:
        return False

    required = {"reaction_smiles", "Reaction_Key", "Reaction_Events", "Reactant_A_Type", "Reactant_B_Type", "Reactant_C_Type"}
    if not required.issubset(set(df.columns)):
        return False

    reaction_smiles = df["reaction_smiles"].fillna("").astype(str).str.strip()
    reaction_keys = df["Reaction_Key"].fillna("").astype(str).str.strip()
    reaction_events = df["Reaction_Events"].fillna("").astype(str).str.strip()
    has_rxn = reaction_smiles.ne("")
    if has_rxn.any():
        if (reaction_keys[has_rxn] == "").any():
            return False
        if (reaction_events[has_rxn] == "").any():
            return False

    a = df["Reactant_A_Type"].fillna("").astype(str).str.strip()
    b = df["Reactant_B_Type"].fillna("").astype(str).str.strip()
    c = df["Reactant_C_Type"].fillna("").astype(str).str.strip()
    if ((a == "") & ((b != "") | (c != ""))).any():
        return False
    if ((b == "") & (c != "")).any():
        return False

    return True


def _normalize_hte_dataframe(df: pd.DataFrame, source_path: Optional[Path] = None) -> pd.DataFrame:
    df = df.copy()

    column_mapping = {
        "reaction_type": "Reaction_Type_Standardized",
        "detected_reaction_type": "Reaction_Type_Standardized",
        "reaction_category": "Reaction_Category",
        "rule_tier": "Rule_Tier",
        "reactant_1": "Reactant_A_Type",
        "reactant_2": "Reactant_B_Type",
        "reactant_3": "Reactant_C_Type",
        "reactant_3_type": "Reactant_C_Type",
        "yield": "AREA_TOTAL_REDUCED",
        "z_score": "z-Score",
        "catalyst": "Catalyst",
        "ligand": "Ligand",
        "base": "Base",
        "solvent": "Solvent",
        "additive": "Additive",
        "Spectator Groups": "spectator_groups",
        "reaction_events": "Reaction_Events",
        "Reaction Events": "Reaction_Events",
    }

    df = df.rename(columns={k: v for k, v in column_mapping.items() if k in df.columns})
    if "Source_Row" not in df.columns:
        df["Source_Row"] = list(range(len(df)))
    if "Reaction_Type_Standardized" not in df.columns:
        df["Reaction_Type_Standardized"] = ""
    if "reaction_id" in df.columns:
        existing_types = df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
        fallback_types = df["reaction_id"].fillna("").astype(str).str.strip()
        missing_type_mask = (existing_types == "") | existing_types.str.lower().eq("unknown")
        if missing_type_mask.any():
            df.loc[missing_type_mask, "Reaction_Type_Standardized"] = fallback_types[missing_type_mask]

    precanonical_fast_path = _looks_like_precanonical_hte_frame(df, source_path)

    # Generate Reaction_Key from reaction_smiles when missing/invalid.
    # Skip for converter-produced canonical datasets that already carry keys.
    if (not precanonical_fast_path) and "reaction_smiles" in df.columns:
        if "Reaction_Key" not in df.columns:
            df["Reaction_Key"] = ""
        reaction_smiles_series = df["reaction_smiles"].fillna("").astype(str).str.strip()
        reaction_key_series = df["Reaction_Key"].fillna("").astype(str).str.strip()
        needs_key_mask = reaction_smiles_series.ne("") & (
            reaction_key_series.eq("") | reaction_key_series.str.lower().eq("nan")
        )
        if needs_key_mask.any():
            def _gen_rxn_key(smiles_val: Any) -> str:
                s = str(smiles_val).strip()
                if not s or s.lower() == "nan":
                    return ""
                try:
                    context = featurize_reaction(s)
                    return str(context.get("reaction_key") or "")
                except Exception:
                    return ""

            df.loc[needs_key_mask, "Reaction_Key"] = reaction_smiles_series[needs_key_mask].apply(_gen_rxn_key)

    if "Reaction_Key" in df.columns:
        if "Reaction_Events" not in df.columns:
            df["Reaction_Events"] = ""
        raw_key_series = df["Reaction_Key"].fillna("").astype(str).str.strip()
        raw_events_series = df["Reaction_Events"].fillna("").astype(str).str.strip()

        if precanonical_fast_path:
            # Trust converter output and avoid per-row canonicalization of keys/events.
            # We still normalize trivial whitespace and derive match keys.
            normalized_keys = raw_key_series.tolist()
            normalized_events = raw_events_series.tolist()
        else:
            normalized_keys = []
            normalized_events = []
            for raw_key, raw_events in zip(raw_key_series.tolist(), raw_events_series.tolist()):
                minimal_key = canonicalize_reaction_key_minimal(raw_key)
                normalized_keys.append(minimal_key)
                if raw_events and raw_events.lower() != "nan":
                    normalized_text = normalize_reaction_events_text(raw_events)
                    normalized_events.append(normalized_text or raw_events)
                else:
                    events_payload = build_reaction_events_payload(raw_key)
                    normalized_events.append(
                        serialize_reaction_events_payload(events_payload)
                    )
        df["Reaction_Key"] = normalized_keys
        df["Reaction_Events"] = normalized_events
        if "Reaction_Events_Key" not in df.columns or not precanonical_fast_path:
            df["Reaction_Events_Key"] = [
                _reaction_events_to_match_key(value) for value in normalized_events
            ]
        else:
            df["Reaction_Events_Key"] = (
                df["Reaction_Events_Key"].fillna("").astype(str).str.strip()
            )
        valid_mask = df["Reaction_Key"].fillna("").astype(str).str.contains("->", regex=False)
        if (~valid_mask).any():
            df.loc[~valid_mask, "Reaction_Key"] = ""
        if "Reaction_Type_Standardized" not in df.columns or not any(df["Reaction_Type_Standardized"]):
            df["Reaction_Type_Standardized"] = df["Reaction_Key"]
        if "Reactant_Types_Key" not in df.columns or not any(df["Reactant_Types_Key"]):
            df["Reactant_Types_Key"] = df["Reaction_Key"]
    if "Reaction_Events" in df.columns and "Reaction_Events_Key" not in df.columns:
        events_series = df["Reaction_Events"].fillna("").astype(str).str.strip()
        df["Reaction_Events_Key"] = [
            _reaction_events_to_match_key(value) for value in events_series.tolist()
        ]

    required_cols = [
        "Reaction_Type_Standardized", "Reactant_A_Type", "Reactant_B_Type",
        "Reactant_C_Type",
        "Catalyst", "Ligand", "Base", "Solvent", "Additive",
        "Secondary Solvent", "Coupling Reagent", "AREA_TOTAL_REDUCED", "z-Score",
        "Reactant_A_Category", "Reactant_B_Category", "Reaction_Category",
        "Source_File", "Source_Group", "Source_Row", "spectator_groups",
        "Reactant_Signature_Core", "Reactant_Signature_Ext", "Rule_Tier", "Reaction_Events",
        "Reaction_Events_Key",
    ]
    for col in required_cols:
        if col not in df.columns:
            if col in ("Source_File", "Source_Group"):
                df[col] = ""
            elif col == "Rule_Tier":
                df[col] = 0.0
            else:
                df[col] = "" if col not in ["AREA_TOTAL_REDUCED", "z-Score"] else 0.0

    # Enforce numeric dtypes for score/yield fields to prevent runtime math errors
    # (e.g., accidental string values from malformed CSV rows).
    df["AREA_TOTAL_REDUCED"] = pd.to_numeric(df["AREA_TOTAL_REDUCED"], errors="coerce").fillna(0.0)
    df["z-Score"] = pd.to_numeric(df["z-Score"], errors="coerce").fillna(0.0)
    df["Rule_Tier"] = df["Rule_Tier"].apply(_coerce_rule_tier_value).astype(float)
    for optional_numeric in ("temperature_C", "time_h", "pressure_bar", "concentration_M"):
        if optional_numeric in df.columns:
            df[optional_numeric] = pd.to_numeric(df[optional_numeric], errors="coerce")

    df["Reaction_Type_Standardized"] = (
        df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
    )

    # Keep reactant columns as object/string so downstream normalization can
    # safely write empty strings without dtype warnings on float columns.
    for reactant_col in ("Reactant_A_Type", "Reactant_B_Type", "Reactant_C_Type"):
        if reactant_col in df.columns:
            df[reactant_col] = df[reactant_col].astype(object)

    def _normalize_reactants_row(row: pd.Series) -> pd.Series:
        # `DataFrame.apply(axis=1)` can hand us a float-typed Series for some rows
        # (e.g. when reactant cells are NaN-heavy). Cast to object before writing
        # normalized string values to avoid pandas dtype-assignment warnings.
        row = row.copy().astype(object)
        raw = [row.get("Reactant_A_Type"), row.get("Reactant_B_Type"), row.get("Reactant_C_Type")]
        cleaned: List[str] = []
        for value in raw:
            if value is None:
                continue
            if isinstance(value, float) and pd.isna(value):
                continue
            text = str(value).strip()
            if not text or text.lower() == "nan":
                continue
            cleaned.append(text)
        a_val = cleaned[0] if len(cleaned) > 0 else ""
        b_val = cleaned[1] if len(cleaned) > 1 else ""
        c_val = cleaned[2] if len(cleaned) > 2 else ""
        row["Reactant_A_Type"] = a_val
        row["Reactant_B_Type"] = b_val
        row["Reactant_C_Type"] = c_val
        row["_reactant_count"] = len(cleaned)
        return row

    for _col in ["Reactant_A_Type", "Reactant_B_Type", "Reactant_C_Type"]:
        if _col in df.columns:
            df[_col] = df[_col].fillna("").astype(str).replace("nan", "").astype(object)
            df[_col] = df[_col].astype(str).str.strip().astype(object)
    if not precanonical_fast_path:
        df = df.apply(_normalize_reactants_row, axis=1)
    df["Intramolecular_Likely"] = df.apply(
        lambda row: _intramolecular_likely_from_fields(
            row.get("Reactant_A_Type"),
            row.get("Reactant_B_Type"),
            row.get("Reactant_C_Type"),
        ),
        axis=1,
    )
    if "Is_Intramolecular" not in df.columns:
        df["Is_Intramolecular"] = df["Intramolecular_Likely"]
    df = df.drop(columns=["_reactant_count"], errors="ignore")

    existing_reactant_types_key = (
        df["Reactant_Types_Key"].fillna("").astype(str).str.strip()
        if "Reactant_Types_Key" in df.columns
        else pd.Series([""] * len(df), index=df.index)
    )
    if (existing_reactant_types_key == "").any():
        computed_reactant_types_key = df.apply(
            lambda row: _reactant_key(
                [row.get("Reactant_A_Type"), row.get("Reactant_B_Type"), row.get("Reactant_C_Type")]
            ),
            axis=1,
        )
        if "Reactant_Types_Key" not in df.columns:
            df["Reactant_Types_Key"] = computed_reactant_types_key
        else:
            missing_rtk_mask = existing_reactant_types_key == ""
            if missing_rtk_mask.any():
                df.loc[missing_rtk_mask, "Reactant_Types_Key"] = computed_reactant_types_key[missing_rtk_mask]
    else:
        df["Reactant_Types_Key"] = existing_reactant_types_key

    if "Reaction_Key" in df.columns:
        keys = df["Reaction_Key"].fillna("").astype(str).str.strip()
        signatures = keys.apply(_reaction_key_to_signatures)
        df["Reactant_Signature_Core"] = signatures.apply(lambda pair: pair[0])
        df["Reactant_Signature_Ext"] = signatures.apply(lambda pair: pair[1])

    df["Reactant_Signature_Core"] = df["Reactant_Signature_Core"].fillna("")
    df["Reactant_Signature_Ext"] = df["Reactant_Signature_Ext"].fillna("")

    core_from_reactants = df.apply(
        lambda row: _reactant_types_to_signature(
            [row.get("Reactant_A_Type"), row.get("Reactant_B_Type"), row.get("Reactant_C_Type")]
        ),
        axis=1,
    )
    core_mask = df["Reactant_Signature_Core"].astype(str).str.strip() == ""
    if core_mask.any():
        df.loc[core_mask, "Reactant_Signature_Core"] = core_from_reactants[core_mask]
    ext_mask = df["Reactant_Signature_Ext"].astype(str).str.strip() == ""
    if ext_mask.any():
        df.loc[ext_mask, "Reactant_Signature_Ext"] = df.loc[ext_mask, "Reactant_Signature_Core"]

    if source_path is not None:
        df["Source_File"] = _format_source_path(source_path)
        df["Source_Group"] = _infer_source_group(source_path)

    return df


@lru_cache(maxsize=4)
def _load_hte_database_cached(
    hte_db_path: str,
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray], Dict[str, Counter], Dict[str, np.ndarray]]:
    """Load and index the HTE database once per path (cached)."""
    db_path = Path(hte_db_path)
    if not db_path.exists():
        raise FileNotFoundError(f"HTE database not found: {db_path}")

    file_paths = _collect_hte_files(db_path)
    if not file_paths:
        raise FileNotFoundError(f"No HTE CSV files found under: {db_path}")

    manifest = _compute_hte_manifest(file_paths)
    cache_dir = _hte_cache_dir(db_path)
    cached = _load_hte_cache(cache_dir, manifest)
    if cached is not None:
        _HTE_CACHE_LAST_LOAD_SOURCE[str(db_path)] = "disk"
        return cached

    frames: List[pd.DataFrame] = []
    for path in file_paths:
        frame = _read_hte_csv(path)
        frame = _normalize_hte_dataframe(frame, source_path=path)
        frames.append(frame)

    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    df = _ensure_rule_tier_column(df)
    LOGGER.info("Loaded HTE database: %s experiments from %s files", len(df), len(file_paths))

    indexed_data: Dict[str, np.ndarray] = {}
    reaction_type_patterns: Dict[str, Counter] = {}
    transformation_indices: Dict[str, np.ndarray] = {}

    LOGGER.info("Building reactant type indices...")

    motif_sets = _load_motif_sets()
    scope_map = _load_scope_map()
    key_to_indices: Dict[str, Set[int]] = defaultdict(set)
    for row in df.itertuples(index=True):
        keys = _expand_reactant_keys(
            getattr(row, "Reactant_A_Type", ""),
            getattr(row, "Reactant_B_Type", ""),
            motif_sets,
            scope_map,
        )
        if not keys:
            key = _reactant_key([getattr(row, "Reactant_A_Type", ""), getattr(row, "Reactant_B_Type", "")])
            if key:
                keys = [key]
        for key in keys:
            key_to_indices[key].add(row.Index)

    for key, indices in key_to_indices.items():
        if not indices:
            continue
        group_df = df.loc[sorted(indices)]
        indexed_data[key] = group_df.index.to_numpy(dtype=np.int64, copy=True)

        rxn_types = group_df["Reaction_Type_Standardized"].value_counts()
        reaction_type_patterns[key] = Counter(rxn_types.to_dict())
    
    # Build transformation-aware index (Reaction_Key preferred when available)
    df["Reaction_Type_Standardized"] = df["Reaction_Type_Standardized"].fillna("Unknown")
    if "Reaction_Key" in df.columns:
        df["Reaction_Key"] = df["Reaction_Key"].fillna("").astype(str).str.strip()
        keyed_df = df[df["Reaction_Key"] != ""]
        for key, group_df in keyed_df.groupby("Reaction_Key"):
            transformation_indices[key] = group_df.index.to_numpy(dtype=np.int64, copy=True)
        unkeyed_df = df[df["Reaction_Key"] == ""]
        if not unkeyed_df.empty:
            if "Reactant_Signature_Core" in unkeyed_df.columns:
                sig_series = unkeyed_df["Reactant_Signature_Core"].fillna("").astype(str).str.strip()
                sig_df = unkeyed_df[sig_series != ""]
                for key, group_df in sig_df.groupby("Reactant_Signature_Core"):
                    transformation_indices[key] = group_df.index.to_numpy(dtype=np.int64, copy=True)
                unkeyed_df = unkeyed_df[sig_series == ""]
            if not unkeyed_df.empty:
                for key, group_df in unkeyed_df.groupby("Reaction_Type_Standardized"):
                    transformation_indices[key] = group_df.index.to_numpy(dtype=np.int64, copy=True)
    else:
        if "Reactant_Signature_Core" in df.columns:
            sig_series = df["Reactant_Signature_Core"].fillna("").astype(str).str.strip()
            sig_df = df[sig_series != ""]
            for key, group_df in sig_df.groupby("Reactant_Signature_Core"):
                transformation_indices[key] = group_df.index.to_numpy(dtype=np.int64, copy=True)
            remaining = df[sig_series == ""]
            if not remaining.empty:
                for key, group_df in remaining.groupby("Reaction_Type_Standardized"):
                    transformation_indices[key] = group_df.index.to_numpy(dtype=np.int64, copy=True)
        else:
            grouped_rxn = df.groupby("Reaction_Type_Standardized")
            for key, group_df in grouped_rxn:
                transformation_indices[key] = group_df.index.to_numpy(dtype=np.int64, copy=True)

    # Add event-aware indices so noisy/missing Reaction_Key values can still match.
    if "Reaction_Events_Key" in df.columns:
        event_key_series = df["Reaction_Events_Key"].fillna("").astype(str).str.strip()
        event_df = df[event_key_series != ""]
        for key, group_df in event_df.groupby("Reaction_Events_Key"):
            if key and key not in transformation_indices:
                transformation_indices[key] = group_df.index.to_numpy(dtype=np.int64, copy=True)

    LOGGER.info(
        "Indexed %s reactant combinations and %s transformation types",
        len(indexed_data),
        len(transformation_indices),
    )
    _save_hte_cache(cache_dir, manifest, df, indexed_data, reaction_type_patterns, transformation_indices)
    _HTE_CACHE_LAST_LOAD_SOURCE[str(db_path)] = "rebuilt"
    return df, indexed_data, reaction_type_patterns, transformation_indices


def _resolve_warm_cache_targets(db_path: Path, source_group: Optional[str]) -> List[Path]:
    if not source_group:
        return [db_path]
    label = str(source_group).strip().lower()
    if not label or label == "all":
        return [db_path]

    if label in ("datasets", "dataset", "lit"):
        label = "literature"
    elif label in ("protocol", "protocols"):
        label = "literature"
    elif label in ("motif", "motifs", "experiment", "experiements"):
        label = "motif"

    if not db_path.is_dir():
        return [db_path]

    if label == "motif":
        for subdir in ("motif", "experiments"):
            sub_path = db_path / subdir
            canonical = sub_path / "HTE_canonical.csv"
            if canonical.exists():
                return [canonical]
            if sub_path.exists():
                return [sub_path]
        return [db_path]
    sub_path = db_path / label
    if sub_path.exists():
        return [sub_path]
    return [db_path]


def check_hte_cache_status(
    hte_db_path: str,
    *,
    source_group: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Quick validity check for the on-disk HTE index cache — no data is loaded.

    Returns a dict with:
        ``valid``    – True if all targets have a current disk (or memory) cache.
        ``targets``  – Per-target status list; each entry has ``status`` which is
                       one of ``"memory"``, ``"disk"``, ``"missing"``, ``"stale"``,
                       or ``"no_files"``.
    """
    base_path = Path(hte_db_path)
    if not base_path.exists():
        return {"valid": False, "reason": "path_not_found", "targets": []}

    normalized_sg = str(source_group or "all").strip().lower()
    targets = _resolve_warm_cache_targets(base_path, normalized_sg)
    results: List[Dict[str, Any]] = []

    for target in targets:
        target_str = str(target)
        target_key = str(Path(target))

        # --- memory hit? ---
        last_source = _HTE_CACHE_LAST_LOAD_SOURCE.get(target_key, "")
        if last_source == "memory":
            results.append({"target": target_str, "status": "memory"})
            continue

        # --- disk cache check (manifest only, no pickle load) ---
        try:
            file_paths = _collect_hte_files(target)
        except Exception:
            results.append({"target": target_str, "status": "no_files"})
            continue
        if not file_paths:
            results.append({"target": target_str, "status": "no_files"})
            continue

        try:
            manifest = _compute_hte_manifest(file_paths)
        except Exception:
            results.append({"target": target_str, "status": "no_files"})
            continue

        cache_dir = _hte_cache_dir(Path(target))
        manifest_path = cache_dir / "manifest.json"
        payload_path = cache_dir / "hte_cache.pkl"

        if not manifest_path.exists() or not payload_path.exists():
            results.append({"target": target_str, "status": "missing", "cache_dir": str(cache_dir)})
            continue

        try:
            stored_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            results.append({"target": target_str, "status": "stale", "cache_dir": str(cache_dir)})
            continue

        if stored_manifest != manifest:
            results.append({"target": target_str, "status": "stale", "cache_dir": str(cache_dir)})
            continue

        # If we previously loaded it (disk or rebuilt), treat as memory now
        if last_source in {"disk", "rebuilt"}:
            results.append({"target": target_str, "status": "memory", "cache_dir": str(cache_dir)})
        else:
            results.append({"target": target_str, "status": "disk", "cache_dir": str(cache_dir)})

    all_valid = all(r["status"] in {"memory", "disk"} for r in results) and bool(results)
    return {"valid": all_valid, "targets": results}


def warm_hte_cache(
    hte_db_path: str = "data/HTE_db",
    *,
    source_group: Optional[str] = None,
    clear_memory_cache: bool = False,
) -> Dict[str, Any]:
    """
    Prebuild on-disk and in-memory HTE indices so first recommendation calls are fast.

    Returns:
        Dict summary containing warmed targets, elapsed time, and index sizes.
    """
    if clear_memory_cache:
        _load_hte_database_cached.cache_clear()

    base_path = Path(hte_db_path)
    if not base_path.exists():
        raise FileNotFoundError(f"HTE database not found: {base_path}")

    targets = _resolve_warm_cache_targets(base_path, source_group)
    summary_targets: List[Dict[str, Any]] = []

    for target in targets:
        target_key = str(Path(target))
        cache_info_before = _load_hte_database_cached.cache_info()
        start = pd.Timestamp.utcnow()
        df, indexed_data, reaction_type_patterns, transformation_indices = _load_hte_database_cached(str(target))
        elapsed_s = (pd.Timestamp.utcnow() - start).total_seconds()
        cache_info_after = _load_hte_database_cached.cache_info()
        if cache_info_after.hits > cache_info_before.hits:
            cache_source = "memory"
        else:
            cache_source = _HTE_CACHE_LAST_LOAD_SOURCE.get(target_key, "unknown")
        summary_targets.append(
            {
                "target": str(target),
                "cache_dir": str(_hte_cache_dir(Path(target))),
                "cache_source": cache_source,
                "cache_reused": cache_source in {"memory", "disk"},
                "num_rows": int(len(df)),
                "reactant_index_keys": int(len(indexed_data)),
                "reaction_type_pattern_keys": int(len(reaction_type_patterns)),
                "transformation_index_keys": int(len(transformation_indices)),
                "elapsed_s": round(float(elapsed_s), 3),
            }
        )

    return {
        "base_path": str(base_path),
        "source_group": source_group or "all",
        "targets": summary_targets,
    }


def _ensure_list(values: Any) -> List[str]:
    if values is None:
        return []
    if isinstance(values, list):
        return [str(v).strip() for v in values if str(v).strip()]
    if isinstance(values, str):
        text = values.strip()
        return [text] if text else []
    text = str(values).strip()
    return [text] if text else []


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


def _format_list(values: Any) -> str:
    items = _dedupe_list(_ensure_list(values))
    return " / ".join(items)


def _first_nonempty_text(values: Iterable[Any]) -> Optional[str]:
    for value in values:
        if value is None:
            continue
        if isinstance(value, float) and pd.isna(value):
            continue
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        return text
    return None


def _normalize_reaction_key(value: Any) -> Optional[str]:
    """Normalize reaction key text; return None for empty placeholders."""
    if value is None:
        return None
    text = canonicalize_reaction_key_minimal(value)
    if not text or text.lower() == "none":
        return None
    compact = text.replace(" ", "")
    if compact in {"[]->[]||[]", "CRK-v1|[]->[]", "|[]->[]"}:
        return None
    return text


def _aggregate_spectator_groups(
    values: pd.Series,
    *,
    query_groups: Optional[Set[str]] = None,
    max_items: int = 4,
) -> str:
    if values is None:
        return ""
    token_counter: Counter[str] = Counter()
    for raw in values:
        if raw is None:
            continue
        if isinstance(raw, float) and pd.isna(raw):
            continue
        text = str(raw).strip()
        if not text:
            continue
        tokens = _split_group_tokens(text)
        for token in tokens:
            token_counter[token] += 1
    if not token_counter:
        return ""

    def _sort_key(group_id: str) -> Tuple[float, int, str]:
        return (
            -spectator_group_weight(group_id),
            -int(token_counter.get(group_id, 0)),
            str(group_id),
        )

    if query_groups:
        overlap = [gid for gid in token_counter if gid in query_groups]
        if overlap:
            picked = sorted(overlap, key=_sort_key)[:max_items]
            return " / ".join(picked)

    weighted = [gid for gid in token_counter if spectator_group_weight(gid) > 0]
    if weighted:
        picked = sorted(weighted, key=_sort_key)[:max_items]
        return " / ".join(picked)

    most_common = token_counter.most_common(max_items)
    return " / ".join(group for group, _ in most_common)


def _normalize_reaction_type_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _clean_reaction_label(value: Optional[str]) -> str:
    text = _normalize_reaction_type_value(value)
    if not text or text.lower() == "unknown":
        return ""
    return text


@lru_cache(maxsize=1024)
def _resolve_reaction_type_label(label: Optional[str]) -> str:
    text = _normalize_reaction_type_value(label)
    if not text:
        return ""
    if _reaction_catalog is None:
        return text
    candidates: List[str] = [text, text.lower()]

    cleaned = re.sub(r"\([^)]*\)", "", text).strip()
    if cleaned:
        candidates.extend([cleaned, cleaned.lower()])

        trimmed = re.sub(r"[_\-\s]*sample\s*\d+\s*$", "", cleaned, flags=re.IGNORECASE).strip()
        trimmed = re.sub(
            r"[_\-\s]*dataset[_\-\s]*converted(?:[_\-\s]*\d+)?\s*$",
            "",
            trimmed,
            flags=re.IGNORECASE,
        ).strip()
        if trimmed:
            candidates.extend([trimmed, trimmed.lower()])
            candidates.extend(
                [
                    trimmed.replace(" ", "_"),
                    trimmed.replace(" ", "-"),
                    trimmed.replace("-", "_"),
                    trimmed.replace("_", "-"),
                ]
            )

    seen: Set[str] = set()
    for candidate in candidates:
        cand = str(candidate).strip()
        if not cand:
            continue
        key = cand.lower()
        if key in seen:
            continue
        seen.add(key)
        resolved = _reaction_catalog.resolve_reaction_type(cand)
        if resolved:
            return resolved

    low = text.lower()
    if "suzuki" in low:
        resolved = _reaction_catalog.resolve_reaction_type("suzuki")
        if resolved:
            return resolved
    return text


@lru_cache(maxsize=256)
def _required_catalyst_classes_for_reaction_type(label: Optional[str]) -> Tuple[str, ...]:
    resolved = _resolve_reaction_type_label(label)
    if not resolved or _reaction_catalog is None:
        return ()
    definition = _reaction_catalog.get_reaction_type(resolved)
    if definition is None:
        return ()
    catalysts = [str(value).strip() for value in getattr(definition, "catalysts", []) or []]
    return tuple(value for value in catalysts if value)


def _series_has_text(values: pd.Series) -> pd.Series:
    text = values.fillna("").astype(str).str.strip()
    return (text != "") & (text.str.lower() != "nan")


def _row_has_catalyst_evidence(df: pd.DataFrame) -> pd.Series:
    mask = pd.Series([False] * len(df), index=df.index)
    if "Catalyst" in df.columns:
        mask = mask | _series_has_text(df["Catalyst"])
    if "Catalyst_Type" in df.columns:
        mask = mask | _series_has_text(df["Catalyst_Type"])
    return mask


def _enforce_required_catalyst_quality(
    matched_df: pd.DataFrame,
    *,
    reaction_type_filter: Optional[str],
    predicted_reaction_type: Optional[str],
    reaction_type_confidence: float,
) -> Tuple[pd.DataFrame, Optional[str], Tuple[str, ...], int, int]:
    if matched_df.empty or "Reaction_Type_Standardized" not in matched_df.columns:
        return matched_df, None, (), 0, 0

    target_reaction = reaction_type_filter
    if (
        not target_reaction
        and predicted_reaction_type
        and predicted_reaction_type != "Unknown"
        and reaction_type_confidence >= 0.5
    ):
        target_reaction = predicted_reaction_type
    target_reaction = _resolve_reaction_type_label(target_reaction)
    required_catalysts = _required_catalyst_classes_for_reaction_type(target_reaction)
    if not target_reaction or not required_catalysts:
        return matched_df, None, (), 0, 0

    type_series = matched_df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
    family_mask = type_series.apply(
        lambda value: _resolve_reaction_type_label(value) == target_reaction if value else False
    )
    if not family_mask.any():
        return matched_df, target_reaction, required_catalysts, 0, 0

    catalyst_evidence_mask = _row_has_catalyst_evidence(matched_df)
    missing_required_catalyst_mask = family_mask & ~catalyst_evidence_mask
    if not missing_required_catalyst_mask.any():
        return matched_df, target_reaction, required_catalysts, 0, 0

    complete_family_mask = family_mask & catalyst_evidence_mask
    if complete_family_mask.any():
        filtered_rows = int(missing_required_catalyst_mask.sum())
        return (
            matched_df.loc[~missing_required_catalyst_mask].copy(),
            target_reaction,
            required_catalysts,
            filtered_rows,
            0,
        )

    return matched_df, target_reaction, required_catalysts, 0, int(missing_required_catalyst_mask.sum())


_CONDITION_FIELD_COLUMNS: Dict[str, Tuple[str, ...]] = {
    "catalyst": ("Catalyst", "Catalyst_Type"),
    "ligand": ("Ligand",),
    "base": ("Base",),
    "solvent": ("Solvent",),
}

_FAMILY_CONDITION_COMPLETENESS_RULES: Dict[str, Dict[str, Any]] = {
    "Suzuki_miyaura": {
        "default": ("catalyst", "base", "solvent"),
        "penalties": {"catalyst": 0.15, "base": 0.55, "solvent": 0.70},
    },
    "Miyaura_borylation": {
        "default": ("catalyst", "base", "solvent"),
        "penalties": {"catalyst": 0.15, "base": 0.60, "solvent": 0.70},
    },
    "C_N_Coupling": {
        "default": ("catalyst", "base", "solvent"),
        "by_catalyst": {
            "Pd": ("catalyst", "ligand", "base", "solvent"),
            "Ni": ("catalyst", "ligand", "base", "solvent"),
            "Cu": ("catalyst", "base", "solvent"),
        },
        "penalties": {"catalyst": 0.15, "ligand": 0.50, "base": 0.60, "solvent": 0.75},
    },
    "C_O_Coupling": {
        "default": ("catalyst", "base", "solvent"),
        "by_catalyst": {
            "Pd": ("catalyst", "ligand", "base", "solvent"),
            "Ni": ("catalyst", "ligand", "base", "solvent"),
            "Cu": ("catalyst", "base", "solvent"),
        },
        "penalties": {"catalyst": 0.15, "ligand": 0.55, "base": 0.60, "solvent": 0.75},
    },
    "C_S_Coupling": {
        "default": ("catalyst", "base", "solvent"),
        "by_catalyst": {
            "Pd": ("catalyst", "ligand", "base", "solvent"),
            "Ni": ("catalyst", "ligand", "base", "solvent"),
            "Cu": ("catalyst", "base", "solvent"),
        },
        "penalties": {"catalyst": 0.15, "ligand": 0.55, "base": 0.60, "solvent": 0.75},
    },
    "Chan_Lam_C_N_Coupling": {
        "default": ("catalyst", "solvent"),
        "penalties": {"catalyst": 0.20, "solvent": 0.75},
    },
}

_METAL_SYMBOL_ALIASES: Dict[str, Tuple[str, ...]] = {
    "Pd": ("palladium",),
    "Ni": ("nickel",),
    "Cu": ("copper",),
    "Ru": ("ruthenium",),
    "Rh": ("rhodium",),
    "Ir": ("iridium",),
    "Fe": ("iron",),
    "Co": ("cobalt",),
    "Ag": ("silver",),
    "Au": ("gold",),
    "Pt": ("platinum",),
    "Zn": ("zinc",),
    "Sc": ("scandium",),
}

_PLAUSIBILITY_MECHANISM_PENALTY = 0.35
_PLAUSIBILITY_SOURCE_EXEMPTIONS = {"rules"}


def _condition_text_present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, float) and pd.isna(value):
        return False
    text = str(value).strip()
    if not text:
        return False
    return text.lower() != "nan"


def _extract_catalyst_metal_tokens(value: Any) -> Tuple[str, ...]:
    if not _condition_text_present(value):
        return ()

    text = str(value)
    lower_text = text.lower()
    found: List[str] = []
    for symbol, aliases in _METAL_SYMBOL_ALIASES.items():
        symbol_pattern = re.compile(rf"(?<![A-Za-z]){re.escape(symbol)}(?=[^a-z]|$)")
        if symbol_pattern.search(text) or any(alias in lower_text for alias in aliases):
            found.append(symbol)
    return tuple(sorted(set(found)))


def _row_contains_unexpected_catalyst_metals(
    row: Mapping[str, Any],
    *,
    allowed_classes: Tuple[str, ...],
) -> bool:
    allowed = {str(value).strip() for value in allowed_classes if str(value).strip()}
    if not allowed:
        return False

    tokens = set(_extract_catalyst_metal_tokens(row.get("Catalyst")))
    tokens.update(_extract_catalyst_metal_tokens(row.get("Catalyst_Type")))
    if not tokens:
        return False

    return bool(tokens - allowed)


def _row_has_implausible_snar_annotation(row: Mapping[str, Any]) -> bool:
    events_text = str(row.get("Reaction_Events") or "").strip()
    if "mech=snar" not in events_text.lower():
        return False

    evidence_text = " ".join(
        str(row.get(field) or "")
        for field in (
            "Reaction_Key",
            "Reaction_Events",
            "reactant_1",
            "reactant_2",
            "reactant_3",
            "formed_motifs",
            "spectator_groups",
        )
    ).upper()

    has_fluoride_site = any(token in evidence_text for token in ("AR-F", "HETEROAR-F", "LG=F"))
    has_iodide_or_bromide_site = any(
        token in evidence_text
        for token in ("AR-I", "HETEROAR-I", "AR-BR", "HETEROAR-BR", "LG=I", "LG=BR")
    )
    return has_iodide_or_bromide_site and not has_fluoride_site


def _row_condition_field_present(row: Mapping[str, Any], field_name: str) -> bool:
    for column in _CONDITION_FIELD_COLUMNS.get(field_name, ()):
        if column in row and _condition_text_present(row[column]):
            return True
    return False


def _infer_catalyst_class(
    catalyst_value: Any,
    catalyst_type_value: Any,
    allowed_classes: Tuple[str, ...],
) -> Optional[str]:
    candidates = [catalyst_type_value, catalyst_value]
    normalized_allowed = [str(value).strip() for value in allowed_classes if str(value).strip()]
    if not normalized_allowed:
        return None

    for candidate in candidates:
        if not _condition_text_present(candidate):
            continue
        text = str(candidate).strip()
        text_upper = text.upper()
        text_tokens = set(re.findall(r"[A-Z][a-z]?|[A-Z]{2,}", text))
        compact_tokens = set(re.findall(r"[A-Z]{1,3}", re.sub(r"[^A-Za-z]", " ", text_upper)))
        for catalyst_class in normalized_allowed:
            cls_upper = catalyst_class.upper()
            if (
                text_upper.startswith(cls_upper)
                or f" {cls_upper}" in f" {text_upper}"
                or cls_upper in text_tokens
                or cls_upper in compact_tokens
            ):
                return catalyst_class
    return None


def _required_condition_fields_for_reaction_type(
    reaction_type: Optional[str],
    catalyst_class: Optional[str] = None,
) -> Tuple[Tuple[str, ...], Dict[str, float]]:
    resolved = _resolve_reaction_type_label(reaction_type)
    profile = _FAMILY_CONDITION_COMPLETENESS_RULES.get(resolved)
    if not profile:
        return (), {}
    penalties = {
        str(key): float(value)
        for key, value in (profile.get("penalties") or {}).items()
    }
    required = tuple(str(value) for value in (profile.get("default") or ()) if str(value))
    by_catalyst = profile.get("by_catalyst") or {}
    if catalyst_class:
        required = tuple(str(value) for value in (by_catalyst.get(catalyst_class) or required) if str(value))
    return required, penalties


def _apply_required_condition_quality_penalties(
    matched_df: pd.DataFrame,
    *,
    reaction_type_filter: Optional[str],
    predicted_reaction_type: Optional[str],
    reaction_type_confidence: float,
) -> Tuple[pd.DataFrame, Optional[str], Dict[str, int], int]:
    if matched_df.empty or "Reaction_Type_Standardized" not in matched_df.columns:
        return matched_df, None, {}, 0

    target_reaction = reaction_type_filter
    if (
        not target_reaction
        and predicted_reaction_type
        and predicted_reaction_type != "Unknown"
        and reaction_type_confidence >= 0.5
    ):
        target_reaction = predicted_reaction_type
    target_reaction = _resolve_reaction_type_label(target_reaction)
    required_catalysts = _required_catalyst_classes_for_reaction_type(target_reaction)
    default_required_fields, _ = _required_condition_fields_for_reaction_type(target_reaction)
    if not target_reaction or not default_required_fields:
        return matched_df, None, {}, 0

    type_series = matched_df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
    family_mask = type_series.apply(
        lambda value: _resolve_reaction_type_label(value) == target_reaction if value else False
    )
    if not family_mask.any():
        return matched_df, target_reaction, {}, 0

    working_df = matched_df.copy()
    existing_match_scores = (
        pd.to_numeric(working_df["match_score"], errors="coerce").fillna(1.0)
        if "match_score" in working_df.columns
        else pd.Series([1.0] * len(working_df), index=working_df.index, dtype=float)
    )
    quality_multiplier = pd.Series([1.0] * len(working_df), index=working_df.index, dtype=float)
    missing_field_labels = pd.Series([""] * len(working_df), index=working_df.index, dtype=object)
    missing_field_counts: Counter[str] = Counter()
    penalized_rows = 0

    source_series = (
        working_df["Source_Group"].apply(_normalize_source_group)
        if "Source_Group" in working_df.columns
        else pd.Series([""] * len(working_df), index=working_df.index)
    )
    catalyst_series = (
        working_df["Catalyst"] if "Catalyst" in working_df.columns else pd.Series([""] * len(working_df), index=working_df.index)
    )
    catalyst_type_series = (
        working_df["Catalyst_Type"] if "Catalyst_Type" in working_df.columns else pd.Series([""] * len(working_df), index=working_df.index)
    )

    for idx in working_df.index[family_mask]:
        if source_series.loc[idx] == "rules":
            continue
        catalyst_class = _infer_catalyst_class(
            catalyst_series.loc[idx],
            catalyst_type_series.loc[idx],
            required_catalysts,
        )
        required_fields, penalty_map = _required_condition_fields_for_reaction_type(target_reaction, catalyst_class)
        if not required_fields:
            continue

        row = working_df.loc[idx]
        missing_fields = [field for field in required_fields if not _row_condition_field_present(row, field)]
        if not missing_fields:
            continue

        penalized_rows += 1
        multiplier = 1.0
        for field in missing_fields:
            missing_field_counts[field] += 1
            multiplier *= float(penalty_map.get(field, 0.75))
        quality_multiplier.loc[idx] = max(0.0, min(1.0, multiplier))
        missing_field_labels.loc[idx] = "|".join(sorted(missing_fields))

    if penalized_rows == 0:
        return matched_df, target_reaction, {}, 0

    working_df["match_score"] = existing_match_scores * quality_multiplier
    working_df["_condition_quality_multiplier"] = quality_multiplier
    working_df["_missing_required_condition_fields"] = missing_field_labels
    return working_df, target_reaction, dict(sorted(missing_field_counts.items())), penalized_rows


def _apply_condition_plausibility_penalties(
    matched_df: pd.DataFrame,
    *,
    reaction_type_filter: Optional[str],
    predicted_reaction_type: Optional[str],
    reaction_type_confidence: float,
) -> Tuple[pd.DataFrame, Optional[str], int, int, Dict[str, int]]:
    if matched_df.empty or "Reaction_Type_Standardized" not in matched_df.columns:
        return matched_df, None, 0, 0, {}

    target_reaction = reaction_type_filter
    if (
        not target_reaction
        and predicted_reaction_type
        and predicted_reaction_type != "Unknown"
        and reaction_type_confidence >= 0.5
    ):
        target_reaction = predicted_reaction_type
    target_reaction = _resolve_reaction_type_label(target_reaction)
    if not target_reaction:
        return matched_df, None, 0, 0, {}

    required_catalysts = _required_catalyst_classes_for_reaction_type(target_reaction)
    if not required_catalysts:
        return matched_df, target_reaction, 0, 0, {}

    type_series = matched_df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
    family_mask = type_series.apply(
        lambda value: _resolve_reaction_type_label(value) == target_reaction if value else False
    )
    if not family_mask.any():
        return matched_df, target_reaction, 0, 0, {}

    working_df = matched_df.copy()
    source_series = (
        working_df["Source_Group"].apply(_normalize_source_group)
        if "Source_Group" in working_df.columns
        else pd.Series([""] * len(working_df), index=working_df.index)
    )
    existing_match_scores = (
        pd.to_numeric(working_df["match_score"], errors="coerce").fillna(1.0)
        if "match_score" in working_df.columns
        else pd.Series([1.0] * len(working_df), index=working_df.index, dtype=float)
    )

    filtered_rows = 0
    penalized_rows = 0
    issue_counts: Counter[str] = Counter()
    rows_to_drop: List[Any] = []

    for idx in working_df.index[family_mask]:
        source_group = str(source_series.loc[idx] or "").strip().lower()
        if source_group in _PLAUSIBILITY_SOURCE_EXEMPTIONS:
            continue

        row = working_df.loc[idx]
        if _row_contains_unexpected_catalyst_metals(row, allowed_classes=required_catalysts):
            rows_to_drop.append(idx)
            filtered_rows += 1
            issue_counts["unexpected_catalyst_metals"] += 1
            continue

        if _row_has_implausible_snar_annotation(row):
            penalized_rows += 1
            issue_counts["implausible_snar_annotation"] += 1
            working_df.at[idx, "match_score"] = float(existing_match_scores.loc[idx]) * _PLAUSIBILITY_MECHANISM_PENALTY

    if rows_to_drop:
        working_df = working_df.drop(index=rows_to_drop)

    return working_df, target_reaction, filtered_rows, penalized_rows, dict(sorted(issue_counts.items()))


def _format_source_reaction_ids(
    group_df: pd.DataFrame,
    *,
    max_items: int = 5,
) -> str:
    if group_df is None or group_df.empty:
        return ""

    ids: List[str] = []
    has_source_file = "Source_File" in group_df.columns
    has_source_row = "Source_Row" in group_df.columns
    has_reaction_id = "reaction_id" in group_df.columns

    source_files = group_df["Source_File"].tolist() if has_source_file else [None] * len(group_df)
    source_rows = group_df["Source_Row"].tolist() if has_source_row else [None] * len(group_df)
    for file_val, row_val in zip(source_files, source_rows):
        file_text = str(file_val or "").strip()
        if not file_text:
            continue
        file_name = Path(file_text).name if file_text else ""
        row_text = ""
        if row_val is not None and not (isinstance(row_val, float) and pd.isna(row_val)):
            try:
                row_text = str(int(row_val))
            except (TypeError, ValueError):
                row_text = str(row_val).strip()
        if row_text:
            ids.append(f"{file_name}:{row_text}")
        else:
            ids.append(file_name or file_text)

    ids = _dedupe_list(ids)
    if not ids and has_reaction_id:
        fallback: List[str] = []
        for value in group_df["reaction_id"]:
            text = str(value or "").strip()
            if text and text.lower() != "nan":
                fallback.append(text)
        ids = _dedupe_list(fallback)
    if not ids:
        return ""
    if max_items > 0 and len(ids) > max_items:
        return ", ".join(ids[:max_items]) + f" (+{len(ids) - max_items} more)"
    return ", ".join(ids)


def _reactant_key(values: Iterable[Optional[str]]) -> str:
    items = _dedupe_list([str(v).strip() for v in values if v])
    return "|".join(sorted(items))


def _normalize_source_group(value: Optional[str]) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    label = str(value).strip().lower()
    if not label or label == "nan":
        return ""
    if label in ("literature", "datasets", "dataset", "lit"):
        return "literature"
    if label in ("protocols", "protocol"):
        return "literature"
    if label in ("motif", "motifs", "experiments", "experiment", "experiements"):
        return "motif"
    if label == "rules":
        return "rules"
    return label


def _coerce_rule_tier_value(value: Any) -> float:
    """Normalize rule tier values to numeric (fallback=0, default=1, preferred=2)."""
    if value is None:
        return 0.0
    if isinstance(value, (int, float)) and not pd.isna(value):
        return float(value)
    text = str(value).strip().lower()
    if not text or text == "nan":
        return 0.0
    try:
        return float(text)
    except ValueError:
        pass
    if text in {"fallback", "backup", "low"}:
        return 0.0
    if text in {"default", "normal", "medium"}:
        return 1.0
    if text in {"preferred", "high", "primary"}:
        return 2.0
    return 0.0


def _ensure_rule_tier_column(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure Rule_Tier exists as a numeric column (cache/backward compatible)."""
    if df is None:
        return df
    # Fast path for already-normalized cached frames.
    if "Rule_Tier" in df.columns:
        try:
            if pd.api.types.is_numeric_dtype(df["Rule_Tier"]):
                return df
        except Exception:
            pass
    frame = df.copy()
    if "Rule_Tier" in frame.columns:
        frame["Rule_Tier"] = frame["Rule_Tier"].apply(_coerce_rule_tier_value).astype(float)
        return frame
    if "rule_tier" in frame.columns:
        frame["Rule_Tier"] = frame["rule_tier"].apply(_coerce_rule_tier_value).astype(float)
        return frame
    frame["Rule_Tier"] = 0.0
    return frame


def _filter_source_group(
    df: pd.DataFrame, source_group: Optional[str]
) -> pd.DataFrame:
    if not source_group or "Source_Group" not in df.columns:
        return df
    label = _normalize_source_group(source_group)
    if not label:
        return df
    series = df["Source_Group"].apply(_normalize_source_group)
    return df[series == label]


def _iter_motif_ids(motifs: Iterable[Any]) -> Iterable[str]:
    for motif in motifs:
        if isinstance(motif, dict):
            cid = motif.get("compound_id") or motif.get("id")
            if cid:
                yield str(cid)
            for alt in motif.get("alt_compound_ids") or []:
                if alt:
                    yield str(alt)
        elif motif:
            yield str(motif)


def _is_aryl_halide_motif(cid: str) -> bool:
    if not cid:
        return False
    if cid in _ARYL_HALIDE_MOTIFS:
        return True
    if cid.startswith("Ar-"):
        suffix = cid.split("-", 1)[-1]
        return suffix in _HALIDE_SUFFIXES or suffix == "X"
    return False


def _is_aryl_boronate_motif(cid: str) -> bool:
    if not cid:
        return False
    return cid.startswith("Ar-B")


def _extract_score_values(result: Any) -> List[float]:
    values: List[float] = []
    if isinstance(result, dict):
        score = result.get("score_0_10")
        if isinstance(score, (int, float)):
            values.append(float(score))
    elif isinstance(result, list):
        for entry in result:
            if isinstance(entry, dict):
                score = entry.get("score_0_10")
                if isinstance(score, (int, float)):
                    values.append(float(score))
    return values


def _extract_aryl_scores(analysis: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
    steric_scores: List[float] = []
    electronic_scores: List[float] = []

    for entry in analysis.get("steric", {}).get("aryl", []) or []:
        if isinstance(entry, dict):
            steric_scores.extend(_extract_score_values(entry.get("result")))

    for entry in analysis.get("electronics", {}).get("aryl", []) or []:
        if isinstance(entry, dict):
            electronic_scores.extend(_extract_score_values(entry.get("result")))

    steric = round(sum(steric_scores) / len(steric_scores), 2) if steric_scores else None
    electronic = (
        round(sum(electronic_scores) / len(electronic_scores), 2)
        if electronic_scores
        else None
    )
    return steric, electronic


@lru_cache(maxsize=4096)
def _aryl_role_scores_for_smiles(smiles: str) -> Dict[str, Dict[str, Optional[float]]]:
    if not smiles:
        return {}
    analysis = featurize_molecule(smiles)
    motifs = analysis.get("motifs", []) or []
    roles: List[str] = []
    for cid in _iter_motif_ids(motifs):
        if _is_aryl_halide_motif(cid):
            roles.append("aryl_halide")
        if _is_aryl_boronate_motif(cid):
            roles.append("aryl_boronate")
    if not roles:
        return {}
    steric, electronic = _extract_aryl_scores(analysis)
    if steric is None and electronic is None:
        return {}
    return {
        role: {"steric": steric, "electronic": electronic}
        for role in sorted(set(roles))
    }


def _merge_role_scores(
    role_maps: Iterable[Dict[str, Dict[str, Optional[float]]]]
) -> Dict[str, Dict[str, Optional[float]]]:
    buckets: Dict[str, Dict[str, List[float]]] = {}
    for role_map in role_maps:
        for role, scores in role_map.items():
            bucket = buckets.setdefault(role, {"steric": [], "electronic": []})
            steric = scores.get("steric")
            electronic = scores.get("electronic")
            if isinstance(steric, (int, float)):
                bucket["steric"].append(float(steric))
            if isinstance(electronic, (int, float)):
                bucket["electronic"].append(float(electronic))

    merged: Dict[str, Dict[str, Optional[float]]] = {}
    for role, values in buckets.items():
        steric_values = values["steric"]
        electronic_values = values["electronic"]
        merged[role] = {
            "steric": round(sum(steric_values) / len(steric_values), 2)
            if steric_values
            else None,
            "electronic": round(sum(electronic_values) / len(electronic_values), 2)
            if electronic_values
            else None,
        }
    return merged


@lru_cache(maxsize=2048)
def _aryl_role_scores_for_reaction(reaction_smiles: str) -> Dict[str, Dict[str, Optional[float]]]:
    if not reaction_smiles:
        return {}
    record = ReactionRecord.from_payload(normalize_reaction(reaction_smiles))
    role_maps: List[Dict[str, Dict[str, Optional[float]]]] = []
    for smi in record.reactant_smiles:
        if not smi:
            continue
        role_maps.append(_aryl_role_scores_for_smiles(smi))
    return _merge_role_scores(role_maps)


def _similarity_from_scores(
    query_scores: Dict[str, Optional[float]],
    row_scores: Dict[str, Optional[float]],
) -> Optional[float]:
    parts: List[float] = []
    for key in ("steric", "electronic"):
        q_val = query_scores.get(key)
        r_val = row_scores.get(key)
        if isinstance(q_val, (int, float)) and isinstance(r_val, (int, float)):
            delta = abs(float(q_val) - float(r_val))
            parts.append(max(0.0, min(1.0, 1.0 - (delta / 10.0))))
    if not parts:
        return None
    return sum(parts) / len(parts)


def _aryl_similarity_weight(
    query_roles: Dict[str, Dict[str, Optional[float]]],
    row_roles: Dict[str, Dict[str, Optional[float]]],
) -> float:
    if not query_roles or not row_roles:
        return 1.0
    scores: List[float] = []
    for role, q_scores in query_roles.items():
        r_scores = row_roles.get(role)
        if not r_scores:
            continue
        similarity = _similarity_from_scores(q_scores, r_scores)
        if similarity is None:
            continue
        scores.append(similarity)
    if not scores:
        return 1.0
    combined = sum(scores) / len(scores)
    return 0.7 + (0.3 * combined)


def _iter_smiles_parts(smiles: str) -> List[str]:
    return [part.strip() for part in str(smiles or "").split(".") if part.strip()]


def _build_reaction_smiles(
    reactant_a_smiles: str,
    reactant_b_smiles: Optional[str],
    product_smiles: Optional[str],
) -> str:
    reactants = [reactant_a_smiles]
    if reactant_b_smiles:
        reactants.append(reactant_b_smiles)
    reactant_text = ".".join([part for part in reactants if part])
    if product_smiles:
        return f"{reactant_text}>>{product_smiles}"
    return reactant_text


def _filter_by_reactant_types(
    df: pd.DataFrame,
    reactant_types: Optional[Iterable[str]],
    match_all: bool = False,
) -> pd.DataFrame:
    if not reactant_types:
        return df
    cols = [col for col in ("Reactant_A_Type", "Reactant_B_Type", "Reactant_C_Type") if col in df.columns]
    if not cols:
        return df

    tokens = [str(t).strip() for t in reactant_types if str(t).strip()]
    if not tokens:
        return df

    masks: List[pd.Series] = []
    for token in tokens:
        pattern = re.escape(token)
        col_mask = None
        for col in cols:
            series = df[col].fillna("").astype(str)
            mask = series.str.contains(pattern, case=False, na=False)
            col_mask = mask if col_mask is None else (col_mask | mask)
        if col_mask is not None:
            masks.append(col_mask)

    if not masks:
        return df

    combined = masks[0]
    for mask in masks[1:]:
        combined = combined & mask if match_all else combined | mask
    return df[combined]


def _clean_reactant_value(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if text.lower() == "inorganic":
        return ""
    return _collapse_motif_field(text)


def _normalize_reactant_values(values: Iterable[Any]) -> Tuple[str, str, str, List[str]]:
    cleaned = [_clean_reactant_value(v) for v in values]
    cleaned = [v for v in cleaned if v]
    a_val = cleaned[0] if len(cleaned) > 0 else ""
    b_val = cleaned[1] if len(cleaned) > 1 else ""
    c_val = cleaned[2] if len(cleaned) > 2 else ""
    return a_val, b_val, c_val, cleaned


_MOTIF_SPLIT_RE = re.compile(r"[|,]")
_COMPOUND_GROUPS_FILE = Path(__file__).resolve().parents[1] / "taxonomy" / "data" / "organic_groups.v1.3.json"
_SCAFFOLD_MOTIFS_FILE = Path(__file__).resolve().parents[1] / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"

_MOTIF_TAG_WEIGHTS = {
    "leaving_group": 40,
    "leaving_group_weak": 25,
    "sulfonate_leaving_group": 30,
    "acyl_electrophile": 30,
    "sulfonyl_halide": 30,
    "nucleophile": 20,
    "n_h_nucleophile": 22,
    "o_h_nucleophile": 22,
    "s_h_nucleophile": 22,
}

_SPECTATOR_GROUP_STOPLIST = {
    "Ar",
    "R",
    "Any",
    "Alkyl",
    "Alkenyl",
    "Alkynyl",
    "H",
}
_GENERIC_FALLBACK_MOTIFS = {
    "Ar-R",
    "Ar-H",
    "Alkyl-H",
}
_RULE_NON_CORE_QUERY_MOTIFS = {
    "Ar-R",
    "Ar-H",
    "Alkyl-H",
    "R_acidic-H",
    "Inorganic",
}


@lru_cache(maxsize=1)
def _load_compound_scope_payload() -> Dict[str, Any]:
    try:
        payload = taxonomy_loader.load_organic_compounds()
    except Exception:
        payload = {}
    return payload if isinstance(payload, dict) else {}
_HALIDE_SUFFIXES = {"Cl", "Br", "I", "F"}
_ARYL_HALIDE_MOTIFS = {"Ar-F", "Ar-Cl", "Ar-Br", "Ar-I", "Ar-X"}
_SPECTATOR_MATCH_WEIGHT = 0.7
_INTRAMOLECULAR_MATCH_BOOST = 1.2
_RULE_TIER_BOOST_STEP = 0.05
_RULE_TIER_MAX_BOOST = 1.2
_AROMATIC_SCAFFOLD_FALLBACK = {"Ar", "HeteroAr", "AromN"}
_HETERO_C_H_SCAFFOLD_GROUPS = {"HeteroAr", "AromN"}
_BEST_SELLER_TOP_N = 5

_NH_HETEROCYCLE_TAG = "nh-heteroaromatic"


def _compound_entry_id(entry: Dict[str, Any]) -> str:
    motif_id = str(entry.get("id") or "").strip()
    if motif_id:
        return motif_id
    scaffold = str(entry.get("A") or "").strip()
    substituent = str(entry.get("B") or "").strip()
    if not scaffold or not substituent:
        return ""
    if substituent.startswith("-"):
        substituent = substituent[1:]
    if not substituent:
        return ""
    return f"{scaffold}-{substituent}"


@lru_cache(maxsize=1)
def _load_motif_sets() -> Dict[str, List[str]]:
    return taxonomy_loader.load_compound_logic_sets()


@lru_cache(maxsize=1)
def _load_scope_map() -> Dict[str, List[str]]:
    try:
        from ..taxonomy import loader as taxonomy_loader

        payload = taxonomy_loader.load_motif_scope_index()
        raw_scope = payload.get("scope_map", {}) if isinstance(payload, dict) else {}
        if isinstance(raw_scope, dict):
            scope_map: Dict[str, List[str]] = {}
            for key, value in raw_scope.items():
                parent = str(key).strip()
                if not parent:
                    continue
                children = [
                    str(item).strip()
                    for item in (value or [])
                    if str(item).strip() and str(item).strip() != parent
                ]
                if children:
                    scope_map[parent] = sorted(set(children))
            if scope_map:
                return scope_map
    except Exception:
        pass

    payload = _load_compound_scope_payload()
    if not payload:
        return {}

    compounds = payload.get("compounds") or []
    any_entries: List[Dict[str, str]] = []
    by_b: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        cid = _compound_entry_id(entry)
        a_val = str(entry.get("A") or "").strip()
        b_val = str(entry.get("B") or "").strip()
        template = str(entry.get("template") or "").strip()
        if not cid or not b_val:
            continue
        record = {"id": cid, "A": a_val, "B": b_val, "template": template}
        if a_val == "Any_Scaffold":
            any_entries.append(record)
        else:
            by_b[b_val].append(record)

    scope_map: Dict[str, List[str]] = {}
    for entry in any_entries:
        children: List[str] = []
        for child in by_b.get(entry["B"], []):
            if child["id"] == entry["id"]:
                continue
            if entry["template"] and child["template"] and entry["template"] != child["template"]:
                continue
            children.append(child["id"])
        if children:
            scope_map[entry["id"]] = sorted(set(children))
    return scope_map


@lru_cache(maxsize=1)
def _load_scope_parent_map() -> Dict[str, Set[str]]:
    parent_map: Dict[str, Set[str]] = defaultdict(set)
    for parent, children in _load_scope_map().items():
        p = str(parent).strip()
        if not p:
            continue
        for child in children:
            c = str(child).strip()
            if c and c != p:
                parent_map[c].add(p)
    return dict(parent_map)


@lru_cache(maxsize=1)
def _load_compound_groups_payload() -> Dict[str, Any]:
    if not _COMPOUND_GROUPS_FILE.exists():
        return {}
    try:
        return load_organic_groups_with_compositions(_COMPOUND_GROUPS_FILE)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def _load_group_tags() -> Dict[str, Set[str]]:
    payload = _load_compound_groups_payload()
    if not payload:
        return {}
    groups = payload.get("groups") or []
    tag_map: Dict[str, Set[str]] = {}
    for entry in groups:
        if not isinstance(entry, dict):
            continue
        gid = str(entry.get("id") or "").strip()
        if not gid:
            continue
        tags = {str(tag).strip() for tag in (entry.get("tags") or []) if str(tag).strip()}
        if tags:
            tag_map[gid] = tags
    return tag_map


@lru_cache(maxsize=1)
def _load_scaffold_motif_ids() -> Set[str]:
    if not _SCAFFOLD_MOTIFS_FILE.exists():
        return set()
    try:
        with _SCAFFOLD_MOTIFS_FILE.open("r", encoding="utf-8") as handle:
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
def _load_nh_heterocycle_scaffolds() -> Set[str]:
    if not _SCAFFOLD_MOTIFS_FILE.exists():
        return set()
    try:
        with _SCAFFOLD_MOTIFS_FILE.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return set()
    motifs: Set[str] = set()
    for entry in payload.get("compounds", []) or []:
        if not isinstance(entry, dict):
            continue
        motif_id = str(entry.get("id") or "").strip()
        description = str(entry.get("description") or "").lower()
        if motif_id and _NH_HETEROCYCLE_TAG in description:
            motifs.add(motif_id)
    return motifs


@lru_cache(maxsize=1)
def _load_scaffold_alias_map() -> Dict[str, List[str]]:
    return {cid: ["AromN-H"] for cid in _load_nh_heterocycle_scaffolds()}


@lru_cache(maxsize=1)
def _load_aromatic_scaffold_ids() -> Set[str]:
    payload = _load_compound_groups_payload()
    if not payload:
        return set(_AROMATIC_SCAFFOLD_FALLBACK)

    aromatic_ids: Set[str] = set()
    for entry in payload.get("groups", []) or []:
        if not isinstance(entry, dict):
            continue
        if str(entry.get("kind") or "").strip() != "scaffold":
            continue
        group_id = str(entry.get("id") or "").strip()
        description = str(entry.get("description") or "").lower()
        if group_id and "aromatic" in description:
            aromatic_ids.add(group_id)

    aromatic_ids.update(_AROMATIC_SCAFFOLD_FALLBACK)
    return aromatic_ids


@lru_cache(maxsize=1)
def _load_motif_compatibility_map() -> Dict[str, Set[str]]:
    payload = _load_compound_scope_payload()
    if not payload:
        return {}

    aromatic_scaffolds = _load_aromatic_scaffold_ids()
    by_substituent: Dict[str, Set[str]] = defaultdict(set)
    for entry in payload.get("compounds", []) or []:
        if not isinstance(entry, dict):
            continue
        motif_id = _compound_entry_id(entry)
        scaffold = str(entry.get("A") or "").strip()
        substituent = str(entry.get("B") or "").strip()
        if not motif_id or not scaffold or not substituent:
            continue
        if scaffold in aromatic_scaffolds:
            by_substituent[substituent].add(motif_id)

    compatibility: Dict[str, Set[str]] = {}
    for members in by_substituent.values():
        if len(members) < 2:
            continue
        for motif_id in members:
            others = members - {motif_id}
            if others:
                compatibility[motif_id] = set(others)
    return compatibility


@lru_cache(maxsize=1)
def _load_compound_tags() -> Dict[str, Set[str]]:
    payload = _load_compound_scope_payload()
    if not payload:
        return {}
    group_tags = _load_group_tags()
    compounds = payload.get("compounds") or []
    tag_map: Dict[str, Set[str]] = {}
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        cid = _compound_entry_id(entry)
        if not cid:
            continue
        tags: Set[str] = set()
        group_a = str(entry.get("A") or "").strip()
        group_b = str(entry.get("B") or "").strip()
        tags.update(group_tags.get(group_a, set()))
        tags.update(group_tags.get(group_b, set()))
        if tags:
            tag_map[cid] = tags
    return tag_map


@lru_cache(maxsize=1)
def _load_compound_ids() -> Set[str]:
    payload = _load_compound_scope_payload()
    if not payload:
        return set()
    compounds = payload.get("compounds") or []
    ids = set()
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        cid = _compound_entry_id(entry)
        if cid:
            ids.add(cid)
    return ids


def _motif_tag_score(motif_id: str) -> int:
    tags = _load_compound_tags().get(motif_id)
    if not tags:
        return 0
    score = 0
    for tag, weight in _MOTIF_TAG_WEIGHTS.items():
        if tag in tags:
            score += weight
    return score


def _filter_generic_motifs(motifs: Iterable[str]) -> List[str]:
    return [m for m in motifs if m and m not in _GENERIC_FALLBACK_MOTIFS]


def _expand_parent_motifs(motifs: Iterable[str]) -> List[str]:
    compound_ids = _load_compound_ids()
    compatibility = _load_motif_compatibility_map()
    expanded = list(motifs)
    alias_map = _load_scaffold_alias_map()
    for motif in motifs:
        text = str(motif).strip()
        aliases = alias_map.get(text) or []
        if aliases:
            expanded.extend(aliases)
        compat_members = compatibility.get(text) or set()
        if compat_members:
            expanded.extend(sorted(compat_members))
        if "-" not in text:
            continue
        prefix, suffix = text.rsplit("-", 1)
        if suffix in _HALIDE_SUFFIXES:
            candidate = f"{prefix}-X"
            if candidate in compound_ids:
                expanded.append(candidate)
            # Expand aromatic halides across leaving groups to enable chemistry-aware
            # fallback (e.g., Ar-Cl query can borrow Ar-Br literature precedent).
            if prefix in _AROMATIC_SCAFFOLD_FALLBACK:
                expanded.extend(f"{prefix}-{hal}" for hal in sorted(_HALIDE_SUFFIXES))
                expanded.append(candidate)
        elif suffix == "X" and prefix in _AROMATIC_SCAFFOLD_FALLBACK:
            expanded.extend(f"{prefix}-{hal}" for hal in sorted(_HALIDE_SUFFIXES))
    return _dedupe_list(expanded)


def _build_fallback_motif_list(
    motifs: Iterable[str],
    reacted: Optional[Set[str]],
    spectators: Optional[Set[str]],
    allow_generic: bool = False,
) -> List[str]:
    base = _dedupe_list([m for m in motifs if m])
    if base and not allow_generic:
        filtered = _filter_generic_motifs(base)
        if filtered:
            base = filtered
        else:
            base = []
    expanded = _expand_parent_motifs(base) if base else []
    prioritized = _prioritize_motifs(expanded or base, reacted, spectators)
    return prioritized or [""]


def _build_fallback_tiers(
    motifs: Iterable[str],
    reacted: Optional[Set[str]],
    spectators: Optional[Set[str]],
) -> List[List[str]]:
    base = _dedupe_list([m for m in motifs if m])
    if not base:
        return [[""]]
    tier_specific = _build_fallback_motif_list(base, reacted, spectators, allow_generic=False)
    tier_generic = _build_fallback_motif_list(base, reacted, spectators, allow_generic=True)
    tiers: List[List[str]] = []
    if tier_specific == [""] and tier_generic != [""]:
        tiers.extend([tier_generic, tier_specific])
    else:
        tiers.append(tier_specific)
        if tier_generic != tier_specific:
            tiers.append(tier_generic)
    deduped: List[List[str]] = []
    seen: Set[Tuple[str, ...]] = set()
    for tier in tiers:
        key = tuple(tier)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(tier)
    return deduped or [[""]]


def _split_motif_tokens(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, float) and pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    return [
        token
        for token in (
            normalize_motif_id(raw_token.strip()) for raw_token in _MOTIF_SPLIT_RE.split(text)
        )
        if token
    ]


@lru_cache(maxsize=4096)
def _expanded_match_tokens(token: str) -> Set[str]:
    text = str(token).strip()
    if not text:
        return set()
    expanded: Set[str] = {text}
    alias_map = _load_scaffold_alias_map()
    expanded.update(alias_map.get(text) or [])
    expanded.update(_load_motif_compatibility_map().get(text) or set())
    expanded.update(_load_scope_map().get(text) or [])
    expanded.update(_load_scope_parent_map().get(text) or set())

    # Treat boronic acid (-B(OH)2) and boronate ester (-B(OR)2) as equivalent
    # coupling partners — the dataset stores both under the B(OR)2 label.
    for member in list(expanded):
        if member.endswith("-B(OH)2"):
            expanded.add(member[: -len("-B(OH)2")] + "-B(OR)2")
        elif member.endswith("-B(OR)2"):
            expanded.add(member[: -len("-B(OR)2")] + "-B(OH)2")

    compound_ids = _load_compound_ids()
    for member in list(expanded):
        if "-" not in member:
            continue
        prefix, suffix = member.rsplit("-", 1)
        if suffix in _HALIDE_SUFFIXES:
            generic = f"{prefix}-X"
            if generic in compound_ids:
                expanded.add(generic)
    return expanded


def _motif_tokens_compatible(left: str, right: str) -> bool:
    left_text = str(left).strip()
    right_text = str(right).strip()
    if not left_text or not right_text:
        return False
    if left_text == right_text:
        return True
    left_expanded = _expanded_match_tokens(left_text)
    right_expanded = _expanded_match_tokens(right_text)
    if not left_expanded or not right_expanded:
        return False
    if right_text in left_expanded or left_text in right_expanded:
        return True
    return bool(left_expanded & right_expanded)


def _token_set_subset_compatible(required: Set[str], available: Set[str]) -> bool:
    if not required:
        return True
    if not available:
        return False
    for req in required:
        if not any(_motif_tokens_compatible(req, cand) for cand in available):
            return False
    return True


def _token_set_overlap_compatible(left: Set[str], right: Set[str]) -> int:
    if not left or not right:
        return 0
    used_right: Set[str] = set()
    overlap = 0
    for left_token in sorted(left):
        for right_token in sorted(right):
            if right_token in used_right:
                continue
            if _motif_tokens_compatible(left_token, right_token):
                used_right.add(right_token)
                overlap += 1
                break
    return overlap


def _signature_matches(query_sets: Iterable[Set[str]], value: Any) -> bool:
    tokens = set(_split_motif_tokens(value))
    if not tokens:
        return False
    for query_tokens in query_sets:
        if not query_tokens:
            continue
        if _token_set_subset_compatible(tokens, query_tokens) or _token_set_subset_compatible(
            query_tokens, tokens
        ):
            return True
    return False


def _select_rule_required_core_tokens(query_tokens: Set[str]) -> Set[str]:
    """Keep only essential motifs for rule-path enforcement."""
    cleaned = {token for token in query_tokens if token}
    if not cleaned:
        return set()

    required = {token for token in cleaned if token not in _RULE_NON_CORE_QUERY_MOTIFS}
    if len(required) >= 2:
        return required

    informative_required = {token for token in required if _motif_tag_score(token) > 0}
    if informative_required:
        return informative_required
    if required:
        return required

    informative_all = {token for token in cleaned if _motif_tag_score(token) > 0}
    return informative_all or cleaned


def _build_rule_required_core_query_sets(query_sets: Iterable[Set[str]]) -> List[Set[str]]:
    narrowed: List[Set[str]] = []
    seen: Set[Tuple[str, ...]] = set()
    for query_tokens in query_sets:
        selected = _select_rule_required_core_tokens(set(query_tokens))
        if not selected:
            continue
        key = tuple(sorted(selected))
        if key in seen:
            continue
        seen.add(key)
        narrowed.append(selected)
    return narrowed


def _encode_signature(tokens: Iterable[str]) -> str:
    """Encode motif tokens into a stable, sorted signature string."""
    if not tokens:
        return ""
    cleaned = {token for token in tokens if token}
    if not cleaned:
        return ""
    return "|".join(sorted(cleaned))


def _reactant_types_to_signature(values: Iterable[Any]) -> str:
    """Build a reactant signature from reactant type fields."""
    tokens: List[str] = []
    for value in values:
        tokens.extend(_split_motif_tokens(value))
    if not tokens:
        return ""
    expanded = _expand_motif_tokens(tokens, _load_motif_sets(), _load_scope_map())
    return _encode_signature(expanded or tokens)


def _reaction_key_to_signatures(key: str) -> Tuple[str, str]:
    """Project a CRK Reaction_Key into core/ext reactant signatures."""
    if not key:
        return "", ""
    reacted, _, spectators = _parse_transformation_key(key)
    core = _encode_signature(reacted)
    ext = _encode_signature(set(reacted) | set(spectators))
    return core, ext


@lru_cache(maxsize=4096)
def _signature_lookup_candidates(signature: str) -> Tuple[str, ...]:
    """
    Build lookup candidates for a reactant signature using scope-parent fallback.

    Example:
    - query token `HeteroAr-B(OH)2` can produce parent `Ar-B(OH)2`
    - query token `HeteroAr-Br` can produce parent `Ar-Br`
    """
    tokens = sorted(set(_split_motif_tokens(signature)))
    if not tokens:
        return tuple()

    parent_map = _load_scope_parent_map()
    options: List[List[str]] = []
    for token in tokens:
        parents = sorted(parent_map.get(token) or set())
        expanded = [token]
        expanded.extend(parent for parent in parents if parent and parent != token)
        options.append(_dedupe_list(expanded))

    candidates: List[str] = []
    seen: Set[str] = set()
    for combo in itertools.product(*options):
        candidate = _encode_signature(combo)
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        candidates.append(candidate)
        # Keep this bounded to avoid combinatorial blow-up on noisy signatures.
        if len(candidates) >= 64:
            break

    return tuple(candidates)


@lru_cache(maxsize=1024)
def _featurize_reaction_for_recommendation(
    reaction_smiles: str,
    *,
    skip_bond_analysis: bool,
) -> Dict[str, Any]:
    payload = featurize_reaction(
        reaction_smiles,
        options={"skip_bond_analysis": bool(skip_bond_analysis)},
    )
    return payload if isinstance(payload, dict) else {}


def _intramolecular_likely_from_fields(
    reactant_a: Any,
    reactant_b: Any,
    reactant_c: Any,
) -> bool:
    tokens_a = _split_motif_tokens(reactant_a)
    tokens_b = _split_motif_tokens(reactant_b)
    tokens_c = _split_motif_tokens(reactant_c)
    reactant_tokens = [tokens for tokens in (tokens_a, tokens_b, tokens_c) if tokens]
    if len(reactant_tokens) != 1:
        return False
    return len(reactant_tokens[0]) > 1


def _apply_intramolecular_boost(
    df: pd.DataFrame,
    query_intramolecular_likely: bool,
) -> None:
    if not query_intramolecular_likely or df is None or df.empty:
        return
    if "Intramolecular_Likely" in df.columns:
        col = "Intramolecular_Likely"
    elif "Is_Intramolecular" in df.columns:
        col = "Is_Intramolecular"
    else:
        return
    series = df[col]
    # Avoid fillna() downcast warnings on object dtype while keeping NaN -> False behavior.
    mask = series.notna() & series.astype(bool)
    if "match_score" in df.columns:
        df.loc[mask, "match_score"] = df.loc[mask, "match_score"] * _INTRAMOLECULAR_MATCH_BOOST


def _collapse_motif_tokens(tokens: Iterable[str]) -> str:
    cleaned = _dedupe_list([str(token).strip() for token in tokens if str(token).strip()])
    if not cleaned:
        return ""
    return ",".join(sorted(cleaned))


def _collapse_motif_field(value: Any) -> str:
    return _collapse_motif_tokens(_split_motif_tokens(value))


def _split_group_tokens(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, float) and pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    return [token.strip() for token in re.split(r"[|/,]", text) if token.strip()]


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


def _spectator_groups_from_motifs(motif_ids: Iterable[str]) -> Set[str]:
    groups: Set[str] = set()
    scaffold_ids = _load_scaffold_motif_ids()
    for motif_id in motif_ids:
        text = str(motif_id).strip()
        if not text:
            continue
        if text in scaffold_ids:
            groups.add(text)
            continue
        group_id = _group_id_from_motif_id(text)
        if group_id and group_id not in _SPECTATOR_GROUP_STOPLIST:
            groups.add(group_id)
    return groups


def _spectator_similarity(query_groups: Set[str], row_groups: Set[str]) -> float:
    return weighted_spectator_similarity(query_groups, row_groups)


def _expand_macro_token(
    token: str,
    motif_sets: Dict[str, List[str]],
    scope_map: Dict[str, List[str]],
) -> List[str]:
    token = token.strip()
    if token.startswith("@"):
        set_name = token[1:]
        members = motif_sets.get(set_name) or []
        if members:
            return members
    expanded = [token]
    if token in scope_map:
        expanded = [token] + scope_map[token]
    alias_map = _load_scaffold_alias_map()
    aliases = alias_map.get(token) or []
    if aliases:
        expanded = _dedupe_list(expanded + aliases)
    return expanded


def _expand_motif_tokens(
    tokens: Iterable[str],
    motif_sets: Dict[str, List[str]],
    scope_map: Dict[str, List[str]],
) -> List[str]:
    expanded: List[str] = []
    for token in tokens:
        expanded.extend(_expand_macro_token(token, motif_sets, scope_map))
    return _dedupe_list(expanded)


def _expand_reactant_field(
    value: Any,
    motif_sets: Dict[str, List[str]],
    scope_map: Dict[str, List[str]],
) -> List[str]:
    tokens = _split_motif_tokens(value)
    if not tokens:
        return [""]
    options = [_expand_macro_token(token, motif_sets, scope_map) for token in tokens]
    expanded_values: List[str] = []
    for combo in itertools.product(*options):
        cleaned = _dedupe_list([item.strip() for item in combo if item and str(item).strip()])
        expanded_values.append(",".join(cleaned))
    return _dedupe_list([item for item in expanded_values if item or item == ""])


def _expand_reactant_keys(
    reactant_a: Any,
    reactant_b: Any,
    motif_sets: Dict[str, List[str]],
    scope_map: Dict[str, List[str]],
) -> List[str]:
    expanded_a = _expand_reactant_field(reactant_a, motif_sets, scope_map)
    expanded_b = _expand_reactant_field(reactant_b, motif_sets, scope_map)
    keys: List[str] = []
    for a_value in expanded_a:
        for b_value in expanded_b:
            key = _reactant_key([a_value, b_value])
            if key:
                keys.append(key)
    return _dedupe_list(keys)


def _parse_transformation_key(key: str) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Parse CRK-v1 format: CRK-v1 |Reactants -> Product_Broad | products_broad: ... | spectators: ...
    
    Returns:
        (reacted_set, formed_set, spectators_set)
    """
    if not key:
        return set(), set(), set()

    text = str(key).strip()
    if "->" not in text:
        return set(), set(), set()

    sections = [section.strip() for section in text.split(" | ") if section.strip()]
    summary = sections[0].strip()
    if summary.startswith("CRK-v1"):
        summary = summary[len("CRK-v1"):].strip()
    if summary.startswith("|"):
        summary = summary[1:].strip()

    def _clean(tokens: Iterable[str]) -> List[str]:
        return [t for t in tokens if t and t != "[]"]

    reacted: Set[str] = set()
    formed: Set[str] = set()
    if "->" in summary:
        reactant_part, product_part = summary.split("->", 1)
        reactant_part = reactant_part.strip()
        product_part = product_part.strip()
        if "(" in product_part:
            product_part = product_part.split("(", 1)[0].strip()
        tokens = _clean(_split_motif_tokens(reactant_part))
        reacted = set(_expand_motif_tokens(tokens, _load_motif_sets(), _load_scope_map()))
        if product_part and product_part != "[]":
            tokens = _clean(_split_motif_tokens(product_part))
            formed = set(_expand_motif_tokens(tokens, _load_motif_sets(), _load_scope_map()))
    spectators: Set[str] = set()
    for section in sections[1:]:
        lower = section.lower()
        if lower.startswith("products_reactive:"):
            if formed:
                continue
            payload = section.split(":", 1)[1].strip()
            tokens = _clean(_split_motif_tokens(payload))
            formed = set(_expand_motif_tokens(tokens, _load_motif_sets(), _load_scope_map()))
        elif lower.startswith("products_broad:"):
            if formed:
                continue
            payload = section.split(":", 1)[1].strip()
            tokens = _clean(_split_motif_tokens(payload))
            formed = set(_expand_motif_tokens(tokens, _load_motif_sets(), _load_scope_map()))
        elif lower.startswith("spectators:"):
            payload = section.split(":", 1)[1].strip()
            tokens = _clean(_split_motif_tokens(payload))
            spectators = set(_expand_motif_tokens(tokens, _load_motif_sets(), _load_scope_map()))

    return reacted, formed, spectators


def _normalize_bond_token_text(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    return text.replace(" ", "")


def _canonicalize_bond_token_for_match(value: str) -> str:
    token = _normalize_bond_token_text(value)
    if not token or "-" not in token:
        return token
    left, right = token.split("-", 1)

    def _elem(part: str) -> str:
        text = str(part or "").strip()
        if "(" in text:
            text = text.split("(", 1)[0].strip()
        return text

    left_el = _elem(left)
    right_el = _elem(right)
    if not left_el or not right_el:
        return token
    return "-".join(sorted((left_el, right_el)))


def _event_signature_from_kinds(kinds: Iterable[str]) -> str:
    kind_set = {str(kind).strip() for kind in (kinds or []) if str(kind).strip()}
    if not kind_set:
        return ""
    ordered_codes: List[str] = []
    for kind in _EVENT_SIGNATURE_PRIORITY:
        if kind not in kind_set:
            continue
        code = _EVENT_SIGNATURE_CODE.get(kind)
        if code and code not in ordered_codes:
            ordered_codes.append(code)
    if not ordered_codes:
        return ""
    return "+".join(ordered_codes[:4])


def _reaction_events_to_match_key(value: Any) -> str:
    """
    Canonical match key from Reaction_Events payload/text.

    Supports:
    - Compact text format (sig/form/break/redox...)
    - Legacy JSON text payload
    - Featurizer-native dict payload (`events`, `bond_pairs`, `redox_assessment`)
    """
    payload: Dict[str, Any] = {}
    if isinstance(value, dict):
        payload = dict(value)
    else:
        payload = deserialize_reaction_events_text(value)
        if not payload:
            return ""

    token_chunks: List[str] = []

    event_signature = str(payload.get("event_signature") or "").strip()
    event_kinds = payload.get("event_kinds")
    if not event_kinds and isinstance(payload.get("events"), list):
        event_kinds = [
            str(event.get("kind") or "").strip()
            for event in payload.get("events") or []
            if isinstance(event, dict) and str(event.get("kind") or "").strip()
        ]
    if not event_signature and isinstance(event_kinds, (list, tuple, set)):
        event_signature = _event_signature_from_kinds(event_kinds)
    if event_signature:
        token_chunks.append("sig=" + event_signature)
    if isinstance(event_kinds, (list, tuple, set)):
        cleaned_kinds = sorted({str(kind).strip() for kind in event_kinds if str(kind).strip()})
        if cleaned_kinds and not event_signature:
            token_chunks.append("kinds=" + "+".join(cleaned_kinds))
    event_families = payload.get("event_families")
    if isinstance(event_families, (list, tuple, set)):
        cleaned_families = sorted(
            {str(family).strip() for family in event_families if str(family).strip()}
        )
        if cleaned_families and not event_signature and not token_chunks:
            token_chunks.append("fam=" + "+".join(cleaned_families))
    mechanism_shortlist = payload.get("mechanism_shortlist")
    if isinstance(mechanism_shortlist, (list, tuple, set)):
        cleaned_mechanisms = sorted(
            {str(name).strip() for name in mechanism_shortlist if str(name).strip()}
        )
        if cleaned_mechanisms and not event_signature and not token_chunks:
            token_chunks.append("mech=" + "+".join(cleaned_mechanisms))

    formed_tokens: List[str] = []
    raw_formed = payload.get("bond_formed")
    if isinstance(raw_formed, (list, tuple, set)):
        formed_tokens.extend(_canonicalize_bond_token_for_match(str(token)) for token in raw_formed)
    bond_pairs = payload.get("bond_pairs")
    if isinstance(bond_pairs, dict):
        for pair in bond_pairs.get("formed") or []:
            if isinstance(pair, (list, tuple)) and len(pair) == 2:
                left = str(pair[0]).strip()
                right = str(pair[1]).strip()
                if left and right:
                    formed_tokens.append("-".join(sorted((left, right))))
    formed_tokens = sorted({token for token in formed_tokens if token})
    if not formed_tokens:
        formed_classes = payload.get("formed_bond_classes")
        if isinstance(formed_classes, (list, tuple, set)):
            formed_tokens = sorted(
                {
                    _canonicalize_bond_token_for_match(str(token))
                    for token in formed_classes
                    if str(token).strip()
                }
            )
            formed_tokens = [token for token in formed_tokens if token]
    if formed_tokens:
        token_chunks.append("form=" + ";".join(formed_tokens))

    broken_tokens: List[str] = []
    raw_broken = payload.get("bond_broken")
    if isinstance(raw_broken, (list, tuple, set)):
        broken_tokens.extend(_canonicalize_bond_token_for_match(str(token)) for token in raw_broken)
    if isinstance(bond_pairs, dict):
        for pair in bond_pairs.get("broken") or []:
            if isinstance(pair, (list, tuple)) and len(pair) == 2:
                left = str(pair[0]).strip()
                right = str(pair[1]).strip()
                if left and right:
                    broken_tokens.append("-".join(sorted((left, right))))
    broken_tokens = sorted({token for token in broken_tokens if token})
    if not broken_tokens:
        broken_classes = payload.get("broken_bond_classes")
        if isinstance(broken_classes, (list, tuple, set)):
            broken_tokens = sorted(
                {
                    _canonicalize_bond_token_for_match(str(token))
                    for token in broken_classes
                    if str(token).strip()
                }
            )
            broken_tokens = [token for token in broken_tokens if token]
    if broken_tokens:
        token_chunks.append("break=" + ";".join(broken_tokens))

    redox = str(payload.get("redox_classification") or "").strip()
    if not redox and isinstance(payload.get("redox_assessment"), dict):
        redox = str((payload.get("redox_assessment") or {}).get("classification") or "").strip()
    if redox:
        token_chunks.append("redox=" + redox)

    # Last-resort fallback when structured event details are sparse.
    if not token_chunks:
        ehyb = str(payload.get("electrophile_hybridization") or "").strip()
        if ehyb:
            token_chunks.append("ehyb=" + ehyb)

    if not token_chunks:
        return ""
    return _EVENT_MATCH_PREFIX + "|".join(token_chunks)


def _reaction_event_match_parts(match_key: Any) -> Dict[str, Any]:
    text = str(match_key or "").strip()
    if not text:
        return {
            "signature": set(),
            "formed": set(),
            "broken": set(),
            "event_kinds": set(),
            "redox": "",
            "reacted_context": set(),
            "formed_context": set(),
        }
    if text.startswith(_EVENT_MATCH_PREFIX):
        text = text[len(_EVENT_MATCH_PREFIX) :]
    chunks = [chunk.strip() for chunk in text.split("|") if chunk.strip()]
    out: Dict[str, Any] = {
        "signature": set(),
        "formed": set(),
        "broken": set(),
        "event_kinds": set(),
        "redox": "",
        "reacted_context": set(),
        "formed_context": set(),
    }
    for chunk in chunks:
        if "=" not in chunk:
            continue
        key, value = chunk.split("=", 1)
        key = key.strip().lower()
        value = value.strip()
        if not value:
            continue
        if key == "sig":
            out["signature"] = {token for token in value.split("+") if token}
        elif key == "kinds":
            out["event_kinds"] = {token for token in value.split("+") if token}
        elif key == "form":
            out["formed"] = {token for token in value.split(";") if token}
        elif key == "break":
            out["broken"] = {token for token in value.split(";") if token}
        elif key == "redox":
            out["redox"] = value
        elif key == "ctx_reacted":
            out["reacted_context"] = {token for token in value.split(";") if token}
        elif key == "ctx_formed":
            out["formed_context"] = {token for token in value.split(";") if token}
    return out


def _jaccard_similarity(left: Set[str], right: Set[str]) -> float:
    if not left and not right:
        return 1.0
    if not left or not right:
        return 0.0
    union = left | right
    if not union:
        return 0.0
    return len(left & right) / len(union)


def _reaction_event_similarity(query_match_key: Any, row_match_key: Any) -> Optional[float]:
    query_parts = _reaction_event_match_parts(query_match_key)
    row_parts = _reaction_event_match_parts(row_match_key)
    total_weight = 0.0
    weighted_score = 0.0

    query_formed = set(query_parts.get("formed") or set())
    if query_formed:
        weight = float(_REACTION_EVENT_COMPONENT_WEIGHTS["formed"])
        total_weight += weight
        weighted_score += weight * _jaccard_similarity(query_formed, set(row_parts.get("formed") or set()))

    query_broken = set(query_parts.get("broken") or set())
    if query_broken:
        weight = float(_REACTION_EVENT_COMPONENT_WEIGHTS["broken"])
        total_weight += weight
        weighted_score += weight * _jaccard_similarity(query_broken, set(row_parts.get("broken") or set()))

    query_signature = set(query_parts.get("signature") or set())
    if query_signature:
        weight = float(_REACTION_EVENT_COMPONENT_WEIGHTS["signature"])
        total_weight += weight
        weighted_score += weight * _jaccard_similarity(
            query_signature,
            set(row_parts.get("signature") or set()),
        )

    query_event_kinds = set(query_parts.get("event_kinds") or set())
    if query_event_kinds:
        weight = float(_REACTION_EVENT_COMPONENT_WEIGHTS["event_kinds"])
        total_weight += weight
        weighted_score += weight * _jaccard_similarity(
            query_event_kinds,
            set(row_parts.get("event_kinds") or set()),
        )

    query_redox = str(query_parts.get("redox") or "").strip()
    if query_redox:
        weight = float(_REACTION_EVENT_COMPONENT_WEIGHTS["redox"])
        total_weight += weight
        row_redox = str(row_parts.get("redox") or "").strip()
        weighted_score += weight * (1.0 if row_redox and row_redox == query_redox else 0.0)

    query_reacted_context = set(query_parts.get("reacted_context") or set())
    if query_reacted_context:
        weight = float(_REACTION_EVENT_COMPONENT_WEIGHTS["reacted_context"])
        total_weight += weight
        weighted_score += weight * _jaccard_similarity(
            query_reacted_context,
            set(row_parts.get("reacted_context") or set()),
        )

    query_formed_context = set(query_parts.get("formed_context") or set())
    if query_formed_context:
        weight = float(_REACTION_EVENT_COMPONENT_WEIGHTS["formed_context"])
        total_weight += weight
        weighted_score += weight * _jaccard_similarity(
            query_formed_context,
            set(row_parts.get("formed_context") or set()),
        )

    if total_weight <= 0.0:
        return None
    return max(0.0, min(1.0, weighted_score / total_weight))


@lru_cache(maxsize=4096)
def _reaction_event_key_from_reaction_smiles(reaction_smiles: str) -> str:
    text = str(reaction_smiles or "").strip()
    if not text or ">" not in text:
        return ""
    try:
        rxn_features = _featurize_reaction_for_recommendation(
            text,
            skip_bond_analysis=True,
        )
    except Exception:
        return ""
    if not isinstance(rxn_features, dict):
        return ""
    reaction_data = rxn_features.get("reaction")
    if not isinstance(reaction_data, dict):
        reaction_data = rxn_features
    if not isinstance(reaction_data, dict):
        return ""
    query_events_payload = (
        reaction_data.get("reaction_events")
        if isinstance(reaction_data.get("reaction_events"), dict)
        else None
    )
    match_key = _reaction_events_to_match_key(query_events_payload or {})
    if not match_key:
        match_key = _reaction_events_to_match_key(
            build_reaction_events_payload(
                reaction_data.get("reaction_key") or "",
                query_events_payload,
            )
        )
    if not match_key:
        return ""

    base_text = match_key[len(_EVENT_MATCH_PREFIX) :] if match_key.startswith(_EVENT_MATCH_PREFIX) else match_key
    chunks = [chunk for chunk in base_text.split("|") if chunk]

    aggregates = reaction_data.get("aggregates") if isinstance(reaction_data, dict) else None
    if isinstance(aggregates, dict):
        reacted_context = sorted(
            {
                normalize_motif_id(str(token).strip())
                for token in (aggregates.get("reacted_motifs") or [])
                if str(token).strip()
            }
        )
        formed_context = sorted(
            {
                normalize_motif_id(str(token).strip())
                for token in (aggregates.get("formed_motifs") or [])
                if str(token).strip()
            }
        )
        reacted_context = [token for token in reacted_context if token]
        formed_context = [token for token in formed_context if token]
        if reacted_context:
            chunks.append("ctx_reacted=" + ";".join(reacted_context))
        if formed_context:
            chunks.append("ctx_formed=" + ";".join(formed_context))

    if not chunks:
        return ""
    return _EVENT_MATCH_PREFIX + "|".join(chunks)


def _derive_query_sets(
    reactant_motifs: Set[str],
    product_motifs: Set[str],
) -> Tuple[Set[str], Set[str], Set[str]]:
    """Derive reacted/formed/spectator motifs from reactant/product sets."""
    reacted = reactant_motifs - product_motifs
    formed = product_motifs - reactant_motifs
    spectators = reactant_motifs & product_motifs
    return reacted, formed, spectators


def _prioritize_motifs(
    motifs: Iterable[str],
    reacted: Optional[Set[str]],
    spectators: Optional[Set[str]],
) -> List[str]:
    ordered = [m for m in motifs if m]
    if not ordered:
        return ordered
    index_map = {m: idx for idx, m in enumerate(ordered)}
    scores = {m: _motif_tag_score(m) for m in ordered}
    def _sort_key(motif: str) -> Tuple[int, int]:
        return (-scores.get(motif, 0), index_map.get(motif, 0))
    if not reacted and not spectators:
        return sorted(ordered, key=_sort_key)
    reacted = reacted or set()
    spectators = spectators or set()
    buckets = {0: [], 1: [], 2: []}
    for motif in ordered:
        if motif in reacted:
            buckets[0].append(motif)
        elif motif in spectators:
            buckets[2].append(motif)
        else:
            buckets[1].append(motif)
    return (
        sorted(buckets[0], key=_sort_key)
        + sorted(buckets[1], key=_sort_key)
        + sorted(buckets[2], key=_sort_key)
    )


@dataclass
class ConditionRecommendation:
    """Single condition recommendation with metadata
    
    Recommendations are ranked primarily by avg_z_score, which measures
    the success of a condition relative to all experiments in the database.
    """
    catalyst: str
    ligand: str
    base: str
    solvent: str
    secondary_solvent: Optional[str] = None
    additive: Optional[str] = None
    coupling_reagent: Optional[str] = None
    temperature: Optional[float] = None
    atmosphere: Optional[str] = None
    spectator_groups: Optional[str] = None
    spectator_score: float = 0.0
    
    # Statistics
    success_rate: float = 0.0  # % of experiments with yield > 50
    avg_yield: float = 0.0
    median_yield: float = 0.0
    num_experiments: int = 0
    avg_z_score: float = 0.0  # Average z-score (PRIMARY ranking metric for condition success)
    confidence_score: float = 0.0  # Secondary score considering z-score and sample size
    match_score: float = 1.0  # How well the transformation matched the query
    condition_quality_score: float = 1.0
    missing_required_fields: Tuple[str, ...] = ()
    
    # Metadata
    reaction_type: Optional[str] = None
    reaction_category: Optional[str] = None
    reaction_id: Optional[str] = None
    reaction_key: Optional[str] = None
    reaction_events: Optional[str] = None
    reactant_types: Tuple[str, str] = ("", "")
    z_score_range: Tuple[float, float] = (0.0, 0.0)


@dataclass
class HTERecommendationResult:
    """Complete recommendation result with ranked conditions"""
    reactant_a_smiles: str
    reactant_b_smiles: Optional[str]
    product_smiles: Optional[str] = None
    
    # Detected types
    reactant_a_type: Optional[str] = None
    reactant_b_type: Optional[str] = None
    reactant_a_category: Optional[str] = None
    reactant_b_category: Optional[str] = None
    product_type: Optional[str] = None
    
    # Predicted reaction type
    predicted_reaction_type: Optional[str] = None
    reaction_type_confidence: float = 0.0

    # Pre-eval motif splits (when product is provided)
    reacted_motifs: Optional[Tuple[str, ...]] = None
    formed_motifs: Optional[Tuple[str, ...]] = None
    spectator_motifs: Optional[Tuple[str, ...]] = None
    query_spectator_groups: Optional[Tuple[str, ...]] = None
    spectator_scoring_applied: bool = False
    spectator_rows_with_groups: int = 0
    spectator_rows_total: int = 0
    spectator_similarity_avg: float = 0.0
    spectator_similarity_range: Tuple[float, float] = (0.0, 0.0)
    query_reaction_key: Optional[str] = None  # Formatted Reaction_Key for the query
    query_reaction_events_key: Optional[str] = None
    
    # Recommendations
    recommendations: List[ConditionRecommendation] = field(default_factory=list)
    recommendations_by_source: Dict[str, List[ConditionRecommendation]] = field(default_factory=dict)
    
    # Metadata
    total_matching_experiments: int = 0
    database_coverage: float = 0.0  # % of database that matches this query
    is_fallback_match: bool = False
    is_filtered_by_detected_type: bool = False  # Whether results were filtered by detected reaction type
    catalyst_requirement_enforced: bool = False
    required_catalyst_family: Optional[str] = None
    required_catalyst_classes: Tuple[str, ...] = ()
    filtered_missing_catalyst_rows: int = 0
    retained_missing_catalyst_rows: int = 0
    condition_quality_family: Optional[str] = None
    penalized_incomplete_condition_rows: int = 0
    missing_required_condition_fields: Dict[str, int] = field(default_factory=dict)
    plausibility_family: Optional[str] = None
    filtered_implausible_catalyst_rows: int = 0
    penalized_implausible_mechanism_rows: int = 0
    plausibility_issue_counts: Dict[str, int] = field(default_factory=dict)
    matched_motifs: Optional[Tuple[str, str]] = None
    timing_ms: Dict[str, float] = field(default_factory=dict)


def _run_precedent_knn(
    reactant_a_smiles: str,
    reactant_b_smiles: Optional[str],
    product_smiles: Optional[str],
    reaction_type: Optional[str],
    top_k: int,
    source_group: Optional[str] = None,
    *,
    prefer_mixfp_for_similarity: bool = True,
    similarity_mixfp_weight: float = 0.75,
) -> List["ConditionRecommendation"]:
    """Standalone precedent KNN search — no HTE data needed.

    Calls ``precedent.knn()`` directly (disk-cached featurized CSV rows).
    Used by the SIMILARITY fast path in ``api.py`` to avoid loading the
    large HTE pkl files.
    """
    try:
        from chemtools import precedent
        from chemtools.featurizers import reaction_pair as feat_pair
        from chemtools.recommend.utils import pick_electrophile_nucleophile
    except Exception:
        return []

    reaction_smiles = _build_reaction_smiles(
        reactant_a_smiles,
        reactant_b_smiles,
        product_smiles,
    )
    has_full_reaction = bool(product_smiles and reaction_smiles)
    use_drfp = has_full_reaction

    reactant_pool: List[str] = []
    if ">" in reaction_smiles:
        record = ReactionRecord.from_payload(normalize_reaction(reaction_smiles))
        reactant_pool = record.reactant_smiles
    if not reactant_pool:
        reactant_pool = [reactant_a_smiles]
        if reactant_b_smiles:
            reactant_pool.append(reactant_b_smiles)

    elec, nuc = pick_electrophile_nucleophile(reactant_pool)
    features = feat_pair.featurize_pair(elec, nuc).get("flat", {}) if (elec or nuc) else {}

    relax = {
        "use_drfp": use_drfp,
        "use_mixfp": bool(prefer_mixfp_for_similarity and has_full_reaction),
        "mixfp_weight": float(similarity_mixfp_weight),
        "reaction_smiles": reaction_smiles if has_full_reaction else "",
        "filter_by_reagent_database": False,
        "precedent_limit": max(top_k, 10),
    }
    knn_k = max(top_k, 10)
    if has_full_reaction:
        # Reaction-centric reranking needs a deeper candidate pool than the
        # final UI list size.
        knn_k = max(knn_k, 200)
        relax["precedent_limit"] = max(int(relax.get("precedent_limit", 10)), 200)

    family = reaction_type or None
    pack = precedent.knn(family=family, features=features, k=knn_k, relax=relax)
    precedents = list(pack.get("precedents", []) or [])
    normalized_source = _normalize_source_group(source_group) if source_group else ""
    precedent_source_filter = ""
    if normalized_source in {"literature", "datasets", "dataset"}:
        precedent_source_filter = "literature"
    elif normalized_source in {"rules", "motif"}:
        precedent_source_filter = normalized_source
    if precedent_source_filter:
        precedents = [
            prec
            for prec in precedents
            if _normalize_source_group(prec.get("source_group")) == precedent_source_filter
        ]

    reaction_match_cache: Dict[str, Set[str]] = {}

    def _reaction_text_key(rsmi: Any) -> str:
        text = str(rsmi or "").strip()
        if not text or ">" not in text:
            return ""
        parts = text.split(">")
        if len(parts) == 2 and ">>" in text:
            reactants, products = parts[0], parts[1]
            agents = ""
        elif len(parts) == 3:
            reactants, agents, products = parts
        else:
            return text

        def _canon_side(side: str) -> str:
            tokens = [tok.strip() for tok in side.split(".") if tok.strip()]
            tokens.sort()
            return ".".join(tokens)

        return ">".join(
            [
                _canon_side(reactants),
                _canon_side(agents),
                _canon_side(products),
            ]
        )

    def _reaction_match_keys(rsmi: Any) -> Set[str]:
        text = str(rsmi or "").strip()
        if not text or ">" not in text:
            return set()
        cached = reaction_match_cache.get(text)
        if cached is not None:
            return cached

        keys: Set[str] = {text}
        try:
            record = ReactionRecord.from_payload(normalize_reaction(text))
        except Exception:
            reaction_match_cache[text] = keys
            return keys

        normalized_text = str(record.normalized or "").strip()
        if normalized_text:
            keys.add(normalized_text)

        reactants = sorted(record.reactant_smiles)
        agents = sorted(
            component.preferred_smiles
            for component in record.agents
            if component.preferred_smiles
        )
        products = sorted(record.product_smiles)
        canonical = ">".join(
            [
                ".".join(reactants),
                ".".join(agents),
                ".".join(products),
            ]
        )
        if canonical:
            keys.add(canonical)

        reaction_match_cache[text] = keys
        return keys

    def _row_to_precedent(row: Dict[str, Any], *, similarity: float = 1.0) -> Dict[str, Any]:
        return {
            "reaction_id": row.get("reaction_id"),
            "dataset_reaction_id": row.get("dataset_reaction_id"),
            "reaction_smiles": row.get("reaction_smiles") or "",
            "similarity": similarity,
            "yield": row.get("yield_value"),
            "base_uid": row.get("base_uid"),
            "solvent_uid": row.get("solvent_uid"),
            "rxn_type": row.get("rxn_type"),
            "source_file": row.get("source_file"),
            "source_group": row.get("source_group"),
            "reference": row.get("reference"),
            "conditions": row.get("conditions"),
        }

    query_reaction_keys = _reaction_match_keys(reaction_smiles)
    query_text_key = _reaction_text_key(reaction_smiles)
    if query_text_key:
        query_reaction_keys.add(query_text_key)
    query_event_key = _reaction_event_key_from_reaction_smiles(reaction_smiles) if has_full_reaction else ""
    require_reaction_smiles = bool(has_full_reaction and query_event_key)
    has_exact_precedent = any(
        query_reaction_keys.intersection(_reaction_match_keys(prec.get("reaction_smiles")))
        for prec in precedents
    ) if query_reaction_keys else False

    if query_reaction_keys and not has_exact_precedent:
        try:
            from chemtools.precedent.loader import (
                _file_family_from_name,
                _infer_source_group_from_path,
                _iter_literature_files,
                _make_row_from_csv,
                _read_csv_records,
            )
        except Exception:
            _iter_literature_files = None  # type: ignore[assignment]

        exact_precedents: List[Dict[str, Any]] = []
        seen_exact_ids: Set[str] = set()
        if _iter_literature_files is not None:
            for path in _iter_literature_files():
                source_label = _normalize_source_group(_infer_source_group_from_path(path))
                if precedent_source_filter and source_label != precedent_source_filter:
                    continue
                try:
                    records = _read_csv_records(path)
                except Exception:
                    continue
                file_family = _file_family_from_name(path)
                for row_index, record in enumerate(records):
                    row_reaction = str(record.get("reaction_smiles") or "").strip()
                    if not row_reaction or ">" not in row_reaction:
                        continue
                    row_text_key = _reaction_text_key(row_reaction)
                    if (
                        row_reaction not in query_reaction_keys
                        and row_text_key not in query_reaction_keys
                    ):
                        continue
                    row = _make_row_from_csv(
                        record,
                        row_index=row_index,
                        file_family=file_family,
                        source_group=source_label,
                    )
                    if row is None:
                        continue
                    row_id = str(row.get("reaction_id") or "").strip()
                    if row_id and row_id in seen_exact_ids:
                        continue
                    if row_id:
                        seen_exact_ids.add(row_id)
                    exact_precedents.append(_row_to_precedent(row, similarity=1.0))

        if exact_precedents:
            existing_ids = {
                str(prec.get("reaction_id") or "").strip()
                for prec in precedents
                if prec.get("reaction_id")
            }

            def _yield_value(prec: Dict[str, Any]) -> float:
                value = prec.get("yield")
                if isinstance(value, (int, float)):
                    return float(value)
                try:
                    return float(str(value).strip())
                except Exception:
                    return -1.0

            prepend = [
                prec
                for prec in exact_precedents
                if str(prec.get("reaction_id") or "").strip() not in existing_ids
            ]
            prepend.sort(
                key=lambda prec: (_yield_value(prec), str(prec.get("reaction_id") or "")),
                reverse=True,
            )
            if prepend:
                precedents = prepend + precedents

    if require_reaction_smiles:
        precedents = [
            prec
            for prec in precedents
            if str(prec.get("reaction_smiles") or "").strip()
        ]

    reaction_event_similarity_cache: Dict[str, Optional[float]] = {}

    def _blend_similarity(row: Dict[str, Any], base_similarity: float) -> float:
        base = max(0.0, min(1.0, float(base_similarity)))
        if not query_event_key:
            return base
        row_reaction_smiles = str(row.get("reaction_smiles") or "").strip()
        if not row_reaction_smiles:
            return base
        cached = reaction_event_similarity_cache.get(row_reaction_smiles)
        if cached is None and row_reaction_smiles in reaction_event_similarity_cache:
            return base
        if row_reaction_smiles not in reaction_event_similarity_cache:
            row_event_key = _reaction_event_key_from_reaction_smiles(row_reaction_smiles)
            reaction_event_similarity_cache[row_reaction_smiles] = _reaction_event_similarity(
                query_event_key,
                row_event_key,
            )
        event_similarity = reaction_event_similarity_cache.get(row_reaction_smiles)
        if event_similarity is None:
            return base

        blend_weight = float(_PRECEDENT_REACTION_EVENT_BLEND_WEIGHT)
        if blend_weight < 0.0:
            blend_weight = 0.0
        elif blend_weight > 1.0:
            blend_weight = 1.0
        blended = ((1.0 - blend_weight) * base) + (blend_weight * float(event_similarity))
        # Keep exact reaction-smiles matches from being penalized by event fallback noise.
        if query_reaction_keys and query_reaction_keys.intersection(_reaction_match_keys(row_reaction_smiles)):
            blended = max(blended, base)
        return max(0.0, min(1.0, blended))

    def _condition_value(conditions: Dict[str, Any], key: str, fallback: Optional[str]) -> str:
        raw = conditions.get(key) if conditions else None
        if not raw:
            return fallback or ""
        return _format_list(raw)

    deduped: Dict[Tuple[str, str, str, str, str], ConditionRecommendation] = {}
    for prec in precedents:
        conditions = prec.get("conditions") or {}
        catalyst = _condition_value(conditions, "catalyst", None)
        ligand = _condition_value(conditions, "ligand", None)
        base = _condition_value(conditions, "base", prec.get("base_uid"))
        solvent = _condition_value(conditions, "solvent", prec.get("solvent_uid"))
        additive = _condition_value(conditions, "additive", None) or None
        coupling_reagent = _condition_value(conditions, "condensation_agent", None) or None

        similarity = _blend_similarity(prec, float(prec.get("similarity") or 0.0))
        yield_val = prec.get("yield")
        avg_yield = float(yield_val) if isinstance(yield_val, (int, float)) else 0.0
        success_rate = 100.0 if avg_yield >= 50.0 else 0.0

        rec = ConditionRecommendation(
            catalyst=catalyst,
            ligand=ligand,
            base=base,
            solvent=solvent,
            additive=additive,
            coupling_reagent=coupling_reagent,
            spectator_groups=prec.get("spectator_groups") or None,
            success_rate=success_rate,
            avg_yield=avg_yield,
            median_yield=avg_yield,
            num_experiments=1,
            avg_z_score=similarity,
            confidence_score=similarity * 100.0,
            match_score=similarity,
            reaction_type=prec.get("dataset_reaction_id") or prec.get("rxn_type"),
            reaction_id=prec.get("reaction_id"),
            reactant_types=("", ""),
            z_score_range=(similarity, similarity),
        )

        key = (rec.catalyst, rec.ligand, rec.base, rec.solvent, rec.additive or "")
        existing = deduped.get(key)
        if existing is None or rec.match_score > existing.match_score:
            deduped[key] = rec

    results = list(deduped.values())
    results.sort(key=lambda item: item.match_score, reverse=True)
    if top_k <= 0:
        return results
    return results[:top_k]


class HTERecommender:
    """
    HTE-based condition recommender using reactant type matching.

    Architecture:
    1. Load and index HTE database by reactant motifs (V2 Taxonomy)
    2. Classify input reactant SMILES to get motifs
    3. Match against database using motif combinations (with hierarchical fallback)
    4. Rank conditions by success rate and confidence
    5. Return top-k recommendations
    """
    
    def __init__(self, hte_db_path: str = "data/HTE_db"):
        """Initialize recommender with HTE database.

        Rule-derived screens should live under `data/HTE_db/rules` so they are
        loaded alongside literature and experimental sources.
        """
        self.db_path = Path(hte_db_path)
        self.df: Optional[pd.DataFrame] = None
        self.indexed_data: Any = {}
        self.reaction_type_patterns: Dict[Tuple[str, str], Counter] = {}
        self.transformation_indices: Any = {}

        df, indexed_data, patterns, trans_indices = _load_hte_database_cached(str(self.db_path))
        self.df = _ensure_rule_tier_column(df)
        self.indexed_data = _coerce_frame_lookup(self.df, indexed_data, cache_size=256)
        self.reaction_type_patterns = dict(patterns)
        self.transformation_indices = _coerce_frame_lookup(self.df, trans_indices, cache_size=512)
    
    def _detect_reactant_types(self, smiles: str) -> Tuple[List[str], Optional[str]]:
        """
        Detect reactant motifs and categories from SMILES using V2 Taxonomy.
        
        Returns:
            (motifs, category) e.g., (["Ar-Br", "Ar-OH"], "Aryl Halide")
        """
        if not smiles:
            return [], None
            
        analysis = featurize_molecule(smiles)
        
        # Extract motif IDs (e.g., "Ar-Br") and include alternate IDs when available.
        motifs: List[str] = []
        for hit in analysis.get("motifs", []):
            compound_id = hit.get("compound_id") or hit.get("id", "")
            if compound_id:
                motifs.append(compound_id)
            for alt_id in hit.get("alt_compound_ids", []) or []:
                if alt_id and alt_id != compound_id:
                    motifs.append(alt_id)
        motifs = _dedupe_list(motifs)

        context_scaffolds: List[str] = []
        nh_scaffolds = _load_nh_heterocycle_scaffolds()
        for entry in analysis.get("context_motifs", []) or []:
            cid = entry.get("compound_id") or entry.get("id")
            if cid and cid in nh_scaffolds:
                context_scaffolds.append(cid)

        if "AromN-H" in motifs and context_scaffolds:
            scaffold = context_scaffolds[0]
            idx = motifs.index("AromN-H")
            motifs = [m for m in motifs if m != "AromN-H"]
            if scaffold not in motifs:
                motifs.insert(min(idx, len(motifs)), scaffold)

        if motifs:
            alias_map = {
                "RCH2-NH2": ["Alkyl-NH2", "R2CH-NH2"],
                "RCH2-NHR": ["Alkyl-NHR", "R2CH-NHR"],
                "RCH2-NR2": ["Alkyl-NR2", "R2CH-NR2"],
            }
            compound_ids = _load_compound_ids()
            expanded = list(motifs)
            for motif in motifs:
                for alias in alias_map.get(motif, []):
                    if alias in compound_ids:
                        expanded.append(alias)
            motifs = _dedupe_list(expanded)
        
        # Use the first motif's category as a general category
        category = analysis.get("motifs", [{}])[0].get("category", "Unknown") if analysis.get("motifs") else "Unknown"
        
        return motifs, category

    def _detect_motif_set(self, smiles: str) -> Set[str]:
        """Detect motifs for a SMILES string (supports multi-part SMILES)."""
        if not smiles:
            return set()

        motifs: List[str] = []
        for part in smiles.split("."):
            part = part.strip()
            if not part:
                continue
            analysis = featurize_molecule(part)
            motifs.extend(
                [m.get("compound_id") or m.get("id", "") for m in analysis.get("motifs", []) if (m.get("compound_id") or m.get("id"))]
            )

        return set(_dedupe_list(motifs))
    
    def _filter_by_catalyst(self, df: pd.DataFrame, catalyst_filter: str) -> pd.DataFrame:
        """
        Filter dataframe by catalyst metal type.
        
        Args:
            df: DataFrame to filter
            catalyst_filter: Metal type or symbol (e.g., 'Pd', 'Cu', 'Ni', 'palladium', 'copper')
        
        Returns:
            Filtered DataFrame
        """
        # Normalize filter term
        filter_lower = catalyst_filter.lower()

        # Map common names to symbols
        metal_map = {
            "palladium": "Pd",
            "copper": "Cu",
            "nickel": "Ni",
            "iridium": "Ir",
            "rhodium": "Rh",
            "ruthenium": "Ru",
            "platinum": "Pt",
            "gold": "Au",
            "silver": "Ag",
            "iron": "Fe",
            "cobalt": "Co",
            "zinc": "Zn",
            "organocatalyst": "organocatalyst",
            "organic catalyst": "organocatalyst",
        }

        # Get the search term (symbol or name)
        search_term = metal_map.get(filter_lower, catalyst_filter)

        if "Catalyst_Type" in df.columns:
            mask = df["Catalyst_Type"].str.contains(search_term, case=False, na=False)
            if mask.any():
                return df[mask]

        # Fallback to catalyst names (case-insensitive)
        mask = df["Catalyst"].str.contains(search_term, case=False, na=False)
        return df[mask]
    
    def _predict_reaction_type(
        self,
        type_a: str,
        type_b: str,
    ) -> Tuple[Optional[str], float]:
        """
        Predict reaction type based on reactant type combination.

        Returns:
            (reaction_type, confidence_score)
        """
        key = _reactant_key([type_a, type_b])

        if key not in self.reaction_type_patterns:
            return None, 0.0

        rxn_counts = self.reaction_type_patterns[key]
        if not rxn_counts:
            return None, 0.0
        
        # Most common reaction type
        top_rxn = rxn_counts.most_common(1)[0]
        reaction_type = top_rxn[0]
        count = top_rxn[1]
        total = sum(rxn_counts.values())
        confidence = count / total
        
        return reaction_type, confidence
    
    def _score_transformation_match(
        self, 
        query_motifs: Set[str], 
        db_key: str,
        query_reacted: Optional[Set[str]] = None,
        query_formed: Optional[Set[str]] = None,
        query_spectators: Optional[Set[str]] = None
    ) -> float:
        """
        Score how well a database entry matches the query motifs.
        
        Logic:
        1. Must match the 'reacted' core motifs.
        2. Higher score for matching 'spectator' motifs.
        3. If product motifs are available, softly prefer matching 'formed' motifs.
        """
        reacted, formed, spectators = _parse_transformation_key(db_key)
        is_transform_key = "->" in str(db_key)

        # 1. Core match: query must contain all 'reacted' motifs
        core_set = query_reacted if query_reacted is not None else query_motifs
        if not reacted or not core_set:
            return 0.0
        if is_transform_key:
            core_match = _token_set_subset_compatible(reacted, core_set)
        else:
            # Signature-style key: allow query motifs to be a subset of a broader DB signature.
            core_match = _token_set_subset_compatible(reacted, core_set) or _token_set_subset_compatible(
                core_set, reacted
            )
        if not core_match:
            return 0.0
            
        # 2. Spectator match
        if query_spectators is None:
            query_spectators = query_motifs - reacted
        
        if not spectators and not query_spectators:
            spectator_score = 1.0
        elif not spectators:
            # Query has spectators but DB doesn't (clean substrate match)
            spectator_score = 0.3
        elif not query_spectators:
            # DB has spectators but query doesn't
            spectator_score = 0.5
        else:
            overlap = _token_set_overlap_compatible(spectators, query_spectators)
            union_count = len(spectators) + len(query_spectators) - overlap
            spectator_score = (overlap / union_count) if union_count else 1.0

        base_score = 0.5 + (0.5 * spectator_score)

        if query_formed is None:
            return base_score

        if not formed and not query_formed:
            return base_score
        if not formed:
            formed_score = 0.3
        elif not query_formed:
            formed_score = 0.6
        else:
            formed_overlap = _token_set_overlap_compatible(formed, query_formed)
            formed_union_count = len(formed) + len(query_formed) - formed_overlap
            if formed_union_count <= 0:
                formed_score = 0.5
            else:
                formed_score = formed_overlap / formed_union_count

        formed_multiplier = 1.0 + (0.1 * formed_score)
        return base_score * (0.7 + 0.3 * formed_score) * formed_multiplier
    
    def _calculate_confidence_score(
        self,
        avg_z_score: float,
        num_experiments: int,
        avg_yield: float
    ) -> float:
        """
        Calculate confidence score combining multiple factors.
        
        Formula: weighted combination of z-score (primary), sample size, and avg yield
        Z-score is the primary metric as it measures the success of a condition.
        """
        # Normalize sample size (log scale, cap at 100)
        sample_score = min(num_experiments, 100) / 100.0
        
        # Z-score weight (primary factor)
        # Normalize z-score: typical range is -3 to +3, scale to 0-1
        # Values above 3 are capped at 1.0, below -3 at 0.0
        z_score_normalized = max(0.0, min(1.0, (avg_z_score + 3.0) / 6.0))
        
        # Average yield weight (secondary factor)
        yield_weight = avg_yield / 100.0
        
        # Combined score with weights
        confidence = (
            0.6 * z_score_normalized +  # 60% from z-score (primary)
            0.25 * sample_score +        # 25% from sample size
            0.15 * yield_weight          # 15% from avg yield
        ) * 100.0
        
        return confidence
    
    def _aggregate_conditions(
        self,
        matched_df: pd.DataFrame,
        top_k: int = 10,
        min_experiments: int = 1,
        query_spectator_groups: Optional[Set[str]] = None,
    ) -> List[ConditionRecommendation]:
        """
        Aggregate and rank condition combinations from matched experiments.
        
        Strategy:
        1. Group by (catalyst, ligand, base, solvent) combination
        2. Score each condition by robust top-end z-score (best-seller style)
        3. Calculate success rate (yield > 50)
        4. Calculate avg/median yield
        5. Compute confidence score (weighted by robust z-score)
        6. Rank by match relevance, robust z-score, and support
        """
        recommendations = []
        if matched_df is None or matched_df.empty:
            return recommendations

        working_df = matched_df.copy()
        # Remove exact duplicate experimental rows without requiring ELN-like IDs.
        dedup_cols = [
            col
            for col in (
                "Reaction_Type_Standardized",
                "Reactant_A_Type",
                "Reactant_B_Type",
                "Catalyst",
                "Ligand",
                "Base",
                "Solvent",
                "Additive",
                "Secondary Solvent",
                "Coupling Reagent",
                "temperature_C",
                "atmosphere",
                "AREA_TOTAL_REDUCED",
                "z-Score",
                "Source_File",
                "Source_Row",
            )
            if col in working_df.columns
        ]
        if dedup_cols:
            working_df = working_df.drop_duplicates(subset=dedup_cols)
        if working_df.empty:
            return recommendations

        if "AREA_TOTAL_REDUCED" in working_df.columns:
            working_df["_yield_numeric"] = pd.to_numeric(
                working_df["AREA_TOTAL_REDUCED"],
                errors="coerce",
            ).fillna(0.0)
        else:
            working_df["_yield_numeric"] = 0.0
        if "z-Score" in working_df.columns:
            working_df["_z_numeric"] = pd.to_numeric(working_df["z-Score"], errors="coerce")
        else:
            working_df["_z_numeric"] = pd.Series([pd.NA] * len(working_df), index=working_df.index)
        if "match_score" in working_df.columns:
            working_df["_match_numeric"] = pd.to_numeric(working_df["match_score"], errors="coerce").fillna(0.0)
        else:
            working_df["_match_numeric"] = 1.0
        if "_condition_quality_multiplier" in working_df.columns:
            working_df["_quality_numeric"] = pd.to_numeric(
                working_df["_condition_quality_multiplier"],
                errors="coerce",
            ).fillna(1.0)
        else:
            working_df["_quality_numeric"] = 1.0
        if "spectator_score" in working_df.columns:
            working_df["_spectator_numeric"] = pd.to_numeric(
                working_df["spectator_score"],
                errors="coerce",
            ).fillna(0.0)
        else:
            working_df["_spectator_numeric"] = 0.0
        if "Source_Group" in working_df.columns:
            working_df["_source_group_norm"] = working_df["Source_Group"].apply(_normalize_source_group)
        else:
            working_df["_source_group_norm"] = ""
        
        # Define success threshold
        SUCCESS_THRESHOLD = 50.0
        
        # Group by the core condition set, then preserve optional fields when
        # they vary meaningfully so distinct regimes do not collapse together.
        condition_cols = ['Catalyst', 'Ligand', 'Base', 'Solvent']
        optional_condition_specs = [
            (("Secondary Solvent",), "_cond_secondary_solvent", "text"),
            (("Additive",), "_cond_additive", "text"),
            (("Coupling Reagent",), "_cond_coupling_reagent", "text"),
            (("temperature_C", "Temperature"), "_cond_temperature", "numeric"),
            (("atmosphere", "Atmosphere"), "_cond_atmosphere", "text"),
        ]
        resolved_optional_cols: Dict[str, str] = {}
        for candidates, temp_col, value_kind in optional_condition_specs:
            source_col = next((name for name in candidates if name in working_df.columns), "")
            if not source_col:
                continue
            resolved_optional_cols[temp_col] = source_col
            if value_kind == "numeric":
                numeric_series = pd.to_numeric(working_df[source_col], errors="coerce").round(3)
                normalized = numeric_series.astype(object).where(numeric_series.notna(), "")
            else:
                text_series = working_df[source_col].fillna("").astype(str).str.strip()
                normalized = text_series.where(~text_series.str.lower().eq("nan"), "")
            working_df[temp_col] = normalized
            distinct_values = {str(value).strip() for value in normalized.tolist()}
            if len(distinct_values) > 1:
                condition_cols.append(temp_col)
        
        grouped = working_df.groupby(condition_cols, dropna=False)
        
        for condition_tuple, group_df in grouped:
            # Rule-derived rows may encode curated single-condition guidance.
            current_min_exp = min_experiments
            if "_source_group_norm" in group_df.columns:
                source_norm = set(
                    value for value in group_df["_source_group_norm"] if str(value).strip()
                )
                is_rule = "rules" in source_norm
                if is_rule:
                    current_min_exp = 1
            
            if len(group_df) < current_min_exp:
                continue
            
            # Extract condition components from the grouped regime.
            if not isinstance(condition_tuple, tuple):
                condition_tuple = (condition_tuple,)
            group_values = dict(zip(condition_cols, condition_tuple))
            catalyst = group_values.get("Catalyst")
            ligand = group_values.get("Ligand")
            base = group_values.get("Base")
            solvent = group_values.get("Solvent")

            def _group_text(source_col: str) -> Optional[str]:
                if source_col in group_values:
                    value = str(group_values[source_col] or "").strip()
                    return value or None
                if source_col not in group_df.columns:
                    return None
                value = _first_nonempty_text(group_df[source_col])
                return value or None

            def _group_float(source_col: str) -> Optional[float]:
                if source_col in group_values:
                    value = group_values[source_col]
                    if value in ("", None) or pd.isna(value):
                        return None
                    try:
                        return float(value)
                    except Exception:
                        return None
                if source_col not in group_df.columns:
                    return None
                numeric = pd.to_numeric(group_df[source_col], errors="coerce").dropna()
                if numeric.empty:
                    return None
                try:
                    return float(numeric.iloc[0])
                except Exception:
                    return None

            sec_solvent = _group_text("_cond_secondary_solvent")
            if sec_solvent is None and "_cond_secondary_solvent" not in resolved_optional_cols:
                sec_solvent = _group_text("Secondary Solvent")
            additive = _group_text("_cond_additive")
            if additive is None and "_cond_additive" not in resolved_optional_cols:
                additive = _group_text("Additive")
            coupling_reagent = _group_text("_cond_coupling_reagent")
            if coupling_reagent is None and "_cond_coupling_reagent" not in resolved_optional_cols:
                coupling_reagent = _group_text("Coupling Reagent")
            temperature = _group_float("_cond_temperature")
            if temperature is None and "_cond_temperature" not in resolved_optional_cols:
                temperature = _group_float("temperature_C")
            atmosphere = _group_text("_cond_atmosphere")
            if atmosphere is None and "_cond_atmosphere" not in resolved_optional_cols:
                atmosphere = _group_text("atmosphere")
            
            # Calculate statistics
            yields = group_df["_yield_numeric"]
            num_exp = len(group_df)
            success_count = (yields > SUCCESS_THRESHOLD).sum()
            success_rate = (success_count / num_exp) * 100.0
            avg_yield = yields.mean()
            median_yield = yields.median()
            
            # Robust z-score: median of top-N z values to reduce outlier risk while
            # keeping a top-end performance signal.
            z_scores = group_df["_z_numeric"].dropna()
            if z_scores.empty:
                continue
            z_sorted = z_scores.sort_values(ascending=False).tolist()
            top_slice = z_sorted[: max(1, min(_BEST_SELLER_TOP_N, len(z_sorted)))]
            avg_z_score = float(pd.Series(top_slice, dtype=float).median())
            z_min = float(min(z_sorted))
            z_max = float(max(z_sorted))
            
            # Confidence score (uses z-score as primary factor)
            confidence = self._calculate_confidence_score(
                avg_z_score, num_exp, avg_yield
            )
            
            # Reaction type/category (most common, ignoring blanks)
            reaction_type = None
            if "Reaction_Type_Standardized" in group_df.columns:
                reaction_type = _first_nonempty_text(group_df["Reaction_Type_Standardized"])
            reaction_category = None
            if "Reaction_Category" in group_df.columns:
                reaction_category = _first_nonempty_text(group_df["Reaction_Category"])

            reaction_key = None
            if "_transformation_key" in group_df.columns:
                reaction_key = _first_nonempty_text(group_df["_transformation_key"])
            elif "Reaction_Key" in group_df.columns:
                reaction_key = _first_nonempty_text(group_df["Reaction_Key"])
            reaction_events = None
            if "Reaction_Events" in group_df.columns:
                reaction_events = _first_nonempty_text(group_df["Reaction_Events"])
            reaction_id = _format_source_reaction_ids(group_df)
            
            # Reactant types (from first row)
            reactant_types = (
                group_df.iloc[0]['Reactant_A_Type'],
                group_df.iloc[0]['Reactant_B_Type']
            )

            spectator_groups = ""
            if "spectator_groups" in group_df.columns:
                spectator_groups = _aggregate_spectator_groups(
                    group_df["spectator_groups"],
                    query_groups=query_spectator_groups,
                )
            spectator_score = 0.0
            if "_spectator_numeric" in group_df.columns:
                spectator_score = float(group_df["_spectator_numeric"].mean())
            
            match_score = float(group_df["_match_numeric"].mean()) if "_match_numeric" in group_df.columns else 1.0
            condition_quality_score = float(group_df["_quality_numeric"].mean()) if "_quality_numeric" in group_df.columns else 1.0
            missing_required_fields: Tuple[str, ...] = ()
            if "_missing_required_condition_fields" in group_df.columns:
                tokens = []
                for value in group_df["_missing_required_condition_fields"]:
                    text = str(value or "").strip()
                    if not text:
                        continue
                    tokens.extend(token.strip() for token in text.split("|") if token.strip())
                if tokens:
                    missing_required_fields = tuple(sorted(set(tokens)))
                    confidence = confidence * condition_quality_score
            
            rec = ConditionRecommendation(
                catalyst=catalyst if pd.notna(catalyst) else "",
                ligand=ligand if pd.notna(ligand) else "",
                base=base if pd.notna(base) else "",
                solvent=solvent if pd.notna(solvent) else "",
                secondary_solvent=sec_solvent if pd.notna(sec_solvent) else None,
                additive=additive if pd.notna(additive) else None,
                coupling_reagent=coupling_reagent if pd.notna(coupling_reagent) else None,
                temperature=temperature,
                atmosphere=atmosphere if pd.notna(atmosphere) else None,
                success_rate=success_rate,
                avg_yield=avg_yield,
                median_yield=median_yield,
                num_experiments=num_exp,
                avg_z_score=avg_z_score,
                confidence_score=confidence,
                match_score=match_score,
                condition_quality_score=condition_quality_score,
                missing_required_fields=missing_required_fields,
                reaction_type=reaction_type,
                reaction_category=reaction_category,
                reaction_id=reaction_id,
                reaction_key=reaction_key,
                reaction_events=reaction_events,
                reactant_types=reactant_types,
                z_score_range=(z_min, z_max),
                spectator_groups=spectator_groups,
                spectator_score=spectator_score,
            )
            
            recommendations.append(rec)
        
        # Sort by match relevance first, then robust z-score and support.
        recommendations.sort(
            key=lambda x: (
                x.match_score,
                x.condition_quality_score,
                x.avg_z_score,
                x.num_experiments,
                x.confidence_score,
                x.spectator_score,
            ),
            reverse=True,
        )

        if top_k <= 0:
            return recommendations
        return recommendations[:top_k]

    def _combine_recommendations_by_source(
        self,
        recommendations_by_source: Dict[str, List[ConditionRecommendation]],
        top_k: int,
    ) -> List[ConditionRecommendation]:
        """Interleave per-source recommendations into a combined list."""
        if not recommendations_by_source:
            return []

        by_source = {k: v for k, v in recommendations_by_source.items() if v}
        if not by_source:
            return []

        total_available = sum(len(items) for items in by_source.values())
        limit = total_available if top_k <= 0 else min(top_k, total_available)

        weight_map = {"motif": 3, "literature": 2, "datasets": 2, "rules": 1}
        priority_map = {"motif": 0, "literature": 1, "datasets": 1, "rules": 2, "other": 3, "unknown": 4}
        sources = sorted(by_source.keys(), key=lambda s: (priority_map.get(str(s), 5), str(s)))

        schedule: List[str] = []
        for source in sources:
            weight = weight_map.get(str(source), 1)
            schedule.extend([source] * max(weight, 1))
        if not schedule:
            schedule = sources

        indices = {source: 0 for source in sources}
        combined: List[ConditionRecommendation] = []
        while len(combined) < limit:
            progressed = False
            for source in schedule:
                idx = indices[source]
                items = by_source[source]
                if idx < len(items):
                    combined.append(items[idx])
                    indices[source] = idx + 1
                    progressed = True
                    if len(combined) >= limit:
                        break
            if not progressed:
                break
        return combined

    def _build_precedent_recommendations(
        self,
        reactant_a_smiles: str,
        reactant_b_smiles: Optional[str],
        product_smiles: Optional[str],
        reaction_type: Optional[str],
        top_k: int,
        source_group: Optional[str] = None,
        *,
        prefer_mixfp_for_similarity: bool = True,
        similarity_mixfp_weight: float = 0.75,
    ) -> List[ConditionRecommendation]:
        return _run_precedent_knn(
            reactant_a_smiles,
            reactant_b_smiles,
            product_smiles,
            reaction_type,
            top_k,
            source_group,
            prefer_mixfp_for_similarity=prefer_mixfp_for_similarity,
            similarity_mixfp_weight=similarity_mixfp_weight,
        )
    
    def recommend(
        self,
        reactant_a_smiles: str,
        reactant_b_smiles: Optional[str] = None,
        product_smiles: Optional[str] = None,
        top_k: int = 10,
        min_experiments: int = 2,
        reaction_type_filter: Optional[str] = None,
        catalyst_filter: Optional[str] = None,
        source_group: Optional[str] = None,
        reaction_key_only: bool = False,
        use_aryl_steric_electronic_weighting: bool = False,
        use_spectator_groups: bool = True,
        prefer_mixfp_for_similarity: bool = True,
        similarity_mixfp_weight: float = 0.75,
        force_precedent_search: bool = False,
    ) -> HTERecommendationResult:
        """
        Recommend conditions based on reactant SMILES.
        
            Args:
                reactant_a_smiles: SMILES of first reactant
                reactant_b_smiles: SMILES of second reactant (optional)
                product_smiles: Optional product SMILES to refine reacted/formed motifs in scoring
                top_k: Number of recommendations to return
                min_experiments: Minimum experiments for a condition to be recommended
                reaction_type_filter: Optional filter for specific reaction type
                catalyst_filter: Optional filter by metal type (e.g., 'Pd', 'Cu', 'Ni', 'palladium', 'copper')
                source_group: Optional source group filter (literature, rules, motif)
                reaction_key_only: Only match using reaction_key/signatures; no reactant-type fallback
                use_aryl_steric_electronic_weighting: Apply aryl steric/electronic weighting when available
                use_spectator_groups: Whether to apply spectator group weighting when available
                prefer_mixfp_for_similarity: Prefer MixFP routing/blending for precedent similarity when available
                similarity_mixfp_weight: MixFP share in similarity blend (0-1) when enabled
        
        Returns:
            HTERecommendationResult with ranked condition recommendations
        """
        result = HTERecommendationResult(
            reactant_a_smiles=reactant_a_smiles,
            reactant_b_smiles=reactant_b_smiles,
            product_smiles=product_smiles
        )
        recommend_started_at = time.perf_counter()
        stage_timings_ms: Dict[str, float] = {}
        _stage_order = (
            "query_prep_ms",
            "match_retrieval_ms",
            "filtering_and_enforcement_ms",
            "aggregation_ms",
            "precedent_merge_ms",
        )

        def _add_stage_timing(stage_name: str, started_at: float) -> None:
            elapsed_ms = round((time.perf_counter() - started_at) * 1000.0, 2)
            stage_timings_ms[stage_name] = round(
                float(stage_timings_ms.get(stage_name, 0.0)) + elapsed_ms, 2
            )

        def _finalize_timing_profile() -> None:
            profile = {name: float(stage_timings_ms.get(name, 0.0)) for name in _stage_order}
            profile["total_ms"] = round((time.perf_counter() - recommend_started_at) * 1000.0, 2)
            result.timing_ms = profile

        # Normalize source group and adjust min_experiments for rule-derived guidance.
        _t_query_prep = time.perf_counter()
        if source_group:
            source_group = _normalize_source_group(source_group)
            if source_group in {"rules"} and min_experiments > 1:
                min_experiments = 1
        normalized_source_group = _normalize_source_group(source_group) if source_group else ""
        fast_motif_mode = normalized_source_group == "motif"

        if reaction_type_filter:
            resolved_filter = _resolve_reaction_type_label(reaction_type_filter)
            if resolved_filter:
                reaction_type_filter = resolved_filter

        target_reaction_for_matching = ""
        if reaction_type_filter:
            target_reaction_for_matching = reaction_type_filter
        elif (
            result.predicted_reaction_type
            and result.predicted_reaction_type != "Unknown"
            and result.reaction_type_confidence >= 0.5
        ):
            target_reaction_for_matching = _resolve_reaction_type_label(
                result.predicted_reaction_type
            )

        def _filter_target_reaction(df: pd.DataFrame) -> pd.DataFrame:
            if (
                not target_reaction_for_matching
                or "Reaction_Type_Standardized" not in df.columns
            ):
                return df
            type_series = df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
            type_mask = type_series.apply(
                lambda value: _resolve_reaction_type_label(value) == target_reaction_for_matching if value else False
            )
            return df[type_mask]
        
        # Step 1: Detect reactant types (used for fallback only)
        type_a, cat_a = self._detect_reactant_types(reactant_a_smiles)
        collapsed_a = _collapse_motif_tokens(type_a) if type_a else ""
        result.reactant_a_type = collapsed_a
        result.reactant_a_category = cat_a

        if reactant_b_smiles:
            type_b, cat_b = self._detect_reactant_types(reactant_b_smiles)
            collapsed_b = _collapse_motif_tokens(type_b) if type_b else ""
            result.reactant_b_type = collapsed_b
            result.reactant_b_category = cat_b
        else:
            type_b, cat_b = [], ""
            result.reactant_b_type = ""
            result.reactant_b_category = ""
            collapsed_b = ""

        # Step 1.5: Detect reaction type using full reaction SMILES (if available)
        # This provides more accurate reaction type detection than pattern-based prediction
        reaction_smiles = None
        if reactant_b_smiles and product_smiles:
            reaction_smiles = f"{reactant_a_smiles}.{reactant_b_smiles}>>{product_smiles}"
        elif product_smiles:
            reaction_smiles = f"{reactant_a_smiles}>>{product_smiles}"
        elif reactant_b_smiles:
            reaction_smiles = f"{reactant_a_smiles}.{reactant_b_smiles}"
        else:
            reaction_smiles = reactant_a_smiles
        
        # Pre-eval: use product motifs to identify reacted vs spectator motifs
        query_reacted = None
        query_formed = None
        query_spectators = None
        query_spectator_groups: Set[str] = set()
        product_motifs_from_rxn = set()
        
        reaction_data: Dict[str, Any] = {}
        has_product_guided_detection = False
        if reaction_smiles and (">" in reaction_smiles or "." in reaction_smiles):
            try:
                rxn_features = _featurize_reaction_for_recommendation(
                    reaction_smiles,
                    skip_bond_analysis=True,
                )
                if isinstance(rxn_features, dict):
                    nested = rxn_features.get("reaction")
                    reaction_data = nested if isinstance(nested, dict) else rxn_features
                if product_smiles and reaction_data:
                    has_product_guided_detection = True
                rxn_type_data = reaction_data.get("reaction_type", {})
                detected_type = None
                detected_confidence = 0.0
                if isinstance(rxn_type_data, dict):
                    detected_type = rxn_type_data.get("reaction_type") or rxn_type_data.get("name")
                    detected_confidence = rxn_type_data.get("confidence", 0.0) or 0.0
                else:
                    detected_type = str(rxn_type_data).strip() if rxn_type_data else None
                    detected_confidence = reaction_data.get("confidence", 0.0) or 0.0
                
                # When product-guided featurization is available, keep that
                # signal as authoritative and avoid reactant-only fallback
                # overriding it downstream.
                if has_product_guided_detection and detected_type:
                    result.predicted_reaction_type = detected_type
                    result.reaction_type_confidence = float(detected_confidence)
                elif detected_type and detected_type != "Unknown" and detected_confidence > 0.5:
                    result.predicted_reaction_type = detected_type
                    result.reaction_type_confidence = detected_confidence
                
                # Extract product motifs from reaction featurization
                if product_smiles:
                    products = reaction_data.get("products", [])
                    for prod in products:
                        if isinstance(prod, dict):
                            for motif in prod.get("motifs", []):
                                if isinstance(motif, dict):
                                    compound_id = motif.get("compound_id") or motif.get("id")
                                    if compound_id:
                                        product_motifs_from_rxn.add(compound_id)
                                else:
                                    text = str(motif).strip()
                                    if text:
                                        product_motifs_from_rxn.add(text)
            except Exception:
                # Fallback to pattern-based prediction if featurization fails
                pass

        if product_smiles:
            # Use product motifs from reaction featurization if available, else detect separately
            if product_motifs_from_rxn:
                product_motifs = product_motifs_from_rxn
            else:
                product_motifs = self._detect_motif_set(product_smiles)
            result.product_type = ",".join(sorted(product_motifs))

            reactant_set = set(type_a) | set(type_b)
            reacted_set: Set[str] = set()
            formed_set: Set[str] = set()
            spectator_set: Set[str] = set()
            aggregates = reaction_data.get("aggregates") if isinstance(reaction_data, dict) else None
            if isinstance(aggregates, dict):
                reacted_set = set(aggregates.get("reacted_motifs") or [])
                formed_set = set(aggregates.get("formed_motifs") or [])
                spectator_set = set(aggregates.get("spectator_motifs") or [])
            if not (reacted_set or formed_set or spectator_set):
                reacted_set, formed_set, spectator_set = _derive_query_sets(reactant_set, product_motifs)
            query_reacted = reacted_set
            query_formed = formed_set
            query_spectators = spectator_set
            query_spectator_groups = _spectator_groups_from_motifs(spectator_set)

            result.reacted_motifs = tuple(sorted(reacted_set))
            result.formed_motifs = tuple(sorted(formed_set))
            result.spectator_motifs = tuple(sorted(spectator_set))
            result.query_spectator_groups = tuple(sorted(query_spectator_groups))

            result.query_reaction_key = _normalize_reaction_key(reaction_data.get("reaction_key"))

        query_events_payload = (
            reaction_data.get("reaction_events")
            if isinstance(reaction_data.get("reaction_events"), dict)
            else None
        )
        query_reaction_events_key = _reaction_events_to_match_key(query_events_payload or {})
        if not query_reaction_events_key:
            query_reaction_events_key = _reaction_events_to_match_key(
                build_reaction_events_payload(
                    reaction_data.get("reaction_key") or "",
                    query_events_payload,
                )
            )
        result.query_reaction_events_key = query_reaction_events_key or None

        # If no reaction key and no reactant types, return empty
        if not result.query_reaction_key and not result.query_reaction_events_key and not type_a:
            _add_stage_timing("query_prep_ms", _t_query_prep)
            _finalize_timing_profile()
            return result

        query_core_signature = ""
        query_ext_signature = ""
        if result.query_reaction_key:
            query_core_signature, query_ext_signature = _reaction_key_to_signatures(
                result.query_reaction_key
            )
        elif query_reacted is not None:
            query_core_signature = _encode_signature(query_reacted)
            query_ext_signature = _encode_signature(
                set(query_reacted) | (query_spectators or set())
            )
        if reaction_key_only and not (
            result.query_reaction_key
            or query_core_signature
            or query_ext_signature
            or result.query_reaction_events_key
        ):
            _add_stage_timing("query_prep_ms", _t_query_prep)
            _finalize_timing_profile()
            return result
        
        lookup_type_a = list(type_a)
        lookup_type_b = list(type_b)
        lookup_collapsed_a = collapsed_a
        lookup_collapsed_b = collapsed_b

        # Expand B(OH)2 ↔ B(OR)2 in lookup types so index lookups find
        # dataset entries that store only one boronate convention.
        def _expand_boronate(tokens: List[str]) -> List[str]:
            extra: List[str] = []
            for t in tokens:
                if t.endswith("-B(OH)2"):
                    alt = t[: -len("-B(OH)2")] + "-B(OR)2"
                    if alt not in tokens:
                        extra.append(alt)
                elif t.endswith("-B(OR)2"):
                    alt = t[: -len("-B(OR)2")] + "-B(OH)2"
                    if alt not in tokens:
                        extra.append(alt)
            return _dedupe_list(tokens + extra)

        lookup_type_a = _expand_boronate(lookup_type_a)
        lookup_type_b = _expand_boronate(lookup_type_b)
        if not lookup_type_a and not lookup_type_b and query_reacted:
            reacted_seed = _prioritize_motifs(sorted(query_reacted), query_reacted, query_spectators or set())
            if reacted_seed:
                if len(reacted_seed) > 1:
                    lookup_type_a = [reacted_seed[0]]
                    lookup_type_b = reacted_seed[1:]
                else:
                    lookup_type_a = [reacted_seed[0]]
                    lookup_type_b = []
                lookup_collapsed_a = _collapse_motif_tokens(lookup_type_a)
                lookup_collapsed_b = _collapse_motif_tokens(lookup_type_b)

        if (
            not target_reaction_for_matching
            and result.predicted_reaction_type
            and result.predicted_reaction_type != "Unknown"
            and result.reaction_type_confidence >= 0.5
        ):
            target_reaction_for_matching = _resolve_reaction_type_label(
                result.predicted_reaction_type
            )

        _add_stage_timing("query_prep_ms", _t_query_prep)

        # Step 3: Match against database (reaction-key first)
        _t_match = time.perf_counter()
        query_motifs = set(type_a) | set(type_b)
        if query_core_signature:
            query_motifs = set(_split_motif_tokens(query_core_signature))
        query_intramolecular_likely = _intramolecular_likely_from_fields(
            collapsed_a,
            collapsed_b,
            "",
        )
        
        key_match_df: Optional[pd.DataFrame] = None
        key_match_label = ""
        key_match_score = 0.0
        key_match_priority = 0
        skip_transformation_scan = normalized_source_group == "motif" and not reaction_key_only
        if result.query_reaction_key and result.query_reaction_key in self.transformation_indices:
            key_match_df = self.transformation_indices[result.query_reaction_key].copy()
            key_match_label = result.query_reaction_key
            key_match_score = 1.0
            key_match_priority = 0
        if key_match_df is None and query_core_signature:
            core_candidates = _signature_lookup_candidates(query_core_signature)
            for idx, candidate in enumerate(core_candidates):
                if candidate not in self.transformation_indices:
                    continue
                key_match_df = self.transformation_indices[candidate].copy()
                key_match_label = candidate
                key_match_score = 0.75 if idx == 0 else 0.7
                key_match_priority = 1 if idx == 0 else 2
                if idx > 0 or candidate != query_core_signature:
                    result.is_fallback_match = True
                break
        if key_match_df is None and query_ext_signature:
            ext_candidates = _signature_lookup_candidates(query_ext_signature)
            for idx, candidate in enumerate(ext_candidates):
                if candidate not in self.transformation_indices:
                    continue
                key_match_df = self.transformation_indices[candidate].copy()
                key_match_label = candidate
                key_match_score = 0.55 if idx == 0 else 0.5
                key_match_priority = 3 if idx == 0 else 4
                result.is_fallback_match = True
                break
        if (
            key_match_df is None
            and (
            result.query_reaction_events_key
            and result.query_reaction_events_key in self.transformation_indices
            )
        ):
            key_match_df = self.transformation_indices[result.query_reaction_events_key].copy()
            key_match_label = result.query_reaction_events_key
            key_match_score = 0.5
            key_match_priority = 5
            result.is_fallback_match = True

        scored_matches = []
        if key_match_df is not None and source_group:
            key_match_df = _filter_source_group(key_match_df, source_group)
            if key_match_df.empty:
                key_match_df = None
        if key_match_df is not None:
            key_match_df = _filter_target_reaction(key_match_df)
            if key_match_df.empty:
                key_match_df = None
        if key_match_df is None and not skip_transformation_scan:
            for db_key, group_df in self.transformation_indices.items():
                score = self._score_transformation_match(
                    query_motifs,
                    db_key,
                    query_reacted=query_reacted,
                    query_formed=query_formed,
                    query_spectators=query_spectators,
                )
                if score > 0:
                    # Weight the z-score by the match score
                    temp_df = group_df.copy()
                    if source_group:
                        temp_df = _filter_source_group(temp_df, source_group)
                        if temp_df.empty:
                            continue
                    temp_df = _filter_target_reaction(temp_df)
                    if temp_df.empty:
                        continue
                    
                    temp_df['match_score'] = score
                    temp_df['_transformation_key'] = db_key
                    _apply_intramolecular_boost(temp_df, query_intramolecular_likely)

                    temp_df['match_priority'] = 1
                        
                    scored_matches.append(temp_df)
        
        matched_df = None
        match_dfs: List[pd.DataFrame] = []
        if key_match_df is not None:
            key_match_df["match_score"] = key_match_score
            key_match_df["_transformation_key"] = key_match_label
            key_match_df["match_priority"] = key_match_priority
            _apply_intramolecular_boost(key_match_df, query_intramolecular_likely)
            match_dfs.append(key_match_df)
            result.matched_motifs = (key_match_label, "")
        if scored_matches:
            match_dfs.extend(scored_matches)
            # Use the highest scoring transformation as the "matched motifs" for display
            best_match = max(scored_matches, key=lambda x: x['match_score'].iloc[0])
            if "_transformation_key" in best_match.columns:
                best_key = best_match["_transformation_key"].iloc[0]
            else:
                best_key = best_match["Reaction_Type_Standardized"].iloc[0]
            result.matched_motifs = (best_key, "")

        key = _reactant_key([lookup_collapsed_a, lookup_collapsed_b])
        direct_match: Optional[pd.DataFrame] = None
        direct_key: Optional[str] = None
        fallback_used = False
        reacted_set = query_reacted or set()
        spectator_set = query_spectators or set()

        has_query_key = bool(
            result.query_reaction_key
            or query_core_signature
            or query_ext_signature
            or result.query_reaction_events_key
        )
        allow_direct_backfill = not reaction_key_only

        if not has_query_key or allow_direct_backfill:
            if key in self.indexed_data:
                direct_match = self.indexed_data[key].copy()
                direct_key = key
                if source_group:
                    direct_match = _filter_source_group(direct_match, source_group)
                    if direct_match.empty:
                        direct_match = None
                        direct_key = None
                if direct_match is not None:
                    direct_match = _filter_target_reaction(direct_match)
                    if direct_match.empty:
                        direct_match = None
                        direct_key = None
                if direct_match is not None:
                    direct_match['match_score'] = 1.0
                    direct_match['match_priority'] = 0
                    _apply_intramolecular_boost(direct_match, query_intramolecular_likely)
                    if not result.matched_motifs:
                        pick_a = _prioritize_motifs(lookup_type_a, reacted_set, spectator_set)
                        pick_b = _prioritize_motifs(lookup_type_b, reacted_set, spectator_set)
                        if pick_a or pick_b:
                            result.matched_motifs = (
                                pick_a[0] if pick_a else "",
                                pick_b[0] if pick_b else "",
                            )
                        else:
                            result.matched_motifs = (result.reactant_a_type, result.reactant_b_type)
            if direct_match is None:
                tiers_a = _build_fallback_tiers(lookup_type_a, reacted_set, spectator_set)
                tiers_b = _build_fallback_tiers(lookup_type_b, reacted_set, spectator_set)
                for tier_idx_a, list_a in enumerate(tiers_a):
                    for tier_idx_b, list_b in enumerate(tiers_b):
                        for ma in list_a:
                            for mb in list_b:
                                if not ma and not mb:
                                    continue
                                candidate = _reactant_key([ma, mb])
                                if candidate in self.indexed_data:
                                    direct_match = self.indexed_data[candidate].copy()
                                    direct_key = candidate
                                    if source_group:
                                        direct_match = _filter_source_group(direct_match, source_group)
                                        if direct_match.empty:
                                            direct_match = None
                                            direct_key = None
                                            continue
                                    direct_match = _filter_target_reaction(direct_match)
                                    if direct_match.empty:
                                        direct_match = None
                                        direct_key = None
                                        continue
                                    direct_match['match_score'] = 1.0
                                    direct_match['match_priority'] = 0
                                    _apply_intramolecular_boost(direct_match, query_intramolecular_likely)
                                    fallback_used = (
                                        fallback_used
                                        or candidate != key
                                        or tier_idx_a > 0
                                        or tier_idx_b > 0
                                    )
                                    if not result.matched_motifs:
                                        result.matched_motifs = (ma, mb)
                                    break
                            if direct_match is not None:
                                break
                        if direct_match is not None:
                            break
                    if direct_match is not None:
                        break

        fallback_match: Optional[pd.DataFrame] = None

        if direct_match is not None:
            match_dfs.append(direct_match)
        if fallback_match is not None:
            match_dfs.append(fallback_match)

        if match_dfs:
            matched_df = pd.concat(match_dfs, axis=0)
            if 'match_priority' in matched_df.columns:
                matched_df = matched_df.sort_values(['match_priority', 'match_score'], ascending=[True, False])
            elif 'match_score' in matched_df.columns:
                matched_df = matched_df.sort_values('match_score', ascending=False)
            matched_df = matched_df[~matched_df.index.duplicated(keep="first")]
            result.total_matching_experiments = len(matched_df)
            if not scored_matches and fallback_used:
                result.is_fallback_match = True

        _add_stage_timing("match_retrieval_ms", _t_match)

        if matched_df is None:
            _t_precedent = time.perf_counter()
            if (
                (source_group or "").lower() in {"", "literature", "datasets", "dataset"}
                and (force_precedent_search or not result.recommendations)
            ):
                precedent_recs = self._build_precedent_recommendations(
                    reactant_a_smiles,
                    reactant_b_smiles,
                    product_smiles,
                    reaction_type_filter or result.predicted_reaction_type,
                    top_k,
                    source_group=source_group,
                    prefer_mixfp_for_similarity=prefer_mixfp_for_similarity,
                    similarity_mixfp_weight=similarity_mixfp_weight,
                )
                if precedent_recs:
                    result.recommendations_by_source["precedent"] = precedent_recs
                    if not result.recommendations:
                        result.recommendations = (
                            precedent_recs[:top_k] if top_k > 0 else list(precedent_recs)
                        )
            _add_stage_timing("precedent_merge_ms", _t_precedent)
            _finalize_timing_profile()
            return result

        _t_filtering = time.perf_counter()
        if "Source_Group" in matched_df.columns:
            matched_df = matched_df.copy()
            matched_df["Source_Group"] = matched_df["Source_Group"].apply(_normalize_source_group)

        # Apply mild priority boost for rule rows using Rule_Tier (fallback=0, default=1, preferred=2+).
        # Keep this bounded so rule priors do not overpower reaction/experiment evidence.
        tier_col = "Rule_Tier" if "Rule_Tier" in matched_df.columns else (
            "rule_tier" if "rule_tier" in matched_df.columns else ""
        )
        if tier_col and "Source_Group" in matched_df.columns:
            rules_mask = matched_df["Source_Group"] == "rules"
            if rules_mask.any():
                rule_tiers = (
                    matched_df.loc[rules_mask, tier_col]
                    .apply(_coerce_rule_tier_value)
                    .astype(float)
                )
                rule_boost = (1.0 + (_RULE_TIER_BOOST_STEP * rule_tiers)).clip(
                    lower=1.0,
                    upper=_RULE_TIER_MAX_BOOST,
                )
                if "match_score" in matched_df.columns:
                    base_scores = pd.to_numeric(
                        matched_df.loc[rules_mask, "match_score"],
                        errors="coerce",
                    ).fillna(0.0)
                    matched_df.loc[rules_mask, "match_score"] = base_scores * rule_boost
                else:
                    matched_df.loc[rules_mask, "match_score"] = rule_boost
        
        # Apply reaction type filter if specified
        if reaction_type_filter:
            if "Reaction_Type_Standardized" in matched_df.columns:
                type_series = matched_df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
                target_reaction = _resolve_reaction_type_label(reaction_type_filter)
                type_mask = type_series.apply(
                    lambda value: _resolve_reaction_type_label(value) == target_reaction if value else False
                )
                matched_df = matched_df[type_mask]
            else:
                matched_df = matched_df.iloc[0:0]
        
        # Apply catalyst filter if specified
        if catalyst_filter:
            matched_df = self._filter_by_catalyst(matched_df, catalyst_filter)

        # Apply spectator group weighting if available
        if use_spectator_groups and query_spectator_groups and "spectator_groups" in matched_df.columns:
            row_group_sets = matched_df["spectator_groups"].apply(
                lambda value: set(_split_group_tokens(value))
            )
            spectator_scores = row_group_sets.apply(
                lambda row_groups: _spectator_similarity(query_spectator_groups, row_groups)
            )
            row_group_counts = row_group_sets.apply(len)
            rows_with_groups = int((row_group_counts > 0).sum())
            rows_total = int(len(row_group_sets))
            if len(spectator_scores) > 0:
                score_min = float(spectator_scores.min())
                score_max = float(spectator_scores.max())
                score_avg = float(spectator_scores.mean())
            else:
                score_min = 0.0
                score_max = 0.0
                score_avg = 0.0
            result.spectator_scoring_applied = True
            result.spectator_rows_with_groups = rows_with_groups
            result.spectator_rows_total = rows_total
            result.spectator_similarity_avg = score_avg
            result.spectator_similarity_range = (score_min, score_max)
            matched_df["spectator_score"] = spectator_scores
            weight = (1.0 - _SPECTATOR_MATCH_WEIGHT) + (_SPECTATOR_MATCH_WEIGHT * spectator_scores)
            if "match_score" in matched_df.columns:
                matched_df["match_score"] = matched_df["match_score"] * weight
            else:
                matched_df["match_score"] = weight

        if use_aryl_steric_electronic_weighting:
            rsmi_col = None
            if "reaction_smiles" in matched_df.columns:
                rsmi_col = "reaction_smiles"
            elif "Reaction_SMILES" in matched_df.columns:
                rsmi_col = "Reaction_SMILES"
            if rsmi_col:
                query_parts = _iter_smiles_parts(reactant_a_smiles)
                query_parts.extend(_iter_smiles_parts(reactant_b_smiles or ""))
                query_role_maps = [_aryl_role_scores_for_smiles(part) for part in query_parts]
                query_roles = _merge_role_scores(query_role_maps)
                if query_roles:
                    weights = [
                        _aryl_similarity_weight(query_roles, _aryl_role_scores_for_reaction(rsmi))
                        if rsmi
                        else 1.0
                        for rsmi in matched_df[rsmi_col].fillna("").astype(str)
                    ]
                    weight_series = pd.Series(weights, index=matched_df.index)
                    if "match_score" in matched_df.columns:
                        matched_df["match_score"] = matched_df["match_score"] * weight_series
                    else:
                        matched_df["match_score"] = weight_series

        # Step 2: Predict reaction type (using reactant patterns; fallback to match frequency)
        if (
            result.matched_motifs
            and not result.predicted_reaction_type
            and not has_product_guided_detection
        ):
            pred_rxn, rxn_conf = self._predict_reaction_type(
                result.reactant_a_type, result.reactant_b_type
            )
            result.predicted_reaction_type = pred_rxn
            result.reaction_type_confidence = rxn_conf
        if (
            (not result.predicted_reaction_type or result.predicted_reaction_type == "Unknown")
            and not has_product_guided_detection
        ):
            if "Reaction_Type_Standardized" in matched_df.columns:
                type_series = matched_df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
                type_series = type_series[type_series != ""]
                if not type_series.empty:
                    counts = type_series.value_counts()
                    result.predicted_reaction_type = counts.index[0]
                    result.reaction_type_confidence = float(counts.iloc[0] / counts.sum())

        # Enforce reaction-type + motif matching for rules/experiments.
        # Rules use a required-core-motif path to avoid over-constraining on
        # helper motifs (e.g., R_acidic-H) that are common in query keys.
        if "Source_Group" in matched_df.columns:
            enforce_rows_mask = matched_df["Source_Group"].isin({"rules", "motif"})
            rules_rows_mask = matched_df["Source_Group"] == "rules"
            motif_rows_mask = matched_df["Source_Group"] == "motif"
            if enforce_rows_mask.any():
                target_reaction = reaction_type_filter
                if (
                    not target_reaction
                    and result.predicted_reaction_type
                    and result.predicted_reaction_type != "Unknown"
                    and result.reaction_type_confidence >= 0.5
                ):
                    target_reaction = result.predicted_reaction_type
                target_reaction = _resolve_reaction_type_label(target_reaction)
                if target_reaction:
                    type_series = matched_df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
                    type_mask = type_series.apply(
                        lambda value: _resolve_reaction_type_label(value) == target_reaction if value else False
                    )
                else:
                    type_mask = pd.Series([True] * len(matched_df), index=matched_df.index)

                query_signatures: List[str] = []
                if query_core_signature:
                    query_signatures.append(query_core_signature)
                if query_ext_signature and query_ext_signature not in query_signatures:
                    query_signatures.append(query_ext_signature)
                reactant_signature = _reactant_types_to_signature([collapsed_a, collapsed_b])
                if reactant_signature and reactant_signature not in query_signatures:
                    query_signatures.append(reactant_signature)
                query_sets = [set(_split_motif_tokens(sig)) for sig in query_signatures if sig]

                if query_sets:
                    core_series = (
                        matched_df["Reactant_Signature_Core"]
                        if "Reactant_Signature_Core" in matched_df.columns
                        else pd.Series([""] * len(matched_df), index=matched_df.index)
                    )
                    ext_series = (
                        matched_df["Reactant_Signature_Ext"]
                        if "Reactant_Signature_Ext" in matched_df.columns
                        else pd.Series([""] * len(matched_df), index=matched_df.index)
                    )
                    rule_query_sets = _build_rule_required_core_query_sets(query_sets) or query_sets
                    rule_core_mask = core_series.apply(
                        lambda value: _signature_matches(rule_query_sets, value)
                    )
                    rule_ext_mask = ext_series.apply(
                        lambda value: _signature_matches(rule_query_sets, value)
                    )
                    exp_core_mask = core_series.apply(lambda value: _signature_matches(query_sets, value))
                    exp_ext_mask = ext_series.apply(lambda value: _signature_matches(query_sets, value))

                    motif_mask = pd.Series([False] * len(matched_df), index=matched_df.index)
                    if rules_rows_mask.any():
                        motif_mask.loc[rules_rows_mask] = (
                            rule_core_mask.loc[rules_rows_mask] | rule_ext_mask.loc[rules_rows_mask]
                        )
                    if motif_rows_mask.any():
                        motif_mask.loc[motif_rows_mask] = (
                            exp_core_mask.loc[motif_rows_mask] | exp_ext_mask.loc[motif_rows_mask]
                        )
                else:
                    motif_mask = pd.Series([False] * len(matched_df), index=matched_df.index)

                enforce_mask = enforce_rows_mask & type_mask & motif_mask
                matched_df = matched_df[~enforce_rows_mask | enforce_mask]

        (
            matched_df,
            required_catalyst_family,
            required_catalyst_classes,
            filtered_missing_catalyst_rows,
            retained_missing_catalyst_rows,
        ) = _enforce_required_catalyst_quality(
            matched_df,
            reaction_type_filter=reaction_type_filter,
            predicted_reaction_type=result.predicted_reaction_type,
            reaction_type_confidence=result.reaction_type_confidence,
        )
        result.required_catalyst_family = required_catalyst_family
        result.required_catalyst_classes = required_catalyst_classes
        result.catalyst_requirement_enforced = bool(required_catalyst_classes)
        result.filtered_missing_catalyst_rows = filtered_missing_catalyst_rows
        result.retained_missing_catalyst_rows = retained_missing_catalyst_rows
        (
            matched_df,
            condition_quality_family,
            missing_required_condition_fields,
            penalized_incomplete_condition_rows,
        ) = _apply_required_condition_quality_penalties(
            matched_df,
            reaction_type_filter=reaction_type_filter,
            predicted_reaction_type=result.predicted_reaction_type,
            reaction_type_confidence=result.reaction_type_confidence,
        )
        result.condition_quality_family = condition_quality_family
        result.missing_required_condition_fields = missing_required_condition_fields
        result.penalized_incomplete_condition_rows = penalized_incomplete_condition_rows
        (
            matched_df,
            plausibility_family,
            filtered_implausible_catalyst_rows,
            penalized_implausible_mechanism_rows,
            plausibility_issue_counts,
        ) = _apply_condition_plausibility_penalties(
            matched_df,
            reaction_type_filter=reaction_type_filter,
            predicted_reaction_type=result.predicted_reaction_type,
            reaction_type_confidence=result.reaction_type_confidence,
        )
        result.plausibility_family = plausibility_family
        result.filtered_implausible_catalyst_rows = filtered_implausible_catalyst_rows
        result.penalized_implausible_mechanism_rows = penalized_implausible_mechanism_rows
        result.plausibility_issue_counts = plausibility_issue_counts

        result.total_matching_experiments = len(matched_df)
        coverage_df = self.df
        if source_group and self.df is not None and "Source_Group" in self.df.columns:
            coverage_df = _filter_source_group(self.df, source_group)
        if coverage_df is not None and len(coverage_df) > 0:
            result.database_coverage = (len(matched_df) / len(coverage_df)) * 100.0
        else:
            result.database_coverage = 0.0

        _add_stage_timing("filtering_and_enforcement_ms", _t_filtering)

        # Step 4: Aggregate and rank conditions
        _t_aggregate = time.perf_counter()
        if len(matched_df) > 0:
            # Filter by detected reaction type if available and confidence is high
            if (
                result.predicted_reaction_type 
                and result.predicted_reaction_type != "Unknown"
                and result.reaction_type_confidence >= 0.5
                and not reaction_type_filter  # Don't override explicit filter
            ):
                # Apply detected reaction type as filter
                if "Reaction_Type_Standardized" in matched_df.columns:
                    type_series = matched_df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
                    detected_normalized = _resolve_reaction_type_label(result.predicted_reaction_type)
                    if detected_normalized:
                        type_mask = type_series.apply(
                            lambda x: _resolve_reaction_type_label(x) == detected_normalized if x else False
                        )
                        matched_df = matched_df[type_mask]
                        result.is_filtered_by_detected_type = True
            
            pool_size = max(top_k * 3, top_k)
            if "Source_Group" in matched_df.columns:
                by_source: Dict[str, List[ConditionRecommendation]] = {}
                for group_name, group_df in matched_df.groupby("Source_Group"):
                    label = str(group_name).strip()
                    if not label or label.lower() == "nan":
                        label = "unknown"
                    group_candidates = self._aggregate_conditions(
                        group_df,
                        pool_size,
                        min_experiments,
                        query_spectator_groups=query_spectator_groups,
                    )
                    if top_k > 0 and len(group_candidates) > top_k:
                        by_source[label] = self._select_diverse_conditions(
                            group_candidates, top_k, prioritize_performance=True
                        )
                    else:
                        by_source[label] = group_candidates
                result.recommendations_by_source = by_source
                result.recommendations = self._combine_recommendations_by_source(by_source, top_k)
            else:
                candidates = self._aggregate_conditions(
                    matched_df,
                    pool_size,
                    min_experiments,
                    query_spectator_groups=query_spectator_groups,
                )
                if top_k > 0 and len(candidates) > top_k:
                    result.recommendations = self._select_diverse_conditions(
                        candidates, top_k, prioritize_performance=True
                    )
                else:
                    result.recommendations = candidates
        _add_stage_timing("aggregation_ms", _t_aggregate)

        _t_precedent = time.perf_counter()
        if (
            (source_group or "").lower() in {"", "literature", "datasets", "dataset"}
            and (force_precedent_search or not result.recommendations)
        ):
            precedent_recs = self._build_precedent_recommendations(
                reactant_a_smiles,
                reactant_b_smiles,
                product_smiles,
                reaction_type_filter or result.predicted_reaction_type,
                top_k,
                source_group=source_group,
                prefer_mixfp_for_similarity=prefer_mixfp_for_similarity,
                similarity_mixfp_weight=similarity_mixfp_weight,
            )
            if precedent_recs:
                result.recommendations_by_source["precedent"] = precedent_recs
                if not result.recommendations:
                    result.recommendations = (
                        precedent_recs[:top_k] if top_k > 0 else list(precedent_recs)
                    )

        _add_stage_timing("precedent_merge_ms", _t_precedent)
        _finalize_timing_profile()
        return result

    def summarize_conditions(
        self,
        *,
        reaction_type_filter: Optional[str] = None,
        reactant_type_filters: Optional[Iterable[str]] = None,
        match_all_reactants: bool = False,
        source_group: Optional[str] = None,
        top_k: int = 10,
        min_experiments: int = 2,
    ) -> Dict[str, Any]:
        """Summarize top conditions for filtered dataset slices."""
        if self.df is None:
            return {
                "total_matching_experiments": 0,
                "database_coverage": 0.0,
                "recommendations": [],
            }

        if reaction_type_filter:
            resolved_filter = _resolve_reaction_type_label(reaction_type_filter)
            if resolved_filter:
                reaction_type_filter = resolved_filter

        filtered = self.df.copy()
        if source_group:
            filtered = _filter_source_group(filtered, source_group)

        if reaction_type_filter:
            filtered = filtered[
                filtered["Reaction_Type_Standardized"] == reaction_type_filter
            ]

        filtered = _filter_by_reactant_types(
            filtered, reactant_type_filters, match_all=match_all_reactants
        )

        coverage_df = self.df
        if source_group:
            coverage_df = _filter_source_group(self.df, source_group)

        total_matches = len(filtered)
        coverage = (
            (total_matches / len(coverage_df)) * 100.0
            if coverage_df is not None and len(coverage_df) > 0
            else 0.0
        )

        if total_matches == 0:
            return {
                "total_matching_experiments": 0,
                "database_coverage": coverage,
                "recommendations": [],
            }

        candidates = self._aggregate_conditions(
            filtered, max(top_k * 3, top_k), min_experiments
        )
        recommendations = (
            self._select_diverse_conditions(candidates, top_k, prioritize_performance=True)
            if top_k > 0 and len(candidates) > top_k
            else candidates
        )
        return {
            "total_matching_experiments": total_matches,
            "database_coverage": coverage,
            "recommendations": recommendations,
        }
    
    def generate_screening_set(
        self,
        reactant_a_smiles: str,
        reactant_b_smiles: Optional[str] = None,
        num_conditions: int = 24,
        min_experiments: int = 1,
        reaction_type_filter: Optional[str] = None,
        catalyst_filter: Optional[str] = None,
        diversity_strategy: str = "balanced"
    ) -> HTERecommendationResult:
        """
        Generate a diverse set of conditions for HTE screening (up to 24 for standard plate).
        
        This is the PRIMARY use case for HTE systems - generating a group of diverse conditions
        to test in parallel on a screening plate.
        
        Args:
            reactant_a_smiles: SMILES of first reactant
            reactant_b_smiles: SMILES of second reactant (optional)
            num_conditions: Number of conditions to generate (default 24 for 4x6 plate)
            min_experiments: Minimum experiments for a condition to be included
            reaction_type_filter: Optional filter for specific reaction type
            catalyst_filter: Optional filter by metal type (e.g., 'Pd', 'Cu', 'Ni')
            diversity_strategy: Strategy for condition selection:
                - "balanced": Mix of top performers + diverse alternatives
                - "top_performers": Focus on highest z-score conditions
                - "diverse": Maximize diversity across reagent space
        
        Returns:
            HTERecommendationResult with up to num_conditions diverse conditions
        """
        # Get initial recommendations with larger pool
        initial_top_k = max(num_conditions * 3, 50)  # Get 3x more for diversity selection
        
        result = self.recommend(
            reactant_a_smiles=reactant_a_smiles,
            reactant_b_smiles=reactant_b_smiles,
            top_k=initial_top_k,
            min_experiments=min_experiments,
            reaction_type_filter=reaction_type_filter,
            catalyst_filter=catalyst_filter
        )
        
        if len(result.recommendations) == 0:
            return result
        
        # Apply diversity strategy
        if diversity_strategy == "top_performers":
            # Simply take top N by z-score
            selected = result.recommendations[:num_conditions]
        
        elif diversity_strategy == "diverse":
            # Maximize diversity across all reagent dimensions
            selected = self._select_diverse_conditions(
                result.recommendations, 
                num_conditions,
                prioritize_performance=False
            )
        
        else:  # "balanced" (default)
            # Take top performers + diverse alternatives
            num_top = min(num_conditions // 3, len(result.recommendations))  # ~1/3 top performers
            num_diverse = num_conditions - num_top
            
            top_performers = result.recommendations[:num_top]
            remaining = result.recommendations[num_top:]
            
            if remaining:
                diverse_picks = self._select_diverse_conditions(
                    remaining, 
                    num_diverse,
                    prioritize_performance=True
                )
                selected = top_performers + diverse_picks
            else:
                selected = top_performers
        
        result.recommendations = selected[:num_conditions]
        return result
    
    def _select_diverse_conditions(
        self, 
        conditions: List[ConditionRecommendation],
        num_to_select: int,
        prioritize_performance: bool = True
    ) -> List[ConditionRecommendation]:
        """
        Select diverse conditions maximizing reagent variation.
        
        Args:
            conditions: Pool of conditions to select from
            num_to_select: Number of conditions to select
            prioritize_performance: If True, weight selection by z-score
        
        Returns:
            List of diverse conditions
        """
        if len(conditions) <= num_to_select:
            return conditions
        
        selected = []
        remaining = list(conditions)
        
        # Track used reagents for diversity
        used_catalysts = set()
        used_ligands = set()
        used_bases = set()
        used_solvents = set()
        
        while len(selected) < num_to_select and remaining:
            best_score = -1
            best_idx = 0
            
            for i, cond in enumerate(remaining):
                # Calculate diversity score (how many new reagents does this add?)
                diversity_score = 0
                if cond.catalyst not in used_catalysts:
                    diversity_score += 1
                if cond.ligand not in used_ligands:
                    diversity_score += 1
                if cond.base not in used_bases:
                    diversity_score += 1
                if cond.solvent not in used_solvents:
                    diversity_score += 1
                
                # Combine diversity with performance if requested
                if prioritize_performance:
                    # Normalize z-score to 0-1 range (assuming typical range -3 to +3)
                    normalized_zscore = (cond.avg_z_score + 3) / 6.0
                    normalized_zscore = max(0, min(1, normalized_zscore))
                    # Include relevance terms so spectator-aware matching remains
                    # effective even after diversity selection.
                    normalized_match = max(0.0, min(1.0, float(cond.match_score)))
                    normalized_spectator = max(0.0, min(1.0, float(cond.spectator_score)))

                    combined_score = (
                        (diversity_score / 4.0) * 0.45
                        + normalized_zscore * 0.30
                        + normalized_match * 0.20
                        + normalized_spectator * 0.05
                    )
                else:
                    combined_score = diversity_score
                
                if combined_score > best_score:
                    best_score = combined_score
                    best_idx = i
            
            # Select best condition
            selected_cond = remaining.pop(best_idx)
            selected.append(selected_cond)
            
            # Update used reagents
            used_catalysts.add(selected_cond.catalyst)
            used_ligands.add(selected_cond.ligand)
            used_bases.add(selected_cond.base)
            used_solvents.add(selected_cond.solvent)
        
        return selected
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics"""
        if self.df is None:
            return {}
        
        return {
            'total_experiments': len(self.df),
            'reaction_types': self.df['Reaction_Type_Standardized'].nunique(),
            'unique_type_combinations': len(self.indexed_data),
            'success_rate_overall': (self.df['AREA_TOTAL_REDUCED'] > 50).mean() * 100,
            'avg_yield': self.df['AREA_TOTAL_REDUCED'].mean(),
            'catalysts': self.df['Catalyst'].nunique(),
            'ligands': self.df['Ligand'].nunique(),
            'bases': self.df['Base'].nunique(),
            'solvents': self.df['Solvent'].nunique()
        }


def format_recommendation(rec: ConditionRecommendation, rank: int = 1) -> str:
    """Format a single recommendation for display"""
    lines = [
        f"\n{'='*80}",
        f"Recommendation #{rank}",
        f"{'='*80}",
        f"⭐ Avg Z-Score: {rec.avg_z_score:.2f} (Primary Ranking Metric)",
        f"Confidence Score: {rec.confidence_score:.1f}/100",
        f"Success Rate: {rec.success_rate:.1f}% ({rec.num_experiments} experiments)",
        f"Avg Yield: {rec.avg_yield:.1f}% | Median: {rec.median_yield:.1f}%",
        f"Match Score: {rec.match_score:.2f}",
        "",
        "🧪 CONDITIONS:",
        f"  Catalyst: {rec.catalyst}",
        f"  Ligand: {rec.ligand}",
        f"  Base: {rec.base}",
        f"  Solvent: {rec.solvent}"
    ]
    
    if rec.secondary_solvent:
        lines.append(f"  Secondary Solvent: {rec.secondary_solvent}")
    if rec.additive:
        lines.append(f"  Additive: {rec.additive}")
    if rec.coupling_reagent:
        lines.append(f"  Coupling Reagent: {rec.coupling_reagent}")
    if rec.temperature is not None:
        lines.append(f"  Temperature (C): {rec.temperature:g}")
    if rec.atmosphere:
        lines.append(f"  Atmosphere: {rec.atmosphere}")
    
    lines.extend([
        "",
        "📊 STATISTICS:",
        f"  Reaction Type: {rec.reaction_type}",
        f"  Z-Score: Avg={rec.avg_z_score:.2f}, Range=[{rec.z_score_range[0]:.2f}, {rec.z_score_range[1]:.2f}]",
        f"  Reactant Types: {rec.reactant_types[0]} + {rec.reactant_types[1]}"
    ])
    
    return "\n".join(lines)


def format_result(result: HTERecommendationResult) -> str:
    """Format complete recommendation result for display"""
    lines = [
        "\n" + "="*80,
        "HTE-BASED CONDITION RECOMMENDATION",
        "="*80,
        f"Reactant A: {result.reactant_a_smiles}",
        f"  Type: {result.reactant_a_type} ({result.reactant_a_category})"
    ]
    
    if result.reactant_b_smiles:
        lines.extend([
            f"Reactant B: {result.reactant_b_smiles}",
            f"  Type: {result.reactant_b_type} ({result.reactant_b_category})"
        ])

    if result.product_smiles:
        lines.append(f"Product: {result.product_smiles}")

    lines.extend([
        "",
        f"🎯 PREDICTED REACTION TYPE: {result.predicted_reaction_type}",
        f"   Confidence: {result.reaction_type_confidence*100:.1f}%",
    ])

    if result.matched_motifs and result.matched_motifs[0]:
        lines.append(f"   Matched Transformation: {result.matched_motifs[0]}")

    lines.extend([
        "",
        "📊 DATABASE MATCH:",
        f"   {result.total_matching_experiments} matching experiments",
        f"   ({result.database_coverage:.2f}% of database)",
        "",
        f"🏆 TOP RECOMMENDATIONS: {len(result.recommendations)} conditions found",
        "   (Ranked by Average Z-Score)",
        "="*80
    ])
    
    # Add individual recommendations
    for i, rec in enumerate(result.recommendations, 1):
        lines.append(format_recommendation(rec, i))
    
    return "\n".join(lines)
