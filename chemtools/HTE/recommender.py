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
from typing import List, Dict, Optional, Tuple, Any, Iterable, Set
from collections import defaultdict, Counter
from functools import lru_cache
import itertools
import re
import pandas as pd
from pathlib import Path
import json

from chemtools.featurizers.structural import featurize_molecule
from chemtools.smiles import normalize_reaction

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _infer_source_group(source_path: Optional[Path]) -> str:
    if not source_path:
        return "unknown"
    parts = [part.lower() for part in source_path.parts]
    for part in parts:
        if part in ("literature", "datasets", "dataset"):
            return "literature"
        if part == "rules":
            return "rules"
        if part in ("experiments", "experiment", "experiements"):
            return "experiments"
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
    candidates.extend(db_path.glob("*.jsonl"))

    for subdir in ("literature", "datasets", "rules", "experiments", "experiment", "experiements"):
        sub_path = db_path / subdir
        if not sub_path.exists():
            continue
        candidates.extend(sub_path.glob("*.csv"))
        candidates.extend(sub_path.glob("*.jsonl"))

    seen = set()
    ordered: List[Path] = []
    for path in sorted(candidates, key=lambda p: str(p)):
        key = str(path.resolve())
        if key in seen:
            continue
        seen.add(key)
        ordered.append(path)
    return ordered


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


def _normalize_hte_dataframe(df: pd.DataFrame, source_path: Optional[Path] = None) -> pd.DataFrame:
    df = df.copy()

    column_mapping = {
        "reaction_type": "Reaction_Type_Standardized",
        "reaction_category": "Reaction_Category",
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
    }

    df = df.rename(columns={k: v for k, v in column_mapping.items() if k in df.columns})
    if "Reaction_Type_Standardized" not in df.columns and "reaction_id" in df.columns:
        df = df.rename(columns={"reaction_id": "Reaction_Type_Standardized"})

    if "Reaction_Key" in df.columns:
        if "Reaction_Type_Standardized" not in df.columns or not any(df["Reaction_Type_Standardized"]):
            df["Reaction_Type_Standardized"] = df["Reaction_Key"]
        if "Reactant_Types_Key" not in df.columns or not any(df["Reactant_Types_Key"]):
            df["Reactant_Types_Key"] = df["Reaction_Key"]

    required_cols = [
        "Reaction_Type_Standardized", "Reactant_A_Type", "Reactant_B_Type",
        "Reactant_C_Type",
        "Catalyst", "Ligand", "Base", "Solvent", "Additive",
        "Secondary Solvent", "Coupling Reagent", "AREA_TOTAL_REDUCED", "z-Score",
        "Reactant_A_Category", "Reactant_B_Category", "Reaction_Category", "Is_Intramolecular",
        "Source_File", "Source_Group", "spectator_groups",
    ]
    for col in required_cols:
        if col not in df.columns:
            if col == "Is_Intramolecular":
                df[col] = df["Reactant_B_Type"].isna() | (df["Reactant_B_Type"] == "")
            elif col in ("Source_File", "Source_Group"):
                df[col] = ""
            else:
                df[col] = "" if col not in ["AREA_TOTAL_REDUCED", "z-Score"] else 0.0

    def _normalize_reactants_row(row: pd.Series) -> pd.Series:
        a_val, b_val, c_val, cleaned = _normalize_reactant_values(
            [row.get("Reactant_A_Type"), row.get("Reactant_B_Type"), row.get("Reactant_C_Type")]
        )
        row["Reactant_A_Type"] = a_val
        row["Reactant_B_Type"] = b_val
        row["Reactant_C_Type"] = c_val
        row["_reactant_count"] = len(cleaned)
        return row

    df = df.apply(_normalize_reactants_row, axis=1)
    df["Is_Intramolecular"] = df["_reactant_count"] <= 1
    df = df.drop(columns=["_reactant_count"])

    df["Reactant_Types_Key"] = df.apply(
        lambda row: _reactant_key(
            [row.get("Reactant_A_Type"), row.get("Reactant_B_Type"), row.get("Reactant_C_Type")]
        ),
        axis=1,
    )

    if source_path is not None:
        df["Source_File"] = _format_source_path(source_path)
        df["Source_Group"] = _infer_source_group(source_path)

    return df


@lru_cache(maxsize=4)
def _load_hte_database_cached(
    hte_db_path: str,
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame], Dict[str, Counter], Dict[str, pd.DataFrame]]:
    """Load and index the HTE database once per path (cached)."""
    db_path = Path(hte_db_path)
    if not db_path.exists():
        raise FileNotFoundError(f"HTE database not found: {db_path}")

    file_paths = _collect_hte_files(db_path)
    if not file_paths:
        raise FileNotFoundError(f"No HTE CSV/JSONL files found under: {db_path}")

    frames: List[pd.DataFrame] = []
    for path in file_paths:
        if path.suffix.lower() == ".jsonl":
            frame = _load_hte_jsonl(path)
        else:
            frame = _read_hte_csv(path)
        frame = _normalize_hte_dataframe(frame, source_path=path)
        frames.append(frame)

    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    print(f"Loaded HTE database: {len(df)} experiments from {len(file_paths)} files")

    indexed_data: Dict[str, pd.DataFrame] = {}
    reaction_type_patterns: Dict[str, Counter] = {}
    transformation_indices: Dict[str, pd.DataFrame] = {}

    print("Building reactant type indices...")

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
        indexed_data[key] = group_df

        rxn_types = group_df["Reaction_Type_Standardized"].value_counts()
        reaction_type_patterns[key] = Counter(rxn_types.to_dict())
    
    # Build transformation-aware index
    df["Reaction_Type_Standardized"] = df["Reaction_Type_Standardized"].fillna("Unknown")
    grouped_rxn = df.groupby("Reaction_Type_Standardized")
    for key, group_df in grouped_rxn:
        transformation_indices[key] = group_df

    print(f"Indexed {len(indexed_data)} reactant combinations and {len(transformation_indices)} transformation types")
    return df, indexed_data, reaction_type_patterns, transformation_indices


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


def _clean_reaction_label(value: Optional[str]) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text or text.lower() in {"nan", "unknown"}:
        return ""
    return text


def _format_reaction_id(reaction_type: Optional[str], reaction_category: Optional[str]) -> str:
    type_text = _clean_reaction_label(reaction_type)
    category_text = _clean_reaction_label(reaction_category)
    if type_text and category_text:
        if type_text.lower() == category_text.lower():
            return type_text
        return f"{type_text} / {category_text}"
    return type_text or category_text


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
    if label in ("experiments", "experiment", "experiements"):
        return "experiments"
    if label == "rules":
        return "rules"
    return label


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
            cid = motif.get("compound_id")
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
    normalized = normalize_reaction(reaction_smiles)
    reactants = normalized.get("reactants") or []
    role_maps: List[Dict[str, Dict[str, Optional[float]]]] = []
    for entry in reactants:
        if not isinstance(entry, dict):
            continue
        smi = entry.get("smiles_norm") or entry.get("largest_smiles") or entry.get("input")
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
    return text


def _normalize_reactant_values(values: Iterable[Any]) -> Tuple[str, str, str, List[str]]:
    cleaned = [_clean_reactant_value(v) for v in values]
    cleaned = [v for v in cleaned if v]
    a_val = cleaned[0] if len(cleaned) > 0 else ""
    b_val = cleaned[1] if len(cleaned) > 1 else ""
    c_val = cleaned[2] if len(cleaned) > 2 else ""
    return a_val, b_val, c_val, cleaned


_MOTIF_SPLIT_RE = re.compile(r"[|,]")
_COMPOUND_LOGIC_FILE = Path(__file__).resolve().parents[1] / "taxonomy" / "data" / "compound_logic.json"
_COMPOUND_SCOPE_FILE = Path(__file__).resolve().parents[1] / "taxonomy" / "data" / "organic_compounds.v1.3.json"
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
_HALIDE_SUFFIXES = {"Cl", "Br", "I", "F"}
_ARYL_HALIDE_MOTIFS = {"Ar-F", "Ar-Cl", "Ar-Br", "Ar-I", "Ar-X"}


@lru_cache(maxsize=1)
def _load_motif_sets() -> Dict[str, List[str]]:
    if not _COMPOUND_LOGIC_FILE.exists():
        return {}
    try:
        with _COMPOUND_LOGIC_FILE.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return {}
    raw_sets = payload.get("motif_sets") or {}
    motif_sets: Dict[str, List[str]] = {}
    for name, entry in raw_sets.items():
        members: List[str] = []
        if isinstance(entry, dict):
            members = entry.get("members") or []
        elif isinstance(entry, list):
            members = entry
        motif_sets[name] = [str(m).strip() for m in members if str(m).strip()]
    return motif_sets


@lru_cache(maxsize=1)
def _load_scope_map() -> Dict[str, List[str]]:
    if not _COMPOUND_SCOPE_FILE.exists():
        return {}
    try:
        with _COMPOUND_SCOPE_FILE.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return {}

    compounds = payload.get("compounds") or []
    any_entries: List[Dict[str, str]] = []
    by_b: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        cid = str(entry.get("id") or "").strip()
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
def _load_group_tags() -> Dict[str, Set[str]]:
    if not _COMPOUND_GROUPS_FILE.exists():
        return {}
    try:
        with _COMPOUND_GROUPS_FILE.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
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
def _load_compound_tags() -> Dict[str, Set[str]]:
    if not _COMPOUND_SCOPE_FILE.exists():
        return {}
    try:
        with _COMPOUND_SCOPE_FILE.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return {}
    group_tags = _load_group_tags()
    compounds = payload.get("compounds") or []
    tag_map: Dict[str, Set[str]] = {}
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        cid = str(entry.get("id") or "").strip()
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
    if not _COMPOUND_SCOPE_FILE.exists():
        return set()
    try:
        with _COMPOUND_SCOPE_FILE.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return set()
    compounds = payload.get("compounds") or []
    ids = set()
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        cid = str(entry.get("id") or "").strip()
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
    expanded = list(motifs)
    for motif in motifs:
        text = str(motif).strip()
        if "-" not in text:
            continue
        prefix, suffix = text.rsplit("-", 1)
        if suffix in _HALIDE_SUFFIXES:
            candidate = f"{prefix}-X"
            if candidate in compound_ids:
                expanded.append(candidate)
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
    return [token.strip() for token in _MOTIF_SPLIT_RE.split(text) if token.strip()]


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
        return text.split("-")[-1].strip()
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
    if not query_groups and not row_groups:
        return 1.0
    if not query_groups:
        return 0.7
    if not row_groups:
        return 0.3
    intersection = query_groups & row_groups
    union = query_groups | row_groups
    if not union:
        return 0.5
    return len(intersection) / len(union)


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
    if token in scope_map:
        return [token] + scope_map[token]
    return [token]


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
    Parse transformation key format: [Reacted] -> [Formed] || [Spectators]
    
    Returns:
        (reacted_set, formed_set, spectators_set)
    """
    if " -> " not in key or " || " not in key:
        # Fallback for old flat keys (treat all as reacted)
        tokens = _split_motif_tokens(key)
        expanded = _expand_motif_tokens(tokens, _load_motif_sets(), _load_scope_map())
        return set(expanded), set(), set()
    
    try:
        parts = key.split(" || ")
        
        def parse_part(p):
            p = p.strip()
            if p == "[]" or p == "None" or not p:
                return set()
            # Remove brackets if present
            if p.startswith("[") and p.endswith("]"):
                p = p[1:-1]
            tokens = _split_motif_tokens(p)
            expanded = _expand_motif_tokens(tokens, _load_motif_sets(), _load_scope_map())
            return set(expanded)

        spectators = parse_part(parts[1])
        
        rxn_parts = parts[0].split(" -> ")
        reacted = parse_part(rxn_parts[0])
        formed = parse_part(rxn_parts[1])
        
        return reacted, formed, spectators
    except:
        # Robust fallback
        tokens = _split_motif_tokens(key)
        expanded = _expand_motif_tokens(tokens, _load_motif_sets(), _load_scope_map())
        return set(expanded), set(), set()


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


def _load_hte_jsonl(path: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue

            reactant_types = _ensure_list(record.get("reactant_types"))
            a_type, b_type, c_type, cleaned = _normalize_reactant_values(reactant_types)
            reactant_categories = _ensure_list(record.get("reactant_categories"))
            catalyst_types = _ensure_list(record.get("catalyst_type"))
            conditions = record.get("conditions") or {}
            metrics = record.get("metrics") or {}

            row = {
                "Reaction_Type_Standardized": record.get("reaction_type") or "Unknown",
                "Reactant_A_Type": a_type,
                "Reactant_B_Type": b_type,
                "Reactant_C_Type": c_type,
                "Reactant_A_Category": reactant_categories[0] if len(reactant_categories) > 0 else "",
                "Reactant_B_Category": reactant_categories[1] if len(reactant_categories) > 1 else "",
                "Reactant_Types_Key": _reactant_key(cleaned),
                "Catalyst_Type": _format_list(catalyst_types),
                "Catalyst": _format_list(conditions.get("catalyst")),
                "Ligand": _format_list(conditions.get("ligand")),
                "Base": _format_list(conditions.get("base")),
                "Solvent": _format_list(conditions.get("solvent")),
                "Secondary Solvent": _format_list(conditions.get("secondary_solvent")),
                "Additive": _format_list(conditions.get("additive")),
                "Coupling Reagent": _format_list(conditions.get("coupling_reagent")),
                "AREA_TOTAL_REDUCED": metrics.get("area_total_reduced"),
                "z-Score": metrics.get("z_score"),
                "spectator_groups": record.get("spectator_groups", ""),
                "Is_Intramolecular": len(cleaned) <= 1,
            }
            rows.append(row)

    return pd.DataFrame(rows)


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
    
    # Statistics
    success_rate: float = 0.0  # % of experiments with yield > 50
    avg_yield: float = 0.0
    median_yield: float = 0.0
    num_experiments: int = 0
    avg_z_score: float = 0.0  # Average z-score (PRIMARY ranking metric for condition success)
    confidence_score: float = 0.0  # Secondary score considering z-score and sample size
    match_score: float = 1.0  # How well the transformation matched the query
    
    # Metadata
    reaction_type: Optional[str] = None
    reaction_category: Optional[str] = None
    reaction_id: Optional[str] = None
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
    
    # Recommendations
    recommendations: List[ConditionRecommendation] = field(default_factory=list)
    recommendations_by_source: Dict[str, List[ConditionRecommendation]] = field(default_factory=dict)
    
    # Metadata
    total_matching_experiments: int = 0
    database_coverage: float = 0.0  # % of database that matches this query
    is_fallback_match: bool = False
    matched_motifs: Optional[Tuple[str, str]] = None


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
        self.indexed_data: Dict[Tuple[str, str], pd.DataFrame] = {}
        self.reaction_type_patterns: Dict[Tuple[str, str], Counter] = {}
        self.transformation_indices: Dict[str, pd.DataFrame] = {}

        df, indexed_data, patterns, trans_indices = _load_hte_database_cached(str(self.db_path))
        self.df = df
        self.indexed_data = dict(indexed_data)
        self.reaction_type_patterns = dict(patterns)
        self.transformation_indices = dict(trans_indices)
    
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
            compound_id = hit.get("compound_id", "")
            if compound_id:
                motifs.append(compound_id)
            for alt_id in hit.get("alt_compound_ids", []) or []:
                if alt_id and alt_id != compound_id:
                    motifs.append(alt_id)
        motifs = _dedupe_list(motifs)

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
                [m.get("compound_id", "") for m in analysis.get("motifs", []) if m.get("compound_id")]
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
        
        # 1. Core match: query must contain all 'reacted' motifs
        core_set = query_reacted if query_reacted is not None else query_motifs
        if not reacted or not reacted.issubset(core_set):
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
            intersection = spectators & query_spectators
            union = spectators | query_spectators
            spectator_score = len(intersection) / len(union)

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
            formed_union = formed | query_formed
            if not formed_union:
                formed_score = 0.5
            else:
                formed_score = len(formed & query_formed) / len(formed_union)

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
        min_experiments: int = 1
    ) -> List[ConditionRecommendation]:
        """
        Aggregate and rank condition combinations from matched experiments.
        
        Strategy:
        1. Group by (catalyst, ligand, base, solvent) combination
        2. Calculate z-score statistics (primary ranking metric)
        3. Calculate success rate (yield > 50)
        4. Calculate avg/median yield
        5. Compute confidence score (weighted by z-score)
        6. Rank by average z-score (primary), then confidence score
        """
        recommendations = []
        
        # Define success threshold
        SUCCESS_THRESHOLD = 50.0
        
        # Group by condition combination
        condition_cols = ['Catalyst', 'Ligand', 'Base', 'Solvent']
        optional_cols = ['Secondary Solvent', 'Additive', 'Coupling Reagent']
        
        grouped = matched_df.groupby(condition_cols, dropna=False)
        
        for condition_tuple, group_df in grouped:
            if len(group_df) < min_experiments:
                continue
            
            # Extract condition components
            catalyst, ligand, base, solvent = condition_tuple
            
            # Get optional components (most common values)
            sec_solvent = group_df['Secondary Solvent'].mode().iloc[0] if not group_df['Secondary Solvent'].isna().all() else None
            additive = group_df['Additive'].mode().iloc[0] if not group_df['Additive'].isna().all() else None
            coupling_reagent = group_df['Coupling Reagent'].mode().iloc[0] if not group_df['Coupling Reagent'].isna().all() else None
            
            # Calculate statistics
            yields = group_df['AREA_TOTAL_REDUCED']
            num_exp = len(group_df)
            success_count = (yields > SUCCESS_THRESHOLD).sum()
            success_rate = (success_count / num_exp) * 100.0
            avg_yield = yields.mean()
            median_yield = yields.median()
            
            # Z-score statistics (primary ranking metric)
            z_scores = group_df['z-Score']
            
            # If match_score is present, weight the z-score
            if 'match_score' in group_df.columns:
                # Weighted average z-score
                avg_z_score = (z_scores * group_df['match_score']).sum() / group_df['match_score'].sum()
            else:
                avg_z_score = z_scores.mean()
                
            z_min = z_scores.min()
            z_max = z_scores.max()
            
            # Confidence score (uses z-score as primary factor)
            confidence = self._calculate_confidence_score(
                avg_z_score, num_exp, avg_yield
            )
            
            # Reaction type/category (most common)
            reaction_type = group_df['Reaction_Type_Standardized'].mode().iloc[0] if not group_df['Reaction_Type_Standardized'].isna().all() else None
            reaction_category = None
            if "Reaction_Category" in group_df.columns and not group_df["Reaction_Category"].isna().all():
                reaction_category = group_df["Reaction_Category"].mode().iloc[0]
            reaction_id = _format_reaction_id(reaction_type, reaction_category)
            
            # Reactant types (from first row)
            reactant_types = (
                group_df.iloc[0]['Reactant_A_Type'],
                group_df.iloc[0]['Reactant_B_Type']
            )
            
            match_score = group_df['match_score'].mean() if 'match_score' in group_df.columns else 1.0
            
            rec = ConditionRecommendation(
                catalyst=catalyst if pd.notna(catalyst) else "",
                ligand=ligand if pd.notna(ligand) else "",
                base=base if pd.notna(base) else "",
                solvent=solvent if pd.notna(solvent) else "",
                secondary_solvent=sec_solvent if pd.notna(sec_solvent) else None,
                additive=additive if pd.notna(additive) else None,
                coupling_reagent=coupling_reagent if pd.notna(coupling_reagent) else None,
                success_rate=success_rate,
                avg_yield=avg_yield,
                median_yield=median_yield,
                num_experiments=num_exp,
                avg_z_score=avg_z_score,
                confidence_score=confidence,
                match_score=match_score,
                reaction_type=reaction_type,
                reaction_category=reaction_category,
                reaction_id=reaction_id,
                reactant_types=reactant_types,
                z_score_range=(z_min, z_max)
            )
            
            recommendations.append(rec)
        
        # Sort by match score (prioritizing intramolecular matches if query is intramolecular),
        # then average z-score (primary performance metric), then confidence score.
        recommendations.sort(key=lambda x: (x.match_score, x.avg_z_score, x.confidence_score), reverse=True)
        
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

        weight_map = {"experiments": 3, "literature": 2, "datasets": 2, "rules": 1}
        priority_map = {"experiments": 0, "literature": 1, "datasets": 1, "rules": 2, "other": 3, "unknown": 4}
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
    ) -> List[ConditionRecommendation]:
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
        use_drfp = bool(product_smiles and reaction_smiles)

        reactant_pool: List[str] = []
        if ">" in reaction_smiles:
            normalized = normalize_reaction(reaction_smiles)
            for entry in normalized.get("reactants", []) or []:
                if not isinstance(entry, dict):
                    continue
                smi = entry.get("smiles_norm") or entry.get("largest_smiles") or entry.get("input")
                if smi:
                    reactant_pool.append(smi)
        if not reactant_pool:
            reactant_pool = [reactant_a_smiles]
            if reactant_b_smiles:
                reactant_pool.append(reactant_b_smiles)

        elec, nuc = pick_electrophile_nucleophile(reactant_pool)
        features = feat_pair.featurize_pair(elec, nuc).get("flat", {}) if (elec or nuc) else {}

        relax = {
            "use_drfp": use_drfp,
            "reaction_smiles": reaction_smiles if use_drfp else "",
            "filter_by_reagent_database": False,
        }

        family = reaction_type or None
        pack = precedent.knn(family=family, features=features, k=max(top_k, 10), relax=relax)
        precedents = list(pack.get("precedents", []) or [])

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

            similarity = float(prec.get("similarity") or 0.0)
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
        return results[:top_k]
    
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
        use_aryl_steric_electronic_weighting: bool = False,
        use_spectator_groups: bool = True,
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
            source_group: Optional source group filter (literature, rules, experiments)
            use_aryl_steric_electronic_weighting: Apply aryl steric/electronic weighting when available
            use_spectator_groups: Whether to apply spectator group weighting when available
        
        Returns:
            HTERecommendationResult with ranked condition recommendations
        """
        result = HTERecommendationResult(
            reactant_a_smiles=reactant_a_smiles,
            reactant_b_smiles=reactant_b_smiles,
            product_smiles=product_smiles
        )
        
        # Step 1: Detect reactant types
        type_a, cat_a = self._detect_reactant_types(reactant_a_smiles)
        result.reactant_a_type = ",".join(type_a) if type_a else ""
        result.reactant_a_category = cat_a
        
        if reactant_b_smiles:
            type_b, cat_b = self._detect_reactant_types(reactant_b_smiles)
            result.reactant_b_type = ",".join(type_b) if type_b else ""
            result.reactant_b_category = cat_b
        else:
            type_b, cat_b = [], ""
            result.reactant_b_type = ""
            result.reactant_b_category = ""
        
        # If no type detected, return empty
        if not type_a:
            return result

        # Pre-eval: use product motifs to identify reacted vs spectator motifs
        query_reacted = None
        query_formed = None
        query_spectators = None
        query_spectator_groups: Set[str] = set()
        if product_smiles:
            product_motifs = self._detect_motif_set(product_smiles)
            result.product_type = ",".join(sorted(product_motifs))

            reactant_set = set(type_a) | set(type_b)
            reacted_set, formed_set, spectator_set = _derive_query_sets(reactant_set, product_motifs)
            query_reacted = reacted_set
            query_formed = formed_set
            query_spectators = spectator_set
            query_spectator_groups = _spectator_groups_from_motifs(spectator_set)

            result.reacted_motifs = tuple(sorted(reacted_set))
            result.formed_motifs = tuple(sorted(formed_set))
            result.spectator_motifs = tuple(sorted(spectator_set))
        
        # Step 3: Match against database
        query_motifs = set(type_a) | set(type_b)
        is_query_intramolecular = (reactant_b_smiles is None)
        
        scored_matches = []
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
                
                # Boost if intramolecular status matches
                if 'Is_Intramolecular' in temp_df.columns:
                    mask = (temp_df['Is_Intramolecular'] == is_query_intramolecular)
                    temp_df.loc[mask, 'match_score'] = score * 1.2
                    temp_df.loc[~mask, 'match_score'] = score
                else:
                    temp_df['match_score'] = score

                temp_df['match_priority'] = 1
                    
                scored_matches.append(temp_df)
        
        matched_df = None
        match_dfs: List[pd.DataFrame] = []
        if scored_matches:
            match_dfs.extend(scored_matches)
            # Use the highest scoring transformation as the "matched motifs" for display
            best_key = max(scored_matches, key=lambda x: x['match_score'].iloc[0])['Reaction_Type_Standardized'].iloc[0]
            result.matched_motifs = (best_key, "")

        key = _reactant_key(list(type_a) + list(type_b))
        direct_match: Optional[pd.DataFrame] = None
        direct_key: Optional[str] = None
        fallback_used = False
        reacted_set = query_reacted or set()
        spectator_set = query_spectators or set()

        if key in self.indexed_data:
            direct_match = self.indexed_data[key].copy()
            direct_key = key
            if source_group:
                direct_match = _filter_source_group(direct_match, source_group)
                if direct_match.empty:
                    direct_match = None
                    direct_key = None
            if direct_match is not None:
                direct_match['match_score'] = 1.0
                direct_match['match_priority'] = 0
                if not result.matched_motifs:
                    pick_a = _prioritize_motifs(type_a, reacted_set, spectator_set)
                    pick_b = _prioritize_motifs(type_b, reacted_set, spectator_set)
                    if pick_a or pick_b:
                        result.matched_motifs = (
                            pick_a[0] if pick_a else "",
                            pick_b[0] if pick_b else "",
                        )
                    else:
                        result.matched_motifs = (result.reactant_a_type, result.reactant_b_type)
        if direct_match is None:
            list_a = _prioritize_motifs(type_a, reacted_set, spectator_set) or [""]
            list_b = _prioritize_motifs(type_b, reacted_set, spectator_set) or [""]
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
                        direct_match['match_score'] = 1.0
                        direct_match['match_priority'] = 0
                        if not result.matched_motifs:
                            result.matched_motifs = (ma, mb)
                        fallback_used = True
                        break
                if direct_match is not None:
                    break
        
        fallback_match: Optional[pd.DataFrame] = None
        expand_for_coverage = False
        if (
            direct_match is not None
            and self.df is not None
            and "Source_Group" in direct_match.columns
            and not source_group
        ):
            direct_groups = {g for g in direct_match["Source_Group"].unique() if str(g).strip()}
            available_groups = {g for g in self.df["Source_Group"].unique() if str(g).strip()}
            if available_groups and direct_groups and not (available_groups <= direct_groups):
                expand_for_coverage = True

        should_expand = direct_match is None or (
            direct_match is not None and len(direct_match) < min_experiments
        ) or expand_for_coverage
        if should_expand:
            tiers_a = _build_fallback_tiers(type_a, reacted_set, spectator_set)
            tiers_b = _build_fallback_tiers(type_b, reacted_set, spectator_set)
            tier_pairs: List[Tuple[int, int, int, List[str], List[str]]] = []
            for idx_a, list_a in enumerate(tiers_a):
                for idx_b, list_b in enumerate(tiers_b):
                    tier_pairs.append((idx_a + idx_b, idx_a, idx_b, list_a, list_b))
            tier_pairs.sort(key=lambda item: (item[0], item[1], item[2]))

            for _, idx_a, idx_b, list_a, list_b in tier_pairs:
                for ma in list_a:
                    for mb in list_b:
                        if not ma and not mb:
                            continue
                        candidate = _reactant_key([ma, mb])
                        if direct_key and candidate == direct_key:
                            continue
                        if candidate in self.indexed_data:
                            fallback_match = self.indexed_data[candidate].copy()
                            if source_group:
                                fallback_match = _filter_source_group(fallback_match, source_group)
                                if fallback_match.empty:
                                    fallback_match = None
                                    continue
                            if (
                                expand_for_coverage
                                and "Source_Group" in fallback_match.columns
                                and direct_match is not None
                            ):
                                direct_groups = {
                                    g
                                    for g in direct_match["Source_Group"].unique()
                                    if str(g).strip()
                                }
                                candidate_groups = {
                                    g
                                    for g in fallback_match["Source_Group"].unique()
                                    if str(g).strip()
                                }
                                if not (candidate_groups - direct_groups):
                                    fallback_match = None
                                    continue
                            if direct_match is None:
                                fallback_match['match_score'] = 1.0
                                fallback_match['match_priority'] = idx_a + idx_b
                            else:
                                fallback_match['match_score'] = 0.85
                                fallback_match['match_priority'] = 1 + idx_a + idx_b
                            if not result.matched_motifs:
                                result.matched_motifs = (ma, mb)
                            fallback_used = True
                            break
                    if fallback_match is not None:
                        break
                if fallback_match is not None:
                    break

        if direct_match is not None:
            match_dfs.append(direct_match)
        if fallback_match is not None:
            match_dfs.append(fallback_match)

        if match_dfs:
            matched_df = pd.concat(match_dfs, axis=0)
            if 'match_priority' in matched_df.columns:
                matched_df = matched_df.sort_values(['match_priority', 'match_score'], ascending=False)
            elif 'match_score' in matched_df.columns:
                matched_df = matched_df.sort_values('match_score', ascending=False)
            matched_df = matched_df[~matched_df.index.duplicated(keep="first")]
            result.total_matching_experiments = len(matched_df)
            if not scored_matches and fallback_used:
                result.is_fallback_match = True

        if matched_df is None:
            return result

        if "Source_Group" in matched_df.columns:
            matched_df = matched_df.copy()
            matched_df["Source_Group"] = matched_df["Source_Group"].apply(_normalize_source_group)
        
        # Apply reaction type filter if specified
        if reaction_type_filter:
            type_filtered = matched_df[matched_df['Reaction_Type_Standardized'] == reaction_type_filter]
            if type_filtered.empty and "Reaction_Category" in matched_df.columns:
                type_filtered = matched_df[matched_df["Reaction_Category"] == reaction_type_filter]
            matched_df = type_filtered
        
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
            if "match_score" in matched_df.columns:
                matched_df["match_score"] = matched_df["match_score"] * (0.7 + 0.3 * spectator_scores)
            else:
                matched_df["match_score"] = 0.7 + 0.3 * spectator_scores

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
        if result.matched_motifs and not result.predicted_reaction_type:
            pred_rxn, rxn_conf = self._predict_reaction_type(
                result.reactant_a_type, result.reactant_b_type
            )
            result.predicted_reaction_type = pred_rxn
            result.reaction_type_confidence = rxn_conf
        if not result.predicted_reaction_type or result.predicted_reaction_type == "Unknown":
            if "Reaction_Type_Standardized" in matched_df.columns:
                type_series = matched_df["Reaction_Type_Standardized"].fillna("").astype(str).str.strip()
                type_series = type_series[type_series != ""]
                if not type_series.empty:
                    counts = type_series.value_counts()
                    result.predicted_reaction_type = counts.index[0]
                    result.reaction_type_confidence = float(counts.iloc[0] / counts.sum())
        
        result.total_matching_experiments = len(matched_df)
        coverage_df = self.df
        if source_group and self.df is not None and "Source_Group" in self.df.columns:
            coverage_df = _filter_source_group(self.df, source_group)
        if coverage_df is not None and len(coverage_df) > 0:
            result.database_coverage = (len(matched_df) / len(coverage_df)) * 100.0
        else:
            result.database_coverage = 0.0
        
        # Step 4: Aggregate and rank conditions
        if len(matched_df) > 0:
            pool_size = max(top_k * 3, top_k)
            if "Source_Group" in matched_df.columns:
                by_source: Dict[str, List[ConditionRecommendation]] = {}
                for group_name, group_df in matched_df.groupby("Source_Group"):
                    label = str(group_name).strip()
                    if not label or label.lower() == "nan":
                        label = "unknown"
                    group_candidates = self._aggregate_conditions(
                        group_df, pool_size, min_experiments
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
                    matched_df, pool_size, min_experiments
                )
                if top_k > 0 and len(candidates) > top_k:
                    result.recommendations = self._select_diverse_conditions(
                        candidates, top_k, prioritize_performance=True
                    )
                else:
                    result.recommendations = candidates
        if (source_group or "").lower() in {"", "literature", "datasets", "dataset"}:
            precedent_recs = self._build_precedent_recommendations(
                reactant_a_smiles,
                reactant_b_smiles,
                product_smiles,
                reaction_type_filter or result.predicted_reaction_type,
                top_k,
            )
            if precedent_recs:
                result.recommendations_by_source["precedent"] = precedent_recs

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

        filtered = self.df.copy()
        if source_group:
            filtered = _filter_source_group(filtered, source_group)

        if reaction_type_filter:
            type_filtered = filtered[
                filtered["Reaction_Type_Standardized"] == reaction_type_filter
            ]
            if type_filtered.empty and "Reaction_Category" in filtered.columns:
                type_filtered = filtered[
                    filtered["Reaction_Category"] == reaction_type_filter
                ]
            filtered = type_filtered

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
                    
                    # 60% diversity, 40% performance
                    combined_score = (diversity_score / 4.0) * 0.6 + normalized_zscore * 0.4
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
        f"",
        f"🧪 CONDITIONS:",
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
    
    lines.extend([
        f"",
        f"📊 STATISTICS:",
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
        f"",
        f"🎯 PREDICTED REACTION TYPE: {result.predicted_reaction_type}",
        f"   Confidence: {result.reaction_type_confidence*100:.1f}%",
    ])

    if result.matched_motifs and result.matched_motifs[0]:
        lines.append(f"   Matched Transformation: {result.matched_motifs[0]}")

    lines.extend([
        f"",
        f"📊 DATABASE MATCH:",
        f"   {result.total_matching_experiments} matching experiments",
        f"   ({result.database_coverage:.2f}% of database)",
        f"",
        f"🏆 TOP RECOMMENDATIONS: {len(result.recommendations)} conditions found",
        f"   (Ranked by Average Z-Score)",
        "="*80
    ])
    
    # Add individual recommendations
    for i, rec in enumerate(result.recommendations, 1):
        lines.append(format_recommendation(rec, i))
    
    return "\n".join(lines)
