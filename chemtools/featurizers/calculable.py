from __future__ import annotations

import logging
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from ..archive.taxonomy import load_registry
from ..archive.taxonomy.calculable_spec import load_calculable_feature_spec
from ..util import rdkit_helpers
from ..util.boolean_expr import evaluate as evaluate_rule
from ..util.smarts_cache import compile_smarts

logger = logging.getLogger(__name__)

_GENERAL_REACTANT_CATEGORIES = {"Alkyl-C-H", "ArH"}
_HEURISTIC_NORMALIZE = str.maketrans(
    {
        "\u2265": ">=",
        "\u2264": "<=",
        "\u2013": "-",
        "\u03b2": "beta",
    }
)
_PROP_CMP_RE = re.compile(r"^(molecular_weight|logP|TPSA|HBA|HBD|rotatable_bonds)\\s*(<=|>=|<|>)\\s*(\\d+(?:\\.\\d+)?)$")
_FSP3_RE = re.compile(r"^fraction sp3 carbons\\s*(<=|>=|<|>)\\s*(\\d+(?:\\.\\d+)?)$")


def _normalize_heuristic(text: str) -> str:
    return text.translate(_HEURISTIC_NORMALIZE).strip()


def _load_reactant_types_from_registry() -> List[dict]:
    try:
        registry = load_registry()
    except Exception:
        return []

    entries: List[dict] = []
    for reactant in registry.reactant_types.values():
        members = []
        for member in reactant.members:
            meta = member.metadata or {}
            members.append(
                {
                    "id": member.id,
                    "name": member.name,
                    "feature_token": meta.get("feature_token"),
                    "metadata": dict(meta),
                }
            )
        entries.append(
            {
                "id": reactant.id,
                "name": reactant.name,
                "description": reactant.description,
                "category": reactant.category,
                "feature_token": (reactant.metadata or {}).get("feature_token"),
                "metadata": dict(reactant.metadata or {}),
                "members": members,
            }
        )
    return entries


def _load_reactant_types_from_file(path: Path) -> List[dict]:
    if not path.exists():
        return []
    try:
        import json
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return []

    if isinstance(payload, dict):
        entries = payload.get("entries") or []
    else:
        entries = payload
    return [entry for entry in entries if isinstance(entry, dict)]


@lru_cache(maxsize=1)
def _load_reactant_types() -> List[dict]:
    entries = _load_reactant_types_from_registry()
    if entries:
        return entries
    data_dir = Path(__file__).resolve().parent.parent / "archive" / "taxonomy" / "data"
    return _load_reactant_types_from_file(data_dir / "reactant_types.json")


def _extract_smarts_any(feature: dict) -> List[str]:
    detect = feature.get("detect") or {}
    smarts_any = detect.get("smarts_any") or []
    if isinstance(smarts_any, str):
        smarts_any = [smarts_any]
    return [s for s in smarts_any if isinstance(s, str)]


@lru_cache(maxsize=1)
def _load_feature_spec() -> Dict[str, Any]:
    spec = load_calculable_feature_spec()
    features = list(spec.get("features", []) or [])
    features_by_token = {
        item.get("token"): item
        for item in features
        if isinstance(item, dict) and item.get("token")
    }

    for entry in _load_reactant_types():
        category_id = entry.get("id")
        category_name = entry.get("name", category_id)
        category_token = entry.get("feature_token")
        category_smarts = ""
        if category_token:
            category_smarts_list = _extract_smarts_any(features_by_token.get(category_token, {}))
            if category_smarts_list:
                category_smarts = category_smarts_list[0]

        for member in entry.get("members", []) or []:
            member_id = member.get("id")
            if not member_id:
                continue
            member_token = member.get("feature_token")
            detect = {}
            if member_token:
                smarts_any = _extract_smarts_any(features_by_token.get(member_token, {}))
                if smarts_any:
                    detect["smarts_any"] = smarts_any
            metadata = {
                "reactant_category": category_id,
                "category_name": category_name,
                "member_name": member.get("name", member_id),
                "group": (entry.get("metadata") or {}).get("group", ""),
                "description": entry.get("description", ""),
                "category_smarts": category_smarts,
            }
            features.append(
                {
                    "token": f"{member_id}_reactant",
                    "type": "bool",
                    "detect": detect,
                    "metadata": metadata,
                }
            )

    return {**spec, "features": features}


def _get_property_values(mol: Any) -> Dict[str, float]:
    try:
        from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    except Exception:
        return {}

    return {
        "molecular_weight": float(Descriptors.MolWt(mol)),
        "logP": float(Descriptors.MolLogP(mol)),
        "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
        "HBA": float(Lipinski.NumHAcceptors(mol)),
        "HBD": float(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": float(Lipinski.NumRotatableBonds(mol)),
    }


def _compare_numeric(value: float, op: str, target: float) -> bool:
    if op == "<":
        return value < target
    if op == "<=":
        return value <= target
    if op == ">":
        return value > target
    if op == ">=":
        return value >= target
    return False


def _count_aromatic_rings(mol: Any) -> int:
    try:
        from rdkit.Chem import rdMolDescriptors
    except Exception:
        return 0
    return int(rdMolDescriptors.CalcNumAromaticRings(mol))


def _count_stereocenters(mol: Any) -> int:
    try:
        from rdkit import Chem
    except Exception:
        return 0
    return len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))


def _has_small_ring(mol: Any) -> bool:
    try:
        rings = mol.GetRingInfo().AtomRings()
    except Exception:
        return False
    return any(len(ring) in {3, 4} for ring in rings)


def _has_fused_ring(mol: Any) -> bool:
    try:
        ring_info = mol.GetRingInfo()
        if len(ring_info.AtomRings()) < 2:
            return False
        return any(
            atom.IsInRing() and ring_info.NumAtomRings(atom.GetIdx()) > 1
            for atom in mol.GetAtoms()
        )
    except Exception:
        return False


def _fraction_sp3(mol: Any) -> Optional[float]:
    try:
        from rdkit.Chem import rdMolDescriptors
    except Exception:
        return None
    return float(rdMolDescriptors.CalcFractionCSP3(mol))


def _heuristic_value(heuristic: str, mol: Any, feature_type: str) -> Optional[Any]:
    normalized = _normalize_heuristic(heuristic)

    if normalized == "count aromatic rings":
        return _count_aromatic_rings(mol)
    if normalized == "count stereogenic centers":
        return _count_stereocenters(mol)
    if normalized == "has stereogenic center":
        return _count_stereocenters(mol) > 0
    if normalized == "contains 3-membered or 4-membered ring":
        return _has_small_ring(mol)
    if normalized == "has fused ring system":
        return _has_fused_ring(mol)
    if normalized == "TPSA>=50 or (HBA+HBD)>=4":
        props = _get_property_values(mol)
        return (props.get("TPSA", 0.0) >= 50.0) or ((props.get("HBA", 0.0) + props.get("HBD", 0.0)) >= 4.0)
    if normalized == "TPSA<=20 and (HBA+HBD)<=2":
        props = _get_property_values(mol)
        return (props.get("TPSA", 0.0) <= 20.0) and ((props.get("HBA", 0.0) + props.get("HBD", 0.0)) <= 2.0)

    match = _PROP_CMP_RE.match(normalized)
    if match:
        props = _get_property_values(mol)
        key, op, value = match.groups()
        if key not in props:
            return None
        return _compare_numeric(props[key], op, float(value))

    match = _FSP3_RE.match(normalized)
    if match:
        value = _fraction_sp3(mol)
        if value is None:
            return None
        op, target = match.groups()
        return _compare_numeric(value, op, float(target))

    if normalized == "count [cX3]-[Cl,Br,I,F] plus [C;X2]=[C;X2]-[Cl,Br,I,F]":
        patterns = ["[cX3][Cl,Br,I,F]", "[C;X2]=[C;X2][Cl,Br,I,F]"]
        count = 0
        for smarts in patterns:
            pattern = compile_smarts(smarts, validate=False)
            if pattern:
                count += len(mol.GetSubstructMatches(pattern))
        return count if feature_type == "int" else count > 0

    return None


def _detect_from_smarts(mol: Any, detect: dict, feature_type: str) -> Optional[Any]:
    if "smarts_any" in detect:
        smarts_any = detect.get("smarts_any") or []
        if isinstance(smarts_any, str):
            smarts_any = [smarts_any]
        for smarts in smarts_any:
            pattern = compile_smarts(smarts, validate=False)
            if pattern and mol.HasSubstructMatch(pattern):
                return True
        return False

    if "smarts_count" in detect:
        smarts_count = detect.get("smarts_count")
        if isinstance(smarts_count, str):
            smarts_count = [smarts_count]
        count = 0
        for smarts in smarts_count or []:
            pattern = compile_smarts(smarts, validate=False)
            if pattern:
                count += len(mol.GetSubstructMatches(pattern))
        if feature_type == "int":
            return count
        return count > 0

    if "heuristic" in detect:
        return _heuristic_value(detect.get("heuristic", ""), mol, feature_type)

    return None


def detect_all_features(smiles: str) -> Dict[str, Any]:
    mol = rdkit_helpers.parse_smiles(smiles)
    if mol is None:
        return {}

    spec = _load_feature_spec()
    features: Dict[str, Any] = {}

    for feature in spec.get("features", []) or []:
        if not isinstance(feature, dict):
            continue
        token = feature.get("token")
        detect = feature.get("detect")
        feature_type = feature.get("type", "bool")
        if not token or not isinstance(detect, dict):
            continue
        value = _detect_from_smarts(mol, detect, feature_type)
        if value is None:
            value = 0 if feature_type == "int" else False
        features[token] = value

    for feature in spec.get("features", []) or []:
        if not isinstance(feature, dict):
            continue
        token = feature.get("token")
        expr = feature.get("derive") or feature.get("derived")
        if token and expr:
            features[token] = bool(evaluate_rule(expr, features))

    for shortcut in spec.get("derived_shortcuts", []) or []:
        if not isinstance(shortcut, dict):
            continue
        token = shortcut.get("token")
        expr = shortcut.get("derive") or shortcut.get("derived")
        if token and expr:
            features[token] = bool(evaluate_rule(expr, features))

    return features


def get_reactant_type_features(smiles: str) -> Dict[str, Any]:
    features = detect_all_features(smiles)
    if not features:
        return {}

    out = dict(features)
    member_types: List[str] = []
    categories: List[str] = []

    for entry in _load_reactant_types():
        category_id = entry.get("id")
        if not category_id:
            continue
        category_token = entry.get("feature_token")
        category_present = bool(features.get(category_token, False)) if category_token else False
        member_present = False
        for member in entry.get("members", []) or []:
            member_id = member.get("id")
            member_token = member.get("feature_token")
            if not member_id or not member_token:
                continue
            if features.get(member_token, False):
                member_present = True
                if member_id not in member_types:
                    member_types.append(member_id)
                out[f"{member_id}_reactant"] = True
        if category_present or member_present:
            categories.append(category_id)

    out["member_types"] = member_types
    out["categories"] = categories
    return out


def _member_feature_lookup() -> Dict[str, dict]:
    spec = _load_feature_spec()
    features = spec.get("features", []) or []
    return {item.get("token"): item for item in features if isinstance(item, dict)}


def classify_reactant_smiles(smiles: str) -> Optional[Dict[str, Any]]:
    reactant_features = get_reactant_type_features(smiles)
    if not reactant_features:
        return None

    features_by_token = _member_feature_lookup()
    best: Optional[Dict[str, Any]] = None
    best_key: Optional[tuple] = None

    for entry in _load_reactant_types():
        category_id = entry.get("id")
        if not category_id:
            continue
        category_name = entry.get("name", category_id)
        category_token = entry.get("feature_token")
        category_present = bool(reactant_features.get(category_token, False)) if category_token else False
        for member in entry.get("members", []) or []:
            member_id = member.get("id")
            member_token = member.get("feature_token")
            if not member_id or not member_token:
                continue
            if not reactant_features.get(member_token, False):
                continue
            member_feature = features_by_token.get(member_token, {})
            smarts_list = _extract_smarts_any(member_feature)
            smarts = smarts_list[0] if smarts_list else ""
            specificity = len(smarts) if smarts else 0
            is_general = category_id in _GENERAL_REACTANT_CATEGORIES
            match = {
                "category": category_id,
                "member_type": member_id,
                "name": member.get("name", member_id),
                "group": (entry.get("metadata") or {}).get("group", ""),
                "smarts": smarts,
                "category_smarts": "",
                "description": entry.get("description", ""),
                "specificity": specificity,
                "is_general": is_general,
                "category_name": category_name,
                "category_present": category_present,
            }
            sort_key = (is_general, -specificity, member_id)
            if best is None or sort_key < best_key:
                best = match
                best_key = sort_key

        if best is None and category_present:
            is_general = category_id in _GENERAL_REACTANT_CATEGORIES
            match = {
                "category": category_id,
                "member_type": category_id,
                "name": category_name,
                "group": (entry.get("metadata") or {}).get("group", ""),
                "smarts": "",
                "category_smarts": "",
                "description": entry.get("description", ""),
                "specificity": 0,
                "is_general": is_general,
                "category_name": category_name,
                "category_present": category_present,
            }
            sort_key = (is_general, 0, category_id)
            if best is None or sort_key < best_key:
                best = match
                best_key = sort_key

    return best


__all__ = [
    "detect_all_features",
    "get_reactant_type_features",
    "classify_reactant_smiles",
    "_load_feature_spec",
]
