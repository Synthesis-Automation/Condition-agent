from __future__ import annotations

import logging
import re
from functools import lru_cache
from typing import Any, Dict, List, Optional

from ..taxonomy.calculable_spec import load_calculable_feature_spec
from ..util import rdkit_helpers
from ..util.boolean_expr import evaluate as evaluate_rule
from ..util.smarts_cache import compile_smarts

logger = logging.getLogger(__name__)

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
    spec = dict(spec)
    spec["features"] = features
    return spec


@lru_cache(maxsize=1)
def _reactant_feature_index() -> Dict[str, Any]:
    spec = _load_feature_spec()
    features = spec.get("features", []) or []
    derived_shortcuts = spec.get("derived_shortcuts", []) or []

    members: List[Dict[str, Any]] = []
    categories: Dict[str, Dict[str, Any]] = {}

    for feature in features:
        if not isinstance(feature, dict):
            continue
        token = feature.get("token")
        if not token:
            continue
        meta = feature.get("metadata") or {}
        member_id = meta.get("reactant_member")
        category_id = meta.get("reactant_category")
        if not member_id or not category_id:
            continue
        members.append(
            {
                "token": token,
                "member_id": member_id,
                "member_name": meta.get("reactant_name") or meta.get("member_name") or member_id,
                "category_id": category_id,
                "category_name": meta.get("category_name") or category_id,
                "group": meta.get("group", ""),
                "description": meta.get("category_description") or meta.get("description") or "",
                "smarts": (_extract_smarts_any(feature) or [""])[0],
                "priority": meta.get("priority", 1),
            }
        )
        if category_id not in categories:
            categories[category_id] = {
                "category_id": category_id,
                "category_name": meta.get("category_name") or category_id,
                "token": None,
            }

    for shortcut in derived_shortcuts:
        if not isinstance(shortcut, dict):
            continue
        token = shortcut.get("token")
        meta = shortcut.get("metadata") or {}
        category_id = meta.get("reactant_category")
        if not token or not category_id:
            continue
        entry = categories.setdefault(
            category_id,
            {"category_id": category_id, "category_name": meta.get("category_name") or category_id},
        )
        entry["token"] = token

    return {"members": members, "categories": categories}


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
    categories_seen = set()

    index = _reactant_feature_index()
    for member in index.get("members", []):
        token = member.get("token")
        member_id = member.get("member_id")
        category_id = member.get("category_id")
        if not token or not member_id:
            continue
        if features.get(token, False):
            member_types.append(member_id)
            if category_id and category_id not in categories_seen:
                categories.append(category_id)
                categories_seen.add(category_id)

    for category_id, entry in index.get("categories", {}).items():
        token = entry.get("token")
        if not token:
            continue
        if features.get(token, False) and category_id not in categories_seen:
            categories.append(category_id)
            categories_seen.add(category_id)

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

    index = _reactant_feature_index()
    categories = index.get("categories", {})
    best: Optional[Dict[str, Any]] = None
    best_key: Optional[tuple] = None

    for member in index.get("members", []):
        token = member.get("token")
        member_id = member.get("member_id")
        category_id = member.get("category_id")
        if not token or not member_id or not category_id:
            continue
        if not reactant_features.get(token, False):
            continue
        category_entry = categories.get(category_id, {})
        category_token = category_entry.get("token")
        category_present = bool(reactant_features.get(category_token, False)) if category_token else False
        
        priority = member.get("priority", 1)
        
        match = {
            "category": category_id,
            "member_type": member_id,
            "name": member.get("member_name", member_id),
            "group": member.get("group", ""),
            "smarts": member.get("smarts", ""),
            "category_smarts": "",
            "description": member.get("description", ""),
            "priority": priority,
            "category_name": category_entry.get("category_name") or member.get("category_name", category_id),
            "category_present": category_present,
        }
        # Higher priority first, then longer ID as tie-breaker
        sort_key = (-priority, -len(member_id), member_id)
        if best is None or sort_key < best_key:
            best = match
            best_key = sort_key

    if best is None:
        for category_id, entry in categories.items():
            token = entry.get("token")
            if not token or not reactant_features.get(token, False):
                continue
            match = {
                "category": category_id,
                "member_type": category_id,
                "name": entry.get("category_name", category_id),
                "group": "",
                "smarts": "",
                "category_smarts": "",
                "description": "",
                "priority": 0,
                "category_name": entry.get("category_name", category_id),
                "category_present": True,
            }
            sort_key = (0, -len(category_id), category_id)
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
