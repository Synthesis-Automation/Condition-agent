#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ar- electronics analysis (POC v1): Gasteiger charge delta at Ar ipso -> 0–10 score.

Implements `CODEX_PLAN_Ar_Electronics_POC_v1.md`:
- detect Ar-* anchor sites from templated atomic features (atom-map `:1`)
- choose smallest aromatic ring containing ipso (optional, used for ring mean)
- compute Gasteiger charges and delta_q = q_ipso - q_ring_mean
- map to electron-poor score (0–10): round(5 + k*delta_q), clamped
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple


def _ensure_repo_root_on_syspath() -> None:
    root = Path(__file__).resolve().parents[3]
    root_str = str(root)
    if root_str not in sys.path:
        sys.path.insert(0, root_str)


try:
    from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available
    from chemtools.util.smarts_cache import compile_smarts
except Exception:
    _ensure_repo_root_on_syspath()
    from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available
    from chemtools.util.smarts_cache import compile_smarts


DEFAULT_DATA_DIR = Path(__file__).resolve().parent


@dataclass(frozen=True)
class AnchorSite:
    label: str
    feature_token: str
    ipso_atom_idx: int


def resolve_data_path(path: Path, *, base_dir: Path = DEFAULT_DATA_DIR) -> Path:
    candidate = path
    if not candidate.is_absolute() and not candidate.exists():
        candidate = base_dir / candidate
    return candidate


def _read_json_object(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object in {path}, got {type(payload).__name__}")
    return payload


def load_compiled_features(path: Path) -> Dict[str, Any]:
    return _read_json_object(path)


def load_groups(path: Path) -> Dict[str, Any]:
    return _read_json_object(path)


def _compose_label(left: str, right: str) -> str:
    if left.endswith("-") and right.startswith("-"):
        return left + right[1:]
    return left + right


def _find_anchor_qidx_by_mapnum(query_mol: Any, map_num: int) -> int:
    for i, atom in enumerate(query_mol.GetAtoms()):
        if atom.GetAtomMapNum() == map_num:
            return i
    raise ValueError(f"No atom with AtomMapNum=={map_num} in query SMARTS")


def find_ar_anchor_sites(
    mol: Any,
    *,
    compiled_features: Dict[str, Any],
    groups_doc: Dict[str, Any],
    context_group_id: str = "Ar",
    context_atom_map_num: int = 1,
) -> List[AnchorSite]:
    g_by_id = {g["id"]: g for g in groups_doc.get("groups", [])}
    left_label = g_by_id.get(context_group_id, {}).get("name", f"{context_group_id}-")

    atomic = compiled_features.get("atomic", [])
    if not isinstance(atomic, list):
        raise ValueError("compiled.atomic must be a list")

    sites: List[AnchorSite] = []
    for feat in atomic:
        if not isinstance(feat, dict):
            continue
        src = feat.get("source", {})
        if not isinstance(src, dict) or src.get("kind") != "templated":
            continue

        gids = src.get("groups", []) or []
        if not isinstance(gids, list) or not gids or gids[0] != context_group_id:
            continue

        token = str(feat.get("token", "")).strip()
        if not token:
            continue

        detect = feat.get("detect", {}) or {}
        smarts_any = detect.get("smarts_any", []) or []
        if not isinstance(smarts_any, list) or not smarts_any:
            continue

        smarts = smarts_any[0]
        if not isinstance(smarts, str) or not smarts.strip():
            continue

        qmol = compile_smarts(smarts, validate=True)
        if qmol is None:
            raise ValueError(f"SMARTS compiled to None for {token}: {smarts}")

        anchor_qidx = _find_anchor_qidx_by_mapnum(qmol, context_atom_map_num)

        right_id = str(gids[1]) if len(gids) > 1 else ""
        right_label = g_by_id.get(right_id, {}).get("name", right_id)
        label = _compose_label(left_label, right_label) if right_label else left_label

        for match in mol.GetSubstructMatches(qmol):
            ipso_idx = int(match[int(anchor_qidx)])
            sites.append(AnchorSite(label=label, feature_token=token, ipso_atom_idx=ipso_idx))

    seen: set[Tuple[int, str]] = set()
    dedup: List[AnchorSite] = []
    for s in sites:
        key = (int(s.ipso_atom_idx), str(s.feature_token))
        if key in seen:
            continue
        seen.add(key)
        dedup.append(s)
    return dedup


def _is_aromatic_ring(mol: Any, ring: Sequence[int]) -> bool:
    return all(mol.GetAtomWithIdx(int(i)).GetIsAromatic() for i in ring)


def choose_smallest_aromatic_ring_containing(mol: Any, atom_idx: int) -> Optional[Tuple[int, ...]]:
    rings = mol.GetRingInfo().AtomRings()
    candidates = [r for r in rings if atom_idx in r and _is_aromatic_ring(mol, r)]
    if not candidates:
        return None
    candidates.sort(key=len)
    return tuple(int(i) for i in candidates[0])


def _round_half_up(x: float) -> int:
    return int(x + 0.5)


def _clamp_int(x: int, lo: int, hi: int) -> int:
    return int(max(lo, min(hi, int(x))))


def electron_poor_score_0_10(delta_q: Optional[float], *, k: float) -> int:
    if delta_q is None:
        return 5
    raw = 5.0 + float(k) * float(delta_q)
    return _clamp_int(_round_half_up(raw), 0, 10)


def _safe_float_prop(atom: Any, prop: str) -> Optional[float]:
    try:
        if not atom.HasProp(prop):
            return None
        v = float(atom.GetProp(prop))
        if v != v:  # NaN
            return None
        return float(v)
    except Exception:
        return None


def analyze_smiles_reactivity_electronics_v1(
    smiles: str,
    *,
    compiled_features_path: Path = Path("calculable_features.compiled.v1.json"),
    groups_path: Path = Path("organic_groups.v1.json"),
    k: float = 25.0,
) -> Dict[str, Any]:
    if not rdkit_available():
        return {"smiles": smiles, "error": "RDKit is not available"}

    mol = parse_smiles(smiles)
    if mol is None:
        return {"smiles": smiles, "error": "Invalid SMILES"}

    compiled_path = resolve_data_path(Path(compiled_features_path))
    groups_path_resolved = resolve_data_path(Path(groups_path))

    compiled = load_compiled_features(compiled_path)
    groups_doc = load_groups(groups_path_resolved)

    sites = find_ar_anchor_sites(mol, compiled_features=compiled, groups_doc=groups_doc)

    try:
        from rdkit.Chem import rdPartialCharges  # type: ignore
    except Exception as e:
        return {"smiles": smiles, "error": f"RDKit partial charges unavailable: {e}"}

    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
    except Exception as e:
        return {"smiles": smiles, "error": f"Gasteiger charge computation failed: {e}"}

    q_ipso_list: List[Optional[float]] = []
    q_ring_mean_list: List[Optional[float]] = []
    delta_q_list: List[Optional[float]] = []
    score_list: List[int] = []
    site_payloads: List[Dict[str, Any]] = []

    for site in sites:
        ipso_atom = mol.GetAtomWithIdx(int(site.ipso_atom_idx))
        q_ipso = _safe_float_prop(ipso_atom, "_GasteigerCharge")

        ring = choose_smallest_aromatic_ring_containing(mol, int(site.ipso_atom_idx))
        q_ring_mean: Optional[float] = None
        if ring is not None and len(ring) > 0:
            charges: List[float] = []
            for idx in ring:
                a = mol.GetAtomWithIdx(int(idx))
                q = _safe_float_prop(a, "_GasteigerCharge")
                if q is not None:
                    charges.append(float(q))
            if charges:
                q_ring_mean = float(sum(charges) / len(charges))

        delta_q: Optional[float]
        if q_ipso is None:
            delta_q = None
        elif q_ring_mean is None:
            delta_q = float(q_ipso)
        else:
            delta_q = float(q_ipso - q_ring_mean)

        score = electron_poor_score_0_10(delta_q, k=float(k))

        q_ipso_list.append(q_ipso)
        q_ring_mean_list.append(q_ring_mean)
        delta_q_list.append(delta_q)
        score_list.append(int(score))

        site_payloads.append(
            {
                "label": site.label,
                "feature_token": site.feature_token,
                "ipso_atom_idx": int(site.ipso_atom_idx),
                "q_ipso": q_ipso,
                "q_ring_mean": q_ring_mean,
                "delta_q": delta_q,
                "electron_poor_score_0_10": int(score),
            }
        )

    return {
        "smiles": smiles,
        "aryl_anchor_site_count": int(len(sites)),
        "aryl_electron_poor_score_0_10_list": score_list,
        "aryl_electron_poor_score_0_10_max": int(max(score_list) if score_list else 0),
        "aryl_ipso_gasteiger_charge_list": q_ipso_list,
        "aryl_ring_mean_gasteiger_charge_list": q_ring_mean_list,
        "aryl_delta_gasteiger_charge_list": delta_q_list,
        "sites": site_payloads,
        "params": {"k": float(k)},
    }

