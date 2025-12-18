#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ar- steric analysis (v2): ortho substitution + topological bulk proxy -> 0–10 score.

Implements the behavior described in `CODEX_IMPL_ReactivitySterics_v2.md`.

Key idea:
- Detect Ar-* "anchor sites" from templated atomic features whose source.groups[0] == "Ar".
- Each templated SMARTS contains atom-map label `:1` on the Ar ipso atom.
- For each anchor site, pick smallest aromatic ring containing ipso, find ortho atoms, then:
  - ortho substitution count: number of ortho atoms with total H count == 0
  - ortho bulk proxy: count heavy atoms in non-ring branches (BFS) within 4 bonds
  - steric score mapping to 0..10
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import json
from pathlib import Path
import sys
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


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


def _read_json_object(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object in {path}, got {type(payload).__name__}")
    return payload


def resolve_data_path(path: Path, *, base_dir: Path = DEFAULT_DATA_DIR) -> Path:
    candidate = path
    if not candidate.is_absolute() and not candidate.exists():
        candidate = base_dir / candidate
    return candidate


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


def _is_aromatic_ring(mol: Any, ring: Sequence[int]) -> bool:
    return all(mol.GetAtomWithIdx(int(i)).GetIsAromatic() for i in ring)


def _choose_smallest_aromatic_ring_containing(mol: Any, atom_idx: int) -> Optional[Tuple[int, ...]]:
    rings = mol.GetRingInfo().AtomRings()
    candidates = [r for r in rings if atom_idx in r and _is_aromatic_ring(mol, r)]
    if not candidates:
        return None
    candidates.sort(key=len)
    return tuple(int(i) for i in candidates[0])


def _ortho_atoms_in_ring(mol: Any, ring: Sequence[int], ipso_idx: int) -> List[int]:
    ring_set = set(int(i) for i in ring)
    ipso = mol.GetAtomWithIdx(int(ipso_idx))
    ortho: List[int] = []
    for nbr in ipso.GetNeighbors():
        j = int(nbr.GetIdx())
        if j in ring_set:
            ortho.append(j)
    ortho.sort()
    return ortho


def ortho_substitution_count(mol: Any, ipso_idx: int) -> int:
    ring = _choose_smallest_aromatic_ring_containing(mol, ipso_idx)
    if ring is None:
        return 0
    ortho_atoms = _ortho_atoms_in_ring(mol, ring, ipso_idx)
    count = 0
    for j in ortho_atoms:
        a = mol.GetAtomWithIdx(int(j))
        if a.GetTotalNumHs() == 0:
            count += 1
    return int(count)


def _count_branch_heavy_atoms(
    mol: Any,
    *,
    branch_start_idx: int,
    ring_set: set[int],
    max_depth: int,
    max_atoms: int,
) -> int:
    visited: set[int] = set()
    q: deque[Tuple[int, int]] = deque()

    visited.add(int(branch_start_idx))
    q.append((int(branch_start_idx), 0))

    heavy = 0
    if mol.GetAtomWithIdx(int(branch_start_idx)).GetAtomicNum() > 1:
        heavy += 1

    while q:
        idx, depth = q.popleft()
        if depth >= max_depth:
            continue

        atom = mol.GetAtomWithIdx(int(idx))
        for nbr in atom.GetNeighbors():
            j = int(nbr.GetIdx())
            if j in ring_set:
                continue
            if j in visited:
                continue

            visited.add(j)
            if len(visited) >= max_atoms:
                return int(heavy)

            if nbr.GetAtomicNum() > 1:
                heavy += 1
            q.append((j, depth + 1))

    return int(heavy)


def ortho_bulk_heavy_pair(
    mol: Any,
    *,
    ipso_idx: int,
    max_depth: int = 4,
    max_atoms: int = 50,
    blocked_value: int = 2,
) -> List[int]:
    ring = _choose_smallest_aromatic_ring_containing(mol, ipso_idx)
    if ring is None:
        return [0, 0]

    ring_set = set(int(i) for i in ring)
    ortho_atoms = _ortho_atoms_in_ring(mol, ring, ipso_idx)

    values: List[int] = []
    for o_idx in ortho_atoms:
        o_atom = mol.GetAtomWithIdx(int(o_idx))
        branch_starts = [int(n.GetIdx()) for n in o_atom.GetNeighbors() if int(n.GetIdx()) not in ring_set]
        if branch_starts:
            counts = [
                _count_branch_heavy_atoms(
                    mol,
                    branch_start_idx=s,
                    ring_set=ring_set,
                    max_depth=max_depth,
                    max_atoms=max_atoms,
                )
                for s in branch_starts
            ]
            values.append(int(max(counts) if counts else 0))
        else:
            values.append(int(blocked_value if o_atom.GetTotalNumHs() == 0 else 0))

    if len(values) >= 2:
        return [int(values[0]), int(values[1])]
    if len(values) == 1:
        return [int(values[0]), 0]
    return [0, 0]


def ortho_bulk_score(ortho_bulk_heavy: Sequence[int], *, cap_each: int = 4) -> int:
    h1 = int(ortho_bulk_heavy[0]) if len(ortho_bulk_heavy) > 0 else 0
    h2 = int(ortho_bulk_heavy[1]) if len(ortho_bulk_heavy) > 1 else 0
    return int(min(cap_each, h1) + min(cap_each, h2))


def _round_half_up(x: float) -> int:
    return int(x + 0.5)


def steric_score_0_10_from_bulk_score(bulk_score: int) -> int:
    b = int(bulk_score)
    if b <= 0:
        return 0
    v = _round_half_up(1 + 9 * b / 8)
    return int(max(0, min(10, v)))


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


def analyze_smiles_reactivity_sterics_v2(
    smiles: str,
    *,
    compiled_features_path: Path = Path("calculable_features.compiled.v1.json"),
    groups_path: Path = Path("organic_groups.v1.json"),
    max_depth: int = 4,
    max_atoms: int = 50,
    blocked_value: int = 2,
    cap_each: int = 4,
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

    ortho_sub_list: List[int] = []
    ortho_bulk_heavy_list: List[List[int]] = []
    ortho_bulk_score_list: List[int] = []
    steric_score_list: List[int] = []
    site_payloads: List[Dict[str, Any]] = []

    for site in sites:
        sub_count = ortho_substitution_count(mol, site.ipso_atom_idx)
        heavy_pair = ortho_bulk_heavy_pair(
            mol,
            ipso_idx=site.ipso_atom_idx,
            max_depth=max_depth,
            max_atoms=max_atoms,
            blocked_value=blocked_value,
        )
        bulk = ortho_bulk_score(heavy_pair, cap_each=cap_each)
        steric = steric_score_0_10_from_bulk_score(bulk)

        ortho_sub_list.append(int(sub_count))
        ortho_bulk_heavy_list.append([int(heavy_pair[0]), int(heavy_pair[1])])
        ortho_bulk_score_list.append(int(bulk))
        steric_score_list.append(int(steric))

        site_payloads.append(
            {
                "label": site.label,
                "feature_token": site.feature_token,
                "ipso_atom_idx": int(site.ipso_atom_idx),
                "ortho_sub_count": int(sub_count),
                "ortho_bulk_heavy": [int(heavy_pair[0]), int(heavy_pair[1])],
                "ortho_bulk_score": int(bulk),
                "steric_score_0_10": int(steric),
            }
        )

    return {
        "smiles": smiles,
        "aryl_anchor_site_count": int(len(sites)),
        "aryl_ortho_sub_count_list": ortho_sub_list,
        "aryl_ortho_bulk_heavy_list": ortho_bulk_heavy_list,
        "aryl_ortho_bulk_score_list": ortho_bulk_score_list,
        "aryl_steric_score_0_10_list": steric_score_list,
        "aryl_steric_score_0_10_max": int(max(steric_score_list) if steric_score_list else 0),
        "sites": site_payloads,
    }

