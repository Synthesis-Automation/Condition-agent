#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ar- steric analysis (POC): ortho substitution counts around Ar-* anchor sites.

Key idea:
- We detect Ar-* "anchor sites" from **templated atomic features** whose source.groups[0] == "Ar".
- Each such feature SMARTS contains an atom-map label :1 on the Ar- ipso atom (e.g., [c:1][Br], [c:1][CX3H1](=O)).
- Using RDKit substructure matches, we map query atom-map :1 to the molecule atom index ("ipso").
- Then we compute ortho substitution count (0/1/2) around that ipso atom.

Definition (POC):
- Identify an aromatic ring containing the ipso atom (choose the smallest aromatic ring containing it).
- The two neighboring atoms of ipso within that ring are the ortho positions.
- An ortho position is counted as "substituted" if its total H count == 0.
  (This counts actual substituents, fused-ring blocking, and aromatic hetero atoms as "blocked".)

This keeps code reusable across Ar-Br / Ar-CHO / Ar-OTf / Ar-CN / etc.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple
import json
from pathlib import Path
import sys


def _ensure_repo_root_on_syspath() -> None:
    """
    Allow running this module as a script from within `chemtools/taxonomy/v1/`.

    When executed directly (e.g. `python ar_context_sterics_v1.py`), Python's
    `sys.path[0]` is this folder, so `import chemtools...` would fail unless the
    repo root is also on `sys.path`.
    """
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
    token: str
    label: str
    ipso_idx: int
    match_atoms: Tuple[int, ...]
    smarts: str

def _read_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object in {path}, got {type(payload).__name__}")
    return payload

def load_compiled_features(path: Path) -> Dict[str, Any]:
    return _read_json(path)

def load_groups(path: Path) -> Dict[str, Any]:
    return _read_json(path)


def resolve_data_path(path: Path, *, base_dir: Path = DEFAULT_DATA_DIR) -> Path:
    """
    Resolve a data file path robustly.

    If `path` is relative and doesn't exist in the current working directory,
    try relative to `base_dir` (this module's folder) and `base_dir/specs/`.
    """
    if path.is_absolute() or path.exists():
        return path

    candidates = [
        base_dir / path,
        base_dir / "specs" / path,
        base_dir / "specs" / path.name,
    ]
    for c in candidates:
        if c.exists():
            return c
    return candidates[-1]

def _compose_label(left: str, right: str) -> str:
    if left.endswith("-") and right.startswith("-"):
        return left + right[1:]
    return left + right

def _find_anchor_qidx_by_mapnum(query_mol, map_num: int) -> int:
    for i, a in enumerate(query_mol.GetAtoms()):
        if a.GetAtomMapNum() == map_num:
            return i
    raise ValueError(f"No atom with AtomMapNum=={map_num} in query SMARTS")

def _is_aromatic_ring(mol, ring: Tuple[int, ...]) -> bool:
    # Consider ring aromatic if all atoms are aromatic
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)

def _choose_smallest_aromatic_ring_containing(mol, atom_idx: int) -> Optional[Tuple[int, ...]]:
    rings = mol.GetRingInfo().AtomRings()
    candidates = [r for r in rings if atom_idx in r and _is_aromatic_ring(mol, r)]
    if not candidates:
        return None
    candidates.sort(key=len)
    return tuple(candidates[0])

def _ortho_atoms_in_ring(mol, ring: Tuple[int, ...], ipso_idx: int) -> List[int]:
    ring_set = set(ring)
    ipso = mol.GetAtomWithIdx(ipso_idx)
    ortho = []
    for nbr in ipso.GetNeighbors():
        j = nbr.GetIdx()
        if j in ring_set:
            ortho.append(j)
    # For simple aromatic rings, should be 2
    return ortho

def ortho_substitution_count(mol, ipso_idx: int) -> Optional[int]:
    ring = _choose_smallest_aromatic_ring_containing(mol, ipso_idx)
    if ring is None:
        return None
    ortho_atoms = _ortho_atoms_in_ring(mol, ring, ipso_idx)
    if not ortho_atoms:
        return None
    count = 0
    for j in ortho_atoms:
        a = mol.GetAtomWithIdx(j)
        # POC definition: ortho is blocked if it has no hydrogens
        if a.GetTotalNumHs() == 0:
            count += 1
    # Clamp to 0..2 for typical rings; still return actual count if unusual
    return count

def find_ar_anchor_sites(
    mol,
    compiled_features: Dict[str, Any],
    groups_doc: Dict[str, Any],
    context_group_id: str = "Ar",
    context_atom_map_num: int = 1,
) -> List[AnchorSite]:
    # Build group label index for composing label
    g_by_id = {g["id"]: g for g in groups_doc.get("groups", [])}
    left_label = g_by_id.get(context_group_id, {}).get("name", f"{context_group_id}-")

    sites: List[AnchorSite] = []

    atomic = compiled_features.get("atomic", [])
    for feat in atomic:
        src = feat.get("source", {})
        if src.get("kind") != "templated":
            continue
        gids = src.get("groups", []) or []
        if not gids or gids[0] != context_group_id:
            continue

        token = feat["token"]
        detect = feat.get("detect", {})
        smarts_list = detect.get("smarts_any", []) or []
        if not smarts_list:
            continue

        right_id = gids[1] if len(gids) > 1 else ""
        right_label = g_by_id.get(right_id, {}).get("name", right_id)
        label = _compose_label(left_label, right_label) if right_label else left_label

        for smarts in smarts_list:
            if not isinstance(smarts, str) or not smarts.strip():
                continue
            qmol = compile_smarts(smarts, validate=True)
            if qmol is None:
                raise ValueError(f"SMARTS compiled to None for {token}: {smarts}")

            anchor_qidx = _find_anchor_qidx_by_mapnum(qmol, context_atom_map_num)

            for match in mol.GetSubstructMatches(qmol):
                ipso_idx = match[anchor_qidx]
                sites.append(
                    AnchorSite(
                        token=token,
                        label=label,
                        ipso_idx=ipso_idx,
                        match_atoms=tuple(match),
                        smarts=smarts,
                    )
                )

    # De-duplicate sites by (ipso_idx, token) to avoid double counting if multiple matches identical
    seen = set()
    dedup = []
    for s in sites:
        key = (s.ipso_idx, s.token)
        if key in seen:
            continue
        seen.add(key)
        dedup.append(s)
    return dedup

def analyze_smiles_ortho_counts(
    smiles: str,
    compiled_features_path: Path,
    groups_path: Path,
) -> Dict[str, Any]:
    if not rdkit_available():
        return {"smiles": smiles, "error": "RDKit is not available"}

    mol = parse_smiles(smiles)
    if mol is None:
        return {"smiles": smiles, "error": "Invalid SMILES"}

    compiled_path = resolve_data_path(compiled_features_path)
    groups_path_resolved = resolve_data_path(groups_path)

    compiled = load_compiled_features(compiled_path)
    groups_doc = load_groups(groups_path_resolved)

    sites = find_ar_anchor_sites(mol, compiled, groups_doc, context_group_id="Ar", context_atom_map_num=1)

    site_results = []
    ortho_list: List[int] = []
    for s in sites:
        c = ortho_substitution_count(mol, s.ipso_idx)
        ortho_list.append(int(c) if c is not None else 0)
        site_results.append({
            "label": s.label,
            "feature_token": s.token,
            "ipso_atom_idx": s.ipso_idx,
            "ortho_sub_count": int(c) if c is not None else 0,
        })

    out: Dict[str, Any] = {
        "smiles": smiles,
        "aryl_anchor_site_count": len(site_results),
        "aryl_ortho_sub_count_list": ortho_list,
        "aryl_ortho_sub_count_max": max(ortho_list) if ortho_list else 0,
        "aryl_ortho_sub_count_sum": sum(ortho_list) if ortho_list else 0,
        "sites": site_results,
    }
    return out
