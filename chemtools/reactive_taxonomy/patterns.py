"""Taxonomy-owned SMARTS loading and mapped candidate extraction."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set

from chemtools.core.smarts import compile_smarts


_HANDLES_PATH = Path(__file__).with_name("data") / "handles.v1.json"


@lru_cache(maxsize=1)
def load_handle_patterns() -> List[Dict[str, Any]]:
    """Load independent reactive-handle SMARTS definitions."""
    with _HANDLES_PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return sorted(
        [dict(item) for item in payload.get("patterns", [])],
        key=lambda item: (-int(item.get("priority", 0)), str(item.get("id", ""))),
    )


def patterns_for(site_type: str) -> List[Dict[str, Any]]:
    return [item for item in load_handle_patterns() if item.get("site_type") == site_type]


def _query_map_positions(pattern: Any) -> Dict[int, int]:
    return {
        int(atom.GetAtomMapNum()): atom.GetIdx()
        for atom in pattern.GetAtoms()
        if atom.GetAtomMapNum()
    }


def matched_role_atoms(mol: Any, site_type: str, role: str) -> Set[int]:
    """Return molecule atom indices assigned to a mapped role by SMARTS."""
    indices: Set[int] = set()
    for definition in patterns_for(site_type):
        raw_roles = definition.get("atom_roles") or {}
        role_maps = raw_roles.get(role)
        if role_maps is None:
            continue
        if not isinstance(role_maps, list):
            role_maps = [role_maps]
        query = compile_smarts(str(definition.get("smarts") or ""), validate=False)
        if query is None:
            continue
        positions = _query_map_positions(query)
        for match in mol.GetSubstructMatches(query, uniquify=True):
            for map_number in role_maps:
                position = positions.get(int(map_number))
                if position is not None:
                    indices.add(int(match[position]))
    return indices


def matched_pattern_ids(mol: Any, site_type: str) -> Set[str]:
    """Return IDs of taxonomy patterns that match a molecule."""
    matched: Set[str] = set()
    for definition in patterns_for(site_type):
        query = compile_smarts(str(definition.get("smarts") or ""), validate=False)
        if query is not None and mol.HasSubstructMatch(query):
            matched.add(str(definition["id"]))
    return matched


def matched_patterns_for_atom(
    mol: Any,
    site_type: str,
    role: str,
    atom_index: int,
) -> List[Dict[str, Any]]:
    """Return ordered pattern definitions assigning ``role`` to an atom."""
    matched: List[Dict[str, Any]] = []
    for definition in patterns_for(site_type):
        raw_role = (definition.get("atom_roles") or {}).get(role)
        if raw_role is None:
            continue
        role_maps = raw_role if isinstance(raw_role, list) else [raw_role]
        query = compile_smarts(str(definition.get("smarts") or ""), validate=False)
        if query is None:
            continue
        positions = _query_map_positions(query)
        found = False
        for match in mol.GetSubstructMatches(query, uniquify=True):
            if any(
                positions.get(int(map_number)) is not None
                and int(match[positions[int(map_number)]]) == int(atom_index)
                for map_number in role_maps
            ):
                found = True
                break
        if found:
            matched.append(definition)
    return matched


__all__ = ["load_handle_patterns", "matched_pattern_ids", "matched_patterns_for_atom", "matched_role_atoms", "patterns_for"]
