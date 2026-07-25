"""Taxonomy-owned SMARTS loading and mapped candidate extraction."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set, Tuple

from .chemistry.smarts_cache import compile_smarts


_HANDLES_PATH = Path(__file__).with_name("definitions") / "handles.v1.json"
_CONTEXTS_PATH = Path(__file__).with_name("definitions") / "contexts.v1.json"


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


class MatchIndex:
    """All taxonomy SMARTS matches for one molecule, calculated once."""

    def __init__(self, mol: Any) -> None:
        self.mol = mol
        self._handle_matches: List[Tuple[Dict[str, Any], Dict[int, int], Tuple[Tuple[int, ...], ...]]] = []
        self._context_matches: Dict[str, Tuple[Dict[int, int], Tuple[Tuple[int, ...], ...]]] = {}
        for definition in load_handle_patterns():
            query = compile_smarts(str(definition.get("smarts") or ""), validate=False)
            if query is None:
                continue
            matches = tuple(
                tuple(int(i) for i in match)
                for match in mol.GetSubstructMatches(
                    query,
                    uniquify=not bool(definition.get("activation_token")),
                )
            )
            self._handle_matches.append((definition, _query_map_positions(query), matches))
        with _CONTEXTS_PATH.open("r", encoding="utf-8") as handle:
            contexts = json.load(handle).get("contexts") or []
        for definition in contexts:
            if definition.get("classification_method") != "mapped_smarts":
                continue
            query = compile_smarts(str(definition.get("smarts") or ""), validate=False)
            if query is None:
                continue
            matches = tuple(tuple(int(i) for i in match) for match in mol.GetSubstructMatches(query, uniquify=True))
            self._context_matches[str(definition["id"])] = (_query_map_positions(query), matches)

    def role_atoms(self, site_type: str, role: str) -> Set[int]:
        indices: Set[int] = set()
        for definition, positions, matches in self._handle_matches:
            if definition.get("site_type") != site_type:
                continue
            raw_maps = (definition.get("atom_roles") or {}).get(role)
            if raw_maps is None:
                continue
            role_maps = raw_maps if isinstance(raw_maps, list) else [raw_maps]
            for match in matches:
                for map_number in role_maps:
                    position = positions.get(int(map_number))
                    if position is not None:
                        indices.add(match[position])
        return indices

    def patterns_for_atom(self, site_type: str, role: str, atom_index: int) -> List[Dict[str, Any]]:
        found: List[Dict[str, Any]] = []
        for definition, positions, matches in self._handle_matches:
            if definition.get("site_type") != site_type:
                continue
            raw_maps = (definition.get("atom_roles") or {}).get(role)
            if raw_maps is None:
                continue
            role_maps = raw_maps if isinstance(raw_maps, list) else [raw_maps]
            role_positions = [positions[int(value)] for value in role_maps if int(value) in positions]
            if any(any(match[position] == atom_index for position in role_positions) for match in matches):
                found.append(definition)
        return found

    def role_matches_for_atom(
        self,
        site_type: str,
        role: str,
        atom_index: int,
    ) -> List[Tuple[Dict[str, Any], Dict[str, Tuple[int, ...]]]]:
        """Return each pattern match assigning ``role`` to ``atom_index``.

        Unlike :meth:`patterns_for_atom`, this preserves distinct matches of the
        same definition. That is required when one carbon is activated by two
        separate electron-withdrawing groups.
        """
        found: List[Tuple[Dict[str, Any], Dict[str, Tuple[int, ...]]]] = []
        for definition, positions, matches in self._handle_matches:
            if definition.get("site_type") != site_type:
                continue
            raw_target_maps = (definition.get("atom_roles") or {}).get(role)
            if raw_target_maps is None:
                continue
            target_maps = (
                raw_target_maps
                if isinstance(raw_target_maps, list)
                else [raw_target_maps]
            )
            target_positions = [
                positions[int(value)]
                for value in target_maps
                if int(value) in positions
            ]
            for match in matches:
                if not any(match[position] == atom_index for position in target_positions):
                    continue
                assignments: Dict[str, Tuple[int, ...]] = {}
                for role_name, raw_maps in (
                    definition.get("atom_roles") or {}
                ).items():
                    role_maps = raw_maps if isinstance(raw_maps, list) else [raw_maps]
                    assignments[str(role_name)] = tuple(
                        int(match[positions[int(value)]])
                        for value in role_maps
                        if int(value) in positions
                    )
                found.append((definition, assignments))
        return found

    def context_match(
        self,
        definition: Dict[str, Any],
        atom_index: int,
        excluded: Set[int],
    ) -> Tuple[int, ...] | None:
        indexed = self._context_matches.get(str(definition["id"]))
        if indexed is None:
            return None
        positions, matches = indexed
        roles = definition.get("atom_roles") or {}
        anchor_position = positions.get(int(roles.get("context_anchor", -1)))
        if anchor_position is None:
            return None
        raw_substituents = roles.get("substituent") or []
        if not isinstance(raw_substituents, list):
            raw_substituents = [raw_substituents]
        substituent_positions = [positions[int(value)] for value in raw_substituents if int(value) in positions]
        return next((
            match for match in matches
            if match[anchor_position] == atom_index
            and not any(match[position] in excluded for position in substituent_positions)
        ), None)


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


__all__ = ["MatchIndex", "load_handle_patterns", "matched_pattern_ids", "matched_patterns_for_atom", "matched_role_atoms", "patterns_for"]
