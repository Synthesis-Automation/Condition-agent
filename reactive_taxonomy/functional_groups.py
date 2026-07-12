"""Taxonomy-driven detection of general substrate functional groups."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .chemistry.smarts_cache import compile_smarts
from .models import FunctionalGroup

_PATH = Path(__file__).with_name("definitions") / "functional_groups.v1.json"


@lru_cache(maxsize=1)
def load_functional_group_definitions() -> Tuple[Dict[str, Any], ...]:
    """Load ordered functional-group definitions."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return tuple(sorted(
        (dict(item) for item in payload.get("groups", [])),
        key=lambda item: (-int(item.get("priority", 0)), str(item.get("id", ""))),
    ))


def detect_functional_groups(mol: Any, component_index: int) -> List[FunctionalGroup]:
    """Detect deduplicated functional groups in one molecular component."""
    found: List[FunctionalGroup] = []
    seen = set()
    for definition in load_functional_group_definitions():
        query = compile_smarts(str(definition.get("smarts") or ""), validate=False)
        if query is None:
            continue
        map_positions = {
            atom.GetAtomMapNum(): atom.GetIdx()
            for atom in query.GetAtoms() if atom.GetAtomMapNum()
        }
        for match in mol.GetSubstructMatches(query, uniquify=True):
            atoms = tuple(sorted(int(index) for index in match))
            key = (str(definition["id"]), atoms)
            if key in seen:
                continue
            seen.add(key)
            anchor_position = map_positions.get(1, 0)
            found.append(FunctionalGroup(
                group_id=str(definition["id"]),
                chemist_label=str(definition.get("label") or definition["id"]),
                component_index=component_index,
                atom_indices=atoms,
                anchor_atom_index=int(match[anchor_position]),
                tags=tuple(str(tag) for tag in definition.get("tags", [])),
                matched_pattern=str(definition.get("smarts") or ""),
            ))
    return found


__all__ = ["detect_functional_groups", "load_functional_group_definitions"]
