"""Definition-driven, non-exclusive molecular motif annotations."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .chemistry.smarts_cache import compile_smarts
from .models import MolecularMotifMatch

_PATH = Path(__file__).with_name("definitions") / "molecular_motifs.v1.json"


@lru_cache(maxsize=1)
def load_molecular_motif_definitions() -> Tuple[Dict[str, Any], ...]:
    """Load ordered molecular-motif definitions."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return tuple(sorted(
        (dict(item) for item in payload.get("groups", [])),
        key=lambda item: (-int(item.get("priority", 0)), str(item.get("id", ""))),
    ))


def detect_molecular_motifs(mol: Any, component_index: int) -> List[MolecularMotifMatch]:
    """Detect deduplicated molecular motifs in one component."""
    candidates: List[Tuple[Dict[str, Any], MolecularMotifMatch]] = []
    seen = set()
    for definition in load_molecular_motif_definitions():
        raw_patterns = definition.get("smarts") or ""
        patterns = raw_patterns if isinstance(raw_patterns, list) else [raw_patterns]
        for smarts in patterns:
            query = compile_smarts(str(smarts), validate=False)
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
                candidates.append((definition, MolecularMotifMatch(
                    motif_id=str(definition["id"]),
                    chemist_label=str(definition.get("label") or definition["id"]),
                    component_index=component_index,
                    atom_indices=atoms,
                    anchor_atom_index=int(match[anchor_position]),
                    tags=tuple(str(tag) for tag in definition.get("tags", [])),
                    matched_pattern=str(smarts),
                )))
    retained: List[Tuple[Dict[str, Any], MolecularMotifMatch]] = []
    for definition, candidate in candidates:
        candidate_atoms = set(candidate.atom_indices)
        suppressed = any(
            candidate.motif_id in set(owner_definition.get("suppresses_on_overlap") or [])
            and bool(candidate_atoms.intersection(owner.atom_indices))
            for owner_definition, owner in retained
        )
        if not suppressed:
            retained.append((definition, candidate))
    return [group for _, group in retained]


__all__ = ["detect_molecular_motifs", "load_molecular_motif_definitions"]
