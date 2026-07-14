"""In-memory indices for canonical generic conversion records."""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Tuple


@dataclass(frozen=True)
class GenericIndexedReaction:
    """A verified precedent reduced to fields required for retrieval."""

    reaction_id: str
    observation_id: str
    reaction_smiles: str
    yield_pct: float
    source_dataset: str
    signature: Dict[str, Any]
    recipe_id: str
    resolved_recipe: Dict[str, Any]
    condition_uncertain: bool

    @property
    def named_family(self) -> str:
        return str(self.signature.get("named_family") or "")

    @property
    def family_confidence(self) -> float:
        return float(self.signature.get("family_confidence") or 0.0)

    @property
    def transformation_class(self) -> str:
        return str(self.signature.get("transformation_class") or "")


@dataclass(frozen=True)
class GenericReactionIndex:
    """Immutable precedent collection with deterministic signature lookup maps."""

    rows: Tuple[GenericIndexedReaction, ...]
    exact: Mapping[str, Tuple[int, ...]]
    handles: Mapping[str, Tuple[int, ...]]
    transformations: Mapping[str, Tuple[int, ...]]
    bond_edits: Mapping[str, Tuple[int, ...]]
    environments: Mapping[str, Tuple[int, ...]]
    families: Mapping[str, Tuple[int, ...]]

    def select(self, positions: Iterable[int]) -> Tuple[GenericIndexedReaction, ...]:
        return tuple(self.rows[position] for position in positions)


_KEY_FIELDS = {
    "exact": "exact_signature_key",
    "handles": "handle_signature_key",
    "transformations": "transformation_signature_key",
    "bond_edits": "bond_edit_signature_key",
    "environments": "environment_signature_key",
}


def _freeze(mapping: Dict[str, list[int]]) -> Dict[str, Tuple[int, ...]]:
    return {key: tuple(values) for key, values in sorted(mapping.items())}


def build_generic_index(
    records: Iterable[Mapping[str, Any]],
    *,
    include_review: bool = False,
) -> GenericReactionIndex:
    """Build lookup maps, admitting only records with signatures and recipes."""
    rows = []
    for record in records:
        tier = str(record.get("admission_tier") or "")
        if tier != "verified" and not (include_review and tier == "review"):
            continue
        signature = record.get("reaction_signature")
        recipe = record.get("resolved_recipe")
        recipe_id = str(record.get("resolved_recipe_id") or "")
        outcome = record.get("yield_pct")
        if not isinstance(signature, Mapping) or not isinstance(recipe, Mapping):
            continue
        if not recipe_id or outcome is None:
            continue
        yield_pct = float(outcome)
        if not 0.0 <= yield_pct <= 100.0:
            continue
        if str(recipe.get("recipe_id") or "") != recipe_id:
            continue
        rows.append(
            GenericIndexedReaction(
                reaction_id=str(record.get("reaction_id") or ""),
                observation_id=str(record.get("observation_id") or ""),
                reaction_smiles=str(record.get("reaction_smiles") or ""),
                yield_pct=yield_pct,
                source_dataset=str(record.get("source_dataset") or ""),
                signature=dict(signature),
                recipe_id=recipe_id,
                resolved_recipe=dict(recipe),
                condition_uncertain=bool(
                    (record.get("condition_resolution") or {}).get("has_uncertainty")
                ),
            )
        )
    rows.sort(key=lambda row: (row.reaction_id, row.observation_id, row.recipe_id))
    maps: Dict[str, Dict[str, list[int]]] = {
        name: defaultdict(list) for name in _KEY_FIELDS
    }
    families: Dict[str, list[int]] = defaultdict(list)
    for position, row in enumerate(rows):
        for name, field in _KEY_FIELDS.items():
            key = str(row.signature.get(field) or "")
            if key:
                maps[name][key].append(position)
        if row.named_family:
            families[row.named_family].append(position)
    return GenericReactionIndex(
        rows=tuple(rows),
        exact=_freeze(maps["exact"]),
        handles=_freeze(maps["handles"]),
        transformations=_freeze(maps["transformations"]),
        bond_edits=_freeze(maps["bond_edits"]),
        environments=_freeze(maps["environments"]),
        families=_freeze(families),
    )


def load_generic_index(
    path: str | Path,
    *,
    include_review: bool = False,
) -> GenericReactionIndex:
    """Load canonical JSONL output from the generic conversion engine."""
    records = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSONL at line {line_number}: {exc.msg}") from exc
            if not isinstance(value, dict):
                raise ValueError(f"JSONL line {line_number} is not an object")
            records.append(value)
    return build_generic_index(records, include_review=include_review)


__all__ = [
    "GenericIndexedReaction",
    "GenericReactionIndex",
    "build_generic_index",
    "load_generic_index",
]
