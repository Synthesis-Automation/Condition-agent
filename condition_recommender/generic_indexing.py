"""In-memory indices for canonical generic conversion records."""

from __future__ import annotations

import json
import hashlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Tuple


@dataclass(frozen=True)
class GenericIndexedReaction:
    """A verified precedent reduced to fields required for retrieval."""

    reaction_id: str
    observation_id: str
    canonical_reaction_id: str
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


def build_generic_index_from_rows(
    rows: Iterable[GenericIndexedReaction],
) -> GenericReactionIndex:
    """Build deterministic lookup maps from already validated index rows."""
    ordered = sorted(
        rows,
        key=lambda row: (
            row.canonical_reaction_id,
            row.reaction_id,
            row.observation_id,
            row.recipe_id,
        ),
    )
    maps: Dict[str, Dict[str, list[int]]] = {
        name: defaultdict(list) for name in _KEY_FIELDS
    }
    families: Dict[str, list[int]] = defaultdict(list)
    for position, row in enumerate(ordered):
        for name, field in _KEY_FIELDS.items():
            key = str(row.signature.get(field) or "")
            if key:
                maps[name][key].append(position)
        if row.named_family:
            families[row.named_family].append(position)
    return GenericReactionIndex(
        rows=tuple(ordered),
        exact=_freeze(maps["exact"]),
        handles=_freeze(maps["handles"]),
        transformations=_freeze(maps["transformations"]),
        bond_edits=_freeze(maps["bond_edits"]),
        environments=_freeze(maps["environments"]),
        families=_freeze(families),
    )


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
                canonical_reaction_id=str(
                    record.get("canonical_reaction_id")
                    or record.get("reaction_id")
                    or record.get("observation_id")
                    or ""
                ),
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
    return build_generic_index_from_rows(rows)


def _index_payload(index: GenericReactionIndex) -> Dict[str, Any]:
    rows = [
        {
            "reaction_id": row.reaction_id,
            "observation_id": row.observation_id,
            "canonical_reaction_id": row.canonical_reaction_id,
            "reaction_smiles": row.reaction_smiles,
            "yield_pct": row.yield_pct,
            "source_dataset": row.source_dataset,
            "signature": row.signature,
            "recipe_id": row.recipe_id,
            "resolved_recipe": row.resolved_recipe,
            "condition_uncertain": row.condition_uncertain,
        }
        for row in index.rows
    ]
    maps = {
        "exact": dict(index.exact),
        "handles": dict(index.handles),
        "transformations": dict(index.transformations),
        "bond_edits": dict(index.bond_edits),
        "environments": dict(index.environments),
        "families": dict(index.families),
    }
    identity = json.dumps(
        {"rows": rows, "maps": maps},
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    )
    return {
        "schema_version": "1.0",
        "artifact_type": "generic_reaction_index",
        "index_id": "GRI1:" + hashlib.sha256(identity.encode("utf-8")).hexdigest(),
        "row_count": len(rows),
        "rows": rows,
        "maps": maps,
    }


def save_generic_index(index: GenericReactionIndex, path: str | Path) -> Dict[str, Any]:
    """Write a deterministic, versioned generic index artifact."""
    payload = _index_payload(index)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        + "\n",
        encoding="utf-8",
    )
    return {
        "schema_version": payload["schema_version"],
        "index_id": payload["index_id"],
        "row_count": payload["row_count"],
        "path": str(destination),
    }


def load_persisted_generic_index(path: str | Path) -> GenericReactionIndex:
    """Load and validate a persisted generic index without rebuilding maps."""
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("artifact_type") != "generic_reaction_index":
        raise ValueError("Not a generic reaction index artifact")
    if payload.get("schema_version") != "1.0":
        raise ValueError("Unsupported generic reaction index schema")
    rows = tuple(
        GenericIndexedReaction(
            reaction_id=str(row["reaction_id"]),
            observation_id=str(row["observation_id"]),
            canonical_reaction_id=str(row["canonical_reaction_id"]),
            reaction_smiles=str(row["reaction_smiles"]),
            yield_pct=float(row["yield_pct"]),
            source_dataset=str(row["source_dataset"]),
            signature=dict(row["signature"]),
            recipe_id=str(row["recipe_id"]),
            resolved_recipe=dict(row["resolved_recipe"]),
            condition_uncertain=bool(row["condition_uncertain"]),
        )
        for row in payload.get("rows") or ()
    )
    maps = payload.get("maps") or {}
    index = GenericReactionIndex(
        rows=rows,
        exact={key: tuple(value) for key, value in (maps.get("exact") or {}).items()},
        handles={key: tuple(value) for key, value in (maps.get("handles") or {}).items()},
        transformations={key: tuple(value) for key, value in (maps.get("transformations") or {}).items()},
        bond_edits={key: tuple(value) for key, value in (maps.get("bond_edits") or {}).items()},
        environments={key: tuple(value) for key, value in (maps.get("environments") or {}).items()},
        families={key: tuple(value) for key, value in (maps.get("families") or {}).items()},
    )
    expected = _index_payload(index)
    if expected["index_id"] != payload.get("index_id"):
        raise ValueError("Generic reaction index integrity check failed")
    if int(payload.get("row_count", -1)) != len(rows):
        raise ValueError("Generic reaction index row count mismatch")
    return index


def load_generic_index(
    path: str | Path,
    *,
    include_review: bool = False,
) -> GenericReactionIndex:
    """Load canonical JSONL output from the generic conversion engine."""
    source = Path(path)
    if source.suffix.casefold() == ".json":
        return load_persisted_generic_index(source)
    records = []
    with source.open("r", encoding="utf-8") as handle:
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
    "build_generic_index_from_rows",
    "load_generic_index",
    "load_persisted_generic_index",
    "save_generic_index",
]
