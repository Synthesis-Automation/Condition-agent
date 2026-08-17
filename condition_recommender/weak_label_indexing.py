"""Separate immutable index for unverified weak-label condition observations."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple


@dataclass(frozen=True)
class WeakLabelParticipant:
    """One normalized source-asserted reactive participant."""

    signature: str
    center_class: Optional[str]
    attachment_class: Optional[str]
    alpha_branched: Optional[bool]


@dataclass(frozen=True)
class WeakLabelIndexedObservation:
    """One label-only observation linked to a canonical condition recipe."""

    source_row_number: int
    reaction_type: str
    participants: Tuple[WeakLabelParticipant, WeakLabelParticipant]
    yield_pct: float
    z_score: Optional[float]
    recipe_id: str
    resolved_recipe: Dict[str, Any]


@dataclass(frozen=True)
class WeakLabelIndex:
    """Rows and deterministic source-reaction-type lookup positions."""

    rows: Tuple[WeakLabelIndexedObservation, ...]
    by_reaction_type: Mapping[str, Tuple[int, ...]]
    source_path: str

    def select_types(
        self,
        reaction_types: Tuple[str, ...],
    ) -> Tuple[WeakLabelIndexedObservation, ...]:
        positions = sorted(
            {
                position
                for reaction_type in reaction_types
                for position in self.by_reaction_type.get(reaction_type, ())
            }
        )
        return tuple(self.rows[position] for position in positions)


def weak_label_recipe_catalog_path(path: str | Path) -> Path:
    """Return the canonical recipe catalog paired with a cleaned CSV."""
    csv_path = Path(path)
    return csv_path.with_name(f"{csv_path.stem}.condition_recipes.jsonl.gz")


def _optional_bool(value: Any) -> Optional[bool]:
    text = str(value or "").strip().casefold()
    if text == "true":
        return True
    if text == "false":
        return False
    return None


def _participant(row: Mapping[str, str], prefix: str) -> WeakLabelParticipant:
    return WeakLabelParticipant(
        signature=str(row.get(f"{prefix}_signature") or "").strip(),
        center_class=(
            str(row.get(f"{prefix}_center_class") or "").strip() or None
        ),
        attachment_class=(
            str(row.get(f"{prefix}_attachment_class") or "").strip() or None
        ),
        alpha_branched=_optional_bool(row.get(f"{prefix}_alpha_branched")),
    )


def load_weak_label_index(path: str | Path) -> WeakLabelIndex:
    """Load a file-version-aware weak-label index and recipe catalog."""
    csv_path = Path(path).resolve()
    catalog_path = weak_label_recipe_catalog_path(csv_path)
    if not csv_path.is_file():
        raise ValueError(f"weak-label CSV does not exist: {csv_path}")
    if not catalog_path.is_file():
        raise ValueError(f"weak-label recipe catalog does not exist: {catalog_path}")
    return _load_weak_label_index_cached(
        str(csv_path),
        csv_path.stat().st_mtime_ns,
        str(catalog_path),
        catalog_path.stat().st_mtime_ns,
    )


@lru_cache(maxsize=4)
def _load_weak_label_index_cached(
    csv_path_text: str,
    csv_mtime_ns: int,
    catalog_path_text: str,
    catalog_mtime_ns: int,
) -> WeakLabelIndex:
    del csv_mtime_ns, catalog_mtime_ns
    recipes: Dict[str, Dict[str, Any]] = {}
    with gzip.open(catalog_path_text, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                recipe = dict(json.loads(line))
            except (TypeError, json.JSONDecodeError) as exc:
                raise ValueError(
                    f"invalid weak-label recipe JSON at line {line_number}"
                ) from exc
            recipe_id = str(recipe.get("recipe_id") or "")
            if not recipe_id.startswith("RCR1:") or recipe_id in recipes:
                raise ValueError(
                    f"invalid or duplicate weak-label recipe ID: {recipe_id}"
                )
            recipes[recipe_id] = recipe

    rows = []
    positions: Dict[str, list[int]] = {}
    with Path(csv_path_text).open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        reader = csv.DictReader(handle)
        required = {
            "source_reaction_type",
            "reactive_site_1_signature",
            "reactive_site_2_signature",
            "yield_pct",
            "condition_recipe_id",
        }
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"weak-label CSV is missing columns: {missing}")
        for source_row_number, raw in enumerate(reader, start=2):
            try:
                yield_pct = float(raw.get("yield_pct") or "")
            except (TypeError, ValueError):
                continue
            if not 0.0 <= yield_pct <= 100.0:
                continue
            try:
                z_score = float(raw.get("z_score") or "")
            except (TypeError, ValueError):
                z_score = None
            recipe_id = str(raw.get("condition_recipe_id") or "").strip()
            recipe = recipes.get(recipe_id)
            if recipe is None:
                raise ValueError(
                    f"weak-label recipe missing at CSV row {source_row_number}: "
                    f"{recipe_id}"
                )
            reaction_type = str(raw.get("source_reaction_type") or "").strip()
            row = WeakLabelIndexedObservation(
                source_row_number=source_row_number,
                reaction_type=reaction_type,
                participants=(
                    _participant(raw, "reactive_site_1"),
                    _participant(raw, "reactive_site_2"),
                ),
                yield_pct=yield_pct,
                z_score=z_score,
                recipe_id=recipe_id,
                resolved_recipe=recipe,
            )
            positions.setdefault(reaction_type, []).append(len(rows))
            rows.append(row)
    return WeakLabelIndex(
        rows=tuple(rows),
        by_reaction_type={
            key: tuple(value) for key, value in sorted(positions.items())
        },
        source_path=csv_path_text,
    )


__all__ = [
    "WeakLabelIndex",
    "WeakLabelIndexedObservation",
    "WeakLabelParticipant",
    "load_weak_label_index",
    "weak_label_recipe_catalog_path",
]
