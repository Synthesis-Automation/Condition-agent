"""Index the cleaned weak-label dataset for condition retrieval."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from .label_conditions import condition_recipe_catalog_path

CONDITION_COLUMNS = (
    "condition_display",
    "temperature_c",
    "time_h",
    "concentration_m",
    "atmosphere",
    "condition_identity_uncertainty",
    "procedure_text",
)


@dataclass(frozen=True)
class LabelParticipant:
    signature: str
    center_class: Optional[str]
    attachment_class: Optional[str]
    alpha_branched: Optional[bool]


@dataclass(frozen=True)
class LabelIndexedReaction:
    source_row_number: int
    reaction_type: str
    participants: Tuple[LabelParticipant, LabelParticipant]
    yield_pct: float
    z_score: float
    conditions: Dict[str, str]
    recipe_id: str
    resolved_recipe: Dict[str, Any]


def load_condition_recipe_catalog(
    path: str | Path,
) -> Dict[str, Dict[str, Any]]:
    """Load and validate the deduplicated nested condition recipes."""
    catalog_path = Path(path)
    if not catalog_path.exists():
        raise ValueError(f"Missing condition recipe catalog: {catalog_path}")
    recipes: Dict[str, Dict[str, Any]] = {}
    with gzip.open(catalog_path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                recipe = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid condition recipe JSONL line {line_number}"
                ) from exc
            recipe_id = str(recipe.get("recipe_id") or "")
            if not recipe_id.startswith("RCR1:"):
                raise ValueError(
                    f"Invalid condition recipe ID at JSONL line {line_number}"
                )
            if recipe_id in recipes:
                raise ValueError(f"Duplicate condition recipe ID: {recipe_id}")
            recipes[recipe_id] = recipe
    return recipes


def _optional_bool(value: str) -> Optional[bool]:
    text = str(value or "").strip().lower()
    if text == "true":
        return True
    if text == "false":
        return False
    return None


def _participant(row: Dict[str, str], prefix: str) -> LabelParticipant:
    return LabelParticipant(
        signature=str(row.get(f"{prefix}_signature") or "").strip(),
        center_class=str(row.get(f"{prefix}_center_class") or "").strip() or None,
        attachment_class=str(row.get(f"{prefix}_attachment_class") or "").strip() or None,
        alpha_branched=_optional_bool(str(row.get(f"{prefix}_alpha_branched") or "")),
    )


def load_label_index(path: str | Path) -> Tuple[LabelIndexedReaction, ...]:
    """Load a file-version-aware cached weak-label index."""
    csv_path = Path(path).resolve()
    catalog_path = condition_recipe_catalog_path(csv_path)
    return _load_label_index_cached(
        str(csv_path),
        csv_path.stat().st_mtime_ns,
        str(catalog_path),
        catalog_path.stat().st_mtime_ns,
    )


@lru_cache(maxsize=4)
def _load_label_index_cached(
    csv_path_text: str,
    csv_mtime_ns: int,
    catalog_path_text: str,
    catalog_mtime_ns: int,
) -> Tuple[LabelIndexedReaction, ...]:
    """Load one immutable dataset version; mtimes form part of the cache key."""
    del csv_mtime_ns, catalog_mtime_ns
    csv_path = Path(csv_path_text)
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "yield_pct",
            "source_reaction_type",
            "reactive_site_1_signature",
            "reactive_site_2_signature",
            "condition_recipe_id",
        }
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"Missing cleaned label columns: {missing}")
        rows = list(reader)
    recipe_catalog = load_condition_recipe_catalog(
        catalog_path_text
    )

    indexed = []
    for row_number, row in enumerate(rows, start=2):
        try:
            yield_pct = float(row.get("yield_pct") or "")
        except (TypeError, ValueError):
            continue
        if not 0.0 <= yield_pct <= 100.0:
            continue
        try:
            z_score = float(row.get("z_score") or 0.0)
        except (TypeError, ValueError):
            z_score = 0.0
        conditions = {name: str(row.get(name) or "").strip() for name in CONDITION_COLUMNS}
        recipe_id = str(row.get("condition_recipe_id") or "").strip()
        if not recipe_id.startswith("RCR1:"):
            raise ValueError(
                f"Invalid canonical condition recipe ID at CSV row {row_number}"
            )
        resolved_recipe = recipe_catalog.get(recipe_id)
        if resolved_recipe is None:
            raise ValueError(
                f"Condition recipe ID missing from catalog at CSV row {row_number}"
            )
        indexed.append(
            LabelIndexedReaction(
                source_row_number=row_number,
                reaction_type=str(row.get("source_reaction_type") or "").strip(),
                participants=(
                    _participant(row, "reactive_site_1"),
                    _participant(row, "reactive_site_2"),
                ),
                yield_pct=yield_pct,
                z_score=z_score,
                conditions=conditions,
                recipe_id=recipe_id,
                resolved_recipe=resolved_recipe,
            )
        )
    return tuple(indexed)


__all__ = [
    "CONDITION_COLUMNS",
    "LabelIndexedReaction",
    "LabelParticipant",
    "condition_recipe_catalog_path",
    "load_condition_recipe_catalog",
    "load_label_index",
]
