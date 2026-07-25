"""Index the cleaned weak-label dataset for condition retrieval."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple


CONDITION_COLUMNS = (
    "base",
    "catalyst",
    "primary_solvent",
    "ligand",
    "additive",
    "coupling_reagent",
    "secondary_solvent",
    "tertiary_solvent",
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


def _recipe_id(conditions: Dict[str, str]) -> str:
    payload = json.dumps(conditions, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]
    return f"LABELCOND1:{digest}"


def load_label_index(path: str | Path) -> Tuple[LabelIndexedReaction, ...]:
    """Load cleaned rows with valid yield and the structured FG schema."""
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "yield_pct",
            "source_reaction_type",
            "reactive_site_1_signature",
            "reactive_site_2_signature",
        }
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"Missing cleaned label columns: {missing}")
        rows = list(reader)

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
                recipe_id=_recipe_id(conditions),
            )
        )
    return tuple(indexed)


__all__ = [
    "CONDITION_COLUMNS",
    "LabelIndexedReaction",
    "LabelParticipant",
    "load_label_index",
]
