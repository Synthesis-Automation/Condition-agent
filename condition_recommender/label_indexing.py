"""Index the cleaned weak-label dataset for condition retrieval."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple


CONDITION_COLUMNS = (
    "Base",
    "Catalyst",
    "Solvent",
    "Ligand",
    "Additive",
    "Coupling Reagent",
    "Secondary Solvent",
    "Tertiary Solvent",
    "conditions",
)


@dataclass(frozen=True)
class LabelParticipant:
    source_label: str
    base_label: str
    display_label: str
    signature: str
    center_class: Optional[str]
    attachment_class: Optional[str]
    alpha_branched: Optional[bool]
    qualifier_scope: Optional[str]
    mapping_status: str


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
        source_label=str(row.get(f"{prefix} Source") or "").strip(),
        base_label=str(row.get(prefix) or "").strip(),
        display_label=str(row.get(f"{prefix} Display") or "").strip(),
        signature=str(row.get(f"{prefix} Signature") or "").strip(),
        center_class=str(row.get(f"{prefix} Center Class") or "").strip() or None,
        attachment_class=str(row.get(f"{prefix} Attachment Class") or "").strip() or None,
        alpha_branched=_optional_bool(str(row.get(f"{prefix} Alpha Branched") or "")),
        qualifier_scope=str(row.get(f"{prefix} Qualifier Scope") or "").strip() or None,
        mapping_status=str(row.get(f"{prefix} Mapping Status") or "unresolved").strip(),
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
            "yield%",
            "Reaction Type",
            "FG A Signature",
            "FG B Signature",
            "FG A Mapping Status",
            "FG B Mapping Status",
        }
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"Missing cleaned label columns: {missing}")
        rows = list(reader)

    indexed = []
    for row_number, row in enumerate(rows, start=2):
        try:
            yield_pct = float(row.get("yield%") or "")
        except (TypeError, ValueError):
            continue
        if not 0.0 <= yield_pct <= 100.0:
            continue
        try:
            z_score = float(row.get("z-Score") or 0.0)
        except (TypeError, ValueError):
            z_score = 0.0
        conditions = {name: str(row.get(name) or "").strip() for name in CONDITION_COLUMNS}
        indexed.append(
            LabelIndexedReaction(
                source_row_number=row_number,
                reaction_type=str(row.get("Reaction Type") or "").strip(),
                participants=(_participant(row, "FG A"), _participant(row, "FG B")),
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
