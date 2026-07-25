"""Convert flat weak-label condition fields into registry-owned recipes."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Tuple

from condition_registry import (
    ConditionComponentInput,
    ConditionProcessStage,
    ResolvedConditionRecipe,
    build_resolved_recipe_from_inputs,
    resolve_substance,
)


SOURCE_CONDITION_FIELDS: Tuple[Tuple[str, str], ...] = (
    ("catalyst", "catalyst"),
    ("ligand", "ligand"),
    ("base", "base"),
    ("primary_solvent", "solvent"),
    ("secondary_solvent", "solvent"),
    ("tertiary_solvent", "solvent"),
    ("additive", "additive"),
    ("coupling_reagent", "coupling_reagent"),
)

ABSENCE_MARKERS = {
    "no catalyst": "catalyst",
    "no ligand": "ligand",
    "no base": "base",
    "no additive": "additive",
    "neat": "solvent",
    "no solvent": "solvent",
    "solvent free": "solvent",
    "solvent-free": "solvent",
}


def condition_recipe_catalog_path(csv_path: str | Path) -> Path:
    """Return the unique canonical recipe catalog paired with a label CSV."""
    path = Path(csv_path)
    return path.with_name(f"{path.stem}.condition_recipes.jsonl.gz")

_STAGE_PATTERNS = (
    re.compile(
        r"(?P<time>\d+(?:\.\d+)?)\s*(?P<unit>h|hr|hrs|hour|hours|min|mins|minute|minutes)"
        r"\s+at\s+(?P<temp>-?\d+(?:\.\d+)?)\s*°?\s*C\b",
        re.IGNORECASE,
    ),
    re.compile(
        r"(?P<temp>-?\d+(?:\.\d+)?)\s*°?\s*C\s+for\s+"
        r"(?P<time>\d+(?:\.\d+)?)\s*(?P<unit>h|hr|hrs|hour|hours|min|mins|minute|minutes)\b",
        re.IGNORECASE,
    ),
)


@dataclass(frozen=True)
class ConvertedLabelConditions:
    """Canonical condition recipe plus concise CSV review fields."""

    recipe: ResolvedConditionRecipe
    condition_display: str
    procedure_text: str

    def to_columns(self) -> Dict[str, str]:
        return {
            "condition_recipe_id": self.recipe.recipe_id,
            "condition_display": self.condition_display,
            "temperature_c": _optional_number(self.recipe.temperature_c),
            "time_h": _optional_number(self.recipe.time_h),
            "concentration_m": _optional_number(self.recipe.concentration_m),
            "atmosphere": self.recipe.atmosphere or "",
            "condition_identity_uncertainty": (
                "true"
                if "CONDITION_IDENTITY_UNCERTAINTY" in self.recipe.warnings
                else "false"
            ),
            "procedure_text": self.procedure_text,
        }


def _optional_number(value: Optional[float]) -> str:
    if value is None:
        return ""
    return f"{value:g}"


def _absence_role(value: str) -> Optional[str]:
    return ABSENCE_MARKERS.get(value.strip().lower())


def _all_parts_resolve(parts: Iterable[str]) -> bool:
    values = tuple(part.strip() for part in parts if part.strip())
    return bool(values) and all(
        _absence_role(value) is not None
        or resolve_substance(name=value).status == "resolved"
        for value in values
    )


def _split_component_cell(value: str) -> Tuple[str, ...]:
    """Split only when delimiters provide safe, deterministic evidence."""
    raw = value.strip()
    if not raw:
        return ()
    if _absence_role(raw) is not None:
        return (raw,)
    if resolve_substance(name=raw).status == "resolved":
        return (raw,)
    if ";" in raw or "|" in raw:
        return tuple(
            part.strip() for part in re.split(r"[;|]", raw) if part.strip()
        )
    if "," in raw:
        parts = tuple(part.strip() for part in raw.split(",") if part.strip())
        if len(parts) > 1 and _all_parts_resolve(parts):
            return parts
    return (raw,)


def parse_process_stages(procedure_text: str) -> Tuple[ConditionProcessStage, ...]:
    """Extract explicit temperature/time stages without interpreting free text."""
    matches = []
    for pattern in _STAGE_PATTERNS:
        matches.extend(pattern.finditer(procedure_text))
    matches.sort(key=lambda item: item.start())
    stages = []
    occupied = set()
    for match in matches:
        span = match.span()
        if any(start < span[1] and span[0] < end for start, end in occupied):
            continue
        occupied.add(span)
        time_value = float(match.group("time"))
        if match.group("unit").lower().startswith("min"):
            time_value /= 60.0
        stages.append(
            ConditionProcessStage(
                stage_index=len(stages) + 1,
                temperature_c=float(match.group("temp")),
                time_h=time_value,
                provenance={"source": "procedure_text", "text": match.group(0)},
            )
        )
    return tuple(stages)


def convert_label_conditions(
    row: Mapping[str, str],
) -> ConvertedLabelConditions:
    """Normalize one flat weak-label condition row into a canonical recipe."""
    inputs = []
    absences = set()
    display_parts = []
    for source_position, (source_field, role_hint) in enumerate(
        SOURCE_CONDITION_FIELDS, start=1
    ):
        raw_cell = str(row.get(source_field) or "").strip()
        if not raw_cell:
            continue
        for raw_identifier in _split_component_cell(raw_cell):
            absence = _absence_role(raw_identifier)
            if absence is not None:
                absences.add(absence)
                display_parts.append(f"No {absence}")
                continue
            inputs.append(
                ConditionComponentInput(
                    raw_identifier=raw_identifier,
                    source_field=source_field,
                    identifier_type="auto",
                    source_role_hint=role_hint,
                    provenance={"source_position": source_position},
                )
            )
            display_parts.append(f"{raw_identifier} [{role_hint}]")

    procedure_text = str(row.get("procedure_text") or "").strip()
    stages = parse_process_stages(procedure_text)
    temperature_c = stages[0].temperature_c if len(stages) == 1 else None
    time_h = stages[0].time_h if len(stages) == 1 else None
    recipe = build_resolved_recipe_from_inputs(
        inputs,
        temperature_c=temperature_c,
        time_h=time_h,
        stages=stages,
        declared_absences=absences,
    )
    return ConvertedLabelConditions(
        recipe=recipe,
        condition_display="; ".join(display_parts),
        procedure_text=procedure_text,
    )


__all__ = [
    "ABSENCE_MARKERS",
    "ConvertedLabelConditions",
    "SOURCE_CONDITION_FIELDS",
    "convert_label_conditions",
    "condition_recipe_catalog_path",
    "parse_process_stages",
]
