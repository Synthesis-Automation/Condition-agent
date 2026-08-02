"""Adapter for the original v2.1 weak-label condition dataset."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterator, List, Tuple

from ..models import (
    CanonicalSourceObservation,
    ConditionComponentClaim,
    ConditionInput,
    ConditionStageInput,
    OutcomeInput,
    ReactionEvidenceInput,
    SourceIdentifier,
)
from .base import (
    clean_text,
    observation_id,
    optional_float,
    raw_fields,
    source_provenance,
    validate_headers,
)


_STAGE_PATTERN = re.compile(
    r"(?P<time>\d+(?:\.\d+)?)\s*"
    r"(?P<time_unit>h|hr|hrs|hour|hours|min|mins|minute|minutes|day|days)"
    r"\s*(?:at|@)\s*"
    r"(?P<temperature>-?\d+(?:\.\d+)?)\s*°?\s*C",
    re.IGNORECASE,
)
_ABSENCE_PATTERN = re.compile(r"^no\s+(.+)$", re.IGNORECASE)


def _time_h(value: str, unit: str) -> float:
    number = float(value)
    lowered = unit.casefold()
    if lowered.startswith("min"):
        return number / 60.0
    if lowered.startswith("day"):
        return number * 24.0
    return number


def _stages(
    procedure: str, component_keys: Tuple[str, ...]
) -> Tuple[ConditionStageInput, ...]:
    stages = []
    for stage_index, match in enumerate(_STAGE_PATTERN.finditer(procedure), start=1):
        stages.append(
            ConditionStageInput(
                stage_index=stage_index,
                component_keys=component_keys,
                temperature_c=float(match.group("temperature")),
                time_h=_time_h(match.group("time"), match.group("time_unit")),
                source_text=match.group(0),
                provenance={"source": "conditions"},
            )
        )
    if stages:
        return tuple(stages)
    return (
        ConditionStageInput(
            stage_index=1,
            component_keys=component_keys,
            source_text=procedure,
            provenance={"source": "conditions"},
        ),
    )


class WeakLabelCsvAdapter:
    """Retain label-only rows without claiming molecular evidence."""

    adapter_id = "weak_label_v2_1.v1"
    adapter_version = "1.0"
    corpus_id = "weak_label"
    required_columns = (
        "yield%",
        "Base",
        "Catalyst",
        "Solvent",
        "Ligand",
        "Additive",
        "Coupling Reagent",
        "Secondary Solvent",
        "Tertiary Solvent",
        "Reaction Type",
        "FG A",
        "FG B",
        "z-Score",
        "conditions",
    )

    @staticmethod
    def _conditions(row: dict[str, str]) -> ConditionInput:
        components: List[ConditionComponentClaim] = []
        absences = []
        for field_name, role_hint in (
            ("Base", "base"),
            ("Catalyst", "catalyst"),
            ("Solvent", "solvent"),
            ("Ligand", "ligand"),
            ("Additive", "additive"),
            ("Coupling Reagent", "coupling_reagent"),
            ("Secondary Solvent", "solvent"),
            ("Tertiary Solvent", "solvent"),
        ):
            value = clean_text(row.get(field_name))
            if not value:
                continue
            absence = _ABSENCE_PATTERN.match(value)
            if absence:
                absences.append(absence.group(1).strip().casefold())
                continue
            key = re.sub(r"[^a-z0-9]+", "_", field_name.casefold()).strip("_")
            components.append(
                ConditionComponentClaim(
                    component_key=key,
                    source_slot=field_name,
                    source_role_hint=role_hint,
                    identifiers=(SourceIdentifier("name", value, field_name),),
                )
            )
        procedure = clean_text(row.get("conditions"))
        component_keys = tuple(component.component_key for component in components)
        stages = _stages(procedure, component_keys)
        warnings = []
        if procedure and not _STAGE_PATTERN.search(procedure):
            warnings.append("UNPARSED_PROCEDURE_OPERATING_CONDITIONS")
        return ConditionInput(
            components=tuple(components),
            stages=stages,
            declared_stage_count=len(stages),
            procedure_text=procedure,
            declared_absences=tuple(sorted(set(absences))),
            warnings=tuple(warnings),
        )

    def iter_observations(
        self, path: Path, *, source_sha256: str
    ) -> Iterator[CanonicalSourceObservation]:
        """Stream all weak-label rows, including formerly filtered rows."""
        validate_headers(path, self.required_columns)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row_number, row in enumerate(csv.DictReader(handle), start=2):
                warnings: List[str] = ["REACTION_STRUCTURE_NOT_AVAILABLE"]
                reaction_type = clean_text(row.get("Reaction Type"))
                fg_a = clean_text(row.get("FG A"))
                fg_b = clean_text(row.get("FG B"))
                if fg_a == fg_b:
                    warnings.append("WEAK_LABEL_IDENTICAL_SITE_LABELS")
                if not fg_a and not fg_b:
                    warnings.append("WEAK_LABEL_SITES_MISSING")
                if "Protecting Group" in {fg_a, fg_b}:
                    warnings.append("WEAK_LABEL_PROTECTING_GROUP_CLAIM")
                conditions = self._conditions(row)
                warnings.extend(conditions.warnings)
                yield_value, warning = optional_float(row.get("yield%"))
                if warning:
                    warnings.append(f"{warning}:yield%")
                z_score, warning = optional_float(row.get("z-Score"))
                if warning:
                    warnings.append(f"{warning}:z-Score")
                record_id = f"{path.stem}:row-{row_number}"
                outcomes = (
                    OutcomeInput(
                        outcome_type="reported_yield_pct",
                        value=yield_value,
                        unit="percent",
                        raw_value=clean_text(row.get("yield%")),
                        source_field="yield%",
                        metadata={"measurement_basis": "source_unspecified"},
                    ),
                    OutcomeInput(
                        outcome_type="source_z_score",
                        value=z_score,
                        unit="dimensionless",
                        raw_value=clean_text(row.get("z-Score")),
                        source_field="z-Score",
                    ),
                )
                yield CanonicalSourceObservation(
                    observation_id=observation_id(
                        adapter_id=self.adapter_id,
                        source_sha256=source_sha256,
                        row_number=row_number,
                        record_id=record_id,
                    ),
                    observation_kind="label_only",
                    source=source_provenance(
                        adapter=self,
                        path=path,
                        source_sha256=source_sha256,
                        row_number=row_number,
                        record_id=record_id,
                        source_groups={"source_reaction_type": reaction_type},
                    ),
                    reaction=ReactionEvidenceInput(
                        evidence_kind="source_labels_unverified",
                        source_reaction_type=reaction_type,
                        source_labels={
                            "reactive_site_1": fg_a,
                            "reactive_site_2": fg_b,
                        },
                        structure_available=False,
                    ),
                    conditions=conditions,
                    outcomes=outcomes,
                    ingestion_status="accepted",
                    warnings=tuple(sorted(set(warnings))),
                    raw_fields=raw_fields(row),
                )


__all__ = ["WeakLabelCsvAdapter"]
