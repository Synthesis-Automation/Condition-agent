"""Read chemistry-free intermediate observations for generic conversion."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, Iterator, Mapping, Optional, Tuple

from condition_registry import ConditionComponentInput, ConditionProcessStage

from ..condition_normalization import normalize_cas_list
from ..ingestion.models import INTERMEDIATE_OBSERVATION_SCHEMA_VERSION
from .input_schema import RawReactionRecord


_IDENTIFIER_PRIORITY = {
    "substance_id": 0,
    "cas": 1,
    "name": 2,
    "mfcd": 3,
    "smiles": 4,
}


def _float(value: Any) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _identifier_input(
    component: Mapping[str, Any],
    groups: Mapping[str, Mapping[str, Any]],
) -> Optional[ConditionComponentInput]:
    identifiers = tuple(
        item
        for item in component.get("identifiers") or ()
        if isinstance(item, Mapping) and str(item.get("value") or "").strip()
    )
    if not identifiers:
        return None
    selected = sorted(
        identifiers,
        key=lambda item: (
            _IDENTIFIER_PRIORITY.get(str(item.get("identifier_type") or ""), 99),
            str(item.get("value") or ""),
        ),
    )[0]
    source_type = str(selected.get("identifier_type") or "auto")
    registry_type = (
        source_type if source_type in {"cas", "name", "substance_id"} else "name"
    )
    group_key = str(component.get("group_key") or "")
    amount = component.get("amount")
    amount_mapping = amount if isinstance(amount, Mapping) else {}
    group = groups.get(group_key) or {}
    return ConditionComponentInput(
        raw_identifier=str(selected.get("value") or "").strip(),
        source_field=str(component.get("source_slot") or "condition_component"),
        identifier_type=registry_type,
        source_role_hint=(
            str(component.get("source_role_hint"))
            if component.get("source_role_hint")
            else None
        ),
        amount=_float(amount_mapping.get("value")),
        amount_unit=str(amount_mapping.get("unit") or "") or None,
        provenance={
            "component_key": str(component.get("component_key") or ""),
            "group_key": group_key or None,
            "source_group": dict(group),
            "source_identifier_type": source_type,
            "source_identifiers": tuple(dict(item) for item in identifiers),
            **dict(component.get("provenance") or {}),
        },
    )


def _condition_inputs(
    conditions: Mapping[str, Any],
) -> Tuple[ConditionComponentInput, ...]:
    groups = {
        str(group.get("group_key") or ""): group
        for group in conditions.get("component_groups") or ()
        if isinstance(group, Mapping)
    }
    values = (
        _identifier_input(component, groups)
        for component in conditions.get("components") or ()
        if isinstance(component, Mapping)
    )
    return tuple(value for value in values if value is not None)


def _process_stages(
    conditions: Mapping[str, Any],
) -> Tuple[ConditionProcessStage, ...]:
    stages = []
    for value in conditions.get("stages") or ():
        if not isinstance(value, Mapping):
            continue
        stages.append(
            ConditionProcessStage(
                stage_index=int(value.get("stage_index") or len(stages) + 1),
                temperature_c=_float(value.get("temperature_c")),
                time_h=_float(value.get("time_h")),
                atmosphere=(
                    str(value.get("atmosphere")) if value.get("atmosphere") else None
                ),
                provenance={
                    "source_text": str(value.get("source_text") or ""),
                    **dict(value.get("provenance") or {}),
                },
            )
        )
    return tuple(stages)


def _primary_yield(
    outcomes: Tuple[Mapping[str, Any], ...],
) -> tuple[Optional[float], str]:
    for outcome in outcomes:
        outcome_type = str(outcome.get("outcome_type") or "")
        if outcome_type.endswith("yield_pct"):
            return _float(outcome.get("value")), outcome_type
    return None, ""


def raw_record_from_intermediate(payload: Mapping[str, Any]) -> RawReactionRecord:
    """Project one intermediate observation into the generic converter input."""
    if (
        str(payload.get("schema_version") or "")
        != INTERMEDIATE_OBSERVATION_SCHEMA_VERSION
    ):
        raise ValueError("Unsupported intermediate observation schema")
    source = payload.get("source") or {}
    reaction = payload.get("reaction") or {}
    conditions = payload.get("conditions") or {}
    if not all(isinstance(value, Mapping) for value in (source, reaction, conditions)):
        raise ValueError("Intermediate observation contains invalid nested objects")
    outcome_values = tuple(
        value for value in payload.get("outcomes") or () if isinstance(value, Mapping)
    )
    yield_pct, primary_outcome_type = _primary_yield(outcome_values)
    inputs = _condition_inputs(conditions)
    stages = _process_stages(conditions)
    single_stage = stages[0] if len(stages) == 1 else None
    catalyst_cas = []
    reagent_cas = []
    solvent_cas = []
    for component in conditions.get("components") or ():
        if not isinstance(component, Mapping):
            continue
        role = str(component.get("source_role_hint") or "")
        for identifier in component.get("identifiers") or ():
            if not isinstance(identifier, Mapping):
                continue
            if str(identifier.get("identifier_type") or "") != "cas":
                continue
            value = str(identifier.get("value") or "")
            if role == "catalyst":
                catalyst_cas.append(value)
            elif role == "solvent":
                solvent_cas.append(value)
            else:
                reagent_cas.append(value)
    raw_source_fields = payload.get("raw_fields") or {}
    return RawReactionRecord(
        source_dataset=str(source.get("corpus_id") or ""),
        source_path=str(source.get("source_file") or ""),
        source_row_number=int(source.get("source_row_number") or 0),
        reaction_id=str(
            source.get("source_record_id") or payload.get("observation_id") or ""
        ),
        source_declared_family=str(reaction.get("source_reaction_type") or ""),
        reaction_smiles=str(reaction.get("reaction_smiles") or ""),
        yield_pct=yield_pct,
        temperature_c=single_stage.temperature_c if single_stage else None,
        time_h=single_stage.time_h if single_stage else None,
        reference=str(source.get("reference") or ""),
        reactant_cas=(),
        product_cas=(),
        reagent_cas=normalize_cas_list(",".join(reagent_cas)),
        catalyst_cas=normalize_cas_list(",".join(catalyst_cas)),
        solvent_cas=normalize_cas_list(",".join(solvent_cas)),
        experimental_procedure=str(conditions.get("procedure_text") or ""),
        stages=str(conditions.get("declared_stage_count") or len(stages) or ""),
        steps="",
        notes="",
        condition_component_inputs=inputs,
        condition_process_stages=stages,
        condition_declared_absences=tuple(
            str(value) for value in conditions.get("declared_absences") or ()
        ),
        primary_outcome_type=primary_outcome_type,
        upstream_observation_id=str(payload.get("observation_id") or ""),
        raw_fields={
            "source_observation": dict(payload),
            "source_raw_fields": dict(raw_source_fields)
            if isinstance(raw_source_fields, Mapping)
            else {},
        },
        schema_version=INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
    )


def iter_intermediate_records(path: str | Path) -> Iterator[RawReactionRecord]:
    """Stream generic converter inputs from a compressed intermediate file."""
    source = Path(path)
    opener = gzip.open if source.suffix.casefold() == ".gz" else open
    with opener(source, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            try:
                payload = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid intermediate JSON at {source}:{line_number}"
                ) from exc
            if not isinstance(payload, Mapping):
                raise ValueError(
                    f"Intermediate row is not an object at {source}:{line_number}"
                )
            yield raw_record_from_intermediate(payload)


__all__ = ["iter_intermediate_records", "raw_record_from_intermediate"]
