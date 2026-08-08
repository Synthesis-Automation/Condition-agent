"""Canonical machine-readable protocol drafts built from resolved recipes."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from typing import Any, Dict, Iterable, Literal, Mapping, Optional, Tuple

from .api import resolve_substance_id
from .models import CONDITION_RECIPE_COMPONENT_BUCKETS, ResolvedConditionRecipe


CONDITION_PROTOCOL_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class ProtocolReactionMaterialInput:
    """Narrow reaction-material observation supplied by the recommender."""

    side: Literal["reactant", "agent", "product"]
    component_index: int
    smiles: str
    canonical_smiles: Optional[str] = None
    cas: Optional[str] = None
    amount: Optional[float] = None
    amount_unit: Optional[str] = None


@dataclass(frozen=True)
class ProtocolMaterial:
    """One identified or unresolved material in a protocol draft."""

    material_id: str
    category: Literal["reaction_input", "condition", "reaction_output"]
    identity_status: str
    substance_id: Optional[str]
    canonical_name: Optional[str]
    cas: Optional[str]
    smiles: Optional[str]
    role: Optional[str]
    role_status: str
    amount: Optional[float]
    amount_unit: Optional[str]
    source_field: Optional[str] = None
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class ProtocolOperatingConditions:
    """Unit-explicit physical setpoints observed for a reaction."""

    temperature_c: Optional[float] = None
    time_h: Optional[float] = None
    concentration_m: Optional[float] = None
    atmosphere: Optional[str] = None


@dataclass(frozen=True)
class ProtocolOperation:
    """One observed process stage; it is not an inferred addition step."""

    operation_id: str
    operation_type: Literal["maintain_conditions"]
    stage_index: int
    temperature_c: Optional[float] = None
    time_h: Optional[float] = None
    atmosphere: Optional[str] = None


@dataclass(frozen=True)
class SynthesisProtocolDraft:
    """Review-gated protocol view that never invents missing operations."""

    protocol_id: str
    recipe_id: str
    recipe_core_id: str
    reaction_smiles: Optional[str]
    materials: Tuple[ProtocolMaterial, ...]
    operating_conditions: ProtocolOperatingConditions
    operations: Tuple[ProtocolOperation, ...]
    execution_readiness: Literal["review_required"]
    missing_required_fields: Tuple[str, ...]
    warnings: Tuple[str, ...] = ()
    schema_version: str = CONDITION_PROTOCOL_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the protocol draft as canonical nested JSON data."""
        return asdict(self)


def _value(source: Any, key: str, default: Any = None) -> Any:
    if isinstance(source, Mapping):
        return source.get(key, default)
    return getattr(source, key, default)


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def _condition_materials(recipe: Any) -> Tuple[ProtocolMaterial, ...]:
    materials = []
    index = 0
    for bucket in CONDITION_RECIPE_COMPONENT_BUCKETS:
        for component in _value(recipe, bucket, ()) or ():
            index += 1
            substance_id = _value(component, "substance_id")
            canonical_name = _value(component, "canonical_name")
            cas = _value(component, "cas")
            if substance_id and (not canonical_name or not cas):
                resolution = resolve_substance_id(str(substance_id))
                if resolution.status == "resolved" and resolution.substance is not None:
                    canonical_name = canonical_name or resolution.substance.canonical_name
                    cas = cas or resolution.substance.cas
            materials.append(
                ProtocolMaterial(
                    material_id=f"condition_{index:03d}",
                    category="condition",
                    identity_status=str(_value(component, "identity_status", "unresolved")),
                    substance_id=substance_id,
                    canonical_name=canonical_name,
                    cas=cas,
                    smiles=None,
                    role=_value(component, "primary_role"),
                    role_status=str(_value(component, "role_status", "unassigned")),
                    amount=_value(component, "amount"),
                    amount_unit=_value(component, "amount_unit"),
                    source_field=_value(component, "source_field"),
                    warnings=tuple(_value(component, "warnings", ()) or ()),
                )
            )
    return tuple(materials)


def _reaction_materials(
    values: Iterable[ProtocolReactionMaterialInput],
) -> Tuple[ProtocolMaterial, ...]:
    materials = []
    side_order = {"reactant": 0, "agent": 1, "product": 2}
    for item in sorted(
        values,
        key=lambda value: (side_order[value.side], value.component_index),
    ):
        category = (
            "reaction_output" if item.side == "product" else "reaction_input"
        )
        materials.append(
            ProtocolMaterial(
                material_id=f"{item.side}_{item.component_index + 1:03d}",
                category=category,
                identity_status="structure_only",
                substance_id=None,
                canonical_name=None,
                cas=item.cas,
                smiles=item.canonical_smiles or item.smiles,
                role=item.side,
                role_status="assigned",
                amount=item.amount,
                amount_unit=item.amount_unit,
            )
        )
    return tuple(materials)


def _operations(recipe: Any) -> Tuple[ProtocolOperation, ...]:
    stages = tuple(_value(recipe, "stages", ()) or ())
    if stages:
        return tuple(
            ProtocolOperation(
                operation_id=f"maintain_{int(_value(stage, 'stage_index', index)):03d}",
                operation_type="maintain_conditions",
                stage_index=int(_value(stage, "stage_index", index)),
                temperature_c=_value(stage, "temperature_c"),
                time_h=_value(stage, "time_h"),
                atmosphere=_value(stage, "atmosphere"),
            )
            for index, stage in enumerate(stages, start=1)
        )
    if any(
        _value(recipe, key) is not None
        for key in ("temperature_c", "time_h", "atmosphere")
    ):
        return (
            ProtocolOperation(
                operation_id="maintain_001",
                operation_type="maintain_conditions",
                stage_index=1,
                temperature_c=_value(recipe, "temperature_c"),
                time_h=_value(recipe, "time_h"),
                atmosphere=_value(recipe, "atmosphere"),
            ),
        )
    return ()


def build_synthesis_protocol_draft(
    recipe: ResolvedConditionRecipe | Mapping[str, Any],
    *,
    reaction_smiles: Optional[str] = None,
    reaction_materials: Iterable[ProtocolReactionMaterialInput] = (),
) -> SynthesisProtocolDraft:
    """Build a deterministic protocol draft while preserving missing data."""
    condition_materials = _condition_materials(recipe)
    observed_reaction_materials = _reaction_materials(reaction_materials)
    materials = (*observed_reaction_materials, *condition_materials)
    operations = _operations(recipe)
    conditions = ProtocolOperatingConditions(
        temperature_c=_value(recipe, "temperature_c"),
        time_h=_value(recipe, "time_h"),
        concentration_m=_value(recipe, "concentration_m"),
        atmosphere=_value(recipe, "atmosphere"),
    )
    missing = {"ordered_operations", "vessel", "mixing", "quench", "workup"}
    if not any(item.category == "reaction_input" for item in materials):
        missing.add("reaction_inputs")
    if conditions.temperature_c is None and not any(
        operation.temperature_c is not None for operation in operations
    ):
        missing.add("operating_conditions.temperature_c")
    if conditions.time_h is None and not any(
        operation.time_h is not None for operation in operations
    ):
        missing.add("operating_conditions.time_h")
    for material in materials:
        if material.category != "reaction_output" and material.amount is None:
            missing.add(f"materials.{material.material_id}.amount")
        if material.category == "condition" and material.identity_status != "resolved":
            missing.add(f"materials.{material.material_id}.identity")
    missing_fields = tuple(sorted(missing))
    payload = {
        "recipe_id": str(_value(recipe, "recipe_id", "")),
        "recipe_core_id": str(_value(recipe, "recipe_core_id", "")),
        "reaction_smiles": reaction_smiles,
        "materials": tuple(asdict(item) for item in materials),
        "operating_conditions": asdict(conditions),
        "operations": tuple(asdict(item) for item in operations),
        "missing_required_fields": missing_fields,
        "schema_version": CONDITION_PROTOCOL_SCHEMA_VERSION,
    }
    protocol_id = "CP1:" + hashlib.sha256(
        _canonical_json(payload).encode("utf-8")
    ).hexdigest()
    warnings = tuple(
        sorted(
            {
                "PROTOCOL_REQUIRES_HUMAN_REVIEW",
                *(
                    "CONDITION_IDENTITY_UNCERTAINTY"
                    for item in condition_materials
                    if item.identity_status != "resolved"
                ),
            }
        )
    )
    return SynthesisProtocolDraft(
        protocol_id=protocol_id,
        recipe_id=payload["recipe_id"],
        recipe_core_id=payload["recipe_core_id"],
        reaction_smiles=reaction_smiles,
        materials=materials,
        operating_conditions=conditions,
        operations=operations,
        execution_readiness="review_required",
        missing_required_fields=missing_fields,
        warnings=warnings,
    )


__all__ = [
    "CONDITION_PROTOCOL_SCHEMA_VERSION",
    "ProtocolMaterial",
    "ProtocolOperatingConditions",
    "ProtocolOperation",
    "ProtocolReactionMaterialInput",
    "SynthesisProtocolDraft",
    "build_synthesis_protocol_draft",
]
