"""Generic conversion of one raw observation into a recommendation record."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any, Dict

from condition_registry import ResolvedConditionRecipe, build_resolved_recipe
from reactive_taxonomy import featurize_reaction

from ..condition_normalization import normalize_cas_list
from ..models import ConditionIdentity, RecommendationRecord
from .admission import decide_admission
from .identities import canonical_reaction_identity, observation_id, raw_recipe_id
from .input_schema import RawReactionRecord
from .signature_serialization import signature_record_fields


def _normalized_conditions(record: RawReactionRecord) -> ConditionIdentity:
    return ConditionIdentity(
        catalyst_cas=normalize_cas_list(",".join(record.catalyst_cas)),
        reagent_cas=normalize_cas_list(",".join(record.reagent_cas)),
        solvent_cas=normalize_cas_list(",".join(record.solvent_cas)),
    )


def _condition_resolution(recipe: ResolvedConditionRecipe) -> Dict[str, Any]:
    components = []
    statuses: Dict[str, int] = {}
    for component in recipe.components:
        statuses[component.identity_status] = (
            statuses.get(component.identity_status, 0) + 1
        )
        components.append(
            {
                "source_field": component.source_field,
                "raw_identifier": component.raw_identifier,
                "status": component.identity_status,
                "substance_id": component.substance_id,
                "canonical_name": component.canonical_name,
                "possible_roles": tuple(role.role_id for role in component.roles),
                "primary_role": component.primary_role,
                "primary_role_confidence": component.primary_role_confidence,
                "warnings": component.warnings,
            }
        )
    return {
        "component_count": len(components),
        "status_counts": dict(sorted(statuses.items())),
        "components": tuple(components),
        "has_uncertainty": bool(recipe.warnings)
        or any(component["status"] != "resolved" for component in components),
        "recipe_warnings": recipe.warnings,
        "schema_version": "1.1",
    }


def convert_record(record: RawReactionRecord) -> RecommendationRecord:
    """Convert one source-faithful row without a declared-family requirement."""
    analysis = featurize_reaction(record.reaction_smiles)
    canonical_identity = canonical_reaction_identity(record.reaction_smiles)
    conditions = _normalized_conditions(record)
    resolved_recipe = build_resolved_recipe(
        {
            "catalyst_cas": record.catalyst_cas,
            "reagent_cas": record.reagent_cas,
            "solvent_cas": record.solvent_cas,
        },
        transformation_class=analysis.transformation_class,
        named_family=analysis.named_family,
        temperature_c=record.temperature_c,
        time_h=record.time_h,
    )
    decision = decide_admission(
        record=record,
        analysis=analysis,
        canonical_identity=canonical_identity,
        conditions=conditions,
    )
    source = {
        "source_dataset": record.source_dataset,
        "source_path": record.source_path,
        "source_row_number": record.source_row_number,
        "source_declared_family": record.source_declared_family,
        "reference": record.reference,
        "reactant_cas": record.reactant_cas,
        "product_cas": record.product_cas,
        "raw_condition_identifiers": {
            "catalyst_cas": record.catalyst_cas,
            "reagent_cas": record.reagent_cas,
            "solvent_cas": record.solvent_cas,
        },
        "experimental_procedure": record.experimental_procedure,
        "stages": record.stages,
        "steps": record.steps,
        "notes": record.notes,
        "raw_fields": record.raw_fields,
        "input_schema_version": record.schema_version,
        "admission_policy_version": decision.policy_version,
    }
    return RecommendationRecord(
        reaction_id=record.reaction_id,
        source_row_number=record.source_row_number,
        reaction_smiles=record.reaction_smiles,
        admission_tier=decision.tier,
        admission_reasons=decision.reasons,
        evidence_quality=analysis.evidence_quality,
        named_family=analysis.named_family,
        reaction_label=analysis.reaction_label,
        reaction_label_status=analysis.reaction_label_status,
        yield_pct=record.yield_pct,
        temperature_c=record.temperature_c,
        time_h=record.time_h,
        conditions=conditions,
        family_environment=asdict(analysis.family_environment)
        if analysis.family_environment
        else None,
        product_connection=asdict(analysis.product_connection)
        if analysis.product_connection
        else None,
        spectator_groups=tuple(asdict(group) for group in analysis.spectator_groups),
        **signature_record_fields(analysis),
        source_dataset=record.source_dataset,
        source_path=record.source_path,
        source_declared_family=record.source_declared_family,
        observation_id=observation_id(record),
        canonical_reaction_id=canonical_identity.reaction_id
        if canonical_identity
        else None,
        canonical_reaction_smiles=canonical_identity.canonical_reaction_smiles
        if canonical_identity
        else None,
        raw_recipe_id=raw_recipe_id(record),
        condition_resolution=_condition_resolution(resolved_recipe),
        resolved_recipe_id=resolved_recipe.recipe_id,
        resolved_recipe=resolved_recipe.to_dict(),
        source=source,
    )


__all__ = ["convert_record"]
