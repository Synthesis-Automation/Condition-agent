"""Generic conversion of one raw observation into a recommendation record."""

from __future__ import annotations

from dataclasses import asdict
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Hashable, TypeVar

from condition_registry import (
    ResolvedConditionRecipe,
    build_resolved_recipe,
    build_resolved_recipe_from_inputs,
)

from ..condition_roles import preferred_roles_for_reaction
from reactive_taxonomy import (
    AtomMappingProvider,
    ExternalMappingAssessment,
    analyze_reaction_with_external_mapping,
    featurize_reaction,
)

from ..condition_normalization import normalize_cas_list
from ..core_eligibility import assess_core_eligibility
from ..fragment_source_support import assess_fragment_source_support
from ..models import ConditionIdentity, PrecedentTier, RecommendationRecord
from .admission import decide_admission
from .identities import canonical_reaction_identity, observation_id, raw_recipe_id
from .input_schema import RawReactionRecord
from .references import normalize_reference
from .reference_series import reference_condition_series_id
from .signature_serialization import signature_record_fields

_T = TypeVar("_T")


@dataclass
class GenericConversionCache:
    """Bounded, invocation-owned cache for deterministic conversion work."""

    max_entries: int = 100_000
    analyses: Dict[Hashable, Any] = field(default_factory=dict)
    external_assessments: Dict[Hashable, Any] = field(default_factory=dict)
    canonical_reactions: Dict[Hashable, Any] = field(default_factory=dict)
    references: Dict[Hashable, Any] = field(default_factory=dict)
    recipes: Dict[Hashable, Any] = field(default_factory=dict)

    def get(
        self,
        values: Dict[Hashable, _T],
        key: Hashable,
        factory: Callable[[], _T],
    ) -> _T:
        """Return a cached deterministic result without unbounded growth."""
        if key in values:
            return values[key]
        value = factory()
        if len(values) >= self.max_entries:
            values.clear()
        values[key] = value
        return value


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
                "role_status": component.role_status,
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
        or any(
            component["status"] != "resolved"
            or component["role_status"] != "assigned"
            for component in components
        ),
        "recipe_warnings": recipe.warnings,
        "schema_version": "2.0",
    }


def convert_record(
    record: RawReactionRecord,
    *,
    cache: GenericConversionCache | None = None,
    mapping_provider: AtomMappingProvider | None = None,
) -> RecommendationRecord:
    """Convert one source-faithful row without a declared-family requirement."""
    assessment: ExternalMappingAssessment | None = None
    if mapping_provider is not None:
        provider_metadata = mapping_provider.metadata
        base_analysis = (
            featurize_reaction(record.reaction_smiles)
            if cache is None
            else cache.get(
                cache.analyses,
                record.reaction_smiles,
                lambda: featurize_reaction(record.reaction_smiles),
            )
        )
        assessment_key = (
            record.reaction_smiles,
            provider_metadata.provider_id,
            provider_metadata.provider_version,
            provider_metadata.model_id,
            provider_metadata.model_sha256,
            "shadow_missing_reaction_core",
        )
        force_resolved_shadow = bool(
            base_analysis.reaction_signature is not None
            and base_analysis.reaction_core is None
        )
        assessment = (
            analyze_reaction_with_external_mapping(
                record.reaction_smiles,
                mapping_provider,
                base_analysis=base_analysis,
                force_resolved_shadow=force_resolved_shadow,
            )
            if cache is None
            else cache.get(
                cache.external_assessments,
                assessment_key,
                lambda: analyze_reaction_with_external_mapping(
                    record.reaction_smiles,
                    mapping_provider,
                    base_analysis=base_analysis,
                    force_resolved_shadow=force_resolved_shadow,
                ),
            )
        )
        # External mapping is supplemental evidence. Canonical signatures,
        # fallback descriptors, and admission remain grounded in the internal
        # observation; a mapper-derived core is attached separately below.
        analysis = base_analysis
    elif cache is None:
        analysis = featurize_reaction(record.reaction_smiles)
    else:
        analysis = cache.get(
            cache.analyses,
            record.reaction_smiles,
            lambda: featurize_reaction(record.reaction_smiles),
        )
    if cache is None:
        canonical_identity = canonical_reaction_identity(record.reaction_smiles)
        reference_identity = normalize_reference(record.reference)
    else:
        canonical_identity = cache.get(
            cache.canonical_reactions,
            record.reaction_smiles,
            lambda: canonical_reaction_identity(record.reaction_smiles),
        )
        reference_identity = cache.get(
            cache.references,
            record.reference,
            lambda: normalize_reference(record.reference),
        )
    conditions = _normalized_conditions(record)
    recipe_key = (
        record.catalyst_cas,
        record.reagent_cas,
        record.solvent_cas,
        tuple(
            (
                item.raw_identifier,
                item.source_field,
                item.identifier_type,
                item.source_role_hint,
                item.amount,
                item.amount_unit,
                repr(sorted(item.provenance.items())),
            )
            for item in record.condition_component_inputs
        ),
        analysis.transformation_class,
        analysis.named_family,
        record.temperature_c,
        record.time_h,
        tuple(
            (
                item.stage_index,
                item.temperature_c,
                item.time_h,
                item.atmosphere,
                repr(sorted(item.provenance.items())),
            )
            for item in record.condition_process_stages
        ),
        record.condition_declared_absences,
    )

    def build_recipe() -> ResolvedConditionRecipe:
        preferred_roles = preferred_roles_for_reaction(
            analysis.transformation_class
        )
        if record.condition_component_inputs:
            return build_resolved_recipe_from_inputs(
                record.condition_component_inputs,
                preferred_roles=preferred_roles,
                temperature_c=record.temperature_c,
                time_h=record.time_h,
                stages=record.condition_process_stages,
                declared_absences=record.condition_declared_absences,
            )
        return build_resolved_recipe(
            {
                "catalyst_cas": record.catalyst_cas,
                "reagent_cas": record.reagent_cas,
                "solvent_cas": record.solvent_cas,
            },
            preferred_roles=preferred_roles,
            temperature_c=record.temperature_c,
            time_h=record.time_h,
        )

    resolved_recipe = (
        build_recipe()
        if cache is None
        else cache.get(cache.recipes, recipe_key, build_recipe)
    )
    source_requirements = (
        analysis.fallback_descriptor.source_requirements
        if analysis.fallback_descriptor is not None
        else ()
    )
    fragment_source_support = assess_fragment_source_support(
        source_requirements,
        resolved_recipe,
    )
    decision = decide_admission(
        record=record,
        analysis=analysis,
        canonical_identity=canonical_identity,
        conditions=conditions,
        resolved_recipe=resolved_recipe,
        fragment_source_support=fragment_source_support,
    )
    signature_fields = signature_record_fields(analysis)
    if (
        assessment is not None
        and analysis.reaction_core is None
        and assessment.analysis.reaction_core is not None
    ):
        signature_fields["reaction_core"] = asdict(
            assessment.analysis.reaction_core
        )
    completeness_payload = (
        asdict(analysis.reaction_completeness)
        if analysis.reaction_completeness
        else None
    )
    partial_payload = (
        asdict(analysis.partial_product_transformation)
        if analysis.partial_product_transformation
        else None
    )
    external_mapping_payload = (
        assessment.to_provenance_dict() if assessment is not None else None
    )
    core_eligibility = assess_core_eligibility(
        {
            "reaction_core": signature_fields.get("reaction_core"),
            "warnings": analysis.warnings,
            "index_eligibility": decision.index_eligibility,
            "chemistry_status": decision.chemistry_status,
            "condition_status": decision.condition_status,
            "condition_stage_status": decision.condition_stage_status,
            "reaction_completeness": completeness_payload,
            "partial_product_transformation": partial_payload,
            "external_atom_mapping": external_mapping_payload,
        }
    )
    condition_series_id = reference_condition_series_id(
        reference_id=reference_identity.reference_id,
        recipe_core_id=resolved_recipe.recipe_core_id,
        chemistry_key=(
            analysis.reaction_signature.transformation_signature_key
            if analysis.reaction_signature
            else ""
        ),
        stages=record.stages,
        steps=record.steps,
        temperature_c=record.temperature_c,
        time_h=record.time_h,
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
    if (
        record.primary_outcome_type
        or record.condition_component_inputs
        or record.condition_process_stages
        or record.condition_declared_absences
    ):
        source.update(
            {
                "primary_outcome_type": record.primary_outcome_type,
                "condition_component_inputs": tuple(
                    asdict(item) for item in record.condition_component_inputs
                ),
                "condition_process_stages": tuple(
                    asdict(item) for item in record.condition_process_stages
                ),
                "condition_declared_absences": (record.condition_declared_absences),
            }
        )
    return RecommendationRecord(
        reaction_id=record.reaction_id,
        source_row_number=record.source_row_number,
        reaction_smiles=record.reaction_smiles,
        admission_tier=decision.tier,
        admission_reasons=decision.reasons,
        evidence_quality=analysis.evidence_quality,
        named_family=analysis.named_family,
        reaction_label=(
            asdict(analysis.reaction_label)
            if analysis.reaction_label is not None
            else None
        ),
        yield_pct=record.yield_pct,
        temperature_c=record.temperature_c,
        time_h=record.time_h,
        conditions=conditions,
        chemistry_status=decision.chemistry_status,
        condition_status=decision.condition_status,
        condition_stage_status=decision.condition_stage_status,
        outcome_status=decision.outcome_status,
        index_eligibility=decision.index_eligibility,
        precedent_tier=(
            PrecedentTier.TRUSTED
            if decision.index_eligibility.value == "eligible"
            else (
                PrecedentTier.REVIEW_CORE
                if decision.index_eligibility.value == "review_only"
                and core_eligibility.tier.value == "review_core"
                else None
            )
        ),
        core_eligibility=core_eligibility.tier,
        core_eligibility_reasons=core_eligibility.reasons,
        core_eligibility_definition_version=(
            core_eligibility.definition_version
        ),
        family_environment=asdict(analysis.family_environment)
        if analysis.family_environment
        else None,
        product_connection=asdict(analysis.product_connection)
        if analysis.product_connection
        else None,
        spectator_groups=tuple(asdict(group) for group in analysis.spectator_groups),
        partial_product_transformation=partial_payload,
        reaction_completeness=completeness_payload,
        reaction_observation=(
            asdict(analysis.observation) if analysis.observation is not None else None
        ),
        reaction_interpretation=(
            asdict(analysis.interpretation)
            if analysis.interpretation is not None
            else None
        ),
        fragment_source_support=tuple(
            support.to_dict() for support in fragment_source_support
        ),
        **signature_fields,
        external_atom_mapping=external_mapping_payload,
        source_dataset=record.source_dataset,
        source_path=record.source_path,
        source_declared_family=record.source_declared_family,
        reference_id=reference_identity.reference_id,
        reference_identity=reference_identity.to_dict(),
        observation_id=record.upstream_observation_id or observation_id(record),
        canonical_reaction_id=canonical_identity.reaction_id
        if canonical_identity
        else None,
        canonical_reaction_smiles=canonical_identity.canonical_reaction_smiles
        if canonical_identity
        else None,
        raw_recipe_id=raw_recipe_id(record),
        resolved_recipe_core_id=resolved_recipe.recipe_core_id,
        condition_resolution=_condition_resolution(resolved_recipe),
        resolved_recipe_id=resolved_recipe.recipe_id,
        resolved_recipe=resolved_recipe.to_dict(),
        reference_condition_series_id=condition_series_id,
        source=source,
    )


__all__ = ["GenericConversionCache", "convert_record"]
