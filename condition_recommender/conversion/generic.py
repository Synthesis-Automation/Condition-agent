"""Generic conversion of one raw observation into a recommendation record."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict
from functools import lru_cache
from typing import Any, Dict

from condition_registry.resolver import ConditionRegistry
from reactive_taxonomy import featurize_reaction

from ..condition_normalization import normalize_cas_list
from ..models import ConditionIdentity, RecommendationRecord
from .admission import decide_admission
from .identities import canonical_reaction_identity, observation_id, raw_recipe_id
from .input_schema import RawReactionRecord
from .signature_serialization import signature_record_fields


@lru_cache(maxsize=1)
def _registry() -> ConditionRegistry:
    return ConditionRegistry()


def _normalized_conditions(record: RawReactionRecord) -> ConditionIdentity:
    return ConditionIdentity(
        catalyst_cas=normalize_cas_list(",".join(record.catalyst_cas)),
        reagent_cas=normalize_cas_list(",".join(record.reagent_cas)),
        solvent_cas=normalize_cas_list(",".join(record.solvent_cas)),
    )


def _condition_resolution(record: RawReactionRecord) -> Dict[str, Any]:
    components = []
    statuses = Counter()
    for source_field in ("catalyst_cas", "reagent_cas", "solvent_cas"):
        for identifier in getattr(record, source_field):
            result = _registry().resolve(cas=identifier)
            statuses[result.status] += 1
            substance = result.substance
            components.append(
                {
                    "source_field": source_field,
                    "raw_identifier": identifier,
                    "status": result.status,
                    "match_kind": result.match_kind,
                    "substance_id": substance.substance_id if substance else None,
                    "canonical_name": substance.canonical_name if substance else None,
                    "possible_roles": tuple(
                        sorted(
                            {
                                assignment.role_id
                                for assignment in (substance.roles if substance else ())
                            }
                        )
                    ),
                }
            )
    return {
        "component_count": len(components),
        "status_counts": dict(sorted(statuses.items())),
        "components": tuple(components),
        "has_uncertainty": any(
            component["status"] != "resolved" for component in components
        ),
        "schema_version": "1.0",
    }


def convert_record(record: RawReactionRecord) -> RecommendationRecord:
    """Convert one source-faithful row without a declared-family requirement."""
    analysis = featurize_reaction(record.reaction_smiles)
    canonical_identity = canonical_reaction_identity(record.reaction_smiles)
    conditions = _normalized_conditions(record)
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
        condition_resolution=_condition_resolution(record),
        source=source,
    )


__all__ = ["convert_record"]
