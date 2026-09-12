"""Refresh canonical records through the normal converter using frozen observations."""

from __future__ import annotations

from typing import Any, Mapping

from condition_registry.models import ConditionComponentInput, ConditionProcessStage
from reactive_taxonomy.observation_serialization import restore_reaction_observation
from reactive_taxonomy.reaction_api import reanalyze_reaction_observation
from reactive_taxonomy.reaction_signatures import reaction_signature_definition_versions
from reactive_taxonomy.reaction_models import (
    REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
    REACTION_CORE_PROJECTION_SCHEMA_VERSION,
)

from .generic import GenericConversionCache, convert_record
from .input_schema import RawReactionRecord


def raw_record_from_canonical(record: Mapping[str, Any]) -> RawReactionRecord:
    """Restore original source inputs, including their existing observation identity."""
    source = record.get("source") or {}
    identifiers = source.get("raw_condition_identifiers")
    if not isinstance(identifiers, Mapping):
        raise ValueError("Canonical refresh requires original condition identifiers")
    return RawReactionRecord(
        source_dataset=str(record["source_dataset"]),
        source_path=str(record["source_path"]),
        source_row_number=int(record["source_row_number"]),
        reaction_id=str(record["reaction_id"]),
        source_declared_family=str(source.get("source_declared_family") or ""),
        reaction_smiles=str(record["reaction_smiles"]),
        yield_pct=record.get("yield_pct"),
        temperature_c=record.get("temperature_c"),
        time_h=record.get("time_h"),
        reference=str(source.get("reference") or ""),
        reactant_cas=tuple(source.get("reactant_cas") or ()),
        product_cas=tuple(source.get("product_cas") or ()),
        catalyst_cas=tuple(identifiers.get("catalyst_cas") or ()),
        reagent_cas=tuple(identifiers.get("reagent_cas") or ()),
        solvent_cas=tuple(identifiers.get("solvent_cas") or ()),
        experimental_procedure=str(source.get("experimental_procedure") or ""),
        stages=str(source.get("stages") or ""),
        steps=str(source.get("steps") or ""),
        notes=str(source.get("notes") or ""),
        condition_component_inputs=tuple(
            ConditionComponentInput(**item)
            for item in source.get("condition_component_inputs") or ()
        ),
        condition_process_stages=tuple(
            ConditionProcessStage(**item)
            for item in source.get("condition_process_stages") or ()
        ),
        condition_declared_absences=tuple(
            source.get("condition_declared_absences") or ()
        ),
        primary_outcome_type=str(source.get("primary_outcome_type") or ""),
        upstream_observation_id=str(record.get("observation_id") or ""),
        raw_fields=dict(source.get("raw_fields") or {}),
    )


def refresh_canonical_record(
    record: Mapping[str, Any], *, cache: GenericConversionCache | None = None
) -> dict[str, Any]:
    """Recompute conditions and annotations without re-inferring current graph evidence.

    The source checksum and exact structural contract are release prerequisites.
    Incompatible observations and externally enriched records require ordinary
    source conversion instead of silently rewriting their structural provenance.
    """
    if record.get("schema_version") not in {"10.1", "10.2", "10.3"}:
        raise ValueError("Unsupported canonical refresh input schema")
    if record.get("external_atom_mapping"):
        raise ValueError(
            "Externally mapped records require their original mapping provider"
        )
    signature = record.get("reaction_signature") or {}
    if (
        signature
        and signature.get("definition_versions")
        != reaction_signature_definition_versions()
    ):
        raise ValueError(
            "Structural definitions changed; full source conversion is required"
        )
    core = record.get("reaction_core") or {}
    if core and (
        core.get("schema_version") != REACTION_CORE_PROJECTION_SCHEMA_VERSION
        or core.get("algorithm_version") != REACTION_CORE_PROJECTION_ALGORITHM_VERSION
    ):
        raise ValueError(
            "Reaction core contract changed; full source conversion is required"
        )
    raw = raw_record_from_canonical(record)
    conversion_cache = cache or GenericConversionCache(max_entries=64)
    value = record.get("reaction_observation")
    if value:
        observation = restore_reaction_observation(value)
        if observation.input_reaction_smiles != raw.reaction_smiles:
            raise ValueError("Stored observation belongs to a different reaction")
        analysis = conversion_cache.get(
            conversion_cache.analyses,
            raw.reaction_smiles,
            lambda: reanalyze_reaction_observation(observation),
        )
        if signature and (
            analysis.reaction_signature is None
            or analysis.reaction_signature.signature_id != signature.get("signature_id")
        ):
            raise ValueError("Reannotation changed frozen reaction signature identity")
    refreshed = convert_record(raw, cache=conversion_cache).to_dict()
    if refreshed.get("observation_id") != record.get("observation_id"):
        raise ValueError("Canonical refresh changed source observation identity")
    return refreshed
