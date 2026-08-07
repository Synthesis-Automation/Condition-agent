"""In-memory indices for canonical generic conversion records."""

from __future__ import annotations

import json
import gzip
from collections import defaultdict
from dataclasses import dataclass, field as dataclass_field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

from reactive_taxonomy import (
    REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
    REACTION_CORE_PROJECTION_SCHEMA_VERSION,
    REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION,
    REACTION_SIGNATURE_SCHEMA_VERSION,
    REACTION_TOPOLOGY_SCHEMA_VERSION,
    reaction_fallback_definition_versions,
    reaction_signature_definition_versions,
)

from .models import (
    CORE_ELIGIBILITY_DEFINITION_VERSION,
    GENERIC_CONVERTER_DEFINITION_VERSION,
    RECOMMENDATION_RECORD_SCHEMA_VERSION,
    CoreEligibility,
    PrecedentIndexScope,
    PrecedentTier,
)
from .evaluation_features import reaction_scaffold_key, reaction_scaffold_tokens
from .signature_features import environment_tokens
from .fallback_similarity import fallback_index_tokens
from .reaction_facets import reaction_facet_keys


GENERIC_INDEX_SCHEMA_VERSION = "6.1"


@dataclass(frozen=True)
class GenericIndexedReaction:
    """A verified precedent reduced to fields required for retrieval."""

    reaction_id: str
    observation_id: str
    canonical_reaction_id: str
    reaction_smiles: str
    yield_pct: Optional[float]
    source_dataset: str
    reference_id: str
    publication_year: Optional[int]
    reference_condition_series_id: str
    scaffold_key: str
    scaffold_tokens: Tuple[str, ...]
    signature: Dict[str, Any]
    reaction_core: Dict[str, Any]
    recipe_id: str
    recipe_core_id: str
    resolved_recipe: Dict[str, Any]
    condition_uncertain: bool
    chemistry_status: str
    condition_status: str
    condition_stage_status: str
    outcome_status: str
    record_schema_version: str
    converter_definition_version: str
    precedent_tier: PrecedentTier
    core_eligibility_definition_version: str
    reaction_label: Dict[str, Any] = dataclass_field(default_factory=dict)
    fallback_descriptor: Dict[str, Any] = dataclass_field(default_factory=dict)
    fragment_source_support: Tuple[Dict[str, Any], ...] = ()

    @property
    def named_family(self) -> str:
        return str(self.signature.get("named_family") or "")

    @property
    def family_confidence(self) -> float:
        return float(self.signature.get("family_confidence") or 0.0)

    @property
    def transformation_class(self) -> str:
        return str(self.signature.get("transformation_class") or "")


@dataclass(frozen=True)
class GenericReactionIndex:
    """Immutable precedent collection with deterministic signature lookup maps."""

    rows: Sequence[GenericIndexedReaction]
    exact: Mapping[str, Tuple[int, ...]]
    handles: Mapping[str, Tuple[int, ...]]
    transformations: Mapping[str, Tuple[int, ...]]
    bond_edits: Mapping[str, Tuple[int, ...]]
    environments: Mapping[str, Tuple[int, ...]]
    facet_exact: Mapping[str, Tuple[int, ...]]
    facet_attachments: Mapping[str, Tuple[int, ...]]
    core_exact: Mapping[str, Tuple[int, ...]]
    core_typed: Mapping[str, Tuple[int, ...]]
    core_shapes: Mapping[str, Tuple[int, ...]]
    core_centers: Mapping[str, Tuple[int, ...]]
    environment_features: Mapping[str, Tuple[int, ...]]
    fallback_features: Mapping[str, Tuple[int, ...]]
    partial_transformations: Mapping[str, Tuple[int, ...]]
    families: Mapping[str, Tuple[int, ...]]
    reaction_signature_schema_version: str
    reaction_core_schema_version: str
    reaction_core_algorithm_version: str
    taxonomy_definition_versions: Tuple[Tuple[str, str], ...]
    fallback_descriptor_schema_version: str
    fallback_definition_versions: Tuple[Tuple[str, str], ...]
    record_schema_versions: Tuple[str, ...]
    converter_definition_versions: Tuple[str, ...]
    precedent_scope: PrecedentIndexScope
    core_eligibility_definition_version: str

    def select(self, positions: Iterable[int]) -> Tuple[GenericIndexedReaction, ...]:
        bulk_select = getattr(self.rows, "select", None)
        if callable(bulk_select):
            return tuple(bulk_select(positions))
        return tuple(self.rows[position] for position in positions)


_KEY_FIELDS = {
    "exact": "exact_signature_key",
    "handles": "handle_signature_key",
    "transformations": "transformation_signature_key",
    "bond_edits": "bond_edit_signature_key",
    "environments": "environment_signature_key",
}


def _freeze(mapping: Dict[str, list[int]]) -> Dict[str, Tuple[int, ...]]:
    return {key: tuple(values) for key, values in sorted(mapping.items())}


def _definition_version_tuple(
    signature: Mapping[str, Any],
) -> Tuple[Tuple[str, str], ...]:
    values = signature.get("definition_versions") or {}
    if not isinstance(values, Mapping):
        raise ValueError("Reaction signature definition_versions must be an object")
    return tuple(sorted((str(key), str(value)) for key, value in values.items()))


def _validate_index_rows(
    rows: Iterable[GenericIndexedReaction],
) -> tuple[
    str,
    Tuple[Tuple[str, str], ...],
    str,
    str,
    Tuple[str, ...],
    Tuple[str, ...],
]:
    values = tuple(rows)
    signature_schemas = {
        str(row.signature.get("schema_version") or "")
        for row in values
        if row.signature
    }
    if signature_schemas and signature_schemas != {REACTION_SIGNATURE_SCHEMA_VERSION}:
        raise ValueError(
            "Incompatible reaction signature schema; regenerate converted records"
        )
    definition_sets = {
        _definition_version_tuple(row.signature) for row in values if row.signature
    }
    current_definitions = tuple(
        sorted(reaction_signature_definition_versions().items())
    )
    if definition_sets and definition_sets != {current_definitions}:
        raise ValueError(
            "Incompatible reaction taxonomy definitions; regenerate converted records"
        )
    topology_schemas = {
        str((row.signature.get("topology") or {}).get("schema_version") or "")
        for row in values
        if row.signature and row.signature.get("topology")
    }
    if topology_schemas and topology_schemas != {
        REACTION_TOPOLOGY_SCHEMA_VERSION
    }:
        raise ValueError(
            "Incompatible reaction topology schema; regenerate converted records"
        )
    core_schemas = {
        str(row.reaction_core.get("schema_version") or "")
        for row in values
        if row.reaction_core
    }
    if core_schemas and core_schemas != {
        REACTION_CORE_PROJECTION_SCHEMA_VERSION
    }:
        raise ValueError(
            "Incompatible reaction core schema; regenerate converted records"
        )
    core_algorithms = {
        str(row.reaction_core.get("algorithm_version") or "")
        for row in values
        if row.reaction_core
    }
    if core_algorithms and core_algorithms != {
        REACTION_CORE_PROJECTION_ALGORITHM_VERSION
    }:
        raise ValueError(
            "Incompatible reaction core algorithm; regenerate converted records"
        )
    fallback_descriptors = tuple(
        row.fallback_descriptor for row in values if row.fallback_descriptor
    )
    fallback_schemas = {
        str(descriptor.get("schema_version") or "")
        for descriptor in fallback_descriptors
    }
    if fallback_schemas and fallback_schemas != {
        REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION
    }:
        raise ValueError(
            "Incompatible fallback descriptor schema; regenerate converted records"
        )
    current_fallback_definitions = tuple(
        sorted(reaction_fallback_definition_versions().items())
    )
    fallback_definition_sets = {
        tuple(
            sorted(
                (str(key), str(value))
                for key, value in (descriptor.get("definition_versions") or {}).items()
            )
        )
        for descriptor in fallback_descriptors
    }
    if fallback_definition_sets and fallback_definition_sets != {
        current_fallback_definitions
    }:
        raise ValueError(
            "Incompatible fallback descriptor definitions; regenerate converted records"
        )
    record_schemas = tuple(sorted({row.record_schema_version for row in values}))
    if record_schemas and record_schemas != (
        RECOMMENDATION_RECORD_SCHEMA_VERSION,
    ):
        raise ValueError(
            "Incompatible recommendation record schema; regenerate converted records"
        )
    converter_versions = tuple(
        sorted({row.converter_definition_version for row in values})
    )
    if converter_versions and converter_versions != (
        GENERIC_CONVERTER_DEFINITION_VERSION,
    ):
        raise ValueError(
            "Incompatible generic converter version; regenerate converted records"
        )
    core_eligibility_versions = {
        row.core_eligibility_definition_version for row in values
    }
    if core_eligibility_versions and core_eligibility_versions != {
        CORE_ELIGIBILITY_DEFINITION_VERSION
    }:
        raise ValueError(
            "Incompatible core eligibility definition; regenerate converted records"
        )
    return (
        REACTION_SIGNATURE_SCHEMA_VERSION,
        current_definitions,
        REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
        record_schemas or (RECOMMENDATION_RECORD_SCHEMA_VERSION,),
        converter_versions or (GENERIC_CONVERTER_DEFINITION_VERSION,),
    )


def build_generic_index_from_rows(
    rows: Iterable[GenericIndexedReaction],
    *,
    precedent_scope: PrecedentIndexScope = PrecedentIndexScope.TRUSTED,
) -> GenericReactionIndex:
    """Build deterministic lookup maps from already validated index rows."""
    ordered = sorted(
        rows,
        key=lambda row: (
            row.canonical_reaction_id,
            row.reaction_id,
            row.observation_id,
            row.recipe_id,
        ),
    )
    (
        signature_schema,
        definition_versions,
        core_schema,
        core_algorithm,
        record_schemas,
        converter_versions,
    ) = _validate_index_rows(ordered)
    row_tiers = {row.precedent_tier for row in ordered}
    if precedent_scope == PrecedentIndexScope.TRUSTED:
        if row_tiers - {PrecedentTier.TRUSTED}:
            raise ValueError("Trusted index cannot contain review-core precedents")
    elif row_tiers - {PrecedentTier.TRUSTED, PrecedentTier.REVIEW_CORE}:
        raise ValueError("Unsupported precedent tier")
    maps: Dict[str, Dict[str, list[int]]] = {
        name: defaultdict(list) for name in _KEY_FIELDS
    }
    families: Dict[str, list[int]] = defaultdict(list)
    core_maps: Dict[str, Dict[str, list[int]]] = {
        name: defaultdict(list)
        for name in ("exact", "typed", "shapes", "centers")
    }
    environment_features: Dict[str, list[int]] = defaultdict(list)
    fallback_features: Dict[str, list[int]] = defaultdict(list)
    partial_transformations: Dict[str, list[int]] = defaultdict(list)
    facet_maps: Dict[str, Dict[str, list[int]]] = {
        "exact": defaultdict(list),
        "attachments": defaultdict(list),
    }
    for position, row in enumerate(ordered):
        for name, field in _KEY_FIELDS.items():
            key = str(row.signature.get(field) or "")
            if key:
                maps[name][key].append(position)
        core_fields = {
            "exact": "exact_core_key",
            "typed": "typed_core_key",
            "shapes": "shape_core_key",
            "centers": "center_transition_key",
        }
        for name, field in core_fields.items():
            key = str(row.reaction_core.get(field) or "")
            if key:
                core_maps[name][key].append(position)
        if row.named_family:
            families[row.named_family].append(position)
        for token in set(environment_tokens(row.signature)):
            environment_features[token].append(position)
        for token in fallback_index_tokens(row.fallback_descriptor):
            fallback_features[token].append(position)
        facet_keys = reaction_facet_keys(
            row.signature,
            row.reaction_core,
            row.fallback_descriptor,
        )
        exact_facet = facet_keys.get("reaction_facet_exact")
        if exact_facet:
            facet_maps["exact"][exact_facet].append(position)
        attachment_facet = facet_keys.get(
            "reaction_facet_attachment_relaxed"
        )
        if attachment_facet:
            facet_maps["attachments"][attachment_facet].append(position)
        partial_key = str(
            row.fallback_descriptor.get("partial_transformation_key") or ""
        )
        if partial_key:
            partial_transformations[partial_key].append(position)
    return GenericReactionIndex(
        rows=tuple(ordered),
        exact=_freeze(maps["exact"]),
        handles=_freeze(maps["handles"]),
        transformations=_freeze(maps["transformations"]),
        bond_edits=_freeze(maps["bond_edits"]),
        environments=_freeze(maps["environments"]),
        facet_exact=_freeze(facet_maps["exact"]),
        facet_attachments=_freeze(facet_maps["attachments"]),
        core_exact=_freeze(core_maps["exact"]),
        core_typed=_freeze(core_maps["typed"]),
        core_shapes=_freeze(core_maps["shapes"]),
        core_centers=_freeze(core_maps["centers"]),
        environment_features=_freeze(environment_features),
        fallback_features=_freeze(fallback_features),
        partial_transformations=_freeze(partial_transformations),
        families=_freeze(families),
        reaction_signature_schema_version=signature_schema,
        reaction_core_schema_version=core_schema,
        reaction_core_algorithm_version=core_algorithm,
        taxonomy_definition_versions=definition_versions,
        fallback_descriptor_schema_version=(
            REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION
        ),
        fallback_definition_versions=tuple(
            sorted(reaction_fallback_definition_versions().items())
        ),
        record_schema_versions=record_schemas,
        converter_definition_versions=converter_versions,
        precedent_scope=precedent_scope,
        core_eligibility_definition_version=(
            CORE_ELIGIBILITY_DEFINITION_VERSION
        ),
    )


def build_generic_index(
    records: Iterable[Mapping[str, Any]],
    *,
    include_review: bool = False,
) -> GenericReactionIndex:
    """Build lookup maps, admitting only records with signatures and recipes."""
    return build_generic_index_from_rows(
        _iter_generic_index_rows(records, include_review=include_review),
        precedent_scope=(
            PrecedentIndexScope.TRUSTED_AND_REVIEW_CORE
            if include_review
            else PrecedentIndexScope.TRUSTED
        ),
    )


def _iter_generic_index_rows(
    records: Iterable[Mapping[str, Any]],
    *,
    include_review: bool = False,
) -> Iterable[GenericIndexedReaction]:
    """Yield admitted retrieval rows without retaining the input corpus."""
    for record in records:
        fragment_source_support = tuple(
            dict(value)
            for value in record.get("fragment_source_support") or ()
            if isinstance(value, Mapping)
        )
        eligibility = _enum_value(record.get("index_eligibility"))
        core_eligibility = _enum_value(record.get("core_eligibility"))
        core_eligibility_definition_version = str(
            record.get("core_eligibility_definition_version") or ""
        )
        if (
            core_eligibility_definition_version
            != CORE_ELIGIBILITY_DEFINITION_VERSION
        ):
            raise ValueError(
                "Converted record lacks current core eligibility definition; "
                "regenerate converted records"
            )
        persisted_precedent_tier = _enum_value(record.get("precedent_tier"))
        if core_eligibility not in {value.value for value in CoreEligibility}:
            raise ValueError(
                "Converted record lacks current core eligibility; regenerate "
                "converted records"
            )
        if (
            eligibility == "eligible"
            and persisted_precedent_tier == PrecedentTier.TRUSTED.value
        ):
            precedent_tier = PrecedentTier.TRUSTED
        elif (
            include_review
            and eligibility == "review_only"
            and core_eligibility == CoreEligibility.REVIEW_CORE.value
            and persisted_precedent_tier == PrecedentTier.REVIEW_CORE.value
        ):
            precedent_tier = PrecedentTier.REVIEW_CORE
        else:
            continue
        signature = record.get("reaction_signature")
        reaction_core = record.get("reaction_core")
        fallback_descriptor = record.get("fallback_descriptor")
        recipe = record.get("resolved_recipe")
        recipe_id = str(record.get("resolved_recipe_id") or "")
        recipe_core_id = str(
            record.get("resolved_recipe_core_id")
            or (recipe or {}).get("recipe_core_id")
            or recipe_id
        )
        outcome = record.get("yield_pct")
        if (
            not isinstance(signature, Mapping)
            and not isinstance(fallback_descriptor, Mapping)
        ) or not isinstance(recipe, Mapping):
            continue
        if not recipe_id or not recipe_core_id:
            continue
        yield_pct = _valid_yield(outcome)
        if str(recipe.get("recipe_id") or "") != recipe_id:
            continue
        embedded_core_id = str(recipe.get("recipe_core_id") or recipe_core_id)
        if embedded_core_id != recipe_core_id:
            continue
        yield GenericIndexedReaction(
            reaction_id=str(record.get("reaction_id") or ""),
            observation_id=str(record.get("observation_id") or ""),
            canonical_reaction_id=str(
                record.get("canonical_reaction_id")
                or record.get("reaction_id")
                or record.get("observation_id")
                or ""
            ),
            reaction_smiles=str(record.get("reaction_smiles") or ""),
            yield_pct=yield_pct,
            source_dataset=str(record.get("source_dataset") or ""),
            reference_id=str(record.get("reference_id") or ""),
            publication_year=(
                int((record.get("reference_identity") or {})["publication_year"])
                if (record.get("reference_identity") or {}).get("publication_year")
                is not None
                else None
            ),
            reference_condition_series_id=str(
                record.get("reference_condition_series_id") or ""
            ),
            scaffold_key=reaction_scaffold_key(
                str(record.get("reaction_smiles") or ""),
                signature if isinstance(signature, Mapping) else {},
            ),
            scaffold_tokens=reaction_scaffold_tokens(
                str(record.get("reaction_smiles") or ""),
                signature if isinstance(signature, Mapping) else {},
            ),
            signature=(dict(signature) if isinstance(signature, Mapping) else {}),
            # A blocked core must never enter core lookup maps. A separately
            # verified signature may remain a trusted precedent through its
            # signature keys, so only the unsafe core view is removed.
            reaction_core=(
                dict(reaction_core)
                if isinstance(reaction_core, Mapping)
                and (
                    core_eligibility == CoreEligibility.TRUSTED_CORE.value
                    or (
                        include_review
                        and core_eligibility == CoreEligibility.REVIEW_CORE.value
                    )
                )
                else {}
            ),
            recipe_id=recipe_id,
            recipe_core_id=recipe_core_id,
            resolved_recipe=dict(recipe),
            condition_uncertain=bool(
                (record.get("condition_resolution") or {}).get("has_uncertainty")
            )
            or _enum_value(record.get("condition_stage_status"))
            == "unassigned_multistage",
            chemistry_status=_enum_value(record.get("chemistry_status")),
            condition_status=_enum_value(record.get("condition_status")),
            condition_stage_status=(
                _enum_value(record.get("condition_stage_status")) or "single_stage"
            ),
            outcome_status=_enum_value(record.get("outcome_status")),
            record_schema_version=str(record.get("schema_version") or ""),
            converter_definition_version=str(
                record.get("converter_definition_version") or ""
            ),
            precedent_tier=precedent_tier,
            core_eligibility_definition_version=(
                core_eligibility_definition_version
            ),
            reaction_label=(
                dict(record.get("reaction_label") or {})
                if isinstance(record.get("reaction_label"), Mapping)
                else {}
            ),
            fallback_descriptor=(
                dict(fallback_descriptor or {})
                if isinstance(fallback_descriptor, Mapping)
                else {}
            ),
            fragment_source_support=fragment_source_support,
        )


def _enum_value(value: Any) -> str:
    return str(value.value if isinstance(value, Enum) else value or "")


def _valid_yield(value: Any) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        outcome = float(value)
    except (TypeError, ValueError):
        return None
    return outcome if 0.0 <= outcome <= 100.0 else None


def _indexed_reaction_payload(row: GenericIndexedReaction) -> Dict[str, Any]:
    """Serialize one retrieval row without storage-format assumptions."""
    return {
        "reaction_id": row.reaction_id,
        "observation_id": row.observation_id,
        "canonical_reaction_id": row.canonical_reaction_id,
        "reaction_smiles": row.reaction_smiles,
        "yield_pct": row.yield_pct,
        "source_dataset": row.source_dataset,
        "reference_id": row.reference_id,
        "publication_year": row.publication_year,
        "reference_condition_series_id": row.reference_condition_series_id,
        "scaffold_key": row.scaffold_key,
        "scaffold_tokens": row.scaffold_tokens,
        "signature": row.signature,
        "reaction_core": row.reaction_core,
        "recipe_id": row.recipe_id,
        "recipe_core_id": row.recipe_core_id,
        "resolved_recipe": row.resolved_recipe,
        "condition_uncertain": row.condition_uncertain,
        "chemistry_status": row.chemistry_status,
        "condition_status": row.condition_status,
        "condition_stage_status": row.condition_stage_status,
        "outcome_status": row.outcome_status,
        "record_schema_version": row.record_schema_version,
        "converter_definition_version": row.converter_definition_version,
        "precedent_tier": row.precedent_tier.value,
        "core_eligibility_definition_version": (
            row.core_eligibility_definition_version
        ),
        "reaction_label": row.reaction_label,
        "fallback_descriptor": row.fallback_descriptor,
        "fragment_source_support": row.fragment_source_support,
    }


def _indexed_reaction_from_payload(
    row: Mapping[str, Any],
) -> GenericIndexedReaction:
    """Deserialize one retrieval row shared by JSON and SQLite storage."""
    return GenericIndexedReaction(
        reaction_id=str(row["reaction_id"]),
        observation_id=str(row["observation_id"]),
        canonical_reaction_id=str(row["canonical_reaction_id"]),
        reaction_smiles=str(row["reaction_smiles"]),
        yield_pct=(
            float(row["yield_pct"]) if row.get("yield_pct") is not None else None
        ),
        source_dataset=str(row["source_dataset"]),
        reference_id=str(row.get("reference_id") or ""),
        publication_year=(
            int(row["publication_year"])
            if row.get("publication_year") is not None
            else None
        ),
        reference_condition_series_id=str(
            row.get("reference_condition_series_id") or ""
        ),
        scaffold_key=str(row.get("scaffold_key") or ""),
        scaffold_tokens=tuple(
            str(value) for value in row.get("scaffold_tokens") or ()
        ),
        signature=dict(row["signature"]),
        reaction_core=dict(row.get("reaction_core") or {}),
        recipe_id=str(row["recipe_id"]),
        recipe_core_id=str(row.get("recipe_core_id") or row["recipe_id"]),
        resolved_recipe=dict(row["resolved_recipe"]),
        condition_uncertain=bool(row["condition_uncertain"]),
        chemistry_status=str(row.get("chemistry_status") or ""),
        condition_status=str(row.get("condition_status") or ""),
        condition_stage_status=str(
            row.get("condition_stage_status") or "single_stage"
        ),
        outcome_status=str(row.get("outcome_status") or ""),
        record_schema_version=str(row["record_schema_version"]),
        converter_definition_version=str(row["converter_definition_version"]),
        precedent_tier=PrecedentTier(str(row["precedent_tier"])),
        core_eligibility_definition_version=str(
            row["core_eligibility_definition_version"]
        ),
        reaction_label=dict(row.get("reaction_label") or {}),
        fallback_descriptor=dict(row.get("fallback_descriptor") or {}),
        fragment_source_support=tuple(
            dict(value)
            for value in row.get("fragment_source_support") or ()
            if isinstance(value, Mapping)
        ),
    )


def _index_maps(
    index: GenericReactionIndex,
) -> Dict[str, Mapping[str, Tuple[int, ...]]]:
    """Return all deterministic lookup maps by their persisted names."""
    return {
        "exact": index.exact,
        "handles": index.handles,
        "transformations": index.transformations,
        "bond_edits": index.bond_edits,
        "environments": index.environments,
        "facet_exact": index.facet_exact,
        "facet_attachments": index.facet_attachments,
        "core_exact": index.core_exact,
        "core_typed": index.core_typed,
        "core_shapes": index.core_shapes,
        "core_centers": index.core_centers,
        "environment_features": index.environment_features,
        "fallback_features": index.fallback_features,
        "partial_transformations": index.partial_transformations,
        "families": index.families,
    }


def _validate_index_metadata(payload: Mapping[str, Any]) -> PrecedentIndexScope:
    """Validate the chemistry and schema contract of a persisted index."""
    if payload.get("artifact_type") != "generic_reaction_index":
        raise ValueError("Not a generic reaction index artifact")
    if payload.get("schema_version") != GENERIC_INDEX_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported generic reaction index schema; rebuild recommendation "
            "artifacts and rebuild the index to create clean trusted and "
            "review-core indexes"
        )
    try:
        precedent_scope = PrecedentIndexScope(str(payload["precedent_scope"]))
    except (KeyError, ValueError) as exc:
        raise ValueError("Invalid precedent index scope; rebuild the index") from exc
    if (
        payload.get("core_eligibility_definition_version")
        != CORE_ELIGIBILITY_DEFINITION_VERSION
    ):
        raise ValueError(
            "Incompatible core eligibility definition; rebuild the index"
        )
    if (
        payload.get("reaction_signature_schema_version")
        != REACTION_SIGNATURE_SCHEMA_VERSION
    ):
        raise ValueError("Incompatible reaction signature schema; rebuild the index")
    if (
        payload.get("reaction_core_schema_version")
        != REACTION_CORE_PROJECTION_SCHEMA_VERSION
    ):
        raise ValueError("Incompatible reaction core schema; rebuild the index")
    if (
        payload.get("reaction_core_algorithm_version")
        != REACTION_CORE_PROJECTION_ALGORITHM_VERSION
    ):
        raise ValueError("Incompatible reaction core algorithm; rebuild the index")
    current_definitions = reaction_signature_definition_versions()
    if payload.get("taxonomy_definition_versions") != current_definitions:
        raise ValueError(
            "Incompatible reaction taxonomy definitions; rebuild the index"
        )
    if (
        payload.get("fallback_descriptor_schema_version")
        != REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION
    ):
        raise ValueError("Incompatible fallback descriptor schema; rebuild the index")
    if payload.get("fallback_definition_versions") != (
        reaction_fallback_definition_versions()
    ):
        raise ValueError(
            "Incompatible fallback descriptor definitions; rebuild the index"
        )
    record_schema_versions = tuple(payload.get("record_schema_versions") or ())
    if record_schema_versions != (RECOMMENDATION_RECORD_SCHEMA_VERSION,):
        raise ValueError("Incompatible recommendation record schema; rebuild the index")
    converter_definition_versions = tuple(
        payload.get("converter_definition_versions") or ()
    )
    if converter_definition_versions != (GENERIC_CONVERTER_DEFINITION_VERSION,):
        raise ValueError("Incompatible generic converter version; rebuild the index")
    return precedent_scope


def load_generic_index(
    path: str | Path,
    *,
    include_review: bool = False,
) -> GenericReactionIndex:
    """Load canonical JSONL output from the generic conversion engine."""
    source = Path(path)
    if source.suffix.casefold() in {".sqlite", ".sqlite3", ".db"}:
        from .sqlite_indexing import load_sqlite_generic_index

        return load_sqlite_generic_index(source)
    if source.name.casefold().endswith(".json.gz"):
        raise ValueError(
            "Persisted JSON recommendation indexes are retired; rebuild or "
            "select the corresponding SQLite index"
        )
    if source.suffix.casefold() == ".json":
        if source.name == "shard_manifest.json":
            from .conversion.sharded import (
                iter_gzip_jsonl,
                validate_sharded_conversion,
            )

            integrity = validate_sharded_conversion(source, verify_rows=False)
            if not integrity["valid"]:
                raise ValueError("Sharded conversion integrity check failed")
            manifest = json.loads(source.read_text(encoding="utf-8"))
            records = (
                record
                for entry in manifest.get("shards") or ()
                for record in iter_gzip_jsonl(source.parent / entry["output_path"])
            )
            return build_generic_index(records, include_review=include_review)
        raise ValueError(
            "Only shard_manifest.json is accepted as canonical JSON input; "
            "persisted runtime indexes must use SQLite"
        )
    opener = gzip.open if source.suffix.casefold() == ".gz" else Path.open
    open_arguments = (
        {"mode": "rt", "encoding": "utf-8"}
        if source.suffix.casefold() == ".gz"
        else {"mode": "r", "encoding": "utf-8"}
    )

    def records() -> Iterable[Dict[str, Any]]:
        with opener(source, **open_arguments) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    value = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ValueError(
                        f"Invalid JSONL at line {line_number}: {exc.msg}"
                    ) from exc
                if not isinstance(value, dict):
                    raise ValueError(f"JSONL line {line_number} is not an object")
                yield value

    return build_generic_index(records(), include_review=include_review)


def validate_generic_index_artifact(path: str | Path) -> Dict[str, Any]:
    """Validate the supported SQLite runtime-index artifact."""
    source = Path(path)
    if source.suffix.casefold() not in {".sqlite", ".sqlite3", ".db"}:
        raise ValueError(
            "Runtime index integrity validation supports SQLite only; canonical "
            "JSONL is validated through the conversion manifest"
        )
    from .sqlite_indexing import validate_sqlite_generic_index

    return validate_sqlite_generic_index(source)


__all__ = [
    "GenericIndexedReaction",
    "GenericReactionIndex",
    "build_generic_index",
    "build_generic_index_from_rows",
    "load_generic_index",
    "validate_generic_index_artifact",
]
