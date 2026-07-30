"""In-memory indices for canonical generic conversion records."""

from __future__ import annotations

import json
import hashlib
import gzip
import io
from collections import Counter, defaultdict
from dataclasses import dataclass, field as dataclass_field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple

from reactive_taxonomy import (
    REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION,
    REACTION_SIGNATURE_SCHEMA_VERSION,
    reaction_fallback_definition_versions,
    reaction_signature_definition_versions,
)

from .models import (
    GENERIC_CONVERTER_DEFINITION_VERSION,
    RECOMMENDATION_RECORD_SCHEMA_VERSION,
)
from .evaluation_features import reaction_scaffold_key, reaction_scaffold_tokens
from .signature_features import environment_tokens
from .fallback_similarity import fallback_index_tokens


GENERIC_INDEX_SCHEMA_VERSION = "2.3"

# Record schema 3.6 / converter v2.6 add optional external-mapping review
# provenance only. They do not change verified signature, fallback descriptor,
# recipe, or index-row contracts, so persisted v3.4/v2.4 and v3.5/v2.5 indices
# remain chemistry-compatible.
_INDEX_COMPATIBLE_RECORD_SCHEMAS = frozenset(
    {"3.4", "3.5", RECOMMENDATION_RECORD_SCHEMA_VERSION}
)
_INDEX_COMPATIBLE_CONVERTER_VERSIONS = frozenset(
    {
        "generic_conversion.v2.4",
        "generic_conversion.v2.5",
        GENERIC_CONVERTER_DEFINITION_VERSION,
    }
)


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
    reaction_label: str = ""
    reaction_label_status: str = "unavailable"
    fallback_descriptor: Dict[str, Any] = dataclass_field(default_factory=dict)

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

    rows: Tuple[GenericIndexedReaction, ...]
    exact: Mapping[str, Tuple[int, ...]]
    handles: Mapping[str, Tuple[int, ...]]
    transformations: Mapping[str, Tuple[int, ...]]
    bond_edits: Mapping[str, Tuple[int, ...]]
    environments: Mapping[str, Tuple[int, ...]]
    environment_features: Mapping[str, Tuple[int, ...]]
    fallback_features: Mapping[str, Tuple[int, ...]]
    families: Mapping[str, Tuple[int, ...]]
    reaction_signature_schema_version: str
    taxonomy_definition_versions: Tuple[Tuple[str, str], ...]
    fallback_descriptor_schema_version: str
    fallback_definition_versions: Tuple[Tuple[str, str], ...]
    record_schema_versions: Tuple[str, ...]
    converter_definition_versions: Tuple[str, ...]

    def select(self, positions: Iterable[int]) -> Tuple[GenericIndexedReaction, ...]:
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
    if record_schemas and not set(record_schemas).issubset(
        _INDEX_COMPATIBLE_RECORD_SCHEMAS
    ):
        raise ValueError(
            "Incompatible recommendation record schema; regenerate converted records"
        )
    converter_versions = tuple(
        sorted({row.converter_definition_version for row in values})
    )
    if converter_versions and not set(converter_versions).issubset(
        _INDEX_COMPATIBLE_CONVERTER_VERSIONS
    ):
        raise ValueError(
            "Incompatible generic converter version; regenerate converted records"
        )
    return (
        REACTION_SIGNATURE_SCHEMA_VERSION,
        current_definitions,
        record_schemas or (RECOMMENDATION_RECORD_SCHEMA_VERSION,),
        converter_versions or (GENERIC_CONVERTER_DEFINITION_VERSION,),
    )


def build_generic_index_from_rows(
    rows: Iterable[GenericIndexedReaction],
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
        record_schemas,
        converter_versions,
    ) = _validate_index_rows(ordered)
    maps: Dict[str, Dict[str, list[int]]] = {
        name: defaultdict(list) for name in _KEY_FIELDS
    }
    families: Dict[str, list[int]] = defaultdict(list)
    environment_features: Dict[str, list[int]] = defaultdict(list)
    fallback_features: Dict[str, list[int]] = defaultdict(list)
    for position, row in enumerate(ordered):
        for name, field in _KEY_FIELDS.items():
            key = str(row.signature.get(field) or "")
            if key:
                maps[name][key].append(position)
        if row.named_family:
            families[row.named_family].append(position)
        for token in set(environment_tokens(row.signature)):
            environment_features[token].append(position)
        for token in fallback_index_tokens(row.fallback_descriptor):
            fallback_features[token].append(position)
    return GenericReactionIndex(
        rows=tuple(ordered),
        exact=_freeze(maps["exact"]),
        handles=_freeze(maps["handles"]),
        transformations=_freeze(maps["transformations"]),
        bond_edits=_freeze(maps["bond_edits"]),
        environments=_freeze(maps["environments"]),
        environment_features=_freeze(environment_features),
        fallback_features=_freeze(fallback_features),
        families=_freeze(families),
        reaction_signature_schema_version=signature_schema,
        taxonomy_definition_versions=definition_versions,
        fallback_descriptor_schema_version=(
            REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION
        ),
        fallback_definition_versions=tuple(
            sorted(reaction_fallback_definition_versions().items())
        ),
        record_schema_versions=record_schemas,
        converter_definition_versions=converter_versions,
    )


def build_generic_index(
    records: Iterable[Mapping[str, Any]],
    *,
    include_review: bool = False,
) -> GenericReactionIndex:
    """Build lookup maps, admitting only records with signatures and recipes."""
    rows = []
    for record in records:
        external_mapping = record.get("external_atom_mapping")
        external_mapping_status = (
            str(external_mapping.get("status") or "")
            if isinstance(external_mapping, Mapping)
            else ""
        )
        if (
            not include_review
            and external_mapping_status
            and not external_mapping_status.startswith("not_requested_")
        ):
            continue
        eligibility = _enum_value(record.get("index_eligibility"))
        if eligibility != "eligible" and not (
            include_review and eligibility == "review_only"
        ):
            continue
        signature = record.get("reaction_signature")
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
        rows.append(
            GenericIndexedReaction(
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
                reaction_label=str(record.get("reaction_label") or ""),
                reaction_label_status=str(
                    record.get("reaction_label_status") or "unavailable"
                ),
                fallback_descriptor=(
                    dict(fallback_descriptor or {})
                    if isinstance(fallback_descriptor, Mapping)
                    else {}
                ),
            )
        )
    return build_generic_index_from_rows(rows)


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


def _index_payload(index: GenericReactionIndex) -> Dict[str, Any]:
    rows = [
        {
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
            "reaction_label": row.reaction_label,
            "reaction_label_status": row.reaction_label_status,
            "fallback_descriptor": row.fallback_descriptor,
        }
        for row in index.rows
    ]
    maps = {
        "exact": dict(index.exact),
        "handles": dict(index.handles),
        "transformations": dict(index.transformations),
        "bond_edits": dict(index.bond_edits),
        "environments": dict(index.environments),
        "environment_features": dict(index.environment_features),
        "fallback_features": dict(index.fallback_features),
        "families": dict(index.families),
    }
    identity = json.dumps(
        {
            "rows": rows,
            "maps": maps,
            "reaction_signature_schema_version": (
                index.reaction_signature_schema_version
            ),
            "taxonomy_definition_versions": dict(index.taxonomy_definition_versions),
            "fallback_descriptor_schema_version": (
                index.fallback_descriptor_schema_version
            ),
            "fallback_definition_versions": dict(index.fallback_definition_versions),
            "record_schema_versions": index.record_schema_versions,
            "converter_definition_versions": index.converter_definition_versions,
        },
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    )
    return {
        "schema_version": GENERIC_INDEX_SCHEMA_VERSION,
        "artifact_type": "generic_reaction_index",
        "reaction_signature_schema_version": (index.reaction_signature_schema_version),
        "taxonomy_definition_versions": dict(index.taxonomy_definition_versions),
        "fallback_descriptor_schema_version": (
            index.fallback_descriptor_schema_version
        ),
        "fallback_definition_versions": dict(index.fallback_definition_versions),
        "record_schema_versions": index.record_schema_versions,
        "converter_definition_versions": index.converter_definition_versions,
        "index_id": "GRI1:" + hashlib.sha256(identity.encode("utf-8")).hexdigest(),
        "row_count": len(rows),
        "rows": rows,
        "maps": maps,
    }


def save_generic_index(index: GenericReactionIndex, path: str | Path) -> Dict[str, Any]:
    """Write a deterministic, versioned generic index artifact."""
    payload = _index_payload(index)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    completed = False
    try:
        if destination.suffix.casefold() == ".gz":
            with temporary.open("wb") as raw:
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    fileobj=raw,
                    compresslevel=6,
                    mtime=0,
                ) as compressed:
                    with io.TextIOWrapper(
                        compressed,
                        encoding="utf-8",
                    ) as handle:
                        json.dump(
                            payload,
                            handle,
                            ensure_ascii=False,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                        handle.write("\n")
        else:
            with temporary.open("w", encoding="utf-8") as handle:
                json.dump(
                    payload,
                    handle,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                )
                handle.write("\n")
        temporary.replace(destination)
        completed = True
    finally:
        if not completed and temporary.is_file():
            temporary.unlink()
    return {
        "schema_version": payload["schema_version"],
        "reaction_signature_schema_version": payload[
            "reaction_signature_schema_version"
        ],
        "taxonomy_definition_versions": payload["taxonomy_definition_versions"],
        "fallback_descriptor_schema_version": payload[
            "fallback_descriptor_schema_version"
        ],
        "fallback_definition_versions": payload["fallback_definition_versions"],
        "record_schema_versions": payload["record_schema_versions"],
        "converter_definition_versions": payload["converter_definition_versions"],
        "index_id": payload["index_id"],
        "row_count": payload["row_count"],
        "path": str(destination),
    }


def load_persisted_generic_index(path: str | Path) -> GenericReactionIndex:
    """Load and validate a persisted generic index without rebuilding maps."""
    source = Path(path)
    opener = gzip.open if source.suffix.casefold() == ".gz" else Path.open
    arguments = (
        {"mode": "rt", "encoding": "utf-8"}
        if source.suffix.casefold() == ".gz"
        else {"mode": "r", "encoding": "utf-8"}
    )
    with opener(source, **arguments) as handle:
        payload = json.load(handle)
    if payload.get("artifact_type") != "generic_reaction_index":
        raise ValueError("Not a generic reaction index artifact")
    if payload.get("schema_version") != GENERIC_INDEX_SCHEMA_VERSION:
        raise ValueError("Unsupported generic reaction index schema; rebuild the index")
    if (
        payload.get("reaction_signature_schema_version")
        != REACTION_SIGNATURE_SCHEMA_VERSION
    ):
        raise ValueError("Incompatible reaction signature schema; rebuild the index")
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
    if not record_schema_versions or not set(record_schema_versions).issubset(
        _INDEX_COMPATIBLE_RECORD_SCHEMAS
    ):
        raise ValueError("Incompatible recommendation record schema; rebuild the index")
    converter_definition_versions = tuple(
        payload.get("converter_definition_versions") or ()
    )
    if not converter_definition_versions or not set(
        converter_definition_versions
    ).issubset(_INDEX_COMPATIBLE_CONVERTER_VERSIONS):
        raise ValueError("Incompatible generic converter version; rebuild the index")
    rows = tuple(
        GenericIndexedReaction(
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
            reaction_label=str(row.get("reaction_label") or ""),
            reaction_label_status=str(
                row.get("reaction_label_status") or "unavailable"
            ),
            fallback_descriptor=dict(row.get("fallback_descriptor") or {}),
        )
        for row in payload.get("rows") or ()
    )
    maps = payload.get("maps") or {}
    index = GenericReactionIndex(
        rows=rows,
        exact={key: tuple(value) for key, value in (maps.get("exact") or {}).items()},
        handles={
            key: tuple(value) for key, value in (maps.get("handles") or {}).items()
        },
        transformations={
            key: tuple(value)
            for key, value in (maps.get("transformations") or {}).items()
        },
        bond_edits={
            key: tuple(value) for key, value in (maps.get("bond_edits") or {}).items()
        },
        environments={
            key: tuple(value) for key, value in (maps.get("environments") or {}).items()
        },
        environment_features={
            key: tuple(value)
            for key, value in (maps.get("environment_features") or {}).items()
        },
        fallback_features={
            key: tuple(value)
            for key, value in (maps.get("fallback_features") or {}).items()
        },
        families={
            key: tuple(value) for key, value in (maps.get("families") or {}).items()
        },
        reaction_signature_schema_version=str(
            payload["reaction_signature_schema_version"]
        ),
        taxonomy_definition_versions=tuple(
            sorted(
                (str(key), str(value))
                for key, value in payload["taxonomy_definition_versions"].items()
            )
        ),
        fallback_descriptor_schema_version=str(
            payload["fallback_descriptor_schema_version"]
        ),
        fallback_definition_versions=tuple(
            sorted(
                (str(key), str(value))
                for key, value in payload["fallback_definition_versions"].items()
            )
        ),
        record_schema_versions=tuple(payload["record_schema_versions"]),
        converter_definition_versions=tuple(payload["converter_definition_versions"]),
    )
    expected = _index_payload(index)
    if expected["index_id"] != payload.get("index_id"):
        raise ValueError("Generic reaction index integrity check failed")
    if int(payload.get("row_count", -1)) != len(rows):
        raise ValueError("Generic reaction index row count mismatch")
    return index


def load_generic_index(
    path: str | Path,
    *,
    include_review: bool = False,
) -> GenericReactionIndex:
    """Load canonical JSONL output from the generic conversion engine."""
    source = Path(path)
    if source.name.casefold().endswith(".json.gz"):
        return load_persisted_generic_index(source)
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
        return load_persisted_generic_index(source)
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
    """Validate persisted map coverage, row identity, and recipe contracts."""
    source = Path(path)
    index = load_persisted_generic_index(source)
    issues = []
    row_count = len(index.rows)
    observation_ids = [row.observation_id for row in index.rows]
    duplicate_observation_count = len(observation_ids) - len(set(observation_ids))
    if duplicate_observation_count:
        issues.append("duplicate_observation_ids")
    maps = {
        "exact": (index.exact, "exact_signature_key"),
        "handles": (index.handles, "handle_signature_key"),
        "transformations": (
            index.transformations,
            "transformation_signature_key",
        ),
        "bond_edits": (index.bond_edits, "bond_edit_signature_key"),
        "environments": (index.environments, "environment_signature_key"),
    }
    for name, (mapping, field) in maps.items():
        for positions in mapping.values():
            if any(position < 0 or position >= row_count for position in positions):
                issues.append(f"out_of_range_position:{name}")
        for position, row in enumerate(index.rows):
            key = str(row.signature.get(field) or "")
            if key and position not in mapping.get(key, ()):
                issues.append(f"missing_reverse_mapping:{name}")
    for position, row in enumerate(index.rows):
        for token in set(environment_tokens(row.signature)):
            if position not in index.environment_features.get(token, ()):
                issues.append("missing_environment_feature_mapping")
        for token in fallback_index_tokens(row.fallback_descriptor):
            if position not in index.fallback_features.get(token, ()):
                issues.append("missing_fallback_feature_mapping")
        if row.named_family and position not in index.families.get(
            row.named_family, ()
        ):
            issues.append("missing_family_mapping")
        if row.resolved_recipe.get("recipe_id") != row.recipe_id:
            issues.append("recipe_identity_mismatch")
        if (
            row.resolved_recipe.get("recipe_core_id") or row.recipe_core_id
        ) != row.recipe_core_id:
            issues.append("recipe_core_identity_mismatch")
    opener = gzip.open if source.suffix.casefold() == ".gz" else Path.open
    arguments = (
        {"mode": "rt", "encoding": "utf-8"}
        if source.suffix.casefold() == ".gz"
        else {"mode": "r", "encoding": "utf-8"}
    )
    with opener(source, **arguments) as handle:
        payload = json.load(handle)
    report = {
        "schema_version": "1.0",
        "artifact_type": "generic_index_integrity",
        "path": str(source),
        "valid": not issues,
        "issues": sorted(set(issues)),
        "index_id": payload["index_id"],
        "index_schema_version": payload["schema_version"],
        "row_count": row_count,
        "file_size_bytes": source.stat().st_size,
        "duplicate_observation_count": duplicate_observation_count,
        "key_counts": {
            "exact": len(index.exact),
            "handles": len(index.handles),
            "transformations": len(index.transformations),
            "bond_edits": len(index.bond_edits),
            "environments": len(index.environments),
            "environment_features": len(index.environment_features),
            "fallback_features": len(index.fallback_features),
            "families": len(index.families),
        },
        "transformation_class_counts": dict(
            sorted(Counter(row.transformation_class for row in index.rows).items())
        ),
        "named_family_counts": dict(
            sorted(Counter(row.named_family or "unnamed" for row in index.rows).items())
        ),
        "reaction_scope_counts": dict(
            sorted(
                Counter(
                    str(
                        (row.signature.get("topology") or {}).get("reaction_scope")
                        or "unknown"
                    )
                    for row in index.rows
                ).items()
            )
        ),
        "unique_reference_count": len(
            {row.reference_id for row in index.rows if row.reference_id}
        ),
        "unique_recipe_core_count": len({row.recipe_core_id for row in index.rows}),
    }
    return report


__all__ = [
    "GenericIndexedReaction",
    "GenericReactionIndex",
    "build_generic_index",
    "build_generic_index_from_rows",
    "load_generic_index",
    "load_persisted_generic_index",
    "save_generic_index",
    "validate_generic_index_artifact",
]
