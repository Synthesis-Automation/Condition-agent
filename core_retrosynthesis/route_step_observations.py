"""Extract patent-split one-step observations from canonical multistep routes."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from condition_recommender.ingestion.adapters.base import (
    observation_id,
    supplied_mapping_status,
)
from condition_recommender.ingestion.models import (
    INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
    CanonicalSourceObservation,
    ConditionComponentClaim,
    ConditionInput,
    ConditionStageInput,
    ReactionEvidenceInput,
    SourceIdentifier,
    SourceProvenance,
)

from .chemistry import canonical_smiles
from .generic_compiler import generic_operator_identity_from_observation
from .route_contract import ReactionRouteTree, iter_route_reactions
from .route_conversion import iter_route_trees
from .route_core import RouteCoreProjection, RouteCoreStep
from .route_core_conversion import iter_route_core_projections


ROUTE_STEP_OBSERVATION_SCHEMA_VERSION = "1.0"
ROUTE_STEP_OBSERVATION_ALGORITHM_VERSION = "route_step_observation.v1"
ROUTE_STEP_OBSERVATION_ADAPTER_ID = "canonical_route_step.v1"
ROUTE_STEP_OBSERVATION_ADAPTER_VERSION = "1.0"
ROUTE_STEP_OBSERVATION_CORPUS_ID = "canonical_multistep_route_steps"
ROUTE_STEP_SPLITS = ("train", "validation", "test")


@dataclass(frozen=True)
class _RouteMembership:
    tree_id: str
    route_core_id: str
    source_route_id: str
    reaction_node_id: str
    step_id: str
    retrosynthetic_depth: int
    observed_remaining_steps: int
    product_occurrence_id: str
    precursor_occurrence_ids: tuple[str, ...]
    internal_precursor_occurrence_ids: tuple[str, ...]
    terminal_precursor_occurrence_ids: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "tree_id": self.tree_id,
            "route_core_id": self.route_core_id,
            "source_route_id": self.source_route_id,
            "reaction_node_id": self.reaction_node_id,
            "step_id": self.step_id,
            "retrosynthetic_depth": self.retrosynthetic_depth,
            "observed_remaining_steps": self.observed_remaining_steps,
            "product_occurrence_id": self.product_occurrence_id,
            "precursor_occurrence_ids": list(self.precursor_occurrence_ids),
            "internal_precursor_occurrence_ids": list(
                self.internal_precursor_occurrence_ids
            ),
            "terminal_precursor_occurrence_ids": list(
                self.terminal_precursor_occurrence_ids
            ),
        }


@dataclass
class _AccumulatedStep:
    split: str
    patent_id: str
    source_reaction_id: str
    reaction_smiles: str
    original_reaction_smiles: str
    reagents_smiles: str
    tree: ReactionRouteTree
    projection: RouteCoreProjection
    step: RouteCoreStep
    memberships: list[_RouteMembership]


def _sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def _two_part_reaction(reaction_smiles: str) -> tuple[str, str, str] | None:
    """Return mapped transformation and the source middle field separately."""

    if reaction_smiles.count(">>") == 1:
        reactants, products = reaction_smiles.split(">>", 1)
        reagents = ""
    elif reaction_smiles.count(">") == 2:
        reactants, reagents, products = reaction_smiles.split(">", 2)
    else:
        return None
    if not reactants or not products:
        return None
    return f"{reactants}>>{products}", reagents, products


def _condition_input(reagents_smiles: str) -> ConditionInput:
    """Preserve middle-field structures without inventing chemical roles."""

    components = tuple(
        ConditionComponentClaim(
            component_key=f"route_agent_{index}",
            source_slot=f"reaction_middle_component_{index}",
            source_role_hint=None,
            identifiers=(
                SourceIdentifier(
                    identifier_type="smiles",
                    value=component,
                    source_field="reaction_middle_field",
                ),
            ),
            introduced_in_stage=1,
            provenance={
                "source_declared_role": "unresolved_reaction_middle_field",
                "source_position": index,
            },
        )
        for index, component in enumerate(
            (item for item in reagents_smiles.split(".") if item),
            start=1,
        )
    )
    stages = (
        (
            ConditionStageInput(
                stage_index=1,
                component_keys=tuple(item.component_key for item in components),
                provenance={"source": "route_reaction_middle_field"},
            ),
        )
        if components
        else ()
    )
    warnings = (
        ("ROUTE_MIDDLE_FIELD_ROLES_UNRESOLVED",)
        if components
        else ("ROUTE_CONDITION_DETAILS_UNAVAILABLE",)
    )
    return ConditionInput(
        components=components,
        stages=stages,
        warnings=warnings,
    )


def _membership(
    tree: ReactionRouteTree,
    projection: RouteCoreProjection,
    step: RouteCoreStep,
) -> _RouteMembership:
    return _RouteMembership(
        tree_id=tree.tree_id,
        route_core_id=projection.route_core_id,
        source_route_id=str(tree.source_route_id or ""),
        reaction_node_id=step.reaction_node_id,
        step_id=step.step_id,
        retrosynthetic_depth=step.retrosynthetic_depth,
        observed_remaining_steps=step.observed_remaining_steps,
        product_occurrence_id=step.product_occurrence_id,
        precursor_occurrence_ids=step.precursor_occurrence_ids,
        internal_precursor_occurrence_ids=step.internal_precursor_occurrence_ids,
        terminal_precursor_occurrence_ids=step.terminal_precursor_occurrence_ids,
    )


def _validate_pair(
    tree: ReactionRouteTree,
    projection: RouteCoreProjection,
) -> None:
    if projection.source_tree_id != tree.tree_id:
        raise ValueError("route tree and route-core identities do not match")
    if projection.reaction_count != tree.reaction_count:
        raise ValueError("route tree and route-core reaction counts differ")
    if projection.patent_id != tree.patent_id:
        raise ValueError("route tree and route-core patent IDs differ")
    if projection.split != tree.split:
        raise ValueError("route tree and route-core splits differ")


def _accumulate_steps(
    route_tree_source: Path,
    route_core_source: Path,
    *,
    minimum_reaction_count: int,
) -> tuple[dict[tuple[str, str], _AccumulatedStep], dict[str, Any]]:
    projections: dict[str, RouteCoreProjection] = {}
    for projection in iter_route_core_projections(route_core_source):
        if projection.source_tree_id in projections:
            raise ValueError("duplicate route-core source tree ID")
        projections[projection.source_tree_id] = projection

    accumulated: dict[tuple[str, str], _AccumulatedStep] = {}
    source_splits: dict[str, str] = {}
    counts: Counter[str] = Counter()
    used_projection_ids: set[str] = set()
    for tree in iter_route_trees(route_tree_source):
        counts["route_count"] += 1
        if tree.reaction_count < minimum_reaction_count:
            counts["excluded_short_route_count"] += 1
            continue
        counts["multistep_route_count"] += 1
        projection = projections.get(tree.tree_id)
        if projection is None:
            raise ValueError(f"route core is missing for tree {tree.tree_id}")
        _validate_pair(tree, projection)
        used_projection_ids.add(projection.route_core_id)
        split = str(tree.split or "")
        patent_id = str(tree.patent_id or "")
        if split not in ROUTE_STEP_SPLITS:
            raise ValueError(f"route tree has unsupported split: {split}")
        if not patent_id:
            raise ValueError("route tree requires a patent ID before extraction")
        core_steps = {item.reaction_node_id: item for item in projection.steps}
        for reaction in iter_route_reactions(tree):
            counts["step_occurrence_count"] += 1
            step = core_steps.get(reaction.reaction_node_id)
            if step is None:
                raise ValueError("route reaction is missing from route core")
            parsed = _two_part_reaction(reaction.reaction_smiles)
            if parsed is None:
                raise ValueError("route step has invalid reaction SMILES")
            reaction_smiles, middle_field, _ = parsed
            source_reaction_id = str(
                reaction.evidence.source_reaction_id
                or step.source_reaction_id
                or reaction.reaction_node_id
            )
            prior_split = source_splits.setdefault(source_reaction_id, split)
            if prior_split != split:
                raise ValueError(
                    "source reaction occurs across patent data splits"
                )
            key = (split, source_reaction_id)
            route_membership = _membership(tree, projection, step)
            current = accumulated.get(key)
            if current is None:
                accumulated[key] = _AccumulatedStep(
                    split=split,
                    patent_id=patent_id,
                    source_reaction_id=source_reaction_id,
                    reaction_smiles=reaction_smiles,
                    original_reaction_smiles=reaction.reaction_smiles,
                    reagents_smiles=(
                        reaction.evidence.reagents_smiles or middle_field
                    ),
                    tree=tree,
                    projection=projection,
                    step=step,
                    memberships=[route_membership],
                )
                continue
            if (
                current.reaction_smiles != reaction_smiles
                or current.patent_id != patent_id
            ):
                raise ValueError(
                    "source reaction ID contradicts structural or patent evidence"
                )
            current.memberships.append(route_membership)
            counts["duplicate_source_occurrence_count"] += 1

    unused = {
        item.route_core_id for item in projections.values()
    }.difference(used_projection_ids)
    if unused and counts["excluded_short_route_count"] == 0:
        raise ValueError("route-core source contains unmatched projections")
    counts["unique_step_count"] = len(accumulated)
    return accumulated, dict(counts)


def _observation(
    item: _AccumulatedStep,
    *,
    tree_source: Path,
    tree_source_sha256: str,
    core_source_sha256: str,
    combined_source_sha256: str,
    row_number: int,
) -> CanonicalSourceObservation:
    step = item.step
    reactants, products = item.reaction_smiles.split(">>", 1)
    canonical_reactants = canonical_smiles(reactants)
    canonical_products = canonical_smiles(products)
    canonical_reaction = (
        f"{canonical_reactants}>>{canonical_products}"
        if canonical_reactants and canonical_products
        else ""
    )
    operator = generic_operator_identity_from_observation(
        step.reaction_smiles,
        step.reaction_signature or {},
    )
    warnings = {
        "ROUTE_CONNECTIVITY_INFERRED_FROM_MOLECULAR_IDENTITY",
        *step.warnings,
    }
    if not step.chemistry_valid or step.reaction_signature is None:
        warnings.add("ROUTE_CORE_NOT_COMPILER_READY")
    record_id = f"{item.patent_id}:{item.source_reaction_id}"
    memberships = sorted(
        item.memberships,
        key=lambda value: (
            value.source_route_id,
            value.tree_id,
            value.reaction_node_id,
        ),
    )
    return CanonicalSourceObservation(
        observation_id=observation_id(
            adapter_id=ROUTE_STEP_OBSERVATION_ADAPTER_ID,
            source_sha256=combined_source_sha256,
            row_number=row_number,
            record_id=record_id,
        ),
        observation_kind="structure_backed",
        source=SourceProvenance(
            corpus_id=ROUTE_STEP_OBSERVATION_CORPUS_ID,
            release_id=tree_source.name,
            adapter_id=ROUTE_STEP_OBSERVATION_ADAPTER_ID,
            adapter_version=ROUTE_STEP_OBSERVATION_ADAPTER_VERSION,
            source_file=tree_source.name,
            source_file_sha256=tree_source_sha256,
            source_row_number=row_number,
            source_record_id=record_id,
            source_groups={
                "patent_id": item.patent_id,
                "split": item.split,
                "source_reaction_id": item.source_reaction_id,
                "source_route_id": str(item.tree.source_route_id or ""),
                "source_tree_id": item.tree.tree_id,
                "route_core_id": item.projection.route_core_id,
                "route_core_sha256": core_source_sha256,
            },
            reference=item.patent_id,
        ),
        reaction=ReactionEvidenceInput(
            evidence_kind="source_structure",
            reaction_smiles=item.reaction_smiles,
            supplied_mapping_status=supplied_mapping_status(
                item.reaction_smiles
            ),
            source_labels={
                "canonical_reaction_smiles": canonical_reaction,
                "reaction_smiles_source_field": "route_reactants_and_product",
                "route_core_quality_status": step.quality_status,
                "route_core_evidence_quality": step.evidence_quality,
                "reaction_signature_id": step.reaction_signature_id,
                "reaction_core_id": step.reaction_core_id,
                "transformation_class": step.transformation_class,
                "named_family": step.named_family,
                "route_operator_id": (
                    operator.operator_id if operator is not None else None
                ),
            },
            structure_available=True,
        ),
        conditions=_condition_input(item.reagents_smiles),
        ingestion_status=(
            "accepted"
            if step.chemistry_valid and step.reaction_signature is not None
            else "review"
        ),
        warnings=tuple(sorted(warnings)),
        raw_fields={
            "original_reaction_smiles": item.original_reaction_smiles,
            "reagents_smiles": item.reagents_smiles,
            "source_reaction_id": item.source_reaction_id,
            "route_memberships": [value.to_dict() for value in memberships],
            "route_core_algorithm_version": item.projection.algorithm_version,
            "route_core_schema_version": item.projection.schema_version,
        },
    )


def _write_observations(
    destination: Path,
    observations: Sequence[CanonicalSourceObservation],
) -> None:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
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
                newline="\n",
            ) as stream:
                for observation in observations:
                    stream.write(
                        json.dumps(
                            observation.to_dict(),
                            ensure_ascii=False,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                    )
                    stream.write("\n")
    temporary.replace(destination)


def extract_route_step_observations(
    route_tree_source: str | Path,
    route_core_source: str | Path,
    output_directory: str | Path,
    *,
    minimum_reaction_count: int = 2,
) -> dict[str, Any]:
    """Write deterministic train/validation/test physical-step observations."""

    if minimum_reaction_count < 2:
        raise ValueError("route-step extraction requires multistep routes")
    tree_source = Path(route_tree_source)
    core_source = Path(route_core_source)
    if not tree_source.is_file() or not core_source.is_file():
        raise ValueError("route tree and route-core sources must exist")
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    tree_sha256 = _sha256(tree_source)
    core_sha256 = _sha256(core_source)
    combined_sha256 = hashlib.sha256(
        f"{tree_sha256}\0{core_sha256}".encode("utf-8")
    ).hexdigest()
    accumulated, counts = _accumulate_steps(
        tree_source,
        core_source,
        minimum_reaction_count=minimum_reaction_count,
    )
    grouped: dict[str, list[_AccumulatedStep]] = {
        split: [] for split in ROUTE_STEP_SPLITS
    }
    for item in accumulated.values():
        grouped[item.split].append(item)

    artifacts = {}
    quality_counts: Counter[str] = Counter()
    for split in ROUTE_STEP_SPLITS:
        items = sorted(
            grouped[split],
            key=lambda item: (
                item.patent_id,
                item.source_reaction_id,
                item.reaction_smiles,
            ),
        )
        observations = tuple(
            _observation(
                item,
                tree_source=tree_source,
                tree_source_sha256=tree_sha256,
                core_source_sha256=core_sha256,
                combined_source_sha256=combined_sha256,
                row_number=row_number,
            )
            for row_number, item in enumerate(items, start=1)
        )
        for item in items:
            quality_counts[item.step.quality_status] += 1
        destination = output / f"route_steps.{split}.observations.jsonl.gz"
        _write_observations(destination, observations)
        artifacts[split] = {
            "path": str(destination.resolve()),
            "row_count": len(observations),
            "sha256": _sha256(destination),
            "size_bytes": destination.stat().st_size,
            "patent_count": len({item.patent_id for item in items}),
        }
    report = {
        "artifact_type": "route_step_observation_extraction",
        "schema_version": ROUTE_STEP_OBSERVATION_SCHEMA_VERSION,
        "algorithm_version": ROUTE_STEP_OBSERVATION_ALGORITHM_VERSION,
        "intermediate_schema_version": INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
        "route_tree_source": str(tree_source.resolve()),
        "route_tree_sha256": tree_sha256,
        "route_core_source": str(core_source.resolve()),
        "route_core_sha256": core_sha256,
        "combined_source_sha256": combined_sha256,
        "minimum_reaction_count": minimum_reaction_count,
        **counts,
        "quality_counts": dict(sorted(quality_counts.items())),
        "artifacts": artifacts,
        "warnings": [
            "Route connectivity is inferred from patent-local molecular identity.",
            "Reaction middle-field components retain unresolved condition roles.",
            "Only training-split observations may be used to compile operators.",
        ],
    }
    manifest = output / "route_step_observation_manifest.json"
    manifest.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return {**report, "manifest_path": str(manifest.resolve())}


__all__ = [
    "ROUTE_STEP_OBSERVATION_ADAPTER_ID",
    "ROUTE_STEP_OBSERVATION_ADAPTER_VERSION",
    "ROUTE_STEP_OBSERVATION_ALGORITHM_VERSION",
    "ROUTE_STEP_OBSERVATION_CORPUS_ID",
    "ROUTE_STEP_OBSERVATION_SCHEMA_VERSION",
    "ROUTE_STEP_SPLITS",
    "extract_route_step_observations",
]
