"""Mine causally coupled two-event strategies from observed route cores."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

from rdkit import Chem

from .chemistry import canonical_smiles, digest
from .route_core import (
    RouteAtomLineageCandidate,
    RouteCoreLineageLink,
    RouteCoreProjection,
    RouteCoreStep,
)
from .route_core_conversion import iter_route_core_projections


COUPLED_ROUTE_STRATEGY_SCHEMA_VERSION = "1.1"
COUPLED_ROUTE_STRATEGY_ALGORITHM_VERSION = "coupled_route_strategy.v2"
COUPLED_ROUTE_STRATEGY_REVIEW_VERSION = "1.1"
RELATIONSHIP_CLASSES = frozenset(
    {
        "handle_progression",
        "same_site_coupled",
        "shared_local_environment",
        "independent_sites",
        "lineage_ambiguous",
        "unresolved",
    }
)
ADMISSION_CLASSES = frozenset({"strict", "review", "rejected"})
DEPENDENCY_CLASSES = frozenset(
    {
        "created_handle_consumed",
        "activation_then_conversion",
        "temporary_group_removed",
        "continued_site_transformation",
        "shared_local_environment",
        "independent_sites",
        "lineage_ambiguous",
        "unresolved",
    }
)
_STRICT_DEPENDENCIES = frozenset(
    {
        "created_handle_consumed",
        "activation_then_conversion",
        "continued_site_transformation",
    }
)


@dataclass(frozen=True)
class CouplingEvidence:
    """Inspectable cross-step site-continuity evidence."""

    producer_active_maps: tuple[int, ...]
    producer_persistent_active_maps: tuple[int, ...]
    consumer_active_maps: tuple[int, ...]
    consumer_intermediate_active_maps: tuple[int, ...]
    overlap_counts: tuple[int, ...]
    formed_overlap_counts: tuple[int, ...]
    transient_bond_counts: tuple[int, ...]
    replacement_bond_counts: tuple[int, ...]
    minimum_distances: tuple[Optional[int], ...]
    ambiguity_invariant: bool
    candidate_limit_reached: bool
    rationale: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible evidence."""

        return {
            "producer_active_maps": list(self.producer_active_maps),
            "producer_persistent_active_maps": list(
                self.producer_persistent_active_maps
            ),
            "consumer_active_maps": list(self.consumer_active_maps),
            "consumer_intermediate_active_maps": list(
                self.consumer_intermediate_active_maps
            ),
            "overlap_counts": list(self.overlap_counts),
            "formed_overlap_counts": list(self.formed_overlap_counts),
            "transient_bond_counts": list(self.transient_bond_counts),
            "replacement_bond_counts": list(self.replacement_bond_counts),
            "minimum_distances": list(self.minimum_distances),
            "ambiguity_invariant": self.ambiguity_invariant,
            "candidate_limit_reached": self.candidate_limit_reached,
            "rationale": list(self.rationale),
        }


@dataclass(frozen=True)
class CoupledRouteStrategyOccurrence:
    """One structurally classified adjacent pair from an observed route."""

    occurrence_id: str
    exact_strategy_id: str
    typed_strategy_id: str
    shape_strategy_id: str
    overall_reaction_id: str
    route_core_id: str
    source_tree_id: str
    source_route_id: Optional[str]
    patent_id: Optional[str]
    split: Optional[str]
    lineage_id: str
    lineage_status: str
    relationship_class: str
    dependency_class: str
    admission_class: str
    coupling_score: float
    first_reaction_node_id: str
    second_reaction_node_id: str
    first_source_reaction_id: Optional[str]
    second_source_reaction_id: Optional[str]
    first_reaction_smiles: str
    second_reaction_smiles: str
    intermediate_smiles: str
    final_product_smiles: str
    overall_reaction_smiles: str
    first_transformation_class: Optional[str]
    second_transformation_class: Optional[str]
    first_named_family: Optional[str]
    second_named_family: Optional[str]
    first_edit_tokens: tuple[str, ...]
    second_edit_tokens: tuple[str, ...]
    combined_edit_tokens: tuple[str, ...]
    evidence: CouplingEvidence
    warnings: tuple[str, ...] = ()
    schema_version: str = COUPLED_ROUTE_STRATEGY_SCHEMA_VERSION
    algorithm_version: str = COUPLED_ROUTE_STRATEGY_ALGORITHM_VERSION

    def __post_init__(self) -> None:
        if self.relationship_class not in RELATIONSHIP_CLASSES:
            raise ValueError("unsupported coupled-strategy relationship")
        if self.dependency_class not in DEPENDENCY_CLASSES:
            raise ValueError("unsupported coupled-strategy dependency")
        if self.admission_class not in ADMISSION_CLASSES:
            raise ValueError("unsupported coupled-strategy admission class")
        if not 0.0 <= self.coupling_score <= 1.0:
            raise ValueError("coupling score must be between zero and one")
        if self.dependency_class in _STRICT_DEPENDENCIES:
            if self.admission_class != "strict":
                raise ValueError("strict dependency requires strict admission")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible strategy occurrence."""

        return {
            "occurrence_id": self.occurrence_id,
            "exact_strategy_id": self.exact_strategy_id,
            "typed_strategy_id": self.typed_strategy_id,
            "shape_strategy_id": self.shape_strategy_id,
            "overall_reaction_id": self.overall_reaction_id,
            "route_core_id": self.route_core_id,
            "source_tree_id": self.source_tree_id,
            "source_route_id": self.source_route_id,
            "patent_id": self.patent_id,
            "split": self.split,
            "lineage_id": self.lineage_id,
            "lineage_status": self.lineage_status,
            "relationship_class": self.relationship_class,
            "dependency_class": self.dependency_class,
            "admission_class": self.admission_class,
            "coupling_score": self.coupling_score,
            "first_reaction_node_id": self.first_reaction_node_id,
            "second_reaction_node_id": self.second_reaction_node_id,
            "first_source_reaction_id": self.first_source_reaction_id,
            "second_source_reaction_id": self.second_source_reaction_id,
            "first_reaction_smiles": self.first_reaction_smiles,
            "second_reaction_smiles": self.second_reaction_smiles,
            "intermediate_smiles": self.intermediate_smiles,
            "final_product_smiles": self.final_product_smiles,
            "overall_reaction_smiles": self.overall_reaction_smiles,
            "first_transformation_class": self.first_transformation_class,
            "second_transformation_class": self.second_transformation_class,
            "first_named_family": self.first_named_family,
            "second_named_family": self.second_named_family,
            "first_edit_tokens": list(self.first_edit_tokens),
            "second_edit_tokens": list(self.second_edit_tokens),
            "combined_edit_tokens": list(self.combined_edit_tokens),
            "evidence": self.evidence.to_dict(),
            "warnings": list(self.warnings),
            "schema_version": self.schema_version,
            "algorithm_version": self.algorithm_version,
        }


def _reaction_sides(reaction_smiles: str) -> Optional[tuple[str, str]]:
    if reaction_smiles.count(">>") == 1:
        reactants, products = reaction_smiles.split(">>", 1)
    elif reaction_smiles.count(">") == 2:
        reactants, _, products = reaction_smiles.split(">", 2)
    else:
        return None
    return (reactants, products) if reactants and products else None


def _active_maps(step: RouteCoreStep) -> set[int]:
    signature = step.reaction_signature or {}
    values: set[int] = set()
    for edit in signature.get("edits") or ():
        for field in ("atom_1", "atom_2"):
            atom = edit.get(field)
            if isinstance(atom, Mapping) and atom.get("atom_map_number"):
                values.add(int(atom["atom_map_number"]))
    return values


def _formed_maps(step: RouteCoreStep) -> set[int]:
    signature = step.reaction_signature or {}
    values: set[int] = set()
    for edit in signature.get("edits") or ():
        if edit.get("edit_type") != "formed":
            continue
        for field in ("atom_1", "atom_2"):
            atom = edit.get(field)
            if isinstance(atom, Mapping) and atom.get("atom_map_number"):
                values.add(int(atom["atom_map_number"]))
    return values


def _component_map_distances(
    step: RouteCoreStep,
    component_index: Optional[int],
    candidates: Sequence[RouteAtomLineageCandidate],
    producer_active: set[int],
    consumer_active: set[int],
) -> tuple[tuple[int, ...], tuple[Optional[int], ...], set[int]]:
    sides = _reaction_sides(step.reaction_smiles)
    if sides is None or component_index is None:
        count = len(candidates)
        return ((0,) * count, (None,) * count, set())
    components = sides[0].split(".")
    if component_index < 0 or component_index >= len(components):
        count = len(candidates)
        return ((0,) * count, (None,) * count, set())
    molecule = Chem.MolFromSmiles(components[component_index])
    if molecule is None:
        count = len(candidates)
        return ((0,) * count, (None,) * count, set())
    map_to_index = {
        int(atom.GetAtomMapNum()): int(atom.GetIdx())
        for atom in molecule.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    intermediate_active = consumer_active.intersection(map_to_index)
    overlaps = []
    distances: list[Optional[int]] = []
    for candidate in candidates:
        translated = {
            consumer_map
            for producer_map, consumer_map in candidate.atom_map_pairs
            if producer_map in producer_active and consumer_map in map_to_index
        }
        overlaps.append(len(translated.intersection(intermediate_active)))
        path_lengths = []
        for first_map in translated:
            for second_map in intermediate_active:
                if first_map == second_map:
                    path_lengths.append(0)
                    continue
                path = Chem.GetShortestPath(
                    molecule,
                    map_to_index[first_map],
                    map_to_index[second_map],
                )
                if path:
                    path_lengths.append(len(path) - 1)
        distances.append(min(path_lengths) if path_lengths else None)
    return tuple(overlaps), tuple(distances), intermediate_active


def _formed_overlap_counts(
    candidates: Sequence[RouteAtomLineageCandidate],
    producer_formed: set[int],
    consumer_intermediate_active: set[int],
) -> tuple[int, ...]:
    return tuple(
        len(
            {
                consumer_map
                for producer_map, consumer_map in candidate.atom_map_pairs
                if producer_map in producer_formed
            }.intersection(consumer_intermediate_active)
        )
        for candidate in candidates
    )


def _edit_bonds(
    step: RouteCoreStep, edit_types: set[str]
) -> set[tuple[int, int]]:
    """Return normalized mapped bonds for selected signature edit types."""

    bonds: set[tuple[int, int]] = set()
    signature = step.reaction_signature or {}
    for edit in signature.get("edits") or ():
        if str(edit.get("edit_type") or "") not in edit_types:
            continue
        first = edit.get("atom_1")
        second = edit.get("atom_2")
        if not isinstance(first, Mapping) or not isinstance(second, Mapping):
            continue
        first_map = int(first.get("atom_map_number") or 0)
        second_map = int(second.get("atom_map_number") or 0)
        if first_map > 0 and second_map > 0 and first_map != second_map:
            bonds.add(tuple(sorted((first_map, second_map))))
    return bonds


def _bond_dependency_counts(
    candidates: Sequence[RouteAtomLineageCandidate],
    producer_formed_bonds: set[tuple[int, int]],
    consumer_broken_bonds: set[tuple[int, int]],
    consumer_constructive_bonds: set[tuple[int, int]],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Measure bonds installed by step 1 and consumed or replaced in step 2."""

    transient_counts: list[int] = []
    replacement_counts: list[int] = []
    for candidate in candidates:
        correspondence = dict(candidate.atom_map_pairs)
        translated = {
            tuple(sorted((correspondence[first], correspondence[second])))
            for first, second in producer_formed_bonds
            if first in correspondence and second in correspondence
        }
        transient = translated.intersection(consumer_broken_bonds)
        replacements = {
            constructive
            for constructive in consumer_constructive_bonds
            if any(
                constructive != removed
                and bool(set(constructive).intersection(removed))
                for removed in transient
            )
        }
        transient_counts.append(len(transient))
        replacement_counts.append(len(replacements))
    return tuple(transient_counts), tuple(replacement_counts)


def _classify(
    link: RouteCoreLineageLink,
    first: RouteCoreStep,
    second: RouteCoreStep,
) -> tuple[str, str, str, float, CouplingEvidence, tuple[str, ...]]:
    producer_active = _active_maps(first)
    consumer_active = _active_maps(second)
    producer_maps = {
        producer_map
        for candidate in link.candidates
        for producer_map, _ in candidate.atom_map_pairs
    }
    persistent_active = producer_active.intersection(producer_maps)
    overlaps, distances, intermediate_active = _component_map_distances(
        second,
        link.consumer_reactant_component_index,
        link.candidates,
        persistent_active,
        consumer_active,
    )
    formed_overlaps = _formed_overlap_counts(
        link.candidates,
        _formed_maps(first).intersection(producer_maps),
        intermediate_active,
    )
    transient_counts, replacement_counts = _bond_dependency_counts(
        link.candidates,
        _edit_bonds(first, {"formed"}),
        _edit_bonds(second, {"broken"}),
        _edit_bonds(second, {"formed", "order_changed"}),
    )
    warnings: list[str] = []
    rationale: list[str] = []
    complete_candidates = (
        bool(link.candidates)
        and not link.candidate_limit_reached
        and link.candidate_count == len(link.candidates)
    )
    chemistry_resolved = all(
        step.chemistry_valid
        and step.reaction_signature is not None
        and step.quality_status not in {"blocked", "unavailable"}
        for step in (first, second)
    )
    if not chemistry_resolved:
        relationship = "unresolved"
        dependency, admission, score = "unresolved", "review", 0.20
        rationale.append("one or both step cores lack admissible chemistry evidence")
    elif not complete_candidates:
        relationship = "lineage_ambiguous"
        dependency, admission, score = "lineage_ambiguous", "review", 0.25
        rationale.append("cross-step atom lineage is incomplete or bounded")
    elif not persistent_active or not intermediate_active:
        relationship = "independent_sites"
        dependency, admission, score = "independent_sites", "rejected", 0.0
        rationale.append("one step does not edit the carried intermediate site")
    elif overlaps and min(overlaps) > 0:
        if formed_overlaps and min(formed_overlaps) > 0:
            relationship = "handle_progression"
            if transient_counts and min(transient_counts) > 0:
                if replacement_counts and min(replacement_counts) > 0:
                    dependency, admission, score = (
                        "activation_then_conversion",
                        "strict",
                        1.0,
                    )
                    rationale.append(
                        "step 2 breaks a step-1-installed bond and forms or "
                        "changes another bond at that handle"
                    )
                elif replacement_counts and max(replacement_counts) > 0:
                    dependency, admission, score = (
                        "lineage_ambiguous",
                        "review",
                        0.40,
                    )
                    rationale.append(
                        "whether the consumed handle is replaced depends on "
                        "the atom-lineage candidate"
                    )
                else:
                    dependency, admission, score = (
                        "temporary_group_removed",
                        "review",
                        0.65,
                    )
                    rationale.append(
                        "step 2 removes a step-1-installed bond without a "
                        "mapped heavy-atom replacement; strategic value needs review"
                    )
            elif transient_counts and max(transient_counts) > 0:
                dependency, admission, score = (
                    "lineage_ambiguous",
                    "review",
                    0.40,
                )
                rationale.append(
                    "whether step 2 breaks the step-1-installed bond depends "
                    "on the atom-lineage candidate"
                )
            else:
                dependency, admission, score = (
                    "created_handle_consumed",
                    "strict",
                    1.0,
                )
                rationale.append(
                    "step 2 transforms an atom in a bond installed by step 1"
                )
        elif formed_overlaps and max(formed_overlaps) > 0:
            relationship = "lineage_ambiguous"
            dependency, admission, score = (
                "lineage_ambiguous",
                "review",
                0.35,
            )
            rationale.append(
                "whether step 2 consumes a step-1-created handle depends on "
                "the atom-lineage candidate"
            )
        else:
            relationship = "same_site_coupled"
            dependency, admission, score = (
                "continued_site_transformation",
                "strict",
                0.90,
            )
            rationale.append(
                "both steps transform the same lineage-traced site without "
                "consuming a newly installed bond"
            )
    elif overlaps and max(overlaps) > 0:
        relationship = "lineage_ambiguous"
        dependency, admission, score = "lineage_ambiguous", "review", 0.35
        rationale.append("site overlap depends on the lineage candidate")
    elif distances and all(value is not None and value <= 2 for value in distances):
        relationship = "shared_local_environment"
        dependency, admission, score = (
            "shared_local_environment",
            "review",
            0.55,
        )
        rationale.append("active sites are lineage-invariant and within two bonds")
    else:
        relationship = "independent_sites"
        dependency, admission, score = "independent_sites", "rejected", 0.0
        rationale.append("active sites are separated by more than two bonds")
    ambiguity_invariant = bool(
        complete_candidates
        and len(set(overlaps)) <= 1
        and len(set(formed_overlaps)) <= 1
        and len(set(transient_counts)) <= 1
        and len(set(replacement_counts)) <= 1
        and len(set(distances)) <= 1
    )
    if link.status != "unique" and not ambiguity_invariant:
        warnings.append("LINEAGE_RELATIONSHIP_NOT_INVARIANT")
    if link.candidate_limit_reached:
        warnings.append("LINEAGE_CANDIDATE_LIMIT_REACHED")
    evidence = CouplingEvidence(
        producer_active_maps=tuple(sorted(producer_active)),
        producer_persistent_active_maps=tuple(sorted(persistent_active)),
        consumer_active_maps=tuple(sorted(consumer_active)),
        consumer_intermediate_active_maps=tuple(sorted(intermediate_active)),
        overlap_counts=overlaps,
        formed_overlap_counts=formed_overlaps,
        transient_bond_counts=transient_counts,
        replacement_bond_counts=replacement_counts,
        minimum_distances=distances,
        ambiguity_invariant=ambiguity_invariant,
        candidate_limit_reached=link.candidate_limit_reached,
        rationale=tuple(rationale),
    )
    return relationship, dependency, admission, score, evidence, tuple(warnings)


def _overall_reaction(
    first: RouteCoreStep,
    second: RouteCoreStep,
    consumer_component_index: Optional[int],
) -> tuple[str, str]:
    first_sides = _reaction_sides(first.reaction_smiles)
    second_sides = _reaction_sides(second.reaction_smiles)
    if first_sides is None or second_sides is None:
        return "", digest("CRSNET1", first.reaction_smiles, second.reaction_smiles)
    inputs = list(first_sides[0].split("."))
    for index, component in enumerate(second_sides[0].split(".")):
        if index != consumer_component_index:
            inputs.append(component)
    canonical_inputs = canonical_smiles(".".join(inputs))
    canonical_product = canonical_smiles(second_sides[1])
    if canonical_inputs is None or canonical_product is None:
        return "", digest("CRSNET1", first.reaction_smiles, second.reaction_smiles)
    reaction = f"{canonical_inputs}>>{canonical_product}"
    return reaction, digest("CRSNET1", reaction)


def extract_coupled_route_strategies(
    projection: RouteCoreProjection,
) -> tuple[CoupledRouteStrategyOccurrence, ...]:
    """Classify every producer-consumer pair in one route-core projection."""

    steps = {step.reaction_node_id: step for step in projection.steps}
    values = []
    for link in projection.lineage_links:
        first = steps.get(link.producer_reaction_node_id)
        second = steps.get(link.consumer_reaction_node_id)
        if first is None or second is None:
            continue
        (
            relationship,
            dependency,
            admission,
            score,
            evidence,
            warnings,
        ) = _classify(link, first, second)
        exact_seed = "|".join(
            (
                str(first.exact_core_key or first.reaction_signature_id or ""),
                str(second.exact_core_key or second.reaction_signature_id or ""),
                dependency,
            )
        )
        typed_seed = "|".join(
            (
                str(first.typed_core_key or first.transformation_class or ""),
                str(second.typed_core_key or second.transformation_class or ""),
                dependency,
            )
        )
        shape_seed = "|".join(
            (
                str(first.shape_core_key or ""),
                str(second.shape_core_key or ""),
                dependency,
            )
        )
        exact_id = digest("CRSE2", exact_seed)
        typed_id = digest("CRST2", typed_seed)
        shape_id = digest("CRSS2", shape_seed)
        occurrence_id = digest(
            "CRSO2",
            projection.route_core_id,
            link.lineage_id,
            exact_id,
        )
        overall_reaction, overall_id = _overall_reaction(
            first, second, link.consumer_reactant_component_index
        )
        second_sides = _reaction_sides(second.reaction_smiles)
        final_product = (
            canonical_smiles(second_sides[1])
            if second_sides is not None
            else None
        )
        values.append(
            CoupledRouteStrategyOccurrence(
                occurrence_id=occurrence_id,
                exact_strategy_id=exact_id,
                typed_strategy_id=typed_id,
                shape_strategy_id=shape_id,
                overall_reaction_id=overall_id,
                route_core_id=projection.route_core_id,
                source_tree_id=projection.source_tree_id,
                source_route_id=projection.source_route_id,
                patent_id=projection.patent_id,
                split=projection.split,
                lineage_id=link.lineage_id,
                lineage_status=link.status,
                relationship_class=relationship,
                dependency_class=dependency,
                admission_class=admission,
                coupling_score=score,
                first_reaction_node_id=first.reaction_node_id,
                second_reaction_node_id=second.reaction_node_id,
                first_source_reaction_id=first.source_reaction_id,
                second_source_reaction_id=second.source_reaction_id,
                first_reaction_smiles=first.reaction_smiles,
                second_reaction_smiles=second.reaction_smiles,
                intermediate_smiles=link.intermediate_smiles,
                final_product_smiles=final_product or "",
                overall_reaction_smiles=overall_reaction,
                first_transformation_class=first.transformation_class,
                second_transformation_class=second.transformation_class,
                first_named_family=first.named_family,
                second_named_family=second.named_family,
                first_edit_tokens=first.edit_tokens,
                second_edit_tokens=second.edit_tokens,
                combined_edit_tokens=tuple(
                    [f"step1:{token}" for token in first.edit_tokens]
                    + [f"step2:{token}" for token in second.edit_tokens]
                ),
                evidence=evidence,
                warnings=tuple(sorted(set(warnings + link.warnings))),
            )
        )
    return tuple(sorted(values, key=lambda item: item.occurrence_id))


def _sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def _diverse_sample(
    occurrences: Iterable[CoupledRouteStrategyOccurrence],
    size: int,
    *,
    seed: int,
    required_route_ids: Sequence[str] = (),
    required_reaction_pairs: Sequence[tuple[str, str]] = (),
) -> tuple[CoupledRouteStrategyOccurrence, ...]:
    values = tuple(occurrences)
    if size < 1 or not values:
        return ()
    ranked = sorted(
        values,
        key=lambda item: digest("CRSSAMPLE1", str(seed), item.occurrence_id),
    )
    selected: list[CoupledRouteStrategyOccurrence] = []
    selected_ids: set[str] = set()
    used_strategies: set[str] = set()
    for first_reaction_id, second_reaction_id in required_reaction_pairs:
        match = next(
            (
                item
                for item in ranked
                if item.first_source_reaction_id == first_reaction_id
                and item.second_source_reaction_id == second_reaction_id
                and item.occurrence_id not in selected_ids
            ),
            None,
        )
        if match is not None and len(selected) < size:
            selected.append(match)
            selected_ids.add(match.occurrence_id)
            used_strategies.add(match.typed_strategy_id)
    for route_id in required_route_ids:
        match = next(
            (
                item
                for item in ranked
                if item.source_route_id == route_id
                and item.occurrence_id not in selected_ids
            ),
            None,
        )
        if match is not None and len(selected) < size:
            selected.append(match)
            selected_ids.add(match.occurrence_id)
            used_strategies.add(match.typed_strategy_id)
    for item in ranked:
        if len(selected) >= size:
            break
        if (
            item.occurrence_id not in selected_ids
            and item.typed_strategy_id not in used_strategies
        ):
            selected.append(item)
            selected_ids.add(item.occurrence_id)
            used_strategies.add(item.typed_strategy_id)
    for item in ranked:
        if len(selected) >= size:
            break
        if item.occurrence_id not in selected_ids:
            selected.append(item)
            selected_ids.add(item.occurrence_id)
    return tuple(selected)


def mine_coupled_route_strategy_poc(
    source: str | Path,
    *,
    strict_sample_size: int = 30,
    review_sample_size: int = 20,
    rejected_sample_size: int = 20,
    seed: int = 20260818,
    required_route_ids: Sequence[str] = (),
    required_reaction_pairs: Sequence[tuple[str, str]] = (),
    max_routes: Optional[int] = None,
) -> dict[str, Any]:
    """Mine, aggregate, and deterministically sample coupled route strategies."""

    source_path = Path(source)
    occurrences: list[CoupledRouteStrategyOccurrence] = []
    route_count = 0
    for projection in iter_route_core_projections(source_path):
        if max_routes is not None and route_count >= max_routes:
            break
        route_count += 1
        occurrences.extend(extract_coupled_route_strategies(projection))
    relationship_counts = Counter(item.relationship_class for item in occurrences)
    dependency_counts = Counter(item.dependency_class for item in occurrences)
    admission_counts = Counter(item.admission_class for item in occurrences)
    comparison_counts = Counter(
        (item.relationship_class, item.dependency_class) for item in occurrences
    )
    grouped: dict[str, list[CoupledRouteStrategyOccurrence]] = defaultdict(list)
    for item in occurrences:
        if item.admission_class == "strict":
            grouped[item.typed_strategy_id].append(item)
    strategy_summaries = []
    for strategy_id, items in grouped.items():
        patents = sorted({item.patent_id for item in items if item.patent_id})
        routes = sorted(
            {item.source_route_id for item in items if item.source_route_id}
        )
        strategy_summaries.append(
            {
                "typed_strategy_id": strategy_id,
                "shape_strategy_id": items[0].shape_strategy_id,
                "relationship_class": items[0].relationship_class,
                "dependency_class": items[0].dependency_class,
                "occurrence_count": len(items),
                "patent_count": len(patents),
                "route_count": len(routes),
                "patent_ids": patents,
                "example_occurrence_id": items[0].occurrence_id,
                "first_edit_tokens": list(items[0].first_edit_tokens),
                "second_edit_tokens": list(items[0].second_edit_tokens),
            }
        )
    strategy_summaries.sort(
        key=lambda item: (
            -int(item["patent_count"]),
            -int(item["occurrence_count"]),
            str(item["typed_strategy_id"]),
        )
    )
    strict = [item for item in occurrences if item.admission_class == "strict"]
    review = [item for item in occurrences if item.admission_class == "review"]
    rejected = [
        item for item in occurrences if item.admission_class == "rejected"
    ]
    sample = (
        _diverse_sample(
            strict,
            strict_sample_size,
            seed=seed,
            required_route_ids=required_route_ids,
            required_reaction_pairs=required_reaction_pairs,
        )
        + _diverse_sample(review, review_sample_size, seed=seed + 1)
        + _diverse_sample(rejected, rejected_sample_size, seed=seed + 2)
    )
    recurring = sum(
        1 for item in strategy_summaries if int(item["patent_count"]) >= 2
    )
    return {
        "artifact_type": "coupled_route_strategy_poc",
        "schema_version": COUPLED_ROUTE_STRATEGY_SCHEMA_VERSION,
        "algorithm_version": COUPLED_ROUTE_STRATEGY_ALGORITHM_VERSION,
        "source": str(source_path),
        "source_sha256": _sha256(source_path),
        "sample_seed": seed,
        "route_count": route_count,
        "lineage_pair_count": len(occurrences),
        "relationship_counts": dict(sorted(relationship_counts.items())),
        "dependency_counts": dict(sorted(dependency_counts.items())),
        "baseline_to_dependency_counts": [
            {
                "relationship_class": relationship,
                "dependency_class": dependency,
                "count": count,
            }
            for (relationship, dependency), count in sorted(
                comparison_counts.items()
            )
        ],
        "admission_counts": dict(sorted(admission_counts.items())),
        "strict_typed_strategy_count": len(strategy_summaries),
        "recurrent_strict_strategy_count": recurring,
        "strategy_summaries": strategy_summaries,
        "sample_counts": dict(Counter(item.admission_class for item in sample)),
        "sample": [item.to_dict() for item in sample],
    }


def write_coupled_route_strategy_report(
    report: Mapping[str, Any], output_path: str | Path
) -> dict[str, Any]:
    """Write the deterministic mining report."""

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(dict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return {
        "output_json": str(destination.resolve()),
        "json_bytes": destination.stat().st_size,
        "route_count": int(report.get("route_count") or 0),
        "lineage_pair_count": int(report.get("lineage_pair_count") or 0),
        "sample_count": len(report.get("sample") or ()),
    }


__all__ = [
    "ADMISSION_CLASSES",
    "COUPLED_ROUTE_STRATEGY_ALGORITHM_VERSION",
    "COUPLED_ROUTE_STRATEGY_REVIEW_VERSION",
    "COUPLED_ROUTE_STRATEGY_SCHEMA_VERSION",
    "DEPENDENCY_CLASSES",
    "RELATIONSHIP_CLASSES",
    "CoupledRouteStrategyOccurrence",
    "CouplingEvidence",
    "extract_coupled_route_strategies",
    "mine_coupled_route_strategy_poc",
    "write_coupled_route_strategy_report",
]
