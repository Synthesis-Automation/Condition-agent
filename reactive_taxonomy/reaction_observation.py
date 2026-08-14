"""Interpretation-independent assembly of structural reaction observations."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Iterable, Tuple

from .reaction_bond_changes import supplied_map_bond_changes
from .reaction_completeness import build_reaction_completeness
from .reaction_core import build_reaction_core_projection
from .reaction_edits import (
    EditNormalizationResult,
    normalize_mapped_edits,
    resolve_structural_evidence,
)
from .reaction_models import (
    ReactionComponent,
    ReactionObservation,
)
from .reaction_topology import build_reaction_topology
from .reaction_parser import parse_reaction_smiles
from .reaction_stoichiometry import infer_reactant_multiplicity


def _synthetic_atom_maps(
    correspondence: Tuple[Tuple[int, int, int, int], ...],
) -> dict[tuple[str, int, int], int]:
    """Give one correspondence hypothesis deterministic internal atom IDs."""
    overrides: dict[tuple[str, int, int], int] = {}
    for number, (reactant_component, reactant_atom, product_component, product_atom) in enumerate(
        sorted(correspondence),
        start=1,
    ):
        overrides[("reactant", reactant_component, reactant_atom)] = number
        overrides[("product", product_component, product_atom)] = number
    return overrides


def _observation_hypotheses(
    *,
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    edit_result: EditNormalizationResult,
) -> tuple[tuple[Any, ...], tuple[Any, ...]]:
    hypotheses = tuple(
        replace(
            hypothesis,
            topology=build_reaction_topology(
                reactants=reactants,
                products=products,
                edit_result=replace(
                    edit_result,
                    edits=hypothesis.edits,
                    evidence=hypothesis.evidence,
                    confidence=hypothesis.confidence,
                    valid=False,
                    stereo_changes=hypothesis.stereo_changes,
                    edit_hypotheses=(),
                    atom_correspondence=hypothesis.atom_correspondence,
                ),
            ),
        )
        for hypothesis in edit_result.edit_hypotheses
    )
    by_id = {hypothesis.hypothesis_id: hypothesis for hypothesis in hypotheses}
    evidence_candidates = tuple(
        replace(
            candidate,
            edit_hypotheses=tuple(
                by_id.get(hypothesis.hypothesis_id, hypothesis)
                for hypothesis in candidate.edit_hypotheses
            ),
        )
        for candidate in edit_result.evidence_candidates
    )
    return hypotheses, evidence_candidates


def build_reaction_observation(
    *,
    input_reaction_smiles: str,
    reactants: Tuple[ReactionComponent, ...],
    agents: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    edit_result: EditNormalizationResult,
    mapped_bond_changes: Tuple[dict[str, Any], ...] = (),
    warnings: Iterable[str] = (),
) -> ReactionObservation:
    """Finalize graph facts after correspondence providers are reconciled."""
    structural_reactants = tuple(component.structure_only() for component in reactants)
    structural_agents = tuple(component.structure_only() for component in agents)
    structural_products = tuple(component.structure_only() for component in products)
    topology = build_reaction_topology(
        reactants=structural_reactants,
        products=structural_products,
        edit_result=edit_result,
    )
    core = build_reaction_core_projection(
        reactants=structural_reactants,
        products=structural_products,
        edits=edit_result.edits,
        stereo_changes=edit_result.stereo_changes,
        evidence=edit_result.evidence,
        confidence=edit_result.confidence,
        atom_map_overrides=(
            _synthetic_atom_maps(edit_result.atom_correspondence)
            if edit_result.atom_correspondence
            else None
        ),
    )
    completeness = build_reaction_completeness(
        reactants=structural_reactants,
        products=structural_products,
        edit_result=edit_result,
    )
    hypotheses, evidence_candidates = _observation_hypotheses(
        reactants=structural_reactants,
        products=structural_products,
        edit_result=edit_result,
    )
    if core is None and hypotheses:
        hypothesis_cores = tuple(
            candidate
            for hypothesis in hypotheses
            for candidate in (
                build_reaction_core_projection(
                    reactants=structural_reactants,
                    products=structural_products,
                    edits=hypothesis.edits,
                    stereo_changes=hypothesis.stereo_changes,
                    evidence=hypothesis.evidence,
                    confidence=hypothesis.confidence,
                    atom_map_overrides=_synthetic_atom_maps(
                        hypothesis.atom_correspondence
                    ),
                ),
            )
            if candidate is not None
        )
        shape_keys = {candidate.shape_core_key for candidate in hypothesis_cores}
        typed_keys = {candidate.typed_core_key for candidate in hypothesis_cores}
        if hypothesis_cores and len(shape_keys) == 1 and len(typed_keys) == 1:
            core = replace(
                min(hypothesis_cores, key=lambda candidate: candidate.core_id),
                evidence="ambiguous_correspondence_core_consensus",
                evidence_status="hypothesis",
                warnings=tuple(
                    sorted(
                        set(
                            min(
                                hypothesis_cores,
                                key=lambda candidate: candidate.core_id,
                            ).warnings
                        ).union({"CORE_SHARED_BY_ALL_EDIT_HYPOTHESES"})
                    )
                ),
            )
    combined_warnings = tuple(
        sorted(
            set(warnings)
            .union(edit_result.warnings)
            .union(completeness.warnings)
        )
    )
    return ReactionObservation(
        input_reaction_smiles=input_reaction_smiles,
        valid=True,
        reactants=structural_reactants,
        agents=structural_agents,
        products=structural_products,
        edits=edit_result.edits,
        stereo_changes=edit_result.stereo_changes,
        evidence_quality=edit_result.evidence,
        evidence_confidence=edit_result.confidence,
        connectivity_edit_graph=edit_result.connectivity_edit_graph,
        evidence_candidates=evidence_candidates,
        edit_hypotheses=hypotheses,
        mapped_bond_changes=mapped_bond_changes,
        topology=topology,
        completeness=completeness,
        core=core,
        warnings=combined_warnings,
    )


def observe_reaction(reaction_smiles: str) -> ReactionObservation:
    """Build the type-agnostic foundation without interpretation metadata."""
    parsed = parse_reaction_smiles(
        reaction_smiles,
        include_molecular_interpretation=False,
    )
    if not parsed.valid:
        return ReactionObservation(
            input_reaction_smiles=reaction_smiles,
            valid=False,
            reactants=tuple(component.structure_only() for component in parsed.reactants),
            agents=tuple(component.structure_only() for component in parsed.agents),
            products=tuple(component.structure_only() for component in parsed.products),
            warnings=parsed.warnings,
            error=parsed.error,
        )
    multiplicity = infer_reactant_multiplicity(
        parsed.reactants,
        parsed.products,
    )
    if multiplicity.inferred_copy_count:
        parsed = replace(
            parsed,
            reactants=multiplicity.reactants,
            warnings=tuple(
                sorted(set(parsed.warnings + multiplicity.warnings))
            ),
        )
    supplied_mapping = normalize_mapped_edits(
        parsed.reactants, parsed.products
    )
    edit_result = resolve_structural_evidence(
        parsed.reactants,
        parsed.products,
        mapped_override=supplied_mapping,
        additional_warnings=multiplicity.warnings,
    )
    return build_reaction_observation(
        input_reaction_smiles=reaction_smiles,
        reactants=parsed.reactants,
        agents=parsed.agents,
        products=parsed.products,
        edit_result=edit_result,
        mapped_bond_changes=tuple(
            supplied_map_bond_changes(reaction_smiles)
        ),
        warnings=parsed.warnings,
    )


__all__ = ["build_reaction_observation", "observe_reaction"]
