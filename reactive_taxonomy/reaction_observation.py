"""Interpretation-independent assembly of structural reaction observations."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Iterable, Mapping, Sequence, Tuple

from .reaction_bond_changes import supplied_map_bond_changes
from .reaction_completeness import build_reaction_completeness
from .reaction_core import build_reaction_core_projection
from .reaction_edits import (
    EditNormalizationResult,
    normalize_mapped_edits,
    resolve_reaction_evidence,
)
from .reaction_models import (
    ReactionComponent,
    ReactionObservation,
    ReactionReconstructionCandidate,
    ReactionSiteReference,
)
from .reaction_spectators import derive_observed_spectator_groups
from .reaction_topology import build_reaction_topology
from .reaction_parser import parse_reaction_smiles
from .reaction_reconstruction import (
    build_reaction_reconstruction_candidates,
    canonical_without_maps,
)


RawCandidate = Tuple[
    Mapping[str, Any],
    Mapping[str, ReactionSiteReference],
]


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
                selected=None,
                edit_result=replace(
                    edit_result,
                    edits=hypothesis.edits,
                    evidence=hypothesis.evidence,
                    confidence=hypothesis.confidence,
                    valid=False,
                    stereo_changes=hypothesis.stereo_changes,
                    edit_hypotheses=(),
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
    reconstruction_sources: Sequence[RawCandidate] = (),
    reconstruction_candidates: Tuple[ReactionReconstructionCandidate, ...] = (),
    selected_reconstruction: ReactionReconstructionCandidate | None = None,
    selected_reconstruction_events: Tuple[ReactionReconstructionCandidate, ...] = (),
    warnings: Iterable[str] = (),
) -> ReactionObservation:
    """Finalize generic facts after evidence providers have been reconciled.

    Reconstruction candidates are structural operator evidence. They contain no
    interpretation, family, or display semantics.
    """
    topology = build_reaction_topology(
        reactants=reactants,
        products=products,
        selected=None,
        edit_result=edit_result,
    )
    core = build_reaction_core_projection(
        reactants=reactants,
        products=products,
        edits=edit_result.edits,
        stereo_changes=edit_result.stereo_changes,
        evidence=edit_result.evidence,
        confidence=edit_result.confidence,
        topology=topology,
    )
    completeness = build_reaction_completeness(
        reactants=reactants,
        products=products,
        raw_candidates=reconstruction_sources,
        selected=selected_reconstruction,
        selected_events=selected_reconstruction_events,
        edit_result=edit_result,
    )
    spectators = derive_observed_spectator_groups(
        reactants,
        edit_result.edits,
        edit_result.evidence,
    )
    hypotheses, evidence_candidates = _observation_hypotheses(
        reactants=reactants,
        products=products,
        edit_result=edit_result,
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
        reactants=reactants,
        agents=agents,
        products=products,
        edits=edit_result.edits,
        stereo_changes=edit_result.stereo_changes,
        evidence_quality=edit_result.evidence,
        evidence_confidence=edit_result.confidence,
        connectivity_edit_graph=edit_result.connectivity_edit_graph,
        evidence_candidates=evidence_candidates,
        edit_hypotheses=hypotheses,
        mapped_bond_changes=mapped_bond_changes,
        reconstruction_candidates=reconstruction_candidates,
        selected_reconstruction=selected_reconstruction,
        selected_reconstruction_events=selected_reconstruction_events,
        spectator_groups=spectators,
        topology=topology,
        completeness=completeness,
        core=core,
        warnings=combined_warnings,
    )


def observe_reaction(reaction_smiles: str) -> ReactionObservation:
    """Build the type-agnostic foundation without interpretation metadata."""
    parsed = parse_reaction_smiles(reaction_smiles)
    if not parsed.valid:
        return ReactionObservation(
            input_reaction_smiles=reaction_smiles,
            valid=False,
            reactants=parsed.reactants,
            agents=parsed.agents,
            products=parsed.products,
            warnings=parsed.warnings,
            error=parsed.error,
        )
    observed_products = {
        canonical
        for component in parsed.products
        for canonical in (canonical_without_maps(component.input_smiles),)
        if canonical is not None
    }
    supplied_mapping = normalize_mapped_edits(
        parsed.reactants, parsed.products
    )
    reconstruction = build_reaction_reconstruction_candidates(
        reactants=parsed.reactants,
        observed_products=observed_products,
        invalid_supplied_mapping=(
            supplied_mapping.evidence == "invalid_atom_mapping"
        ),
        max_candidates=100,
    )
    edit_result = resolve_reaction_evidence(
        parsed.reactants,
        parsed.products,
        selected=reconstruction.selected_candidate,
        selected_events=reconstruction.selected_events,
        candidates=reconstruction.candidates,
        mapped_override=supplied_mapping,
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
        reconstruction_sources=reconstruction.raw_candidates,
        reconstruction_candidates=reconstruction.candidates,
        selected_reconstruction=reconstruction.selected_candidate,
        selected_reconstruction_events=reconstruction.selected_events,
        warnings=parsed.warnings + reconstruction.warnings,
    )


__all__ = ["build_reaction_observation", "observe_reaction"]
