"""Grammar-independent assembly of structural reaction observations."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Iterable, Mapping, Sequence, Tuple

from .reaction_bond_changes import supplied_map_bond_changes
from .reaction_completeness import build_reaction_completeness
from .reaction_core import build_reaction_core_projection
from .reaction_edits import EditNormalizationResult, resolve_reaction_evidence
from .reaction_models import (
    ReactionCandidate,
    ReactionComponent,
    ReactionObservation,
    ReactionSiteReference,
)
from .reaction_spectators import derive_observed_spectator_groups
from .reaction_topology import build_reaction_topology
from .reaction_parser import parse_reaction_smiles


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
    operator_candidates: Sequence[RawCandidate] = (),
    selected_operator_candidate: ReactionCandidate | None = None,
    selected_operator_events: Tuple[ReactionCandidate, ...] = (),
    warnings: Iterable[str] = (),
) -> ReactionObservation:
    """Finalize generic facts after evidence providers have been reconciled.

    Operator candidates are accepted only for product-provenance accounting.
    Their grammar or family semantics are never copied into the observation.
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
        raw_candidates=operator_candidates,
        selected=selected_operator_candidate,
        selected_events=selected_operator_events,
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
        spectator_groups=spectators,
        topology=topology,
        completeness=completeness,
        core=core,
        warnings=combined_warnings,
    )


def observe_reaction(reaction_smiles: str) -> ReactionObservation:
    """Build the type-agnostic reaction foundation without loading grammars."""
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
    edit_result = resolve_reaction_evidence(
        parsed.reactants,
        parsed.products,
        selected=None,
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
