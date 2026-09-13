"""Explain rejected core matches without granting retrieval eligibility.

Counterfactual graph comparisons are review aids, not chemistry relaxations.
They retain graph connectivity and expose which label changes would be needed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

from rdkit import Chem

from .shared_core_graph import canonical_graph
from .shared_reaction_core import SharedReactionCore, compare_reaction_cores


@dataclass(frozen=True)
class SharedCoreDiagnostic:
    """Observed comparison and explicitly non-authorizing counterfactuals."""

    eligible: bool
    classification: str
    matching_counterfactuals: tuple[str, ...]
    reasons: tuple[str, ...]
    requires_review: bool = True


def _decode_graph(payload: str) -> tuple[list[Any], list[tuple[int, int, Any]]]:
    vocabulary, smiles = json.loads(payload)
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid canonical core graph")
    nodes = {}
    labels = []
    edge_nodes = []
    for atom in molecule.GetAtoms():
        kind, label = json.loads(vocabulary[atom.GetIsotope() - 1])
        if kind == "atom":
            nodes[atom.GetIdx()] = len(labels)
            labels.append(label)
        elif kind == "edge":
            edge_nodes.append((atom, label))
        else:
            raise ValueError("Unknown canonical graph node kind")
    edges = []
    for atom, label in edge_nodes:
        neighbors = tuple(atom.GetNeighbors())
        if len(neighbors) != 2:
            raise ValueError("Invalid canonical graph incidence edge")
        edges.append((nodes[neighbors[0].GetIdx()], nodes[neighbors[1].GetIdx()], label))
    return labels, edges


def _counterfactual_key(payload: str, profile: str) -> str:
    labels, edges = _decode_graph(payload)
    for label in labels:
        if not isinstance(label, list) or label[0] != "center":
            continue
        before, after = label[1:3]
        if before is None or after is None:
            continue
        fields = {
            "unchanged_ring_membership": (8,),
            "unchanged_stereochemistry": (7,),
            "hydrogen_substitution_with_same_delta": (6,),
            "combined_static_context": (6, 7, 8),
        }[profile]
        for field in fields:
            if field == 6:
                if isinstance(before[field], int) and isinstance(after[field], int):
                    delta = after[field] - before[field]
                    before[field], after[field] = "initial_hydrogens", ["delta", delta]
            elif before[field] == after[field]:
                before[field] = after[field] = "unchanged_context"
    return canonical_graph(labels, edges)


def diagnose_shared_core_difference(
    query: SharedReactionCore, precedent: SharedReactionCore
) -> SharedCoreDiagnostic:
    """Distinguish unavailable evidence and selected static-context barriers.

    No counterfactual changes edits, element, charge, aromaticity, hybridization,
    isotope, radical state, event connectivity or center distances. A match here
    must never be used as an eligibility decision in recommendation.
    """
    comparison = compare_reaction_cores(query, precedent)
    if comparison.eligible:
        return SharedCoreDiagnostic(True, "qualified_shared_core", (), comparison.reasons)
    if query.definition_hash != precedent.definition_hash:
        return SharedCoreDiagnostic(False, "definition_mismatch", (), comparison.reasons)
    if not query.levels or not precedent.levels:
        return SharedCoreDiagnostic(
            False, "observation_unavailable", (),
            tuple(sorted(set(query.unavailable_reasons + precedent.unavailable_reasons))),
        )
    profiles = (
        "unchanged_ring_membership", "unchanged_stereochemistry",
        "hydrogen_substitution_with_same_delta", "combined_static_context",
    )
    left = {item.level: item.graph_payload for item in query.levels}
    right = {item.level: item.graph_payload for item in precedent.levels}
    if "retained_typed" not in left or "retained_typed" not in right:
        return SharedCoreDiagnostic(False, "broad_projection_unavailable", (), comparison.reasons)
    matches = tuple(
        profile for profile in profiles
        if _counterfactual_key(left["retained_typed"], profile)
        == _counterfactual_key(right["retained_typed"], profile)
    )
    return SharedCoreDiagnostic(
        False, "static_context_review_candidate" if matches else "protected_graph_or_other_context",
        matches, comparison.reasons,
    )
