"""Canonical molecule-to-reaction route trees and deterministic distances."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from typing import Any, Iterable, Optional, Protocol

from .chemistry import digest


ROUTE_TREE_SCHEMA_VERSION = "1.0"


class RouteStepLike(Protocol):
    """Structural fields needed to build a route tree."""

    step_id: str
    depth: int
    product_smiles: str
    product_node_id: str
    precursor_smiles: tuple[str, ...]
    precursor_node_ids: tuple[str, ...]
    candidate: Any


class RouteLeafLike(Protocol):
    """Terminal fields needed to build a route tree."""

    route_node_id: str
    canonical_smiles: str
    depth: int
    terminal: bool
    terminal_evidence: str
    unresolved_reason: Optional[str]


@dataclass(frozen=True)
class RouteTreeReactionNode:
    """A validated retrosynthetic action and its precursor molecule nodes."""

    reaction_node_id: str
    step_id: str
    depth: int
    proposed_reaction_smiles: str
    operator_id: str
    disconnection_site_key: str
    children: tuple["RouteTreeMoleculeNode", ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible reaction subtree."""

        return {
            "reaction_node_id": self.reaction_node_id,
            "step_id": self.step_id,
            "depth": self.depth,
            "proposed_reaction_smiles": self.proposed_reaction_smiles,
            "operator_id": self.operator_id,
            "disconnection_site_key": self.disconnection_site_key,
            "children": [child.to_dict() for child in self.children],
        }


@dataclass(frozen=True)
class RouteTreeMoleculeNode:
    """One molecule occurrence in a route, distinct from molecular identity."""

    molecule_node_id: str
    smiles: str
    depth: int
    terminal: bool
    terminal_evidence: str
    unresolved_reason: Optional[str]
    reaction: Optional[RouteTreeReactionNode] = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible molecule subtree."""

        value = asdict(self)
        value["reaction"] = self.reaction.to_dict() if self.reaction else None
        return value


@dataclass(frozen=True)
class CanonicalRouteTree:
    """A versioned, occurrence-preserving route-tree representation."""

    tree_id: str
    target_smiles: str
    root: RouteTreeMoleculeNode
    reaction_count: int
    maximum_depth: int
    fingerprint_tokens: tuple[str, ...]
    schema_version: str = ROUTE_TREE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible tree."""

        return {
            "tree_id": self.tree_id,
            "target_smiles": self.target_smiles,
            "root": self.root.to_dict(),
            "reaction_count": self.reaction_count,
            "maximum_depth": self.maximum_depth,
            "fingerprint_tokens": list(self.fingerprint_tokens),
            "schema_version": self.schema_version,
        }


def build_canonical_route_tree(
    target_smiles: str,
    root_node_id: str,
    steps: Iterable[RouteStepLike],
    leaves: Iterable[RouteLeafLike],
) -> CanonicalRouteTree:
    """Build a stable alternating molecule/reaction tree from route edges."""

    ordered_steps = tuple(steps)
    step_by_product = {step.product_node_id: step for step in ordered_steps}
    leaf_by_node = {leaf.route_node_id: leaf for leaf in leaves}

    def molecule_node(node_id: str, smiles: str, depth: int) -> RouteTreeMoleculeNode:
        step = step_by_product.get(node_id)
        leaf = leaf_by_node.get(node_id)
        reaction = None
        if step is not None:
            children = tuple(
                molecule_node(child_id, child_smiles, step.depth)
                for child_id, child_smiles in zip(
                    step.precursor_node_ids,
                    step.precursor_smiles,
                )
            )
            reaction = RouteTreeReactionNode(
                reaction_node_id=digest("RTREACTION1", node_id, step.step_id),
                step_id=step.step_id,
                depth=step.depth,
                proposed_reaction_smiles=step.candidate.proposed_reaction_smiles,
                operator_id=step.candidate.operator_id,
                disconnection_site_key=step.candidate.disconnection_site_key,
                children=children,
            )
        return RouteTreeMoleculeNode(
            molecule_node_id=node_id,
            smiles=smiles,
            depth=depth,
            terminal=bool(leaf and leaf.terminal),
            terminal_evidence=(leaf.terminal_evidence if leaf else "expanded"),
            unresolved_reason=(leaf.unresolved_reason if leaf else None),
            reaction=reaction,
        )

    root = molecule_node(root_node_id, target_smiles, 0)
    tokens = tuple(
        sorted(
            (
                f"reaction:{step.depth}:{step.candidate.operator_id}:"
                f"{step.candidate.disconnection_site_key}:"
                f"{'.'.join(sorted(step.precursor_smiles))}"
            )
            for step in ordered_steps
        )
    )
    maximum_depth = max((step.depth for step in ordered_steps), default=0)
    tree_id = digest(
        "ROUTETREE1",
        target_smiles,
        str(maximum_depth),
        *tokens,
    )
    return CanonicalRouteTree(
        tree_id=tree_id,
        target_smiles=target_smiles,
        root=root,
        reaction_count=len(ordered_steps),
        maximum_depth=maximum_depth,
        fingerprint_tokens=tokens,
    )


def route_tree_distance(
    left: CanonicalRouteTree,
    right: CanonicalRouteTree,
) -> float:
    """Return a deterministic 0..1 multiset-Jaccard route distance."""

    left_tokens = Counter(left.fingerprint_tokens)
    right_tokens = Counter(right.fingerprint_tokens)
    union = sum((left_tokens | right_tokens).values())
    if union == 0:
        return 0.0 if left.target_smiles == right.target_smiles else 1.0
    intersection = sum((left_tokens & right_tokens).values())
    return round(1.0 - (intersection / union), 8)


def route_distance_matrix(
    trees: Iterable[CanonicalRouteTree],
) -> tuple[tuple[float, ...], ...]:
    """Return a symmetric distance matrix in the supplied tree order."""

    values = tuple(trees)
    return tuple(
        tuple(route_tree_distance(left, right) for right in values)
        for left in values
    )


__all__ = [
    "ROUTE_TREE_SCHEMA_VERSION",
    "CanonicalRouteTree",
    "RouteTreeMoleculeNode",
    "RouteTreeReactionNode",
    "build_canonical_route_tree",
    "route_distance_matrix",
    "route_tree_distance",
]
