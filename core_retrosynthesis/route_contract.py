"""Evidence-neutral, occurrence-preserving retrosynthesis route trees."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterator, Optional


ROUTE_TREE_SCHEMA_VERSION = "2.0"
ROUTE_KINDS = frozenset({"observed", "planned"})
ROUTE_EVIDENCE_KINDS = frozenset({"observed", "inferred", "predicted"})


@dataclass(frozen=True)
class RouteStepEvidence:
    """Provenance and uncertainty for one route reaction occurrence."""

    evidence_kind: str
    source_dataset_id: Optional[str] = None
    source_route_id: Optional[str] = None
    source_reaction_id: Optional[str] = None
    patent_id: Optional[str] = None
    connectivity_method: str = "explicit"
    confidence: Optional[float] = None
    reactants_smiles: str = ""
    reagents_smiles: str = ""
    product_smiles_mapped: str = ""
    abstracted_reaction_smiles: str = ""
    abstraction_status: str = "none"
    reaction_signature_id: Optional[str] = None
    reaction_signature_schema_version: Optional[str] = None
    condition_recipe_id: Optional[str] = None
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.evidence_kind not in ROUTE_EVIDENCE_KINDS:
            raise ValueError(
                f"Unsupported route evidence kind: {self.evidence_kind}"
            )
        if self.confidence is not None and not 0.0 <= self.confidence <= 1.0:
            raise ValueError("Route evidence confidence must be in [0, 1]")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible evidence record."""

        value = asdict(self)
        value["warnings"] = list(self.warnings)
        return value

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "RouteStepEvidence":
        """Reconstruct typed evidence from its serialized contract."""

        return cls(
            evidence_kind=str(value["evidence_kind"]),
            source_dataset_id=value.get("source_dataset_id"),
            source_route_id=value.get("source_route_id"),
            source_reaction_id=value.get("source_reaction_id"),
            patent_id=value.get("patent_id"),
            connectivity_method=str(value.get("connectivity_method") or "explicit"),
            confidence=value.get("confidence"),
            reactants_smiles=str(value.get("reactants_smiles") or ""),
            reagents_smiles=str(value.get("reagents_smiles") or ""),
            product_smiles_mapped=str(value.get("product_smiles_mapped") or ""),
            abstracted_reaction_smiles=str(
                value.get("abstracted_reaction_smiles") or ""
            ),
            abstraction_status=str(value.get("abstraction_status") or "none"),
            reaction_signature_id=value.get("reaction_signature_id"),
            reaction_signature_schema_version=value.get(
                "reaction_signature_schema_version"
            ),
            condition_recipe_id=value.get("condition_recipe_id"),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
        )


@dataclass(frozen=True)
class PlannedRouteAction:
    """Optional prediction-specific annotation on a neutral reaction node."""

    operator_id: str
    disconnection_site_key: str
    template_id: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible action annotation."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "PlannedRouteAction":
        """Reconstruct a planned action annotation."""

        return cls(
            operator_id=str(value["operator_id"]),
            disconnection_site_key=str(value["disconnection_site_key"]),
            template_id=value.get("template_id"),
        )


@dataclass(frozen=True)
class RouteReactionNode:
    """One reaction occurrence and all precursor molecule occurrences."""

    reaction_node_id: str
    step_id: str
    depth: int
    reaction_smiles: str
    evidence: RouteStepEvidence
    children: tuple["MoleculeOccurrenceNode", ...]
    planned_action: Optional[PlannedRouteAction] = None

    @property
    def proposed_reaction_smiles(self) -> str:
        """Compatibility view for the former planned-route-only contract."""

        return self.reaction_smiles

    @property
    def operator_id(self) -> Optional[str]:
        """Return the planned operator, when this is a predicted action."""

        return self.planned_action.operator_id if self.planned_action else None

    @property
    def disconnection_site_key(self) -> Optional[str]:
        """Return the planned disconnection site, when available."""

        return (
            self.planned_action.disconnection_site_key
            if self.planned_action
            else None
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible reaction subtree."""

        return {
            "reaction_node_id": self.reaction_node_id,
            "step_id": self.step_id,
            "depth": self.depth,
            "reaction_smiles": self.reaction_smiles,
            # Compatibility field retained through the v2 migration.
            "proposed_reaction_smiles": self.reaction_smiles,
            "operator_id": self.operator_id,
            "disconnection_site_key": self.disconnection_site_key,
            "evidence": self.evidence.to_dict(),
            "planned_action": (
                self.planned_action.to_dict()
                if self.planned_action is not None
                else None
            ),
            "children": [child.to_dict() for child in self.children],
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "RouteReactionNode":
        """Reconstruct a reaction occurrence and its precursor subtree."""

        evidence = value.get("evidence")
        if not isinstance(evidence, dict):
            raise ValueError("Route reaction evidence must be an object")
        planned_action = value.get("planned_action")
        children = value.get("children")
        if not isinstance(children, list):
            raise ValueError("Route reaction children must be a list")
        return cls(
            reaction_node_id=str(value["reaction_node_id"]),
            step_id=str(value["step_id"]),
            depth=int(value["depth"]),
            reaction_smiles=str(
                value.get("reaction_smiles")
                or value.get("proposed_reaction_smiles")
                or ""
            ),
            evidence=RouteStepEvidence.from_dict(evidence),
            children=tuple(
                MoleculeOccurrenceNode.from_dict(child) for child in children
            ),
            planned_action=(
                PlannedRouteAction.from_dict(planned_action)
                if isinstance(planned_action, dict)
                else None
            ),
        )


@dataclass(frozen=True)
class MoleculeOccurrenceNode:
    """One molecule occurrence, distinct from canonical molecular identity."""

    occurrence_id: str
    smiles: str
    depth: int
    terminal: bool
    terminal_evidence: str
    unresolved_reason: Optional[str]
    reaction: Optional[RouteReactionNode] = None

    @property
    def molecule_node_id(self) -> str:
        """Compatibility alias for the former planned-route field name."""

        return self.occurrence_id

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible molecule subtree."""

        return {
            "occurrence_id": self.occurrence_id,
            # Compatibility field retained through the v2 migration.
            "molecule_node_id": self.occurrence_id,
            "smiles": self.smiles,
            "depth": self.depth,
            "terminal": self.terminal,
            "terminal_evidence": self.terminal_evidence,
            "unresolved_reason": self.unresolved_reason,
            "reaction": self.reaction.to_dict() if self.reaction else None,
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "MoleculeOccurrenceNode":
        """Reconstruct a molecule occurrence and optional reaction subtree."""

        reaction = value.get("reaction")
        occurrence_id = value.get("occurrence_id") or value.get(
            "molecule_node_id"
        )
        if not occurrence_id:
            raise ValueError("Molecule occurrence ID is required")
        return cls(
            occurrence_id=str(occurrence_id),
            smiles=str(value["smiles"]),
            depth=int(value["depth"]),
            terminal=bool(value["terminal"]),
            terminal_evidence=str(value.get("terminal_evidence") or "none"),
            unresolved_reason=value.get("unresolved_reason"),
            reaction=(
                RouteReactionNode.from_dict(reaction)
                if isinstance(reaction, dict)
                else None
            ),
        )


@dataclass(frozen=True)
class ReactionRouteTree:
    """One concrete observed or planned alternating retrosynthesis tree."""

    tree_id: str
    route_kind: str
    target_smiles: str
    root: MoleculeOccurrenceNode
    reaction_count: int
    maximum_depth: int
    fingerprint_tokens: tuple[str, ...]
    source_dataset_id: Optional[str] = None
    source_route_id: Optional[str] = None
    patent_id: Optional[str] = None
    split: Optional[str] = None
    source_record_schema_version: Optional[str] = None
    connectivity_method: str = "explicit"
    confidence: Optional[float] = None
    higher_level_reaction_count: Optional[int] = None
    higher_level_depth: Optional[int] = None
    abstraction_reduction: Optional[int] = None
    warnings: tuple[str, ...] = ()
    schema_version: str = ROUTE_TREE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.route_kind not in ROUTE_KINDS:
            raise ValueError(f"Unsupported route kind: {self.route_kind}")
        if self.confidence is not None and not 0.0 <= self.confidence <= 1.0:
            raise ValueError("Route confidence must be in [0, 1]")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible nested route tree."""

        return {
            "tree_id": self.tree_id,
            "route_kind": self.route_kind,
            "target_smiles": self.target_smiles,
            "root": self.root.to_dict(),
            "reaction_count": self.reaction_count,
            "maximum_depth": self.maximum_depth,
            "fingerprint_tokens": list(self.fingerprint_tokens),
            "source_dataset_id": self.source_dataset_id,
            "source_route_id": self.source_route_id,
            "patent_id": self.patent_id,
            "split": self.split,
            "source_record_schema_version": self.source_record_schema_version,
            "connectivity_method": self.connectivity_method,
            "confidence": self.confidence,
            "higher_level_reaction_count": self.higher_level_reaction_count,
            "higher_level_depth": self.higher_level_depth,
            "abstraction_reduction": self.abstraction_reduction,
            "warnings": list(self.warnings),
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "ReactionRouteTree":
        """Reconstruct and validate one serialized v2 route tree."""

        root = value.get("root")
        if not isinstance(root, dict):
            raise ValueError("Route root must be an object")
        tree = cls(
            tree_id=str(value["tree_id"]),
            route_kind=str(value["route_kind"]),
            target_smiles=str(value["target_smiles"]),
            root=MoleculeOccurrenceNode.from_dict(root),
            reaction_count=int(value["reaction_count"]),
            maximum_depth=int(value["maximum_depth"]),
            fingerprint_tokens=tuple(
                str(item) for item in value.get("fingerprint_tokens") or ()
            ),
            source_dataset_id=value.get("source_dataset_id"),
            source_route_id=value.get("source_route_id"),
            patent_id=value.get("patent_id"),
            split=value.get("split"),
            source_record_schema_version=value.get(
                "source_record_schema_version"
            ),
            connectivity_method=str(value.get("connectivity_method") or "explicit"),
            confidence=value.get("confidence"),
            higher_level_reaction_count=value.get(
                "higher_level_reaction_count"
            ),
            higher_level_depth=value.get("higher_level_depth"),
            abstraction_reduction=value.get("abstraction_reduction"),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            schema_version=str(
                value.get("schema_version") or ROUTE_TREE_SCHEMA_VERSION
            ),
        )
        assert_valid_route_tree(tree)
        return tree


@dataclass(frozen=True)
class RouteTreeValidationReport:
    """Structural validation result for one alternating route tree."""

    valid: bool
    molecule_occurrence_count: int
    reaction_count: int
    maximum_depth: int
    issues: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible validation report."""

        value = asdict(self)
        value["issues"] = list(self.issues)
        return value


def iter_molecule_occurrences(
    tree: ReactionRouteTree,
) -> Iterator[MoleculeOccurrenceNode]:
    """Yield molecule occurrences in deterministic depth-first order."""

    def visit(node: MoleculeOccurrenceNode) -> Iterator[MoleculeOccurrenceNode]:
        yield node
        if node.reaction is not None:
            for child in node.reaction.children:
                yield from visit(child)

    yield from visit(tree.root)


def iter_route_reactions(tree: ReactionRouteTree) -> Iterator[RouteReactionNode]:
    """Yield reaction occurrences in deterministic depth-first order."""

    for molecule in iter_molecule_occurrences(tree):
        if molecule.reaction is not None:
            yield molecule.reaction


def validate_route_tree(tree: ReactionRouteTree) -> RouteTreeValidationReport:
    """Validate alternating topology, occurrence identity, and depth fields."""

    issues: list[str] = []
    molecule_ids: set[str] = set()
    reaction_ids: set[str] = set()
    reaction_count = 0
    molecule_count = 0
    maximum_depth = 0

    def issue(code: str, identifier: str) -> None:
        issues.append(f"{code}:{identifier}")

    def visit(node: MoleculeOccurrenceNode, expected_depth: int) -> None:
        nonlocal molecule_count, reaction_count, maximum_depth
        molecule_count += 1
        maximum_depth = max(maximum_depth, node.depth)
        if node.occurrence_id in molecule_ids:
            issue("duplicate_molecule_occurrence_id", node.occurrence_id)
            return
        molecule_ids.add(node.occurrence_id)
        if node.depth != expected_depth:
            issue("molecule_depth_mismatch", node.occurrence_id)
        if node.terminal and node.reaction is not None:
            issue("terminal_molecule_has_reaction", node.occurrence_id)
        reaction = node.reaction
        if reaction is None:
            return
        reaction_count += 1
        if reaction.reaction_node_id in reaction_ids:
            issue("duplicate_reaction_node_id", reaction.reaction_node_id)
            return
        reaction_ids.add(reaction.reaction_node_id)
        if reaction.depth != node.depth + 1:
            issue("reaction_depth_mismatch", reaction.reaction_node_id)
        if not reaction.children:
            issue("reaction_without_precursors", reaction.reaction_node_id)
        for child in reaction.children:
            visit(child, reaction.depth)

    visit(tree.root, 0)
    if tree.target_smiles != tree.root.smiles:
        issue("target_root_identity_mismatch", tree.tree_id)
    if tree.reaction_count != reaction_count:
        issue("reaction_count_mismatch", tree.tree_id)
    if tree.maximum_depth != maximum_depth:
        issue("maximum_depth_mismatch", tree.tree_id)
    if tuple(sorted(tree.fingerprint_tokens)) != tree.fingerprint_tokens:
        issue("fingerprint_tokens_not_sorted", tree.tree_id)
    return RouteTreeValidationReport(
        valid=not issues,
        molecule_occurrence_count=molecule_count,
        reaction_count=reaction_count,
        maximum_depth=maximum_depth,
        issues=tuple(issues),
    )


def assert_valid_route_tree(tree: ReactionRouteTree) -> None:
    """Raise when a route violates the shared tree contract."""

    report = validate_route_tree(tree)
    if not report.valid:
        raise ValueError("; ".join(report.issues))


__all__ = [
    "ROUTE_EVIDENCE_KINDS",
    "ROUTE_KINDS",
    "ROUTE_TREE_SCHEMA_VERSION",
    "MoleculeOccurrenceNode",
    "PlannedRouteAction",
    "ReactionRouteTree",
    "RouteReactionNode",
    "RouteStepEvidence",
    "RouteTreeValidationReport",
    "assert_valid_route_tree",
    "iter_molecule_occurrences",
    "iter_route_reactions",
    "validate_route_tree",
]
