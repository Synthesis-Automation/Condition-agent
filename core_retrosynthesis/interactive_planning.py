"""Deterministic, occurrence-preserving edits for human-selected route plans.

Search evidence is retained as reported, including on import. Session validation
checks molecular identity and topology; it never upgrades a saved prediction to
an observed reaction or independently revalidates its reported chemistry score.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
from typing import Any, Mapping

from rdkit import Chem

from .chemistry import canonical_smiles, digest
from .route_contract import (
    MoleculeOccurrenceNode,
    PlannedRouteAction,
    ReactionRouteTree,
    RouteReactionNode,
    RouteStepEvidence,
    assert_valid_route_tree,
)

SCHEMA_VERSION = "interactive_planning.v1"
MAX_NODES = 200
MAX_DEPTH = 20
MAX_HISTORY = 30
MAX_SEARCHES = 100


def _object(value: Any) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError("Expected a planning object")
    return value


def _text(value: Any) -> str:
    if not isinstance(value, str) or not value or len(value) > 20_000:
        raise ValueError("Expected a nonempty bounded string")
    return value


def _integer(value: Any) -> int:
    if type(value) is not int or value < 0:
        raise ValueError("Candidate indices must be nonnegative integers")
    return value


def _messages(value: Any) -> None:
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        raise ValueError("Saved warnings must be a list of strings")


def _molecule(value: Any) -> str:
    smiles = canonical_smiles(_text(value))
    if not smiles or "." in smiles:
        raise ValueError("Planning requires one connected molecule per node")
    return smiles


def normalize_settings(value: Mapping[str, Any]) -> dict[str, Any]:
    """Normalize the search settings retained with each expansion."""
    defaults = dict(
        library_mode="full",
        top_k=5,
        include_l0=True,
        use_context=True,
        diversify=True,
        use_precursor_realism=True,
        use_forward_validation=True,
    )
    if set(value) - set(defaults):
        raise ValueError("Unknown planning search setting")
    result = {**defaults, **value}
    if result["library_mode"] not in ("full", "compact"):
        raise ValueError("Unknown operator library")
    if type(result["top_k"]) is not int or not 1 <= result["top_k"] <= 50:
        raise ValueError("Top strategies must be between 1 and 50")
    for key in defaults.keys() - {"library_mode", "top_k"}:
        if type(result[key]) is not bool:
            raise ValueError("Search switches must be booleans")
    return result


@dataclass(frozen=True)
class PlanningChoice:
    """A concrete realization within a retained single-step search."""

    search_id: str
    strategy_index: int
    realization_index: int

    @property
    def choice_id(self) -> str:
        """Stable reference for conditions and occurrence identities."""
        return f"{self.search_id}:{self.strategy_index}:{self.realization_index}"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible reference."""
        return dict(
            search_id=self.search_id,
            strategy_index=self.strategy_index,
            realization_index=self.realization_index,
        )


@dataclass(frozen=True)
class PlanningNode:
    """One molecule occurrence with an optional selected precursor set."""

    node_id: str
    smiles: str
    stopped: bool = False
    choice: PlanningChoice | None = None
    children: tuple[PlanningNode, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Serialize this molecule occurrence and its selected branch."""
        return dict(
            node_id=self.node_id,
            smiles=self.smiles,
            stopped=self.stopped,
            choice=self.choice.to_dict() if self.choice else None,
            children=[child.to_dict() for child in self.children],
        )


@dataclass(frozen=True)
class PlanningSearch:
    """Full reported search evidence, separate from route selection."""

    search_id: str
    target_smiles: str
    settings: Mapping[str, Any]
    result: Mapping[str, Any]

    def to_dict(self) -> dict[str, Any]:
        """Serialize settings and the unmodified search evidence."""
        return dict(
            search_id=self.search_id,
            target_smiles=self.target_smiles,
            settings=dict(self.settings),
            result=dict(self.result),
        )


@dataclass(frozen=True)
class PlanningSession:
    """A resumable selected route, alternatives, conditions, and edit history."""

    root: PlanningNode
    past: tuple[PlanningNode, ...] = ()
    future: tuple[PlanningNode, ...] = ()
    searches: tuple[PlanningSearch, ...] = ()
    conditions: Mapping[str, Any] | None = None
    stock: Mapping[str, Any] | None = None
    schema_version: str = SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Export a session without claiming independent evidence revalidation."""
        return dict(
            schema_version=self.schema_version,
            root=self.root.to_dict(),
            past=[root.to_dict() for root in self.past],
            future=[root.to_dict() for root in self.future],
            searches=[search.to_dict() for search in self.searches],
            conditions=dict(self.conditions or {}),
            stock=dict(self.stock or {}),
        )


def start_session(target_smiles: str) -> PlanningSession:
    """Start an unresolved plan from one parsed molecular graph."""
    smiles = _molecule(target_smiles)
    return PlanningSession(PlanningNode("root", smiles))


def candidate_for_choice(
    session: PlanningSession,
    choice: PlanningChoice,
) -> tuple[PlanningSearch, Mapping[str, Any], tuple[str, ...]]:
    """Resolve reported evidence and reconcile all precursor/product graphs."""
    _integer(choice.strategy_index)
    _integer(choice.realization_index)
    search = next(
        (s for s in session.searches if s.search_id == choice.search_id), None
    )
    if search is None:
        raise ValueError("Selected search is missing")
    try:
        strategy = search.result["strategies"][choice.strategy_index]
        variants = [strategy["representative"], *strategy["alternate_realizations"]]
        candidate = _object(variants[choice.realization_index])
    except (IndexError, KeyError, TypeError) as exc:
        raise ValueError("Selected precursor choice is missing") from exc
    if candidate.get("forward_validation_status") != "verified_signature":
        raise ValueError("Only signature-verified search proposals may be selected")
    for field in ("score", "independent_reference_support"):
        value = candidate.get(field, 0)
        if type(value) not in (int, float) or not math.isfinite(value):
            raise ValueError("Saved candidate metrics must be finite numbers")
    for field in ("abstraction_level", "transformation_kind"):
        if candidate.get(field) is not None and not isinstance(candidate[field], str):
            raise ValueError("Saved candidate annotations must be strings")
    if candidate.get("forward_assessment") is not None:
        _object(candidate["forward_assessment"])
    target = _molecule(candidate.get("target_smiles"))
    if target != search.target_smiles:
        raise ValueError("Candidate target contradicts its search target")
    mol = Chem.MolFromSmiles(_text(candidate.get("precursor_smiles")))
    if mol is None or not mol.GetNumAtoms():
        raise ValueError("Invalid precursor structures")
    precursors = tuple(
        sorted(
            _molecule(Chem.MolToSmiles(fragment, isomericSmiles=True))
            for fragment in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        )
    )
    reaction = _text(candidate.get("proposed_reaction_smiles")).split(">")
    if len(reaction) != 3 or _molecule(reaction[2]) != target:
        raise ValueError("Reaction product contradicts the candidate target")
    if canonical_smiles(reaction[0]) != ".".join(precursors):
        raise ValueError("Reaction reactants contradict the complete precursor set")
    return search, candidate, precursors


def find_node(root: PlanningNode, node_id: str) -> PlanningNode:
    """Find one occurrence without conflating identical molecules."""
    if root.node_id == node_id:
        return root
    for child in root.children:
        try:
            return find_node(child, node_id)
        except LookupError:
            pass
    raise LookupError("Selected molecule is no longer in this route")


def _replace_node(
    root: PlanningNode, node_id: str, replacement: PlanningNode
) -> PlanningNode:
    if root.node_id == node_id:
        return replacement
    return replace(
        root,
        children=tuple(
            _replace_node(child, node_id, replacement) for child in root.children
        ),
    )


def _validate_root(session: PlanningSession, root: PlanningNode) -> None:
    seen: set[str] = set()

    def visit(node: PlanningNode, ancestors: tuple[str, ...]) -> None:
        if len(seen) >= MAX_NODES or len(ancestors) > MAX_DEPTH:
            raise ValueError("Plan exceeds 200 molecules or 20 reaction levels")
        if node.node_id in seen:
            raise ValueError("Duplicate molecule occurrence ID")
        seen.add(node.node_id)
        if node.smiles in ancestors:
            raise ValueError("This step creates a cycle back to an ancestor molecule")
        if node.choice:
            if node.stopped:
                raise ValueError(
                    "An expanded molecule cannot also be a starting material"
                )
            search, _, precursors = candidate_for_choice(session, node.choice)
            if search.target_smiles != node.smiles:
                raise ValueError("Chosen reaction does not produce this molecule")
            if tuple(child.smiles for child in node.children) != precursors:
                raise ValueError("A chosen step must retain every precursor occurrence")
        elif node.children:
            raise ValueError("Children require a chosen reaction")
        for child in node.children:
            visit(child, (*ancestors, node.smiles))

    visit(root, ())


def restore_session(value: Mapping[str, Any]) -> PlanningSession:
    """Validate saved structures and topology, retaining evidence as reported."""
    if value.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported planning session version")
    searches = value.get("searches", [])
    past, future = value.get("past", []), value.get("future", [])
    if not all(isinstance(items, list) for items in (searches, past, future)):
        raise ValueError("Searches and history must be lists")
    if len(searches) > MAX_SEARCHES or max(len(past), len(future)) > MAX_HISTORY:
        raise ValueError("Saved session exceeds search or undo history limits")
    records = []
    for item in searches:
        item = _object(item)
        result = _object(item.get("result"))
        _messages(result.get("warnings", []))
        if not isinstance(result.get("strategies"), list):
            raise ValueError("Saved search has no strategy list")
        for strategy in result["strategies"]:
            strategy = _object(strategy)
            _object(strategy.get("representative"))
            if not isinstance(strategy.get("alternate_realizations"), list):
                raise ValueError("Saved strategy has no precursor choices")
            for candidate in strategy["alternate_realizations"]:
                _object(candidate)
        records.append(
            PlanningSearch(
                _text(item.get("search_id")),
                _molecule(item.get("target_smiles")),
                normalize_settings(_object(item.get("settings"))),
                result,
            )
        )
    if len({s.search_id for s in records}) != len(records):
        raise ValueError("Duplicate search ID")

    def read_node(raw: Any, depth: int = 0) -> PlanningNode:
        if depth > MAX_DEPTH:
            raise ValueError("Saved route exceeds maximum depth")
        item = _object(raw)
        choice = item.get("choice")
        if choice is not None:
            choice = _object(choice)
            choice = PlanningChoice(
                _text(choice.get("search_id")),
                _integer(choice.get("strategy_index")),
                _integer(choice.get("realization_index")),
            )
        children = item.get("children")
        if not isinstance(children, list) or len(children) > MAX_NODES:
            raise ValueError("Invalid precursor occurrence list")
        if type(item.get("stopped")) is not bool:
            raise ValueError("Starting-material designation must be a boolean")
        return PlanningNode(
            _text(item.get("node_id")),
            _molecule(item.get("smiles")),
            item["stopped"],
            choice,
            tuple(read_node(c, depth + 1) for c in children),
        )

    session = PlanningSession(
        read_node(value.get("root")),
        tuple(read_node(n) for n in past),
        tuple(read_node(n) for n in future),
        tuple(records),
        _object(value.get("conditions", {})),
        _object(value.get("stock", {})),
    )
    for evidence in (session.conditions or {}).values():
        evidence = _object(evidence)
        _text(evidence.get("status"))
        _messages(evidence.get("warnings", []))
        recipes = evidence.get("recommendations", [])
        if not isinstance(recipes, list):
            raise ValueError("Saved condition recommendations must be a list")
        for recommendation in recipes:
            recipe = _object(_object(recommendation).get("resolved_recipe"))
            for components in recipe.values():
                if isinstance(components, list):
                    for component in components:
                        # Preserve other list metadata; reject null recipe components.
                        if component is None:
                            raise ValueError(
                                "Saved recipe contains a missing component"
                            )
    for evidence in (session.stock or {}).values():
        evidence = _object(evidence)
        _text(evidence.get("status"))
    for search in session.searches:
        strategies = search.result["strategies"]
        if len(strategies) > 50:
            raise ValueError("Saved search exceeds the strategy limit")
        for index, strategy in enumerate(strategies):
            if len(strategy["alternate_realizations"]) > 100:
                raise ValueError("Saved strategy exceeds the precursor-choice limit")
            for variant in range(1 + len(strategy["alternate_realizations"])):
                candidate_for_choice(
                    session, PlanningChoice(search.search_id, index, variant)
                )
    for root in (session.root, *session.past, *session.future):
        if root.node_id != "root" or root.smiles != session.root.smiles:
            raise ValueError("Edit history must preserve the target molecule")
        _validate_root(session, root)
    return session


def record_search(
    session: PlanningSession,
    node_id: str,
    settings: Mapping[str, Any],
    result: Mapping[str, Any],
) -> PlanningSession:
    """Retain a search without treating empty results as route completion."""
    node = find_node(session.root, node_id)
    settings = normalize_settings(settings)
    if _molecule(result.get("target_smiles")) != node.smiles:
        raise ValueError("Search response target does not match the selected molecule")
    search_id = digest(
        "PLANSEARCH1",
        node.smiles,
        json.dumps(settings, sort_keys=True),
        json.dumps(result, sort_keys=True),
    )
    records = tuple(s for s in session.searches if s.search_id != search_id)
    if len(records) >= MAX_SEARCHES:
        raise ValueError("Session has 100 searches; export it and start a new plan")
    search = PlanningSearch(search_id, node.smiles, settings, result)
    return restore_session(replace(session, searches=(*records, search)).to_dict())


def edit_session(
    session: PlanningSession,
    action: str,
    node_id: str = "root",
    choice: PlanningChoice | None = None,
) -> PlanningSession:
    """Choose, remove, stop, reopen, undo, or redo a branch atomically."""
    if action == "undo":
        if not session.past:
            raise ValueError("Nothing to undo")
        return replace(
            session,
            root=session.past[-1],
            past=session.past[:-1],
            future=(*session.future, session.root)[-MAX_HISTORY:],
        )
    if action == "redo":
        if not session.future:
            raise ValueError("Nothing to redo")
        return replace(
            session,
            root=session.future[-1],
            future=session.future[:-1],
            past=(*session.past, session.root)[-MAX_HISTORY:],
        )
    node = find_node(session.root, node_id)
    if action == "select" and choice is not None:
        if node.choice == choice:
            return session
        _, _, precursors = candidate_for_choice(session, choice)
        children = tuple(
            PlanningNode(
                digest("PLANNODE1", node_id, choice.choice_id, str(index)),
                smiles,
            )
            for index, smiles in enumerate(precursors)
        )
        updated = replace(node, choice=choice, children=children, stopped=False)
    elif action == "remove":
        updated = replace(node, choice=None, children=(), stopped=False)
    elif action in ("stop", "reopen"):
        if node.choice:
            raise ValueError(
                "Remove the chosen step before declaring a starting material"
            )
        updated = replace(node, stopped=action == "stop")
    else:
        raise ValueError("Unknown planning edit")
    root = _replace_node(session.root, node_id, updated)
    if root == session.root:
        return session
    _validate_root(session, root)
    return replace(
        session, root=root, past=(*session.past, session.root)[-MAX_HISTORY:], future=()
    )


def record_conditions(
    session: PlanningSession,
    node_id: str,
    evidence: Mapping[str, Any],
) -> PlanningSession:
    """Attach reported conditions to the exact chosen reaction, including history."""
    node = find_node(session.root, node_id)
    if node.choice is None:
        raise ValueError("Choose a step before looking up conditions")
    return replace(
        session,
        conditions={
            **(session.conditions or {}),
            node.choice.choice_id: dict(evidence),
        },
    )


def record_stock(
    session: PlanningSession,
    node_id: str,
    evidence: Mapping[str, Any],
) -> PlanningSession:
    """Retain an exact stock lookup separately from the user's stopping choice."""
    node = find_node(session.root, node_id)
    if _molecule(evidence.get("smiles")) != node.smiles:
        raise ValueError("Stock evidence belongs to another molecule")
    return replace(
        session, stock={**(session.stock or {}), node.smiles: dict(evidence)}
    )


def export_route(session: PlanningSession) -> ReactionRouteTree:
    """Export the selected plan through the canonical predicted-route contract."""
    _validate_root(session, session.root)
    tokens: list[str] = []
    maximum_depth = 0

    def visit(node: PlanningNode, depth: int) -> MoleculeOccurrenceNode:
        nonlocal maximum_depth
        maximum_depth = max(maximum_depth, depth)
        reaction = None
        if node.choice:
            _, candidate, _ = candidate_for_choice(session, node.choice)
            reaction_smiles = str(candidate["proposed_reaction_smiles"])
            tokens.append(f"reaction:{depth + 1}:{reaction_smiles}")
            reaction = RouteReactionNode(
                reaction_node_id=digest(
                    "PLANREACTION1", node.node_id, node.choice.choice_id
                ),
                step_id=f"step:{node.node_id}",
                depth=depth + 1,
                reaction_smiles=reaction_smiles,
                evidence=RouteStepEvidence(
                    evidence_kind="predicted",
                    source_dataset_id="core_retrosynthesis",
                    connectivity_method="user_selected_search_proposal",
                    warnings=(
                        "Search validation is retained as reported; session topology checks do not revalidate chemistry.",
                    ),
                ),
                planned_action=PlannedRouteAction(
                    operator_id=str(candidate.get("operator_id") or ""),
                    disconnection_site_key=str(
                        candidate.get("disconnection_site_key") or ""
                    ),
                    template_id=str(candidate.get("template_id") or ""),
                ),
                children=tuple(visit(child, depth + 1) for child in node.children),
            )
        return MoleculeOccurrenceNode(
            occurrence_id=node.node_id,
            smiles=node.smiles,
            depth=depth,
            terminal=node.stopped,
            terminal_evidence="user_designated_starting_material"
            if node.stopped
            else "expanded"
            if reaction
            else "unresolved",
            unresolved_reason=None if node.stopped or reaction else "not_expanded",
            reaction=reaction,
        )

    root = visit(session.root, 0)

    def chemistry_shape(node: MoleculeOccurrenceNode) -> dict[str, Any]:
        reaction = node.reaction
        return dict(
            smiles=node.smiles,
            terminal=node.terminal,
            reaction=reaction.reaction_smiles if reaction else None,
            operator=reaction.operator_id if reaction else None,
            children=[chemistry_shape(child) for child in reaction.children]
            if reaction
            else [],
        )

    tree = ReactionRouteTree(
        tree_id=digest("PLANROUTE1", json.dumps(chemistry_shape(root), sort_keys=True)),
        route_kind="planned",
        target_smiles=root.smiles,
        root=root,
        reaction_count=len(tokens),
        maximum_depth=maximum_depth,
        fingerprint_tokens=tuple(sorted(tokens)),
        connectivity_method="user_selected_search_proposals",
        warnings=(
            "User-designated starting materials are not verified stock matches.",
            "Chosen steps are predictions; retain the session export for full search and condition evidence.",
        ),
    )
    assert_valid_route_tree(tree)
    return tree


def session_response(session: PlanningSession) -> dict[str, Any]:
    """Return the resumable session, selected route, and honest leaf counts."""
    from .route_contract import iter_molecule_occurrences

    if len(json.dumps(session.to_dict())) > 18_000_000:
        raise ValueError(
            "Session is full. Export it before starting a new plan; this action was not saved."
        )
    tree = export_route(session)
    leaves = [node for node in iter_molecule_occurrences(tree) if node.reaction is None]
    return dict(
        session=session.to_dict(),
        route_tree=tree.to_dict(),
        summary=dict(
            reaction_count=tree.reaction_count,
            unresolved_count=sum(not n.terminal for n in leaves),
            starting_material_count=sum(n.terminal for n in leaves),
            maximum_depth=tree.maximum_depth,
        ),
    )
