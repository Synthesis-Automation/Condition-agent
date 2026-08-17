"""Type-agnostic hierarchical retrieval over generic reaction signatures."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterator, Literal, Mapping, Tuple

from reactive_taxonomy import departing_fragment_tokens

from .compatibility import CompatibilityAssessment, filter_compatible_precedents
from .core_retrieval import reaction_core_query_eligible
from .edit_prototypes import (
    anonymous_edit_center_compatible,
    anonymous_edit_prototype,
    anonymous_edit_similarity,
)
from .generic_indexing import GenericIndexedReaction, GenericReactionIndex
from .models import RetrievalLevelTrace
from .signature_features import (
    environment_profile_similarity,
    environment_tokens,
)
from .reaction_facets import (
    active_atom_states_compatible,
    load_reaction_facet_rules,
    reaction_facet_keys,
)
from .similarity import generic_signature_similarity, reaction_scope
from .support import summarize_evidence_support

_RULES_PATH = Path(__file__).with_name("definitions") / "generic_retrieval.v1.json"
_SUPPORTED_RETRIEVAL_LEVELS = {
    "exact_signature",
    "handle_signature",
    "named_family",
    "transformation_signature",
    "environment_neighbors",
    "bond_edit_signature",
    "reaction_core_exact",
    "reaction_core_local",
    "reaction_core_context",
    "edit_graph_neighbors",
}
RetrievalStrategy = Literal[
    "hybrid",
    "family_only",
    "generic_only",
    "transformation_prior",
    "legacy_pilot",
]
_EVALUATION_STRATEGIES = {
    "family_only",
    "generic_only",
    "transformation_prior",
    "legacy_pilot",
}
RETROSYNTHESIS_PRECEDENT_BRIDGE_VERSION = "1.0"


@lru_cache(maxsize=1)
def load_generic_retrieval_rules() -> Dict[str, Any]:
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "1.8":
        raise ValueError("unsupported generic retrieval definition schema")
    if str(rules.get("definition_id") or "") != "generic_retrieval.v1":
        raise ValueError("unexpected generic retrieval definition ID")
    if not str(rules.get("calibration_status") or "").strip():
        raise ValueError("generic retrieval definition requires calibration status")
    ladder = tuple(str(value) for value in rules.get("retrieval_ladder") or ())
    if (
        not ladder
        or len(set(ladder)) != len(ladder)
        or set(ladder) != _SUPPORTED_RETRIEVAL_LEVELS
    ):
        raise ValueError("generic retrieval ladder is incomplete or invalid")
    evaluation_ladders = rules.get("evaluation_ladders")
    if not isinstance(evaluation_ladders, Mapping) or set(
        evaluation_ladders
    ) != _EVALUATION_STRATEGIES:
        raise ValueError("generic retrieval evaluation ladders are incomplete")
    for strategy, values in evaluation_ladders.items():
        strategy_ladder = tuple(str(value) for value in values or ())
        if (
            not strategy_ladder
            or len(set(strategy_ladder)) != len(strategy_ladder)
            or not set(strategy_ladder) <= _SUPPORTED_RETRIEVAL_LEVELS
        ):
            raise ValueError(
                f"generic retrieval evaluation ladder is invalid: {strategy}"
            )
    if int(rules["environment_neighbor_limit"]) < 1:
        raise ValueError("environment_neighbor_limit must be positive")
    environment_threshold = float(rules["environment_neighbor_min_similarity"])
    if not 0.0 < environment_threshold <= 1.0:
        raise ValueError(
            "environment_neighbor_min_similarity must be in (0, 1]"
        )
    if int(rules["edit_graph_neighbor_limit"]) < 1:
        raise ValueError("edit_graph_neighbor_limit must be positive")
    edit_threshold = float(rules["edit_graph_neighbor_min_similarity"])
    if not 0.0 < edit_threshold <= 1.0:
        raise ValueError(
            "edit_graph_neighbor_min_similarity must be in (0, 1]"
        )
    return rules


def _positions(mapping: Mapping[str, Tuple[int, ...]], key: Any) -> set[int]:
    return set(mapping.get(str(key or ""), ()))


def _compatible_edit_positions(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    reaction_core: Mapping[str, Any] | None = None,
    query_reaction_smiles: str = "",
) -> set[int]:
    """Enforce net-edit and reactive-center compatibility before similarity."""
    positions = _positions(
        index.bond_edits,
        signature.get("bond_edit_signature_key"),
    )
    return _structurally_compatible_positions(
        signature,
        index,
        positions,
        reaction_core=reaction_core,
        query_reaction_smiles=query_reaction_smiles,
    )


def _structurally_compatible_positions(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    positions: set[int],
    *,
    reaction_core: Mapping[str, Any] | None = None,
    query_reaction_smiles: str = "",
) -> set[int]:
    """Filter tiers on contradictory center and departing-fragment graphs."""
    if not positions:
        return positions
    query = anonymous_edit_prototype(signature)
    ordered_positions = tuple(sorted(positions))
    rows = index.select(ordered_positions)
    compatible = set()
    for position, row in zip(ordered_positions, rows):
        candidate = anonymous_edit_prototype(row.signature)
        edit_compatible = (
            query is None
            or candidate is None
            or anonymous_edit_center_compatible(query, candidate)
        )
        if edit_compatible and active_atom_states_compatible(
            reaction_core,
            row.reaction_core,
        ):
            compatible.add(position)
    query_departing = (
        departing_fragment_tokens(query_reaction_smiles, signature)
        if query_reaction_smiles
        else ()
    )
    if not query_departing or not compatible:
        return compatible
    known = set(index.departing_fragments.get("DF1:KNOWN", ()))
    matching = {
        position
        for token in query_departing
        for position in index.departing_fragments.get(token, ())
    }
    return compatible - known | compatible & matching


def _candidate_levels(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
    query_reaction_smiles: str = "",
) -> Iterator[tuple[str, set[int]]]:
    """Yield retrieval tiers in order without computing later fallbacks early."""
    compatible = _compatible_edit_positions(
        signature,
        index,
        reaction_core=reaction_core,
        query_reaction_smiles=query_reaction_smiles,
    )
    rules = load_generic_retrieval_rules()
    ladder = (
        rules["retrieval_ladder"]
        if strategy == "hybrid"
        else (rules["evaluation_ladders"] or {}).get(strategy)
    )
    if not ladder:
        raise ValueError(f"Unsupported generic retrieval strategy: {strategy}")
    core_eligible = bool(
        reaction_core
        and reaction_core_query_eligible(reaction_core, index)[0]
    )
    for raw_level in ladder:
        level = str(raw_level)
        if level == "exact_signature":
            positions = _positions(
                index.exact, signature.get("exact_signature_key")
            ) & compatible
        elif level == "handle_signature":
            positions = _positions(
                index.handles, signature.get("handle_signature_key")
            ) & compatible
        elif level == "named_family":
            family = str(signature.get("named_family") or "")
            family_confidence = float(
                signature.get("family_confidence") or 0.0
            )
            positions = (
                _positions(index.families, family) & compatible
                if family
                and family_confidence
                >= float(rules["high_confidence_family_threshold"])
                else set()
            )
        elif level == "transformation_signature":
            positions = _positions(
                index.transformations,
                signature.get("transformation_signature_key"),
            ) & compatible
        elif level == "environment_neighbors":
            positions = _environment_neighbor_positions(
                signature, index, compatible
            )
        elif level == "bond_edit_signature":
            positions = compatible
        elif level.startswith("reaction_core_"):
            positions = _structurally_compatible_positions(
                signature,
                index,
                _core_level_positions(
                    reaction_core,
                    index,
                    level=level,
                    query_eligible=core_eligible,
                ),
                reaction_core=reaction_core,
                query_reaction_smiles=query_reaction_smiles,
            )
        elif level == "edit_graph_neighbors":
            positions = _structurally_compatible_positions(
                signature,
                index,
                _edit_graph_neighbor_positions(
                    signature,
                    index,
                    exclude=compatible,
                    query_reaction_smiles=query_reaction_smiles,
                ),
                reaction_core=reaction_core,
                query_reaction_smiles=query_reaction_smiles,
            )
        else:  # pragma: no cover - validated definitions prevent this branch.
            raise ValueError(f"Unsupported generic retrieval level: {level}")
        yield level, positions


def _core_level_positions(
    reaction_core: Mapping[str, Any] | None,
    index: GenericReactionIndex,
    *,
    level: str,
    query_eligible: bool | None = None,
) -> set[int]:
    """Return one bulk-loaded, event-compatible reaction-core tier."""
    if not reaction_core or not (
        reaction_core_query_eligible(reaction_core, index)[0]
        if query_eligible is None
        else query_eligible
    ):
        return set()
    event_count = int(reaction_core.get("event_count") or 0)
    definitions = {
        "reaction_core_exact": (
            index.core_exact,
            str(reaction_core.get("exact_core_key") or ""),
        ),
        "reaction_core_local": (
            index.core_typed,
            str(reaction_core.get("typed_core_key") or ""),
        ),
        "reaction_core_context": (
            index.core_shapes,
            str(reaction_core.get("shape_core_key") or ""),
        ),
    }
    if level not in definitions:
        raise ValueError(f"Unsupported reaction-core retrieval level: {level}")
    mapping, key = definitions[level]
    positions = tuple(mapping.get(key, ()))
    rows = index.select(positions)
    return {
        position
        for position, row in zip(positions, rows)
        if row.signature
        and row.reaction_core
        and int(row.reaction_core.get("event_count") or 0) == event_count
    }


def _edit_graph_neighbor_positions(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    exclude: set[int],
    query_reaction_smiles: str = "",
) -> set[int]:
    """Find chemistry-gated approximate edit graphs without family routing."""
    query = anonymous_edit_prototype(signature)
    if query is None:
        return set()
    query = _project_departing_ring_inventory(
        query,
        query_reaction_smiles,
        signature,
    )
    rules = load_generic_retrieval_rules()
    threshold = float(rules["edit_graph_neighbor_min_similarity"])
    limit = int(rules["edit_graph_neighbor_limit"])
    scored = []
    candidate_positions = getattr(
        index.rows,
        "edit_graph_candidate_positions",
        None,
    )
    if callable(candidate_positions):
        positions = tuple(candidate_positions(query))
    else:
        positions = tuple(range(len(index.rows)))
    included_positions = tuple(
        position for position in positions if position not in exclude
    )
    rows = index.select(included_positions)
    for position, row in zip(included_positions, rows):
        if not row.signature:
            continue
        candidate = anonymous_edit_prototype(row.signature)
        if candidate is None:
            continue
        candidate = _project_departing_ring_inventory(
            candidate,
            row.reaction_smiles,
            row.signature,
        )
        score = anonymous_edit_similarity(query, candidate)
        if score >= threshold:
            scored.append(
                (
                    score,
                    row.canonical_reaction_id,
                    row.reaction_id,
                    row.observation_id,
                    position,
                )
            )
    scored.sort(key=lambda item: (-item[0],) + item[1:])
    return {item[-1] for item in scored[:limit]}


def _project_departing_ring_inventory(
    prototype: Any,
    reaction_smiles: str,
    signature: Mapping[str, Any],
) -> Any:
    """Remove cycles belonging wholly to a structurally observed departure."""

    if prototype.ring_count_delta >= 0 or not reaction_smiles:
        return prototype
    ring_count = _departing_ring_count(reaction_smiles, signature)
    if not ring_count:
        return prototype
    return replace(
        prototype,
        ring_count_delta=min(0, prototype.ring_count_delta + ring_count),
    )


def _departing_ring_count(
    reaction_smiles: str,
    signature: Mapping[str, Any],
) -> int:
    """Return cycles wholly contained in observed departing fragments."""

    ring_count = 0
    for token in departing_fragment_tokens(reaction_smiles, signature):
        parts = token.split(":", 3)
        if len(parts) >= 3 and parts[1].startswith("R"):
            ring_count += int(parts[1][1:] or 0)
    return ring_count


def _environment_neighbor_positions(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    compatible: set[int],
) -> set[int]:
    """Narrow edit-compatible rows using interpretable local environments."""
    query_tokens = environment_tokens(signature)
    if not query_tokens:
        return set()
    rules = load_generic_retrieval_rules()
    threshold = float(rules["environment_neighbor_min_similarity"])
    limit = int(rules["environment_neighbor_limit"])
    token_positions = set()
    for token in query_tokens:
        token_positions.update(index.environment_features.get(token, ()))
    positions = tuple(sorted(compatible & token_positions))
    rows = index.select(positions)
    scored = []
    for position, row in zip(positions, rows):
        score = environment_profile_similarity(
            signature,
            row.signature,
        )
        if score >= threshold:
            scored.append(
                (
                    score,
                    row.canonical_reaction_id,
                    row.reaction_id,
                    row.observation_id,
                    position,
                )
            )
    scored.sort(
        key=lambda item: (-item[0],) + item[1:]
    )
    return {item[-1] for item in scored[:limit]}


def _minimum_support(minimum_pool_size: int | None) -> int:
    rules = load_generic_retrieval_rules()
    minimum = (
        int(minimum_pool_size)
        if minimum_pool_size is not None
        else int(rules["minimum_pool_size"])
    )
    if minimum < 1:
        raise ValueError("minimum_pool_size must be positive")
    return minimum


def _no_compatible_bond_edit_trace(minimum: int) -> RetrievalLevelTrace:
    """Return the stable trace used when every chemistry tier is empty."""
    return RetrievalLevelTrace(
        level="bond_edit_gate",
        candidate_count=0,
        independent_candidate_count=0,
        compatible_candidate_count=0,
        independent_compatible_candidate_count=0,
        excluded_candidate_count=0,
        minimum_independent_support=minimum,
        status="no_compatible_bond_edit",
    )


def retrieve_generic_pool_with_trace(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
) -> tuple[
    str,
    Tuple[GenericIndexedReaction, ...],
    Tuple[RetrievalLevelTrace, ...],
]:
    """Select by independent support and retain every attempted tier."""
    minimum = _minimum_support(minimum_pool_size)
    levels = _candidate_levels(
        signature,
        index,
        reaction_core=reaction_core,
        strategy=strategy,
    )
    fallback: tuple[
        str,
        Tuple[GenericIndexedReaction, ...],
        int,
        int,
    ] | None = None
    traces = []
    saw_candidates = False
    for level, positions in levels:
        if not positions:
            traces.append(
                RetrievalLevelTrace(
                    level=level,
                    candidate_count=0,
                    independent_candidate_count=0,
                    compatible_candidate_count=0,
                    independent_compatible_candidate_count=0,
                    excluded_candidate_count=0,
                    minimum_independent_support=minimum,
                    status="empty",
                )
            )
            continue
        saw_candidates = True
        rows = index.select(sorted(positions))
        support = summarize_evidence_support(rows)
        trace = RetrievalLevelTrace(
            level=level,
            candidate_count=len(rows),
            independent_candidate_count=support.independent_count,
            compatible_candidate_count=len(rows),
            independent_compatible_candidate_count=support.independent_count,
            excluded_candidate_count=0,
            minimum_independent_support=minimum,
            status="insufficient_independent_support",
        )
        traces.append(trace)
        if fallback is None or (
            support.independent_count,
            len(rows),
        ) > (
            fallback[3],
            len(fallback[1]),
        ):
            fallback = (
                level,
                rows,
                len(traces) - 1,
                support.independent_count,
            )
        if support.independent_count >= minimum:
            traces[-1] = replace(trace, status="selected")
            return level, rows, tuple(traces)
    if fallback is None:
        if not saw_candidates:
            trace = _no_compatible_bond_edit_trace(minimum)
            return "no_compatible_bond_edit", (), (trace,)
        return "no_compatible_precedent", (), tuple(traces)
    level, rows, trace_index, _ = fallback
    traces[trace_index] = replace(
        traces[trace_index], status="selected_limited_support"
    )
    return (
        f"{level}_limited_support",
        rows,
        tuple(traces),
    )


def retrieve_generic_pool(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
) -> Tuple[str, Tuple[GenericIndexedReaction, ...]]:
    """Compatibility wrapper returning the historical two-value result."""
    level, rows, _ = retrieve_generic_pool_with_trace(
        signature,
        index,
        minimum_pool_size=minimum_pool_size,
        reaction_core=reaction_core,
        strategy=strategy,
    )
    return level, rows


@dataclass(frozen=True)
class CompatibleRetrievalResult:
    """Selected compatibility-filtered pool and its complete retrieval trace."""

    level: str
    pool: tuple[tuple[GenericIndexedReaction, CompatibilityAssessment], ...]
    candidate_count: int
    independent_candidate_count: int
    excluded_candidate_count: int
    independent_compatible_candidate_count: int
    trace: Tuple[RetrievalLevelTrace, ...]


def retrieve_preferred_reaction_pool_with_trace(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    reaction_ids: Tuple[str, ...],
    *,
    reaction_core: Mapping[str, Any] | None = None,
    minimum_pool_size: int | None = None,
    query_reaction_smiles: str = "",
) -> CompatibleRetrievalResult:
    """Resolve exact source IDs, then apply every structural and recipe gate."""

    minimum = _minimum_support(minimum_pool_size)
    requested = tuple(dict.fromkeys(str(value) for value in reaction_ids if value))
    raw_positions = {
        position
        for reaction_id in requested
        for position in index.reaction_ids.get(reaction_id, ())
    }
    structural_positions = _structurally_compatible_positions(
        signature,
        index,
        raw_positions,
        reaction_core=reaction_core,
        query_reaction_smiles=query_reaction_smiles,
    )
    raw_rows = index.select(sorted(raw_positions))
    structural_rows = index.select(sorted(structural_positions))
    accepted, recipe_excluded = filter_compatible_precedents(
        signature,
        structural_rows,
    )
    raw_support = summarize_evidence_support(raw_rows)
    accepted_support = summarize_evidence_support(
        tuple(row for row, _ in accepted)
    )
    excluded_count = (
        len(raw_positions) - len(structural_positions) + len(recipe_excluded)
    )
    status = (
        "selected"
        if accepted_support.independent_count >= minimum
        else "selected_limited_support"
        if accepted
        else "no_compatible_precedent"
        if raw_positions
        else "empty"
    )
    base_level = "retrosynthesis_precedent"
    level = (
        f"{base_level}_limited_support"
        if accepted and accepted_support.independent_count < minimum
        else base_level
    )
    trace = RetrievalLevelTrace(
        level=base_level,
        candidate_count=len(raw_rows),
        independent_candidate_count=raw_support.independent_count,
        compatible_candidate_count=len(accepted),
        independent_compatible_candidate_count=accepted_support.independent_count,
        excluded_candidate_count=excluded_count,
        minimum_independent_support=minimum,
        status=status,
    )
    return CompatibleRetrievalResult(
        level=level,
        pool=accepted,
        candidate_count=len(raw_rows),
        independent_candidate_count=raw_support.independent_count,
        excluded_candidate_count=excluded_count,
        independent_compatible_candidate_count=accepted_support.independent_count,
        trace=(trace,),
    )


@dataclass(frozen=True)
class ProgressiveCompatibleRetrievalResult:
    """Ordered chemistry tiers collected until enough recipe cores exist."""

    tiers: tuple[
        tuple[
            str,
            tuple[tuple[GenericIndexedReaction, CompatibilityAssessment], ...],
        ],
        ...,
    ]
    candidate_count: int
    independent_candidate_count: int
    compatible_candidate_count: int
    independent_compatible_candidate_count: int
    excluded_candidate_count: int
    trace: Tuple[RetrievalLevelTrace, ...]

    @property
    def level(self) -> str:
        if not self.tiers:
            return "no_compatible_condition_precedent"
        if len(self.tiers) == 1:
            return self.tiers[0][0]
        return "progressive_reaction_facets"


def retrieve_progressive_compatible_pools_with_trace(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    reaction_core: Mapping[str, Any] | None,
    fallback_descriptor: Mapping[str, Any] | None = None,
    target_recipe_count: int,
    minimum_pool_size: int | None = None,
    query_reaction_smiles: str = "",
) -> ProgressiveCompatibleRetrievalResult | None:
    """Collect ordered facet tiers without allowing broad rows to displace exact ones.

    ``None`` means a trustworthy graph-facet projection was unavailable and the
    historical retrieval policy should be used unchanged.
    """
    if target_recipe_count < 1:
        raise ValueError("target_recipe_count must be positive")
    keys = reaction_facet_keys(
        signature,
        reaction_core,
        fallback_descriptor,
    )
    rules = load_reaction_facet_rules()
    if not keys or any(level not in keys for level in rules["retrieval_ladder"]):
        return None
    minimum = _minimum_support(minimum_pool_size)
    lookup_maps = {
        "reaction_facet_exact": index.facet_exact,
        "reaction_facet_attachment_relaxed": index.facet_attachments,
    }
    seen_positions: set[int] = set()
    seen_recipe_cores: set[str] = set()
    raw_rows_by_position: dict[int, GenericIndexedReaction] = {}
    accepted_rows_by_position: dict[int, GenericIndexedReaction] = {}
    tiers = []
    traces = []
    excluded_total = 0

    def process_level(level: str, raw_positions: set[int]) -> bool:
        nonlocal excluded_total
        positions = raw_positions - seen_positions
        seen_positions.update(positions)
        if not positions:
            traces.append(
                RetrievalLevelTrace(
                    level=level,
                    candidate_count=0,
                    independent_candidate_count=0,
                    compatible_candidate_count=0,
                    independent_compatible_candidate_count=0,
                    excluded_candidate_count=0,
                    minimum_independent_support=minimum,
                    status="empty",
                )
            )
            return False
        ordered_positions = tuple(sorted(positions))
        rows = index.select(ordered_positions)
        raw_rows_by_position.update(zip(ordered_positions, rows))
        accepted, excluded = filter_compatible_precedents(signature, rows)
        accepted_ids = {id(row) for row, _ in accepted}
        accepted_positions = {
            position
            for position, row in zip(ordered_positions, rows)
            if id(row) in accepted_ids
        }
        accepted_rows_by_position.update(
            (position, row)
            for position, row in zip(ordered_positions, rows)
            if position in accepted_positions
        )
        excluded_total += len(excluded)
        raw_support = summarize_evidence_support(rows)
        accepted_rows = tuple(row for row, _ in accepted)
        accepted_support = summarize_evidence_support(accepted_rows)
        new_recipe_cores = {
            row.recipe_core_id for row in accepted_rows
        } - seen_recipe_cores
        if accepted and new_recipe_cores:
            tiers.append((level, accepted))
            seen_recipe_cores.update(new_recipe_cores)
            status = (
                "selected_target_reached"
                if len(seen_recipe_cores) >= target_recipe_count
                else "selected_progressive"
                if accepted_support.independent_count >= minimum
                else "selected_limited_support"
            )
        elif accepted:
            status = "no_new_recipe_core"
        else:
            status = "no_compatible_recipe"
        traces.append(
            RetrievalLevelTrace(
                level=level,
                candidate_count=len(rows),
                independent_candidate_count=raw_support.independent_count,
                compatible_candidate_count=len(accepted),
                independent_compatible_candidate_count=(
                    accepted_support.independent_count
                ),
                excluded_candidate_count=len(excluded),
                minimum_independent_support=minimum,
                status=status,
            )
        )
        return len(seen_recipe_cores) >= target_recipe_count

    target_reached = False
    for level in rules["retrieval_ladder"]:
        target_reached = process_level(
            level,
            _structurally_compatible_positions(
                signature,
                index,
                _positions(lookup_maps[level], keys[level]),
                reaction_core=reaction_core,
                query_reaction_smiles=query_reaction_smiles,
            ),
        )
        if target_reached:
            break
    if not target_reached:
        for level, positions in _candidate_levels(
            signature,
            index,
            reaction_core=reaction_core,
            strategy="hybrid",
            query_reaction_smiles=query_reaction_smiles,
        ):
            target_reached = process_level(level, positions)
            if target_reached:
                break
    raw_rows = tuple(raw_rows_by_position.values())
    compatible_rows = tuple(accepted_rows_by_position.values())
    raw_support = summarize_evidence_support(raw_rows)
    compatible_support = summarize_evidence_support(compatible_rows)
    return ProgressiveCompatibleRetrievalResult(
        tiers=tuple(tiers),
        candidate_count=len(raw_rows),
        independent_candidate_count=raw_support.independent_count,
        compatible_candidate_count=len(compatible_rows),
        independent_compatible_candidate_count=compatible_support.independent_count,
        excluded_candidate_count=excluded_total,
        trace=tuple(traces),
    )


def retrieve_compatible_generic_pool_with_trace(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
    query_reaction_smiles: str = "",
) -> CompatibleRetrievalResult:
    """Apply compatibility before independent-support checks at every tier."""
    minimum = _minimum_support(minimum_pool_size)
    levels = _candidate_levels(
        signature,
        index,
        reaction_core=reaction_core,
        strategy=strategy,
        query_reaction_smiles=query_reaction_smiles,
    )
    fallback = None
    traces = []
    saw_candidates = False
    for level, positions in levels:
        if not positions:
            traces.append(
                RetrievalLevelTrace(
                    level=level,
                    candidate_count=0,
                    independent_candidate_count=0,
                    compatible_candidate_count=0,
                    independent_compatible_candidate_count=0,
                    excluded_candidate_count=0,
                    minimum_independent_support=minimum,
                    status="empty",
                )
            )
            continue
        saw_candidates = True
        rows = index.select(sorted(positions))
        accepted, excluded = filter_compatible_precedents(signature, rows)
        raw_support = summarize_evidence_support(rows)
        accepted_rows = tuple(row for row, _ in accepted)
        accepted_support = summarize_evidence_support(accepted_rows)
        status = (
            "insufficient_independent_support"
            if accepted
            else "no_compatible_recipe"
        )
        trace = RetrievalLevelTrace(
            level=level,
            candidate_count=len(rows),
            independent_candidate_count=raw_support.independent_count,
            compatible_candidate_count=len(accepted),
            independent_compatible_candidate_count=(
                accepted_support.independent_count
            ),
            excluded_candidate_count=len(excluded),
            minimum_independent_support=minimum,
            status=status,
        )
        traces.append(trace)
        if not accepted:
            continue
        candidate = (
            level,
            accepted,
            len(rows),
            raw_support.independent_count,
            len(excluded),
            accepted_support.independent_count,
            len(traces) - 1,
        )
        if fallback is None or (
            accepted_support.independent_count,
            len(accepted),
        ) > (
            fallback[5],
            len(fallback[1]),
        ):
            fallback = candidate
        if accepted_support.independent_count >= minimum:
            traces[-1] = replace(trace, status="selected")
            return CompatibleRetrievalResult(
                level=level,
                pool=accepted,
                candidate_count=len(rows),
                independent_candidate_count=raw_support.independent_count,
                excluded_candidate_count=len(excluded),
                independent_compatible_candidate_count=(
                    accepted_support.independent_count
                ),
                trace=tuple(traces),
            )
    if fallback is None:
        if not saw_candidates:
            trace = _no_compatible_bond_edit_trace(minimum)
            return CompatibleRetrievalResult(
                "no_compatible_bond_edit", (), 0, 0, 0, 0, (trace,)
            )
        bond_rows = index.select(
            sorted(
                _compatible_edit_positions(
                    signature,
                    index,
                    reaction_core=reaction_core,
                    query_reaction_smiles=query_reaction_smiles,
                )
            )
        )
        _, excluded = filter_compatible_precedents(signature, bond_rows)
        support = summarize_evidence_support(bond_rows)
        return CompatibleRetrievalResult(
            "no_compatible_condition_precedent",
            (),
            len(bond_rows),
            support.independent_count,
            len(excluded),
            0,
            tuple(traces),
        )
    (
        level,
        accepted,
        raw_count,
        independent_raw_count,
        excluded_count,
        independent_accepted_count,
        trace_index,
    ) = fallback
    traces[trace_index] = replace(
        traces[trace_index], status="selected_limited_support"
    )
    return CompatibleRetrievalResult(
        f"{level}_limited_support",
        accepted,
        raw_count,
        independent_raw_count,
        excluded_count,
        independent_accepted_count,
        tuple(traces),
    )


def retrieve_compatible_generic_pool(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
    query_reaction_smiles: str = "",
) -> tuple[
    str,
    tuple[tuple[GenericIndexedReaction, CompatibilityAssessment], ...],
    int,
    int,
]:
    """Compatibility wrapper returning the historical four-value result."""
    result = retrieve_compatible_generic_pool_with_trace(
        signature,
        index,
        minimum_pool_size=minimum_pool_size,
        reaction_core=reaction_core,
        strategy=strategy,
        query_reaction_smiles=query_reaction_smiles,
    )
    return (
        result.level,
        result.pool,
        result.candidate_count,
        result.excluded_candidate_count,
    )


__all__ = [
    "generic_signature_similarity",
    "CompatibleRetrievalResult",
    "ProgressiveCompatibleRetrievalResult",
    "RetrievalStrategy",
    "RETROSYNTHESIS_PRECEDENT_BRIDGE_VERSION",
    "load_generic_retrieval_rules",
    "reaction_scope",
    "retrieve_compatible_generic_pool",
    "retrieve_compatible_generic_pool_with_trace",
    "retrieve_progressive_compatible_pools_with_trace",
    "retrieve_preferred_reaction_pool_with_trace",
    "retrieve_generic_pool",
    "retrieve_generic_pool_with_trace",
]
