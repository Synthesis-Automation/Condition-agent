"""Assess the strength of indexed precedents for a proposed reaction step.

The grade describes structural precedent under the source conditions. A valid
yield below the policy threshold excludes an observation from transferable
support. Source groups are conservative proxies for repeated reports.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Callable
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
import json
from math import isfinite
from pathlib import Path
from typing import Any

from reactive_taxonomy.fragment_search import indexed_product
from reactive_taxonomy.step_precedent_context import (
    StepContextComparison,
    StepPrecedentContext,
    build_step_context,
    compare_step_contexts,
)

from .generic_indexing import GenericReactionIndex
from .reaction_precedent_support import (
    ReactionPrecedentEvidence,
    collect_reaction_precedent_support,
)
from .shared_core_index import SharedCoreIndex


_POLICY_PATH = Path(__file__).with_name("definitions") / "step_plausibility.v2.json"
_LEVELS = ("whole_reaction", "observed_local", "retained_local", "retained_typed")
_GRADES = (
    "no_qualified_precedent", "structural_analogue", "context_supported",
    "transferable_precedent", "repeated_transferable_precedent",
)


@dataclass(frozen=True)
class StepPlausibilityPolicy:
    """Versioned ordinal rules and explicitly heuristic cutoffs."""

    definition_id: str
    definition_version: str
    score_semantics: str
    close_levels: tuple[str, ...]
    transferable_levels: tuple[str, ...]
    transferable_analogue_context: str
    minimum_partner_similarity: float
    low_yield_below_pct: float
    repeated_support_minimum_source_groups: int
    ranks: tuple[tuple[str, int | None], ...]


def load_step_plausibility_policy() -> StepPlausibilityPolicy:
    """Load and validate the experimental standard without mutating shared state."""
    raw = json.loads(_POLICY_PATH.read_text(encoding="utf-8"))
    close, transferable = tuple(raw["close_levels"]), tuple(raw["transferable_levels"])
    sim = raw["minimum_partner_similarity"]
    pct = raw["low_yield_below_pct"]
    repeats = raw["repeated_support_minimum_source_groups"]
    expected_ranks = {"unresolved": None, **{grade: i for i, grade in enumerate(_GRADES)}}
    if (
        raw["definition_id"] != "step_plausibility.v2"
        or raw["schema_version"] != "2.0"
        or not isinstance(raw["definition_version"], str) or not raw["definition_version"]
        or raw["status"] != "development_pending_independent_review"
        or raw["score_semantics"] != "ordinal_strength_of_documented_structural_support"
        or raw["transferable_analogue_context"] != "bidirectional"
        or raw["ranks"] != expected_ranks
        or any(type(raw["ranks"][grade]) is not int for grade in _GRADES)
        or not transferable or not close or transferable != close[:len(transferable)]
        or close != _LEVELS[:len(close)] or "retained_local" in transferable
        or type(sim) not in (int, float) or not 0 < sim <= 1
        or type(pct) not in (int, float) or not 0 < pct <= 100
        or type(repeats) is not int or repeats < 2
    ):
        raise ValueError("Invalid step plausibility policy")
    return StepPlausibilityPolicy(
        raw["definition_id"], raw["definition_version"], raw["score_semantics"],
        close, transferable, raw["transferable_analogue_context"], sim, pct, repeats,
        tuple(raw["ranks"].items()),
    )


@dataclass(frozen=True)
class StepPlausibilityMatch:
    """One audited source observation, its chemical context and measurement."""

    precedent: ReactionPrecedentEvidence
    context: StepContextComparison
    outcome: str
    context_covered_close: bool
    structurally_comparable: bool
    source_group_id: str | None = None
    excluded_by_low_yield: bool = False
    contributes_to_support: bool = False


@dataclass(frozen=True)
class StepSupportCounts:
    """Counts in the full visited pool; outcome counts cover comparable rows only."""

    distinct_reference_count: int
    context_covered_close_observation_count: int
    context_covered_close_reference_count: int
    comparable_observation_count: int
    supporting_observation_count: int
    supporting_reference_count: int
    supporting_distinct_reaction_count: int
    supporting_source_group_count: int
    at_or_above_threshold_yield_count: int
    low_yield_count: int
    zero_yield_count: int
    unreported_yield_count: int
    unusable_yield_count: int
    unknown_yield_count: int
    missing_reference_supporting_count: int


@dataclass(frozen=True)
class StepPlausibilityAssessment:
    """Deterministic, bounded evidence assessment for one precursor/product pair."""

    grade: str
    evidence_rank: int | None
    measured_evidence_status: str
    query_reaction_smiles: str
    strongest_level: str | None
    candidate_count: int
    qualified_count: int
    candidate_limit: int
    match_limit: int
    candidate_truncated: bool
    matches_truncated: bool
    level_counts: tuple[tuple[str, int], ...]
    query_context: StepPrecedentContext
    support_counts: StepSupportCounts
    supporting_source_groups: tuple[tuple[str, ...], ...]
    matches: tuple[StepPlausibilityMatch, ...]
    exclusions: tuple[tuple[str, int], ...]
    warnings: tuple[str, ...]
    policy: StepPlausibilityPolicy
    index_source_identity: str
    shared_core_definition_hash: str | None
    schema_version: str = "step_plausibility.v2"

    def to_dict(self) -> dict[str, Any]:
        """Return source-linked evidence and the rules that assigned its grade."""
        return asdict(self)


def _classify_outcome(row: ReactionPrecedentEvidence, threshold: float) -> str:
    value = row.yield_pct
    if value is None:
        return "unusable" if row.outcome_status in ("invalid", "conflicting", "ambiguous") else "unreported"
    if (row.outcome_status != "usable" or type(value) not in (int, float)
            or not isfinite(value) or not 0 <= value <= 100):
        return "unusable"
    if value == 0:
        return "zero_yield"
    return "at_or_above_threshold" if value >= threshold else "low_yield"


def _normalized_reaction(reaction_smiles: str) -> tuple[str, str]:
    """Canonical map-free reactants and products retaining charge, stereo and multiplicity."""
    left, _, right = reaction_smiles.split(">")
    return indexed_product(left)[0], indexed_product(right)[0]


def _source_groups(
    matches: list[StepPlausibilityMatch],
    normalized_reaction: Callable[[str], tuple[str, str]],
) -> dict[str, str]:
    """Collapse reference IDs connected by a common observation or probable copy.

    A probable copy has the same map-free stereospecific substrate and product
    graphs and the same known canonical recipe ID. Missing recipe IDs cannot
    connect different reports on chemical identity alone. Group membership is
    transitive and independent of yield magnitude or availability; it can undercount genuine
    replication and cannot resolve all patent families or shared laboratories.
    """
    refs = {m.precedent.reference_id for m in matches if m.precedent.reference_id}
    parents = {ref: ref for ref in refs}

    def representative(key: str) -> str:
        while parents[key] != key:
            parents[key] = parents[parents[key]]
            key = parents[key]
        return key

    def connect(ids: set[str]) -> None:
        roots = {representative(key) for key in ids}
        leader = min(roots)
        for root in roots:
            parents[root] = leader

    same_observation: dict[tuple[str, tuple[str, str]], set[str]] = defaultdict(set)
    experiments: dict[tuple[tuple[str, str], str], set[str]] = defaultdict(set)
    for match in matches:
        row = match.precedent
        if not row.reference_id:
            continue
        graph = normalized_reaction(row.reaction_smiles)
        if row.observation_id:
            same_observation[row.observation_id, graph].add(row.reference_id)
        recipe_id = str((row.resolved_recipe or {}).get("recipe_id") or "")
        if recipe_id and match.context.status != "unresolved":
            experiments[graph, recipe_id].add(row.reference_id)
    for refs in same_observation.values():
        connect(refs)
    for refs in experiments.values():
        connect(refs)
    return {reference: representative(reference) for reference in parents}


def _display_sample(
    matches: list[StepPlausibilityMatch], limit: int,
) -> tuple[StepPlausibilityMatch, ...]:
    """Show a matched counterexample early and avoid filling the page with one report."""
    def priority(m: StepPlausibilityMatch) -> tuple:
        row = m.precedent
        similarity = m.context.minimum_partner_similarity
        return (not m.structurally_comparable, not m.context_covered_close,
                _LEVELS.index(row.level), -(similarity if similarity is not None else -1),
                row.reference_id, row.reaction_id, row.observation_id or "", row.reaction_smiles)

    ordered = sorted(range(len(matches)), key=lambda i: priority(matches[i]))
    first_support: dict[str, int] = {}
    adverse = []
    for i in ordered:
        m = matches[i]
        if not m.structurally_comparable:
            continue
        if m.source_group_id and m.contributes_to_support:
            first_support.setdefault(m.source_group_id, i)
        elif m.outcome in ("low_yield", "zero_yield"):
            adverse.append(i)
    support = list(first_support.values())
    priority_ids = dict.fromkeys(support[:1] + adverse[:1] + support[1:] + ordered)
    return tuple(matches[i] for i in list(priority_ids)[:limit])


def _assess_structural_match(
    query: StepPrecedentContext, source: StepPrecedentContext,
    row: ReactionPrecedentEvidence, policy: StepPlausibilityPolicy,
) -> StepPlausibilityMatch:
    context = compare_step_contexts(query, source)
    covered = context.status == "covered" and row.level in policy.close_levels
    min_similarity = context.minimum_partner_similarity
    transferable_analogue = (
        min_similarity is not None and min_similarity >= policy.minimum_partner_similarity
        and context.reverse_status == "covered" and not query.has_specified_stereo
    )
    outcome = _classify_outcome(row, policy.low_yield_below_pct)
    return StepPlausibilityMatch(
        row, context, outcome, covered, covered and row.level in policy.transferable_levels
        and (row.level == "whole_reaction" or transferable_analogue),
        excluded_by_low_yield=outcome in ("low_yield", "zero_yield"),
    )


def assess_step_plausibility(
    reaction_smiles: str, index: GenericReactionIndex, shared_index: SharedCoreIndex,
    *, candidate_limit: int = 128, match_limit: int = 20,
) -> StepPlausibilityAssessment:
    """Grade graph-qualified precedents with transferable context and source attribution.

    Both indexes are the bound inputs used for condition recommendation. Every
    qualified visited row enters the assessment before the evidence display cap.
    A numeric yield describes the source experiment under its own recipe.
    """
    assessment, _ = _assess_step_and_pool(
        reaction_smiles, index, shared_index, candidate_limit=candidate_limit, match_limit=match_limit,
    )
    return assessment


def _assess_step_and_pool(
    reaction_smiles: str, index: GenericReactionIndex, shared_index: SharedCoreIndex,
    *, candidate_limit: int, match_limit: int,
) -> tuple[StepPlausibilityAssessment, list[StepPlausibilityMatch]]:
    """Return the display and the same full pool for downstream evidence views."""
    if type(match_limit) is not int or not 1 <= match_limit <= 50:
        raise ValueError("match_limit must be an integer between 1 and 50")
    policy = load_step_plausibility_policy()
    pool = collect_reaction_precedent_support(
        reaction_smiles, index, shared_index, candidate_limit=candidate_limit,
    )
    context_for = lru_cache(maxsize=None)(build_step_context)
    query = context_for(reaction_smiles)
    matches = [_assess_structural_match(query, context_for(row.reaction_smiles), row, policy)
               for row in pool.matches]
    normalized_reaction = lru_cache(maxsize=None)(_normalized_reaction)
    groups = _source_groups(matches, normalized_reaction)
    # The same observation may have a sparse copy. A reported low target yield
    # on that observation applies to its copies, without identifying distinct
    # experiments as failures solely because they used the same recipe.
    low_observations = {
        (m.precedent.observation_id, normalized_reaction(m.precedent.reaction_smiles))
        for m in matches if m.excluded_by_low_yield and m.precedent.observation_id
    }
    for i, match in enumerate(matches):
        row = match.precedent
        group = groups.get(row.reference_id)
        low = match.excluded_by_low_yield or (
            bool(low_observations) and row.observation_id is not None
            and (row.observation_id, normalized_reaction(row.reaction_smiles)) in low_observations
        )
        matches[i] = replace(
            match, source_group_id=group, excluded_by_low_yield=low,
            contributes_to_support=match.structurally_comparable and group is not None and not low,
        )
    close = [m for m in matches if m.context_covered_close]
    comparable = [m for m in close if m.structurally_comparable]
    outcomes = Counter(m.outcome for m in comparable)
    supporting = [m for m in comparable if not m.excluded_by_low_yield]
    supporting_references = {m.precedent.reference_id for m in supporting if m.contributes_to_support}
    supporting_groups: dict[str, list[str]] = defaultdict(list)
    for reference in sorted(supporting_references):
        supporting_groups[groups[reference]].append(reference)
    counts = StepSupportCounts(
        distinct_reference_count=pool.distinct_reference_count,
        context_covered_close_observation_count=len(close),
        context_covered_close_reference_count=len({
            m.precedent.reference_id for m in close if m.precedent.reference_id
        }),
        comparable_observation_count=len(comparable),
        supporting_observation_count=len(supporting),
        supporting_reference_count=len(supporting_references),
        supporting_distinct_reaction_count=len({
            normalized_reaction(m.precedent.reaction_smiles)
            for m in supporting if m.contributes_to_support
        }),
        supporting_source_group_count=len(supporting_groups),
        at_or_above_threshold_yield_count=outcomes["at_or_above_threshold"],
        low_yield_count=outcomes["low_yield"],
        zero_yield_count=outcomes["zero_yield"],
        unreported_yield_count=outcomes["unreported"],
        unusable_yield_count=outcomes["unusable"],
        unknown_yield_count=outcomes["unreported"] + outcomes["unusable"],
        missing_reference_supporting_count=sum(not m.precedent.reference_id for m in supporting),
    )
    has_adverse = counts.low_yield_count + counts.zero_yield_count > 0
    measured_status = {
        (False, False): "no_comparable_numeric_outcome",
        (True, False): "only_at_or_above_threshold",
        (False, True): "only_low_or_zero",
        (True, True): "mixed_at_or_above_threshold_and_low_or_zero",
    }[counts.at_or_above_threshold_yield_count > 0, has_adverse]
    grade = "no_qualified_precedent"
    if pool.status == "unresolved":
        grade = "unresolved"
    elif counts.supporting_source_group_count >= policy.repeated_support_minimum_source_groups:
        grade = "repeated_transferable_precedent"
    elif counts.supporting_source_group_count:
        grade = "transferable_precedent"
    elif close:
        grade = "context_supported"
    elif matches:
        grade = "structural_analogue"
    warnings = set(pool.warnings) | set(query.warnings)
    warnings.update(("REFERENCE_GROUPS_ARE_PROXIES_FOR_INDEPENDENT_REPORTS",
                     "SOURCE_OUTCOMES_APPLY_TO_RECORDED_RECIPES"))
    if has_adverse:
        warnings.add("COMPARABLE_LOW_OR_ZERO_OUTCOMES_INSPECT_SOURCE_CONDITIONS")
    if any(m.excluded_by_low_yield and m.outcome not in ("low_yield", "zero_yield") for m in comparable):
        warnings.add("LOW_YIELD_REPORTED_FOR_SAME_OBSERVATION")
    if any(m.context.status == "missing" for m in matches):
        warnings.add("SOME_QUALIFIED_PRECEDENTS_MISS_REQUIRED_PARTNER_CONTEXT")
    if any(m.context.status == "unresolved" for m in matches):
        warnings.add("SOME_QUALIFIED_PRECEDENT_CONTEXT_UNRESOLVED")
    if any(m.context.source_only for m in matches):
        warnings.add("SOME_SOURCES_HAVE_ADDITIONAL_FUNCTIONAL_CONTEXT")
    if any(m.context.source_unresolved for m in matches):
        warnings.add("SOME_SOURCE_TO_QUERY_CONTEXT_UNRESOLVED")
    if query.status == "resolved" and not query.requirements:
        warnings.add("NO_TRACKED_SPECTATOR_REQUIREMENTS")
    if any(
        m.precedent.level in policy.transferable_levels
        and m.precedent.level != "whole_reaction"
        and m.context.minimum_partner_similarity is not None
        and m.context.minimum_partner_similarity < policy.minimum_partner_similarity
        for m in close
    ):
        warnings.add("SOME_LOCAL_MATCHES_HAVE_DISSIMILAR_REACTING_PARTNERS")
    if counts.supporting_reference_count > counts.supporting_source_group_count:
        warnings.add("POSSIBLE_DUPLICATE_SUPPORTING_REPORTS_GROUPED")
    if query.has_specified_stereo and any(m.precedent.level != "whole_reaction" for m in close):
        warnings.add("ANALOGUE_STEREOCHEMICAL_TRANSFER_UNVERIFIED")
    if query.status != "resolved":
        warnings.add("QUERY_FUNCTIONAL_CONTEXT_UNRESOLVED")
    assessment = StepPlausibilityAssessment(
        grade=grade, evidence_rank=dict(policy.ranks)[grade], strongest_level=pool.strongest_level,
        measured_evidence_status=measured_status, query_reaction_smiles=reaction_smiles,
        candidate_count=pool.candidate_count, qualified_count=pool.qualified_count,
        candidate_limit=candidate_limit, match_limit=match_limit,
        candidate_truncated=pool.candidate_truncated, matches_truncated=len(matches) > match_limit,
        level_counts=pool.level_counts, query_context=query, support_counts=counts,
        supporting_source_groups=tuple(tuple(refs) for _, refs in sorted(supporting_groups.items())),
        matches=_display_sample(matches, match_limit), exclusions=pool.exclusions,
        warnings=tuple(sorted(warnings)), policy=policy, index_source_identity=shared_index.source_identity,
        shared_core_definition_hash=pool.definition_hash,
    )
    return assessment, matches
