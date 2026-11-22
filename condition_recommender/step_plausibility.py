"""Assess the strength of indexed precedents for a proposed reaction step.

The grade describes documented examples under their own conditions. A missing
outcome contributes structural evidence only. Source groups are conservative
proxies for repeated reports, not verified independent replications.
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


_POLICY_PATH = Path(__file__).with_name("definitions") / "step_plausibility.v1.json"
_LEVELS = ("whole_reaction", "observed_local", "retained_local", "retained_typed")
_GRADES = (
    "no_qualified_precedent", "structural_analogue", "context_supported",
    "measured_support", "repeated_measured_support",
)


@dataclass(frozen=True)
class StepPlausibilityPolicy:
    """Versioned ordinal rules and explicitly heuristic cutoffs."""

    definition_id: str
    definition_version: str
    score_semantics: str
    close_levels: tuple[str, ...]
    measured_levels: tuple[str, ...]
    measured_analogue_context: str
    minimum_partner_similarity: float
    positive_yield_strictly_above_pct: float
    repeated_support_minimum_source_groups: int
    ranks: tuple[tuple[str, int | None], ...]


def load_step_plausibility_policy() -> StepPlausibilityPolicy:
    """Load and validate the experimental standard without mutating shared state."""
    raw = json.loads(_POLICY_PATH.read_text(encoding="utf-8"))
    close, measured = tuple(raw["close_levels"]), tuple(raw["measured_levels"])
    sim = raw["minimum_partner_similarity"]
    pct = raw["positive_yield_strictly_above_pct"]
    repeats = raw["repeated_support_minimum_source_groups"]
    expected_ranks = {"unresolved": None, **{grade: i for i, grade in enumerate(_GRADES)}}
    if (
        raw["definition_id"] != "step_plausibility.v1"
        or raw["schema_version"] != "1.0"
        or not isinstance(raw["definition_version"], str) or not raw["definition_version"]
        or raw["status"] != "development_pending_independent_review"
        or raw["score_semantics"] != "ordinal_strength_of_documented_structural_support"
        or raw["measured_analogue_context"] != "bidirectional"
        or raw["ranks"] != expected_ranks
        or any(type(raw["ranks"][grade]) is not int for grade in _GRADES)
        or not measured or not close or measured != close[:len(measured)]
        or close != _LEVELS[:len(close)] or "retained_local" in measured
        or type(sim) not in (int, float) or not 0 < sim <= 1
        or type(pct) not in (int, float) or not 0 <= pct < 100
        or type(repeats) is not int or repeats < 2
    ):
        raise ValueError("Invalid step plausibility policy")
    return StepPlausibilityPolicy(
        raw["definition_id"], raw["definition_version"], raw["score_semantics"],
        close, measured, raw["measured_analogue_context"], sim, pct, repeats,
        tuple(raw["ranks"].items()),
    )


@dataclass(frozen=True)
class StepPlausibilityMatch:
    """One audited source observation, its chemical context and measurement."""

    precedent: ReactionPrecedentEvidence
    context: StepContextComparison
    outcome: str
    context_covered_close: bool
    comparable_for_measured_support: bool
    source_group_id: str | None = None


@dataclass(frozen=True)
class StepSupportCounts:
    """Counts in the full visited pool; outcome counts cover comparable rows only."""

    distinct_reference_count: int
    context_covered_close_observation_count: int
    context_covered_close_reference_count: int
    comparable_observation_count: int
    positive_yield_count: int
    positive_reference_count: int
    positive_distinct_reaction_count: int
    positive_source_group_count: int
    low_yield_count: int
    zero_yield_count: int
    unreported_yield_count: int
    unusable_yield_count: int
    unknown_yield_count: int
    missing_reference_positive_count: int


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
    positive_source_groups: tuple[tuple[str, ...], ...]
    matches: tuple[StepPlausibilityMatch, ...]
    exclusions: tuple[tuple[str, int], ...]
    warnings: tuple[str, ...]
    policy: StepPlausibilityPolicy
    index_source_identity: str
    shared_core_definition_hash: str | None
    schema_version: str = "step_plausibility.v1"

    def to_dict(self) -> dict[str, Any]:
        """Return source-linked evidence and the rules that assigned its grade."""
        return asdict(self)


def _classify_outcome(row: ReactionPrecedentEvidence, threshold: float) -> str:
    value = row.yield_pct
    if value is None:
        return "unusable" if row.outcome_status == "invalid" else "unreported"
    if (row.outcome_status != "usable" or type(value) not in (int, float)
            or not isfinite(value) or not 0 <= value <= 100):
        return "unusable"
    if value == 0:
        return "zero_yield"
    return "above_threshold" if value > threshold else "low_yield"


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
    graphs, canonical recipe ID, and exact numeric outcome. With a missing recipe
    identity on a matching graph/outcome, group the compatible reports
    conservatively. Group membership is transitive; this can undercount genuine
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

    same_observation: dict[str, set[str]] = defaultdict(set)
    experiments: dict[tuple[str, str, float], dict[str, set[str]]] = {}
    for match in matches:
        row = match.precedent
        if not row.reference_id:
            continue
        if row.observation_id:
            same_observation[row.observation_id].add(row.reference_id)
        if match.outcome in ("unreported", "unusable") or match.context.status == "unresolved":
            continue
        left, right = normalized_reaction(row.reaction_smiles)
        key = (left, right, float(row.yield_pct))
        recipe_id = str((row.resolved_recipe or {}).get("recipe_id") or "")
        experiments.setdefault(key, {}).setdefault(recipe_id, set()).add(row.reference_id)
    for refs in same_observation.values():
        connect(refs)
    for by_recipe in experiments.values():
        if "" in by_recipe:
            connect(set.union(*by_recipe.values()))
        else:
            for refs in by_recipe.values():
                connect(refs)
    return {reference: representative(reference) for reference in parents}


def _display_sample(
    matches: list[StepPlausibilityMatch], limit: int,
) -> tuple[StepPlausibilityMatch, ...]:
    """Show a matched counterexample early and avoid filling the page with one report."""
    def priority(m: StepPlausibilityMatch) -> tuple:
        row = m.precedent
        similarity = m.context.minimum_partner_similarity
        return (not m.comparable_for_measured_support, not m.context_covered_close,
                _LEVELS.index(row.level), -(similarity if similarity is not None else -1),
                row.reference_id, row.reaction_id, row.observation_id or "", row.reaction_smiles)

    ordered = sorted(range(len(matches)), key=lambda i: priority(matches[i]))
    first_positive: dict[str, int] = {}
    adverse = []
    for i in ordered:
        m = matches[i]
        if not m.comparable_for_measured_support:
            continue
        if m.source_group_id and m.outcome == "above_threshold":
            first_positive.setdefault(m.source_group_id, i)
        elif m.outcome in ("low_yield", "zero_yield"):
            adverse.append(i)
    positives = list(first_positive.values())
    priority_ids = dict.fromkeys(positives[:1] + adverse[:1] + positives[1:] + ordered)
    return tuple(matches[i] for i in list(priority_ids)[:limit])


def assess_step_plausibility(
    reaction_smiles: str, index: GenericReactionIndex, shared_index: SharedCoreIndex,
    *, candidate_limit: int = 128, match_limit: int = 20,
) -> StepPlausibilityAssessment:
    """Grade graph-qualified precedents with transferable context and measured support.

    Both indexes are the bound inputs used for condition recommendation. Every
    qualified visited row enters the assessment before the evidence display cap.
    A numeric yield describes the source experiment under its own recipe.
    """
    if type(match_limit) is not int or not 1 <= match_limit <= 50:
        raise ValueError("match_limit must be an integer between 1 and 50")
    policy = load_step_plausibility_policy()
    pool = collect_reaction_precedent_support(
        reaction_smiles, index, shared_index, candidate_limit=candidate_limit,
    )
    context_for = lru_cache(maxsize=None)(build_step_context)
    query = context_for(reaction_smiles)
    matches = []
    for row in pool.matches:
        context = compare_step_contexts(query, context_for(row.reaction_smiles))
        covered = context.status == "covered" and row.level in policy.close_levels
        min_similarity = context.minimum_partner_similarity
        similar_partners = (
            min_similarity is not None and min_similarity >= policy.minimum_partner_similarity
        )
        transferable_analogue = (
            similar_partners and context.reverse_status == "covered"
            and not query.has_specified_stereo
        )
        matches.append(StepPlausibilityMatch(
            row, context, _classify_outcome(row, policy.positive_yield_strictly_above_pct),
            covered, covered and row.level in policy.measured_levels
            and (row.level == "whole_reaction" or transferable_analogue),
        ))
    normalized_reaction = lru_cache(maxsize=None)(_normalized_reaction)
    groups = _source_groups(matches, normalized_reaction)
    matches = [replace(m, source_group_id=groups.get(m.precedent.reference_id)) for m in matches]
    close = [m for m in matches if m.context_covered_close]
    comparable = [m for m in close if m.comparable_for_measured_support]
    outcomes = Counter(m.outcome for m in comparable)
    positive = [m for m in comparable if m.outcome == "above_threshold"]
    positive_references = {m.precedent.reference_id for m in positive if m.precedent.reference_id}
    positive_groups: dict[str, list[str]] = defaultdict(list)
    for reference in sorted(positive_references):
        positive_groups[groups[reference]].append(reference)
    counts = StepSupportCounts(
        distinct_reference_count=pool.distinct_reference_count,
        context_covered_close_observation_count=len(close),
        context_covered_close_reference_count=len({
            m.precedent.reference_id for m in close if m.precedent.reference_id
        }),
        comparable_observation_count=len(comparable),
        positive_yield_count=len(positive),
        positive_reference_count=len(positive_references),
        positive_distinct_reaction_count=len({
            normalized_reaction(m.precedent.reaction_smiles)
            for m in positive if m.precedent.reference_id
        }),
        positive_source_group_count=len(positive_groups),
        low_yield_count=outcomes["low_yield"],
        zero_yield_count=outcomes["zero_yield"],
        unreported_yield_count=outcomes["unreported"],
        unusable_yield_count=outcomes["unusable"],
        unknown_yield_count=outcomes["unreported"] + outcomes["unusable"],
        missing_reference_positive_count=sum(not m.precedent.reference_id for m in positive),
    )
    has_adverse = counts.low_yield_count + counts.zero_yield_count > 0
    measured_status = {
        (False, False): "no_comparable_numeric_outcome",
        (True, False): "only_above_threshold",
        (False, True): "only_low_or_zero",
        (True, True): "mixed_above_threshold_and_low_or_zero",
    }[bool(positive), has_adverse]
    grade = "no_qualified_precedent"
    if pool.status == "unresolved":
        grade = "unresolved"
    elif counts.positive_source_group_count >= policy.repeated_support_minimum_source_groups:
        grade = "repeated_measured_support"
    elif counts.positive_source_group_count:
        grade = "measured_support"
    elif close:
        grade = "context_supported"
    elif matches:
        grade = "structural_analogue"
    warnings = set(pool.warnings) | set(query.warnings)
    warnings.update(("REFERENCE_GROUPS_ARE_PROXIES_FOR_INDEPENDENT_REPORTS",
                     "SOURCE_OUTCOMES_APPLY_TO_RECORDED_RECIPES"))
    if has_adverse:
        warnings.add("COMPARABLE_LOW_OR_ZERO_OUTCOMES_INSPECT_SOURCE_CONDITIONS")
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
        m.precedent.level in policy.measured_levels
        and m.precedent.level != "whole_reaction"
        and m.context.minimum_partner_similarity is not None
        and m.context.minimum_partner_similarity < policy.minimum_partner_similarity
        for m in close
    ):
        warnings.add("SOME_LOCAL_MATCHES_HAVE_DISSIMILAR_REACTING_PARTNERS")
    if counts.positive_reference_count > counts.positive_source_group_count:
        warnings.add("POSSIBLE_DUPLICATE_POSITIVE_REPORTS_GROUPED")
    if query.has_specified_stereo and any(m.precedent.level != "whole_reaction" for m in close):
        warnings.add("ANALOGUE_STEREOCHEMICAL_TRANSFER_UNVERIFIED")
    if query.status != "resolved":
        warnings.add("QUERY_FUNCTIONAL_CONTEXT_UNRESOLVED")
    return StepPlausibilityAssessment(
        grade=grade, evidence_rank=dict(policy.ranks)[grade], strongest_level=pool.strongest_level,
        measured_evidence_status=measured_status, query_reaction_smiles=reaction_smiles,
        candidate_count=pool.candidate_count, qualified_count=pool.qualified_count,
        candidate_limit=candidate_limit, match_limit=match_limit,
        candidate_truncated=pool.candidate_truncated, matches_truncated=len(matches) > match_limit,
        level_counts=pool.level_counts, query_context=query, support_counts=counts,
        positive_source_groups=tuple(tuple(refs) for _, refs in sorted(positive_groups.items())),
        matches=_display_sample(matches, match_limit), exclusions=pool.exclusions,
        warnings=tuple(sorted(warnings)), policy=policy, index_source_identity=shared_index.source_identity,
        shared_core_definition_hash=pool.definition_hash,
    )
