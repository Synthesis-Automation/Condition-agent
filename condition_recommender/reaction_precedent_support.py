"""Graph-qualified reaction evidence without condition-recipe aggregation.

This reads the canonical index and shared projections used by condition
retrieval. It does not recommend conditions or infer experimental success.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any, Mapping

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.shared_reaction_core import (
    LEVELS,
    SharedReactionCore,
    build_shared_reaction_core,
    compare_reaction_cores,
)

from .generic_indexing import GenericReactionIndex
from .shared_core_index import SharedCoreIndex


@dataclass(frozen=True)
class ReactionPrecedentEvidence:
    """One qualified source observation, with its original recipe and outcome."""

    reaction_id: str
    reference_id: str
    reaction_smiles: str
    level: str
    differences: tuple[str, ...]
    warnings: tuple[str, ...]
    source: str
    observation_id: str | None = None
    template_id: str | None = None
    yield_pct: float | None = None
    outcome_status: str = "unavailable"
    resolved_recipe: Mapping[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize source evidence without transferring its operating values."""
        return asdict(self)


@dataclass(frozen=True)
class ReactionPrecedentSupport:
    """Bounded support counts over qualified observations, before display limits."""

    status: str
    strongest_level: str | None
    matches: tuple[ReactionPrecedentEvidence, ...]
    level_counts: tuple[tuple[str, int], ...]
    distinct_reference_count: int
    candidate_count: int
    qualified_count: int
    candidate_truncated: bool
    matches_truncated: bool
    exclusions: tuple[tuple[str, int], ...]
    warnings: tuple[str, ...]
    definition_hash: str | None
    scope: str = "bounded_canonical_index_lookup"
    schema_version: str = "reaction_precedent_support.v1"

    def to_dict(self) -> dict[str, Any]:
        """Return deterministic nested JSON."""
        return asdict(self)


def reaction_evidence_core(reaction_smiles: str) -> SharedReactionCore | None:
    """Project only reactions with verified product accounting and graph evidence."""
    analysis = featurize_reaction(reaction_smiles)
    if (
        not analysis.valid
        or analysis.reaction_signature is None
        or analysis.reaction_core is None
        or analysis.reaction_completeness is None
        or analysis.reaction_completeness.status != "verified"
    ):
        return None
    return build_shared_reaction_core(
        reaction_smiles, asdict(analysis.reaction_signature),
        asdict(analysis.reaction_core),
    )


def qualify_reaction_precedent(
    query: SharedReactionCore, precedent: SharedReactionCore, *,
    reaction_id: str, reference_id: str, reaction_smiles: str,
    source: str, observation_id: str | None = None,
    template_id: str | None = None, yield_pct: float | None = None,
    outcome_status: str = "unavailable",
    resolved_recipe: Mapping[str, Any] | None = None,
) -> ReactionPrecedentEvidence | None:
    """Use the same original-query graph qualification as condition retrieval.

    Exactness comes from canonical reaction identity plus qualified edits,
    never from fingerprint equality, product-only identity, or family names.
    """
    comparison = compare_reaction_cores(query, precedent)
    if not comparison.eligible or comparison.level is None:
        return None
    return ReactionPrecedentEvidence(
        reaction_id=reaction_id, reference_id=reference_id,
        reaction_smiles=reaction_smiles, level=comparison.level,
        differences=comparison.differences, warnings=comparison.reasons,
        source=source, observation_id=observation_id, template_id=template_id,
        yield_pct=yield_pct, outcome_status=outcome_status,
        resolved_recipe=resolved_recipe,
    )


def summarize_reaction_support(
    matches: tuple[ReactionPrecedentEvidence, ...], *, match_limit: int,
    candidate_count: int, candidate_truncated: bool,
    exclusions: tuple[tuple[str, int], ...] = (),
    warnings: tuple[str, ...] = (), definition_hash: str | None = None,
    scope: str = "bounded_canonical_index_lookup",
) -> ReactionPrecedentSupport:
    """Deduplicate repeated observations; do not count templates as publications."""
    levels = ("whole_reaction", *LEVELS)
    unique: dict[tuple[str, str, str, str], ReactionPrecedentEvidence] = {}
    for match in matches:
        key = (match.reference_id, match.reaction_id,
               match.observation_id or "", match.reaction_smiles)
        prior = unique.get(key)
        if prior is None or levels.index(match.level) < levels.index(prior.level):
            unique[key] = match
    ordered = tuple(sorted(unique.values(), key=lambda item: (
        levels.index(item.level), item.reference_id, item.reaction_id,
        item.observation_id or "", item.template_id or "",
    )))
    counts = tuple((level, sum(item.level == level for item in ordered)) for level in levels)
    return ReactionPrecedentSupport(
        status="supported" if ordered else "no_qualified_precedent",
        strongest_level=ordered[0].level if ordered else None,
        matches=ordered[:match_limit], level_counts=counts,
        distinct_reference_count=len({item.reference_id for item in ordered if item.reference_id}),
        candidate_count=candidate_count, qualified_count=len(ordered),
        candidate_truncated=candidate_truncated, matches_truncated=len(ordered) > match_limit,
        exclusions=exclusions, warnings=tuple(sorted(set(warnings))),
        definition_hash=definition_hash, scope=scope,
    )


def assess_reaction_precedent_support(
    reaction_smiles: str, index: GenericReactionIndex, shared_index: SharedCoreIndex,
    *, candidate_limit: int = 128, match_limit: int = 20,
) -> ReactionPrecedentSupport:
    """Find the strongest qualified support in a bounded exact-to-broad ladder.

    Counts describe visited observations, not exhaustive corpus coverage or
    independent replication. Source recipes are retained without filtering on a
    different proposed recipe; recipe applicability is a separate assessment.
    """
    if type(candidate_limit) is not int or not 1 <= candidate_limit <= 512:
        raise ValueError("candidate_limit must be an integer between 1 and 512")
    if type(match_limit) is not int or not 1 <= match_limit <= 50:
        raise ValueError("match_limit must be an integer between 1 and 50")
    shared_index.validate_binding(index)
    query = reaction_evidence_core(reaction_smiles)
    if query is None or not query.levels:
        empty = summarize_reaction_support(
            (), match_limit=match_limit, candidate_count=0, candidate_truncated=False,
            warnings=("QUERY_GRAPH_EVIDENCE_UNAVAILABLE",),
        )
        return replace(empty, status="unresolved")
    keys = [("whole_reaction", query.reaction_identity)]
    keys.extend((level.level, level.key) for level in query.levels)
    # Product-side candidates still undergo edit qualification. This also
    # explains why same-product/different-transformation records were excluded.
    keys.append(("product_identity", query.product_identity))
    positions: list[int] = []
    seen: set[int] = set()
    truncated = False
    for kind, key in keys:
        found, capped = shared_index.lookup(kind, key, candidate_limit)
        truncated |= capped
        for position in found:
            if position in seen:
                continue
            if len(positions) == candidate_limit:
                truncated = True
                continue
            positions.append(position)
            seen.add(position)
    matches = []
    exclusions: dict[str, int] = {}
    for position, row in zip(positions, index.select(positions)):
        projection = shared_index.projection(position, row)
        evidence = qualify_reaction_precedent(
            query, projection, reaction_id=row.reaction_id, reference_id=row.reference_id,
            reaction_smiles=row.reaction_smiles, source="canonical_condition_index",
            observation_id=row.observation_id, yield_pct=row.yield_pct,
            outcome_status=row.outcome_status, resolved_recipe=row.resolved_recipe,
        )
        if evidence is not None:
            matches.append(evidence)
        else:
            comparison = compare_reaction_cores(query, projection)
            for reason in comparison.reasons:
                exclusions[reason] = exclusions.get(reason, 0) + 1
    return summarize_reaction_support(
        tuple(matches), match_limit=match_limit, candidate_count=len(positions),
        candidate_truncated=truncated, exclusions=tuple(sorted(exclusions.items())),
        warnings=(*query.warnings, "SHARED_CORE_PENDING_INDEPENDENT_REVIEW",
                  "COUNTS_DESCRIBE_VISITED_OBSERVATIONS_NOT_FULL_CORPUS"),
        definition_hash=query.definition_hash,
    )
