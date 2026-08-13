"""Hierarchical SITE1 to SYN1/REAL1 ranking with completion priors.

The ranker consumes only structurally validated candidates and evidence already
stored in a generic operator library.  It never generates reactions, changes a
graph-validation decision, or uses a named reaction family as a routing key.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Iterable, Mapping, Sequence, Tuple

from .generic_models import (
    GenericDisconnectionCandidate,
    GenericHandleCompletionGroup,
    GenericTemplateLibrary,
)
from .ranking_policy import (
    RetrosynthesisRankingPolicy,
    load_retrosynthesis_ranking_policy,
    structural_score_bands,
)


HIERARCHICAL_RANKING_POLICY_PATH = (
    Path(__file__).with_name("definitions") / "hierarchical_ranking.v1.json"
)
_BACKOFF_ORDER = ("operator_synthon", "operator", "global")


@dataclass(frozen=True)
class HierarchicalRankingPolicy:
    """Validated immutable policy for evidence-scoped hierarchical ranking."""

    definition_id: str
    schema_version: str
    required_forward_validation_status: str
    preserve_abstraction_level_order: bool
    preserve_structural_score_bands: bool
    backoff_order: Tuple[str, ...]
    minimum_context_independent_support: int
    laplace_alpha: float
    independent_support_saturation: int
    site_weights: Tuple[Tuple[str, float], ...]
    synthon_weights: Tuple[Tuple[str, float], ...]
    realization_weights: Tuple[Tuple[str, float], ...]
    interleave_sites_within_partition: bool
    interleave_synthons_within_site: bool

    def weights(self, stage: str) -> Mapping[str, float]:
        """Return the configured component weights for one ranking stage."""

        values = {
            "site": self.site_weights,
            "synthon": self.synthon_weights,
            "realization": self.realization_weights,
        }
        try:
            return dict(values[stage])
        except KeyError as error:
            raise ValueError(f"unsupported hierarchical ranking stage: {stage}") from error


@dataclass(frozen=True)
class CompletionPriorEvidence:
    """One smoothed completion prior with its complete evidence trace."""

    completion_group_id: str
    probability: float | None
    backoff_level: str
    independent_support: int
    total_context_support: int
    alternative_count: int


@dataclass(frozen=True)
class CompletionPriorIndex:
    """Read-only completion evidence indexed by operator and synthon identity."""

    groups: Tuple[GenericHandleCompletionGroup, ...]
    group_by_template_id: Mapping[str, GenericHandleCompletionGroup]
    group_by_realization_id: Mapping[str, GenericHandleCompletionGroup]
    groups_by_operator: Mapping[str, Tuple[GenericHandleCompletionGroup, ...]]
    groups_by_operator_synthon: Mapping[
        Tuple[str, str], Tuple[GenericHandleCompletionGroup, ...]
    ]

    def assess(
        self,
        candidate: GenericDisconnectionCandidate,
        *,
        policy: HierarchicalRankingPolicy | None = None,
    ) -> CompletionPriorEvidence:
        """Return the narrowest supported prior for a candidate completion."""

        resolved = policy or load_hierarchical_ranking_policy()
        group = self.group_by_template_id.get(candidate.template_id)
        if group is None and candidate.realization_id:
            group = self.group_by_realization_id.get(candidate.realization_id)
        if group is None:
            return CompletionPriorEvidence("", None, "unavailable", 0, 0, 0)

        contexts: dict[str, Tuple[GenericHandleCompletionGroup, ...]] = {}
        if candidate.operator_id and candidate.synthon_signature:
            contexts["operator_synthon"] = self.groups_by_operator_synthon.get(
                (candidate.operator_id, candidate.synthon_signature),
                (),
            )
        if candidate.operator_id:
            contexts["operator"] = self.groups_by_operator.get(
                candidate.operator_id,
                (),
            )
        contexts["global"] = self.groups

        for level in resolved.backoff_order:
            alternatives = contexts.get(level, ())
            if group not in alternatives:
                continue
            total_support = sum(
                item.independent_reference_support for item in alternatives
            )
            if total_support < resolved.minimum_context_independent_support:
                continue
            alternative_count = len(alternatives)
            denominator = (
                total_support + resolved.laplace_alpha * alternative_count
            )
            probability = (
                group.independent_reference_support + resolved.laplace_alpha
            ) / denominator
            return CompletionPriorEvidence(
                completion_group_id=group.completion_group_id,
                probability=round(probability, 8),
                backoff_level=level,
                independent_support=group.independent_reference_support,
                total_context_support=total_support,
                alternative_count=alternative_count,
            )
        return CompletionPriorEvidence(
            completion_group_id=group.completion_group_id,
            probability=None,
            backoff_level="unavailable",
            independent_support=group.independent_reference_support,
            total_context_support=0,
            alternative_count=0,
        )


def _positive_number(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must be a positive number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized <= 0.0:
        raise ValueError(f"{field} must be a positive number")
    return normalized


def _positive_integer(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{field} must be a positive integer")
    return value


def _validated_weights(
    value: object,
    *,
    required: frozenset[str],
    field: str,
) -> Tuple[Tuple[str, float], ...]:
    if not isinstance(value, dict) or set(value) != set(required):
        raise ValueError(f"{field} weights must contain exactly {sorted(required)}")
    weights = []
    for name in sorted(required):
        raw = value[name]
        if isinstance(raw, bool) or not isinstance(raw, (int, float)):
            raise ValueError(f"{field} {name} weight must be nonnegative")
        weight = float(raw)
        if not math.isfinite(weight) or weight < 0.0:
            raise ValueError(f"{field} {name} weight must be nonnegative")
        weights.append((name, weight))
    if not math.isclose(sum(weight for _, weight in weights), 1.0, abs_tol=1e-12):
        raise ValueError(f"{field} weights must sum to one")
    return tuple(weights)


@lru_cache(maxsize=1)
def load_hierarchical_ranking_policy() -> HierarchicalRankingPolicy:
    """Load and validate the canonical hierarchical ranking definition."""

    value = json.loads(HIERARCHICAL_RANKING_POLICY_PATH.read_text(encoding="utf-8"))
    if value.get("definition_id") != "hierarchical_retrosynthesis_ranking.v1":
        raise ValueError("unexpected hierarchical ranking definition ID")
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported hierarchical ranking schema")
    applicability = value.get("applicability")
    prior = value.get("completion_prior")
    scores = value.get("scores")
    selection = value.get("selection")
    if not all(
        isinstance(section, dict)
        for section in (applicability, prior, scores, selection)
    ):
        raise ValueError("hierarchical ranking definition has incomplete sections")
    required_status = applicability.get("required_forward_validation_status")
    if not isinstance(required_status, str) or not required_status.strip():
        raise ValueError("hierarchical ranking requires a validation status")
    preserve_levels = applicability.get("preserve_abstraction_level_order")
    preserve_bands = applicability.get("preserve_structural_score_bands")
    if preserve_levels is not True or preserve_bands is not True:
        raise ValueError("hierarchical ranking must preserve levels and score bands")
    backoff = tuple(prior.get("backoff_order") or ())
    if backoff != _BACKOFF_ORDER:
        raise ValueError("completion-prior backoff order is unsupported")
    site = _validated_weights(
        scores.get("site"),
        required=frozenset({"structural_score", "independent_support"}),
        field="site",
    )
    synthon = _validated_weights(
        scores.get("synthon"),
        required=frozenset({"structural_score", "independent_support"}),
        field="synthon",
    )
    realization = _validated_weights(
        scores.get("realization"),
        required=frozenset(
            {"structural_score", "independent_support", "completion_prior"}
        ),
        field="realization",
    )
    interleave_sites = selection.get("interleave_sites_within_partition")
    interleave_synthons = selection.get("interleave_synthons_within_site")
    if not isinstance(interleave_sites, bool) or not isinstance(
        interleave_synthons, bool
    ):
        raise ValueError("hierarchical interleaving settings must be boolean")
    return HierarchicalRankingPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        required_forward_validation_status=required_status,
        preserve_abstraction_level_order=preserve_levels,
        preserve_structural_score_bands=preserve_bands,
        backoff_order=backoff,
        minimum_context_independent_support=_positive_integer(
            prior.get("minimum_context_independent_support"),
            "minimum completion-prior context support",
        ),
        laplace_alpha=_positive_number(
            prior.get("laplace_alpha"),
            "completion-prior Laplace alpha",
        ),
        independent_support_saturation=_positive_integer(
            scores.get("independent_support_saturation"),
            "independent-support saturation",
        ),
        site_weights=site,
        synthon_weights=synthon,
        realization_weights=realization,
        interleave_sites_within_partition=interleave_sites,
        interleave_synthons_within_site=interleave_synthons,
    )


def _unique_membership(
    groups: Sequence[GenericHandleCompletionGroup],
    field: str,
) -> Mapping[str, GenericHandleCompletionGroup]:
    values: dict[str, GenericHandleCompletionGroup] = {}
    for group in groups:
        for identifier in getattr(group, field):
            if not identifier:
                continue
            current = values.get(identifier)
            if current is not None and current.completion_group_id != group.completion_group_id:
                raise ValueError(
                    f"{field} member belongs to multiple completion groups: {identifier}"
                )
            values[identifier] = group
    return MappingProxyType(values)


def build_completion_prior_index(
    library: GenericTemplateLibrary | object,
) -> CompletionPriorIndex:
    """Build deterministic read-only lookup tables from library evidence."""

    groups = tuple(
        sorted(
            getattr(library, "completion_groups", ()) or (),
            key=lambda group: group.completion_group_id,
        )
    )
    by_operator: dict[str, list[GenericHandleCompletionGroup]] = {}
    by_operator_synthon: dict[
        Tuple[str, str], list[GenericHandleCompletionGroup]
    ] = {}
    for group in groups:
        by_operator.setdefault(group.operator_id, []).append(group)
        for synthon_signature in group.synthon_signatures:
            if synthon_signature:
                by_operator_synthon.setdefault(
                    (group.operator_id, synthon_signature),
                    [],
                ).append(group)
    return CompletionPriorIndex(
        groups=groups,
        group_by_template_id=_unique_membership(groups, "template_ids"),
        group_by_realization_id=_unique_membership(groups, "realization_ids"),
        groups_by_operator=MappingProxyType(
            {key: tuple(values) for key, values in sorted(by_operator.items())}
        ),
        groups_by_operator_synthon=MappingProxyType(
            {
                key: tuple(values)
                for key, values in sorted(by_operator_synthon.items())
            }
        ),
    )


def _support_score(support: int, saturation: int) -> float:
    return min(1.0, math.log1p(max(0, support)) / math.log1p(saturation))


def _weighted_score(
    components: Mapping[str, float | None],
    weights: Mapping[str, float],
) -> float:
    available = {
        name: float(value)
        for name, value in components.items()
        if value is not None and weights.get(name, 0.0) > 0.0
    }
    denominator = sum(weights[name] for name in available)
    if denominator <= 0.0:
        return 0.0
    return round(
        sum(weights[name] * value for name, value in available.items())
        / denominator,
        8,
    )


def _round_robin(
    groups: Sequence[Sequence[GenericDisconnectionCandidate]],
    *,
    enabled: bool,
) -> list[GenericDisconnectionCandidate]:
    if not enabled:
        return [candidate for group in groups for candidate in group]
    values: list[GenericDisconnectionCandidate] = []
    for depth in range(max((len(group) for group in groups), default=0)):
        for group in groups:
            if depth < len(group):
                values.append(group[depth])
    return values


def rank_hierarchical_candidates(
    candidates: Iterable[GenericDisconnectionCandidate],
    library: GenericTemplateLibrary | object,
    *,
    policy: HierarchicalRankingPolicy | None = None,
    structural_policy: RetrosynthesisRankingPolicy | None = None,
    prior_index: CompletionPriorIndex | None = None,
) -> Tuple[GenericDisconnectionCandidate, ...]:
    """Rank validated candidates as SITE1, then SYN1, then REAL1 choices.

    Candidates can move only inside the same abstraction-level and effective
    structural-score-band partition.  Sites and synthons are round-robin
    interleaved after evidence ranking, retaining chemistry-distinct choices.
    """

    resolved = policy or load_hierarchical_ranking_policy()
    ranking_policy = structural_policy or load_retrosynthesis_ranking_policy()
    values = tuple(candidates)
    invalid = tuple(
        candidate.forward_validation_status
        for candidate in values
        if candidate.forward_validation_status
        != resolved.required_forward_validation_status
    )
    if invalid:
        raise ValueError(
            "hierarchical ranking requires forward-validated candidates; "
            f"found {sorted(set(invalid))}"
        )
    if not values:
        return ()

    calculated_bands = structural_score_bands(
        values,
        width=ranking_policy.diversity_score_band_width,
        separate_by_level=resolved.preserve_abstraction_level_order,
    )

    def partition(candidate: GenericDisconnectionCandidate) -> Tuple[str, int]:
        level = candidate.abstraction_level if resolved.preserve_abstraction_level_order else "*"
        if candidate.ranking_policy_definition_id == ranking_policy.definition_id:
            band = candidate.effective_structural_score_band
        else:
            band = calculated_bands[id(candidate)]
        return level, band

    structurally_ordered = tuple(
        sorted(
            values,
            key=lambda candidate: (
                ranking_policy.level_rank(partition(candidate)[0]),
                partition(candidate)[1],
                -candidate.score,
                -candidate.independent_reference_support,
                candidate.precursor_smiles,
                candidate.template_id,
            ),
        )
    )
    pre_ranks = {
        id(candidate): rank
        for rank, candidate in enumerate(structurally_ordered, start=1)
    }
    completion_index = prior_index or build_completion_prior_index(library)
    annotated = []
    for candidate in structurally_ordered:
        evidence = completion_index.assess(candidate, policy=resolved)
        annotated.append(
            replace(
                candidate,
                completion_group_id=evidence.completion_group_id,
                completion_prior_probability=evidence.probability,
                completion_prior_backoff_level=evidence.backoff_level,
                completion_prior_independent_support=evidence.independent_support,
                completion_prior_total_support=evidence.total_context_support,
                completion_prior_alternative_count=evidence.alternative_count,
                hierarchical_partition_key=partition(candidate),
                pre_hierarchical_rank=pre_ranks[id(candidate)],
                hierarchical_ranking_definition_id=resolved.definition_id,
            )
        )

    partitions: dict[
        Tuple[str, int], list[GenericDisconnectionCandidate]
    ] = {}
    for candidate in annotated:
        partitions.setdefault(candidate.hierarchical_partition_key, []).append(candidate)

    ranked: list[GenericDisconnectionCandidate] = []
    for partition_key in sorted(
        partitions,
        key=lambda key: (ranking_policy.level_rank(key[0]), key[1]),
    ):
        site_groups: dict[str, list[GenericDisconnectionCandidate]] = {}
        for candidate in partitions[partition_key]:
            site_key = candidate.disconnection_site_key or candidate.precursor_smiles
            site_groups.setdefault(site_key, []).append(candidate)
        site_scores = {
            site_key: _weighted_score(
                {
                    "structural_score": max(value.score for value in group),
                    "independent_support": max(
                        _support_score(
                            value.independent_reference_support,
                            resolved.independent_support_saturation,
                        )
                        for value in group
                    ),
                },
                resolved.weights("site"),
            )
            for site_key, group in site_groups.items()
        }
        ordered_sites = sorted(
            site_groups,
            key=lambda site_key: (
                -site_scores[site_key],
                min(value.pre_hierarchical_rank for value in site_groups[site_key]),
                site_key,
            ),
        )
        ranked_site_groups: list[list[GenericDisconnectionCandidate]] = []
        for site_rank, site_key in enumerate(ordered_sites, start=1):
            synthon_groups: dict[str, list[GenericDisconnectionCandidate]] = {}
            for candidate in site_groups[site_key]:
                synthon_key = (
                    candidate.synthon_signature
                    or candidate.realization_id
                    or candidate.precursor_smiles
                )
                synthon_groups.setdefault(synthon_key, []).append(candidate)
            synthon_scores = {
                synthon_key: _weighted_score(
                    {
                        "structural_score": max(value.score for value in group),
                        "independent_support": max(
                            _support_score(
                                value.independent_reference_support,
                                resolved.independent_support_saturation,
                            )
                            for value in group
                        ),
                    },
                    resolved.weights("synthon"),
                )
                for synthon_key, group in synthon_groups.items()
            }
            ordered_synthons = sorted(
                synthon_groups,
                key=lambda synthon_key: (
                    -synthon_scores[synthon_key],
                    min(
                        value.pre_hierarchical_rank
                        for value in synthon_groups[synthon_key]
                    ),
                    synthon_key,
                ),
            )
            ranked_synthon_groups: list[list[GenericDisconnectionCandidate]] = []
            for synthon_rank, synthon_key in enumerate(ordered_synthons, start=1):
                group = synthon_groups[synthon_key]
                realization_scores = {
                    id(candidate): _weighted_score(
                        {
                            "structural_score": candidate.score,
                            "independent_support": _support_score(
                                candidate.independent_reference_support,
                                resolved.independent_support_saturation,
                            ),
                            "completion_prior": candidate.completion_prior_probability,
                        },
                        resolved.weights("realization"),
                    )
                    for candidate in group
                }
                ordered_realizations = sorted(
                    group,
                    key=lambda candidate: (
                        -realization_scores[id(candidate)],
                        candidate.pre_hierarchical_rank,
                        candidate.realization_id,
                        candidate.precursor_smiles,
                        candidate.template_id,
                    ),
                )
                ranked_synthon_groups.append(
                    [
                        replace(
                            candidate,
                            hierarchical_site_score=site_scores[site_key],
                            hierarchical_synthon_score=synthon_scores[synthon_key],
                            hierarchical_realization_score=realization_scores[
                                id(candidate)
                            ],
                            hierarchical_site_rank=site_rank,
                            hierarchical_synthon_rank=synthon_rank,
                            hierarchical_realization_rank=realization_rank,
                        )
                        for realization_rank, candidate in enumerate(
                            ordered_realizations,
                            start=1,
                        )
                    ]
                )
            ranked_site_groups.append(
                _round_robin(
                    ranked_synthon_groups,
                    enabled=resolved.interleave_synthons_within_site,
                )
            )
        ranked.extend(
            _round_robin(
                ranked_site_groups,
                enabled=resolved.interleave_sites_within_partition,
            )
        )
    return tuple(
        replace(candidate, hierarchical_rank=rank)
        for rank, candidate in enumerate(ranked, start=1)
    )


__all__ = [
    "HIERARCHICAL_RANKING_POLICY_PATH",
    "CompletionPriorEvidence",
    "CompletionPriorIndex",
    "HierarchicalRankingPolicy",
    "build_completion_prior_index",
    "load_hierarchical_ranking_policy",
    "rank_hierarchical_candidates",
]
