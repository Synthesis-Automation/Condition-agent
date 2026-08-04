"""Chemistry-first exploratory retrieval for unknown reaction development."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple

from reactive_taxonomy import (
    AtomMappingProvider,
    analyze_reaction_with_external_mapping,
    featurize_reaction,
)

from .discovery_models import (
    DiscoveryScoreTrace,
    ReactionDiscoveryHit,
    ReactionDiscoveryResult,
)
from .edit_prototypes import (
    AnonymousEditPrototype,
    anonymous_edit_prototype,
    anonymous_edit_prototype_from_hypothesis,
    anonymous_edit_similarity,
)
from .evaluation_features import reaction_scaffold_tokens
from .generic_api import GenericConditionRecommender
from .generic_indexing import GenericIndexedReaction, GenericReactionIndex
from .preference_scoring import assess_partner_category_similarity
from .signature_features import environment_tokens
from .similarity import assess_signature_similarity, jaccard


_RULES_PATH = Path(__file__).with_name("definitions") / "discovery_retrieval.v1.json"
_FACTORS = (
    "edit_similarity",
    "reaction_center",
    "local_environment",
    "partner_category",
    "spectator_groups",
    "reaction_topology",
    "reactive_scaffold",
)
_VIEWS = {
    "closest_chemistry",
    "diverse_strategies",
    "successful_precedents",
    "failure_informed",
}


@lru_cache(maxsize=1)
def load_discovery_rules() -> Dict[str, Any]:
    """Load and validate the versioned exploratory retrieval definition."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if rules.get("schema_version") != "1.0":
        raise ValueError("unsupported discovery retrieval schema")
    if rules.get("definition_id") != "discovery_retrieval.v1":
        raise ValueError("unexpected discovery retrieval definition ID")
    weights = rules.get("weights")
    if not isinstance(weights, Mapping) or set(weights) != set(_FACTORS):
        raise ValueError("discovery weights do not match factor vocabulary")
    normalized = {str(key): float(value) for key, value in weights.items()}
    if any(value < 0.0 for value in normalized.values()):
        raise ValueError("discovery weights must be non-negative")
    if abs(sum(normalized.values()) - 1.0) > 1e-9:
        raise ValueError("discovery weights must sum to one")
    if not 0.0 <= float(rules["edit_graph_min_similarity"]) <= 1.0:
        raise ValueError("discovery edit-graph threshold must be in [0, 1]")
    if int(rules["maximum_candidates"]) < 1:
        raise ValueError("discovery maximum_candidates must be positive")
    rules["weights"] = normalized
    return rules


@dataclass(frozen=True)
class _QueryAlternative:
    hypothesis_id: Optional[str]
    signature: Optional[Dict[str, Any]]
    reaction_core: Optional[Dict[str, Any]]
    prototype: Optional[AnonymousEditPrototype]


def _reaction_center_similarity(
    query: Mapping[str, Any] | None,
    precedent: Mapping[str, Any] | None,
) -> Optional[float]:
    if not query or not precedent:
        return None
    if query.get("exact_core_key") == precedent.get("exact_core_key"):
        return 1.0
    if query.get("typed_core_key") == precedent.get("typed_core_key"):
        return 0.8
    if query.get("shape_core_key") == precedent.get("shape_core_key"):
        return 0.6
    query_transitions = {
        str(value.get("transition_key") or "")
        for value in query.get("atom_transitions") or ()
        if isinstance(value, Mapping) and value.get("transition_key")
    }
    precedent_transitions = {
        str(value.get("transition_key") or "")
        for value in precedent.get("atom_transitions") or ()
        if isinstance(value, Mapping) and value.get("transition_key")
    }
    return jaccard(query_transitions, precedent_transitions)


def _available_weighted_score(
    components: Mapping[str, Optional[float]],
    configured: Mapping[str, float],
) -> tuple[float, Dict[str, float], Dict[str, float]]:
    available = {
        name: float(value) for name, value in components.items() if value is not None
    }
    total = sum(configured[name] for name in available)
    if not available or total <= 0.0:
        return 0.0, {}, {}
    effective = {name: configured[name] / total for name in available}
    contributions = {
        name: round(effective[name] * available[name], 6) for name in available
    }
    return round(sum(contributions.values()), 6), effective, contributions


def _candidate_positions(
    alternative: _QueryAlternative,
    index: GenericReactionIndex,
) -> Dict[int, set[str]]:
    tiers: Dict[int, set[str]] = defaultdict(set)

    def add(level: str, positions: Iterable[int]) -> None:
        for position in positions:
            tiers[int(position)].add(level)

    signature = alternative.signature
    if signature:
        direct = (
            ("exact_signature", index.exact, "exact_signature_key"),
            ("handle_signature", index.handles, "handle_signature_key"),
            (
                "transformation_signature",
                index.transformations,
                "transformation_signature_key",
            ),
            ("bond_edit_signature", index.bond_edits, "bond_edit_signature_key"),
        )
        for level, mapping, field in direct:
            add(level, mapping.get(str(signature.get(field) or ""), ()))
        for token in environment_tokens(signature):
            add("environment_neighbor", index.environment_features.get(token, ()))
    core = alternative.reaction_core
    if core:
        for level, mapping, field in (
            ("reaction_core_exact", index.core_exact, "exact_core_key"),
            ("reaction_core_local", index.core_typed, "typed_core_key"),
            ("reaction_core_context", index.core_shapes, "shape_core_key"),
        ):
            add(level, mapping.get(str(core.get(field) or ""), ()))
    prototype = alternative.prototype
    if prototype is not None:
        candidate_lookup = getattr(index.rows, "edit_graph_candidate_positions", None)
        candidates = (
            candidate_lookup(prototype)
            if callable(candidate_lookup)
            else range(len(index.rows))
        )
        threshold = float(load_discovery_rules()["edit_graph_min_similarity"])
        for position in candidates:
            row = index.rows[position]
            candidate = (
                anonymous_edit_prototype(row.signature) if row.signature else None
            )
            if candidate is None:
                continue
            if anonymous_edit_similarity(prototype, candidate) >= threshold:
                tiers[int(position)].add("edit_graph_neighbor")
    return tiers


def _relation_class(tiers: Iterable[str]) -> str:
    values = set(tiers)
    if "exact_signature" in values:
        return "direct_signature_analogue"
    if values.intersection({"handle_signature", "bond_edit_signature"}):
        return "direct_edit_analogue"
    if values.intersection(
        {"environment_neighbor", "reaction_core_exact", "reaction_core_local"}
    ):
        return "local_environment_analogue"
    if values.intersection({"reaction_core_context", "transformation_signature"}):
        return "reaction_core_analogue"
    return "edit_graph_analogue"


def _score_row(
    query_smiles: str,
    alternative: _QueryAlternative,
    row: GenericIndexedReaction,
    tiers: Iterable[str],
) -> tuple[float, DiscoveryScoreTrace]:
    signature_assessment = (
        assess_signature_similarity(
            alternative.signature,
            row.signature,
            query_reaction_core=alternative.reaction_core,
            precedent_reaction_core=row.reaction_core,
        )
        if alternative.signature and row.signature
        else None
    )
    precedent_prototype = (
        anonymous_edit_prototype(row.signature) if row.signature else None
    )
    edit_score = (
        anonymous_edit_similarity(alternative.prototype, precedent_prototype)
        if alternative.prototype is not None and precedent_prototype is not None
        else None
    )
    partner_score, partner_evidence = assess_partner_category_similarity(
        alternative.reaction_core,
        row.reaction_core,
    )
    query_scaffolds = (
        reaction_scaffold_tokens(query_smiles, alternative.signature)
        if alternative.signature
        else ()
    )
    components: Dict[str, Optional[float]] = {
        "edit_similarity": edit_score,
        "reaction_center": _reaction_center_similarity(
            alternative.reaction_core, row.reaction_core
        ),
        "local_environment": (
            signature_assessment.components["environment"]
            if signature_assessment
            else None
        ),
        "partner_category": partner_score,
        "spectator_groups": (
            signature_assessment.components["spectators"]
            if signature_assessment
            else None
        ),
        "reaction_topology": (
            signature_assessment.components["reaction_topology"]
            if signature_assessment
            else None
        ),
        "reactive_scaffold": (
            jaccard(query_scaffolds, row.scaffold_tokens)
            if query_scaffolds and row.scaffold_tokens
            else None
        ),
    }
    rules = load_discovery_rules()
    score, effective, contributions = _available_weighted_score(
        components, rules["weights"]
    )
    matches = [f"retrieval:{value}" for value in sorted(set(tiers))]
    mismatches = []
    for match in partner_evidence.get("matches") or ():
        if not isinstance(match, Mapping):
            continue
        text = f"reactant category: {match.get('query')} vs {match.get('precedent')}"
        (matches if match.get("status") == "exact_category" else mismatches).append(
            text
        )
    for name, value in components.items():
        if value is not None and value >= 0.8:
            matches.append(f"{name}={value:.3f}")
        elif value is not None and value < 0.4:
            mismatches.append(f"{name}={value:.3f}")
    trace = DiscoveryScoreTrace(
        components={
            name: (round(value, 6) if value is not None else None)
            for name, value in components.items()
        },
        contributions=contributions,
        configured_weights=dict(rules["weights"]),
        effective_weights={name: round(value, 6) for name, value in effective.items()},
        matches=tuple(dict.fromkeys(matches)),
        mismatches=tuple(dict.fromkeys(mismatches)),
        definition_id=str(rules["definition_id"]),
        definition_version=str(rules["schema_version"]),
    )
    return score, trace


def _hit_sort_key(hit: ReactionDiscoveryHit, view: str) -> tuple[Any, ...]:
    if view == "successful_precedents":
        outcome = hit.yield_pct if hit.yield_pct is not None else -1.0
        return (
            -hit.discovery_score,
            -outcome,
            hit.canonical_reaction_id,
            hit.reaction_id,
        )
    if view == "failure_informed":
        missing = hit.yield_pct is None
        outcome = hit.yield_pct if hit.yield_pct is not None else -1.0
        return (
            missing,
            outcome,
            -hit.discovery_score,
            hit.canonical_reaction_id,
            hit.reaction_id,
        )
    return (
        -hit.discovery_score,
        hit.canonical_reaction_id,
        hit.reaction_id,
        hit.observation_id,
    )


def _diversify(
    hits: Iterable[ReactionDiscoveryHit], top_k: int
) -> Tuple[ReactionDiscoveryHit, ...]:
    buckets: Dict[str, list[ReactionDiscoveryHit]] = defaultdict(list)
    for hit in hits:
        buckets[hit.relation_class].append(hit)
    ordered_classes = tuple(load_discovery_rules()["relation_priority"])
    result = []
    while len(result) < top_k and any(buckets.values()):
        for relation in ordered_classes:
            if buckets[relation] and len(result) < top_k:
                result.append(buckets[relation].pop(0))
    return tuple(result)


@dataclass(frozen=True)
class ReactionDiscoveryExplorer:
    """Reusable explorer backed by the same validated reaction SQLite index."""

    index: GenericReactionIndex
    source_path: str = ""
    mapping_provider: AtomMappingProvider | None = None

    @classmethod
    def from_path(
        cls,
        path: str | Path,
        *,
        mapping_provider: AtomMappingProvider | None = None,
        include_review: bool = False,
    ) -> "ReactionDiscoveryExplorer":
        loaded = GenericConditionRecommender.from_path(
            path,
            mapping_provider=mapping_provider,
            include_review=include_review,
        )
        return cls(loaded.index, loaded.source_path, loaded.mapping_provider)

    def discover(
        self,
        reaction_smiles: str,
        *,
        top_k: int = 20,
        view: str = "closest_chemistry",
        include_low_yield: bool = True,
        include_unreported_outcomes: bool = True,
    ) -> ReactionDiscoveryResult:
        """Find structural analogues and expose their observed conditions."""
        rules = load_discovery_rules()
        definition_version = f"{rules['definition_id']}@{rules['schema_version']}"
        if top_k < 1:
            return ReactionDiscoveryResult(
                reaction_smiles, False, error="TOP_K_MUST_BE_POSITIVE"
            )
        if view not in _VIEWS:
            return ReactionDiscoveryResult(
                reaction_smiles, False, error="UNSUPPORTED_DISCOVERY_VIEW"
            )
        if not self.index.rows:
            return ReactionDiscoveryResult(
                reaction_smiles, False, error="EMPTY_GENERIC_INDEX"
            )
        base = featurize_reaction(reaction_smiles)
        assessment = (
            analyze_reaction_with_external_mapping(
                reaction_smiles, self.mapping_provider, base_analysis=base
            )
            if self.mapping_provider is not None and base.valid
            else None
        )
        analysis = assessment.analysis if assessment is not None else base
        if not analysis.valid:
            return ReactionDiscoveryResult(
                reaction_smiles, False, error=analysis.error or "INVALID_REACTION"
            )
        signature = (
            asdict(analysis.reaction_signature) if analysis.reaction_signature else None
        )
        core = asdict(analysis.reaction_core) if analysis.reaction_core else None
        alternatives = []
        if signature:
            alternatives.append(
                _QueryAlternative(
                    None, signature, core, anonymous_edit_prototype(signature)
                )
            )
        else:
            for hypothesis in analysis.edit_hypotheses:
                payload = asdict(hypothesis)
                alternatives.append(
                    _QueryAlternative(
                        str(payload.get("hypothesis_id") or "") or None,
                        None,
                        core,
                        anonymous_edit_prototype_from_hypothesis(payload),
                    )
                )
        if not alternatives and core:
            alternatives.append(_QueryAlternative(None, None, core, None))
        hypothesis_ids = tuple(
            str(value.hypothesis_id) for value in alternatives if value.hypothesis_id
        )
        label = asdict(analysis.reaction_label) if analysis.reaction_label else {}
        context = dict(
            query_signature_id=str((signature or {}).get("signature_id") or "") or None,
            query_reaction_core_id=str((core or {}).get("core_id") or "") or None,
            query_edit_hypothesis_ids=hypothesis_ids,
            reaction_label=label,
            named_family=(signature or {}).get("named_family"),
            transformation_class=(signature or {}).get("transformation_class")
            or analysis.transformation_class,
            spectator_groups=tuple(
                asdict(value) for value in analysis.spectator_groups
            ),
            reaction_partners=tuple((signature or {}).get("partners") or ()),
            discovery_view=view,
            retrieval_definition_version=definition_version,
        )
        if not alternatives:
            return ReactionDiscoveryResult(
                reaction_smiles,
                False,
                **context,
                error="QUERY_HAS_NO_DISCOVERABLE_REACTION_EVIDENCE",
            )
        position_alternatives: Dict[int, list[tuple[_QueryAlternative, set[str]]]] = (
            defaultdict(list)
        )
        for alternative in alternatives:
            for position, tiers in _candidate_positions(
                alternative, self.index
            ).items():
                position_alternatives[position].append((alternative, tiers))
        maximum = int(rules["maximum_candidates"])
        if len(position_alternatives) > maximum:
            prioritized = sorted(
                position_alternatives,
                key=lambda position: (
                    min(
                        rules["relation_priority"].index(_relation_class(tiers))
                        for _, tiers in position_alternatives[position]
                    ),
                    position,
                ),
            )[:maximum]
            position_alternatives = {
                position: position_alternatives[position] for position in prioritized
            }
        rows = dict(
            zip(
                sorted(position_alternatives),
                self.index.select(sorted(position_alternatives)),
            )
        )
        hits = []
        low_threshold = float(rules["low_yield_threshold_pct"])
        for position, choices in position_alternatives.items():
            row = rows[position]
            if not include_unreported_outcomes and row.yield_pct is None:
                continue
            if (
                not include_low_yield
                and row.yield_pct is not None
                and row.yield_pct < low_threshold
            ):
                continue
            scored = []
            for alternative, tiers in choices:
                score, trace = _score_row(reaction_smiles, alternative, row, tiers)
                scored.append(
                    (score, alternative.hypothesis_id or "", alternative, tiers, trace)
                )
            score, _, alternative, tiers, trace = max(
                scored, key=lambda value: (value[0], value[1])
            )
            cautions = []
            if row.yield_pct is None:
                cautions.append(
                    "Outcome was not reported; this is structural evidence only"
                )
            elif row.yield_pct < low_threshold:
                cautions.append(
                    f"Low observed yield ({row.yield_pct:.1f}%) may reveal a failure boundary"
                )
            if alternative.hypothesis_id:
                cautions.append(
                    "Related to one ambiguous query edit hypothesis, not a verified query center"
                )
            if row.precedent_tier.value != "trusted":
                cautions.append("Review-tier precedent; expert review is required")
            insights = [
                f"Observed conditions from {row.source_dataset or 'dataset precedent'}",
                *trace.matches[:4],
            ]
            hits.append(
                ReactionDiscoveryHit(
                    rank=0,
                    reaction_id=row.reaction_id,
                    observation_id=row.observation_id,
                    canonical_reaction_id=row.canonical_reaction_id,
                    reaction_smiles=row.reaction_smiles,
                    reaction_label=dict(row.reaction_label),
                    relation_class=_relation_class(tiers),
                    relation_tiers=tuple(sorted(tiers)),
                    discovery_score=score,
                    yield_pct=row.yield_pct,
                    outcome_status=row.outcome_status,
                    evidence_tier=row.precedent_tier.value,
                    chemistry_status=row.chemistry_status,
                    source_dataset=row.source_dataset,
                    reference_id=row.reference_id,
                    resolved_recipe=dict(row.resolved_recipe),
                    recipe_id=row.recipe_id,
                    recipe_core_id=row.recipe_core_id,
                    hypothesis_id=alternative.hypothesis_id,
                    score_trace=trace,
                    insights=tuple(insights),
                    cautions=tuple(cautions),
                )
            )
        hits.sort(key=lambda value: _hit_sort_key(value, view))
        selected = (
            _diversify(hits, top_k)
            if view == "diverse_strategies"
            else tuple(hits[:top_k])
        )
        selected = tuple(
            replace(hit, rank=index + 1) for index, hit in enumerate(selected)
        )
        relation_counts = dict(Counter(hit.relation_class for hit in hits))
        warnings = ["DISCOVERY_CONDITIONS_ARE_OBSERVED_NOT_RECOMMENDED"]
        if hypothesis_ids:
            warnings.append("AMBIGUOUS_QUERY_HYPOTHESES_SCORED_SEPARATELY")
        if assessment is not None and assessment.mapping_result is not None:
            warnings.append(f"EXTERNAL_MAPPING_STATUS:{assessment.status}")
        if not selected:
            return ReactionDiscoveryResult(
                reaction_smiles,
                False,
                **context,
                candidate_count=len(position_alternatives),
                relation_counts=relation_counts,
                warnings=tuple(warnings),
                error="NO_STRUCTURAL_DISCOVERY_ANALOGUE",
            )
        return ReactionDiscoveryResult(
            reaction_smiles,
            True,
            **context,
            candidate_count=len(position_alternatives),
            relation_counts=relation_counts,
            hits=selected,
            warnings=tuple(warnings),
        )


__all__ = ["ReactionDiscoveryExplorer", "load_discovery_rules"]
