"""Reaction-label condition retrieval over the cleaned HTE dataset."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Tuple

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.reaction_grammars import load_reaction_grammars

from .label_indexing import LabelIndexedReaction, LabelParticipant, load_label_index
from .models import LabelConditionRecommendation, LabelRecommendationResult


_RULES_PATH = Path(__file__).with_name("definitions") / "label_retrieval.v1.json"


@dataclass(frozen=True)
class QueryParticipant:
    role: str
    signature: str
    label: str
    center_class: str | None
    attachment_classes: Tuple[str, ...]
    alpha_branched: bool | None


@lru_cache(maxsize=1)
def load_label_retrieval_rules() -> Dict[str, Any]:
    """Load and validate weak-label retrieval configuration."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    grammar_ids = {str(item["id"]) for item in load_reaction_grammars()}
    configured = set((rules.get("grammar_source_reaction_types") or {}).keys())
    unknown = sorted(configured - grammar_ids)
    if unknown:
        raise ValueError(f"Unknown configured reaction grammars: {unknown}")
    supported_scopes = set(rules.get("supported_reaction_scopes") or ())
    allowed_scopes = {
        "intramolecular",
        "intermolecular",
        "mixed",
        "unimolecular",
        "unresolved",
    }
    if not supported_scopes or not supported_scopes <= allowed_scopes:
        raise ValueError("Invalid supported_reaction_scopes")
    for key in ("participant_similarity_weights", "ranking_weights"):
        weights = rules.get(key) or {}
        if (
            not weights
            or abs(sum(float(value) for value in weights.values()) - 1.0) > 1e-9
        ):
            raise ValueError(f"Weights must sum to one: {key}")
    return rules


def _query_participants(analysis: Any) -> Tuple[QueryParticipant, ...]:
    selected = analysis.selected_candidate
    if selected is None:
        return ()
    environments: Dict[Tuple[int, str], Any] = {}
    for component in analysis.reactants:
        for environment in component.compound_analysis.site_environments:
            environments[(component.component_index, environment.site_id)] = environment

    participants = []
    for role, site in selected.role_assignments.items():
        environment = environments.get((site.component_index, site.site_id))
        steric = dict(environment.steric) if environment is not None else {}
        attached = list(steric.get("attached_groups") or [])
        alpha_values = [
            item.get("alpha_branched")
            for item in attached
            if item.get("alpha_branched") is not None
        ]
        participants.append(
            QueryParticipant(
                role=str(role),
                signature=str(site.canonical_signature),
                label=str(site.chemist_label),
                center_class=str(steric.get("center_substitution_class") or "") or None,
                attachment_classes=tuple(
                    sorted(
                        {
                            str(item["attachment_carbon_class"])
                            for item in attached
                            if item.get("attachment_carbon_class")
                        }
                    )
                ),
                alpha_branched=any(bool(value) for value in alpha_values)
                if alpha_values
                else None,
            )
        )
    return tuple(sorted(participants, key=lambda item: item.role))


def _signature_similarity(query: str, precedent: str) -> float:
    if not query or not precedent:
        return 0.0
    if query == precedent:
        return 1.0
    query_parts, precedent_parts = query.split("|"), precedent.split("|")
    if not query_parts or not precedent_parts or query_parts[0] != precedent_parts[0]:
        return 0.0
    width = max(len(query_parts), len(precedent_parts)) - 1
    if width <= 0:
        return 0.35
    matches = sum(
        left == right for left, right in zip(query_parts[1:], precedent_parts[1:])
    )
    return round(0.35 + 0.65 * matches / width, 6)


def _qualifier_similarity(
    query: QueryParticipant, precedent: LabelParticipant
) -> float:
    comparisons: List[bool] = []
    if precedent.center_class:
        comparisons.append(query.center_class == precedent.center_class)
    if precedent.attachment_class:
        comparisons.append(precedent.attachment_class in query.attachment_classes)
    if precedent.alpha_branched is not None:
        comparisons.append(query.alpha_branched is precedent.alpha_branched)
    if not comparisons:
        return 0.5
    signed_evidence = sum(1.0 if comparison else -1.0 for comparison in comparisons)
    return max(0.0, min(1.0, 0.5 + 0.5 * signed_evidence / 3.0))


def _pair_similarity(
    query: Tuple[QueryParticipant, QueryParticipant],
    precedent: Tuple[LabelParticipant, LabelParticipant],
) -> Tuple[float, float, float]:
    weights = load_label_retrieval_rules()["participant_similarity_weights"]
    alternatives = []
    for order in ((0, 1), (1, 0)):
        signatures = [
            _signature_similarity(
                query[index].signature, precedent[order[index]].signature
            )
            for index in range(2)
        ]
        qualifiers = [
            _qualifier_similarity(query[index], precedent[order[index]])
            for index in range(2)
        ]
        signature_score = sum(signatures) / 2.0
        qualifier_score = sum(qualifiers) / 2.0
        label_score = (
            float(weights["signature"]) * signature_score
            + float(weights["qualifier"]) * qualifier_score
        )
        alternatives.append((label_score, signature_score, qualifier_score))
    return max(alternatives, key=lambda item: (item[0], item[1], item[2]))


def _sigmoid(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-max(-20.0, min(20.0, value))))


def _rank_recipes(
    query: Tuple[QueryParticipant, QueryParticipant],
    rows: Tuple[LabelIndexedReaction, ...],
    *,
    top_k: int,
) -> Tuple[LabelConditionRecommendation, ...]:
    rules = load_label_retrieval_rules()
    scored = []
    for row in rows:
        label_score, signature_score, qualifier_score = _pair_similarity(
            query, row.participants
        )
        if signature_score <= 0.0:
            continue
        scored.append((label_score, signature_score, qualifier_score, row))

    groups: Dict[str, List[Tuple[float, float, float, LabelIndexedReaction]]] = (
        defaultdict(list)
    )
    for item in scored:
        groups[item[3].recipe_id].append(item)

    ranking_weights = rules["ranking_weights"]
    maximum_precedents = int(rules.get("maximum_precedents_per_recipe", 5))
    ranked = []
    for recipe_id, members in groups.items():
        members.sort(
            key=lambda item: (-item[0], -item[3].yield_pct, item[3].source_row_number)
        )
        best_label, best_signature, best_qualifier, _ = members[0]
        similarity_weights = [max(item[0], 0.05) for item in members]
        weight_sum = sum(similarity_weights)
        expected_yield = (
            sum(
                weight * item[3].yield_pct
                for weight, item in zip(similarity_weights, members)
            )
            / weight_sum
        )
        mean_z_score = (
            sum(
                weight * item[3].z_score
                for weight, item in zip(similarity_weights, members)
            )
            / weight_sum
        )
        support_score = min(1.0, math.log1p(len(members)) / math.log1p(10))
        score = (
            float(ranking_weights["label_similarity"]) * best_label
            + float(ranking_weights["yield"]) * expected_yield / 100.0
            + float(ranking_weights["z_score"]) * _sigmoid(mean_z_score)
            + float(ranking_weights["support"]) * support_score
        )
        ranked.append(
            (
                score,
                expected_yield,
                mean_z_score,
                best_label,
                best_signature,
                best_qualifier,
                recipe_id,
                members,
            )
        )

    ranked.sort(key=lambda item: (-item[0], -item[1], item[6]))
    recommendations = []
    for rank, item in enumerate(ranked[:top_k], start=1):
        (
            score,
            expected_yield,
            mean_z_score,
            label_score,
            signature_score,
            qualifier_score,
            recipe_id,
            members,
        ) = item
        best_row = members[0][3]
        explanation = ["Compatible reactive-taxonomy grammar"]
        if signature_score == 1.0:
            explanation.append("Exact unordered FG-signature pair")
        else:
            explanation.append("Partial FG-signature pair match")
        if qualifier_score == 1.0:
            explanation.append("Source steric qualifiers match the query")
        recommendations.append(
            LabelConditionRecommendation(
                rank=rank,
                recipe_id=recipe_id,
                score=round(score, 6),
                label_similarity=round(label_score, 6),
                signature_similarity=round(signature_score, 6),
                qualifier_similarity=round(qualifier_score, 6),
                expected_yield_pct=round(expected_yield, 2),
                mean_z_score=round(mean_z_score, 3),
                support=len(members),
                source_reaction_types=tuple(
                    sorted({member[3].reaction_type for member in members})
                ),
                source_row_numbers=tuple(
                    member[3].source_row_number
                    for member in members[:maximum_precedents]
                ),
                conditions=dict(best_row.conditions),
                explanation=tuple(explanation),
            )
        )
    return tuple(recommendations)


def recommend_conditions_from_labels(
    reaction_smiles: str,
    *,
    records_path: str | Path = "datasets/reaction_label/v2.1_cleaned.csv",
    top_k: int = 5,
) -> LabelRecommendationResult:
    """Return top condition recipes using reaction grammar and FG signatures."""
    if top_k < 1:
        return LabelRecommendationResult(
            reaction_smiles, False, error="TOP_K_MUST_BE_POSITIVE"
        )
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid:
        return LabelRecommendationResult(
            reaction_smiles,
            False,
            error=analysis.error or "INVALID_REACTION",
        )
    if (
        analysis.evidence_quality != "exact_product_reconstruction"
        or analysis.selected_candidate is None
    ):
        return LabelRecommendationResult(
            reaction_smiles,
            False,
            query_label=analysis.reaction_label,
            error="QUERY_REACTION_NOT_EXACTLY_VERIFIED",
        )

    selected = analysis.selected_candidate
    rules = load_label_retrieval_rules()
    topology = analysis.reaction_topology
    reaction_scope = topology.reaction_scope if topology is not None else "unresolved"
    supported_scopes = set(rules.get("supported_reaction_scopes") or ())
    if reaction_scope not in supported_scopes:
        return LabelRecommendationResult(
            reaction_smiles,
            False,
            query_label=analysis.reaction_label,
            grammar_id=selected.grammar_id,
            warnings=("LABEL_DATASET_HAS_NO_REACTION_TOPOLOGY",),
            error="QUERY_TOPOLOGY_NOT_SUPPORTED_BY_LABEL_DATASET",
        )
    source_types = tuple(
        (rules.get("grammar_source_reaction_types") or {}).get(selected.grammar_id, ())
    )
    if not source_types:
        return LabelRecommendationResult(
            reaction_smiles,
            False,
            query_label=analysis.reaction_label,
            grammar_id=selected.grammar_id,
            error="QUERY_GRAMMAR_NOT_SUPPORTED_BY_LABEL_DATASET",
        )
    participants = _query_participants(analysis)
    if len(participants) != 2:
        return LabelRecommendationResult(
            reaction_smiles,
            False,
            query_label=analysis.reaction_label,
            grammar_id=selected.grammar_id,
            error="QUERY_REQUIRES_TWO_REACTIVE_PARTICIPANTS",
        )

    rows = tuple(
        row
        for row in load_label_index(records_path)
        if row.reaction_type in source_types
    )
    if not rows:
        return LabelRecommendationResult(
            reaction_smiles,
            False,
            query_label=analysis.reaction_label,
            grammar_id=selected.grammar_id,
            query_signatures=tuple(item.signature for item in participants),
            error="NO_FAMILY_COMPATIBLE_PRECEDENTS",
        )
    recommendations = _rank_recipes(participants, rows, top_k=top_k)
    if not recommendations:
        return LabelRecommendationResult(
            reaction_smiles,
            False,
            query_label=analysis.reaction_label,
            grammar_id=selected.grammar_id,
            query_signatures=tuple(item.signature for item in participants),
            candidate_count=len(rows),
            error="NO_SIGNATURE_COMPATIBLE_PRECEDENTS",
        )
    recipe_count = len({row.recipe_id for row in rows})
    return LabelRecommendationResult(
        query_reaction_smiles=reaction_smiles,
        valid=True,
        query_label=analysis.reaction_label,
        grammar_id=selected.grammar_id,
        query_signatures=tuple(item.signature for item in participants),
        candidate_count=len(rows),
        recipe_count=recipe_count,
        recommendations=recommendations,
        warnings=(
            "WEAK_LABEL_PRECEDENTS_NOT_STRUCTURE_VERIFIED",
            "REACTIVE_SITES_MATCHED_AS_UNORDERED_PARTICIPANTS",
            "CONDITION_NAMES_NOT_REGISTRY_NORMALIZED",
        ),
    )


__all__ = ["load_label_retrieval_rules", "recommend_conditions_from_labels"]
