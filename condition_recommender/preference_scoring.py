"""Transparent partner-category and functional-group tolerance evidence."""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple

from .ranking_preferences import functional_group_distance_decay


_ALIPHATIC_CLASSES = {"alkyl", "ring_aliphatic", "alkenyl", "alkynyl"}
_AROMATIC_CLASSES = {"aryl", "heteroaryl"}


def _multiset_jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    a = Counter(str(value) for value in left)
    b = Counter(str(value) for value in right)
    if not a or not b:
        return 0.0
    keys = set(a).union(b)
    intersection = sum(min(a[key], b[key]) for key in keys)
    union = sum(max(a[key], b[key]) for key in keys)
    return intersection / union if union else 0.0


def _attachment_family(classes: Sequence[str]) -> str:
    values = {value for value in classes if value}
    if not values:
        return "unresolved"
    if values <= _ALIPHATIC_CLASSES:
        return "aliphatic"
    if values <= _AROMATIC_CLASSES:
        return "aryl"
    if values & _ALIPHATIC_CLASSES and values & _AROMATIC_CLASSES:
        return "mixed_aryl_aliphatic"
    if "acyl" in values:
        return "acyl"
    return "mixed_other"


def _substitution_label(element: str, degree: int, hydrogens: int) -> str:
    if element == "N":
        if degree == 1 and hydrogens >= 1:
            return "primary"
        if degree == 2 and hydrogens >= 1:
            return "secondary"
        if degree == 3 and hydrogens == 0:
            return "tertiary"
        if degree >= 4:
            return "quaternary"
    return f"degree_{degree}"


def _category_label(
    *, element: str, substitution: str, family: str, aromatic: bool
) -> str:
    if element == "N":
        if family == "aliphatic":
            return f"{substitution} aliphatic amine"
        if family == "aryl":
            return f"{substitution} aryl amine"
        if family == "mixed_aryl_aliphatic":
            return f"{substitution} mixed aryl/aliphatic amine"
        if family == "acyl":
            return f"{substitution} amide-like nitrogen"
        return f"{substitution} nitrogen center"
    state = "aromatic" if aromatic else family.replace("_", " ")
    return f"{state} {element} center".strip()


def reacting_center_categories(
    reaction_core: Mapping[str, Any] | None,
) -> Tuple[dict[str, Any], ...]:
    """Derive compositional reactant categories from one serialized core.

    The categories are display/scoring evidence only. They never alter reaction
    identity, atom correspondence, admission, or retrieval keys.
    """
    if not reaction_core:
        return ()
    ports_by_center: dict[tuple[int, int], list[Mapping[str, Any]]] = defaultdict(list)
    for remote in reaction_core.get("remote_subgraphs") or ():
        if not isinstance(remote, Mapping) or remote.get("side") != "reactant":
            continue
        if str(remote.get("continuity") or "") == "departing":
            continue
        for port in remote.get("attachment_ports") or ():
            if not isinstance(port, Mapping) or port.get("side") != "reactant":
                continue
            key = (
                int(port.get("core_component_index") or 0),
                int(port.get("core_atom_index") or 0),
            )
            ports_by_center[key].append(port)

    categories = []
    for transition in reaction_core.get("atom_transitions") or ():
        if not isinstance(transition, Mapping):
            continue
        before = transition.get("before_state") or {}
        if not isinstance(before, Mapping):
            continue
        key = (
            int(before.get("component_index") or 0),
            int(before.get("atom_index") or 0),
        )
        ports = ports_by_center.get(key, [])
        stable_count = int(transition.get("stable_remote_subgraph_count") or 0)
        if not ports and stable_count == 0:
            continue
        element = str(before.get("element") or "")
        degree = int(before.get("heavy_atom_degree") or 0)
        hydrogens = int(before.get("total_hydrogens") or 0)
        aromatic = bool(before.get("aromatic"))
        classes = tuple(
            str((port.get("substituent_profile") or {}).get("base_class") or "")
            for port in ports
        )
        family = _attachment_family(classes)
        substitution = _substitution_label(element, degree, hydrogens)
        tokens = [
            f"element:{element}",
            f"substitution:{substitution}",
            f"attachment_family:{family}",
            f"hydrogens:{hydrogens}",
            f"charge:{int(before.get('formal_charge') or 0)}",
            f"aromatic:{int(aromatic)}",
            f"hybridization:{before.get('hybridization') or 'unknown'}",
        ]
        tokens.extend(f"remote_class:{value}" for value in classes if value)
        tokens.extend(
            "neighbor:" + "|".join(str(token).split("|")[:2])
            for token in before.get("neighbor_tokens") or ()
        )
        categories.append(
            {
                "label": _category_label(
                    element=element,
                    substitution=substitution,
                    family=family,
                    aromatic=aromatic,
                ),
                "element": element,
                "substitution": substitution,
                "attachment_family": family,
                "component_index": key[0],
                "atom_index": key[1],
                "role": str(transition.get("role") or "participant"),
                "tokens": tuple(sorted(tokens)),
            }
        )
    return tuple(
        sorted(
            categories,
            key=lambda value: (
                value["element"],
                value["role"],
                value["component_index"],
                value["atom_index"],
            ),
        )
    )


def assess_partner_category_similarity(
    query_core: Mapping[str, Any] | None,
    precedent_core: Mapping[str, Any] | None,
) -> tuple[Optional[float], dict[str, Any]]:
    """Compare graph-derived reacting-center categories with an audit trace."""
    query = reacting_center_categories(query_core)
    precedent = reacting_center_categories(precedent_core)
    if not query or not precedent:
        return None, {
            "status": "unavailable",
            "query_categories": tuple(item["label"] for item in query),
            "precedent_categories": tuple(item["label"] for item in precedent),
            "matches": (),
        }
    unused = set(range(len(precedent)))
    matches = []
    scores = []
    for query_category in query:
        candidates = [
            (index, candidate)
            for index, candidate in enumerate(precedent)
            if index in unused and candidate["element"] == query_category["element"]
        ]
        if not candidates:
            scores.append(0.0)
            matches.append(
                {
                    "query": query_category["label"],
                    "precedent": None,
                    "score": 0.0,
                    "status": "missing_matching_center",
                }
            )
            continue
        ranked = [
            (
                _multiset_jaccard(
                    query_category["tokens"], candidate["tokens"]
                ),
                index,
                candidate,
            )
            for index, candidate in candidates
        ]
        score, selected, candidate = max(
            ranked, key=lambda value: (value[0], -value[1])
        )
        unused.remove(selected)
        scores.append(score)
        matches.append(
            {
                "query": query_category["label"],
                "precedent": candidate["label"],
                "score": round(score, 6),
                "status": "exact_category"
                if query_category["label"] == candidate["label"]
                else "category_mismatch",
            }
        )
    denominator = max(len(query), len(precedent))
    score = sum(scores) / denominator if denominator else 0.0
    return round(score, 6), {
        "status": "available",
        "query_categories": tuple(item["label"] for item in query),
        "precedent_categories": tuple(item["label"] for item in precedent),
        "matches": tuple(matches),
    }


def aggregate_partner_category_evidence(
    query_core: Mapping[str, Any] | None,
    precedents: Iterable[tuple[str, Mapping[str, Any]]],
) -> tuple[Optional[float], dict[str, Any]]:
    """Average category similarity over independent precedents."""
    assessments = []
    for evidence_id, core in precedents:
        score, evidence = assess_partner_category_similarity(query_core, core)
        assessments.append(
            {"evidence_id": evidence_id, "score": score, **evidence}
        )
    available = [item for item in assessments if item["score"] is not None]
    if not available:
        query_categories = reacting_center_categories(query_core)
        return None, {
            "status": "unavailable",
            "query_categories": tuple(
                item["label"] for item in query_categories
            ),
            "precedent_assessments": tuple(assessments),
        }
    score = sum(float(item["score"]) for item in available) / len(available)
    return round(score, 6), {
        "status": "available",
        "query_categories": tuple(
            item["label"] for item in reacting_center_categories(query_core)
        ),
        "independent_precedent_count": len(available),
        "precedent_assessments": tuple(available[:5]),
    }


def assess_functional_group_tolerance(
    query_signature: Mapping[str, Any],
    precedents: Iterable[tuple[str, Mapping[str, Any]]],
) -> tuple[Optional[float], dict[str, Any]]:
    """Score direct unchanged-group precedent coverage; absence stays unknown."""
    query_groups: dict[str, Mapping[str, Any]] = {}
    for group in query_signature.get("spectator_groups") or ():
        if not isinstance(group, Mapping) or not group.get("group_id"):
            continue
        query_groups.setdefault(str(group["group_id"]), group)
    if not query_groups:
        return None, {"status": "not_applicable", "groups": ()}
    precedent_values = tuple(precedents)
    decay = functional_group_distance_decay()
    evidence = []
    total_weight = 0.0
    matched_weight = 0.0
    for group_id, group in sorted(query_groups.items()):
        distance_value = group.get("graph_distance")
        distance = int(distance_value) if distance_value is not None else None
        weight = 1.0 / (1.0 + decay * distance) if distance is not None else 1.0
        matched_ids = []
        for evidence_id, signature in precedent_values:
            candidate_ids = {
                str(value.get("group_id"))
                for value in signature.get("spectator_groups") or ()
                if isinstance(value, Mapping) and value.get("group_id")
            }
            if group_id in candidate_ids:
                matched_ids.append(evidence_id)
        total_weight += weight
        if matched_ids:
            matched_weight += weight
        evidence.append(
            {
                "group_id": group_id,
                "chemist_label": str(group.get("chemist_label") or group_id),
                "graph_distance": distance,
                "distance_weight": round(weight, 6),
                "matched_independent_count": len(set(matched_ids)),
                "status": "directly_demonstrated"
                if matched_ids
                else "unknown_not_tolerant",
                "evidence_ids": tuple(sorted(set(matched_ids))[:5]),
            }
        )
    score = matched_weight / total_weight if total_weight else 0.0
    return round(score, 6), {
        "status": "available",
        "score_basis": "direct_unchanged_group_precedent_coverage",
        "groups": tuple(evidence),
    }


__all__ = [
    "aggregate_partner_category_evidence",
    "assess_functional_group_tolerance",
    "assess_partner_category_similarity",
    "reacting_center_categories",
]
