"""Declarative, graph-derived intramolecular reactive-pair assessments."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Sequence, Tuple

from rdkit import Chem

from .api import analyze_molecule
from .chemistry.rdkit_utils import parse_smiles
from .models import MoleculeAnalysis, ReactiveSiteHypothesis


REACTIVE_PAIR_INTERACTION_DEFINITION_ID = "reactive_pair_interactions.v1"
REACTIVE_PAIR_INTERACTION_SCHEMA_VERSION = "1.0"
REACTIVE_PAIR_INTERACTION_DEFINITION_PATH = (
    Path(__file__).with_name("definitions")
    / "reactive_pair_interactions.v1.json"
)

_SITE_TYPES = frozenset(
    {
        "leaving_group",
        "pronucleophile_XH",
        "nucleophile_anion",
        "transfer_group",
        "addition_donor",
        "eliminable_pair",
        "electrophilic_center",
        "aromatic_CH",
        "unsaturated_bond",
        "dipolar_group",
        "heteroatom_bond",
    }
)
_OPERATORS = frozenset({"eq", "in", "gte"})
_SEVERITIES = frozenset({"low", "medium", "high", "critical"})
_WARNING_STRENGTHS = frozenset({"advisory", "strong"})
_TOPOLOGY_EVALUATIONS = frozenset(
    {"none", "graph_distance", "potential_bond_closure"}
)


@dataclass(frozen=True)
class ReactivePairSiteReference:
    """Stable reference to one site participating in an interaction."""

    hypothesis_id: str
    component_index: int
    atom_index: int
    site_type: str
    canonical_signature: str
    chemist_label: str
    availability: str


@dataclass(frozen=True)
class ReactivePairInteractionAssessment:
    """One definition-derived pair of conflicting sites in one component."""

    assessment_id: str
    rule_id: str
    interaction_class: str
    scope: str
    component_index: int
    component_smiles: str
    left_site: ReactivePairSiteReference
    right_site: ReactivePairSiteReference
    graph_distance: int | None
    potential_closure_ring_size: int | None
    intrinsic_severity: str
    warning_strength: str
    message: str
    definition_id: str = REACTIVE_PAIR_INTERACTION_DEFINITION_ID
    schema_version: str = REACTIVE_PAIR_INTERACTION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible assessment."""

        return asdict(self)


def _selector_errors(selector: object, location: str) -> list[str]:
    errors: list[str] = []
    if not isinstance(selector, Mapping):
        return [f"{location}_must_be_object"]
    site_type = str(selector.get("site_type") or "")
    if site_type not in _SITE_TYPES:
        errors.append(f"{location}_unknown_site_type:{site_type}")
    atom_field = str(selector.get("interaction_atom_field") or "")
    if not atom_field.startswith("details."):
        errors.append(f"{location}_invalid_interaction_atom_field")
    constraints = selector.get("constraints")
    if not isinstance(constraints, list):
        errors.append(f"{location}_constraints_must_be_list")
        return errors
    for index, constraint in enumerate(constraints):
        prefix = f"{location}_constraint_{index}"
        if not isinstance(constraint, Mapping):
            errors.append(f"{prefix}_must_be_object")
            continue
        if not str(constraint.get("field") or ""):
            errors.append(f"{prefix}_missing_field")
        operator = str(constraint.get("operator") or "")
        if operator not in _OPERATORS:
            errors.append(f"{prefix}_unknown_operator:{operator}")
        if operator == "in":
            values = constraint.get("values")
            if not isinstance(values, list) or not values:
                errors.append(f"{prefix}_requires_values")
        elif operator in {"eq", "gte"} and "value" not in constraint:
            errors.append(f"{prefix}_requires_value")
    return errors


def validate_reactive_pair_interaction_definition(
    payload: Mapping[str, Any],
) -> list[str]:
    """Return deterministic validation errors for a pair-rule definition."""

    errors: list[str] = []
    if payload.get("definition_id") != REACTIVE_PAIR_INTERACTION_DEFINITION_ID:
        errors.append("unexpected_reactive_pair_definition_id")
    if payload.get("schema_version") != REACTIVE_PAIR_INTERACTION_SCHEMA_VERSION:
        errors.append("unsupported_reactive_pair_schema")
    if payload.get("allowed_scopes") != ["same_component"]:
        errors.append("reactive_pair_scope_must_be_same_component")
    rules = payload.get("rules")
    if not isinstance(rules, list) or not rules:
        return [*errors, "reactive_pair_rules_must_be_nonempty"]
    rule_ids: list[str] = []
    for index, rule in enumerate(rules):
        location = f"reactive_pair_rule_{index}"
        if not isinstance(rule, Mapping):
            errors.append(f"{location}_must_be_object")
            continue
        rule_id = str(rule.get("rule_id") or "")
        rule_ids.append(rule_id)
        if not rule_id:
            errors.append(f"{location}_missing_rule_id")
        if not str(rule.get("interaction_class") or ""):
            errors.append(f"{location}_missing_interaction_class")
        if rule.get("scope") != "same_component":
            errors.append(f"{location}_scope_must_be_same_component")
        severity = str(rule.get("intrinsic_severity") or "")
        if severity not in _SEVERITIES:
            errors.append(f"{location}_invalid_severity:{severity}")
        strength = str(rule.get("warning_strength") or "")
        if strength not in _WARNING_STRENGTHS:
            errors.append(f"{location}_invalid_warning_strength:{strength}")
        topology = str(rule.get("topology_evaluation") or "")
        if topology not in _TOPOLOGY_EVALUATIONS:
            errors.append(f"{location}_invalid_topology:{topology}")
        if not str(rule.get("message") or ""):
            errors.append(f"{location}_missing_message")
        errors.extend(_selector_errors(rule.get("left_site"), f"{location}_left"))
        errors.extend(_selector_errors(rule.get("right_site"), f"{location}_right"))
    if len(rule_ids) != len(set(rule_ids)):
        errors.append("duplicate_reactive_pair_rule_id")
    return errors


@lru_cache(maxsize=1)
def load_reactive_pair_interaction_definition() -> Mapping[str, Any]:
    """Load and validate the canonical declarative pair rules."""

    payload = json.loads(
        REACTIVE_PAIR_INTERACTION_DEFINITION_PATH.read_text(encoding="utf-8")
    )
    errors = validate_reactive_pair_interaction_definition(payload)
    if errors:
        raise ValueError("invalid reactive-pair definition: " + ", ".join(errors))
    return payload


def _field_value(site: ReactiveSiteHypothesis, field: str) -> Any:
    value: Any = site
    for part in field.split("."):
        if isinstance(value, Mapping):
            value = value.get(part)
        else:
            value = getattr(value, part, None)
    return value


def _matches_constraint(site: ReactiveSiteHypothesis, constraint: Mapping[str, Any]) -> bool:
    value = _field_value(site, str(constraint["field"]))
    operator = str(constraint["operator"])
    if operator == "eq":
        return value == constraint.get("value")
    if operator == "in":
        return value in constraint.get("values", ())
    if operator == "gte":
        threshold = constraint.get("value")
        return (
            isinstance(value, (int, float))
            and not isinstance(value, bool)
            and isinstance(threshold, (int, float))
            and not isinstance(threshold, bool)
            and value >= threshold
        )
    raise ValueError(f"unsupported reactive-pair constraint operator: {operator}")


def _matches_selector(site: ReactiveSiteHypothesis, selector: Mapping[str, Any]) -> bool:
    return site.site_type == selector["site_type"] and all(
        _matches_constraint(site, constraint)
        for constraint in selector["constraints"]
    )


def _interaction_atom(site: ReactiveSiteHypothesis, selector: Mapping[str, Any]) -> int:
    value = _field_value(site, str(selector["interaction_atom_field"]))
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(
            f"reactive-pair site {site.hypothesis_id} has no integer interaction atom"
        )
    return value


def _site_reference(
    site: ReactiveSiteHypothesis,
    selector: Mapping[str, Any],
) -> ReactivePairSiteReference:
    return ReactivePairSiteReference(
        hypothesis_id=site.hypothesis_id,
        component_index=site.component_index,
        atom_index=_interaction_atom(site, selector),
        site_type=site.site_type,
        canonical_signature=site.canonical_signature,
        chemist_label=site.chemist_label,
        availability=site.availability,
    )


def _assessment_id(value: Sequence[Any]) -> str:
    encoded = json.dumps(value, separators=(",", ":"), ensure_ascii=True)
    return "RPIA1:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:20]


def assess_reactive_pair_interactions(
    molecule: str | MoleculeAnalysis,
) -> Tuple[ReactivePairInteractionAssessment, ...]:
    """Find declaratively defined incompatible sites within one component.

    Sites in different dot-separated components are intentionally not paired:
    their coexistence in a reaction mixture is a condition/selectivity question,
    not evidence that either precursor molecule is intrinsically incompatible.
    """

    analysis = analyze_molecule(molecule) if isinstance(molecule, str) else molecule
    if not analysis.valid:
        raise ValueError(analysis.structure.error or "invalid molecular structure")
    payload = load_reactive_pair_interaction_definition()
    components = {
        component.component_index: component for component in analysis.components
    }
    sites_by_component: dict[int, list[ReactiveSiteHypothesis]] = {}
    for site in analysis.reactive_site_hypotheses:
        sites_by_component.setdefault(site.component_index, []).append(site)

    assessments: list[ReactivePairInteractionAssessment] = []
    seen: set[tuple[str, int, str, str]] = set()
    for rule in payload["rules"]:
        left_selector = rule["left_site"]
        right_selector = rule["right_site"]
        for component_index, sites in sorted(sites_by_component.items()):
            left_sites = tuple(
                site for site in sites if _matches_selector(site, left_selector)
            )
            right_sites = tuple(
                site for site in sites if _matches_selector(site, right_selector)
            )
            if not left_sites or not right_sites:
                continue
            component = components[component_index]
            rd_molecule = parse_smiles(component.input_smiles)
            if rd_molecule is None:
                continue
            for left in left_sites:
                for right in right_sites:
                    if left.hypothesis_id == right.hypothesis_id:
                        continue
                    key = (
                        str(rule["rule_id"]),
                        component_index,
                        left.hypothesis_id,
                        right.hypothesis_id,
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    left_reference = _site_reference(left, left_selector)
                    right_reference = _site_reference(right, right_selector)
                    path = Chem.GetShortestPath(
                        rd_molecule,
                        left_reference.atom_index,
                        right_reference.atom_index,
                    )
                    topology = str(rule["topology_evaluation"])
                    graph_distance = (
                        len(path) - 1 if topology != "none" and path else None
                    )
                    closure_size = (
                        len(path)
                        if topology == "potential_bond_closure" and path
                        else None
                    )
                    assessment_identity = (
                        rule["rule_id"],
                        component.canonical_smiles,
                        left_reference.atom_index,
                        right_reference.atom_index,
                    )
                    assessments.append(
                        ReactivePairInteractionAssessment(
                            assessment_id=_assessment_id(assessment_identity),
                            rule_id=str(rule["rule_id"]),
                            interaction_class=str(rule["interaction_class"]),
                            scope="same_component",
                            component_index=component_index,
                            component_smiles=component.canonical_smiles,
                            left_site=left_reference,
                            right_site=right_reference,
                            graph_distance=graph_distance,
                            potential_closure_ring_size=closure_size,
                            intrinsic_severity=str(rule["intrinsic_severity"]),
                            warning_strength=str(rule["warning_strength"]),
                            message=str(rule["message"]),
                        )
                    )
    return tuple(
        sorted(
            assessments,
            key=lambda item: (
                item.component_index,
                item.rule_id,
                item.left_site.atom_index,
                item.right_site.atom_index,
                item.assessment_id,
            ),
        )
    )


__all__ = [
    "REACTIVE_PAIR_INTERACTION_DEFINITION_ID",
    "REACTIVE_PAIR_INTERACTION_SCHEMA_VERSION",
    "ReactivePairInteractionAssessment",
    "ReactivePairSiteReference",
    "assess_reactive_pair_interactions",
    "load_reactive_pair_interaction_definition",
    "validate_reactive_pair_interaction_definition",
]
