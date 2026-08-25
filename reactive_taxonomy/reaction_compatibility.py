"""Graph-derived reaction-regime compatibility diagnostics.

These diagnostics do not assign a named reaction from reagent strings. They
match validated or inferred edit motifs and unchanged molecular spectators,
then report the regime conflict encoded in a versioned definition.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence, Tuple

from .reaction_api import featurize_reaction
from .molecular_motifs import load_molecular_motif_definitions
from .reaction_models import (
    ReactionAnalysis,
    ReactionAtomReference,
    ReactionEdit,
)


REACTION_COMPATIBILITY_DEFINITION_ID = "reaction_compatibility.v1"
REACTION_COMPATIBILITY_SCHEMA_VERSION = "1.0"
REACTION_COMPATIBILITY_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "reaction_compatibility.v1.json"
)

_EDIT_MOTIFS = frozenset({"halo_carbon_transfer_to_carbonyl"})
_SEVERITIES = frozenset({"low", "medium", "high", "critical"})
_WARNING_STRENGTHS = frozenset({"advisory", "strong"})
_HALOGENS = frozenset({"Cl", "Br", "I"})


@dataclass(frozen=True)
class ReactionCompatibilityAssessment:
    """One declarative conflict between an edit-derived regime and spectator."""

    assessment_id: str
    rule_id: str
    interaction_class: str
    required_edit_motif: str
    inferred_regime: str
    spectator_group_id: str
    spectator_component_index: int
    spectator_tags: Tuple[str, ...]
    matched_edit_indices: Tuple[int, ...]
    intrinsic_severity: str
    warning_strength: str
    evidence_quality: str
    confidence: float
    message: str
    definition_id: str = REACTION_COMPATIBILITY_DEFINITION_ID
    schema_version: str = REACTION_COMPATIBILITY_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible assessment."""

        return asdict(self)


def validate_reaction_compatibility_definition(
    payload: Mapping[str, Any],
) -> list[str]:
    """Return deterministic validation errors for reaction compatibility rules."""

    errors: list[str] = []
    if payload.get("definition_id") != REACTION_COMPATIBILITY_DEFINITION_ID:
        errors.append("unexpected_reaction_compatibility_definition_id")
    if payload.get("schema_version") != REACTION_COMPATIBILITY_SCHEMA_VERSION:
        errors.append("unsupported_reaction_compatibility_schema")
    rules = payload.get("rules")
    if not isinstance(rules, list) or not rules:
        return [*errors, "reaction_compatibility_rules_must_be_nonempty"]
    known_tags = {
        str(tag)
        for definition in load_molecular_motif_definitions()
        for tag in definition.get("tags") or ()
    }
    rule_ids: list[str] = []
    allowed_keys = {
        "rule_id",
        "interaction_class",
        "required_edit_motif",
        "spectator_tags_any",
        "inferred_regime",
        "intrinsic_severity",
        "warning_strength",
        "message",
    }
    for index, rule in enumerate(rules):
        location = f"reaction_compatibility_rule_{index}"
        if not isinstance(rule, Mapping):
            errors.append(f"{location}_must_be_object")
            continue
        unknown = set(rule) - allowed_keys
        if unknown:
            errors.append(f"{location}_unsupported_keys:{sorted(unknown)}")
        rule_id = str(rule.get("rule_id") or "")
        rule_ids.append(rule_id)
        if not rule_id:
            errors.append(f"{location}_missing_rule_id")
        if not str(rule.get("interaction_class") or ""):
            errors.append(f"{location}_missing_interaction_class")
        motif = str(rule.get("required_edit_motif") or "")
        if motif not in _EDIT_MOTIFS:
            errors.append(f"{location}_unsupported_edit_motif:{motif}")
        tags = tuple(str(value) for value in rule.get("spectator_tags_any") or ())
        if not tags:
            errors.append(f"{location}_requires_spectator_tags")
        unknown_tags = set(tags) - known_tags
        if unknown_tags:
            errors.append(f"{location}_unknown_spectator_tags:{sorted(unknown_tags)}")
        if not str(rule.get("inferred_regime") or ""):
            errors.append(f"{location}_missing_inferred_regime")
        severity = str(rule.get("intrinsic_severity") or "")
        if severity not in _SEVERITIES:
            errors.append(f"{location}_invalid_severity:{severity}")
        strength = str(rule.get("warning_strength") or "")
        if strength not in _WARNING_STRENGTHS:
            errors.append(f"{location}_invalid_warning_strength:{strength}")
        if not str(rule.get("message") or "").strip():
            errors.append(f"{location}_missing_message")
    if len(rule_ids) != len(set(rule_ids)):
        errors.append("duplicate_reaction_compatibility_rule_id")
    return errors


@lru_cache(maxsize=1)
def load_reaction_compatibility_definition() -> Mapping[str, Any]:
    """Load and validate the canonical reaction compatibility rules."""

    payload = json.loads(
        REACTION_COMPATIBILITY_DEFINITION_PATH.read_text(encoding="utf-8")
    )
    errors = validate_reaction_compatibility_definition(payload)
    if errors:
        raise ValueError(
            "invalid reaction-compatibility definition: " + ", ".join(errors)
        )
    return payload


def _atom_key(
    atom: ReactionAtomReference,
) -> tuple[str, int, int, int | None]:
    return (
        atom.side,
        atom.component_index,
        atom.atom_index,
        atom.atom_map_number,
    )


def _bond_atoms(
    edit: ReactionEdit,
    elements: frozenset[str],
) -> dict[str, ReactionAtomReference] | None:
    if edit.atom_2 is None:
        return None
    atoms = (edit.atom_1, edit.atom_2)
    if frozenset(atom.element for atom in atoms) != elements:
        return None
    return {atom.element: atom for atom in atoms}


def _halo_carbonyl_transfer_edits(
    edits: Sequence[ReactionEdit],
) -> tuple[tuple[int, int, int], ...]:
    broken: list[tuple[int, ReactionAtomReference]] = []
    carbonyl: list[tuple[int, ReactionAtomReference]] = []
    formed: list[
        tuple[int, tuple[ReactionAtomReference, ReactionAtomReference]]
    ] = []
    for index, edit in enumerate(edits):
        if edit.edit_type == "broken" and edit.old_order == "SINGLE":
            if edit.atom_2 is None:
                continue
            elements = frozenset((edit.atom_1.element, edit.atom_2.element))
            if "C" in elements and elements.intersection(_HALOGENS):
                carbon = edit.atom_1 if edit.atom_1.element == "C" else edit.atom_2
                broken.append((index, carbon))
        elif (
            edit.edit_type == "order_changed"
            and edit.old_order == "DOUBLE"
            and edit.new_order == "SINGLE"
        ):
            atoms = _bond_atoms(edit, frozenset({"C", "O"}))
            if atoms is not None:
                carbonyl.append((index, atoms["C"]))
        elif (
            edit.edit_type == "formed"
            and edit.old_order is None
            and edit.new_order == "SINGLE"
            and edit.atom_2 is not None
            and edit.atom_1.element == "C"
            and edit.atom_2.element == "C"
        ):
            formed.append((index, (edit.atom_1, edit.atom_2)))
    matches = []
    for broken_index, halo_carbon in broken:
        for carbonyl_index, carbonyl_carbon in carbonyl:
            expected = {_atom_key(halo_carbon), _atom_key(carbonyl_carbon)}
            for formed_index, formed_atoms in formed:
                if {_atom_key(atom) for atom in formed_atoms} == expected:
                    matches.append(
                        tuple(sorted((broken_index, carbonyl_index, formed_index)))
                    )
    return tuple(sorted(set(matches)))


def _assessment_id(value: Sequence[Any]) -> str:
    encoded = json.dumps(value, separators=(",", ":"), ensure_ascii=True)
    return "RCIA1:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:20]


def _assess_analysis(
    analysis: ReactionAnalysis,
) -> tuple[ReactionCompatibilityAssessment, ...]:
    if not analysis.valid or analysis.reaction_signature is None:
        return ()
    edit_matches = _halo_carbonyl_transfer_edits(analysis.reaction_signature.edits)
    if not edit_matches:
        return ()
    payload = load_reaction_compatibility_definition()
    assessments: list[ReactionCompatibilityAssessment] = []
    for rule in payload["rules"]:
        if rule["required_edit_motif"] != "halo_carbon_transfer_to_carbonyl":
            continue
        required_tags = set(rule["spectator_tags_any"])
        for spectator in analysis.spectator_groups:
            if not required_tags.intersection(spectator.tags):
                continue
            matched_indices = tuple(
                sorted({index for match in edit_matches for index in match})
            )
            identity = (
                analysis.reaction_signature.signature_id,
                rule["rule_id"],
                spectator.group_id,
                spectator.component_index,
                spectator.atom_indices,
                matched_indices,
            )
            edit_confidences = tuple(
                analysis.reaction_signature.edits[index].confidence
                for index in matched_indices
            )
            assessments.append(
                ReactionCompatibilityAssessment(
                    assessment_id=_assessment_id(identity),
                    rule_id=str(rule["rule_id"]),
                    interaction_class=str(rule["interaction_class"]),
                    required_edit_motif=str(rule["required_edit_motif"]),
                    inferred_regime=str(rule["inferred_regime"]),
                    spectator_group_id=spectator.group_id,
                    spectator_component_index=spectator.component_index,
                    spectator_tags=tuple(spectator.tags),
                    matched_edit_indices=matched_indices,
                    intrinsic_severity=str(rule["intrinsic_severity"]),
                    warning_strength=str(rule["warning_strength"]),
                    evidence_quality=analysis.evidence_quality,
                    confidence=min(edit_confidences) if edit_confidences else 0.0,
                    message=str(rule["message"]),
                )
            )
    return tuple(
        sorted(
            assessments,
            key=lambda item: (
                item.rule_id,
                item.spectator_component_index,
                item.spectator_group_id,
                item.assessment_id,
            ),
        )
    )


@lru_cache(maxsize=20_000)
def _assess_reaction_smiles(
    reaction_smiles: str,
) -> tuple[ReactionCompatibilityAssessment, ...]:
    return _assess_analysis(featurize_reaction(reaction_smiles))


def assess_reaction_compatibility(
    reaction: str | ReactionAnalysis,
) -> tuple[ReactionCompatibilityAssessment, ...]:
    """Assess edit-derived reaction regimes against unchanged spectators."""

    if isinstance(reaction, str):
        return _assess_reaction_smiles(reaction)
    return _assess_analysis(reaction)


__all__ = [
    "REACTION_COMPATIBILITY_DEFINITION_ID",
    "REACTION_COMPATIBILITY_SCHEMA_VERSION",
    "ReactionCompatibilityAssessment",
    "assess_reaction_compatibility",
    "load_reaction_compatibility_definition",
    "validate_reaction_compatibility_definition",
]
