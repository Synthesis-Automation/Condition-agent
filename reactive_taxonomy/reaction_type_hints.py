"""Typed reaction-type hints projected from observed graph edits.

Hints are optional interpretations. They consume finalized observations and
molecular reactive-site annotations, and never participate in reaction identity.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Sequence

from .reaction_models import (
    REACTION_TYPE_HINT_SCHEMA_VERSION,
    ReactionComponent,
    ReactionHintParticipant,
    ReactionObservation,
    ReactionPatternMatch,
    ReactionTypeHint,
)
from .reaction_patterns import load_reaction_pattern_definitions


REACTION_TYPE_HINT_DEFINITION_ID = "reaction_type_hints.v1"
REACTION_TYPE_HINT_DEFINITION_VERSION = "reaction_type_hints.v1@1.0"
_PATH = Path(__file__).with_name("definitions") / "reaction_type_hints.v1.json"
_SITE_PRIORITY = {
    "leaving_group": 100,
    "transfer_group": 95,
    "pronucleophile": 90,
    "aromatic_CH": 88,
    "unsaturated_bond": 85,
    "electrophilic_center": 80,
    "heteroatom_bond": 75,
    "addition_donor": 70,
}


@lru_cache(maxsize=1)
def load_reaction_type_hint_definitions() -> dict[str, dict[str, Any]]:
    """Load and validate graph-pattern to type-hint projections."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = dict(json.load(handle))
    if payload.get("schema_version") != "1.0":
        raise ValueError("unsupported reaction-type-hint schema")
    if payload.get("definition_id") != REACTION_TYPE_HINT_DEFINITION_ID:
        raise ValueError("unexpected reaction-type-hint definition ID")
    known_patterns = {
        str(item["id"]) for item in load_reaction_pattern_definitions()
    }
    definitions: dict[str, dict[str, Any]] = {}
    type_ids: set[str] = set()
    for raw in payload.get("hints") or ():
        definition = dict(raw)
        pattern_id = str(definition.get("pattern_id") or "")
        type_id = str(definition.get("type_id") or "")
        if not pattern_id or pattern_id not in known_patterns:
            raise ValueError(f"unknown reaction-type-hint pattern: {pattern_id}")
        if pattern_id in definitions or not type_id or type_id in type_ids:
            raise ValueError(f"duplicate reaction type hint: {pattern_id}/{type_id}")
        if not str(definition.get("category") or "") or not str(
            definition.get("display_name") or ""
        ):
            raise ValueError(f"incomplete reaction type hint: {pattern_id}")
        roles = []
        for slot in definition.get("participant_slots") or ():
            role = str(slot.get("role") or "")
            prefixes = tuple(str(value) for value in slot.get("signature_prefixes") or ())
            if not role or not prefixes or any(not value for value in prefixes):
                raise ValueError(f"invalid participant slot: {pattern_id}")
            roles.append(role)
        if len(roles) != len(set(roles)):
            raise ValueError(f"duplicate participant role: {pattern_id}")
        definition["participant_slots"] = tuple(
            {
                "role": str(slot["role"]),
                "signature_prefixes": tuple(
                    str(value) for value in slot["signature_prefixes"]
                ),
            }
            for slot in definition.get("participant_slots") or ()
        )
        definition["qualifier_tokens"] = tuple(
            str(value) for value in definition.get("qualifier_tokens") or ()
        )
        definitions[pattern_id] = definition
        type_ids.add(type_id)
    return definitions


def validate_reaction_type_hint_definitions() -> list[str]:
    """Return definition validation errors without hiding their cause."""
    try:
        load_reaction_type_hint_definitions()
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        return [str(exc)]
    return []


def _active_atoms(
    observation: ReactionObservation,
    edit_indices: Sequence[int],
) -> dict[int, set[int]]:
    active: dict[int, set[int]] = {}
    for index in edit_indices:
        edit = observation.edits[int(index)]
        for atom in (edit.atom_1, edit.atom_2):
            if atom is not None and atom.side == "reactant":
                active.setdefault(atom.component_index, set()).add(atom.atom_index)
    return active


def _environment(component: ReactionComponent, site_id: str) -> Any | None:
    return next(
        (
            value
            for value in component.molecule_analysis.reactive_site_environments
            if value.hypothesis_id == site_id
        ),
        None,
    )


def _participant(
    component: ReactionComponent,
    site: Any,
    active: set[int],
    *,
    role: str | None,
    confidence: float,
) -> ReactionHintParticipant:
    environment = _environment(component, str(site.hypothesis_id))
    profile = environment.reactivity_profile if environment is not None else None
    attached = (
        tuple(profile.context.attached_groups)
        if profile is not None and profile.context_kind == "heteroatom"
        else ()
    )
    alpha_values = tuple(
        bool(value.alpha_branched)
        for value in attached
        if value.alpha_branched is not None
    )
    return ReactionHintParticipant(
        component_index=component.component_index,
        site_id=str(site.hypothesis_id),
        site_type=str(site.site_type),
        canonical_signature=str(site.canonical_signature),
        chemist_label=str(site.chemist_label),
        active_atom_indices=tuple(sorted(active.intersection(site.atom_indices))),
        role=role,
        role_confidence=confidence if role else 0.0,
        center_substitution_class=(
            str(profile.reactive_center.substitution_class)
            if profile is not None and profile.reactive_center.substitution_class
            else None
        ),
        attachment_carbon_classes=tuple(
            sorted(
                {
                    str(value.attachment_carbon_class)
                    for value in attached
                    if value.attachment_carbon_class
                }
            )
        ),
        alpha_branched=(any(alpha_values) if alpha_values else None),
    )


def _candidate_sites(
    component: ReactionComponent,
    active: set[int],
) -> tuple[tuple[Any, int], ...]:
    return tuple(
        (site, len(active.intersection(site.atom_indices)))
        for site in component.molecule_analysis.reactive_site_hypotheses
        if active.intersection(site.atom_indices)
    )


def _prefix_width(signature: str, prefixes: Sequence[str]) -> int:
    return max(
        (len(prefix) for prefix in prefixes if signature.startswith(prefix)),
        default=0,
    )


def _slot_participants(
    components: Sequence[ReactionComponent],
    active: Mapping[int, set[int]],
    slots: Sequence[Mapping[str, Any]],
    *,
    confidence: float,
) -> tuple[ReactionHintParticipant, ...]:
    by_index = {
        component.component_index: component
        for component in components
        if component.component_index in active
    }
    if len(by_index) < len(slots):
        return ()
    assignments = []
    for component_order in itertools.permutations(sorted(by_index), len(slots)):
        selected = []
        total_prefix = total_overlap = total_priority = 0
        for slot, component_index in zip(slots, component_order):
            component = by_index[component_index]
            prefixes = tuple(slot.get("signature_prefixes") or ())
            candidates = []
            for site, overlap in _candidate_sites(component, active[component_index]):
                prefix_width = _prefix_width(str(site.canonical_signature), prefixes)
                if prefix_width:
                    candidates.append(
                        (
                            -prefix_width,
                            -overlap,
                            -_SITE_PRIORITY.get(str(site.site_type), 0),
                            len(site.atom_indices),
                            str(site.hypothesis_id),
                            site,
                        )
                    )
            if not candidates:
                selected = []
                break
            choice = min(candidates)
            site = choice[-1]
            total_prefix += -choice[0]
            total_overlap += -choice[1]
            total_priority += -choice[2]
            selected.append(
                _participant(
                    component,
                    site,
                    active[component_index],
                    role=str(slot["role"]),
                    confidence=confidence,
                )
            )
        if selected:
            assignments.append(
                (
                    -total_prefix,
                    -total_overlap,
                    -total_priority,
                    tuple(item.component_index for item in selected),
                    tuple(selected),
                )
            )
    return min(assignments)[-1] if assignments else ()


def _generic_participants(
    components: Sequence[ReactionComponent],
    active: Mapping[int, set[int]],
) -> tuple[ReactionHintParticipant, ...]:
    participants = []
    for component in components:
        atoms = active.get(component.component_index)
        if not atoms:
            continue
        candidates = _candidate_sites(component, atoms)
        if not candidates:
            continue
        site, _ = min(
            candidates,
            key=lambda item: (
                -item[1],
                -_SITE_PRIORITY.get(str(item[0].site_type), 0),
                len(item[0].atom_indices),
                str(item[0].hypothesis_id),
            ),
        )
        participants.append(
            _participant(component, site, atoms, role=None, confidence=0.0)
        )
    return tuple(sorted(participants, key=lambda item: item.component_index))


def _digest(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return f"RTH1:{hashlib.sha256(encoded).hexdigest()[:24]}"


def build_reaction_type_hints(
    observation: ReactionObservation,
    reactants: Sequence[ReactionComponent],
    pattern_matches: Sequence[ReactionPatternMatch],
) -> tuple[ReactionTypeHint, ...]:
    """Project supported pattern matches onto graph-bound active sites."""
    definitions = load_reaction_type_hint_definitions()
    hints = []
    scope = observation.topology.reaction_scope if observation.topology else "unresolved"
    for match in pattern_matches:
        definition = definitions.get(match.pattern_id)
        if definition is None:
            continue
        active = _active_atoms(observation, match.matched_edit_indices)
        slots = tuple(definition.get("participant_slots") or ())
        participants = (
            _slot_participants(
                reactants,
                active,
                slots,
                confidence=match.confidence,
            )
            if slots
            else _generic_participants(reactants, active)
        )
        warnings = (
            ("ACTIVE_HINT_PARTICIPANTS_UNRESOLVED",)
            if slots and len(participants) != len(slots)
            else ()
        )
        qualifiers = tuple(
            sorted(
                set(definition.get("qualifier_tokens") or ())
                | {
                    f"site:{participant.canonical_signature}"
                    for participant in participants
                }
            )
        )
        identity = {
            "definition": REACTION_TYPE_HINT_DEFINITION_VERSION,
            "pattern_id": match.pattern_id,
            "type_id": definition["type_id"],
            "scope": scope,
            "edits": match.matched_edit_indices,
            "participants": [
                {
                    "signature": participant.canonical_signature,
                    "role": participant.role,
                    "site_type": participant.site_type,
                    "active_atoms": participant.active_atom_indices,
                }
                for participant in participants
            ],
        }
        hints.append(
            ReactionTypeHint(
                hint_id=_digest(identity),
                type_id=str(definition["type_id"]),
                category=str(definition["category"]),
                display_name=str(definition["display_name"]),
                pattern_id=match.pattern_id,
                confidence=match.confidence,
                reaction_scope=scope,
                matched_edit_indices=match.matched_edit_indices,
                matched_core_event_ids=match.matched_core_event_ids,
                participants=participants,
                qualifier_tokens=qualifiers,
                compatible_named_families=match.compatible_named_families,
                requires_condition_evidence=match.requires_condition_evidence,
                evidence=tuple(
                    dict.fromkeys(
                        (*match.evidence, "reaction_pattern", "active_site_projection")
                    )
                ),
                warnings=warnings,
            )
        )
    return tuple(hints)


__all__ = [
    "REACTION_TYPE_HINT_DEFINITION_ID",
    "REACTION_TYPE_HINT_DEFINITION_VERSION",
    "build_reaction_type_hints",
    "load_reaction_type_hint_definitions",
    "validate_reaction_type_hint_definitions",
]
