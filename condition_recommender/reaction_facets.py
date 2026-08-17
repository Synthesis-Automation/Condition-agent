"""Generic graph-derived reaction facets for progressive retrieval."""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Tuple


_RULES_PATH = (
    Path(__file__).with_name("definitions")
    / "reaction_facet_retrieval.v1.json"
)
_REMOTE_CLASSES = {
    "aryl",
    "heteroaryl",
    "alkyl",
    "ring_aliphatic",
    "alkenyl",
    "alkynyl",
    "acyl",
    "heteroatom",
    "generic_R",
}


@lru_cache(maxsize=1)
def load_reaction_facet_rules() -> dict[str, Any]:
    """Load and validate the versioned reaction-facet definition."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if rules.get("schema_version") != "1.1":
        raise ValueError("unsupported reaction-facet definition schema")
    if rules.get("definition_id") != "reaction_facet_retrieval.v1":
        raise ValueError("unexpected reaction-facet definition ID")
    if not str(rules.get("calibration_status") or "").strip():
        raise ValueError("reaction-facet definition requires calibration status")
    state_fields = tuple(str(value) for value in rules.get("state_fields") or ())
    hard_state_fields = tuple(
        str(value)
        for value in rules.get("hard_compatibility_state_fields") or ()
    )
    topology_fields = tuple(
        str(value) for value in rules.get("topology_fields") or ()
    )
    if not state_fields or len(set(state_fields)) != len(state_fields):
        raise ValueError("reaction-facet state fields are invalid")
    if (
        not hard_state_fields
        or len(set(hard_state_fields)) != len(hard_state_fields)
        or not set(hard_state_fields) <= set(state_fields)
    ):
        raise ValueError("reaction-facet hard compatibility fields are invalid")
    if rules.get("hard_compatibility_requires_hydrogen_change") is not True:
        raise ValueError(
            "reaction-facet hard compatibility must require hydrogen change"
        )
    if not topology_fields or len(set(topology_fields)) != len(topology_fields):
        raise ValueError("reaction-facet topology fields are invalid")
    parents = rules.get("attachment_parents")
    if not isinstance(parents, Mapping) or set(parents) != _REMOTE_CLASSES:
        raise ValueError("reaction-facet attachment hierarchy is incomplete")
    prefixes = tuple(
        str(value) for value in rules.get("active_site_signature_prefixes") or ()
    )
    if not prefixes or any(not value for value in prefixes):
        raise ValueError("reaction-facet active-site prefixes are invalid")
    ladder = tuple(str(value) for value in rules.get("retrieval_ladder") or ())
    if ladder != (
        "reaction_facet_exact",
        "reaction_facet_attachment_relaxed",
    ):
        raise ValueError("reaction-facet retrieval ladder is invalid")
    rules["state_fields"] = state_fields
    rules["hard_compatibility_state_fields"] = hard_state_fields
    rules["topology_fields"] = topology_fields
    rules["active_site_signature_prefixes"] = prefixes
    rules["retrieval_ladder"] = ladder
    rules["attachment_parents"] = {
        str(key): str(value) for key, value in parents.items()
    }
    return rules


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def _digest(prefix: str, payload: Any) -> str:
    encoded = _canonical_json(payload).encode("utf-8")
    return f"{prefix}:" + hashlib.sha256(encoded).hexdigest()[:32]


def _state_payload(state: Any, fields: Tuple[str, ...]) -> Any:
    if not isinstance(state, Mapping):
        return None
    return tuple((field, state.get(field)) for field in fields)


def active_atom_state_tokens(
    reaction_core: Mapping[str, Any] | None,
) -> Tuple[str, ...]:
    """Return exact before/after states for atoms that lose or gain hydrogen.

    These atoms encode nucleophile substitution level and electronic class.
    Atom-map identifiers and remote substituent identity are deliberately
    excluded so equivalent chemistry remains representation invariant.
    """

    if not reaction_core:
        return ()
    fields = load_reaction_facet_rules()["hard_compatibility_state_fields"]
    tokens = []
    for transition in reaction_core.get("atom_transitions") or ():
        if not isinstance(transition, Mapping):
            continue
        if int(transition.get("incident_edit_count") or 0) < 1:
            continue
        before = transition.get("before_state")
        after = transition.get("after_state")
        if not isinstance(before, Mapping) or not isinstance(after, Mapping):
            continue
        before_h = before.get("total_hydrogens")
        after_h = after.get("total_hydrogens")
        if (
            not isinstance(before_h, int)
            or not isinstance(after_h, int)
            or before_h == after_h
        ):
            continue
        tokens.append(
            _canonical_json(
                {
                    "before": _state_payload(before, fields),
                    "after": _state_payload(after, fields),
                }
            )
        )
    return tuple(sorted(tokens))


def active_atom_states_compatible(
    query_core: Mapping[str, Any] | None,
    precedent_core: Mapping[str, Any] | None,
) -> bool:
    """Reject a precedent only when both cores contradict at active X-H atoms."""

    query = active_atom_state_tokens(query_core)
    precedent = active_atom_state_tokens(precedent_core)
    return not query or not precedent or query == precedent


def _active_site_signatures(
    transition: Mapping[str, Any],
    fallback_descriptor: Mapping[str, Any],
    prefixes: Tuple[str, ...],
) -> Tuple[str, ...]:
    """Return an unambiguous X-H site class consistent with one atom change."""
    before = transition.get("before_state") or {}
    after = transition.get("after_state") or {}
    if not isinstance(before, Mapping) or not isinstance(after, Mapping):
        return ()
    element = str(before.get("element") or "")
    before_h = before.get("total_hydrogens")
    after_h = after.get("total_hydrogens")
    if (
        not element
        or not isinstance(before_h, int)
        or not isinstance(after_h, int)
        or before_h <= after_h
    ):
        return ()
    candidates = set()
    for raw_token in fallback_descriptor.get("reactant_site_tokens") or ():
        token = str(raw_token)
        if not any(token.startswith(prefix) for prefix in prefixes):
            continue
        parts = token.split("|")
        if len(parts) < 4:
            continue
        center = parts[1]
        hydrogen = parts[2]
        if center.startswith(element) and hydrogen == f"H{before_h}":
            candidates.add(token)
    return tuple(candidates) if len(candidates) == 1 else ()


def _attachment_tokens(
    core: Mapping[str, Any],
    transition: Mapping[str, Any],
    *,
    relaxed: bool,
    parents: Mapping[str, str],
) -> Tuple[Tuple[str, str], ...]:
    atom_map_number = transition.get("atom_map_number")
    if atom_map_number is None:
        before = transition.get("before_state") or {}
        if isinstance(before, Mapping):
            atom_map_number = before.get("atom_map_number")
    if atom_map_number is None:
        return ()
    tokens = set()
    for remote in core.get("remote_subgraphs") or ():
        if not isinstance(remote, Mapping) or remote.get("side") != "reactant":
            continue
        continuity = str(remote.get("continuity") or "unresolved")
        for port in remote.get("attachment_ports") or ():
            if not isinstance(port, Mapping):
                continue
            if port.get("core_atom_map_number") != atom_map_number:
                continue
            profile = port.get("substituent_profile") or {}
            base_class = (
                str(profile.get("base_class") or "")
                if isinstance(profile, Mapping)
                else ""
            )
            if not base_class:
                base_class = str(remote.get("remote_class") or "")
            if not base_class:
                continue
            value = parents.get(base_class, base_class) if relaxed else base_class
            tokens.add((continuity, value))
    return tuple(sorted(tokens))


def _facet_payload(
    signature: Mapping[str, Any],
    core: Mapping[str, Any],
    fallback_descriptor: Mapping[str, Any],
    *,
    relaxed: bool,
) -> Mapping[str, Any] | None:
    rules = load_reaction_facet_rules()
    bond_key = str(signature.get("bond_edit_signature_key") or "")
    transitions = tuple(
        value
        for value in core.get("atom_transitions") or ()
        if isinstance(value, Mapping)
    )
    if not bond_key or not transitions:
        return None
    state_fields = rules["state_fields"]
    transition_tokens = []
    for transition in transitions:
        before = _state_payload(transition.get("before_state"), state_fields)
        after = _state_payload(transition.get("after_state"), state_fields)
        if before is None and after is None:
            continue
        transition_tokens.append(
            {
                "role": str(transition.get("role") or "participant"),
                "before": before,
                "after": after,
                "attachments": _attachment_tokens(
                    core,
                    transition,
                    relaxed=relaxed,
                    parents=rules["attachment_parents"],
                ),
                "active_site_signatures": _active_site_signatures(
                    transition,
                    fallback_descriptor,
                    rules["active_site_signature_prefixes"],
                ),
            }
        )
    if not transition_tokens:
        return None
    transition_tokens.sort(key=_canonical_json)
    topology = signature.get("topology") or {}
    topology_payload = tuple(
        (field, topology.get(field))
        for field in rules["topology_fields"]
        if isinstance(topology, Mapping)
    )
    return {
        "definition_id": rules["definition_id"],
        "definition_version": rules["schema_version"],
        "bond_edit_signature_key": bond_key,
        "event_count": int(signature.get("event_count") or 0),
        "event_scope": str(signature.get("event_scope") or "unresolved"),
        "topology": topology_payload,
        "transitions": transition_tokens,
    }


def reaction_facet_keys(
    signature: Mapping[str, Any],
    reaction_core: Mapping[str, Any] | None,
    fallback_descriptor: Mapping[str, Any] | None = None,
) -> dict[str, str]:
    """Return exact-class and attachment-relaxed generic retrieval keys."""
    if not reaction_core:
        return {}
    fallback = fallback_descriptor or {}
    exact = _facet_payload(signature, reaction_core, fallback, relaxed=False)
    relaxed = _facet_payload(signature, reaction_core, fallback, relaxed=True)
    if exact is None or relaxed is None:
        return {}
    return {
        "reaction_facet_exact": _digest("RFC1", exact),
        "reaction_facet_attachment_relaxed": _digest("RFR1", relaxed),
    }


__all__ = [
    "active_atom_state_tokens",
    "active_atom_states_compatible",
    "load_reaction_facet_rules",
    "reaction_facet_keys",
]
