"""Evidence-based contextual condition-role resolution."""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

from .api import resolve_identifier, resolve_substance, resolve_substance_id
from .loader import load_role_definitions
from .models import (
    CONDITION_IDENTIFIER_TYPES,
    ContextualRoleAssignment,
    ResolvedConditionComponent,
)

_RULES_PATH = Path(__file__).with_name("definitions") / "role_resolution.v2.json"


@lru_cache(maxsize=1)
def load_role_resolution_rules() -> Dict[str, Any]:
    """Load role resolution and recipe placement rules."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "2.0":
        raise ValueError("Unsupported role resolution schema")
    known_roles = {
        str(item["id"]) for item in load_role_definitions().get("roles", ())
    }
    referenced = set(rules.get("role_buckets") or ())
    for values in rules.get("source_role_preferences", {}).values():
        referenced.update(str(value) for value in values)
    for values in rules.get("source_role_hint_preferences", {}).values():
        referenced.update(str(value) for value in values)
    referenced.update(
        str(value)
        for value in rules.get("source_role_fallbacks", {}).values()
        if value
    )
    unknown = referenced - known_roles
    if unknown:
        raise ValueError(f"Role resolution references unknown roles: {sorted(unknown)}")
    return rules


@lru_cache(maxsize=1)
def _role_priorities() -> Dict[str, int]:
    return {
        str(item["id"]): int(item.get("priority", 100))
        for item in load_role_definitions().get("roles", ())
    }


def _ranked_roles(
    capabilities: Iterable[Any],
    *,
    source_field: str,
    source_role_hint: Optional[str],
    preferred_roles: Iterable[str],
) -> Tuple[ContextualRoleAssignment, ...]:
    rules = load_role_resolution_rules()
    source_preferences = set(
        rules.get("source_role_preferences", {}).get(source_field, ())
    )
    hint_preferences = set(
        rules.get("source_role_hint_preferences", {}).get(
            source_role_hint or "", ()
        )
    )
    context_preferences = {str(value) for value in preferred_roles}
    ranked = []
    for capability in capabilities:
        role_id = capability.role_id
        hint_match = role_id in hint_preferences
        source_match = role_id in source_preferences
        context_match = role_id in context_preferences
        evidence = ["curated_role_capability"]
        if hint_match:
            confidence = 0.95
            evidence.append("source_role_hint_match")
        elif source_match and context_match:
            confidence = 0.90
            evidence.extend(("source_field_role_match", "reaction_role_preference"))
        elif source_match:
            confidence = 0.85
            evidence.append("source_field_role_match")
        elif context_match:
            confidence = 0.75
            evidence.append("reaction_role_preference")
        else:
            confidence = 0.65
        ranked.append(
            ContextualRoleAssignment(
                role_id=role_id,
                confidence=confidence,
                evidence=tuple(evidence),
            )
        )
    return tuple(
        sorted(
            ranked,
            key=lambda item: (
                -item.confidence,
                _role_priorities().get(item.role_id, 100),
                item.role_id,
            ),
        )
    )


def _role_decision(
    roles: Tuple[ContextualRoleAssignment, ...],
    *,
    source_field: str,
    source_role_hint: Optional[str],
) -> Tuple[str, Optional[str], Optional[float]]:
    if not roles:
        return "unassigned", None, None
    rules = load_role_resolution_rules()
    hint_preferences = set(
        rules.get("source_role_hint_preferences", {}).get(
            source_role_hint or "", ()
        )
    )
    source_preferences = set(
        rules.get("source_role_preferences", {}).get(source_field, ())
    )
    available = {role.role_id for role in roles}
    if hint_preferences and not available.intersection(hint_preferences):
        return "conflicting", None, None
    if source_preferences and not available.intersection(source_preferences):
        return "conflicting", None, None
    top_confidence = roles[0].confidence
    top_roles = tuple(role for role in roles if role.confidence == top_confidence)
    if len(top_roles) != 1:
        return "ambiguous", None, None
    primary = top_roles[0]
    return "assigned", primary.role_id, primary.confidence


def _resolution_provenance(
    *,
    base: Optional[Dict[str, Any]],
    result: Any,
    resolved_identifier_type: str,
    preferred_roles: Iterable[str],
    rules: Dict[str, Any],
) -> Dict[str, Any]:
    matched = result.matched_identifier
    return {
        **dict(base or {}),
        "identity_match_kind": result.match_kind,
        "identity_identifier_id": matched.identifier_id if matched else None,
        "identity_identifier_type": (
            matched.identifier_type if matched else resolved_identifier_type
        ),
        "identifier_type": resolved_identifier_type,
        "preferred_roles": tuple(sorted({str(value) for value in preferred_roles})),
        "definition_id": rules["definition_id"],
    }


def resolve_contextual_component(
    identifier: str,
    *,
    source_field: str,
    identifier_type: str = "auto",
    source_role_hint: Optional[str] = None,
    preferred_roles: Iterable[str] = (),
    amount: Optional[float] = None,
    amount_unit: Optional[str] = None,
    provenance: Optional[Dict[str, Any]] = None,
) -> ResolvedConditionComponent:
    """Resolve identity and select a role only when evidence is decisive."""
    rules = load_role_resolution_rules()
    preferred_roles = tuple(sorted({str(value) for value in preferred_roles}))
    supported_identifier_types = {
        "auto",
        "name",
        "substance_id",
        *CONDITION_IDENTIFIER_TYPES,
    }
    if identifier_type not in supported_identifier_types:
        raise ValueError(f"Unsupported condition identifier type: {identifier_type}")
    resolved_identifier_type = identifier_type
    if resolved_identifier_type == "auto":
        resolved_identifier_type = (
            "cas" if re.fullmatch(r"\d{2,7}-\d{2}-\d", identifier.strip()) else "name"
        )
    if resolved_identifier_type == "cas":
        result = resolve_substance(cas=identifier)
    elif resolved_identifier_type == "substance_id":
        result = resolve_substance_id(identifier)
    elif resolved_identifier_type == "name":
        result = resolve_substance(name=identifier)
    else:
        result = resolve_identifier(identifier, identifier_type=resolved_identifier_type)

    if result.status == "resolved" and result.substance is not None:
        substance = result.substance
        roles = _ranked_roles(
            substance.roles,
            source_field=source_field,
            source_role_hint=source_role_hint,
            preferred_roles=preferred_roles,
        )
        role_status, primary_role, primary_confidence = _role_decision(
            roles,
            source_field=source_field,
            source_role_hint=source_role_hint,
        )
        warnings = []
        if role_status == "unassigned":
            warnings.append("REGISTRY_IDENTITY_WITHOUT_ROLE")
        elif role_status == "ambiguous":
            warnings.append("MULTIPLE_POSSIBLE_ROLES")
        elif role_status == "conflicting":
            warnings.append("SOURCE_ROLE_CONFLICT")
        return ResolvedConditionComponent(
            raw_identifier=identifier,
            source_field=source_field,
            identity_status=result.status,
            substance_id=substance.substance_id,
            canonical_name=substance.canonical_name,
            roles=roles,
            role_status=role_status,
            primary_role=primary_role,
            primary_role_confidence=primary_confidence,
            cas=substance.cas,
            amount=amount,
            amount_unit=amount_unit,
            source_role_hint=source_role_hint,
            warnings=tuple(warnings),
            provenance=_resolution_provenance(
                base=provenance,
                result=result,
                resolved_identifier_type=resolved_identifier_type,
                preferred_roles=preferred_roles,
                rules=rules,
            ),
        )

    hint_roles = tuple(
        str(value)
        for value in rules.get("source_role_hint_preferences", {}).get(
            source_role_hint or "", ()
        )
    )
    fallback_role = rules.get("source_role_fallbacks", {}).get(source_field)
    candidate_ids = hint_roles or ((str(fallback_role),) if fallback_role else ())
    roles = tuple(
        ContextualRoleAssignment(
            role_id=role_id,
            confidence=0.55 if hint_roles else 0.50,
            evidence=(
                "source_role_hint" if hint_roles else "source_field_role",
                f"identity_{result.status}",
            ),
        )
        for role_id in candidate_ids
    )
    if len(roles) == 1:
        role_status = "assigned"
        primary_role = roles[0].role_id
        primary_confidence = roles[0].confidence
    elif roles:
        role_status = "ambiguous"
        primary_role = None
        primary_confidence = None
    else:
        role_status = "unassigned"
        primary_role = None
        primary_confidence = None
    return ResolvedConditionComponent(
        raw_identifier=identifier,
        source_field=source_field,
        identity_status=result.status,
        substance_id=None,
        canonical_name=None,
        roles=roles,
        role_status=role_status,
        primary_role=primary_role,
        primary_role_confidence=primary_confidence,
        amount=amount,
        amount_unit=amount_unit,
        source_role_hint=source_role_hint,
        warnings=("CONDITION_IDENTITY_UNCERTAINTY",),
        provenance=_resolution_provenance(
            base=provenance,
            result=result,
            resolved_identifier_type=resolved_identifier_type,
            preferred_roles=preferred_roles,
            rules=rules,
        ),
    )


__all__ = ["load_role_resolution_rules", "resolve_contextual_component"]
