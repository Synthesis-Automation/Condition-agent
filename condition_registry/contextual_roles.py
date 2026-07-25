"""Contextual condition-role resolution with ranked evidence."""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

from .api import resolve_substance, resolve_substance_id
from .loader import load_taxonomy
from .models import ContextualRoleAssignment, ResolvedConditionComponent

_RULES_PATH = Path(__file__).with_name("definitions") / "role_resolution.v1.json"


@lru_cache(maxsize=1)
def load_role_resolution_rules() -> Dict[str, Any]:
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        return dict(json.load(handle))


@lru_cache(maxsize=1)
def _role_priorities() -> Dict[str, int]:
    return {
        str(item["id"]): int(item.get("priority", 100))
        for item in load_taxonomy().get("roles", [])
    }


def _ranked_roles(
    roles: Iterable[Any],
    *,
    source_field: str,
    source_role_hint: Optional[str],
    transformation_class: Optional[str],
) -> Tuple[ContextualRoleAssignment, ...]:
    rules = load_role_resolution_rules()
    source_preferences = tuple(
        rules.get("source_role_preferences", {}).get(source_field, ())
    )
    transformation_preferences = tuple(
        rules.get("transformation_role_preferences", {}).get(
            transformation_class or "", ()
        )
    )
    hint_preferences = tuple(
        rules.get("source_role_hint_preferences", {}).get(
            source_role_hint or "", ()
        )
    )
    generic_roles = set(rules.get("generic_roles", ()))
    assignments = tuple(roles)
    if any(assignment.role_id not in generic_roles for assignment in assignments):
        assignments = tuple(
            assignment
            for assignment in assignments
            if assignment.role_id not in generic_roles
        )
    ranked = []
    for assignment in assignments:
        role_id = assignment.role_id
        source_match = role_id in source_preferences
        transformation_match = role_id in transformation_preferences
        hint_match = role_id in hint_preferences
        specific = role_id not in generic_roles
        confidence = 0.55
        evidence = ["curated_registry_role"]
        if specific:
            confidence += 0.15
        if source_match:
            confidence += 0.2
            evidence.append("source_field_role_match")
        if transformation_match:
            confidence += 0.08
            evidence.append("transformation_role_preference")
        if hint_match:
            confidence += 0.12
            evidence.append("source_role_hint_match")
        ranked.append(
            ContextualRoleAssignment(
                role_id=role_id,
                family_id=assignment.family_id,
                confidence=round(min(confidence, 0.98), 3),
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
                item.family_id or "",
            ),
        )
    )


def resolve_contextual_component(
    identifier: str,
    *,
    source_field: str,
    identifier_type: str = "auto",
    source_role_hint: Optional[str] = None,
    transformation_class: Optional[str] = None,
    named_family: Optional[str] = None,
    amount: Optional[float] = None,
    amount_unit: Optional[str] = None,
    provenance: Optional[Dict[str, Any]] = None,
) -> ResolvedConditionComponent:
    """Resolve identity and rank possible roles using reaction context."""
    rules = load_role_resolution_rules()
    if identifier_type not in {"auto", "cas", "name", "substance_id"}:
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
    else:
        result = resolve_substance(name=identifier)
    warnings = []
    if result.status == "resolved" and result.substance is not None:
        substance = result.substance
        roles = _ranked_roles(
            substance.roles,
            source_field=source_field,
            source_role_hint=source_role_hint,
            transformation_class=transformation_class,
        )
        if not roles:
            fallback = str(
                rules.get("source_fallback_roles", {}).get(
                    source_field, "other_reagent"
                )
            )
            roles = (
                ContextualRoleAssignment(
                    fallback,
                    None,
                    0.45,
                    ("source_field_fallback", "registry_identity_without_role"),
                ),
            )
            warnings.append("REGISTRY_IDENTITY_WITHOUT_ROLE")
        source_preferences = set(
            rules.get("source_role_preferences", {}).get(source_field, ())
        )
        source_preferences.update(
            rules.get("source_role_hint_preferences", {}).get(
                source_role_hint or "", ()
            )
        )
        if source_preferences and not any(
            role.role_id in source_preferences for role in roles
        ):
            warnings.append("SOURCE_FIELD_ROLE_MISMATCH")
        if len(roles) > 1:
            warnings.append("MULTIPLE_POSSIBLE_ROLES")
        primary = roles[0]
        return ResolvedConditionComponent(
            raw_identifier=identifier,
            source_field=source_field,
            identity_status=result.status,
            substance_id=substance.substance_id,
            canonical_name=substance.canonical_name,
            roles=roles,
            primary_role=primary.role_id,
            primary_role_confidence=primary.confidence,
            amount=amount,
            amount_unit=amount_unit,
            source_role_hint=source_role_hint,
            warnings=tuple(sorted(set(warnings))),
            provenance={
                **dict(provenance or {}),
                "identity_match_kind": result.match_kind,
                "identifier_type": resolved_identifier_type,
                "transformation_class": transformation_class,
                "named_family": named_family,
                "definition_id": rules["definition_id"],
            },
        )
    fallback = str(
        rules.get("source_role_hint_fallbacks", {}).get(
            source_role_hint or "",
            rules.get("source_fallback_roles", {}).get(
                source_field, "other_reagent"
            ),
        )
    )
    confidence = 0.55 if source_role_hint else (
        0.7 if source_field == "solvent_cas" else 0.25
    )
    role = ContextualRoleAssignment(
        fallback,
        None,
        confidence,
        (
            "source_role_hint_fallback" if source_role_hint else "source_field_fallback",
            f"identity_{result.status}",
        ),
    )
    return ResolvedConditionComponent(
        raw_identifier=identifier,
        source_field=source_field,
        identity_status=result.status,
        substance_id=None,
        canonical_name=None,
        roles=(role,),
        primary_role=role.role_id,
        primary_role_confidence=role.confidence,
        amount=amount,
        amount_unit=amount_unit,
        source_role_hint=source_role_hint,
        warnings=("CONDITION_IDENTITY_UNCERTAINTY",),
        provenance={
            **dict(provenance or {}),
            "identifier_type": resolved_identifier_type,
            "transformation_class": transformation_class,
            "named_family": named_family,
            "definition_id": rules["definition_id"],
        },
    )


__all__ = ["load_role_resolution_rules", "resolve_contextual_component"]
