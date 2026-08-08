"""Build canonical resolved condition recipes from contextual components."""

from __future__ import annotations

import hashlib
import json
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple

from .contextual_roles import load_role_resolution_rules, resolve_contextual_component
from .models import (
    ConditionComponentInput,
    ConditionProcessStage,
    ContextualRoleAssignment,
    ResolvedConditionComponent,
    ResolvedConditionRecipe,
)


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def _component_token(component: ResolvedConditionComponent) -> Tuple[Any, ...]:
    identity = component.substance_id or f"raw:{component.raw_identifier.strip().lower()}"
    return (
        identity,
        component.role_status,
        component.primary_role,
        component.amount,
        component.amount_unit,
    )


def _component_core_token(
    component: ResolvedConditionComponent,
) -> Tuple[Any, ...]:
    """Identify a role-aware substance without variant-level quantities."""
    identity = component.substance_id or f"raw:{component.raw_identifier.strip().lower()}"
    return (
        identity,
        component.role_status,
        component.primary_role,
    )


def _sorted_components(
    values: Iterable[ResolvedConditionComponent],
) -> Tuple[ResolvedConditionComponent, ...]:
    return tuple(sorted(values, key=_component_token))


def _merge_duplicate_components(
    values: Iterable[ResolvedConditionComponent],
) -> Tuple[ResolvedConditionComponent, ...]:
    grouped: Dict[str, list[ResolvedConditionComponent]] = {}
    for component in values:
        identity = component.substance_id or (
            f"raw:{component.raw_identifier.strip().lower()}"
        )
        grouped.setdefault(identity, []).append(component)
    merged = []
    for identity in sorted(grouped):
        candidates = grouped[identity]
        primary = sorted(
            candidates,
            key=lambda item: (
                -(item.primary_role_confidence if item.primary_role_confidence is not None else -1.0),
                item.primary_role or "",
                item.source_field,
            ),
        )[0]
        roles: Dict[str, ContextualRoleAssignment] = {}
        for candidate in candidates:
            for role in candidate.roles:
                existing = roles.get(role.role_id)
                if existing is None or role.confidence > existing.confidence:
                    roles[role.role_id] = role
        ranked_roles = tuple(
            sorted(
                roles.values(),
                key=lambda role: (-role.confidence, role.role_id),
            )
        )
        assigned_roles = {
            candidate.primary_role
            for candidate in candidates
            if candidate.role_status == "assigned" and candidate.primary_role
        }
        if len(assigned_roles) == 1:
            selected_role = next(iter(assigned_roles))
            selected_assignment = next(
                role for role in ranked_roles if role.role_id == selected_role
            )
            ranked_roles = (
                selected_assignment,
                *(role for role in ranked_roles if role.role_id != selected_role),
            )
            role_status = "assigned"
            primary_role = selected_role
            primary_confidence = selected_assignment.confidence
        elif len(assigned_roles) > 1:
            role_status = "conflicting"
            primary_role = None
            primary_confidence = None
        elif any(candidate.role_status == "conflicting" for candidate in candidates):
            role_status = "conflicting"
            primary_role = None
            primary_confidence = None
        elif ranked_roles:
            role_status = "ambiguous"
            primary_role = None
            primary_confidence = None
        else:
            role_status = "unassigned"
            primary_role = None
            primary_confidence = None
        source_fields = tuple(sorted({item.source_field for item in candidates}))
        warnings = {warning for item in candidates for warning in item.warnings}
        if len(candidates) > 1:
            warnings.add("DUPLICATE_SOURCE_IDENTITY_MERGED")
        provenance = dict(primary.provenance)
        provenance["source_fields"] = source_fields
        merged.append(
            ResolvedConditionComponent(
                raw_identifier=primary.raw_identifier,
                source_field="|".join(source_fields),
                identity_status=primary.identity_status,
                substance_id=primary.substance_id,
                canonical_name=primary.canonical_name,
                roles=ranked_roles,
                role_status=role_status,
                primary_role=primary_role,
                primary_role_confidence=primary_confidence,
                amount=primary.amount,
                amount_unit=primary.amount_unit,
                source_role_hint=primary.source_role_hint,
                warnings=tuple(sorted(warnings)),
                provenance=provenance,
            )
        )
    return tuple(merged)


def build_resolved_recipe_from_components(
    components: Iterable[ResolvedConditionComponent],
    *,
    temperature_c: Optional[float] = None,
    time_h: Optional[float] = None,
    concentration_m: Optional[float] = None,
    atmosphere: Optional[str] = None,
    stages: Iterable[ConditionProcessStage] = (),
    declared_absences: Iterable[str] = (),
) -> ResolvedConditionRecipe:
    """Build a canonical recipe from already resolved, role-aware components."""
    rules = load_role_resolution_rules()
    bucket_names = tuple(sorted(set(rules["role_buckets"].values())))
    buckets: Dict[str, list[ResolvedConditionComponent]] = {
        name: [] for name in bucket_names
    }
    warnings = []
    for component in _merge_duplicate_components(components):
        bucket = str(
            rules["role_buckets"].get(component.primary_role, "other_components")
        )
        buckets.setdefault(bucket, []).append(component)
        warnings.extend(component.warnings)
    sorted_buckets = {
        name: _sorted_components(values) for name, values in buckets.items()
    }
    definition_versions = {
        "role_resolution.v2.json": str(rules["schema_version"]),
    }
    normalized_stages = tuple(sorted(stages, key=lambda item: item.stage_index))
    normalized_absences = tuple(
        sorted({str(value).strip().lower() for value in declared_absences if str(value).strip()})
    )
    identity_payload = {
        "buckets": {
            name: tuple(_component_token(component) for component in values)
            for name, values in sorted(sorted_buckets.items())
        },
        "temperature_c": temperature_c,
        "time_h": time_h,
        "concentration_m": concentration_m,
        "atmosphere": atmosphere,
        "stages": tuple(
            (
                stage.stage_index,
                stage.temperature_c,
                stage.time_h,
                stage.atmosphere,
            )
            for stage in normalized_stages
        ),
        "declared_absences": normalized_absences,
        "definition_versions": definition_versions,
        "schema_version": "2.0",
    }
    core_identity_payload = {
        "buckets": {
            name: tuple(_component_core_token(component) for component in values)
            for name, values in sorted(sorted_buckets.items())
        },
        "definition_versions": definition_versions,
        "schema_version": "2.0",
    }
    recipe_core_id = "RCORE2:" + hashlib.sha256(
        _canonical_json(core_identity_payload).encode("utf-8")
    ).hexdigest()
    recipe_id = "RCR2:" + hashlib.sha256(
        _canonical_json(identity_payload).encode("utf-8")
    ).hexdigest()
    return ResolvedConditionRecipe(
        recipe_id=recipe_id,
        recipe_core_id=recipe_core_id,
        catalysts=sorted_buckets.get("catalysts", ()),
        ligands=sorted_buckets.get("ligands", ()),
        bases=sorted_buckets.get("bases", ()),
        acids=sorted_buckets.get("acids", ()),
        condensation_agents=sorted_buckets.get("condensation_agents", ()),
        oxidants=sorted_buckets.get("oxidants", ()),
        reductants=sorted_buckets.get("reductants", ()),
        additives=sorted_buckets.get("additives", ()),
        solvents=sorted_buckets.get("solvents", ()),
        other_components=sorted_buckets.get("other_components", ()),
        temperature_c=temperature_c,
        time_h=time_h,
        concentration_m=concentration_m,
        atmosphere=atmosphere,
        stages=normalized_stages,
        declared_absences=normalized_absences,
        warnings=tuple(sorted(set(warnings))),
        definition_versions=definition_versions,
    )


def build_resolved_recipe_from_inputs(
    inputs: Iterable[ConditionComponentInput],
    *,
    preferred_roles: Iterable[str] = (),
    temperature_c: Optional[float] = None,
    time_h: Optional[float] = None,
    concentration_m: Optional[float] = None,
    atmosphere: Optional[str] = None,
    stages: Iterable[ConditionProcessStage] = (),
    declared_absences: Iterable[str] = (),
) -> ResolvedConditionRecipe:
    """Resolve typed raw inputs and build one canonical condition recipe."""
    components = (
        resolve_contextual_component(
            item.raw_identifier,
            source_field=item.source_field,
            identifier_type=item.identifier_type,
            source_role_hint=item.source_role_hint,
            preferred_roles=preferred_roles,
            amount=item.amount,
            amount_unit=item.amount_unit,
            provenance=item.provenance,
        )
        for item in inputs
        if item.raw_identifier.strip()
    )
    return build_resolved_recipe_from_components(
        components,
        temperature_c=temperature_c,
        time_h=time_h,
        concentration_m=concentration_m,
        atmosphere=atmosphere,
        stages=stages,
        declared_absences=declared_absences,
    )


def build_resolved_recipe(
    source_components: Mapping[str, Iterable[str]],
    *,
    preferred_roles: Iterable[str] = (),
    temperature_c: Optional[float] = None,
    time_h: Optional[float] = None,
    concentration_m: Optional[float] = None,
    atmosphere: Optional[str] = None,
) -> ResolvedConditionRecipe:
    """Resolve identities and contextual roles into stable recipe buckets."""
    resolved_components = []
    seen = set()
    for source_field in ("catalyst_cas", "reagent_cas", "solvent_cas"):
        for identifier in source_components.get(source_field, ()):
            raw = str(identifier).strip()
            source_identity = (source_field, raw)
            if not raw or source_identity in seen:
                continue
            seen.add(source_identity)
            component = resolve_contextual_component(
                raw,
                source_field=source_field,
                identifier_type="cas",
                preferred_roles=preferred_roles,
            )
            resolved_components.append(component)
    return build_resolved_recipe_from_components(
        resolved_components,
        temperature_c=temperature_c,
        time_h=time_h,
        concentration_m=concentration_m,
        atmosphere=atmosphere,
    )


__all__ = [
    "build_resolved_recipe",
    "build_resolved_recipe_from_components",
    "build_resolved_recipe_from_inputs",
]
