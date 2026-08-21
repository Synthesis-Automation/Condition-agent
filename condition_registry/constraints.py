"""Typed normalization of user-authorized condition constraints."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from typing import Any, Dict, Literal, Mapping, Optional, Tuple

from .api import get_registry
from .models import CONDITION_RECIPE_COMPONENT_BUCKETS
from .vocabulary import load_condition_vocabulary


CONDITION_CONSTRAINT_SCHEMA_VERSION = "1.0"
ConditionConstraintProvenance = Literal["explicit_user", "confirmed_user"]
ConditionConstraintKind = Literal[
    "excluded_substance",
    "excluded_role",
    "maximum_temperature_c",
    "required_atmosphere",
    "required_solvent",
]

_ATMOSPHERES = {
    "n2": "nitrogen",
    "nitrogen": "nitrogen",
    "ar": "argon",
    "argon": "argon",
    "air": "air",
    "o2": "oxygen",
    "oxygen": "oxygen",
    "h2": "hydrogen",
    "hydrogen": "hydrogen",
}


def _constraint_id(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:20]
    return f"CCON1:{digest}"


@dataclass(frozen=True)
class ConditionConstraint:
    """One normalized hard constraint owned by the condition registry."""

    constraint_id: str
    kind: ConditionConstraintKind
    provenance: ConditionConstraintProvenance
    normalized_value: str
    numeric_value: Optional[float] = None
    source_value: str = ""
    resolution_evidence: Tuple[str, ...] = ()
    schema_version: str = CONDITION_CONSTRAINT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.constraint_id or not self.normalized_value:
            raise ValueError("normalized condition constraint identity is required")
        if self.provenance not in {"explicit_user", "confirmed_user"}:
            raise ValueError("hard condition constraints require user authority")
        if self.kind == "maximum_temperature_c" and self.numeric_value is None:
            raise ValueError("maximum temperature requires a numeric value")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ConditionConstraintSet:
    """Unique normalized condition constraints for a recommendation request."""

    constraints: Tuple[ConditionConstraint, ...] = ()
    schema_version: str = CONDITION_CONSTRAINT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        identities = [(item.kind, item.normalized_value) for item in self.constraints]
        if len(identities) != len(set(identities)):
            raise ValueError("condition constraints must be unique")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def append(self, constraint: ConditionConstraint) -> "ConditionConstraintSet":
        return ConditionConstraintSet(self.constraints + (constraint,))


@dataclass(frozen=True)
class ConditionConstraintResolution:
    """Resolution outcome that preserves ambiguity instead of guessing."""

    status: Literal["resolved", "ambiguous", "unresolved", "invalid"]
    constraint: Optional[ConditionConstraint] = None
    candidates: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()


def normalize_condition_constraint(
    kind: ConditionConstraintKind,
    value: Any,
    *,
    provenance: ConditionConstraintProvenance,
) -> ConditionConstraintResolution:
    """Resolve one user-authorized raw value using registry-owned semantics."""

    source_value = str(value).strip()
    if not source_value:
        return ConditionConstraintResolution(
            status="invalid",
            warnings=("CONDITION_CONSTRAINT_VALUE_EMPTY",),
        )
    normalized: str
    numeric_value = None
    evidence = []
    if kind in {"excluded_substance", "required_solvent"}:
        resolution = get_registry().resolve_identifier(source_value)
        if resolution.status == "ambiguous":
            return ConditionConstraintResolution(
                status="ambiguous",
                candidates=resolution.candidates,
                warnings=resolution.warnings,
            )
        if resolution.status != "resolved" or resolution.substance is None:
            resolution = get_registry().resolve_id(source_value)
        if resolution.status != "resolved" or resolution.substance is None:
            return ConditionConstraintResolution(
                status="unresolved",
                warnings=tuple(resolution.warnings) or (
                    "CONDITION_SUBSTANCE_UNRESOLVED",
                ),
            )
        normalized = resolution.substance.substance_id
        evidence.extend(
            (
                f"registry_substance:{normalized}",
                f"identity_match:{resolution.match_kind or 'exact_substance_id'}",
            )
        )
    elif kind == "excluded_role":
        normalized = source_value.casefold().replace(" ", "_")
        role_ids = load_condition_vocabulary().role_ids
        if normalized.endswith("s") and normalized[:-1] in role_ids:
            normalized = normalized[:-1]
        if normalized not in role_ids:
            return ConditionConstraintResolution(
                status="invalid",
                warnings=("CONDITION_ROLE_UNSUPPORTED",),
            )
        evidence.append(f"condition_role:{normalized}")
    elif kind == "maximum_temperature_c":
        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            return ConditionConstraintResolution(
                status="invalid",
                warnings=("CONDITION_TEMPERATURE_INVALID",),
            )
        if numeric_value < -273.15 or numeric_value > 1000.0:
            return ConditionConstraintResolution(
                status="invalid",
                warnings=("CONDITION_TEMPERATURE_OUT_OF_RANGE",),
            )
        normalized = format(numeric_value, ".12g")
        evidence.append("temperature_unit:celsius")
    elif kind == "required_atmosphere":
        token = source_value.casefold().replace(" ", "")
        normalized = _ATMOSPHERES.get(token, "")
        if not normalized:
            return ConditionConstraintResolution(
                status="unresolved",
                warnings=("CONDITION_ATMOSPHERE_UNRESOLVED",),
            )
        evidence.append(f"atmosphere:{normalized}")
    else:
        raise ValueError(f"unsupported condition constraint kind: {kind}")
    identity = {
        "kind": kind,
        "normalized_value": normalized,
        "numeric_value": numeric_value,
        "provenance": provenance,
        "schema_version": CONDITION_CONSTRAINT_SCHEMA_VERSION,
    }
    return ConditionConstraintResolution(
        status="resolved",
        constraint=ConditionConstraint(
            constraint_id=_constraint_id(identity),
            kind=kind,
            provenance=provenance,
            normalized_value=normalized,
            numeric_value=numeric_value,
            source_value=source_value,
            resolution_evidence=tuple(evidence),
        ),
    )


def condition_constraint_conflicts(
    recipe: Mapping[str, Any],
    constraints: ConditionConstraintSet,
) -> Tuple[str, ...]:
    """Return deterministic conflict codes for a canonical resolved recipe."""

    components = []
    for bucket in CONDITION_RECIPE_COMPONENT_BUCKETS:
        values = recipe.get(bucket) or ()
        if isinstance(values, (list, tuple)):
            components.extend(
                (bucket, item) for item in values if isinstance(item, Mapping)
            )
    flat = recipe.get("components") or ()
    if not components and isinstance(flat, (list, tuple)):
        components.extend(
            (str(item.get("role") or ""), item)
            for item in flat
            if isinstance(item, Mapping)
        )
    conflicts = []
    for constraint in constraints.constraints:
        code = f"condition_constraint:{constraint.constraint_id}"
        if constraint.kind == "excluded_substance" and any(
            str(item.get("substance_id") or "") == constraint.normalized_value
            for _, item in components
        ):
            conflicts.append(f"{code}:excluded_substance")
        elif constraint.kind == "excluded_role" and any(
            (
                bucket.rstrip("s") == constraint.normalized_value
                or str(item.get("primary_role") or item.get("role") or "")
                == constraint.normalized_value
            )
            for bucket, item in components
        ):
            conflicts.append(f"{code}:excluded_role")
        elif constraint.kind == "required_solvent":
            solvents = [
                item
                for bucket, item in components
                if bucket == "solvents"
                or str(item.get("primary_role") or item.get("role") or "")
                == "solvent"
            ]
            if not any(
                str(item.get("substance_id") or "") == constraint.normalized_value
                for item in solvents
            ):
                conflicts.append(f"{code}:required_solvent_missing")
        elif constraint.kind == "maximum_temperature_c":
            temperature = recipe.get("temperature_c")
            if temperature is None:
                protocol = recipe.get("synthesis_protocol") or {}
                if isinstance(protocol, Mapping):
                    temperature = protocol.get("temperature_c")
            if temperature is None:
                conflicts.append(f"{code}:temperature_unknown")
            elif float(temperature) > float(constraint.numeric_value):
                conflicts.append(f"{code}:temperature_exceeded")
        elif constraint.kind == "required_atmosphere":
            atmosphere = str(recipe.get("atmosphere") or "").casefold()
            aliases = {
                "nitrogen": ("nitrogen", "n2"),
                "argon": ("argon", "ar"),
                "oxygen": ("oxygen", "o2"),
                "hydrogen": ("hydrogen", "h2"),
                "air": ("air",),
            }[constraint.normalized_value]
            if not any(token in atmosphere for token in aliases):
                conflicts.append(f"{code}:required_atmosphere_missing")
    return tuple(conflicts)
