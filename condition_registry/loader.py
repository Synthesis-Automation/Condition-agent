"""Load the unified condition-substance registry and role vocabulary."""

from __future__ import annotations

import json
import hashlib
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping, Tuple

from .models import RoleCapability, Substance, SubstanceIdentifier
from .normalization import identifier_normalization_profile, normalize_cas

DEFINITIONS_DIR = Path(__file__).with_name("definitions")
SUBSTANCES_PATH = DEFINITIONS_DIR / "substances.v2.jsonl"
ROLES_PATH = DEFINITIONS_DIR / "roles.v2.json"

_SUBSTANCE_FIELDS = {"id", "name", "cas", "smiles", "aliases", "roles"}
_ALIAS_FIELDS = {"type", "value", "language", "shared"}


@lru_cache(maxsize=1)
def load_role_definitions() -> Dict[str, Any]:
    """Load and validate the package-owned role vocabulary."""
    with ROLES_PATH.open("r", encoding="utf-8") as handle:
        definitions = dict(json.load(handle))
    if str(definitions.get("schema_version") or "") != "roles.v2":
        raise ValueError("Unsupported condition role schema")
    roles = tuple(definitions.get("roles") or ())
    role_ids = tuple(str(item.get("id") or "") for item in roles)
    if not role_ids or any(not role_id for role_id in role_ids):
        raise ValueError("Condition roles require non-empty IDs")
    if len(role_ids) != len(set(role_ids)):
        raise ValueError("Condition role IDs must be unique")
    return definitions


def iter_substance_records(
    path: str | Path = SUBSTANCES_PATH,
) -> Iterator[Tuple[int, Dict[str, Any]]]:
    """Yield one parsed substance record per non-empty JSONL line."""
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            payload = json.loads(line)
            if not isinstance(payload, dict):
                raise ValueError(f"Substance line {line_number} must be an object")
            yield line_number, payload


def _alias_identifier(
    payload: Mapping[str, Any],
    *,
    substance_id: str,
) -> SubstanceIdentifier:
    identifier_type = str(payload.get("type") or "").strip()
    value = str(payload.get("value") or "").strip()
    digest = hashlib.sha256(
        f"{substance_id}|{identifier_type}|{value}".encode("utf-8")
    ).hexdigest()[:16]
    return SubstanceIdentifier(
        identifier_id=f"alias:{substance_id}:{digest}",
        substance_id=substance_id,
        identifier_type=identifier_type,
        value=value,
        language=str(payload.get("language") or "").strip() or None,
        normalization_profile=identifier_normalization_profile(identifier_type),
        allow_ambiguous=bool(payload.get("shared", False)),
    )


def _role_from_id(
    role_id: str,
    *,
    substance_id: str,
) -> RoleCapability:
    return RoleCapability(
        role_id=role_id,
        capability_id=f"role:{substance_id}:{role_id}",
        evidence="curated_registry",
    )


def record_to_substance(payload: Mapping[str, Any], line_number: int) -> Substance:
    """Convert and locally validate one unified substance record."""
    unknown_fields = set(payload) - _SUBSTANCE_FIELDS
    if unknown_fields:
        raise ValueError(
            f"Unsupported fields on substance line {line_number}: "
            + ", ".join(sorted(unknown_fields))
        )
    substance_id = str(payload.get("id") or "").strip()
    canonical_name = str(payload.get("name") or "").strip()
    if not substance_id or not canonical_name:
        raise ValueError(f"Substance line {line_number} requires ID and name")
    cas = str(payload.get("cas") or "").strip() or None
    if cas is not None and normalize_cas(cas) is None:
        raise ValueError(f"Invalid CAS on substance line {line_number}")
    identifiers = [
        SubstanceIdentifier(
            identifier_id=f"identity:{substance_id}:name",
            substance_id=substance_id,
            identifier_type="canonical_name",
            value=canonical_name,
            is_preferred=True,
            normalization_profile="chemical_name_v1",
        )
    ]
    if cas:
        identifiers.append(
            SubstanceIdentifier(
                identifier_id=f"identity:{substance_id}:cas",
                substance_id=substance_id,
                identifier_type="cas",
                value=cas,
                is_preferred=True,
                normalization_profile="cas_v1",
            )
        )
    alias_payloads = payload.get("aliases") or ()
    if not isinstance(alias_payloads, (list, tuple)):
        raise ValueError(f"Aliases on substance line {line_number} must be a list")
    for alias in alias_payloads:
        if not isinstance(alias, Mapping):
            raise ValueError(f"Aliases on substance line {line_number} must be objects")
        unknown_alias_fields = set(alias) - _ALIAS_FIELDS
        if unknown_alias_fields:
            raise ValueError(
                f"Unsupported alias fields on substance line {line_number}: "
                + ", ".join(sorted(unknown_alias_fields))
            )
        identifiers.append(
            _alias_identifier(
                alias,
                substance_id=substance_id,
            )
        )
    role_ids = payload.get("roles") or ()
    if not isinstance(role_ids, (list, tuple)) or any(
        not isinstance(value, str) for value in role_ids
    ):
        raise ValueError(f"Roles on substance line {line_number} must be a string list")
    return Substance(
        substance_id=substance_id,
        canonical_name=canonical_name,
        cas=cas,
        smiles=str(payload.get("smiles") or "").strip() or None,
        identifiers=tuple(identifiers),
        roles=tuple(
            _role_from_id(
                str(value),
                substance_id=substance_id,
            )
            for value in role_ids
        ),
        properties={},
        provenance={},
    )


def load_substances(
    *,
    substances_path: str | Path = SUBSTANCES_PATH,
) -> Tuple[Substance, ...]:
    """Load active substances from the unified JSONL definition."""
    substances = tuple(
        record_to_substance(payload, line_number)
        for line_number, payload in iter_substance_records(substances_path)
    )
    known_roles = {
        str(item["id"]) for item in load_role_definitions().get("roles", ())
    }
    seen_substances: set[str] = set()
    seen_identifiers: set[str] = set()
    seen_capabilities: set[str] = set()
    for substance in substances:
        if substance.substance_id in seen_substances:
            raise ValueError(f"Duplicate substance ID: {substance.substance_id}")
        seen_substances.add(substance.substance_id)
        for identifier in substance.identifiers:
            if identifier.identifier_id in seen_identifiers:
                raise ValueError(
                    f"Duplicate substance identifier ID: {identifier.identifier_id}"
                )
            seen_identifiers.add(identifier.identifier_id)
        for capability in substance.roles:
            if capability.role_id not in known_roles:
                raise ValueError(
                    f"Unknown role capability: {substance.substance_id}:"
                    f"{capability.role_id}"
                )
            capability_id = capability.capability_id or ""
            if not capability_id or capability_id in seen_capabilities:
                raise ValueError(f"Invalid or duplicate role capability: {capability_id}")
            seen_capabilities.add(capability_id)
    return substances


__all__ = [
    "DEFINITIONS_DIR",
    "ROLES_PATH",
    "SUBSTANCES_PATH",
    "iter_substance_records",
    "load_role_definitions",
    "load_substances",
    "record_to_substance",
]
