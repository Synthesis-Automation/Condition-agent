"""Validated atomic curation of unified condition-substance records."""

from __future__ import annotations

import json
import os
import tempfile
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from .loader import SUBSTANCES_PATH, iter_substance_records, load_role_definitions
from .models import (
    CONDITION_IDENTIFIER_TYPES,
    CONDITION_NAME_IDENTIFIER_TYPES,
    RoleCapability,
    Substance,
    SubstanceIdentifier,
)
from .normalization import (
    normalize_cas,
    normalize_identifier,
)
from .resolver import ConditionRegistry

_CURATION_LOCK = threading.RLock()


class CompoundAdditionError(ValueError):
    """Raised when a proposed registry mutation fails validation."""

    def __init__(self, errors: Iterable[str]) -> None:
        self.errors = tuple(errors)
        super().__init__("; ".join(self.errors))


@dataclass(frozen=True)
class CompoundAliasInput:
    """One supplemental identifier supplied during compound curation."""

    identifier_type: str
    value: str
    language: Optional[str] = None
    is_preferred: bool = False
    allow_ambiguous: bool = False


@dataclass(frozen=True)
class CompoundAdditionRequest:
    """Curator-supplied fields for one condition substance."""

    canonical_name: str
    cas: str
    smiles: Optional[str] = None
    abbreviation: Optional[str] = None
    formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    substance_kind: Optional[str] = None
    density: Optional[float] = None
    boiling_point_c: Optional[float] = None
    melting_point_c: Optional[float] = None
    roles: Tuple[RoleCapability, ...] = ()
    aliases: Tuple[CompoundAliasInput, ...] = ()
    curator_notes: Optional[str] = None


@dataclass(frozen=True)
class CompoundAdditionResult:
    substance: Substance
    canonical_smiles: Optional[str]
    formula: Optional[str]
    molecular_weight: Optional[float]
    alias_count: int
    substances_path: str


@dataclass(frozen=True)
class SubstanceAliasAdditionRequest:
    substance_id: str
    identifier_type: str
    value: str
    language: Optional[str] = None
    is_preferred: bool = False
    allow_ambiguous: bool = False


@dataclass(frozen=True)
class SubstanceAliasAdditionResult:
    added: Tuple[SubstanceIdentifier, ...]
    skipped_existing: Tuple[SubstanceIdentifier, ...]
    substances_path: str


@dataclass(frozen=True)
class _PreparedRecord:
    request: CompoundAdditionRequest
    record: Dict[str, Any]
    canonical_smiles: Optional[str]
    formula: Optional[str]
    molecular_weight: Optional[float]


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def _read_records(path: Path) -> list[Dict[str, Any]]:
    return [payload for _, payload in iter_substance_records(path)]


def _write_jsonl_atomic(path: Path, records: Iterable[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = ""
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="\n",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            for record in records:
                handle.write(_canonical_json(record) + "\n")
        os.replace(temporary_name, path)
    finally:
        if temporary_name and Path(temporary_name).exists():
            Path(temporary_name).unlink()


def _restore_bytes(path: Path, content: bytes) -> None:
    temporary_name = ""
    try:
        with tempfile.NamedTemporaryFile(
            "wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".restore",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(content)
        os.replace(temporary_name, path)
    finally:
        if temporary_name and Path(temporary_name).exists():
            Path(temporary_name).unlink()


def _clear_default_registry_caches(substances_path: Path) -> None:
    if substances_path.resolve() != SUBSTANCES_PATH.resolve():
        return
    from .api import get_registry
    from .vocabulary import condition_registry_definition_versions

    get_registry.cache_clear()
    condition_registry_definition_versions.cache_clear()


def _optional_number(value: Optional[float], field_name: str, errors: list[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        errors.append(f"INVALID_NUMBER:{field_name}")
        return None
    if field_name in {"molecular_weight", "density"} and parsed <= 0:
        errors.append(f"NON_POSITIVE_NUMBER:{field_name}")
    return parsed


def _alias_payload(
    *,
    identifier_type: str,
    value: str,
    language: Optional[str] = None,
    shared: bool = False,
) -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "type": identifier_type,
        "value": value,
    }
    if language:
        payload["language"] = language
    if shared:
        payload["shared"] = True
    return payload


def _prepare_record(
    request: CompoundAdditionRequest,
    *,
    registry: ConditionRegistry,
    editing_substance_id: Optional[str] = None,
) -> _PreparedRecord:
    errors: list[str] = []
    canonical_name = str(request.canonical_name or "").strip()
    cas = normalize_cas(str(request.cas or "").strip())
    if not canonical_name:
        errors.append("MISSING_NAME")
    if cas is None:
        errors.append("INVALID_CAS")
    substance_id = editing_substance_id or (f"cas:{cas}" if cas else "")
    if cas and editing_substance_id is None:
        existing = registry.resolve(cas=cas)
        if existing.status != "unresolved":
            errors.append(f"CAS_ALREADY_REGISTERED:{cas}")
    if editing_substance_id is not None:
        existing = registry.resolve_id(editing_substance_id)
        if existing.substance is None:
            errors.append(f"EDIT_TARGET_NOT_FOUND:{editing_substance_id}")
        elif existing.substance.cas != cas:
            errors.append("CAS_CHANGE_NOT_ALLOWED")

    canonical_smiles: Optional[str] = None
    derived_formula: Optional[str] = None
    derived_mw: Optional[float] = None
    if request.smiles:
        molecule = Chem.MolFromSmiles(str(request.smiles))
        if molecule is None:
            errors.append("INVALID_SMILES")
        else:
            canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
            derived_formula = rdMolDescriptors.CalcMolFormula(molecule)
            derived_mw = round(float(Descriptors.MolWt(molecule)), 6)
    supplied_formula = str(request.formula or "").strip() or None
    if supplied_formula and derived_formula and supplied_formula != derived_formula:
        errors.append(f"FORMULA_SMILES_MISMATCH:{supplied_formula}:{derived_formula}")
    molecular_weight = _optional_number(
        request.molecular_weight, "molecular_weight", errors
    )
    if molecular_weight is not None and derived_mw is not None:
        if abs(molecular_weight - derived_mw) > 0.1:
            errors.append(
                f"MOLECULAR_WEIGHT_SMILES_MISMATCH:{molecular_weight}:{derived_mw}"
            )
    density = _optional_number(request.density, "density", errors)
    boiling_point = _optional_number(request.boiling_point_c, "boiling_point_c", errors)
    melting_point = _optional_number(request.melting_point_c, "melting_point_c", errors)

    known_roles = {
        str(item["id"]) for item in load_role_definitions().get("roles", ())
    }
    seen_roles: set[str] = set()
    for role in request.roles:
        if role.role_id not in known_roles:
            errors.append(f"UNKNOWN_ROLE:{role.role_id}")
        if role.role_id in seen_roles:
            errors.append(f"DUPLICATE_ROLE:{role.role_id}")
        seen_roles.add(role.role_id)

    abbreviation = str(request.abbreviation or "").strip()
    proposed_identifiers: list[Tuple[str, str]] = [("canonical_name", canonical_name)]
    if cas:
        proposed_identifiers.append(("cas", cas))
    if abbreviation:
        proposed_identifiers.append(("abbreviation", abbreviation))
    seen_normalized: set[Tuple[str, str]] = set()
    aliases = []
    for alias in request.aliases:
        alias_type = str(alias.identifier_type or "").strip()
        value = str(alias.value or "").strip()
        if alias_type not in CONDITION_IDENTIFIER_TYPES:
            errors.append(f"UNKNOWN_IDENTIFIER_TYPE:{alias_type}")
            continue
        if alias_type in {"canonical_name", "cas", "abbreviation"}:
            errors.append(f"RESERVED_ALIAS_TYPE:{alias_type}")
            continue
        if not value:
            errors.append("EMPTY_ALIAS_VALUE")
            continue
        aliases.append(alias)
        proposed_identifiers.append((alias_type, value))
    for identifier_type, value in proposed_identifiers:
        normalized = normalize_identifier(value, identifier_type)
        if normalized is None:
            errors.append(f"INVALID_IDENTIFIER_VALUE:{identifier_type}:{value}")
            continue
        namespace = (
            "name" if identifier_type in CONDITION_NAME_IDENTIFIER_TYPES else identifier_type
        )
        key = (namespace, normalized)
        if key in seen_normalized:
            errors.append(f"DUPLICATE_PROPOSED_IDENTIFIER:{value}")
        seen_normalized.add(key)

    for alias in aliases:
        result = registry.resolve_identifier(alias.value, identifier_type=alias.identifier_type)
        owned = (
            editing_substance_id is not None
            and result.substance is not None
            and result.substance.substance_id == editing_substance_id
        )
        if result.status != "unresolved" and not owned and not alias.allow_ambiguous:
            errors.append(f"ALIAS_ALREADY_REGISTERED:{alias.identifier_type}:{alias.value}")
    if errors:
        raise CompoundAdditionError(errors)

    alias_payloads: list[Dict[str, Any]] = []
    if abbreviation:
        alias_payloads.append(
            _alias_payload(
                identifier_type="abbreviation",
                value=abbreviation,
            )
        )
    for alias in aliases:
        alias_payloads.append(
            _alias_payload(
                identifier_type=alias.identifier_type,
                value=alias.value,
                language=alias.language,
                shared=alias.allow_ambiguous,
            )
        )
    record: Dict[str, Any] = {
        "id": substance_id,
        "name": canonical_name,
        "cas": cas,
    }
    if canonical_smiles:
        record["smiles"] = canonical_smiles
    if alias_payloads:
        record["aliases"] = alias_payloads
    if request.roles:
        record["roles"] = [role.role_id for role in request.roles]
    return _PreparedRecord(
        request=request,
        record=record,
        canonical_smiles=canonical_smiles,
        formula=supplied_formula or derived_formula,
        molecular_weight=molecular_weight if molecular_weight is not None else derived_mw,
    )


def _mark_shared_aliases(
    records: list[Dict[str, Any]],
    aliases: Iterable[CompoundAliasInput],
) -> None:
    shared = {
        ("name" if alias.identifier_type in CONDITION_NAME_IDENTIFIER_TYPES else alias.identifier_type,
         normalize_identifier(alias.value, alias.identifier_type))
        for alias in aliases
        if alias.allow_ambiguous
    }
    for record in records:
        for identifier in record.get("aliases") or ():
            identifier_type = str(identifier.get("type") or "")
            namespace = "name" if identifier_type in CONDITION_NAME_IDENTIFIER_TYPES else identifier_type
            key = (namespace, normalize_identifier(str(identifier.get("value") or ""), identifier_type))
            if key in shared:
                identifier["shared"] = True


def _verified_registry(path: Path, substance_id: str) -> Tuple[ConditionRegistry, Substance]:
    registry = ConditionRegistry(substances_path=path)
    result = registry.resolve_id(substance_id)
    if result.status != "resolved" or result.substance is None:
        raise RuntimeError(f"Substance could not be resolved after write: {substance_id}")
    return registry, result.substance


def add_compound(
    request: CompoundAdditionRequest,
    *,
    substances_path: str | Path = SUBSTANCES_PATH,
) -> CompoundAdditionResult:
    """Validate and atomically append one unified substance record."""
    path = Path(substances_path)
    with _CURATION_LOCK:
        registry = ConditionRegistry(substances_path=path)
        prepared = _prepare_record(request, registry=registry)
        records = _read_records(path)
        _mark_shared_aliases(records, request.aliases)
        records.append(prepared.record)
        original = path.read_bytes()
        try:
            _write_jsonl_atomic(path, records)
            _, substance = _verified_registry(path, prepared.record["id"])
        except Exception:
            _restore_bytes(path, original)
            raise
        _clear_default_registry_caches(path)
        return CompoundAdditionResult(
            substance=substance,
            canonical_smiles=prepared.canonical_smiles,
            formula=prepared.formula,
            molecular_weight=prepared.molecular_weight,
            alias_count=len(request.aliases),
            substances_path=str(path),
        )


def update_compound(
    substance_id: str,
    request: CompoundAdditionRequest,
    *,
    substances_path: str | Path = SUBSTANCES_PATH,
) -> CompoundAdditionResult:
    """Validate and atomically replace one unified substance record."""
    path = Path(substances_path)
    with _CURATION_LOCK:
        registry = ConditionRegistry(substances_path=path)
        prepared = _prepare_record(
            request,
            registry=registry,
            editing_substance_id=substance_id,
        )
        records = _read_records(path)
        matches = [
            index
            for index, record in enumerate(records)
            if str(record.get("id") or "") == substance_id
        ]
        if len(matches) != 1:
            raise CompoundAdditionError((f"EDIT_TARGET_DEFINITION_AMBIGUOUS:{substance_id}",))
        _mark_shared_aliases(records, request.aliases)
        records[matches[0]] = prepared.record
        original = path.read_bytes()
        try:
            _write_jsonl_atomic(path, records)
            _, substance = _verified_registry(path, substance_id)
        except Exception:
            _restore_bytes(path, original)
            raise
        _clear_default_registry_caches(path)
        return CompoundAdditionResult(
            substance=substance,
            canonical_smiles=prepared.canonical_smiles,
            formula=prepared.formula,
            molecular_weight=prepared.molecular_weight,
            alias_count=len(request.aliases),
            substances_path=str(path),
        )


def add_substance_aliases(
    requests: Iterable[SubstanceAliasAdditionRequest],
    *,
    substances_path: str | Path = SUBSTANCES_PATH,
) -> SubstanceAliasAdditionResult:
    """Validate and atomically add aliases inside unified substance records."""
    requests = tuple(requests)
    path = Path(substances_path)
    with _CURATION_LOCK:
        registry = ConditionRegistry(substances_path=path)
        records = _read_records(path)
        by_id = {str(record.get("id") or ""): record for record in records}
        added_payloads: list[Tuple[str, Dict[str, Any]]] = []
        skipped: list[SubstanceIdentifier] = []
        errors: list[str] = []
        for request in requests:
            target = registry.resolve_id(request.substance_id)
            if target.status != "resolved" or target.substance is None:
                errors.append(f"UNKNOWN_SUBSTANCE_ID:{request.substance_id}")
                continue
            if request.identifier_type not in CONDITION_IDENTIFIER_TYPES:
                errors.append(f"UNKNOWN_IDENTIFIER_TYPE:{request.identifier_type}")
                continue
            normalized = normalize_identifier(request.value, request.identifier_type)
            if normalized is None:
                errors.append(f"INVALID_ALIAS_VALUE:{request.identifier_type}:{request.value}")
                continue
            existing = registry.resolve_identifier(
                request.value, identifier_type=request.identifier_type
            )
            if (
                existing.status == "resolved"
                and existing.substance is not None
                and existing.substance.substance_id == request.substance_id
            ):
                if existing.matched_identifier is not None:
                    skipped.append(existing.matched_identifier)
                continue
            if existing.status != "unresolved" and not request.allow_ambiguous:
                errors.append(
                    f"ALIAS_ALREADY_REGISTERED:{request.identifier_type}:{request.value}"
                )
                continue
            payload = _alias_payload(
                identifier_type=request.identifier_type,
                value=request.value,
                language=request.language,
                shared=request.allow_ambiguous,
            )
            added_payloads.append((request.substance_id, payload))
        if errors:
            raise CompoundAdditionError(errors)
        shared_aliases = tuple(
            CompoundAliasInput(
                identifier_type=request.identifier_type,
                value=request.value,
                allow_ambiguous=True,
            )
            for request in requests
            if request.allow_ambiguous
        )
        _mark_shared_aliases(records, shared_aliases)
        for substance_id, payload in added_payloads:
            by_id[substance_id].setdefault("aliases", []).append(payload)
        original = path.read_bytes()
        try:
            _write_jsonl_atomic(path, records)
            verified = ConditionRegistry(substances_path=path)
            added = []
            for substance_id, payload in added_payloads:
                result = verified.resolve_identifier(
                    str(payload["value"]),
                    identifier_type=str(payload["type"]),
                )
                if result.status not in {"resolved", "ambiguous"}:
                    raise RuntimeError(f"Alias could not be resolved after write: {payload['value']}")
                if result.matched_identifier is None:
                    raise RuntimeError(
                        f"Alias resolution omitted its identifier: {payload['value']}"
                    )
                added.append(result.matched_identifier)
        except Exception:
            _restore_bytes(path, original)
            raise
        _clear_default_registry_caches(path)
        return SubstanceAliasAdditionResult(
            added=tuple(added),
            skipped_existing=tuple(skipped),
            substances_path=str(path),
        )


__all__ = [
    "CompoundAdditionError",
    "CompoundAdditionRequest",
    "CompoundAdditionResult",
    "CompoundAliasInput",
    "SubstanceAliasAdditionRequest",
    "SubstanceAliasAdditionResult",
    "add_compound",
    "add_substance_aliases",
    "update_compound",
]
