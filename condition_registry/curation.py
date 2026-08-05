"""Validated, atomic curation workflow for adding condition substances."""

from __future__ import annotations

import csv
import hashlib
import os
import tempfile
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from .loader import (
    ADDITIONS_PATH,
    IDENTIFIERS_PATH,
    SUBSTANCES_PATH,
    load_taxonomy,
)
from .models import (
    CONDITION_IDENTIFIER_TYPES,
    CONDITION_NAME_IDENTIFIER_TYPES,
    RoleAssignment,
    Substance,
)
from .normalization import (
    identifier_normalization_profile,
    normalize_cas,
    normalize_identifier,
)
from .resolver import ConditionRegistry

ADDITION_FIELDNAMES = (
    "substance_id",
    "name",
    "abbreviation",
    "cas",
    "smiles",
    "formula",
    "type",
    "density",
    "mw",
    "bp",
    "mp",
    "volatile",
    "viscose",
    "role_1",
    "family_1",
    "tag_1",
    "role_2",
    "family_2",
    "tag_2",
    "status",
    "source",
    "curator_notes",
)

IDENTIFIER_FIELDNAMES = (
    "identifier_id",
    "substance_id",
    "identifier_type",
    "value",
    "language",
    "is_preferred",
    "source",
    "confidence",
    "status",
    "normalization_profile",
    "allow_ambiguous",
)

_CURATION_LOCK = threading.Lock()


class CompoundAdditionError(ValueError):
    """Raised when a proposed compound cannot be curated safely."""

    def __init__(self, errors: Iterable[str]) -> None:
        self.errors = tuple(dict.fromkeys(str(error) for error in errors))
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
    """Curator-supplied fields for one new condition substance."""

    canonical_name: str
    cas: str
    source: str
    smiles: Optional[str] = None
    abbreviation: Optional[str] = None
    formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    substance_kind: Optional[str] = None
    density: Optional[float] = None
    boiling_point_c: Optional[float] = None
    melting_point_c: Optional[float] = None
    roles: Tuple[RoleAssignment, ...] = ()
    aliases: Tuple[CompoundAliasInput, ...] = ()
    curator_notes: Optional[str] = None


@dataclass(frozen=True)
class CompoundAdditionResult:
    """Identity and normalized values written by a curation operation."""

    substance: Substance
    canonical_smiles: Optional[str]
    formula: Optional[str]
    molecular_weight: Optional[float]
    alias_count: int
    additions_path: str
    identifiers_path: str


@dataclass(frozen=True)
class _PreparedAddition:
    request: CompoundAdditionRequest
    substance_id: str
    cas: str
    canonical_name: str
    canonical_smiles: Optional[str]
    formula: Optional[str]
    molecular_weight: Optional[float]
    aliases: Tuple[CompoundAliasInput, ...]


def _optional_float(value: Optional[float], field_name: str, errors: list[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        errors.append(f"INVALID_NUMBER:{field_name}")
        return None
    if parsed <= 0 and field_name in {"molecular_weight", "density"}:
        errors.append(f"NON_POSITIVE_NUMBER:{field_name}")
    return parsed


def _prepare_request(
    request: CompoundAdditionRequest,
    *,
    registry: ConditionRegistry,
) -> _PreparedAddition:
    errors: list[str] = []
    canonical_name = str(request.canonical_name or "").strip()
    source = str(request.source or "").strip()
    if not canonical_name:
        errors.append("MISSING_CANONICAL_NAME")
    if not source:
        errors.append("MISSING_SOURCE")
    cas = normalize_cas(request.cas)
    if cas is None:
        errors.append("INVALID_CAS")
        cas = str(request.cas or "").strip()
    elif registry.resolve(cas=cas).status != "unresolved":
        errors.append(f"CAS_ALREADY_REGISTERED:{cas}")
    if canonical_name:
        name_result = registry.resolve(name=canonical_name)
        if name_result.status != "unresolved":
            errors.append(f"NAME_ALREADY_REGISTERED:{canonical_name}")

    canonical_smiles = None
    derived_formula = None
    derived_mw = None
    smiles = str(request.smiles or "").strip()
    if smiles:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            errors.append("INVALID_SMILES")
        else:
            canonical_smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
            derived_formula = rdMolDescriptors.CalcMolFormula(molecule)
            derived_mw = round(float(Descriptors.MolWt(molecule)), 6)

    supplied_formula = str(request.formula or "").strip() or None
    if supplied_formula and derived_formula and supplied_formula != derived_formula:
        errors.append(
            f"FORMULA_SMILES_MISMATCH:{supplied_formula}:{derived_formula}"
        )
    molecular_weight = _optional_float(
        request.molecular_weight, "molecular_weight", errors
    )
    if molecular_weight is not None and derived_mw is not None:
        if abs(molecular_weight - derived_mw) > 0.1:
            errors.append(
                f"MOLECULAR_WEIGHT_SMILES_MISMATCH:{molecular_weight}:{derived_mw}"
            )
    _optional_float(request.density, "density", errors)
    _optional_float(request.boiling_point_c, "boiling_point_c", errors)
    _optional_float(request.melting_point_c, "melting_point_c", errors)

    taxonomy = load_taxonomy()
    known_roles = {str(item["id"]) for item in taxonomy["roles"]}
    family_roles = {
        str(item["id"]): str(item["role_id"])
        for item in taxonomy["families"]
    }
    if len(request.roles) > 2:
        errors.append("TOO_MANY_ROLES_FOR_SUBSTANCE_SCHEMA")
    for assignment in request.roles:
        if assignment.role_id not in known_roles:
            errors.append(f"UNKNOWN_ROLE:{assignment.role_id}")
        if assignment.family_id and assignment.family_id not in family_roles:
            errors.append(f"UNKNOWN_FAMILY:{assignment.family_id}")
        if (
            assignment.family_id
            and family_roles.get(assignment.family_id) != assignment.role_id
        ):
            errors.append(
                "ROLE_FAMILY_MISMATCH:"
                f"{assignment.role_id}:{assignment.family_id}"
            )

    aliases = []
    seen_aliases: set[Tuple[str, str]] = set()
    normalized_canonical_name = normalize_identifier(
        canonical_name, "canonical_name"
    )
    if normalized_canonical_name:
        seen_aliases.add(("name", normalized_canonical_name))
    abbreviation = str(request.abbreviation or "").strip()
    if abbreviation:
        normalized_abbreviation = normalize_identifier(abbreviation, "abbreviation")
        if normalized_abbreviation:
            seen_aliases.add(("name", normalized_abbreviation))
        existing = registry.resolve(name=abbreviation)
        if existing.status != "unresolved":
            errors.append(f"ABBREVIATION_ALREADY_REGISTERED:{abbreviation}")
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
        normalized = normalize_identifier(value, alias_type)
        if normalized is None:
            errors.append(f"INVALID_ALIAS_VALUE:{alias_type}:{value}")
            continue
        namespace = "name" if alias_type in CONDITION_NAME_IDENTIFIER_TYPES else alias_type
        key = (namespace, normalized)
        if key in seen_aliases:
            errors.append(f"DUPLICATE_PROPOSED_ALIAS:{value}")
            continue
        seen_aliases.add(key)
        existing = (
            registry.resolve(name=value)
            if alias_type in CONDITION_NAME_IDENTIFIER_TYPES
            else registry.resolve_identifier(value, identifier_type=alias_type)
        )
        if existing.status != "unresolved" and not alias.allow_ambiguous:
            errors.append(f"ALIAS_ALREADY_REGISTERED:{alias_type}:{value}")
        aliases.append(
            CompoundAliasInput(
                identifier_type=alias_type,
                value=value,
                language=str(alias.language or "").strip() or None,
                is_preferred=bool(alias.is_preferred),
                allow_ambiguous=bool(alias.allow_ambiguous),
            )
        )
    if errors:
        raise CompoundAdditionError(errors)
    return _PreparedAddition(
        request=request,
        substance_id=f"cas:{cas}",
        cas=cas,
        canonical_name=canonical_name,
        canonical_smiles=canonical_smiles,
        formula=supplied_formula or derived_formula,
        molecular_weight=molecular_weight if molecular_weight is not None else derived_mw,
        aliases=tuple(aliases),
    )


def _format_number(value: Optional[float]) -> str:
    return "" if value is None else format(float(value), ".12g")


def _addition_row(prepared: _PreparedAddition) -> Dict[str, str]:
    roles = (*prepared.request.roles, None, None)[:2]
    row = {field: "" for field in ADDITION_FIELDNAMES}
    row.update(
        {
            "substance_id": prepared.substance_id,
            "name": prepared.canonical_name,
            "abbreviation": str(prepared.request.abbreviation or "").strip(),
            "cas": prepared.cas,
            "smiles": prepared.canonical_smiles or "",
            "formula": prepared.formula or "",
            "type": str(prepared.request.substance_kind or "").strip(),
            "density": _format_number(prepared.request.density),
            "mw": _format_number(prepared.molecular_weight),
            "bp": _format_number(prepared.request.boiling_point_c),
            "mp": _format_number(prepared.request.melting_point_c),
            "status": "active",
            "source": str(prepared.request.source).strip(),
            "curator_notes": str(prepared.request.curator_notes or "").strip(),
        }
    )
    for index, assignment in enumerate(roles, start=1):
        if assignment is None:
            continue
        row[f"role_{index}"] = assignment.role_id
        row[f"family_{index}"] = assignment.family_id or ""
        row[f"tag_{index}"] = assignment.tag or ""
    return row


def _identifier_rows(prepared: _PreparedAddition) -> Tuple[Dict[str, str], ...]:
    rows = []
    source = str(prepared.request.source).strip()
    for alias in prepared.aliases:
        normalized = normalize_identifier(alias.value, alias.identifier_type) or alias.value
        digest = hashlib.sha256(
            f"{prepared.substance_id}|{alias.identifier_type}|{normalized}".encode("utf-8")
        ).hexdigest()[:16]
        rows.append(
            {
                "identifier_id": f"sid:{prepared.cas}:{digest}",
                "substance_id": prepared.substance_id,
                "identifier_type": alias.identifier_type,
                "value": alias.value,
                "language": alias.language or "",
                "is_preferred": str(alias.is_preferred).lower(),
                "source": source,
                "confidence": "1.0",
                "status": "active",
                "normalization_profile": (
                    identifier_normalization_profile(alias.identifier_type) or ""
                ),
                "allow_ambiguous": str(alias.allow_ambiguous).lower(),
            }
        )
    return tuple(rows)


def _declare_existing_shared_aliases(
    rows: Iterable[Mapping[str, str]],
    aliases: Iterable[CompoundAliasInput],
) -> Tuple[Dict[str, str], ...]:
    """Mark matching supplemental identifiers shared when explicitly requested."""
    shared_keys = set()
    for alias in aliases:
        if not alias.allow_ambiguous:
            continue
        normalized = normalize_identifier(alias.value, alias.identifier_type)
        if normalized is None:
            continue
        namespace = (
            "name"
            if alias.identifier_type in CONDITION_NAME_IDENTIFIER_TYPES
            else alias.identifier_type
        )
        shared_keys.add((namespace, normalized))
    updated = []
    for raw_row in rows:
        row = dict(raw_row)
        identifier_type = str(row.get("identifier_type") or "")
        value = str(row.get("value") or "")
        normalized = normalize_identifier(value, identifier_type)
        namespace = (
            "name"
            if identifier_type in CONDITION_NAME_IDENTIFIER_TYPES
            else identifier_type
        )
        if normalized is not None and (namespace, normalized) in shared_keys:
            row["allow_ambiguous"] = "true"
        updated.append(row)
    return tuple(updated)


def _read_csv(path: Path, expected_fields: Sequence[str]) -> list[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(expected_fields):
            raise CompoundAdditionError((f"DEFINITION_HEADER_MISMATCH:{path.name}",))
        return [dict(row) for row in reader]


def _write_csv_atomic(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = ""
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8-sig",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
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


def add_compound(
    request: CompoundAdditionRequest,
    *,
    substances_path: str | Path = SUBSTANCES_PATH,
    additions_path: str | Path = ADDITIONS_PATH,
    identifiers_path: str | Path = IDENTIFIERS_PATH,
) -> CompoundAdditionResult:
    """Validate and atomically add one compound plus supplemental aliases."""
    substances_path = Path(substances_path)
    additions_path = Path(additions_path)
    identifiers_path = Path(identifiers_path)
    with _CURATION_LOCK:
        registry = ConditionRegistry(
            substances_path=substances_path,
            additions_path=additions_path,
            identifiers_path=identifiers_path,
        )
        prepared = _prepare_request(request, registry=registry)
        addition_rows = _read_csv(additions_path, ADDITION_FIELDNAMES)
        identifier_rows = _declare_existing_shared_aliases(
            _read_csv(identifiers_path, IDENTIFIER_FIELDNAMES),
            prepared.aliases,
        )
        original_additions = additions_path.read_bytes()
        original_identifiers = identifiers_path.read_bytes()
        try:
            _write_csv_atomic(
                additions_path,
                ADDITION_FIELDNAMES,
                (*addition_rows, _addition_row(prepared)),
            )
            _write_csv_atomic(
                identifiers_path,
                IDENTIFIER_FIELDNAMES,
                (*identifier_rows, *_identifier_rows(prepared)),
            )
            verified_registry = ConditionRegistry(
                substances_path=substances_path,
                additions_path=additions_path,
                identifiers_path=identifiers_path,
            )
            verified = verified_registry.resolve_id(prepared.substance_id)
            if verified.status != "resolved" or verified.substance is None:
                raise RuntimeError("New compound could not be resolved after write")
            for alias in prepared.aliases:
                alias_result = verified_registry.resolve_identifier(
                    alias.value,
                    identifier_type=alias.identifier_type,
                )
                if alias_result.status not in {"resolved", "ambiguous"}:
                    raise RuntimeError(
                        f"New alias could not be resolved after write: {alias.value}"
                    )
        except Exception:
            _restore_bytes(additions_path, original_additions)
            _restore_bytes(identifiers_path, original_identifiers)
            raise
        if (
            substances_path.resolve() == SUBSTANCES_PATH.resolve()
            and additions_path.resolve() == ADDITIONS_PATH.resolve()
            and identifiers_path.resolve() == IDENTIFIERS_PATH.resolve()
        ):
            from .api import get_registry
            from .vocabulary import condition_registry_definition_versions

            get_registry.cache_clear()
            condition_registry_definition_versions.cache_clear()
        return CompoundAdditionResult(
            substance=verified.substance,
            canonical_smiles=prepared.canonical_smiles,
            formula=prepared.formula,
            molecular_weight=prepared.molecular_weight,
            alias_count=len(prepared.aliases),
            additions_path=str(additions_path),
            identifiers_path=str(identifiers_path),
        )


__all__ = [
    "ADDITION_FIELDNAMES",
    "IDENTIFIER_FIELDNAMES",
    "CompoundAdditionError",
    "CompoundAdditionRequest",
    "CompoundAdditionResult",
    "CompoundAliasInput",
    "add_compound",
]
