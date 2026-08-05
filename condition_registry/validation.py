"""Validation and audit of package-owned condition registry data."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, Iterable, List, Mapping, Tuple

from .loader import (
    ADDITIONS_PATH,
    IDENTIFIERS_PATH,
    SUBSTANCES_PATH,
    iter_addition_rows,
    iter_identifier_rows,
    iter_substance_rows,
    load_taxonomy,
    row_to_identifier,
)
from .models import CONDITION_NAME_IDENTIFIER_TYPES
from .normalization import (
    identifier_normalization_profile,
    normalize_cas,
    normalize_identifier,
    normalize_name,
)


def _substance_id(row: Mapping[str, str], row_number: int) -> str:
    explicit = str(row.get("substance_id") or "").strip()
    if explicit:
        return explicit
    cas = normalize_cas(str(row.get("cas") or "").strip())
    return f"cas:{cas}" if cas else f"registry-row:{row_number}"


def _validate_substance_row(
    row: Dict[str, str],
    *,
    row_number: int,
    definition_name: str,
    is_addition: bool,
    roles: set[str],
    families: Dict[str, str],
    known_substance_ids: set[str],
    cas_counts: Counter[str],
    name_counts: Counter[str],
    identifier_owners: DefaultDict[Tuple[str, str], set[str]],
) -> Dict[str, Any] | None:
    row_issues = []
    name = str(row.get("name") or "").strip()
    cas_raw = str(row.get("cas") or "").strip()
    if not name:
        row_issues.append("MISSING_NAME")
    cas = normalize_cas(cas_raw) if cas_raw else None
    if cas_raw and cas is None:
        row_issues.append("INVALID_CAS")
    smiles = str(row.get("smiles") or "").strip()
    if not cas_raw and normalize_cas(smiles):
        row_issues.append("CAS_IN_SMILES_COLUMN")
    substance_id = _substance_id(row, row_number)
    if is_addition:
        explicit_id = str(row.get("substance_id") or "").strip()
        status = str(row.get("status") or "").strip()
        if not explicit_id:
            row_issues.append("MISSING_SUBSTANCE_ID")
        if not cas_raw:
            row_issues.append("MISSING_CAS")
        if cas and explicit_id != f"cas:{cas}":
            row_issues.append("SUBSTANCE_ID_CAS_MISMATCH")
        if not str(row.get("source") or "").strip():
            row_issues.append("MISSING_SUBSTANCE_SOURCE")
        if status not in {"active", "deprecated"}:
            row_issues.append(f"INVALID_SUBSTANCE_STATUS:{status}")
        if substance_id in known_substance_ids:
            row_issues.append(f"DUPLICATE_SUBSTANCE_ID:{substance_id}")
        if cas and cas_counts[cas]:
            row_issues.append(f"DUPLICATE_CAS:{cas}")
    for index in (1, 2):
        role = str(row.get(f"role_{index}") or "").strip()
        family = str(row.get(f"family_{index}") or "").strip()
        if role and role not in roles:
            row_issues.append(f"UNKNOWN_ROLE:{role}")
        if family and family not in families:
            row_issues.append(f"UNKNOWN_FAMILY:{family}")
        if role and family and families.get(family) != role:
            row_issues.append(f"ROLE_FAMILY_MISMATCH:{role}:{family}")
    known_substance_ids.add(substance_id)
    if cas:
        cas_counts[cas] += 1
        identifier_owners[("cas", cas)].add(substance_id)
    normalized_name = normalize_name(name)
    if normalized_name:
        name_counts[normalized_name] += 1
    canonical_key = normalize_identifier(name, "canonical_name")
    if canonical_key:
        identifier_owners[("name", canonical_key)].add(substance_id)
    abbreviation = str(row.get("abbreviation") or "").strip().strip('"')
    abbreviation_key = normalize_identifier(abbreviation, "abbreviation")
    if abbreviation_key:
        identifier_owners[("name", abbreviation_key)].add(substance_id)
    if not row_issues:
        return None
    return {
        "definition": definition_name,
        "row_number": row_number,
        "name": name,
        "cas_raw": cas_raw,
        "issues": row_issues,
    }


def _validate_identifiers(
    identifier_rows: List[Dict[str, str]],
    *,
    known_substance_ids: set[str],
    identifier_owners: DefaultDict[Tuple[str, str], set[str]],
) -> List[Dict[str, Any]]:
    identifier_issues: List[Dict[str, Any]] = []
    identifier_id_counts: Counter[str] = Counter(
        str(row.get("identifier_id") or "").strip()
        for row in identifier_rows
    )
    supplemental_owners: DefaultDict[Tuple[str, str], set[str]] = defaultdict(set)
    supplemental_value_counts: Counter[Tuple[str, str, str]] = Counter()
    for row in identifier_rows:
        try:
            candidate = row_to_identifier(row)
        except (TypeError, ValueError):
            continue
        normalized = normalize_identifier(candidate.value, candidate.identifier_type)
        if normalized is None:
            continue
        namespace = (
            "name"
            if candidate.identifier_type in CONDITION_NAME_IDENTIFIER_TYPES
            else candidate.identifier_type
        )
        key = (namespace, normalized)
        supplemental_owners[key].add(candidate.substance_id)
        supplemental_value_counts[(namespace, normalized, candidate.substance_id)] += 1
    for row_number, row in enumerate(identifier_rows, start=2):
        row_issues = []
        identifier_id = str(row.get("identifier_id") or "").strip()
        substance_id = str(row.get("substance_id") or "").strip()
        if identifier_id and identifier_id_counts[identifier_id] > 1:
            row_issues.append(f"DUPLICATE_IDENTIFIER_ID:{identifier_id}")
        if substance_id not in known_substance_ids:
            row_issues.append(f"UNKNOWN_SUBSTANCE_ID:{substance_id}")
        if not str(row.get("source") or "").strip():
            row_issues.append("MISSING_IDENTIFIER_SOURCE")
        try:
            identifier = row_to_identifier(row)
        except (TypeError, ValueError) as error:
            row_issues.append(f"INVALID_IDENTIFIER:{error}")
            identifier = None
        if identifier is not None:
            expected_profile = identifier_normalization_profile(
                identifier.identifier_type
            )
            if identifier.normalization_profile != expected_profile:
                row_issues.append(
                    "NORMALIZATION_PROFILE_MISMATCH:"
                    f"{identifier.normalization_profile or ''}:"
                    f"{expected_profile or ''}"
                )
            normalized = normalize_identifier(
                identifier.value, identifier.identifier_type
            )
            if normalized is None:
                row_issues.append(
                    f"INVALID_{identifier.identifier_type.upper()}_VALUE"
                )
            else:
                namespace = (
                    "name"
                    if identifier.identifier_type in CONDITION_NAME_IDENTIFIER_TYPES
                    else identifier.identifier_type
                )
                key = (namespace, normalized)
                owners = identifier_owners[key] | supplemental_owners[key]
                if len(owners) > 1 and not identifier.allow_ambiguous:
                    row_issues.append("UNDECLARED_AMBIGUOUS_IDENTIFIER")
                if identifier.substance_id in identifier_owners[key] or (
                    supplemental_value_counts[
                        (namespace, normalized, identifier.substance_id)
                    ]
                    > 1
                ):
                    row_issues.append("DUPLICATE_IDENTIFIER_VALUE")
        if row_issues:
            identifier_issues.append(
                {
                    "row_number": row_number,
                    "identifier_id": identifier_id,
                    "substance_id": substance_id,
                    "value": str(row.get("value") or ""),
                    "issues": row_issues,
                }
            )
    return identifier_issues


def validate_registry(
    *,
    substances_path: str | Path = SUBSTANCES_PATH,
    additions_path: str | Path = ADDITIONS_PATH,
    identifiers_path: str | Path = IDENTIFIERS_PATH,
) -> Dict[str, Any]:
    """Validate substance rows and normalized one-to-many identifiers."""
    taxonomy = load_taxonomy()
    roles = {str(item.get("id")) for item in taxonomy.get("roles", [])}
    families = {
        str(item.get("id")): str(item.get("role_id"))
        for item in taxonomy.get("families", [])
    }
    issues: List[Dict[str, Any]] = []
    cas_counts: Counter[str] = Counter()
    name_counts: Counter[str] = Counter()
    known_substance_ids: set[str] = set()
    identifier_owners: DefaultDict[Tuple[str, str], set[str]] = defaultdict(set)
    accepted = 0
    definitions: Iterable[
        Tuple[str, Iterable[Dict[str, str]], bool]
    ] = (
        (Path(substances_path).name, iter_substance_rows(substances_path), False),
        (Path(additions_path).name, iter_addition_rows(additions_path), True),
    )
    for definition_name, rows, is_addition in definitions:
        for row_number, row in enumerate(rows, start=2):
            issue = _validate_substance_row(
                row,
                row_number=row_number,
                definition_name=definition_name,
                is_addition=is_addition,
                roles=roles,
                families=families,
                known_substance_ids=known_substance_ids,
                cas_counts=cas_counts,
                name_counts=name_counts,
                identifier_owners=identifier_owners,
            )
            if issue is None:
                accepted += 1
            else:
                issues.append(issue)
    identifier_rows = list(iter_identifier_rows(identifiers_path))
    identifier_issues = _validate_identifiers(
        identifier_rows,
        known_substance_ids=known_substance_ids,
        identifier_owners=identifier_owners,
    )
    substance_issue_counts = Counter(
        issue for item in issues for issue in item["issues"]
    )
    identifier_issue_counts = Counter(
        issue for item in identifier_issues for issue in item["issues"]
    )
    return {
        "schema_version": "1.2",
        "total_rows": accepted + len(issues),
        "accepted_rows": accepted,
        "issue_rows": len(issues),
        "duplicate_cas": sum(
            count - 1 for count in cas_counts.values() if count > 1
        ),
        "duplicate_normalized_names": sum(
            count - 1 for count in name_counts.values() if count > 1
        ),
        "issue_counts": dict(sorted(substance_issue_counts.items())),
        "issues": issues,
        "identifier_total_rows": len(identifier_rows),
        "identifier_accepted_rows": len(identifier_rows) - len(identifier_issues),
        "identifier_issue_rows": len(identifier_issues),
        "identifier_issue_counts": dict(sorted(identifier_issue_counts.items())),
        "identifier_issues": identifier_issues,
        "has_errors": bool(issues or identifier_issues),
    }


__all__ = ["validate_registry"]
