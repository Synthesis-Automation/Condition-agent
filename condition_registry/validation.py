"""Validation and audit of the unified condition-substance registry."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, List, Tuple

from .loader import (
    SUBSTANCES_PATH,
    iter_substance_records,
    load_role_definitions,
    record_to_substance,
)
from .models import CONDITION_NAME_IDENTIFIER_TYPES
from .normalization import normalize_identifier


def validate_registry(
    *,
    substances_path: str | Path = SUBSTANCES_PATH,
) -> Dict[str, Any]:
    """Validate unified records, identifiers, and role capabilities."""
    known_roles = {
        str(item["id"]) for item in load_role_definitions().get("roles", ())
    }
    issues: List[Dict[str, Any]] = []
    identifier_issues: List[Dict[str, Any]] = []
    substance_ids: set[str] = set()
    identifier_ids: set[str] = set()
    capability_ids: set[str] = set()
    name_counts: Counter[str] = Counter()
    cas_counts: Counter[str] = Counter()
    identifier_owners: DefaultDict[Tuple[str, str], set[str]] = defaultdict(set)
    alias_shared_flags: DefaultDict[Tuple[str, str], list[bool]] = defaultdict(list)
    identifier_types: DefaultDict[Tuple[str, str], set[str]] = defaultdict(set)
    total_rows = 0

    try:
        records = tuple(iter_substance_records(substances_path))
    except (OSError, ValueError) as error:
        return {
            "schema_version": "2.0",
            "total_rows": 0,
            "accepted_rows": 0,
            "issue_rows": 1,
            "issues": [{"line_number": 0, "issues": [f"INVALID_JSONL:{error}"]}],
            "identifier_total_rows": 0,
            "identifier_accepted_rows": 0,
            "identifier_issue_rows": 0,
            "identifier_issues": [],
            "duplicate_cas": 0,
            "duplicate_normalized_names": 0,
            "issue_counts": {"INVALID_JSONL": 1},
            "has_errors": True,
        }

    for line_number, payload in records:
        total_rows += 1
        row_issues = []
        try:
            substance = record_to_substance(payload, line_number)
        except (TypeError, ValueError) as error:
            issues.append(
                {
                    "line_number": line_number,
                    "substance_id": str(payload.get("id") or ""),
                    "issues": [f"INVALID_RECORD:{error}"],
                }
            )
            continue
        if substance.substance_id in substance_ids:
            row_issues.append(f"DUPLICATE_SUBSTANCE_ID:{substance.substance_id}")
        substance_ids.add(substance.substance_id)
        if substance.cas:
            cas_counts[substance.cas] += 1
        canonical_key = normalize_identifier(substance.canonical_name, "canonical_name")
        if canonical_key:
            name_counts[canonical_key] += 1
        for identifier in substance.identifiers:
            current_issues = []
            if identifier.identifier_id in identifier_ids:
                current_issues.append(f"DUPLICATE_IDENTIFIER_ID:{identifier.identifier_id}")
            identifier_ids.add(identifier.identifier_id)
            normalized = normalize_identifier(identifier.value, identifier.identifier_type)
            if normalized is None:
                current_issues.append(
                    f"INVALID_{identifier.identifier_type.upper()}_VALUE"
                )
            else:
                namespace = (
                    "name"
                    if identifier.identifier_type in CONDITION_NAME_IDENTIFIER_TYPES
                    else identifier.identifier_type
                )
                key = (namespace, normalized)
                identifier_owners[key].add(substance.substance_id)
                if identifier.identifier_type != "canonical_name":
                    alias_shared_flags[key].append(identifier.allow_ambiguous)
                identifier_types[key].add(identifier.identifier_type)
            if current_issues:
                identifier_issues.append(
                    {
                        "line_number": line_number,
                        "identifier_id": identifier.identifier_id,
                        "substance_id": substance.substance_id,
                        "issues": current_issues,
                    }
                )
        for capability in substance.roles:
            capability_id = capability.capability_id or ""
            if capability_id in capability_ids:
                row_issues.append(f"DUPLICATE_ROLE_CAPABILITY:{capability_id}")
            capability_ids.add(capability_id)
            if capability.role_id not in known_roles:
                row_issues.append(f"UNKNOWN_ROLE:{capability.role_id}")
        if row_issues:
            issues.append(
                {
                    "line_number": line_number,
                    "substance_id": substance.substance_id,
                    "issues": row_issues,
                }
            )

    for key, owners in identifier_owners.items():
        if len(owners) <= 1:
            continue
        if identifier_types[key] == {"canonical_name"}:
            continue
        if not all(alias_shared_flags[key]):
            identifier_issues.append(
                {
                    "line_number": 0,
                    "identifier_id": "",
                    "substance_id": "|".join(sorted(owners)),
                    "issues": ["UNDECLARED_AMBIGUOUS_IDENTIFIER"],
                }
            )
    issue_counts = Counter(
        issue.split(":", 1)[0]
        for row in (*issues, *identifier_issues)
        for issue in row["issues"]
    )
    return {
        "schema_version": "2.0",
        "total_rows": total_rows,
        "accepted_rows": total_rows - len(issues),
        "issue_rows": len(issues),
        "issues": issues,
        "identifier_total_rows": len(identifier_ids),
        "identifier_accepted_rows": len(identifier_ids) - len(identifier_issues),
        "identifier_issue_rows": len(identifier_issues),
        "identifier_issues": identifier_issues,
        "role_capability_rows": len(capability_ids),
        "duplicate_cas": sum(
            count - 1 for count in cas_counts.values() if count > 1
        ),
        "duplicate_normalized_names": sum(
            count - 1 for count in name_counts.values() if count > 1
        ),
        "issue_counts": dict(sorted(issue_counts.items())),
        "has_errors": bool(issues or identifier_issues),
    }


__all__ = ["validate_registry"]
