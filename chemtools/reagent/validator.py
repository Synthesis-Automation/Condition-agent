"""
Reagent Database Validator
==========================

Validates CSV entries against the minimal reagent schema.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Dict, List, Optional


REQUIRED_FIELDS = {"name"}
ROLE_FIELDS = ("role_1", "role_2", "role")
FAMILY_FIELDS = ("family_1", "family_2", "family_id")
TAG_FIELDS = ("tag_1", "tag_2", "tag")
SMILES_FIELDS = ("smiles", "smile")
OPTIONAL_FIELDS = {
    "abbreviation",
    "cas",
    *SMILES_FIELDS,
    *FAMILY_FIELDS,
    *TAG_FIELDS,
    *ROLE_FIELDS,
}

CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")
FAMILY_ID_PATTERN = re.compile(r"^[a-z0-9_]+$")


def _load_csv_entries(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [row for row in reader]


def _entry_roles(entry: Dict[str, Any]) -> List[str]:
    roles = []
    for field in ROLE_FIELDS:
        value = str(entry.get(field) or "").strip()
        if value:
            roles.append(value)
    return roles


def _entry_has_role(entry: Dict[str, Any], role: str) -> bool:
    if not role:
        return False
    return role in _entry_roles(entry)


def validate_entry(
    entry: Dict[str, Any],
    role: Optional[str] = None,
    strict: bool = True,
) -> List[Dict[str, str]]:
    """
    Validate a single reagent entry against the CSV schema.
    """
    issues: List[Dict[str, str]] = []

    for field in REQUIRED_FIELDS:
        if not entry.get(field):
            issues.append({
                "severity": "error",
                "field": field,
                "message": f"Missing required field '{field}'",
            })

    roles = [str(entry.get(field) or "").strip() for field in ROLE_FIELDS]
    roles = [value for value in roles if value]
    if not roles:
        issues.append({
            "severity": "error",
            "field": "role_1",
            "message": "Missing required field 'role_1' (or legacy 'role')",
        })

    if role and roles and role not in roles:
        issues.append({
            "severity": "error",
            "field": "role_1",
            "message": f"Entry roles {roles} do not include expected '{role}'",
        })

    if strict:
        for field in REQUIRED_FIELDS | OPTIONAL_FIELDS:
            if field in entry and entry[field] is not None and not isinstance(entry[field], str):
                issues.append({
                    "severity": "error",
                    "field": field,
                    "message": f"Field '{field}' must be string, got {type(entry[field]).__name__}",
                })

    cas = (entry.get("cas") or "").strip()
    if cas and not cas.startswith("$") and not CAS_PATTERN.match(cas):
        issues.append({
            "severity": "warning",
            "field": "cas",
            "message": f"CAS '{cas}' is not in standard format",
        })

    for field in FAMILY_FIELDS:
        family_id = (entry.get(field) or "").strip()
        if family_id and not FAMILY_ID_PATTERN.match(family_id):
            issues.append({
                "severity": "warning",
                "field": field,
                "message": f"Family id '{family_id}' is not snake_case",
            })

    return issues


def _validate_entries(
    entries: List[Dict[str, str]],
    role: Optional[str],
    strict: bool,
    file_path: Path,
) -> Dict[str, Any]:
    result = {
        "role": role or "*",
        "file": str(file_path),
        "total_entries": 0,
        "valid_entries": 0,
        "invalid_entries": 0,
        "errors": [],
        "warnings": [],
        "entry_issues": {},
    }

    filtered = entries
    if role:
        filtered = [entry for entry in entries if _entry_has_role(entry, role)]

    result["total_entries"] = len(filtered)

    for idx, entry in enumerate(filtered):
        issues = validate_entry(entry, role=role, strict=strict)
        if issues:
            result["entry_issues"][f"index_{idx}"] = issues
            errors = [iss for iss in issues if iss["severity"] == "error"]
            warnings = [iss for iss in issues if iss["severity"] == "warning"]
            if errors:
                result["invalid_entries"] += 1
            else:
                result["valid_entries"] += 1
            result["errors"].extend([f"Entry {idx}: {iss['message']}" for iss in errors])
            result["warnings"].extend([f"Entry {idx}: {iss['message']}" for iss in warnings])
        else:
            result["valid_entries"] += 1

    return result


def validate_database(
    registry_dir: str | Path,
    strict: bool = True,
    roles: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Validate entire reagent database (CSV).
    """
    registry_dir = Path(registry_dir)
    csv_path = registry_dir / "reagents.csv"

    if not registry_dir.exists():
        return {
            "registry_dir": str(registry_dir),
            "error": f"Registry directory not found: {registry_dir}",
        }
    if not csv_path.exists():
        return {
            "registry_dir": str(registry_dir),
            "error": f"Reagent CSV not found: {csv_path}",
        }

    entries = _load_csv_entries(csv_path)
    available_roles = sorted({role for entry in entries for role in _entry_roles(entry)})
    roles_to_check = roles or available_roles

    result = {
        "registry_dir": str(registry_dir),
        "roles_checked": list(roles_to_check),
        "total_files": 1,
        "total_entries": 0,
        "valid_entries": 0,
        "invalid_entries": 0,
        "error_count": 0,
        "warning_count": 0,
        "by_role": {},
    }

    for role in roles_to_check:
        role_result = _validate_entries(entries, role, strict, csv_path)
        result["by_role"][role] = role_result
        result["total_entries"] += role_result["total_entries"]
        result["valid_entries"] += role_result["valid_entries"]
        result["invalid_entries"] += role_result["invalid_entries"]
        result["error_count"] += len(role_result["errors"])
        result["warning_count"] += len(role_result["warnings"])

    return result


def validate_role_file(
    registry_dir: str | Path,
    role: str,
    strict: bool = True,
) -> Dict[str, Any]:
    """
    Validate entries for a specific role in the CSV registry.
    """
    results = validate_database(registry_dir, strict=strict, roles=[role])
    if "by_role" in results and role in results["by_role"]:
        return results["by_role"][role]
    return results


def print_validation_summary(results: Dict[str, Any], verbose: bool = False) -> None:
    """
    Print validation results in human-readable format.
    """
    print("\n" + "=" * 70)
    print("REAGENT DATABASE VALIDATION REPORT")
    print("=" * 70)

    if "error" in results:
        print(f"\nERROR: {results['error']}")
        return

    print(f"\nRegistry: {results['registry_dir']}")
    print(f"Roles checked: {', '.join(results['roles_checked'])}")
    print(f"Total files: {results['total_files']}")
    print(f"Total entries: {results['total_entries']}")

    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Valid entries:   {results['valid_entries']}")
    print(f"Invalid entries: {results['invalid_entries']}")
    print(f"Errors:          {results['error_count']}")
    print(f"Warnings:        {results['warning_count']}")

    if verbose:
        print(f"\n{'='*70}")
        print("DETAILS")
        print(f"{'='*70}")
        for role, role_data in results["by_role"].items():
            if role_data.get("errors") or role_data.get("warnings"):
                print(f"\n[{role}]")
                for err in role_data.get("errors", []):
                    print(f"  ERROR: {err}")
                for warn in role_data.get("warnings", []):
                    print(f"  WARNING: {warn}")


def print_critical_errors_summary(results: Dict[str, Any]) -> None:
    """
    Print only critical errors from validation results.
    """
    if "error" in results:
        print(f"\nERROR: {results['error']}")
        return

    errors_found = False
    for role, role_data in results.get("by_role", {}).items():
        errors = role_data.get("errors", [])
        if not errors:
            continue
        errors_found = True
        print(f"\n[{role}]")
        for err in errors:
            print(f"  ERROR: {err}")

    if not errors_found:
        print("\nNo critical errors found.")
