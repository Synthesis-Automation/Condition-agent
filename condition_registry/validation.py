"""Validation and audit of package-owned condition registry data."""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List

from .loader import iter_substance_rows, load_taxonomy
from .normalization import normalize_cas, normalize_name


def validate_registry() -> Dict[str, Any]:
    taxonomy = load_taxonomy()
    roles = {str(item.get("id")) for item in taxonomy.get("roles", [])}
    families = {str(item.get("id")): str(item.get("role_id")) for item in taxonomy.get("families", [])}
    issues: List[Dict[str, Any]] = []
    cas_counts: Counter[str] = Counter(); name_counts: Counter[str] = Counter(); accepted = 0
    for row_number, row in enumerate(iter_substance_rows(), start=2):
        row_issues = []
        name = str(row.get("name") or "").strip(); cas_raw = str(row.get("cas") or "").strip()
        if not name: row_issues.append("MISSING_NAME")
        cas = normalize_cas(cas_raw) if cas_raw else None
        if cas_raw and cas is None: row_issues.append("INVALID_CAS")
        smiles = str(row.get("smiles") or "").strip()
        if not cas_raw and normalize_cas(smiles): row_issues.append("CAS_IN_SMILES_COLUMN")
        for index in (1, 2):
            role = str(row.get(f"role_{index}") or "").strip(); family = str(row.get(f"family_{index}") or "").strip()
            if role and role not in roles: row_issues.append(f"UNKNOWN_ROLE:{role}")
            if family and family not in families: row_issues.append(f"UNKNOWN_FAMILY:{family}")
            if role and family and families.get(family) != role: row_issues.append(f"ROLE_FAMILY_MISMATCH:{role}:{family}")
        if cas: cas_counts[cas] += 1
        normalized_name = normalize_name(name)
        if normalized_name: name_counts[normalized_name] += 1
        if row_issues:
            issues.append({"row_number": row_number, "name": name, "cas_raw": cas_raw, "issues": row_issues})
        else:
            accepted += 1
    return {
        "schema_version": "1.0", "total_rows": accepted + len(issues), "accepted_rows": accepted,
        "issue_rows": len(issues), "duplicate_cas": sum(count - 1 for count in cas_counts.values() if count > 1),
        "duplicate_normalized_names": sum(count - 1 for count in name_counts.values() if count > 1),
        "issue_counts": dict(sorted(Counter(issue for item in issues for issue in item["issues"]).items())),
        "issues": issues,
    }


__all__ = ["validate_registry"]
