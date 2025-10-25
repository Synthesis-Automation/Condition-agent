"""
Reagent Database Validator
==========================

Validates reagent entries against the schema definition.
Checks for required fields, correct data types, and proper structure.

Usage:
    from chemtools.reagent import validate_database, validate_entry
    
    # Validate entire database
    results = validate_database("data/reagent_db")
    print(f"Valid: {results['valid_count']}, Invalid: {results['invalid_count']}")
    
    # Validate single entry
    entry = {...}
    issues = validate_entry(entry, role="ligand")
    if issues:
        print("Validation errors:", issues)
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from .constants import ROLE_FILES


# Required fields for all entries
REQUIRED_FIELDS = {
    "id",
    "name",
    "abbreviation",
    "aliases",
    "cas",
    "inchi_key",
    "smiles",
    "roles",
    "embedding_text",
}


# Role-specific required fields in roles.<role> object
ROLE_REQUIRED_FIELDS = {
    "ligand": {"families", "donors", "denticity"},
    "base": {"families", "basicity", "nucleophilicity", "sterics"},
    "acid": {"families", "acidity"},
    "oxidant": {"families", "strength_band"},
    "reductant": {"families", "strength_band"},
    "solvent": {"families", "proticity", "polarity", "coordination"},
    "additive": {"families"},
    "metal_precursor": {"families", "metal", "oxidation_states"},
    "preformed_metal_catalyst": {"families", "metal", "oxidation_states", "ligand_type"},
    "organo_catalyst": {"families", "activation_mode", "chirality"},
    "enzyme": {"families", "source", "cofactor_requirement"},
    "condensation_agent": {"families", "strength_band"},
    "other_reagent": {"families"},
}


def validate_entry(
    entry: Dict[str, Any],
    role: Optional[str] = None,
    strict: bool = True,
) -> List[Dict[str, str]]:
    """
    Validate a single reagent entry against schema.
    
    Args:
        entry: Entry dictionary to validate
        role: Expected role (if None, extracts from entry)
        strict: If True, check types and formats strictly
        
    Returns:
        List of validation issues. Empty list means valid.
        Each issue is a dict with 'severity', 'field', 'message'.
    """
    issues: List[Dict[str, str]] = []
    
    # Check required top-level fields
    missing_fields = REQUIRED_FIELDS - set(entry.keys())
    for field in missing_fields:
        issues.append({
            "severity": "error",
            "field": field,
            "message": f"Missing required field '{field}'"
        })
    
    # Validate field types
    if strict:
        # id: string
        if "id" in entry and entry["id"] is not None and not isinstance(entry["id"], str):
            issues.append({
                "severity": "error",
                "field": "id",
                "message": f"Field 'id' must be string, got {type(entry['id']).__name__}"
            })
        
        # name: string
        if "name" in entry and entry["name"] is not None and not isinstance(entry["name"], str):
            issues.append({
                "severity": "error",
                "field": "name",
                "message": f"Field 'name' must be string, got {type(entry['name']).__name__}"
            })
        
        # abbreviation: array
        if "abbreviation" in entry and entry["abbreviation"] is not None:
            if not isinstance(entry["abbreviation"], list):
                issues.append({
                    "severity": "error",
                    "field": "abbreviation",
                    "message": f"Field 'abbreviation' must be array, got {type(entry['abbreviation']).__name__}"
                })
            elif not all(isinstance(x, str) for x in entry["abbreviation"]):
                issues.append({
                    "severity": "error",
                    "field": "abbreviation",
                    "message": "All items in 'abbreviation' must be strings"
                })
        
        # aliases: array
        if "aliases" in entry and entry["aliases"] is not None:
            if not isinstance(entry["aliases"], list):
                issues.append({
                    "severity": "error",
                    "field": "aliases",
                    "message": f"Field 'aliases' must be array, got {type(entry['aliases']).__name__}"
                })
            elif not all(isinstance(x, str) for x in entry["aliases"]):
                issues.append({
                    "severity": "error",
                    "field": "aliases",
                    "message": "All items in 'aliases' must be strings"
                })
        
        # cas: string or null
        if "cas" in entry and entry["cas"] is not None and not isinstance(entry["cas"], str):
            issues.append({
                "severity": "error",
                "field": "cas",
                "message": f"Field 'cas' must be string or null, got {type(entry['cas']).__name__}"
            })
        
        # roles: object
        if "roles" in entry:
            if not isinstance(entry["roles"], dict):
                issues.append({
                    "severity": "error",
                    "field": "roles",
                    "message": f"Field 'roles' must be object, got {type(entry['roles']).__name__}"
                })
    
    # Validate roles structure
    if "roles" in entry and isinstance(entry["roles"], dict):
        roles_obj = entry["roles"]
        
        # Determine which role to validate
        if role:
            roles_to_check = [role] if role in roles_obj else []
        else:
            roles_to_check = list(roles_obj.keys())
        
        if not roles_to_check and role:
            issues.append({
                "severity": "error",
                "field": "roles",
                "message": f"Entry missing expected role '{role}'"
            })
        
        for role_key in roles_to_check:
            role_data = roles_obj.get(role_key)
            if not isinstance(role_data, dict):
                issues.append({
                    "severity": "error",
                    "field": f"roles.{role_key}",
                    "message": f"Role '{role_key}' must be object, got {type(role_data).__name__}"
                })
                continue
            
            # Check role-specific required fields
            required = ROLE_REQUIRED_FIELDS.get(role_key, set())
            missing = required - set(role_data.keys())
            for field in missing:
                issues.append({
                    "severity": "error",
                    "field": f"roles.{role_key}.{field}",
                    "message": f"Missing required field '{field}' for role '{role_key}'"
                })
            
            # Validate families field (required for all roles)
            if "families" in role_data:
                families = role_data["families"]
                if not isinstance(families, list):
                    issues.append({
                        "severity": "error",
                        "field": f"roles.{role_key}.families",
                        "message": f"Field 'families' must be array, got {type(families).__name__}"
                    })
                elif not families:
                    issues.append({
                        "severity": "warning",
                        "field": f"roles.{role_key}.families",
                        "message": "Field 'families' is empty array"
                    })
    
    # Validate ID format (warning only)
    if "id" in entry and entry["id"]:
        entry_id = str(entry["id"])
        # Check if it looks like InChIKey (preferred format)
        if len(entry_id) == 27 and entry_id.count("-") == 2:
            # Valid InChIKey format
            pass
        elif "-" in entry_id and len(entry_id.split("-")) == 3:
            # Looks like CAS number (acceptable)
            pass
        else:
            issues.append({
                "severity": "warning",
                "field": "id",
                "message": f"ID format '{entry_id}' is non-standard (prefer InChIKey or CAS)"
            })
    
    return issues


def validate_role_file(
    file_path: Path,
    role: str,
    strict: bool = True,
) -> Dict[str, Any]:
    """
    Validate all entries in a role file.
    
    Args:
        file_path: Path to role JSON file
        role: Role name (ligand, base, etc.)
        strict: If True, check types strictly
        
    Returns:
        {
            "role": "ligand",
            "file": "/path/to/ligand.json",
            "total_entries": 100,
            "valid_entries": 95,
            "invalid_entries": 5,
            "errors": [...],
            "warnings": [...],
            "entry_issues": {
                "index_0": [...],
                "index_5": [...]
            }
        }
    """
    result = {
        "role": role,
        "file": str(file_path),
        "total_entries": 0,
        "valid_entries": 0,
        "invalid_entries": 0,
        "errors": [],
        "warnings": [],
        "entry_issues": {},
    }
    
    if not file_path.exists():
        result["errors"].append(f"File not found: {file_path}")
        return result
    
    try:
        content = file_path.read_text(encoding="utf-8")
        entries = json.loads(content)
    except json.JSONDecodeError as e:
        result["errors"].append(f"Invalid JSON: {e}")
        return result
    except Exception as e:
        result["errors"].append(f"Failed to read file: {e}")
        return result
    
    if not isinstance(entries, list):
        result["errors"].append(f"File must contain array, got {type(entries).__name__}")
        return result
    
    result["total_entries"] = len(entries)
    
    for idx, entry in enumerate(entries):
        if not isinstance(entry, dict):
            result["entry_issues"][f"index_{idx}"] = [{
                "severity": "error",
                "field": "entry",
                "message": f"Entry must be object, got {type(entry).__name__}"
            }]
            result["invalid_entries"] += 1
            continue
        
        issues = validate_entry(entry, role=role, strict=strict)
        
        if issues:
            result["entry_issues"][f"index_{idx}"] = issues
            
            # Count errors and warnings
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
    Validate entire reagent database.
    
    Args:
        registry_dir: Path to registry directory
        strict: If True, check types strictly
        roles: List of roles to validate (None = all)
        
    Returns:
        {
            "registry_dir": "/path/to/reagents",
            "roles_checked": ["ligand", "base", ...],
            "total_files": 5,
            "total_entries": 500,
            "valid_entries": 490,
            "invalid_entries": 10,
            "error_count": 15,
            "warning_count": 25,
            "by_role": {
                "ligand": {...},
                "base": {...}
            }
        }
    """
    registry_dir = Path(registry_dir)
    
    if not registry_dir.exists():
        return {
            "registry_dir": str(registry_dir),
            "error": f"Registry directory not found: {registry_dir}",
        }
    
    # Determine which roles to check
    if roles:
        roles_to_check = {r: ROLE_FILES.get(r) for r in roles if r in ROLE_FILES}
    else:
        roles_to_check = ROLE_FILES.copy()
    
    result = {
        "registry_dir": str(registry_dir),
        "roles_checked": list(roles_to_check.keys()),
        "total_files": 0,
        "total_entries": 0,
        "valid_entries": 0,
        "invalid_entries": 0,
        "error_count": 0,
        "warning_count": 0,
        "by_role": {},
    }
    
    for role, filename in roles_to_check.items():
        file_path = registry_dir / filename
        
        # Try alternate naming convention if primary not found
        if not file_path.exists():
            # Try simple name: "ligand.json" instead of "taxonomy_ligand.json"
            simple_name = f"{role}.json"
            alt_path = registry_dir / simple_name
            if alt_path.exists():
                file_path = alt_path
            else:
                result["by_role"][role] = {
                    "error": f"File not found: {file_path} (also tried {alt_path})",
                    "file": str(file_path),
                }
                continue
        
        result["total_files"] += 1
        
        role_result = validate_role_file(file_path, role, strict=strict)
        result["by_role"][role] = role_result
        
        result["total_entries"] += role_result["total_entries"]
        result["valid_entries"] += role_result["valid_entries"]
        result["invalid_entries"] += role_result["invalid_entries"]
        result["error_count"] += len(role_result["errors"])
        result["warning_count"] += len(role_result["warnings"])
    
    return result


def print_validation_summary(results: Dict[str, Any], verbose: bool = False) -> None:
    """
    Print validation results in human-readable format.
    
    Args:
        results: Results from validate_database()
        verbose: If True, show all issues
    """
    print("\n" + "="*70)
    print("REAGENT DATABASE VALIDATION REPORT")
    print("="*70)
    
    if "error" in results:
        print(f"\n❌ ERROR: {results['error']}")
        return
    
    print(f"\nRegistry: {results['registry_dir']}")
    print(f"Roles checked: {', '.join(results['roles_checked'])}")
    print(f"Total files: {results['total_files']}")
    print(f"Total entries: {results['total_entries']}")
    
    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"{'='*70}")
    
    valid_pct = (results['valid_entries'] / results['total_entries'] * 100) if results['total_entries'] > 0 else 0
    
    print(f"✅ Valid entries:   {results['valid_entries']:5d} ({valid_pct:.1f}%)")
    print(f"❌ Invalid entries: {results['invalid_entries']:5d}")
    print(f"🔴 Errors:          {results['error_count']:5d}")
    print(f"⚠️  Warnings:        {results['warning_count']:5d}")
    
    print(f"\n{'='*70}")
    print(f"BY ROLE")
    print(f"{'='*70}")
    
    for role, role_result in results["by_role"].items():
        if "error" in role_result:
            print(f"\n{role.upper()}: ❌ {role_result['error']}")
            continue
        
        total = role_result["total_entries"]
        valid = role_result["valid_entries"]
        invalid = role_result["invalid_entries"]
        
        status = "✅" if invalid == 0 else "⚠️" if invalid < total * 0.1 else "❌"
        
        print(f"\n{status} {role.upper()}")
        print(f"   Total: {total}, Valid: {valid}, Invalid: {invalid}")
        print(f"   Errors: {len(role_result['errors'])}, Warnings: {len(role_result['warnings'])}")
        
        if verbose and role_result.get("entry_issues"):
            print(f"\n   Issues:")
            for entry_key, issues in role_result["entry_issues"].items():
                print(f"   - {entry_key}:")
                for issue in issues:
                    severity_icon = "🔴" if issue["severity"] == "error" else "⚠️"
                    print(f"     {severity_icon} {issue['field']}: {issue['message']}")
    
    print(f"\n{'='*70}\n")


def print_critical_errors_summary(results: Dict[str, Any]) -> None:
    """
    Print a focused summary of critical errors only (not warnings).
    Shows compound names by re-reading the files.
    
    Args:
        results: Results from validate_database()
    """
    if "error" in results:
        return
    
    # Collect all critical errors with file info
    critical_errors = []
    
    for role, role_result in results["by_role"].items():
        if "error" in role_result or not role_result.get("entry_issues"):
            continue
        
        file_path = Path(role_result["file"])
        
        # Load entries to get names
        entries = []
        if file_path.exists():
            try:
                content = file_path.read_text(encoding="utf-8")
                entries = json.loads(content)
            except:
                pass
        
        for entry_key, issues in role_result["entry_issues"].items():
            for issue in issues:
                if issue["severity"] == "error":
                    # Extract index from "index_16"
                    idx = int(entry_key.split("_")[1]) if "_" in entry_key else -1
                    
                    # Get compound name
                    compound_name = "Unknown"
                    compound_id = "Unknown"
                    if 0 <= idx < len(entries):
                        entry = entries[idx]
                        compound_name = entry.get("name", "Unknown")
                        compound_id = entry.get("id", "Unknown")
                        # Try to get abbreviation if name is too long
                        abbr = entry.get("abbreviation", [])
                        if abbr and isinstance(abbr, list) and len(abbr) > 0:
                            if len(compound_name) > 40:
                                compound_name = f"{abbr[0]} ({compound_name[:37]}...)"
                            else:
                                compound_name = f"{abbr[0]} ({compound_name})"
                    
                    critical_errors.append({
                        "role": role.upper(),
                        "entry_key": entry_key,
                        "field": issue["field"],
                        "message": issue["message"],
                        "compound_name": compound_name,
                        "compound_id": compound_id,
                    })
    
    if not critical_errors:
        print(f"\n{'='*70}")
        print("✅ NO CRITICAL ERRORS - All entries are valid!")
        print(f"{'='*70}\n")
        return
    
    print(f"\n{'='*70}")
    print(f"🔴 CRITICAL ERRORS SUMMARY ({len(critical_errors)} errors)")
    print(f"{'='*70}")
    print("\nThe following compounds have CRITICAL errors that must be fixed:\n")
    
    # Group by role
    errors_by_role = {}
    for error in critical_errors:
        role = error["role"]
        if role not in errors_by_role:
            errors_by_role[role] = []
        errors_by_role[role].append(error)
    
    for role, errors in errors_by_role.items():
        print(f"\n📌 {role} ({len(errors)} error{'s' if len(errors) > 1 else ''}):")
        for i, error in enumerate(errors, 1):
            print(f"   {i}. {error['compound_name']}")
            print(f"      ID: {error['compound_id']}")
            print(f"      ❌ {error['field']}: {error['message']}")
            if i < len(errors):
                print()
    
    print(f"\n{'='*70}")
    print("💡 For detailed information, run with --verbose flag:")
    print(f"   python scripts/validate_reagent_db.py --verbose")
    print(f"\n📄 See REAGENT_ERRORS_REPORT.md for fix instructions")
    print(f"{'='*70}\n")


# CLI interface
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python -m chemtools.reagent.validator <registry_dir> [--strict] [--verbose]")
        print("\nExample:")
        print("  python -m chemtools.reagent.validator data/reagent_db")
        print("  python -m chemtools.reagent.validator data/reagent_db --verbose")
        sys.exit(1)
    
    registry_dir = sys.argv[1]
    strict = "--strict" in sys.argv or "-s" in sys.argv
    verbose = "--verbose" in sys.argv or "-v" in sys.argv
    
    results = validate_database(registry_dir, strict=strict)
    print_validation_summary(results, verbose=verbose)
    
    # Exit code based on results
    if results.get("error") or results.get("invalid_entries", 0) > 0:
        sys.exit(1)
    else:
        sys.exit(0)
