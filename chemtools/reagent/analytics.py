"""
Reagent Database Analytics
===========================

Provides statistical analysis and insights about the reagent database.

Usage:
    from chemtools.reagent import get_database_statistics, print_database_summary
    
    # Get detailed statistics
    stats = get_database_statistics()
    print(f"Total reagents: {stats['total_reagents']}")
    
    # Pretty print summary
    print_database_summary()
"""

from __future__ import annotations

from typing import Dict, List, Any, Optional
from pathlib import Path
from collections import Counter, defaultdict

from .lookup import get_all_reagent_types, load_reagent_database, get_data_dir
from .constants import ROLE_FILES


def get_database_statistics(registry_dir: Optional[str | Path] = None) -> Dict[str, Any]:
    """
    Get comprehensive statistics about the reagent database.
    
    Args:
        registry_dir: Path to reagent registry (default: data/reagent_db)
        
    Returns:
        Dictionary with statistics:
        {
            "total_reagents": 383,
            "by_type": {"ligand": 155, "base": 47, ...},
            "total_with_cas": 350,
            "total_with_inchikey": 200,
            "total_with_smiles": 100,
            "families_by_role": {...},
            "top_families": [...],
            ...
        }
    """
    if registry_dir is None:
        registry_dir = get_data_dir() / "reagents"
    else:
        registry_dir = Path(registry_dir)
    
    stats = {
        "registry_dir": str(registry_dir),
        "total_reagents": 0,
        "by_type": {},
        "total_with_cas": 0,
        "total_with_inchikey": 0,
        "total_with_smiles": 0,
        "total_with_abbreviations": 0,
        "families_by_role": defaultdict(list),
        "family_counts": Counter(),
        "top_families": [],
        "id_format_stats": {
            "inchikey": 0,
            "cas": 0,
            "other": 0
        },
        "multi_role_reagents": 0,
        "roles_per_reagent": Counter(),
    }
    
    # Get all available types
    types = get_all_reagent_types()
    
    for reagent_type in types:
        db = load_reagent_database(reagent_type)
        
        # Skip if db is not a list (e.g., taxonomy files which are role definitions)
        if not isinstance(db, list):
            continue
        
        count = len(db)
        stats["by_type"][reagent_type] = count
        stats["total_reagents"] += count
        
        # Analyze each entry
        for entry in db:
            # Check field presence
            if entry.get("cas"):
                stats["total_with_cas"] += 1
            
            if entry.get("inchi_key"):
                stats["total_with_inchikey"] += 1
            
            if entry.get("smiles"):
                stats["total_with_smiles"] += 1
            
            if entry.get("abbreviation"):
                stats["total_with_abbreviations"] += 1
            
            # ID format analysis
            entry_id = entry.get("id", "")
            if entry_id:
                if len(entry_id) == 27 and entry_id.count("-") == 2:
                    stats["id_format_stats"]["inchikey"] += 1
                elif "-" in entry_id and entry_id.replace("-", "").isdigit():
                    stats["id_format_stats"]["cas"] += 1
                else:
                    stats["id_format_stats"]["other"] += 1
            
            # Role and family analysis
            roles = entry.get("roles", {})
            num_roles = len(roles)
            stats["roles_per_reagent"][num_roles] += 1
            
            if num_roles > 1:
                stats["multi_role_reagents"] += 1
            
            for role, role_data in roles.items():
                if isinstance(role_data, dict):
                    families = role_data.get("families", [])
                    for family in families:
                        stats["families_by_role"][role].append(family)
                        stats["family_counts"][family] += 1
    
    # Convert defaultdict to regular dict and get unique families
    stats["families_by_role"] = {
        role: list(set(families))
        for role, families in stats["families_by_role"].items()
    }
    
    # Get top 20 families
    stats["top_families"] = stats["family_counts"].most_common(20)
    
    # Convert Counter to dict for JSON serialization
    stats["roles_per_reagent"] = dict(stats["roles_per_reagent"])
    
    # Calculate percentages
    total = stats["total_reagents"]
    if total > 0:
        stats["percent_with_cas"] = (stats["total_with_cas"] / total) * 100
        stats["percent_with_inchikey"] = (stats["total_with_inchikey"] / total) * 100
        stats["percent_with_smiles"] = (stats["total_with_smiles"] / total) * 100
        stats["percent_with_abbreviations"] = (stats["total_with_abbreviations"] / total) * 100
    
    return stats


def get_family_statistics(role: str, registry_dir: Optional[str | Path] = None) -> Dict[str, Any]:
    """
    Get statistics about families for a specific role.
    
    Args:
        role: Reagent role (ligand, base, etc.)
        registry_dir: Path to reagent registry
        
    Returns:
        {
            "role": "ligand",
            "total_reagents": 155,
            "total_families": 25,
            "families": [
                {"name": "phosphine", "count": 45},
                {"name": "nho_carbene", "count": 20},
                ...
            ]
        }
    """
    if registry_dir is None:
        registry_dir = get_data_dir() / "reagents"
    else:
        registry_dir = Path(registry_dir)
    
    db = load_reagent_database(role)
    
    family_counts = Counter()
    
    for entry in db:
        roles = entry.get("roles", {})
        role_data = roles.get(role, {})
        if isinstance(role_data, dict):
            families = role_data.get("families", [])
            for family in families:
                family_counts[family] += 1
    
    return {
        "role": role,
        "total_reagents": len(db),
        "total_families": len(family_counts),
        "families": [
            {"name": name, "count": count}
            for name, count in family_counts.most_common()
        ]
    }


def find_reagents_by_family(role: str, family: str) -> List[Dict[str, Any]]:
    """
    Find all reagents in a specific family.
    
    Args:
        role: Reagent role
        family: Family name
        
    Returns:
        List of reagent entries in that family
    """
    db = load_reagent_database(role)
    
    matches = []
    for entry in db:
        roles = entry.get("roles", {})
        role_data = roles.get(role, {})
        if isinstance(role_data, dict):
            families = role_data.get("families", [])
            if family in families:
                matches.append(entry)
    
    return matches


def get_missing_data_report(registry_dir: Optional[str | Path] = None) -> Dict[str, Any]:
    """
    Generate report of missing data in database.
    
    Args:
        registry_dir: Path to reagent registry
        
    Returns:
        {
            "by_field": {
                "cas": {"missing": 33, "percent": 8.6},
                "inchi_key": {"missing": 183, "percent": 47.8},
                ...
            },
            "by_type": {
                "ligand": {"missing_cas": 5, "missing_inchikey": 40, ...},
                ...
            }
        }
    """
    if registry_dir is None:
        registry_dir = get_data_dir() / "reagents"
    else:
        registry_dir = Path(registry_dir)
    
    types = get_all_reagent_types()
    
    total_count = 0
    missing_by_field = {
        "cas": 0,
        "inchi_key": 0,
        "smiles": 0,
        "abbreviation": 0,
        "aliases": 0,
    }
    
    missing_by_type = {}
    
    for reagent_type in types:
        db = load_reagent_database(reagent_type)
        total_count += len(db)
        
        type_missing = {
            "total": len(db),
            "missing_cas": 0,
            "missing_inchikey": 0,
            "missing_smiles": 0,
            "missing_abbreviation": 0,
            "missing_aliases": 0,
        }
        
        for entry in db:
            if not entry.get("cas"):
                missing_by_field["cas"] += 1
                type_missing["missing_cas"] += 1
            
            if not entry.get("inchi_key"):
                missing_by_field["inchi_key"] += 1
                type_missing["missing_inchikey"] += 1
            
            if not entry.get("smiles"):
                missing_by_field["smiles"] += 1
                type_missing["missing_smiles"] += 1
            
            if not entry.get("abbreviation"):
                missing_by_field["abbreviation"] += 1
                type_missing["missing_abbreviation"] += 1
            
            if not entry.get("aliases"):
                missing_by_field["aliases"] += 1
                type_missing["missing_aliases"] += 1
        
        missing_by_type[reagent_type] = type_missing
    
    # Calculate percentages
    by_field = {}
    for field, count in missing_by_field.items():
        by_field[field] = {
            "missing": count,
            "percent": (count / total_count * 100) if total_count > 0 else 0
        }
    
    return {
        "total_reagents": total_count,
        "by_field": by_field,
        "by_type": missing_by_type,
    }


def print_database_summary(registry_dir: Optional[str | Path] = None) -> None:
    """
    Print a formatted summary of the reagent database.
    
    Args:
        registry_dir: Path to reagent registry
    """
    stats = get_database_statistics(registry_dir)
    
    print("\n" + "="*70)
    print("REAGENT DATABASE SUMMARY")
    print("="*70)
    
    print(f"\nRegistry: {stats['registry_dir']}")
    print(f"Total reagents: {stats['total_reagents']}")
    
    print(f"\n{'='*70}")
    print("BY TYPE")
    print("="*70)
    
    for reagent_type, count in sorted(stats["by_type"].items(), key=lambda x: -x[1]):
        print(f"  {reagent_type:20s}: {count:4d} reagents")
    
    print(f"\n{'='*70}")
    print("DATA COMPLETENESS")
    print("="*70)
    
    print(f"  CAS numbers:       {stats['total_with_cas']:4d} ({stats.get('percent_with_cas', 0):.1f}%)")
    print(f"  InChIKeys:         {stats['total_with_inchikey']:4d} ({stats.get('percent_with_inchikey', 0):.1f}%)")
    print(f"  SMILES:            {stats['total_with_smiles']:4d} ({stats.get('percent_with_smiles', 0):.1f}%)")
    print(f"  Abbreviations:     {stats['total_with_abbreviations']:4d} ({stats.get('percent_with_abbreviations', 0):.1f}%)")
    
    print(f"\n{'='*70}")
    print("ID FORMATS")
    print("="*70)
    
    id_stats = stats["id_format_stats"]
    print(f"  InChIKey format:   {id_stats['inchikey']:4d}")
    print(f"  CAS format:        {id_stats['cas']:4d}")
    print(f"  Other format:      {id_stats['other']:4d}")
    
    print(f"\n{'='*70}")
    print("MULTI-ROLE REAGENTS")
    print("="*70)
    
    print(f"  Reagents with multiple roles: {stats['multi_role_reagents']}")
    for num_roles, count in sorted(stats["roles_per_reagent"].items()):
        if num_roles > 0:
            role_label = "role" if num_roles == 1 else "roles"
            print(f"    {num_roles} {role_label}: {count} reagents")
    
    print(f"\n{'='*70}")
    print("TOP 10 FAMILIES (across all roles)")
    print("="*70)
    
    for family, count in stats["top_families"][:10]:
        print(f"  {family:30s}: {count:3d} reagents")
    
    print(f"\n{'='*70}")
    print("FAMILIES BY ROLE")
    print("="*70)
    
    for role in sorted(stats["families_by_role"].keys()):
        families = stats["families_by_role"][role]
        print(f"  {role:20s}: {len(families):2d} families")
    
    print(f"\n{'='*70}\n")


def print_role_summary(role: str, registry_dir: Optional[str | Path] = None) -> None:
    """
    Print detailed summary for a specific role.
    
    Args:
        role: Reagent role
        registry_dir: Path to reagent registry
    """
    stats = get_family_statistics(role, registry_dir)
    
    print("\n" + "="*70)
    print(f"{role.upper()} DATABASE SUMMARY")
    print("="*70)
    
    print(f"\nTotal {role}s: {stats['total_reagents']}")
    print(f"Total families: {stats['total_families']}")
    
    print(f"\n{'='*70}")
    print("FAMILIES")
    print("="*70)
    
    for family_data in stats["families"]:
        name = family_data["name"]
        count = family_data["count"]
        percent = (count / stats['total_reagents'] * 100) if stats['total_reagents'] > 0 else 0
        print(f"  {name:40s}: {count:3d} ({percent:4.1f}%)")
    
    print(f"\n{'='*70}\n")


def print_missing_data_report(registry_dir: Optional[str | Path] = None) -> None:
    """
    Print report of missing data in database.
    
    Args:
        registry_dir: Path to reagent registry
    """
    report = get_missing_data_report(registry_dir)
    
    print("\n" + "="*70)
    print("MISSING DATA REPORT")
    print("="*70)
    
    print(f"\nTotal reagents: {report['total_reagents']}")
    
    print(f"\n{'='*70}")
    print("MISSING BY FIELD")
    print("="*70)
    
    for field, data in report["by_field"].items():
        print(f"  {field:20s}: {data['missing']:4d} missing ({data['percent']:5.1f}%)")
    
    print(f"\n{'='*70}")
    print("MISSING BY TYPE")
    print("="*70)
    
    for reagent_type, data in sorted(report["by_type"].items()):
        print(f"\n  {reagent_type.upper()} (total: {data['total']})")
        print(f"    Missing CAS:          {data['missing_cas']:3d}")
        print(f"    Missing InChIKey:     {data['missing_inchikey']:3d}")
        print(f"    Missing SMILES:       {data['missing_smiles']:3d}")
        print(f"    Missing abbreviation: {data['missing_abbreviation']:3d}")
        print(f"    Missing aliases:      {data['missing_aliases']:3d}")
    
    print(f"\n{'='*70}\n")


# CLI interface
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python -m chemtools.reagent.analytics <command> [options]")
        print("\nCommands:")
        print("  summary              - Show database summary")
        print("  role <role>          - Show summary for specific role")
        print("  families <role>      - Show family statistics for role")
        print("  missing              - Show missing data report")
        print("\nExamples:")
        print("  python -m chemtools.reagent.analytics summary")
        print("  python -m chemtools.reagent.analytics role ligand")
        print("  python -m chemtools.reagent.analytics missing")
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == "summary":
        print_database_summary()
    elif command == "role" and len(sys.argv) > 2:
        print_role_summary(sys.argv[2])
    elif command == "families" and len(sys.argv) > 2:
        print_role_summary(sys.argv[2])
    elif command == "missing":
        print_missing_data_report()
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)
