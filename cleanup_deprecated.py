"""
Cleanup script to remove deprecated v1 data and code references.

This script identifies and optionally removes:
1. Old data directories (rule_db, protocol_db without _v2)
2. Old migration scripts that are no longer needed
3. Old test files referencing v1 paths
4. Deprecated wrapper code

The system now exclusively uses:
- data/rule_db_v2/
- data/protocol_db_v2/
- chemtools/recommend/unified.py (DRFP-based)
- chemtools/schema/* (v2.0 schema)
"""

import os
from pathlib import Path
from typing import List, Tuple
import shutil

# Repository root
REPO_ROOT = Path(__file__).parent


def identify_v1_data_dirs() -> List[Tuple[Path, str]]:
    """Find old v1 data directories."""
    old_dirs = []
    
    data_dir = REPO_ROOT / "data"
    
    # Old rule/protocol databases
    for dir_name in ["rule_db", "protocol_db"]:
        old_path = data_dir / dir_name
        if old_path.exists():
            file_count = len(list(old_path.glob("*.json")))
            old_dirs.append((old_path, f"{dir_name}/ ({file_count} JSON files)"))
    
    return old_dirs


def identify_deprecated_scripts() -> List[Tuple[Path, str]]:
    """Find deprecated migration/helper scripts."""
    deprecated = []
    
    scripts_dir = REPO_ROOT / "scripts"
    
    # Migration scripts (one-time use, no longer needed)
    migration_scripts = [
        "migrate_protocols_to_v2.py",  # Migration complete
        "migrate_rules_to_v2.py",       # Migration complete
        "add_reference_reactions.py",   # Migration utility
    ]
    
    for script_name in migration_scripts:
        script_path = scripts_dir / script_name
        if script_path.exists():
            deprecated.append((script_path, f"Migration script (one-time use)"))
    
    return deprecated


def identify_deprecated_test_files() -> List[Tuple[Path, str]]:
    """Find test files that reference old v1 paths."""
    deprecated = []
    
    # Root-level test files
    root_tests = [
        "test_similarity_rule_selection.py",  # Uses data/rule_db
    ]
    
    for test_name in root_tests:
        test_path = REPO_ROOT / test_name
        if test_path.exists():
            deprecated.append((test_path, "References old rule_db path"))
    
    # Test files in tests/ directory - check if they reference old paths
    tests_dir = REPO_ROOT / "tests"
    if tests_dir.exists():
        for test_file in tests_dir.glob("test_*.py"):
            with open(test_file, 'r', encoding='utf-8') as f:
                content = f.read()
                # Check if file references old paths without _v2
                # Exclude temporary directory usage (Path(tmpdir) / "protocol_db")
                if ('data/rule_db' in content or 'data\\rule_db' in content) and 'rule_db_v2' not in content:
                    deprecated.append((test_file, "References old data/rule_db path"))
                elif ('data/protocol_db' in content or 'data\\protocol_db' in content) and 'protocol_db_v2' not in content:
                    deprecated.append((test_file, "References old data/protocol_db path"))
    
    return deprecated


def identify_deprecated_wrapper_code() -> List[Tuple[Path, str]]:
    """Find wrapper code referencing old paths."""
    deprecated = []
    
    # Check chem_assistant wrapper
    wrapper_path = REPO_ROOT / "chem_assistant" / "chemtools_wrapper.py"
    if wrapper_path.exists():
        with open(wrapper_path, 'r', encoding='utf-8') as f:
            content = f.read()
            if 'rule_db"' in content and 'rule_db_v2' not in content:
                deprecated.append((wrapper_path, "Contains references to rule_db (line 132-140)"))
    
    return deprecated


def print_report():
    """Print cleanup report."""
    print("=" * 80)
    print("DEPRECATED CODE AND DATA CLEANUP REPORT")
    print("=" * 80)
    print()
    
    # V1 data directories
    v1_dirs = identify_v1_data_dirs()
    print("1. OLD DATA DIRECTORIES (v1)")
    print("-" * 80)
    if v1_dirs:
        for path, desc in v1_dirs:
            print(f"  ❌ {path}")
            print(f"     {desc}")
            print(f"     → Replaced by: {path.name}_v2/")
    else:
        print("  ✅ No old v1 data directories found")
    print()
    
    # Deprecated scripts
    scripts = identify_deprecated_scripts()
    print("2. DEPRECATED MIGRATION SCRIPTS")
    print("-" * 80)
    if scripts:
        for path, desc in scripts:
            print(f"  ❌ {path.name}")
            print(f"     {desc}")
    else:
        print("  ✅ No deprecated scripts found")
    print()
    
    # Deprecated tests
    tests = identify_deprecated_test_files()
    print("3. DEPRECATED TEST FILES")
    print("-" * 80)
    if tests:
        for path, desc in tests:
            print(f"  ❌ {path}")
            print(f"     {desc}")
    else:
        print("  ✅ No deprecated test files found")
    print()
    
    # Wrapper code
    wrappers = identify_deprecated_wrapper_code()
    print("4. DEPRECATED WRAPPER CODE")
    print("-" * 80)
    if wrappers:
        for path, desc in wrappers:
            print(f"  ⚠️  {path}")
            print(f"     {desc}")
            print(f"     → Should update to use rule_db_v2")
    else:
        print("  ✅ No deprecated wrapper code found")
    print()
    
    # Summary
    total_issues = len(v1_dirs) + len(scripts) + len(tests) + len(wrappers)
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total items to clean up: {total_issues}")
    print()
    
    if total_issues > 0:
        print("Recommended actions:")
        print()
        if v1_dirs:
            print("  1. Archive old data directories:")
            print("     mkdir -p data/archive")
            for path, _ in v1_dirs:
                print(f"     mv {path} data/archive/")
        
        if scripts:
            print()
            print("  2. Move migration scripts to archive:")
            print("     mkdir -p scripts/archive")
            for path, _ in scripts:
                print(f"     mv {path} scripts/archive/")
        
        if tests:
            print()
            print("  3. Remove or update deprecated tests:")
            for path, _ in tests:
                print(f"     # Review and update: {path}")
        
        if wrappers:
            print()
            print("  4. Update wrapper code to use v2 paths:")
            for path, _ in wrappers:
                print(f"     # Update: {path}")
                print(f"     # Change 'rule_db' → 'rule_db_v2'")
                print(f"     # Change 'protocol_db' → 'protocol_db_v2'")
    else:
        print("✅ No cleanup needed - codebase is clean!")
    print()


def archive_old_data(dry_run=True):
    """Archive old v1 data directories."""
    v1_dirs = identify_v1_data_dirs()
    
    if not v1_dirs:
        print("No old data directories to archive.")
        return
    
    archive_dir = REPO_ROOT / "data" / "archive"
    
    if dry_run:
        print("DRY RUN - Would archive:")
        for path, desc in v1_dirs:
            print(f"  {path} → {archive_dir / path.name}")
        print("\nRun with --execute to actually move files")
    else:
        archive_dir.mkdir(exist_ok=True)
        for path, desc in v1_dirs:
            dest = archive_dir / path.name
            print(f"Archiving: {path} → {dest}")
            shutil.move(str(path), str(dest))
        print(f"\n✅ Archived {len(v1_dirs)} directories to {archive_dir}")


def archive_migration_scripts(dry_run=True):
    """Archive completed migration scripts."""
    scripts = identify_deprecated_scripts()
    
    if not scripts:
        print("No migration scripts to archive.")
        return
    
    archive_dir = REPO_ROOT / "scripts" / "archive"
    
    if dry_run:
        print("DRY RUN - Would archive:")
        for path, desc in scripts:
            print(f"  {path.name}")
        print("\nRun with --execute to actually move files")
    else:
        archive_dir.mkdir(exist_ok=True)
        for path, desc in scripts:
            dest = archive_dir / path.name
            print(f"Archiving: {path.name} → {dest}")
            shutil.move(str(path), str(dest))
        print(f"\n✅ Archived {len(scripts)} scripts to {archive_dir}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Cleanup deprecated v1 data and code"
    )
    parser.add_argument(
        "--action",
        choices=["report", "archive-data", "archive-scripts"],
        default="report",
        help="Action to perform"
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually perform actions (default is dry-run)"
    )
    
    args = parser.parse_args()
    
    if args.action == "report":
        print_report()
    
    elif args.action == "archive-data":
        archive_old_data(dry_run=not args.execute)
    
    elif args.action == "archive-scripts":
        archive_migration_scripts(dry_run=not args.execute)
