#!/usr/bin/env python3
"""
Validate All Suzuki Database Rules
===================================

This script tests that EVERY rule in suzuki_db.json can actually be reached
by the matching system. For each rule, we:

1. Load a specific test reaction from sample_reactions.py
2. Run the matcher
3. Verify the expected rule ID is returned
4. Generate a comprehensive validation report

Usage:
    python scripts/validate_all_suzuki_rules.py
    
Expected Output:
    - Console summary of pass/fail for each rule
    - Detailed report showing which entry matched
    - Priority conflicts (if multiple rules match)
    - Total coverage statistics
"""

import sys
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "tests"))

from chemtools.scdb_matcher.loader import load_db
from chemtools.scdb_matcher.matcher import match
from sample_reactions import get_suzuki_db_tests


def validate_single_rule(rule_id: str, test_data: dict, db) -> dict:
    """
    Validate that a single rule can be matched.
    
    Args:
        rule_id: The rule ID (e.g., "SCDB-SUZ-ARBR-ORTHO-XPhos")
        test_data: Dict with 'smiles', 'description', 'expected_features'
        db: Loaded database instance
    
    Returns:
        Dict with validation results
    """
    smiles = test_data["smiles"]
    description = test_data["description"]
    
    # Run the matcher
    try:
        result = match(db, rxn_smiles=smiles)
        
        # Check if we got a match - entry_id will be set if matched
        if not result or not result.entry_id:
            return {
                "rule_id": rule_id,
                "status": "FAIL",
                "reason": "No match found (matcher returned no entry_id)",
                "description": description,
                "matched_id": None,
                "priority": None
            }
        
        matched_id = result.entry_id
        priority = result.priority
        
        # Check if we matched the expected rule
        if matched_id == rule_id:
            return {
                "rule_id": rule_id,
                "status": "PASS",
                "reason": f"Correctly matched expected rule",
                "description": description,
                "matched_id": matched_id,
                "priority": priority,
                "conditions": result.conditions
            }
        else:
            return {
                "rule_id": rule_id,
                "status": "FAIL",
                "reason": f"Matched different rule: {matched_id} (priority {priority})",
                "description": description,
                "matched_id": matched_id,
                "priority": priority,
                "conditions": result.conditions
            }
            
    except Exception as e:
        import traceback
        return {
            "rule_id": rule_id,
            "status": "ERROR",
            "reason": f"Exception during matching: {str(e)}",
            "description": description,
            "matched_id": None,
            "priority": None,
            "traceback": traceback.format_exc()
        }


def print_result(result: dict):
    """Pretty print a single validation result"""
    status = result["status"]
    rule_id = result["rule_id"]
    
    # Status indicators
    if status == "PASS":
        status_str = f"[PASS]"
    elif status == "FAIL":
        status_str = f"[FAIL]"
    else:
        status_str = f"[ERROR]"
    
    print(f"\n{status_str} {rule_id}")
    print(f"  Description: {result['description']}")
    print(f"  Reason: {result['reason']}")
    
    if result.get("matched_id"):
        print(f"  Matched ID: {result['matched_id']}")
        print(f"  Priority: {result.get('priority', 'N/A')}")
        
        # Show key conditions if matched
        if result.get("conditions"):
            cond = result["conditions"]
            pd = cond.get("pd_source", ["N/A"])[0] if isinstance(cond.get("pd_source"), list) else cond.get("pd_source", "N/A")
            ligand = cond.get("ligands", ["N/A"])[0] if isinstance(cond.get("ligands"), list) else cond.get("ligands", "N/A")
            temp = cond.get("temperature_C", ["N/A"])[0] if isinstance(cond.get("temperature_C"), list) else cond.get("temperature_C", "N/A")
            print(f"  Conditions: {pd} / {ligand} / {temp}C")


def generate_markdown_report(results: List[dict], output_path: Path):
    """Generate a markdown validation report"""
    
    pass_count = sum(1 for r in results if r["status"] == "PASS")
    fail_count = sum(1 for r in results if r["status"] == "FAIL")
    error_count = sum(1 for r in results if r["status"] == "ERROR")
    total_count = len(results)
    
    pass_rate = (pass_count / total_count * 100) if total_count > 0 else 0
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(f"# Suzuki Database Validation Report\n\n")
        f.write(f"**Generated:** {Path(__file__).name}\n\n")
        
        f.write(f"## Summary Statistics\n\n")
        f.write(f"- **Total Rules:** {total_count}\n")
        f.write(f"- **Passed:** {pass_count} ({pass_rate:.1f}%)\n")
        f.write(f"- **Failed:** {fail_count}\n")
        f.write(f"- **Errors:** {error_count}\n\n")
        
        f.write(f"## Detailed Results\n\n")
        
        # Group by status
        for status in ["PASS", "FAIL", "ERROR"]:
            status_results = [r for r in results if r["status"] == status]
            if not status_results:
                continue
            
            if status == "PASS":
                emoji = "âœ?
            elif status == "FAIL":
                emoji = "âœ?
            else:
                emoji = "âš?
            
            f.write(f"### {emoji} {status} ({len(status_results)})\n\n")
            
            for result in status_results:
                f.write(f"#### {result['rule_id']}\n\n")
                f.write(f"- **Description:** {result['description']}\n")
                f.write(f"- **Reason:** {result['reason']}\n")
                
                if result.get("matched_id"):
                    f.write(f"- **Matched ID:** `{result['matched_id']}`\n")
                    f.write(f"- **Priority:** {result.get('priority', 'N/A')}\n")
                    
                    if result.get("conditions"):
                        cond = result["conditions"]
                        pd = cond.get("pd_source", ["N/A"])[0] if isinstance(cond.get("pd_source"), list) else cond.get("pd_source", "N/A")
                        ligand = cond.get("ligands", ["N/A"])[0] if isinstance(cond.get("ligands"), list) else cond.get("ligands", "N/A")
                        temp = cond.get("temperature_C", ["N/A"])[0] if isinstance(cond.get("temperature_C"), list) else cond.get("temperature_C", "N/A")
                        solv = cond.get("solvent", ["N/A"])[0] if isinstance(cond.get("solvent"), list) else cond.get("solvent", "N/A")
                        
                        f.write(f"- **Conditions:**\n")
                        f.write(f"  - Pd source: `{pd}`\n")
                        f.write(f"  - Ligand: `{ligand}`\n")
                        f.write(f"  - Temperature: {temp}Â°C\n")
                        f.write(f"  - Solvent: `{solv}`\n")
                
                f.write(f"\n")
        
        f.write(f"\n## Recommendations\n\n")
        
        if fail_count > 0:
            f.write(f"### Failed Rules Analysis\n\n")
            failed = [r for r in results if r["status"] == "FAIL"]
            for result in failed:
                f.write(f"- **{result['rule_id']}**: Matched `{result['matched_id']}` instead\n")
                f.write(f"  - **Action:** Check SMARTS patterns, feature requirements, and priority\n")
                f.write(f"  - **Conflict:** Priority {result.get('priority', 'N/A')} beat expected rule\n\n")
        
        if error_count > 0:
            f.write(f"### Error Rules Analysis\n\n")
            errors = [r for r in results if r["status"] == "ERROR"]
            for result in errors:
                f.write(f"- **{result['rule_id']}**: {result['reason']}\n")
                f.write(f"  - **Action:** Check SMARTS syntax and feature detector implementation\n\n")
        
        if pass_rate == 100:
            f.write(f"### âœ?All Rules Validated Successfully!\n\n")
            f.write(f"Every rule in suzuki_db.json can be successfully matched by the system.\n")
            f.write(f"The database is complete and all entries are reachable.\n\n")


def main():
    """Main validation workflow"""
    
    print("=" * 70)
    print("Suzuki Database Validation Test")
    print("=" * 70)
    
    # Paths
    db_path = project_root / "data" / "conditionDB" / "suzuki_db.json"
    report_path = project_root / "docs" / "suzuki_db_validation_report.md"
    
    # Load database
    print(f"\nLoading database: {db_path}")
    db = load_db(db_path)
    print(f"  Loaded {len(db.entries)} entries")
    
    # Load test cases
    test_cases = get_suzuki_db_tests()
    
    print(f"  Loaded {len(test_cases)} test cases from sample_reactions.py")
    print(f"\nStarting validation...\n")
    
    # Validate each rule
    results = []
    for rule_id, test_data in test_cases.items():
        result = validate_single_rule(rule_id, test_data, db)
        results.append(result)
        print_result(result)
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    
    pass_count = sum(1 for r in results if r["status"] == "PASS")
    fail_count = sum(1 for r in results if r["status"] == "FAIL")
    error_count = sum(1 for r in results if r["status"] == "ERROR")
    total_count = len(results)
    
    pass_rate = (pass_count / total_count * 100) if total_count > 0 else 0
    
    print(f"\nTotal Rules Tested: {total_count}")
    print(f"  [PASS]: {pass_count} ({pass_rate:.1f}%)")
    print(f"  [FAIL]: {fail_count}")
    print(f"  [ERROR]: {error_count}")
    
    # Show errors if any
    if error_count > 0:
        print(f"\nERROR DETAILS:")
        for r in results:
            if r["status"] == "ERROR" and "traceback" in r:
                print(f"\n{r['rule_id']}:")
                print(r["traceback"])
    
    # Generate markdown report
    print(f"\nGenerating detailed report: {report_path}")
    generate_markdown_report(results, report_path)
    
    if pass_rate == 100:
        print(f"\n[SUCCESS] ALL TESTS PASSED!")
    else:
        print(f"\n[WARNING] SOME TESTS FAILED")
    print("=" * 70)
    
    # Exit code
    sys.exit(0 if fail_count == 0 and error_count == 0 else 1)


if __name__ == "__main__":
    main()
