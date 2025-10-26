"""
Test Sample Reactions
=====================

Tests all reactions in sample_reactions.py using the analyze_reaction function.
Reports any errors or failures.
"""

import sys
import os
from pathlib import Path
from typing import Dict, List, Any

# Set UTF-8 encoding for Windows console
if sys.platform == "win32":
    os.environ["PYTHONIOENCODING"] = "utf-8"

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "tests"))

from chemtools.analysis import analyze_reaction
from sample_reactions import SAMPLE_REACTIONS, BUCHWALD_HARTWIG_REACTIONS

def extract_smiles_from_entry(entry: str) -> str:
    """Extract SMILES from a sample reaction entry."""
    if ">>>" in entry or entry.startswith("Select"):
        return None
    
    # Find the last opening parenthesis (description starts there)
    # This avoids splitting SMILES like B(O)O
    last_paren = entry.rfind(" (")
    if last_paren == -1:
        # No description in parentheses - check if it's still a valid reaction
        return entry.strip() if ">>" in entry else None
    
    # Extract SMILES before the description
    smiles = entry[:last_paren].strip()
    return smiles if ">>" in smiles else None

def test_reaction(smiles: str, description: str, index: int) -> Dict[str, Any]:
    """Test a single reaction and return results."""
    result = {
        "index": index,
        "description": description,
        "smiles": smiles,
        "success": False,
        "error": None,
        "family": None,
        "reactants": None
    }
    
    try:
        analysis = analyze_reaction(smiles)
        result["success"] = True
        result["family"] = analysis.get("family", {}).get("canonical_id", "UNKNOWN")
        result["reactants"] = len(analysis.get("reactants", []))
        result["analysis"] = analysis
    except Exception as e:
        result["error"] = str(e)
    
    return result

def test_sample_reactions():
    """Test all SAMPLE_REACTIONS."""
    print("=" * 80)
    print("TESTING SAMPLE_REACTIONS")
    print("=" * 80)
    
    results = []
    errors = []
    successful = 0
    total = 0
    
    for idx, entry in enumerate(SAMPLE_REACTIONS):
        smiles = extract_smiles_from_entry(entry)
        if smiles is None:
            continue
        
        total += 1
        description = entry.split("(")[1].split(")")[0] if "(" in entry else "No description"
        
        result = test_reaction(smiles, description, idx)
        results.append(result)
        
        if result["success"]:
            successful += 1
            family = result.get('family') or 'UNKNOWN'
            
            # Extract reactant types
            reactant_types = []
            analysis = result.get('analysis', {})
            for r in analysis.get('reactants', []):
                tax = r.get('taxonomy', {})
                best = tax.get('best_match')
                if best:
                    reactant_types.append(best.get('member_type', '?'))
                else:
                    reactant_types.append('?')
            
            reactants_str = '+'.join(reactant_types) if reactant_types else 'N/A'
            print(f"OK  [{idx:3d}] {family:20s} | {reactants_str:30s} | {description[:40]}")
        else:
            errors.append(result)
            print(f"ERR [{idx:3d}] ERROR: {result['error'][:80]}")
            print(f"      SMILES: {smiles[:80]}")
            print(f"      DESC:   {description[:80]}")
    
    print("\n" + "=" * 80)
    print(f"SUMMARY: {successful}/{total} reactions passed ({100*successful/total:.1f}%)")
    print("=" * 80)
    
    return results, errors

def test_buchwald_hartwig_reactions():
    """Test all BUCHWALD_HARTWIG_REACTIONS."""
    print("\n" + "=" * 80)
    print("TESTING BUCHWALD_HARTWIG_REACTIONS")
    print("=" * 80)
    
    results = []
    errors = []
    successful = 0
    total = 0
    
    for idx, entry in enumerate(BUCHWALD_HARTWIG_REACTIONS):
        smiles = extract_smiles_from_entry(entry)
        if smiles is None:
            continue
        
        total += 1
        description = entry.split("(")[1].split(")")[0] if "(" in entry else "No description"
        
        result = test_reaction(smiles, description, idx)
        results.append(result)
        
        if result["success"]:
            successful += 1
            family = result.get('family') or 'UNKNOWN'
            
            # Extract reactant types
            reactant_types = []
            analysis = result.get('analysis', {})
            for r in analysis.get('reactants', []):
                tax = r.get('taxonomy', {})
                best = tax.get('best_match')
                if best:
                    reactant_types.append(best.get('member_type', '?'))
                else:
                    reactant_types.append('?')
            
            reactants_str = '+'.join(reactant_types) if reactant_types else 'N/A'
            print(f"OK  [{idx:3d}] {family:20s} | {reactants_str:30s} | {description[:40]}")
        else:
            errors.append(result)
            print(f"ERR [{idx:3d}] ERROR: {result['error'][:80]}")
            print(f"      SMILES: {smiles[:80]}")
            print(f"      DESC:   {description[:80]}")
    
    print("\n" + "=" * 80)
    print(f"SUMMARY: {successful}/{total} reactions passed ({100*successful/total:.1f}%)")
    print("=" * 80)
    
    return results, errors

def print_error_details(errors: List[Dict[str, Any]]):
    """Print detailed error information."""
    if not errors:
        print("\n[OK] No errors found!")
        return
    
    print("\n" + "=" * 80)
    print(f"DETAILED ERROR REPORT ({len(errors)} errors)")
    print("=" * 80)
    
    # Group errors by type
    error_types = {}
    for err in errors:
        error_msg = err["error"]
        if error_msg not in error_types:
            error_types[error_msg] = []
        error_types[error_msg].append(err)
    
    for error_msg, err_list in error_types.items():
        print(f"\n[ERROR] Error Type: {error_msg}")
        print(f"   Count: {len(err_list)}")
        print(f"   Examples:")
        for err in err_list[:3]:  # Show first 3 examples
            print(f"     - [{err['index']:3d}] {err['description'][:60]}")
            print(f"       SMILES: {err['smiles'][:70]}")

def main():
    """Main test function."""
    # Test SAMPLE_REACTIONS
    sample_results, sample_errors = test_sample_reactions()
    
    # Test BUCHWALD_HARTWIG_REACTIONS
    bh_results, bh_errors = test_buchwald_hartwig_reactions()
    
    # Print error details
    if sample_errors:
        print("\n" + "=" * 80)
        print("SAMPLE_REACTIONS ERRORS")
        print("=" * 80)
        print_error_details(sample_errors)
    
    if bh_errors:
        print("\n" + "=" * 80)
        print("BUCHWALD_HARTWIG_REACTIONS ERRORS")
        print("=" * 80)
        print_error_details(bh_errors)
    
    # Final summary
    total_tested = len(sample_results) + len(bh_results)
    total_errors = len(sample_errors) + len(bh_errors)
    total_passed = total_tested - total_errors
    
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Total reactions tested: {total_tested}")
    print(f"Passed:                 {total_passed} ({100*total_passed/total_tested:.1f}%)")
    print(f"Failed:                 {total_errors} ({100*total_errors/total_tested:.1f}%)")
    print("=" * 80)
    
    # Exit with error code if any failures
    return 0 if total_errors == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
