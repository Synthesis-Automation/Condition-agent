#!/usr/bin/env python3
"""
Test all sample reactions with cross-family recommendation.
Generates a comprehensive summary of results.
"""

import sys
from pathlib import Path
import json
from collections import defaultdict
import time

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "tests"))

from sample_reactions import SAMPLE_REACTIONS
from chemtools import chem


def extract_reaction_type(reaction_str: str) -> str:
    """Extract reaction type from the description in parentheses."""
    if "(" not in reaction_str:
        return "Unknown"
    
    desc = reaction_str.split("(")[-1].split(")")[0]
    
    # Map to major categories
    if any(x in desc for x in ["Suzuki", "Stille", "Sonogashira", "Heck", "Negishi", "Kumada"]):
        return "C-C Coupling"
    elif any(x in desc for x in ["Buchwald-Hartwig", "B-H", "C-N"]):
        return "C-N Coupling"
    elif "C-O" in desc or "Ullmann Ether" in desc:
        return "C-O Coupling"
    elif "C-S" in desc:
        return "C-S Coupling"
    elif any(x in desc for x in ["Amide", "Amidation", "Esterification"]):
        return "Amide/Ester Formation"
    elif any(x in desc for x in ["Hydrogenation", "NaBH4", "LiAlH4", "reduction", "Transfer Hydrogenation", "Birch", "DIBAL"]):
        return "Reduction"
    elif any(x in desc for x in ["Oxidation", "PCC", "Swern", "MnO2", "Jones", "Dess-Martin", "IBX", "TEMPO", "mCPBA", "Baeyer-Villiger"]):
        return "Oxidation"
    elif any(x in desc for x in ["SN1", "SN2", "SNAr", "Allylic", "Finkelstein", "Appel", "Mitsunobu"]):
        return "Substitution"
    elif any(x in desc for x in ["E1", "E2"]):
        return "Elimination"
    elif any(x in desc for x in ["Diels-Alder", "Click", "CuAAC", "Cycloaddition"]):
        return "Cycloaddition"
    elif any(x in desc for x in ["Aldol", "Wittig", "Mannich", "Michael", "Claisen", "Knoevenagel", "Henry", "Reformatsky", "Horner"]):
        return "Carbonyl Chemistry"
    elif any(x in desc for x in ["Grignard"]):
        return "Organometallic"
    elif any(x in desc for x in ["Metathesis"]):
        return "Metathesis"
    elif any(x in desc for x in ["Paal-Knorr", "Hantzsch", "Biginelli", "Fischer Indole", "Pictet-Spengler", 
                                   "Benzimidazole", "Benzoxazole", "Benzothiazole", "Imidazole", "Carbazole", "Quinoline"]):
        return "Heterocycle Synthesis"
    elif any(x in desc for x in ["Boc", "Cbz", "FMOC", "TBS", "PMB", "Acetonide", "Trityl", "Protection", "Deprotection"]):
        return "Protecting Groups"
    else:
        return "Other"


def test_reaction(reaction_smiles: str, k: int = 50) -> dict:
    """
    Test a single reaction with cross-family recommendation.
    
    Returns dict with:
        - success: bool
        - recommendations: list of dicts
        - error: str (if any)
        - execution_time: float
    """
    start_time = time.time()
    
    try:
        result = chem.recommend.conditions(
            reaction=reaction_smiles,
            k=k,
            search_all_families=True,
            relax=None
        )
        
        execution_time = time.time() - start_time
        
        recommendations = result.get("recommendations", [])
        
        return {
            "success": True,
            "recommendations": recommendations,
            "num_recommendations": len(recommendations),
            "reaction_family": result.get("reaction_family", "Unknown"),
            "error": None,
            "execution_time": execution_time
        }
    
    except Exception as e:
        execution_time = time.time() - start_time
        return {
            "success": False,
            "recommendations": [],
            "num_recommendations": 0,
            "reaction_family": None,
            "error": str(e),
            "execution_time": execution_time
        }


def main():
    """Main test execution."""
    print("=" * 80)
    print("CROSS-FAMILY RECOMMENDATION TEST: ALL SAMPLE REACTIONS")
    print("=" * 80)
    print()
    
    # Filter out the placeholder first entry
    reactions_to_test = [r for r in SAMPLE_REACTIONS if not r.startswith("Select a sample")]
    
    print(f"Total reactions to test: {len(reactions_to_test)}")
    print()
    
    # Statistics tracking
    results_by_type = defaultdict(lambda: {
        "total": 0,
        "success": 0,
        "failed": 0,
        "with_recommendations": 0,
        "without_recommendations": 0,
        "total_recommendations": 0,
        "total_time": 0.0
    })
    
    overall_stats = {
        "total": 0,
        "success": 0,
        "failed": 0,
        "with_recommendations": 0,
        "without_recommendations": 0,
        "total_recommendations": 0,
        "total_time": 0.0
    }
    
    detailed_results = []
    
    # Test each reaction
    for i, reaction_str in enumerate(reactions_to_test, 1):
        # Extract SMILES and description
        if ">>" not in reaction_str:
            continue
        
        # Handle both standard (reactants>>products) and non-standard (reactants>>[reagent]>>products) formats
        parts = reaction_str.split(">>")
        if len(parts) == 2:
            # Standard format: reactants>>products (description)
            smiles = parts[0] + ">>" + parts[1].split(" (")[0]
        elif len(parts) >= 3:
            # Non-standard format: reactants>>[reagent]>>products (description)
            # Convert to standard: reactants.reagent>>products
            reactants = parts[0]
            reagent = parts[1]
            products = parts[2].split(" (")[0]
            smiles = f"{reactants}.{reagent}>>{products}"
        else:
            continue
        
        description = reaction_str.split("(")[-1].rstrip(")") if "(" in reaction_str else "No description"
        reaction_type = extract_reaction_type(reaction_str)
        
        print(f"[{i}/{len(reactions_to_test)}] Testing: {description[:60]}...")
        
        # Test reaction
        result = test_reaction(smiles)
        
        # Update statistics
        overall_stats["total"] += 1
        results_by_type[reaction_type]["total"] += 1
        
        if result["success"]:
            overall_stats["success"] += 1
            results_by_type[reaction_type]["success"] += 1
            
            if result["num_recommendations"] > 0:
                overall_stats["with_recommendations"] += 1
                overall_stats["total_recommendations"] += result["num_recommendations"]
                results_by_type[reaction_type]["with_recommendations"] += 1
                results_by_type[reaction_type]["total_recommendations"] += result["num_recommendations"]
                print(f"  ✓ Success: {result['num_recommendations']} recommendations ({result['execution_time']:.2f}s)")
            else:
                overall_stats["without_recommendations"] += 1
                results_by_type[reaction_type]["without_recommendations"] += 1
                print(f"  ⚠ Success but no recommendations ({result['execution_time']:.2f}s)")
        else:
            overall_stats["failed"] += 1
            results_by_type[reaction_type]["failed"] += 1
            print(f"  ✗ Failed: {result['error'][:80]}")
        
        overall_stats["total_time"] += result["execution_time"]
        results_by_type[reaction_type]["total_time"] += result["execution_time"]
        
        # Store detailed result
        detailed_results.append({
            "index": i,
            "smiles": smiles,
            "description": description,
            "reaction_type": reaction_type,
            "result": result
        })
    
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    
    # Overall statistics
    print("Overall Statistics:")
    print(f"  Total reactions tested:        {overall_stats['total']}")
    print(f"  Successful:                    {overall_stats['success']} ({overall_stats['success']/overall_stats['total']*100:.1f}%)")
    print(f"  Failed:                        {overall_stats['failed']} ({overall_stats['failed']/overall_stats['total']*100:.1f}%)")
    print(f"  With recommendations:          {overall_stats['with_recommendations']} ({overall_stats['with_recommendations']/overall_stats['total']*100:.1f}%)")
    print(f"  Without recommendations:       {overall_stats['without_recommendations']} ({overall_stats['without_recommendations']/overall_stats['total']*100:.1f}%)")
    print(f"  Total recommendations found:   {overall_stats['total_recommendations']}")
    if overall_stats['with_recommendations'] > 0:
        print(f"  Average recommendations/rxn:   {overall_stats['total_recommendations']/overall_stats['with_recommendations']:.1f}")
    print(f"  Total execution time:          {overall_stats['total_time']:.2f}s")
    print(f"  Average time per reaction:     {overall_stats['total_time']/overall_stats['total']:.2f}s")
    print()
    
    # Statistics by reaction type
    print("Statistics by Reaction Type:")
    print()
    
    # Sort by total count
    sorted_types = sorted(results_by_type.items(), key=lambda x: x[1]["total"], reverse=True)
    
    for reaction_type, stats in sorted_types:
        print(f"  {reaction_type}:")
        print(f"    Total:                {stats['total']}")
        print(f"    Success:              {stats['success']} ({stats['success']/stats['total']*100:.1f}%)")
        print(f"    With recommendations: {stats['with_recommendations']} ({stats['with_recommendations']/stats['total']*100:.1f}%)")
        if stats['with_recommendations'] > 0:
            print(f"    Avg recommendations:  {stats['total_recommendations']/stats['with_recommendations']:.1f}")
        print(f"    Avg time:             {stats['total_time']/stats['total']:.2f}s")
        print()
    
    # Save detailed results to JSON
    output_file = "test_all_sample_reactions_results.json"
    with open(output_file, "w") as f:
        json.dump({
            "overall_stats": overall_stats,
            "stats_by_type": dict(results_by_type),
            "detailed_results": detailed_results
        }, f, indent=2)
    
    print(f"Detailed results saved to: {output_file}")
    print()
    
    # Show top failures if any
    failures = [r for r in detailed_results if not r["result"]["success"]]
    if failures:
        print("=" * 80)
        print(f"TOP FAILURES ({len(failures)} total):")
        print("=" * 80)
        for i, failure in enumerate(failures[:10], 1):
            print(f"{i}. {failure['description'][:60]}")
            print(f"   Error: {failure['result']['error'][:80]}")
            print()
    
    # Show reactions without recommendations
    no_recs = [r for r in detailed_results if r["result"]["success"] and r["result"]["num_recommendations"] == 0]
    if no_recs:
        print("=" * 80)
        print(f"REACTIONS WITHOUT RECOMMENDATIONS ({len(no_recs)} total):")
        print("=" * 80)
        for i, no_rec in enumerate(no_recs[:10], 1):
            print(f"{i}. {no_rec['description'][:60]}")
            print(f"   Reaction family: {no_rec['result']['reaction_family']}")
            print()


if __name__ == "__main__":
    main()
