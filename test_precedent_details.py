#!/usr/bin/env python
"""Test ML recommendation output with precedent details."""

import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from chemtools.recommend import core as recommend_core

def test_precedent_details():
    """Test that precedent details are included in ML output."""
    
    print("=" * 70)
    print("Testing ML Recommendation - Precedent Details in Output")
    print("=" * 70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    print(f"\nReaction: {reaction}")
    print(f"Family: Buchwald-Hartwig C-N Coupling\n")
    
    # Get recommendation with precedent details
    result = recommend_core.recommend_from_reaction(
        reaction=reaction,
        k=25,
        family_override="C_N_Coupling_Pd",
        use_fusion=False,
    )
    
    # Get formatted output
    formatted = result.get("formatted", {})
    
    # Save full output to JSON
    output_file = Path("ml_recommendation_with_precedents.json")
    with open(output_file, 'w') as f:
        json.dump(formatted, f, indent=2)
    print(f"✅ Full output saved to: {output_file}\n")
    
    # Check for precedents_used section
    precedents_used = formatted.get("precedents_used", {})
    
    if not precedents_used:
        print("❌ ERROR: No precedents_used section found in output!")
        return
    
    print("=" * 70)
    print("PRECEDENTS USED IN RECOMMENDATION")
    print("=" * 70)
    
    # Overall statistics
    total_count = precedents_used.get("total_count", 0)
    stats = precedents_used.get("statistics", {})
    
    print(f"\n📊 Overall Statistics:")
    print(f"  Total precedents: {total_count}")
    print(f"  Average yield: {stats.get('average_yield')}%")
    print(f"  Yield range: {stats.get('yield_range')}")
    print(f"  Temperature range: {stats.get('temperature_range')}°C")
    print(f"  Time range: {stats.get('time_range')}h")
    
    # Top precedents
    top_precs = precedents_used.get("top_precedents", [])
    print(f"\n📚 Top {len(top_precs)} Precedents (Ranked by Similarity):")
    
    for prec in top_precs[:5]:  # Show first 5
        print(f"\n  [{prec.get('rank')}] Reaction ID: {prec.get('reaction_id')}")
        print(f"      Core: {prec.get('core')}")
        print(f"      Yield: {prec.get('yield')}%")
        
        # Show catalysts
        catalysts = prec.get("catalysts", [])
        if catalysts:
            print(f"      Catalysts:")
            for cat in catalysts:
                print(f"        - {cat.get('name')} (CAS: {cat.get('cas')})")
        
        # Show reagents
        reagents = prec.get("reagents", [])
        if reagents:
            print(f"      Reagents:")
            for r in reagents:
                print(f"        - {r.get('name')} ({r.get('role')}, CAS: {r.get('cas')})")
        
        # Show solvents
        solvents = prec.get("solvents", [])
        if solvents:
            print(f"      Solvents:")
            for s in solvents:
                print(f"        - {s.get('name')} (CAS: {s.get('cas')})")
        
        # Show conditions
        conditions = prec.get("conditions", {})
        if conditions:
            temp = conditions.get("temperature_C")
            time = conditions.get("time_h")
            if temp or time:
                print(f"      Conditions: ", end="")
                if temp:
                    print(f"{temp}°C", end="")
                if time:
                    print(f", {time}h", end="")
                print()
    
    # Core-matched precedents
    core_matched = precedents_used.get("core_matched_precedents", {})
    core = core_matched.get("core")
    core_count = core_matched.get("count", 0)
    core_examples = core_matched.get("examples", [])
    
    print(f"\n🎯 Precedents Matching Chosen Core ({core}):")
    print(f"  Count: {core_count}")
    
    if core_examples:
        print(f"  Top {len(core_examples)} examples:")
        for ex in core_examples:
            print(f"    [{ex.get('rank')}] {ex.get('reaction_id')}: ", end="")
            print(f"{ex.get('base')} + {ex.get('solvent')} → {ex.get('yield')}% yield")
    
    print("\n" + "=" * 70)
    print("✅ Precedent details successfully included in output!")
    print("=" * 70)
    print(f"\nFull JSON output: {output_file}")

if __name__ == "__main__":
    test_precedent_details()
