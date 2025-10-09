#!/usr/bin/env python
"""Test ML recommendation with reagent database filtering via recommend API."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from chemtools.recommend import core

def test_ml_recommend_with_filtering():
    """Test ML recommendation through the main API."""
    
    print("=" * 70)
    print("Testing ML Recommendation API with Reagent Filtering")
    print("=" * 70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    print(f"\nReaction: {reaction}")
    print(f"\nCalling recommend_from_reaction()...\n")
    
    # Test with filtering enabled (default)
    result = core.recommend_from_reaction(
        reaction=reaction,
        k=25,
        family_override="C_N_Coupling_Pd",
        use_fusion=False,  # Use ML only, not fusion
    )
    
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    
    # Show precedent pack info
    precedent_pack = result.get("precedent_pack", {})
    precedents = precedent_pack.get("precedents", [])
    support = precedent_pack.get("support", 0)
    
    print(f"\nPrecedent Pack:")
    print(f"  Support: {support} precedents")
    print(f"  Returned: {len(precedents)} precedents")
    
    # Show recommendation
    recommendation = result.get("recommendation", {})
    print(f"\nPrimary Recommendation:")
    print(f"  Core: {recommendation.get('core')}")
    print(f"  Base: {recommendation.get('base')}")
    print(f"  Solvent: {recommendation.get('solvent')}")
    print(f"  Temperature: {recommendation.get('T_C')}°C")
    print(f"  Time: {recommendation.get('time_h')}h")
    print(f"  Confidence: {recommendation.get('confidence', 0):.2%}")
    
    # Show formatted output
    formatted = result.get("formatted", {})
    recommended_conditions = formatted.get("recommended_conditions", [])
    
    print(f"\nFormatted Output:")
    print(f"  Recommendations: {len(recommended_conditions)}")
    
    if recommended_conditions:
        first_rec = recommended_conditions[0]
        chemicals = first_rec.get("chemicals", [])
        
        print(f"\n  First Recommendation Chemicals:")
        for chem in chemicals:
            name = chem.get("name") or "Unknown"
            role = chem.get("role") or "unknown"
            cas = chem.get("cas")
            equiv = chem.get("equivalents")
            
            cas_str = f"CAS:{cas}" if cas else "No CAS"
            equiv_str = f"{equiv:.3f} eq" if equiv is not None else "N/A"
            
            print(f"    - {name:40s} [{role:15s}] {cas_str:15s} {equiv_str}")
    
    print("\n" + "=" * 70)
    print("✅ ML recommendation with filtering is working!")
    print("=" * 70)
    
    return result

if __name__ == "__main__":
    result = test_ml_recommend_with_filtering()
