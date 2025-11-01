#!/usr/bin/env python
"""
Test Step 2: Condition Recommendations
Tests the recommend_from_reaction API with regenerated datasets.
"""
import sys
from pathlib import Path
import time

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chemtools.recommend import recommend_from_reaction


def test_condition_recommendations():
    """Test condition recommendations with regenerated datasets."""
    
    print("\n" + "="*70)
    print("  Test Step 2: Condition Recommendations")
    print("="*70)
    
    print("\n📋 Testing recommend_from_reaction API:")
    print("   - Basic recommendations")
    print("   - Family detection")
    print("   - Condition scoring")
    
    # Test 2a: C-N coupling (should detect Cu or Pd family)
    print(f"\n{'─'*70}")
    print(f"  Test 2a: C-N Coupling Recommendations")
    print(f"{'─'*70}")
    
    rxn = 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'
    print(f"   Query: {rxn}")
    
    print(f"\n   ⏱️  Running recommendation...")
    start_time = time.time()
    
    result = recommend_from_reaction(reaction=rxn, k=10)
    
    elapsed = time.time() - start_time
    print(f"   ⏱️  Completed in {elapsed:.3f} seconds")
    
    # Extract results
    family = result.get('family', 'N/A')
    confidence = result.get('confidence', 0)
    recommendation = result.get('recommendation', {})
    formatted = result.get('formatted', {})
    conditions = formatted.get('recommended_conditions', [])
    support = result.get('precedent_pack', {}).get('support', 0)
    
    print(f"\n   ✅ Results:")
    print(f"      - Detected family: {family}")
    print(f"      - Confidence: {confidence:.2f}")
    print(f"      - Conditions found: {len(conditions)}")
    print(f"      - Dataset support: {support} reactions")
    
    # Show top recommendations
    if conditions:
        print(f"\n   📊 Top 3 recommendations:")
        for i, cond in enumerate(conditions[:3], 1):
            # Extract nested structures
            cond_details = cond.get('conditions', {})
            chemicals = cond.get('chemicals', [])
            summary = cond.get('summary', '')
            
            print(f"\n   [{i}] Rank: {cond.get('rank', 'N/A')}")
            print(f"       → Summary: {summary if summary else 'N/A'}")
            
            # Extract catalyst system
            catalysts = [c for c in chemicals if c.get('role') in ['metal_catalyst', 'ligand']]
            if catalysts:
                cat_names = [c.get('name', 'Unknown') for c in catalysts]
                print(f"       → Catalysts: {', '.join(cat_names)}")
            
            # Extract base
            bases = [c for c in chemicals if c.get('role') == 'base']
            if bases:
                base_names = [c.get('name', c.get('cas', 'Unknown')) for c in bases]
                print(f"       → Base: {', '.join(base_names)}")
            
            # Extract solvent
            solvents = [c for c in chemicals if c.get('role') == 'solvent']
            if solvents:
                solv_names = [c.get('name', c.get('cas', 'Unknown')) for c in solvents]
                print(f"       → Solvent: {', '.join(solv_names)}")
            
            # Temperature/Time from conditions
            temp = cond_details.get('temperature')
            if temp:
                print(f"       → Temperature: {temp}")
            
            time_val = cond_details.get('time')
            if time_val:
                print(f"       → Time: {time_val}")
    else:
        print(f"\n   ⚠️  No conditions found")
        print(f"   💡 This may be expected if:")
        print(f"      - Dataset is very small")
        print(f"      - No matching reactions in dataset")
        print(f"      - Family detection failed")
    
    # Test 2b: With family override (force Pd dataset)
    print(f"\n{'─'*70}")
    print(f"  Test 2b: With Family Override (C_N_Coupling_Pd)")
    print(f"{'─'*70}")
    
    print(f"   Query: {rxn}")
    print(f"   Override: C_N_Coupling_Pd (Buchwald-Hartwig)")
    
    print(f"\n   ⏱️  Running recommendation...")
    start_time_pd = time.time()
    
    result_pd = recommend_from_reaction(
        reaction=rxn,
        k=10,
        family_override='C_N_Coupling_Pd'
    )
    
    elapsed_pd = time.time() - start_time_pd
    print(f"   ⏱️  Completed in {elapsed_pd:.3f} seconds")
    
    # Extract results
    family_pd = result_pd.get('family', 'N/A')
    formatted_pd = result_pd.get('formatted', {})
    conditions_pd = formatted_pd.get('recommended_conditions', [])
    support_pd = result_pd.get('precedent_pack', {}).get('support', 0)
    
    print(f"\n   ✅ Results:")
    print(f"      - Family (overridden): {family_pd}")
    print(f"      - Conditions found: {len(conditions_pd)}")
    print(f"      - Dataset support: {support_pd} reactions")
    
    # Show top Pd recommendations
    if conditions_pd:
        print(f"\n   📊 Top 3 Pd-catalyzed recommendations:")
        for i, cond in enumerate(conditions_pd[:3], 1):
            # Extract nested structures
            cond_details = cond.get('conditions', {})
            chemicals = cond.get('chemicals', [])
            summary = cond.get('summary', '')
            
            print(f"\n   [{i}] Rank: {cond.get('rank', 'N/A')}")
            print(f"       → Summary: {summary if summary else 'N/A'}")
            
            # Extract Pd catalyst system
            catalysts = [c for c in chemicals if c.get('role') in ['metal_catalyst', 'ligand']]
            if catalysts:
                cat_names = [c.get('name', 'Unknown') for c in catalysts]
                print(f"       → Pd System: {', '.join(cat_names)}")
            
            # Base
            bases = [c for c in chemicals if c.get('role') == 'base']
            if bases:
                base_names = [c.get('name', c.get('cas', 'Unknown')) for c in bases]
                print(f"       → Base: {', '.join(base_names)}")
            
            # Solvent
            solvents = [c for c in chemicals if c.get('role') == 'solvent']
            if solvents:
                solv_names = [c.get('name', c.get('cas', 'Unknown')) for c in solvents]
                print(f"       → Solvent: {', '.join(solv_names)}")
    
    # Test 2c: With DRFP similarity
    print(f"\n{'─'*70}")
    print(f"  Test 2c: With DRFP Similarity Search")
    print(f"{'─'*70}")
    
    print(f"   Query: {rxn}")
    print(f"   DRFP threshold: 0.3 (minimum similarity)")
    
    relax = {
        "use_drfp": True,
        "reaction_smiles": rxn,
        "drfp_threshold": 0.3,
        "drfp_weight": 0.6,
        "precompute_drfp": False,
        "selective_loading": True
    }
    
    print(f"\n   ⏱️  Running recommendation with DRFP...")
    start_time_drfp = time.time()
    
    result_drfp = recommend_from_reaction(
        reaction=rxn,
        k=10,
        relax=relax
    )
    
    elapsed_drfp = time.time() - start_time_drfp
    print(f"   ⏱️  Completed in {elapsed_drfp:.3f} seconds")
    
    # Extract results
    formatted_drfp = result_drfp.get('formatted', {})
    conditions_drfp = formatted_drfp.get('recommended_conditions', [])
    
    print(f"\n   ✅ Results:")
    print(f"      - Conditions found: {len(conditions_drfp)}")
    print(f"      - Using DRFP similarity for ranking")
    
    if conditions_drfp:
        print(f"\n   📊 Top recommendation (DRFP-ranked):")
        top = conditions_drfp[0]
        cond_details = top.get('conditions', {})
        chemicals = top.get('chemicals', [])
        summary = top.get('summary', '')
        
        print(f"       → Rank: {top.get('rank', 'N/A')}")
        print(f"       → Summary: {summary if summary else 'N/A'}")
        
        # Catalysts
        catalysts = [c for c in chemicals if c.get('role') in ['metal_catalyst', 'ligand']]
        if catalysts:
            cat_names = [c.get('name', 'Unknown') for c in catalysts]
            print(f"       → Catalysts: {', '.join(cat_names)}")
        
        # Base
        bases = [c for c in chemicals if c.get('role') == 'base']
        if bases:
            base_names = [c.get('name', c.get('cas', 'Unknown')) for c in bases]
            print(f"       → Base: {', '.join(base_names)}")
        
        # Solvent
        solvents = [c for c in chemicals if c.get('role') == 'solvent']
        if solvents:
            solv_names = [c.get('name', c.get('cas', 'Unknown')) for c in solvents]
            print(f"       → Solvent: {', '.join(solv_names)}")
    
    # Summary
    print(f"\n{'─'*70}")
    print(f"  Performance Summary")
    print(f"{'─'*70}")
    
    print(f"\n   ⏱️  Execution Times:")
    print(f"      - Basic recommendation: {elapsed:.3f}s")
    print(f"      - With family override: {elapsed_pd:.3f}s")
    print(f"      - With DRFP similarity: {elapsed_drfp:.3f}s")
    print(f"      - Average: {(elapsed + elapsed_pd + elapsed_drfp)/3:.3f}s")
    
    print(f"\n   📊 Results Summary:")
    print(f"      - Auto-detected family: {family}")
    print(f"      - Conditions (auto): {len(conditions)}")
    print(f"      - Conditions (Pd override): {len(conditions_pd)}")
    print(f"      - Conditions (DRFP): {len(conditions_drfp)}")
    
    # Validation
    print(f"\n   ✅ Validation:")
    
    family_detected = family and family != 'N/A'
    if family_detected:
        print(f"      ✅ Family auto-detection working")
    else:
        print(f"      ⚠️  Family not detected")
    
    override_worked = family_pd == 'C_N_Coupling_Pd'
    if override_worked:
        print(f"      ✅ Family override working")
    else:
        print(f"      ⚠️  Family override not applied correctly")
    
    fast_execution = all(t < 2.0 for t in [elapsed, elapsed_pd, elapsed_drfp])
    if fast_execution:
        print(f"      ✅ All executions completed in < 2 seconds")
    else:
        print(f"      ⚠️  Some executions slower than expected")
    
    has_results = len(conditions) > 0 or len(conditions_pd) > 0 or len(conditions_drfp) > 0
    if has_results:
        print(f"      ✅ Recommendations generated successfully")
    else:
        print(f"      ⚠️  No recommendations generated (may be expected for small datasets)")
    
    print("\n" + "="*70)
    print("  ✅ Test Step 2 Complete!")
    print("="*70 + "\n")
    
    # Return summary
    return {
        "family_detected": family,
        "basic_conditions": len(conditions),
        "pd_conditions": len(conditions_pd),
        "drfp_conditions": len(conditions_drfp),
        "avg_time": (elapsed + elapsed_pd + elapsed_drfp) / 3,
        "success": family_detected and override_worked and fast_execution
    }


if __name__ == "__main__":
    summary = test_condition_recommendations()
    
    # Exit with appropriate code
    if summary["success"]:
        print("🎉 All validation checks passed!\n")
        sys.exit(0)
    else:
        print("⚠️  Some validation checks failed. Review output above.\n")
        sys.exit(1)
