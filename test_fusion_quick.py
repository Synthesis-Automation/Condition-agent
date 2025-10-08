"""Quick test of fusion on a simple Suzuki reaction"""

from chemtools.recommend.core import recommend_from_reaction

# Simple Suzuki reaction from sample_reactions.py
suzuki_smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"

print("=" * 80)
print("TESTING FUSION ON SUZUKI REACTION")
print("=" * 80)
print(f"\nReaction: {suzuki_smiles}")
print("Description: Suzuki - Simple Ph-Ph (bromobenzene + phenylboronic acid)")

print("\n[1] BASELINE k-NN RECOMMENDATIONS")
print("-" * 80)
try:
    baseline = recommend_from_reaction(
        reaction=suzuki_smiles,
        use_fusion=False,
        k=5
    )
    
    if baseline and 'formatted' in baseline:
        recs = baseline['formatted'].get('recommended_conditions', [])
        print(f"✅ Got {len(recs)} baseline recommendations")
        if recs:
            print(f"\nTop recommendation:")
            top = recs[0]
            summary = top.get('summary', {})
            print(f"  Catalyst: {summary.get('core', 'N/A')}")
            print(f"  Base: {summary.get('base', {}).get('name', 'N/A')}")
            print(f"  Solvent: {summary.get('solvent', {}).get('name', 'N/A')}")
    else:
        print(f"❌ Unexpected baseline result structure")
        print(f"   Keys: {list(baseline.keys()) if baseline else 'None'}")
        
except Exception as e:
    print(f"❌ ERROR in baseline: {e}")
    import traceback
    traceback.print_exc()

print("\n\n[2] FUSION MULTI-SOURCE RECOMMENDATIONS")
print("-" * 80)
try:
    fusion = recommend_from_reaction(
        reaction=suzuki_smiles,
        use_fusion=True,
        k=5
    )
    
    if fusion and 'formatted' in fusion:
        recs = fusion['formatted'].get('recommended_conditions', [])
        print(f"✅ Got {len(recs)} fusion recommendations")
        
        if recs:
            print(f"\nTop recommendation:")
            top = recs[0]
            summary = top.get('summary', {})
            print(f"  Catalyst: {summary.get('core', 'N/A')}")
            print(f"  Base: {summary.get('base', {}).get('name', 'N/A')}")
            print(f"  Solvent: {summary.get('solvent', {}).get('name', 'N/A')}")
            
            # Component scores
            if 'component_scores' in top:
                cs = top['component_scores']
                print(f"\n  Component Scores:")
                print(f"    PS (Precedent): {cs.get('PS', 0):.3f}")
                print(f"    AS (Analytics): {cs.get('AS', 0):.3f}")
                print(f"    RS (Rules): {cs.get('RS', 0):.3f}")
                print(f"    MS (ML): {cs.get('MS', 0):.3f}")
        
        # Fusion metadata
        if 'fusion_meta' in fusion:
            meta = fusion['fusion_meta']
            weights = meta.get('adaptive_weights', {})
            print(f"\n  Adaptive Weights:")
            print(f"    α (precedent): {weights.get('alpha', weights.get('α', 0)):.3f}")
            print(f"    β (analytics): {weights.get('beta', weights.get('β', 0)):.3f}")
            print(f"    γ (rules): {weights.get('gamma', weights.get('γ', 0)):.3f}")
            print(f"    δ (ML): {weights.get('delta', weights.get('δ', 0)):.3f}")
            
            evidence = meta.get('evidence_summary', {})
            print(f"\n  Evidence Quality:")
            print(f"    Precedent count: {evidence.get('precedent_count', evidence.get('precedents', 0))}")
            print(f"    Diversity score: {evidence.get('diversity_score', evidence.get('diversity', 0)):.3f}")
    else:
        print(f"❌ Unexpected fusion result structure")
        print(f"   Keys: {list(fusion.keys()) if fusion else 'None'}")
        if fusion and 'formatted' in fusion:
            print(f"   Formatted keys: {list(fusion['formatted'].keys())}")
        
except Exception as e:
    print(f"❌ ERROR in fusion: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 80)
print("TEST COMPLETE")
print("=" * 80)
