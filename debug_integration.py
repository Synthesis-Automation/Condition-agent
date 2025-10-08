"""
Debug script to check fusion evidence structure in core.py integration.
"""

from chemtools.recommend import recommend_from_reaction

reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

print("Running recommendation with fusion=True...")
results = recommend_from_reaction(
    reaction=reaction,
    k=50,
    use_fusion=True,
    max_variants=3
)

print(f"\n✅ Recommendation complete!")
print(f"\nTop-level keys: {results.keys()}")

if 'fusion_meta' in results:
    fusion_meta = results['fusion_meta']
    print(f"\nFusion meta keys: {fusion_meta.keys()}")
    
    evidence_summary = fusion_meta.get('evidence_summary', {})
    print(f"\nEvidence summary:")
    for key, value in evidence_summary.items():
        print(f"  {key}: {value} (type: {type(value).__name__})")
else:
    print("\n⚠️  No fusion_meta found!")
