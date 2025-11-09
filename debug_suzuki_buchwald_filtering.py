"""
Debug: Why are Suzuki and Buchwald-Hartwig both appearing?
===========================================================

Query: Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1

Expected:
- Has aryl chloride (Cl on indole)
- Has boronic acid (B(O)O on pyridine)
- Should match Suzuki (aryl halide + boronic acid) ✅
- Should NOT match Buchwald-Hartwig (needs aryl halide + amine, but no amine!) ❌

Let's check what's happening...
"""

from chemtools.recommend.unified import UnifiedRecommender
from chemtools.rule.analyzer import FeatureAnalyzer

reaction_smiles = "Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1"

print("=" * 80)
print("Debug: Feature Detection and Rule Filtering")
print("=" * 80)
print()

# Step 1: Detect features
print("Step 1: Feature Detection")
print("-" * 80)
analyzer = FeatureAnalyzer()
features = analyzer.analyze_reaction(reaction_smiles, combine_method="union")

print(f"Reaction: {reaction_smiles}")
print()
print("Detected Features:")
for feature, present in sorted(features.items()):
    if present:
        print(f"  ✓ {feature}")

print()
print()

# Step 2: Check Suzuki applies_if
print("Step 2: Check Suzuki applies_if")
print("-" * 80)

recommender = UnifiedRecommender("build/unified_index_complete")
suzuki_data = recommender.get_source_details("suzuki_v2")

if suzuki_data and 'applies_if' in suzuki_data:
    applies_if = suzuki_data['applies_if']
    print(f"Suzuki applies_if: {applies_if}")
    print()
    
    # Check ALL conditions
    if 'all' in applies_if:
        print("ALL conditions (all must be true):")
        for cond in applies_if['all']:
            present = features.get(cond, False)
            status = "✅" if present else "❌"
            print(f"  {status} {cond}: {present}")
    
    # Check ANY conditions
    if 'any' in applies_if:
        print("ANY conditions (at least one must be true):")
        for cond in applies_if['any']:
            present = features.get(cond, False)
            status = "✅" if present else "❌"
            print(f"  {status} {cond}: {present}")
    
    # Final check
    passes = recommender._check_applies_if(features, applies_if)
    print()
    print(f"Result: {'✅ PASSES' if passes else '❌ FAILS'}")

print()
print()

# Step 3: Check Buchwald-Hartwig applies_if
print("Step 3: Check Buchwald-Hartwig applies_if")
print("-" * 80)

bh_data = recommender.get_source_details("c_n_coupling_pd_v2")

if bh_data and 'applies_if' in bh_data:
    applies_if = bh_data['applies_if']
    print(f"Buchwald-Hartwig applies_if: {applies_if}")
    print()
    
    # Check ALL conditions
    if 'all' in applies_if:
        print("ALL conditions (all must be true):")
        for cond in applies_if['all']:
            present = features.get(cond, False)
            status = "✅" if present else "❌"
            print(f"  {status} {cond}: {present}")
    
    # Check ANY conditions
    if 'any' in applies_if:
        print("ANY conditions (at least one must be true):")
        for cond in applies_if['any']:
            present = features.get(cond, False)
            status = "✅" if present else "❌"
            print(f"  {status} {cond}: {present}")
    
    # Final check
    passes = recommender._check_applies_if(features, applies_if)
    print()
    print(f"Result: {'✅ PASSES' if passes else '❌ FAILS'}")

print()
print()

# Step 4: Check actual recommendation results
print("Step 4: Actual Recommendation Results")
print("-" * 80)

results = recommender.recommend(
    reaction_smiles,
    top_k=10,
    source_types=['rule'],
    validate_rules=True
)

print(f"Found {len(results)} rules after applies_if filtering:")
for r in results:
    print(f"  • {r.name} (similarity: {r.similarity:.3f})")

print()
print()

# Step 5: Analysis
print("=" * 80)
print("Analysis")
print("=" * 80)
print()

print("Query substrates:")
print("  1. Clc1cccc2c1cc[nH]2 - 7-chloroindole (aryl chloride)")
print("  2. c1ccc(B(O)O)nc1 - pyridine-3-boronic acid")
print()

print("Expected:")
print("  ✅ Suzuki should match (aryl halide + boronic acid)")
print("  ❌ Buchwald-Hartwig should NOT match (needs amine, but query has indole)")
print()

print("Question: Does indole count as an amine?")
print("  • Indole has N-H in the ring")
print("  • But it's NOT a nucleophilic amine for C-N coupling")
print("  • Feature detector might be incorrectly flagging 'amine_present'")
print()

if features.get('amine_present') or features.get('secondary_amine_present'):
    print("⚠️  WARNING: Amine detected! This may be the issue.")
    print("   Indole's N-H is being treated as an amine.")
    print("   Buchwald-Hartwig filter is passing because it sees 'amine'.")
else:
    print("✅ No amine detected - working correctly")
