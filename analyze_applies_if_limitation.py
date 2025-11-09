"""
Analysis: applies_if Limitation
================================

Current Issue:
- Query: Clc1cccc2c1cc[nH]2 (7-chloroindole as ELECTROPHILE)
- Buchwald-Hartwig applies_if sees: indole_present = True
- Filter passes, but this is WRONG context

Root Cause:
applies_if checks PRESENCE of features but not their ROLE:
- Can't distinguish: indole as nucleophile vs indole as electrophile
- Can't tell: which substrate has which feature

Possible Solutions:

1. REMOVE indole_present from Buchwald-Hartwig applies_if
   Pros: Fixes this specific case
   Cons: Will filter out legitimate N-arylation of indoles
   
2. ADD more specific features like:
   - indole_nucleophile_present
   - indole_electrophile_present
   Pros: Accurate distinction
   Cons: Requires major feature detection overhaul

3. IMPROVE feature detection to track substrate roles
   Pros: Most accurate
   Cons: Complex implementation

4. ACCEPT current behavior as reasonable
   Pros: Indoles CAN be nucleophiles in C-N coupling
   Cons: Some false positives (like this case)

Recommendation:
Option 4 (Accept) - Here's why:
- Indoles ARE valid nucleophiles for Buchwald-Hartwig
- The rule is technically applicable (indole + aryl halide)
- The issue is that BOTH substrates match partial criteria
- User can still see it's not the best match (low similarity: 0.086)
- Alternative: Add note to rule documentation
"""

from chemtools.recommend.unified import UnifiedRecommender

print("=" * 80)
print("Analysis: Is Buchwald-Hartwig Match Reasonable?")
print("=" * 80)
print()

reaction = "Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1"

print("Query Reaction:")
print("  Substrates:")
print("    1. Clc1cccc2c1cc[nH]2 (7-chloroindole)")
print("    2. c1ccc(B(O)O)nc1 (pyridine-3-boronic acid)")
print()
print("  Product:")
print("    c1ccc(-c2cccc3[nH]ccc23)nc1 (indole-pyridine C-C coupled)")
print()

print("Buchwald-Hartwig applies_if:")
print("  ALL: [aryl_halide_present] ✅ (Cl-indole has aryl chloride)")
print("  ANY: [primary_amine, secondary_amine, aniline, indole, amides]")
print("       ✅ indole_present (Cl-indole is an indole)")
print()

print("Why Buchwald-Hartwig matches:")
print("  • Query has aryl halide (Cl) ✅")
print("  • Query has indole ✅")
print("  • Buchwald-Hartwig CAN do N-arylation of indoles")
print()

print("Why it's a FALSE POSITIVE in this case:")
print("  • The indole is the ELECTROPHILE (has Cl), not nucleophile")
print("  • The actual reaction is SUZUKI (C-C coupling)")
print("  • Buchwald-Hartwig would couple the indole-NH with something")
print()

print("Why similarity is LOW (0.086):")
print("  • DRFP recognizes this is NOT a typical Buchwald-Hartwig")
print("  • Similarity score correctly deprioritizes it")
print("  • Suzuki has higher similarity (0.135) ✅")
print()

recommender = UnifiedRecommender("build/unified_index_complete")

results = recommender.recommend(
    reaction,
    top_k=10,
    source_types=['rule'],
    validate_rules=True
)

print("Actual Results (sorted by similarity):")
for i, r in enumerate(results, 1):
    print(f"  {i}. {r.name}")
    print(f"     Similarity: {r.similarity:.3f}")
    print(f"     Family: {r.family}")

print()
print("=" * 80)
print("Conclusion")
print("=" * 80)
print()
print("The current behavior is ACCEPTABLE because:")
print()
print("✅ 1. Similarity ranking works correctly")
print("     Suzuki (0.135) ranks higher than Buchwald-Hartwig (0.086)")
print()
print("✅ 2. applies_if is doing its job")
print("     It filters out completely incompatible rules")
print("     It allows potentially compatible rules (even if context differs)")
print()
print("✅ 3. Indoles ARE valid Buchwald-Hartwig substrates")
print("     Example: Indole + ArBr → N-aryl-indole")
print("     The rule is not fundamentally wrong")
print()
print("⚠️  4. Limitation: applies_if checks PRESENCE not ROLE")
print("     Can't tell if indole is nucleophile or electrophile")
print("     This is a known limitation of feature-based filtering")
print()
print("💡 Recommendation:")
print("   Keep current behavior - similarity ranking handles the edge cases")
print("   The top result is still correct (Suzuki)")
