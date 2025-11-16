"""
Quick test to verify the benzylic halide classification fix is working
"""

from chemtools.HTE.recommender import HTERecommender

print("Testing Benzylic Halide Classification Fix\n")
print("=" * 70)

# Test molecules
test_cases = [
    ("Brc1ccccc1", "Phenyl bromide", "ArBr"),
    ("Brc1ccc(C)cc1", "4-Bromotoluene (USER'S QUERY)", "ArBr"),
    ("BrCc1ccccc1", "Benzyl bromide", "Bn-Br"),
]

recommender = HTERecommender()

print("\n1. Classification Test:")
print("-" * 70)
all_correct = True
for smiles, name, expected in test_cases:
    reactant_type, category = recommender._detect_reactant_types(smiles)
    status = "✅" if reactant_type == expected else "❌"
    if reactant_type != expected:
        all_correct = False
    print(f"{status} {name}")
    print(f"   SMILES: {smiles}")
    print(f"   Expected: {expected}, Got: {reactant_type}")

if all_correct:
    print("\n✅ ALL CLASSIFICATIONS CORRECT!")
else:
    print("\n❌ SOME CLASSIFICATIONS WRONG - Cache issue?")

# Test the actual query
print("\n" + "=" * 70)
print("2. User's Query Test (ArBr + ArNH2 + Cu + C-N Coupling):")
print("-" * 70)

result = recommender.recommend(
    reactant_a_smiles="Brc1ccc(C)cc1",  # 4-bromotoluene
    reactant_b_smiles="Nc1ccccc1",      # aniline
    catalyst_filter="Cu",
    reaction_type_filter="C_N_Coupling",
    top_k=3,
    min_experiments=1
)

print(f"\nReactant A: Brc1ccc(C)cc1 → {result.reactant_a_type}")
print(f"Reactant B: Nc1ccccc1 → {result.reactant_b_type}")
print(f"Predicted Reaction: {result.predicted_reaction_type}")
print(f"Matching Experiments: {result.total_matching_experiments}")

if result.total_matching_experiments > 0:
    print(f"\n✅ SUCCESS! Found {result.total_matching_experiments} experiments")
    print("\nTop 3 Conditions:")
    for i, rec in enumerate(result.recommendations[:3], 1):
        print(f"\n{i}. Z-Score: {rec.avg_z_score:.2f}")
        print(f"   {rec.catalyst} / {rec.ligand} / {rec.base} / {rec.solvent}")
        print(f"   Avg Yield: {rec.avg_yield:.1f}%")
else:
    print("\n❌ FAILED! No experiments found")
    print("\n⚠️  This means either:")
    print("   1. The GUI app has cached the OLD patterns")
    print("   2. Python bytecode cache (.pyc) needs clearing")
    print("   3. The fix didn't apply correctly")

print("\n" + "=" * 70)
print("INSTRUCTIONS:")
print("-" * 70)
print("If test shows ❌, try:")
print("  1. Kill the GUI app completely (Ctrl+C or Task Manager)")
print("  2. Delete Python cache: rm -r __pycache__ chemtools/__pycache__")
print("  3. Restart the GUI app")
print("=" * 70)
