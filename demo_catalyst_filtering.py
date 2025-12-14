"""
Demo: Catalyst Metal Type Filtering in HTE Recommender

This demo shows how to filter recommendations by catalyst metal type (Pd, Cu, Ni, etc.)
"""
from chemtools.HTE import HTERecommender, format_result

print("="*80)
print("HTE CATALYST FILTERING DEMO")
print("="*80)

# Initialize recommender
print("\n🚀 Initializing HTE Recommender...")
recommender = HTERecommender()

# Test Case 1: C-N Coupling - Compare Palladium vs Copper catalysts
print("\n" + "="*80)
print("TEST CASE 1: C-N Coupling (ArCl + ArNH2)")
print("Reactants: Chlorobenzene + Aniline")
print("="*80)

# All catalysts
print("\n📊 ALL CATALYSTS:")
result_all = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",
    reactant_b_smiles="c1ccc(N)cc1",
    top_k=3
)
print(f"  Matching experiments: {result_all.total_matching_experiments}")
if result_all.recommendations:
    for i, rec in enumerate(result_all.recommendations[:3], 1):
        print(f"  {i}. {rec.catalyst}")
        print(f"     Success: {rec.success_rate:.1f}% ({rec.num_experiments} exp, {rec.avg_yield:.1f}% avg)")

# Palladium only
print("\n🔵 PALLADIUM CATALYSTS ONLY:")
result_pd = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",
    reactant_b_smiles="c1ccc(N)cc1",
    top_k=3,
    catalyst_filter="Pd"
)
print(f"  Matching experiments: {result_pd.total_matching_experiments}")
if result_pd.recommendations:
    for i, rec in enumerate(result_pd.recommendations[:3], 1):
        print(f"  {i}. {rec.catalyst}")
        print(f"     Success: {rec.success_rate:.1f}% ({rec.num_experiments} exp, {rec.avg_yield:.1f}% avg)")
else:
    print("  No recommendations found")

# Copper only
print("\n🟠 COPPER CATALYSTS ONLY:")
result_cu = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",
    reactant_b_smiles="c1ccc(N)cc1",
    top_k=3,
    catalyst_filter="Cu"
)
print(f"  Matching experiments: {result_cu.total_matching_experiments}")
if result_cu.recommendations:
    for i, rec in enumerate(result_cu.recommendations[:3], 1):
        print(f"  {i}. {rec.catalyst}")
        print(f"     Success: {rec.success_rate:.1f}% ({rec.num_experiments} exp, {rec.avg_yield:.1f}% avg)")
else:
    print("  No recommendations found")


# Test Case 2: C-O Coupling - Should have copper catalysts
print("\n" + "="*80)
print("TEST CASE 2: C-O Coupling")
print("Reactants: Bromobenzene only (ArBr)")
print("="*80)

# Try to find C-O coupling conditions
result_co = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",
    reactant_b_smiles=None,
    top_k=5,
    reaction_type_filter="CO-Coupling"
)
print(f"\n📊 ALL CATALYSTS (CO-Coupling):")
print(f"  Matching experiments: {result_co.total_matching_experiments}")

if result_co.recommendations:
    for i, rec in enumerate(result_co.recommendations[:5], 1):
        print(f"  {i}. {rec.catalyst}")
        print(f"     Success: {rec.success_rate:.1f}% ({rec.num_experiments} exp, {rec.avg_yield:.1f}% avg)")

# Copper only for C-O coupling
result_co_cu = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",
    reactant_b_smiles=None,
    top_k=5,
    reaction_type_filter="CO-Coupling",
    catalyst_filter="copper"
)
print(f"\n🟠 COPPER CATALYSTS ONLY (CO-Coupling):")
print(f"  Matching experiments: {result_co_cu.total_matching_experiments}")

if result_co_cu.recommendations:
    for i, rec in enumerate(result_co_cu.recommendations[:5], 1):
        print(f"  {i}. {rec.catalyst}")
        print(f"     Success: {rec.success_rate:.1f}% ({rec.num_experiments} exp, {rec.avg_yield:.1f}% avg)")


# Test Case 3: Suzuki Coupling - Mostly palladium
print("\n" + "="*80)
print("TEST CASE 3: Suzuki Coupling")
print("Reactants: Chlorobenzene + Phenylboronic acid")
print("="*80)

result_suzuki_pd = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",
    reactant_b_smiles="c1ccc(B(O)O)cc1",
    top_k=3,
    reaction_type_filter="Suzuki",
    catalyst_filter="palladium"
)
print(f"\n🔵 PALLADIUM CATALYSTS (Suzuki):")
print(f"  Matching experiments: {result_suzuki_pd.total_matching_experiments}")

if result_suzuki_pd.recommendations:
    rec = result_suzuki_pd.recommendations[0]
    print(f"\n  TOP RECOMMENDATION:")
    print(f"  Catalyst: {rec.catalyst}")
    print(f"  Ligand: {rec.ligand}")
    print(f"  Base: {rec.base}")
    print(f"  Solvent: {rec.solvent}")
    print(f"  Success: {rec.success_rate:.1f}% ({rec.num_experiments} exp, {rec.avg_yield:.1f}% avg)")
    print(f"  Confidence: {rec.confidence_score:.1f}/100")


print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print("""
✅ Catalyst filtering is now supported!

Usage:
  Python API: 
    result = recommender.recommend(
        reactant_a_smiles="...",
        reactant_b_smiles="...",
        catalyst_filter="Pd"  # or "Cu", "Ni", "palladium", "copper", etc.
    )
  
  CLI:
    python -m chemtools.HTE.cli -a "..." -b "..." --catalyst Pd
    python -m chemtools.HTE.cli -a "..." -b "..." --catalyst copper

Supported metals:
  - Palladium (Pd, palladium) - 45,205 experiments (68.2%)
  - Copper (Cu, copper) - 5,907 experiments (8.9%)
  - Nickel (Ni, nickel) - 619 experiments (0.9%)
  - Plus: Ir, Rh, Ru, Pt, Au, Ag, Fe, Co, Zn
""")
print("="*80)


