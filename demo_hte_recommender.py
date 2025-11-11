"""
Demo script for HTE-based condition recommendation system
"""
from chemtools.HTE import HTERecommender, format_result

# Initialize recommender
print("Initializing HTE Recommender...")
recommender = HTERecommender()

# Show database statistics
stats = recommender.get_statistics()
print("\n📊 DATABASE STATISTICS:")
for key, value in stats.items():
    print(f"  {key}: {value}")

print("\n" + "="*80)
print("TESTING HTE RECOMMENDATIONS")
print("="*80)

# Test Case 1: C-N Coupling (ArBr + RNH2)
print("\n\n🧪 TEST CASE 1: C-N Coupling")
print("Reactants: Bromobenzene + Ethylamine")
result1 = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",  # Bromobenzene (ArBr)
    reactant_b_smiles="CCN",  # Ethylamine (RNH2)
    top_k=5
)
print(format_result(result1))

# Test Case 2: Suzuki Coupling (ArBr + ArB(OH)2)
print("\n\n🧪 TEST CASE 2: Suzuki Coupling")
print("Reactants: 4-Bromobenzonitrile + Phenylboronic acid")
result2 = recommender.recommend(
    reactant_a_smiles="c1cc(Br)ccc1C#N",  # 4-Bromobenzonitrile
    reactant_b_smiles="c1ccc(B(O)O)cc1",  # Phenylboronic acid
    top_k=5
)
print(format_result(result2))

# Test Case 3: C-N Coupling (ArCl + ArNH2)
print("\n\n🧪 TEST CASE 3: C-N Coupling with Aniline")
print("Reactants: Chlorobenzene + Aniline")
result3 = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",  # Chlorobenzene
    reactant_b_smiles="c1ccc(N)cc1",  # Aniline
    top_k=5
)
print(format_result(result3))

# Test Case 4: Single reactant (CH-Activation)
print("\n\n🧪 TEST CASE 4: Single Reactant (ArBr)")
print("Reactant: Bromobenzene")
result4 = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",  # Bromobenzene
    reactant_b_smiles=None,
    top_k=5
)
print(format_result(result4))

# Test Case 5: Alkyl halide (different reactant type)
print("\n\n🧪 TEST CASE 5: Alkyl Halide")
print("Reactants: Ethyl bromide + Aniline")
result5 = recommender.recommend(
    reactant_a_smiles="CCBr",  # Ethyl bromide
    reactant_b_smiles="c1ccc(N)cc1",  # Aniline
    top_k=3
)
print(format_result(result5))

# Test Case 6: Filter by reaction type
print("\n\n🧪 TEST CASE 6: Suzuki with ArCl")
print("Reactants: Chlorobenzene + Phenylboronic acid (Suzuki only)")
result6 = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",  # Chlorobenzene
    reactant_b_smiles="c1ccc(B(O)O)cc1",  # Phenylboronic acid
    top_k=5,
    reaction_type_filter="Suzuki"
)
print(format_result(result6))

print("\n" + "="*80)
print("DEMO COMPLETE")
print("="*80)
