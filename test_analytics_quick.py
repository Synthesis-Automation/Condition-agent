"""
Quick test of dataset analytics module.
"""

from chemtools import chem

print("Testing Dataset Analytics Module...")
print()

# Test 1: Get all families
print("[1] Get all families...")
families = chem.analytics.get_all_families()
print(f"    Found {len(families)} families: {', '.join(families)}")
print()

# Test 2: Get stats
print("[2] Get basic stats for Suzuki...")
stats = chem.analytics.get_stats("Suzuki")
print(f"    Total reactions: {stats['total_reactions']:,}")
print(f"    Unique bases: {stats['unique_bases']}")
print(f"    Unique solvents: {stats['unique_solvents']}")
if stats['yield_stats']:
    print(f"    Avg yield: {stats['yield_stats']['mean']:.1f}%")
print()

# Test 3: Get common catalysts
print("[3] Get top 5 catalysts for C_N_Coupling_Pd...")
catalysts = chem.analytics.get_common_catalysts("C_N_Coupling_Pd", top_n=5)
for name, count, avg_yield in catalysts:
    yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
    print(f"    {count:>4} reactions | {yield_str:>6} | {name[:50]}")
print()

# Test 4: Get common bases
print("[4] Get top 5 bases for Suzuki...")
bases = chem.analytics.get_common_bases("Suzuki", top_n=5)
for name, count, avg_yield in bases:
    yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
    print(f"    {count:>4} reactions | {yield_str:>6} | {name[:50]}")
print()

# Test 5: Get plate recommendations
print("[5] Get 12 plate recommendations for C_N_Coupling_Pd...")
plate = chem.analytics.get_plate_recommendations("C_N_Coupling_Pd", n_conditions=12, min_yield=60.0)
print(f"    Generated {len(plate)} conditions")
if plate:
    cond = plate[0]
    print(f"    Top condition:")
    print(f"      Catalyst: {cond['catalyst']}")
    print(f"      Ligand: {cond['ligand']}")
    print(f"      Base: {cond['base']}")
    print(f"      Avg yield: {cond['avg_yield']:.1f}%")
    print(f"      Frequency: {cond['frequency']} precedents")
print()

print("[OK] All tests passed!")
