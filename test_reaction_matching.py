#!/usr/bin/env python3
"""Test reaction type matching"""
from chemtools.HTE import HTEAnalytics

a = HTEAnalytics()

print("Test 1: 'C-N'")
df1 = a.list_reactant_pairs(reaction_type='C-N', catalyst_filter='Cu', min_experiments=1)
print(f"Results: {len(df1)} pairs")

print("\nTest 2: 'C_N'")
df2 = a.list_reactant_pairs(reaction_type='C_N', catalyst_filter='Cu', min_experiments=1)
print(f"Results: {len(df2)} pairs")

print("\nTest 3: 'C_N_Coupling'")
df3 = a.list_reactant_pairs(reaction_type='C_N_Coupling', catalyst_filter='Cu', min_experiments=1)
print(f"Results: {len(df3)} pairs")
if len(df3) > 0:
    print("\nTop 5 Cu-catalyzed C_N_Coupling pairs:")
    print(df3.head(5)[['Reactant_A_Type', 'Reactant_B_Type', 'Num_Experiments', 'Success_Rate', 'Top_Catalyst']])

print("\nTest 4: 'coupling' (general)")
df4 = a.list_reactant_pairs(reaction_type='coupling', catalyst_filter='Cu', min_experiments=1)
print(f"Results: {len(df4)} pairs")

print("\nTest 5: No reaction filter, just Cu")
df5 = a.list_reactant_pairs(reaction_type=None, catalyst_filter='Cu', min_experiments=1)
print(f"Results: {len(df5)} pairs")
print("\nReaction types:")
print(df5['Reaction_Type'].value_counts())
