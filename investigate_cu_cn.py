#!/usr/bin/env python3
"""Investigate Cu catalyst data in HTE database"""
from chemtools.HTE import HTEAnalytics
import pandas as pd

a = HTEAnalytics()

print("=" * 80)
print("1. Check all Cu-catalyzed reactions:")
print("=" * 80)
df_cu_all = a.list_reactant_pairs(reaction_type=None, catalyst_filter='Cu', min_experiments=1)
print(f"Total Cu-catalyzed pairs: {len(df_cu_all)}")
print(f"Total Cu experiments: {df_cu_all['Num_Experiments'].sum() if len(df_cu_all) > 0 else 0}")

if len(df_cu_all) > 0:
    print("\nReaction type distribution:")
    print(df_cu_all.groupby('Reaction_Type')['Num_Experiments'].sum().sort_values(ascending=False).head(10))
    
    print("\nTop 10 Cu-catalyzed pairs:")
    print(df_cu_all.head(10)[['Reaction_Type', 'Reactant_A_Type', 'Reactant_B_Type', 'Num_Experiments', 'Success_Rate']])

print("\n" + "=" * 80)
print("2. Check C-N coupling reactions (all catalysts):")
print("=" * 80)
df_cn_all = a.list_reactant_pairs(reaction_type='C-N', catalyst_filter=None, min_experiments=1)
print(f"Total C-N coupling pairs: {len(df_cn_all)}")

if len(df_cn_all) > 0:
    print("\nTop 10 C-N coupling pairs:")
    print(df_cn_all.head(10)[['Reactant_A_Type', 'Reactant_B_Type', 'Num_Experiments', 'Top_Catalyst']])

print("\n" + "=" * 80)
print("3. Check catalyst stats for C-N coupling:")
print("=" * 80)
df_cat = a.get_catalyst_stats(reaction_type='C-N')
print(f"Total catalysts for C-N: {len(df_cat)}")

if len(df_cat) > 0:
    print("\nMetal distribution:")
    metal_dist = df_cat.groupby('Metal')['Num_Experiments'].sum().sort_values(ascending=False)
    print(metal_dist)
    
    print("\nTop Cu catalysts for C-N:")
    cu_cats = df_cat[df_cat['Metal'] == 'Cu']
    if len(cu_cats) > 0:
        print(cu_cats.head(10)[['Catalyst', 'Num_Experiments', 'Success_Rate']])
    else:
        print("No Cu catalysts found for C-N coupling")

print("\n" + "=" * 80)
print("4. Check reaction type naming:")
print("=" * 80)
rxn_summary = a.get_reaction_type_summary()
cn_related = rxn_summary[rxn_summary['Reaction_Type'].str.contains('N', case=False, na=False)]
print("Reaction types containing 'N':")
print(cn_related[['Reaction_Type', 'Num_Experiments', 'Top_Catalyst']])
