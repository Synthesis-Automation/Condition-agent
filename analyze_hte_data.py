"""
Analyze HTE_0.csv dataset structure and statistics
"""
import pandas as pd
from collections import defaultdict

# Load data
df = pd.read_csv('data/HTE_db/HTE_0.csv')

print("="*80)
print("HTE DATABASE ANALYSIS")
print("="*80)

print(f"\n📊 DATASET OVERVIEW")
print(f"Total experiments: {len(df)}")
print(f"Columns: {len(df.columns)}")

print(f"\n📋 COLUMNS:")
for col in df.columns:
    print(f"  - {col}")

print(f"\n⚗️  REACTION TYPES:")
print(df['Reaction_Type_Standardized'].value_counts())

print(f"\n📈 YIELD STATISTICS:")
print(f"  Min: {df['AREA_TOTAL_REDUCED'].min():.2f}")
print(f"  Max: {df['AREA_TOTAL_REDUCED'].max():.2f}")
print(f"  Mean: {df['AREA_TOTAL_REDUCED'].mean():.2f}")
print(f"  Median: {df['AREA_TOTAL_REDUCED'].median():.2f}")

print(f"\n🔬 REACTANT A TYPES (Top 15):")
print(df['Reactant_A_Type'].value_counts().head(15))

print(f"\n🔬 REACTANT B TYPES (Top 15):")
print(df['Reactant_B_Type'].value_counts().head(15))

print(f"\n🏷️  REACTANT A CATEGORIES:")
print(df['Reactant_A_Category'].value_counts().head(15))

print(f"\n🏷️  REACTANT B CATEGORIES:")
print(df['Reactant_B_Category'].value_counts().head(15))

# Analyze condition components
print(f"\n⚙️  CATALYSTS (Top 15):")
print(df['Catalyst'].value_counts().head(15))

print(f"\n🧪 LIGANDS (Top 15):")
print(df['Ligand'].value_counts().head(15))

print(f"\n🧂 BASES (Top 15):")
print(df['Base'].value_counts().head(15))

print(f"\n🧴 SOLVENTS (Top 15):")
print(df['Solvent'].value_counts().head(15))

# Analyze combinations
print(f"\n🔗 REACTANT TYPE COMBINATIONS (Top 20):")
df['type_combo'] = df['Reactant_A_Type'].fillna('') + ' + ' + df['Reactant_B_Type'].fillna('')
print(df['type_combo'].value_counts().head(20))

# Success rate analysis
print(f"\n✅ SUCCESS ANALYSIS (Yield > 50):")
df['success'] = df['AREA_TOTAL_REDUCED'] > 50
success_rate = df['success'].mean() * 100
print(f"Overall success rate: {success_rate:.1f}%")

print(f"\n✅ SUCCESS BY REACTION TYPE:")
for rxn_type in df['Reaction_Type_Standardized'].unique():
    if pd.notna(rxn_type):
        subset = df[df['Reaction_Type_Standardized'] == rxn_type]
        sr = (subset['AREA_TOTAL_REDUCED'] > 50).mean() * 100
        print(f"  {rxn_type}: {sr:.1f}% ({len(subset)} experiments)")

# Analyze reactant type patterns for each reaction
print(f"\n🎯 REACTANT TYPE PATTERNS BY REACTION:")
for rxn_type in ['C_N_Coupling', 'Suzuki', 'CH-Activation', 'CO-Coupling']:
    subset = df[df['Reaction_Type_Standardized'] == rxn_type]
    if len(subset) > 0:
        print(f"\n  {rxn_type} ({len(subset)} experiments):")
        type_combos = subset['type_combo'].value_counts().head(5)
        for combo, count in type_combos.items():
            avg_yield = subset[subset['type_combo'] == combo]['AREA_TOTAL_REDUCED'].mean()
            print(f"    {combo}: {count} experiments, avg yield: {avg_yield:.1f}")

print("\n" + "="*80)
