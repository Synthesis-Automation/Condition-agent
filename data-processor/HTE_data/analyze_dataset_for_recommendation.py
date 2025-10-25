"""Analyze the z-Score dataset to understand recommendation opportunities."""

import pandas as pd
import numpy as np
from collections import defaultdict

def analyze_dataset():
    df = pd.read_csv('data-processor/other_data/z-Score Peaks with FG_STANDARDIZED.csv')
    
    print("=" * 80)
    print("DATASET STRUCTURE ANALYSIS FOR CONDITION RECOMMENDATION")
    print("=" * 80)
    
    print(f"\nTotal reactions: {len(df):,}")
    print(f"Columns: {len(df.columns)}")
    
    print("\nKey columns for recommendation:")
    print("  INPUT features:")
    print("    - Reactant_Type_Electrophile (e.g., ArBr, ArCl)")
    print("    - Reactant_Type_Nucleophile (e.g., RNH2, ArB(OR)2)")
    print("    - Reaction_Type_Standardized (e.g., Buchwald-Hartwig, Suzuki-Miyaura)")
    print("    - FG A, FG B (functional groups)")
    print("  OUTPUT conditions:")
    print("    - Catalyst")
    print("    - Ligand")
    print("    - Base")
    print("    - Solvent (Primary, Secondary, Tertiary)")
    print("    - Additive")
    print("    - Coupling Reagent")
    print("  PERFORMANCE metric:")
    print("    - z-Score (higher = better)")
    print("    - AREA_TOTAL_REDUCED (yield proxy)")
    
    print("\n" + "=" * 80)
    print("z-SCORE DISTRIBUTION")
    print("=" * 80)
    print(df['z-Score'].describe())
    
    # Analyze high performers
    high = df[df['z-Score'] > 1.0]
    medium = df[(df['z-Score'] > 0) & (df['z-Score'] <= 1.0)]
    low = df[df['z-Score'] <= 0]
    
    print(f"\nPerformance tiers:")
    print(f"  High (z > 1.0): {len(high):,} ({len(high)/len(df)*100:.1f}%)")
    print(f"  Medium (0 < z <= 1.0): {len(medium):,} ({len(medium)/len(df)*100:.1f}%)")
    print(f"  Low (z <= 0): {len(low):,} ({len(low)/len(df)*100:.1f}%)")
    
    print("\n" + "=" * 80)
    print("TOP REACTION TYPES (by high z-score count)")
    print("=" * 80)
    for rxn, count in high['Reaction_Type_Standardized'].value_counts().head(10).items():
        avg_z = high[high['Reaction_Type_Standardized']==rxn]['z-Score'].mean()
        total_rxn = len(df[df['Reaction_Type_Standardized']==rxn])
        success_rate = count / total_rxn * 100
        print(f"  {rxn}:")
        print(f"    High performers: {count:,} / {total_rxn:,} ({success_rate:.1f}%)")
        print(f"    Avg z-Score (high): {avg_z:.2f}")
    
    print("\n" + "=" * 80)
    print("CONDITION COVERAGE ANALYSIS")
    print("=" * 80)
    
    # Analyze condition combinations
    print(f"\nUnique catalysts: {df['Catalyst'].nunique()}")
    print(f"Unique ligands: {df['Ligand'].nunique()}")
    print(f"Unique bases: {df['Base'].nunique()}")
    print(f"Unique solvents: {df['Solvent'].nunique()}")
    
    # Most common successful conditions
    print("\n" + "=" * 80)
    print("MOST SUCCESSFUL CONDITIONS (z > 1.0)")
    print("=" * 80)
    
    print("\nTop Catalysts:")
    for cat, count in high['Catalyst'].value_counts().head(10).items():
        avg_z = high[high['Catalyst']==cat]['z-Score'].mean()
        print(f"  {cat}: {count:,} uses (avg z: {avg_z:.2f})")
    
    print("\nTop Ligands:")
    for lig, count in high['Ligand'].value_counts().head(10).items():
        avg_z = high[high['Ligand']==lig]['z-Score'].mean()
        print(f"  {lig}: {count:,} uses (avg z: {avg_z:.2f})")
    
    print("\nTop Bases:")
    for base, count in high['Base'].value_counts().head(10).items():
        avg_z = high[high['Base']==base]['z-Score'].mean()
        print(f"  {base}: {count:,} uses (avg z: {avg_z:.2f})")
    
    print("\nTop Solvents:")
    for solv, count in high['Solvent'].value_counts().head(10).items():
        avg_z = high[high['Solvent']==solv]['z-Score'].mean()
        print(f"  {solv}: {count:,} uses (avg z: {avg_z:.2f})")

if __name__ == "__main__":
    analyze_dataset()
