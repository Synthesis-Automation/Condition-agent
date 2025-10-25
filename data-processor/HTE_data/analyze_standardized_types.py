"""Analyze the standardized types used in the processed CSV."""

import pandas as pd
from pathlib import Path

def analyze_types():
    csv_path = Path(__file__).parent / "z-Score Peaks with FG_STANDARDIZED.csv"
    df = pd.read_csv(csv_path)
    
    print("UNIQUE STANDARDIZED REACTANT TYPES")
    print("=" * 80)
    
    # Electrophile types
    print("\nElectrophile Types (from Aryl Halide):")
    e_types = df['Reactant_Type_Electrophile'].dropna()
    e_types = e_types[e_types != '']
    all_e = set()
    for val in e_types:
        all_e.update([v.strip() for v in str(val).split(',')])
    
    print(f"Total unique: {len(all_e)}")
    for t in sorted(all_e):
        count = sum(t in str(row) for row in e_types)
        print(f"  - {t} ({count:,} occurrences)")
    
    # Nucleophile types
    print("\nNucleophile Types:")
    n_types = df['Reactant_Type_Nucleophile'].dropna()
    n_types = n_types[n_types != '']
    all_n = set()
    for val in n_types:
        all_n.update([v.strip() for v in str(val).split(',')])
    
    print(f"Total unique: {len(all_n)}")
    for t in sorted(all_n):
        count = sum(t in str(row) for row in n_types)
        print(f"  - {t} ({count:,} occurrences)")
    
    # Categories
    print("\n" + "=" * 80)
    print("CATEGORY USAGE")
    print("=" * 80)
    
    print("\nElectrophile Categories:")
    e_cats = df['Reactant_Category_Electrophile'].dropna()
    e_cats = e_cats[e_cats != '']
    all_e_cats = set()
    for val in e_cats:
        all_e_cats.update([v.strip() for v in str(val).split(',')])
    for cat in sorted(all_e_cats):
        count = sum(cat in str(row) for row in e_cats)
        print(f"  - {cat} ({count:,} occurrences)")
    
    print("\nNucleophile Categories:")
    n_cats = df['Reactant_Category_Nucleophile'].dropna()
    n_cats = n_cats[n_cats != '']
    all_n_cats = set()
    for val in n_cats:
        all_n_cats.update([v.strip() for v in str(val).split(',')])
    for cat in sorted(all_n_cats):
        count = sum(cat in str(row) for row in n_cats)
        print(f"  - {cat} ({count:,} occurrences)")

if __name__ == "__main__":
    analyze_types()
