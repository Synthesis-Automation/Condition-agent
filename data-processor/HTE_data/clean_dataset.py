#!/usr/bin/env python3
"""
Clean and simplify the z-Score Peaks dataset by removing redundant and low-value fields.
Uses generic reactant naming to accommodate all reaction types (not all have electrophile/nucleophile).
Reduces from 28 columns to 16 essential columns.
"""

import pandas as pd
import os

# Define input and output file paths
INPUT_FILE = "z-Score Peaks with FG_STANDARDIZED.csv"
OUTPUT_FILE = "z-Score Peaks CLEANED.csv"

# Define the essential fields to keep and their new names (if renamed)
KEEP_COLUMNS = [
    # Reaction Classification
    "Reaction_Type_Standardized",
    
    # Performance Metrics
    "AREA_TOTAL_REDUCED",
    "z-Score",
    
    # Reactant Information (use FG A/B as source, more reliable than Reactant_Type fields)
    "FG A",  # Will rename to Reactant_A
    "FG B",  # Will rename to Reactant_B
    
    # Reactant Classification from standardized fields (backup/additional info)
    "Reactant_Type_Electrophile",      # Will rename to Reactant_A_Type
    "Reactant_Type_Nucleophile",       # Will rename to Reactant_B_Type
    "Reactant_Category_Electrophile",  # Will rename to Reactant_A_Category
    "Reactant_Category_Nucleophile",   # Will rename to Reactant_B_Category
    
    # Reaction Conditions
    "Catalyst",
    "Ligand",
    "Base",
    "Solvent",
    "Secondary Solvent",
    "Additive",
    "Coupling Reagent",
]

# Column renaming map (old_name -> new_name)
RENAME_COLUMNS = {
    "FG A": "Reactant_A",
    "FG B": "Reactant_B",
    "Reactant_Type_Electrophile": "Reactant_A_Type",
    "Reactant_Type_Nucleophile": "Reactant_B_Type",
    "Reactant_Category_Electrophile": "Reactant_A_Category",
    "Reactant_Category_Nucleophile": "Reactant_B_Category",
}

def main():
    print("=" * 80)
    print("Dataset Cleaning Script")
    print("=" * 80)
    
    # Check if input file exists
    if not os.path.exists(INPUT_FILE):
        print(f"❌ ERROR: Input file '{INPUT_FILE}' not found!")
        return
    
    # Load the dataset
    print(f"\n📂 Loading dataset: {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    
    print(f"   ✓ Loaded {len(df):,} rows")
    print(f"   ✓ Original columns: {len(df.columns)}")
    
    # Display original column names
    print(f"\n📋 Original columns ({len(df.columns)}):")
    for i, col in enumerate(df.columns, 1):
        print(f"   {i:2d}. {col}")
    
    # Check which columns to keep exist
    missing_columns = [col for col in KEEP_COLUMNS if col not in df.columns]
    if missing_columns:
        print(f"\n⚠️  WARNING: The following columns to keep are missing:")
        for col in missing_columns:
            print(f"   - {col}")
        print("\n   Proceeding with available columns only...")
        available_columns = [col for col in KEEP_COLUMNS if col in df.columns]
    else:
        available_columns = KEEP_COLUMNS
    
    # Create cleaned dataset
    print(f"\n🔧 Creating cleaned dataset with {len(available_columns)} columns...")
    df_cleaned = df[available_columns].copy()
    
    # Rename columns to generic names
    print(f"\n🏷️  Renaming columns to generic names...")
    df_cleaned = df_cleaned.rename(columns=RENAME_COLUMNS)
    
    # Display renamed columns
    for old_name, new_name in RENAME_COLUMNS.items():
        if old_name in available_columns:
            print(f"   {old_name:<40} → {new_name}")
    
    # Display statistics
    print(f"\n📊 Cleaned Dataset Statistics:")
    print(f"   Rows: {len(df_cleaned):,}")
    print(f"   Columns: {len(df_cleaned.columns)} (reduced from {len(df.columns)})")
    print(f"   Reduction: {(1 - len(df_cleaned.columns)/len(df.columns))*100:.1f}%")
    
    print(f"\n📋 Kept columns ({len(df_cleaned.columns)}):")
    for i, col in enumerate(df_cleaned.columns, 1):
        null_count = df_cleaned[col].isnull().sum()
        null_pct = (null_count / len(df_cleaned)) * 100
        unique_count = df_cleaned[col].nunique()
        print(f"   {i:2d}. {col:<40} | Nulls: {null_count:6,} ({null_pct:5.1f}%) | Unique: {unique_count:,}")
    
    # Display removed columns
    removed_columns = [col for col in df.columns if col not in available_columns]
    print(f"\n🗑️  Removed columns ({len(removed_columns)}):")
    for i, col in enumerate(removed_columns, 1):
        print(f"   {i:2d}. {col}")
    
    # Save cleaned dataset
    print(f"\n💾 Saving cleaned dataset: {OUTPUT_FILE}")
    df_cleaned.to_csv(OUTPUT_FILE, index=False)
    
    # File size comparison
    original_size = os.path.getsize(INPUT_FILE) / (1024 * 1024)  # MB
    cleaned_size = os.path.getsize(OUTPUT_FILE) / (1024 * 1024)  # MB
    size_reduction = (1 - cleaned_size / original_size) * 100
    
    print(f"\n📦 File Size Comparison:")
    print(f"   Original: {original_size:,.2f} MB")
    print(f"   Cleaned:  {cleaned_size:,.2f} MB")
    print(f"   Savings:  {original_size - cleaned_size:,.2f} MB ({size_reduction:.1f}% reduction)")
    
    print(f"\n✅ Dataset cleaning complete!")
    print(f"   Output file: {OUTPUT_FILE}")
    print("=" * 80)
    
    # Quick validation
    print(f"\n🔍 Quick Validation:")
    print(f"   Sample rows from cleaned dataset:")
    print(df_cleaned.head(3).to_string())
    
    print(f"\n✓ All done!")

if __name__ == "__main__":
    main()
