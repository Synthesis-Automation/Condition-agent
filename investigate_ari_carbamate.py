#!/usr/bin/env python3
"""
Investigate ArI + Carbamate + Cu conditions
"""
from chemtools.HTE import HTERecommender
import pandas as pd

# Load database and check raw data
recommender = HTERecommender()

# Filter for C_N_Coupling + Cu + ArI + Carbamate
df = recommender.df[
    (recommender.df['Reaction_Type_Standardized'] == 'C_N_Coupling') &
    (recommender.df['Catalyst_Metal'] == 'Cu') &
    (recommender.df['Reactant_A_Type'] == 'ArI') &
    (recommender.df['Reactant_B_Type'] == 'Carbamate')
]

print(f"Total experiments: {len(df)}")
print(f"\nColumns: {list(df.columns[:20])}")

if len(df) > 0:
    # Group by condition and get statistics
    condition_cols = ['Catalyst', 'Ligand', 'Base', 'Solvent']
    available_cols = [col for col in condition_cols if col in df.columns]
    
    print(f"\nGrouping by: {available_cols}")
    
    grouped = df.groupby(available_cols).agg({
        'Product_Yield_PCT_Area_UV': ['count', 'mean', 'median', 'std']
    }).reset_index()
    
    grouped.columns = available_cols + ['Count', 'Avg_Yield', 'Median_Yield', 'Std_Yield']
    grouped['Success_Rate'] = df.groupby(available_cols)['Product_Yield_PCT_Area_UV'].apply(
        lambda x: (x > 50).sum() / len(x) * 100
    ).values
    
    # Sort by count
    grouped = grouped.sort_values('Count', ascending=False)
    
    print(f"\nTop 10 conditions:")
    print(grouped.head(10).to_string())
else:
    print("\nNo experiments found!")
    
    # Check what's actually there
    print("\nChecking reactant types for Cu C_N_Coupling:")
    cu_cn = recommender.df[
        (recommender.df['Reaction_Type_Standardized'] == 'C_N_Coupling') &
        (recommender.df['Catalyst_Metal'] == 'Cu')
    ]
    
    print(f"Total Cu C-N coupling: {len(cu_cn)}")
    print("\nReactant A types:")
    print(cu_cn['Reactant_A_Type'].value_counts().head(10))
    print("\nReactant B types:")
    print(cu_cn['Reactant_B_Type'].value_counts().head(10))
