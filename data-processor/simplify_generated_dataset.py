"""
Simplify generated dataset from process_reactions.py to only include columns needed for recommendation.

This script takes the full CSV output from process_reactions.py and creates a simplified version
with only the columns needed by simple_condition_recommender.py.

Usage:
    python simplify_generated_dataset.py --input reactions.csv --output reactions_simplified.csv
"""

import argparse
import pandas as pd
import json
from typing import List


def parse_json_list(value: str) -> List[str]:
    """Parse a JSON list string, return empty list if parsing fails."""
    if pd.isna(value) or value == '':
        return []
    try:
        result = json.loads(value)
        return result if isinstance(result, list) else []
    except:
        return []


def simplify_dataset(input_csv: str, output_csv: str):
    """
    Simplify dataset to only columns needed for recommendation.
    
    Columns to keep (16 total):
    1. Reaction_Type_Standardized (from ReactionType)
    2. AREA_TOTAL_REDUCED (placeholder - not in source, will be empty)
    3. z-Score (placeholder - not in source, will be empty)
    4. Reactant_A (from ReactantTypes[0] if available)
    5. Reactant_B (from ReactantTypes[1] if available)
    6. Reactant_A_Type (from ReactantTypes[0])
    7. Reactant_B_Type (from ReactantTypes[1])
    8. Reactant_A_Category (from ReactantCategories[0])
    9. Reactant_B_Category (from ReactantCategories[1])
    10. Catalyst (from CatalystCoreGeneric)
    11. Ligand
    12. Base (from Reagent where ReagentRole='BASE')
    13. Solvent
    14. Secondary Solvent (empty - not in source)
    15. Additive (from Reagent where ReagentRole='ADDITIVE')
    16. Coupling Reagent (from Reagent where ReagentRole contains 'COUPLING')
    """
    
    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)
    print(f"Input: {len(df)} rows, {len(df.columns)} columns")
    
    # Create simplified dataframe
    simplified = pd.DataFrame()
    
    # 1. Reaction type (direct copy)
    simplified['Reaction_Type_Standardized'] = df['ReactionType']
    
    # 2-3. Performance metrics (empty - these come from experimental results)
    simplified['AREA_TOTAL_REDUCED'] = ''
    simplified['z-Score'] = ''
    
    # 4-9. Reactant information from classification
    print("\nExtracting reactant types...")
    reactant_types_list = df['ReactantTypes'].apply(parse_json_list)
    reactant_categories_list = df['ReactantCategories'].apply(parse_json_list)
    
    # Use SMILES as reactant identifiers if types aren't available
    smiles_list = df['ReactantSMILES'].fillna('').str.split('.')
    
    simplified['Reactant_A'] = reactant_types_list.apply(lambda x: x[0] if len(x) > 0 else '')
    simplified['Reactant_B'] = reactant_types_list.apply(lambda x: x[1] if len(x) > 1 else '')
    simplified['Reactant_A_Type'] = reactant_types_list.apply(lambda x: x[0] if len(x) > 0 else '')
    simplified['Reactant_B_Type'] = reactant_types_list.apply(lambda x: x[1] if len(x) > 1 else '')
    simplified['Reactant_A_Category'] = reactant_categories_list.apply(lambda x: x[0] if len(x) > 0 else '')
    simplified['Reactant_B_Category'] = reactant_categories_list.apply(lambda x: x[1] if len(x) > 1 else '')
    
    # 10. Catalyst (from CatalystCoreGeneric - first item)
    print("Extracting catalyst...")
    catalyst_list = df['CatalystCoreGeneric'].apply(parse_json_list)
    simplified['Catalyst'] = catalyst_list.apply(lambda x: x[0] if len(x) > 0 else '')
    
    # 11. Ligand (from Ligand - first item)
    print("Extracting ligand...")
    ligand_list = df['Ligand'].apply(parse_json_list)
    simplified['Ligand'] = ligand_list.apply(lambda x: x[0] if len(x) > 0 else '')
    
    # 12-16. Reagents by role
    print("Extracting reagents by role...")
    reagent_list = df['Reagent'].apply(parse_json_list)
    role_list = df['ReagentRole'].apply(parse_json_list)
    
    def extract_by_role(reagents, roles, target_role):
        """Extract first reagent matching the target role."""
        for i, role in enumerate(roles):
            if target_role.upper() in role.upper():
                if i < len(reagents):
                    return reagents[i]
        return ''
    
    simplified['Base'] = [extract_by_role(r, ro, 'BASE') for r, ro in zip(reagent_list, role_list)]
    simplified['Solvent'] = df['Solvent'].apply(lambda x: parse_json_list(x)[0] if parse_json_list(x) else '')
    simplified['Secondary Solvent'] = [parse_json_list(s)[1] if len(parse_json_list(s)) > 1 else '' 
                                       for s in df['Solvent']]
    simplified['Additive'] = [extract_by_role(r, ro, 'ADDITIVE') for r, ro in zip(reagent_list, role_list)]
    simplified['Coupling Reagent'] = [extract_by_role(r, ro, 'COUPLING') for r, ro in zip(reagent_list, role_list)]
    
    # Save
    print(f"\nSaving to {output_csv}...")
    simplified.to_csv(output_csv, index=False)
    
    print(f"\n✓ Output: {len(simplified)} rows, {len(simplified.columns)} columns")
    print(f"\nColumn list:")
    for i, col in enumerate(simplified.columns, 1):
        non_empty = (simplified[col] != '').sum()
        print(f"  {i:2d}. {col:30s} ({non_empty:,} non-empty)")
    
    print(f"\n✓ Simplified dataset saved to: {output_csv}")
    print("\nNOTE: AREA_TOTAL_REDUCED and z-Score are empty - these must come from experimental data.")
    print("      This dataset shows the reaction conditions structure for the recommender.")


def main():
    parser = argparse.ArgumentParser(
        description="Simplify process_reactions.py output to recommender-ready format"
    )
    parser.add_argument('--input', '-i', required=True, help='Input CSV from process_reactions.py')
    parser.add_argument('--output', '-o', required=True, help='Output simplified CSV')
    
    args = parser.parse_args()
    
    simplify_dataset(args.input, args.output)


if __name__ == '__main__':
    main()
