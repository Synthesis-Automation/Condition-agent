import pandas as pd
from collections import defaultdict

# Read the CSV
df = pd.read_csv(r'c:\Git-softwares\Condition-agent\data-processor\other_data\z-Score Peaks with FG.csv')

# Define reagent columns (not reactant columns)
reagent_columns = {
    'Base': 'Base',
    'Catalyst': 'Catalyst',
    'Ligand': 'Ligand',
    'Additive': 'Additive',
    'Solvent': 'Solvent',
    'Secondary Solvent': 'Solvent',
    'Tertiary Solvent': 'Solvent',
    'Coupling Reagent': 'Coupling Reagent'
}

# Reactant columns to exclude (these are substrates/starting materials)
reactant_columns = ['Aryl Halide', 'N-Nucleophile/Boronate Type', 'FG A', 'FG B']

# Dictionary to store reagents by role
reagents_by_role = defaultdict(set)

# Extract all unique reagents
for column, role in reagent_columns.items():
    if column in df.columns:
        # Get all non-null values
        values = df[column].dropna()
        
        # Process each value
        for val in values:
            if val and str(val).strip() and str(val).lower() != 'nan':
                # Split by comma if multiple reagents listed
                reagents = [r.strip() for r in str(val).split(',')]
                for reagent in reagents:
                    if reagent and reagent.lower() not in ['missing', 'no ligand', '']:
                        reagents_by_role[role].add(reagent)

# Also get reaction types
reaction_types = set(df['Reaction Type'].dropna().unique())

# Sort reagents within each role
for role in reagents_by_role:
    reagents_by_role[role] = sorted(list(reagents_by_role[role]))

# Print summary
print("=" * 80)
print("REAGENT EXTRACTION SUMMARY")
print("=" * 80)
print(f"\nTotal rows in dataset: {len(df)}")
print(f"Reaction types: {len(reaction_types)}")
print("\nReagent counts by role:")
for role, reagents in sorted(reagents_by_role.items()):
    print(f"  {role}: {len(reagents)} unique reagents")

# Save to markdown file
output_file = 'extracted_reagents_registry.md'

with open(output_file, 'w', encoding='utf-8') as f:
    f.write("# Reagent Registry Extraction from z-Score Peaks Dataset\n\n")
    f.write("## Dataset Overview\n\n")
    f.write(f"- **Total Reactions:** {len(df):,}\n")
    f.write(f"- **Reaction Types:** {len(reaction_types)}\n")
    f.write(f"- **Unique Reagent Roles:** {len(reagents_by_role)}\n")
    f.write(f"- **Total Unique Reagents:** {sum(len(r) for r in reagents_by_role.values())}\n\n")
    
    f.write("---\n\n")
    
    # Write reaction types
    f.write("## Reaction Types Covered\n\n")
    reaction_type_counts = df['Reaction Type'].value_counts()
    for rxn_type, count in reaction_type_counts.items():
        f.write(f"- **{rxn_type}** ({count:,} reactions)\n")
    f.write("\n---\n\n")
    
    # Write reagents by role
    for role in sorted(reagents_by_role.keys()):
        reagents = reagents_by_role[role]
        f.write(f"## {role} ({len(reagents)} reagents)\n\n")
        
        for reagent in reagents:
            # Count occurrences
            count = 0
            for col in reagent_columns:
                if reagent_columns[col] == role and col in df.columns:
                    count += df[col].astype(str).str.contains(reagent, regex=False, na=False).sum()
            
            f.write(f"### {reagent}\n")
            f.write(f"- **Role:** {role}\n")
            f.write(f"- **Occurrences:** {count}\n")
            
            # Find reaction types where this reagent is used
            rxn_types_for_reagent = set()
            for col in reagent_columns:
                if reagent_columns[col] == role and col in df.columns:
                    mask = df[col].astype(str).str.contains(reagent, regex=False, na=False)
                    rxn_types_for_reagent.update(df[mask]['Reaction Type'].dropna().unique())
            
            if rxn_types_for_reagent:
                f.write(f"- **Used in reactions:** {', '.join(sorted(rxn_types_for_reagent))}\n")
            
            f.write("\n")
        
        f.write("---\n\n")
    
    # Add footer
    f.write("## Notes\n\n")
    f.write("- Reactant columns (Aryl Halide, N-Nucleophile/Boronate Type, FG A, FG B) were excluded\n")
    f.write("- Multiple reagents in a single field were split and counted separately\n")
    f.write("- 'Missing' and 'No Ligand' entries were excluded\n")
    f.write("- This data can be used to enrich the reagent registry with real experimental usage patterns\n\n")
    f.write(f"**Data source:** `z-Score Peaks with FG.csv`  \n")
    f.write(f"**Total records:** {len(df):,}  \n")
    f.write(f"**Date extracted:** 2025-10-24\n")

print(f"\n✓ Reagent registry saved to: {output_file}")
print(f"\nTop 10 most common reagents by role:")
for role in sorted(reagents_by_role.keys())[:3]:  # Show first 3 roles
    print(f"\n{role}:")
    reagents = reagents_by_role[role]
    counts = []
    for reagent in reagents:
        count = 0
        for col in reagent_columns:
            if reagent_columns[col] == role and col in df.columns:
                count += df[col].astype(str).str.contains(reagent, regex=False, na=False).sum()
        counts.append((reagent, count))
    
    counts.sort(key=lambda x: x[1], reverse=True)
    for reagent, count in counts[:10]:
        print(f"  {reagent}: {count} times")
