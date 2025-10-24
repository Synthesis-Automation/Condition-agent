"""Extract unique reaction types from sample_reactions.py"""

import re
from pathlib import Path

# Read the sample_reactions.py file
sample_file = Path(__file__).parent.parent / "tests" / "sample_reactions.py"
with open(sample_file, 'r', encoding='utf-8') as f:
    content = f.read()

# Extract reaction types from comments in parentheses
# Pattern: text in parentheses at end of SMILES string descriptions
pattern = r'\(([^)]+(?:Coupling|Reaction|Protection|Deprotection|Formation|Substitution|Chemistry|Cycloaddition|Arylation|Amination))\)'
matches = re.findall(pattern, content, re.IGNORECASE)

# Also extract from section headers
section_pattern = r'# ([A-Z\-\s]+REACTIONS?)'
section_matches = re.findall(section_pattern, content)

# Clean and deduplicate
reaction_types = set()

# Process inline matches
for match in matches:
    # Clean up the match
    cleaned = match.strip()
    # Remove specific descriptors to get base reaction type
    cleaned = re.sub(r'\s*-\s*.*$', '', cleaned)  # Remove everything after dash
    cleaned = re.sub(r'\s*\(.*?\)', '', cleaned)  # Remove nested parentheses
    reaction_types.add(cleaned)

# Process section headers
for match in section_matches:
    cleaned = match.strip()
    # Remove "REACTIONS" suffix
    cleaned = re.sub(r'\s+REACTIONS?$', '', cleaned)
    reaction_types.add(cleaned)

print("=" * 70)
print("REACTION TYPES FOUND IN sample_reactions.py")
print("=" * 70)

sorted_types = sorted(reaction_types)
for idx, rtype in enumerate(sorted_types, 1):
    print(f"{idx:2}. {rtype}")

print(f"\nTotal unique reaction types: {len(sorted_types)}")

# Find specific reaction mentions in comments
print("\n" + "=" * 70)
print("SPECIFIC REACTION MENTIONS")
print("=" * 70)

specific_reactions = {
    'Stille', 'Kumada', 'Chan-Lam', 'Click Chemistry', 
    'Ullmann', 'Boc Protection', 'Cbz Protection', 
    'TBS Protection', 'PMB Protection', 'Fmoc Protection',
    'Benzyl Protection', 'Acetyl Protection'
}

found_specific = set()
for reaction in specific_reactions:
    if reaction in content:
        found_specific.add(reaction)
        
for reaction in sorted(found_specific):
    print(f"  - {reaction}")

print(f"\nTotal: {len(found_specific)}")
