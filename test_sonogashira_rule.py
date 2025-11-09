"""
Test script to verify Sonogashira rule database is accessible.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from chem_assistant.chemtools_wrapper import find_conditions_by_rule

# Test reaction from user's query
reaction_smiles = "Clc1ccc(C#N)cc1.C#CC(C)(C)C>>CC(C)(C)C#Cc1ccc(C#N)cc1"

print("Testing Sonogashira rule database access...")
print(f"Reaction: {reaction_smiles}")
print()

try:
    result = find_conditions_by_rule(
        reaction_smiles=reaction_smiles,
        database="sonogashira",
        auto_detect=False,
        explain=True
    )
    
    print("✅ SUCCESS! Found conditions:")
    print()
    print(f"Database used: {result.get('database_used', 'N/A')}")
    print(f"Matched rule: {result.get('matched_rule_id', 'N/A')}")
    print()
    print("Recommended conditions:")
    conditions = result.get('recommended_conditions', {})
    for key, value in conditions.items():
        print(f"  {key}: {value}")
    
    if 'explanation' in result:
        print()
        print("Explanation:")
        print(result['explanation'])
        
except Exception as e:
    print(f"❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
