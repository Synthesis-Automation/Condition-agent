"""
Direct test of Sonogashira rule database without heavy imports.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from chemtools.rule import RuleEngine

# Test reaction from user's query
reaction_smiles = "Clc1ccc(C#N)cc1.C#CC(C)(C)C>>CC(C)(C)C#Cc1ccc(C#N)cc1"
db_path = project_root / "data" / "rule_db_v2" / "sonogashira_v2.json"

print("Testing Sonogashira rule database...")
print(f"Database path: {db_path}")
print(f"File exists: {db_path.exists()}")
print()

if not db_path.exists():
    print(f"❌ ERROR: Database file not found at {db_path}")
    sys.exit(1)

try:
    # Load the rule engine
    print("Loading rule engine...")
    engine = RuleEngine.from_file(db_path)
    print(f"✅ Rule engine loaded successfully")
    print()
    
    # Test with the reaction
    print(f"Testing reaction: {reaction_smiles}")
    print()
    
    result = engine.recommend(reaction_smiles)
    
    print("✅ SUCCESS! Got recommendation")
    print()
    
    # Print the formatted summary
    print(result.format_summary())
        
except Exception as e:
    print(f"❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
