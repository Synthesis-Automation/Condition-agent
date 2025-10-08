"""Debug feature detection for specific test cases."""

import sys
sys.path.insert(0, '.')

from chemtools.scdb_matcher.matcher import match
from chemtools.scdb_matcher.loader import load_db

# Load database
db = load_db("data/conditionDB/suzuki_db.json")

# Test case: BULKY-NUC-XPHOS (PhBr + 2,6-dimethylphenylboronic acid)
test_smiles = "Brc1ccccc1.OB(O)c1c(C)cccc1C>>c1ccccc1-c1c(C)cccc1C"

print(f"Testing: BULKY-NUC-XPHOS")
print(f"SMILES: {test_smiles}")
print()

# Match
result = match(db, rxn_smiles=test_smiles)

print(f"Matched ID: {result.entry_id}")
print(f"Priority: {result.priority}")
print()

# Check features in trace
if hasattr(result, 'trace') and result.trace:
    trace = result.trace
    print("=== Feature Detection ===")
    if 'set_features' in trace:
        print("\nSet Features:")
        for key, value in trace['set_features'].items():
            print(f"  {key}: {value}")
    if 'numeric_features' in trace:
        print("\nNumeric Features:")
        for key, value in trace['numeric_features'].items():
            print(f"  {key}: {value}")
    print()
    
    # Check which entries were considered
    if 'applicable_entries' in trace:
        print(f"\n=== Applicable Entries ({len(trace['applicable_entries'])}) ===")
        for entry_id in trace['applicable_entries']:
            print(f"  - {entry_id}")
    
    if 'filtered_entries' in trace:
        print(f"\n=== Filtered Out ({len(trace['filtered_entries'])}) ===")
        for entry_id, reason in trace['filtered_entries'].items():
            print(f"  - {entry_id}: {reason}")
