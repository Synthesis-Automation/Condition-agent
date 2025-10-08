import sys
sys.path.insert(0, '.')

from chemtools.scdb_matcher.matcher import match
from chemtools.scdb_matcher.loader import load_db

# Load DB
db = load_db("data/conditionDB/suzuki_db.json")

# Test BULKY-NUC-XPHOS
rxn = "Brc1ccccc1.OB(O)c1c(C)cccc1C>>Cc1cccc(C)c1-c1ccccc1"
print(f"Testing: {rxn}\n")

result = match(db, rxn_smiles=rxn)

print(f"Matched: {result.entry_id} (priority {result.priority})")

# Extract features from the database entries to see what BULKY-NUC-XPHOS requires
bulky_entry = next((e for e in db.entries if e.id == "SCDB-SUZ-BULKY-NUC-XPHOS"), None)
if bulky_entry:
    print(f"\nBULKY-NUC-XPHOS requirements:")
    print(f"  feature_requirements: {bulky_entry.feature_requirements}")
    print(f"  priority: {bulky_entry.priority}")
