#!/usr/bin/env python3
"""Debug script to trace full matching logic."""

import sys
sys.path.insert(0, ".")

from chemtools.scdb_matcher import loader, match

# Test reaction with 2-bromotoluene (1 ortho substituent)
rxn_smiles = "Brc1ccccc1C.c1ccc(B(O)O)cc1>>Cc1ccccc1-c1ccccc1"

# Load database and match
db = loader.load_db("data/conditionDB/suzuki_db.json")
result = match(db, rxn_smiles)

print("="*70)
print("MATCH RESULT")
print("="*70)
print(f"Entry ID: {result.entry_id}")
print(f"Match Type: {result.match_type}")
print(f"Priority: {result.priority}")

print("\n" + "="*70)
print("TRACE - ALL ENTRIES")
print("="*70)

if result.trace and 'evaluations' in result.trace:
    for eval_data in result.trace['evaluations']:
        entry_id = eval_data.get('id', 'UNKNOWN')
        matched = eval_data.get('matched', False)
        filtered = eval_data.get('filtered', False)
        entry_type = eval_data.get('type', 'unknown')
        priority = eval_data.get('priority', 0)
        
        status = "âœ?MATCHED" if matched else ("âŠ?FILTERED" if filtered else "âœ?NO MATCH")
        
        print(f"\n{entry_id}")
        print(f"  Type: {entry_type}, Priority: {priority}")
        print(f"  Status: {status}")
        
        if entry_id == "SCDB-SUZ-ARBR-ORTHO-XPhos":
            print(f"  >>> THIS IS THE XPHOS ENTRY <<<")
            print(f"  Full evaluation: {eval_data}")

print("\n" + "="*70)
