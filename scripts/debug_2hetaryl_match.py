#!/usr/bin/env python3
"""Debug why 2-bromopyridine matches ortho rule."""

import sys
sys.path.insert(0, ".")

from chemtools.scdb_matcher import loader, match

db = loader.load_db("data/conditionDB/suzuki_db.json")

# Test 2-bromopyridine
rxn = "Brc1ccccn1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccn2)cc1"
result = match(db, rxn)

print("="*80)
print("TRACE: 2-Bromopyridine Matching")
print("="*80)

if result.trace and 'evaluations' in result.trace:
    # Find the relevant entries
    target_ids = ["SCDB-SUZ-2HETARYL-SPHOS", "SCDB-SUZ-ARBR-ORTHO-XPhos"]
    
    for eval_data in result.trace['evaluations']:
        entry_id = eval_data.get('id', 'UNKNOWN')
        
        if entry_id in target_ids:
            matched = eval_data.get('matched', False)
            filtered = eval_data.get('filtered', False)
            priority = eval_data.get('priority', 0)
            
            status = "MATCHED" if matched else ("FILTERED" if filtered else "NO MATCH")
            
            print(f"\n{entry_id}")
            print(f"  Priority: {priority}")
            print(f"  Status: {status}")
            
            if 'matches' in eval_data:
                print(f"  SMARTS matches: {len(eval_data['matches'])}")
                for m in eval_data['matches']:
                    print(f"    - {m.get('smarts', 'N/A')}")
            
            if 'missing' in eval_data:
                print(f"  Missing patterns: {eval_data['missing']}")

print("\n" + "="*80)
print("ANALYSIS")
print("="*80)

print("\nThe issue:")
print("  2-Bromopyridine has nitrogen adjacent to the C-Br bond.")
print("  The ortho detector counts heteroatoms as 'ortho substituents'!")
print("  So pyridine nitrogen counts as an ortho substituent.")

print("\nSolution:")
print("  The 2-hetaryl rule needs HIGHER priority (65 > 50)")
print("  This ensures heteroaryls are caught before generic ortho rule.")

print("\nAlternatively:")
print("  Modify ortho feature detector to exclude ring heteroatoms")
print("  Only count exocyclic substituents as 'ortho substituents'")

print("\n" + "="*80)
