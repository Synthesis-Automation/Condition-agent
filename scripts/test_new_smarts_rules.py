#!/usr/bin/env python3
"""Test the newly added SMARTS-based rules."""

import sys
sys.path.insert(0, ".")

from chemtools.scdb_matcher import loader, match

# Load database
db = loader.load_db("data/conditionDB/suzuki_db.json")

# Test cases for new rules
test_cases = [
    {
        "name": "Ortho-substituted ArBr (existing rule)",
        "rxn": "Brc1ccccc1C.c1ccc(B(O)O)cc1>>Cc1ccccc1-c1ccccc1",
        "expected_entry": "SCDB-SUZ-ARBR-ORTHO-XPhos",
        "expected_ligand": "XPhos"
    },
    {
        "name": "Simple ArBr (should use default)",
        "rxn": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccccc1-c1ccccc1",
        "expected_entry": "SCDB-SUZ-DEFAULT-ArI_ArBr",
        "expected_ligand": "PPh3"
    },
    {
        "name": "2-Bromopyridine (new 2-hetaryl rule)",
        "rxn": "Brc1ccccn1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccn2)cc1",
        "expected_entry": "SCDB-SUZ-2HETARYL-SPHOS",
        "expected_ligand": "SPhos"
    },
]

print("="*80)
print("TESTING NEW SMARTS-BASED RULES")
print("="*80)

for i, test in enumerate(test_cases, 1):
    print(f"\n{i}. {test['name']}")
    print("-" * 80)
    print(f"Reaction: {test['rxn'][:60]}...")
    
    try:
        result = match(db, test['rxn'])
        
        print(f"\nResult:")
        print(f"  Entry ID: {result.entry_id}")
        print(f"  Match Type: {result.match_type}")
        print(f"  Priority: {result.priority}")
        
        # Check if it matches expected
        if result.entry_id == test['expected_entry']:
            print(f"  Status: OK - Matched expected entry")
        else:
            print(f"  Status: UNEXPECTED")
            print(f"  Expected: {test['expected_entry']}")
            print(f"  Got: {result.entry_id}")
        
        # Show conditions
        if hasattr(result, 'conditions') and result.conditions:
            pd_source = result.conditions.get('pd_source', ['N/A'])
            ligands = result.conditions.get('ligands', ['N/A'])
            temp = result.conditions.get('temperature_C', ['N/A'])
            
            print(f"\nRecommended Conditions:")
            print(f"  Pd Source: {pd_source[0] if isinstance(pd_source, list) else pd_source}")
            print(f"  Ligand: {ligands[0] if isinstance(ligands, list) else ligands}")
            print(f"  Temperature: {temp[0] if isinstance(temp, list) else temp} C")
        
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "="*80)
print("SUMMARY")
print("="*80)

# Count entries by type
scheme_count = sum(1 for e in db.entries if e.type == "scheme")
default_count = sum(1 for e in db.entries if e.type == "default_condition")
total = len(db.entries)

print(f"\nDatabase Statistics:")
print(f"  Total entries: {total}")
print(f"  Scheme entries: {scheme_count}")
print(f"  Default entries: {default_count}")

print(f"\nNew entries added:")
print(f"  - SCDB-SUZ-BULKY-NUC-XPHOS (bulky nucleophiles)")
print(f"  - SCDB-SUZ-2HETARYL-SPHOS (2-hetaryl halides)")
print(f"  - SCDB-SUZ-ARI-META-TBUXPHOS (meta-disubstituted)")

print(f"\nHybrid SMARTS + Features approach:")
print(f"  Simple SMARTS: Easy to read and write")
print(f"  Feature requirements: Precise selectivity")
print(f"  Priority system: Control rule precedence")

print("\n" + "="*80)
