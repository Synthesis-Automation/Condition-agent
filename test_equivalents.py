#!/usr/bin/env python
"""Quick test for equivalents calculation in rule-based recommendations."""

import sys
import json
from pathlib import Path

# Add project root to path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from chemtools import chem, output_formatter

def test_equivalents():
    reaction_smiles = "c1ccc(B(O)O)cc1.Brc1ccc(C(F)(F)F)cc1>>c1ccc(-c2ccc(C(F)(F)F)cc2)cc1"
    
    print("Testing equivalents calculation...")
    print(f"Reaction: {reaction_smiles}\n")
    
    # Load rule database
    db_path = Path("data/conditionDB/suzuki_db.json")
    if not db_path.exists():
        print(f"Database not found: {db_path}")
        return
    
    db = chem.rules.load_database(str(db_path), cache=True)
    matches = chem.rules.match(db, reaction_smiles)
    
    if not matches:
        print("No matches found")
        return
    
    # Convert MatchResult to dict
    match_dict = matches.to_json_dict() if hasattr(matches, 'to_json_dict') else matches
    
    # Format matches
    result = output_formatter.format_rule_match_result(reaction_smiles, match_dict)
    
    if result and "recommended_conditions" in result:
        print(f"Found {len(result['recommended_conditions'])} recommendations\n")
        
        # Save to JSON for inspection
        output_file = Path("test_equivalents_output.json")
        with open(output_file, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"Full output saved to: {output_file}\n")
        
        for i, rec in enumerate(result["recommended_conditions"][:3], 1):
            print(f"=== Recommendation {i} ===")
            print(f"Rule: {rec.get('rule_id', 'N/A')}")
            print(f"Chemicals:")
            
            for chem_entry in rec.get("chemicals", []):
                name = chem_entry.get("name") or "Unknown"
                role = chem_entry.get("role") or "Unknown"
                equiv = chem_entry.get("equivalents")
                
                equiv_str = f"{equiv:.3f}" if equiv is not None else "None"
                print(f"  - {name:30s} [{role:15s}] equivalents: {equiv_str}")
            print()
    else:
        print("No recommendations found")
        print(json.dumps(result, indent=2))

if __name__ == "__main__":
    test_equivalents()
