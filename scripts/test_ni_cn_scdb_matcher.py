"""
Test Ni C-N Coupling SchemeConditionDB Matcher
==============================================

Tests the new cn_coupling_ni_db.json rule database against
sample C-N coupling reactions.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.scdb_matcher import loader, matcher


def get_sample_cn_reactions():
    """Get sample C-N coupling reactions for testing."""
    return [
        {
            "name": "ArBr + aniline",
            "rxn_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
        },
        {
            "name": "ArCl + aniline",
            "rxn_smiles": "Clc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
        },
        {
            "name": "ArI + aniline",
            "rxn_smiles": "Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
        },
        {
            "name": "ArBr + pyrrolidine",
            "rxn_smiles": "Brc1ccccc1.C1CCNC1>>C1CCN(c2ccccc2)C1"
        },
        {
            "name": "ArCl + pyrrolidine",
            "rxn_smiles": "Clc1ccccc1.C1CCNC1>>C1CCN(c2ccccc2)C1"
        },
        {
            "name": "ArBr + morpholine",
            "rxn_smiles": "Brc1ccccc1.C1COCCN1>>C1COCCN1c1ccccc1"
        },
        {
            "name": "ArCl + morpholine",
            "rxn_smiles": "Clc1ccccc1.C1COCCN1>>C1COCCN1c1ccccc1"
        },
        {
            "name": "ArBr + piperidine",
            "rxn_smiles": "Brc1ccccc1.C1CCNCC1>>C1CCN(c2ccccc2)CC1"
        },
        {
            "name": "ArCl + indole",
            "rxn_smiles": "Clc1ccccc1.c1c[nH]c2ccccc12>>c1c(-c2ccccc2)[nH]c2ccccc12"
        },
        {
            "name": "ArBr + carbazole",
            "rxn_smiles": "Brc1ccccc1.c1ccc2c(c1)[nH]c1ccccc12>>c1ccc2c(c1)n(-c1ccccc1)c1ccccc12"
        },
        {
            "name": "Heteroaryl Cl + aniline",
            "rxn_smiles": "Clc1ccccn1.Nc1ccccc1>>c1ccc(Nc2ccccn2)cc1"
        },
        {
            "name": "Heteroaryl Br + pyrrolidine",
            "rxn_smiles": "Brc1ccccn1.C1CCNC1>>C1CCN(c2ccccn2)C1"
        },
        {
            "name": "ArOTf + aniline",
            "rxn_smiles": "c1ccc(OS(=O)(=O)C(F)(F)F)cc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
        },
        {
            "name": "ArCl + sec-amine (dibutylamine)",
            "rxn_smiles": "Clc1ccccc1.CCCCNCCCC>>CCCCN(CCCC)c1ccccc1"
        },
        {
            "name": "ArBr + imidazole",
            "rxn_smiles": "Brc1ccccc1.c1cnc[nH]1>>c1ccc(-n2ccnc2)cc1"
        },
        {
            "name": "ArCl + benzimidazole",
            "rxn_smiles": "Clc1ccccc1.c1ccc2[nH]cnc2c1>>c1ccc2nc(-c3ccccc3)n2c1"
        },
    ]


def test_ni_cn_rules():
    """Test all Ni C-N coupling rules against sample reactions."""
    
    # Load Ni C-N coupling database
    db_path = Path("data/conditionDB/cn_coupling_ni.json")
    
    print("=" * 80)
    print("Testing Ni C-N Coupling SchemeConditionDB")
    print("=" * 80)
    
    if not db_path.exists():
        print(f"â?Database file not found: {db_path}")
        return
    
    # Load database
    try:
        db = loader.load_db(str(db_path))
        print(f"\nâœ?Loaded database: {len(db.entries)} rules")
        print(f"   Reaction type: {db.reaction_type}")
        print(f"   Schema version: {db.schema_version}")
    except Exception as e:
        print(f"â?Failed to load database: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Get C-N coupling sample reactions
    cn_reactions = get_sample_cn_reactions()
    print(f"\nğŸ“‹ Testing against {len(cn_reactions)} C-N coupling reactions")
    
    # Test each reaction
    results = {
        "matched": [],
        "no_match": [],
        "errors": []
    }
    
    print("\n" + "=" * 80)
    print("Testing Individual Reactions")
    print("=" * 80)
    
    for i, rxn in enumerate(cn_reactions, 1):
        rxn_smiles = rxn["rxn_smiles"]
        name = rxn.get("name", f"Reaction {i}")
        
        print(f"\n[{i}/{len(cn_reactions)}] {name}")
        print(f"  SMILES: {rxn_smiles}")
        
        try:
            # Try to match
            result = matcher.match(db, rxn_smiles)
            
            # result is a MatchResult object, not a dict
            if result:
                print(f"  âœ?MATCHED:")
                print(f"     Rule ID: {result.entry_id}")
                print(f"     Name: {result.entry_name}")
                print(f"     Match type: {result.match_type}")
                print(f"     Priority: {result.priority}")
                
                # Show recommended conditions
                conds = result.conditions
                if conds:
                    catalyst = conds.get("catalyst", {})
                    if catalyst:
                        cat_name = catalyst.get("name", "")
                        ligand = catalyst.get("ligand", {})
                        if ligand and ligand.get("name"):
                            print(f"     Catalyst: {cat_name} + {ligand['name']}")
                        else:
                            print(f"     Catalyst: {cat_name}")
                    
                    base = conds.get("base", {})
                    if base and base.get("name"):
                        print(f"     Base: {base['name']}")
                    
                    solvent = conds.get("solvent", {})
                    if solvent and solvent.get("name"):
                        print(f"     Solvent: {solvent['name']}")
                    
                    temp = conds.get("temperature_c")
                    if temp:
                        print(f"     Temperature: {temp}Â°C")
                
                results["matched"].append({
                    "name": name,
                    "smiles": rxn_smiles,
                    "rule_id": result.entry_id,
                    "rule_name": result.entry_name
                })
            else:
                print(f"  âš ï¸  NO MATCH")
                results["no_match"].append({
                    "name": name,
                    "smiles": rxn_smiles
                })
        
        except Exception as e:
            print(f"  â?ERROR: {e}")
            results["errors"].append({
                "name": name,
                "smiles": rxn_smiles,
                "error": str(e)
            })
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    total = len(cn_reactions)
    matched = len(results["matched"])
    no_match = len(results["no_match"])
    errors = len(results["errors"])
    
    print(f"\nTotal reactions tested: {total}")
    print(f"  âœ?Matched: {matched} ({matched/total*100:.1f}%)")
    print(f"  âš ï¸  No match: {no_match} ({no_match/total*100:.1f}%)")
    print(f"  â?Errors: {errors} ({errors/total*100:.1f}%)")
    
    # Show no-match reactions
    if results["no_match"]:
        print(f"\nâš ï¸  Reactions with no match ({len(results['no_match'])}):")
        for item in results["no_match"]:
            print(f"  - {item['name']}")
    
    # Show errors
    if results["errors"]:
        print(f"\nâ?Reactions with errors ({len(results['errors'])}):")
        for item in results["errors"]:
            print(f"  - {item['name']}: {item['error']}")
    
    # Test specific Ni-relevant reactions
    print("\n" + "=" * 80)
    print("Testing Specific Ni-Relevant Reactions")
    print("=" * 80)
    
    ni_specific = [
        {
            "name": "ArCl + aniline (thermal Ni)",
            "rxn_smiles": "Clc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
        },
        {
            "name": "ArBr + aniline (thermal Ni)",
            "rxn_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
        },
        {
            "name": "ArCl + pyrrolidine (thermal Ni)",
            "rxn_smiles": "Clc1ccccc1.C1CCNC1>>C1CCN(c2ccccc2)C1"
        },
        {
            "name": "ArCl + morpholine (thermal Ni)",
            "rxn_smiles": "Clc1ccccc1.C1COCCN1>>C1COCCN1c1ccccc1"
        },
        {
            "name": "Heteroaryl chloride + aniline",
            "rxn_smiles": "Clc1ccccn1.Nc1ccccc1>>c1ccc(Nc2ccccn2)cc1"
        },
        {
            "name": "ArOTf + aniline (Ni catalysis)",
            "rxn_smiles": "c1ccc(OS(=O)(=O)C(F)(F)F)cc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
        }
    ]
    
    for rxn in ni_specific:
        print(f"\nğŸ§ª {rxn['name']}")
        print(f"   {rxn['rxn_smiles']}")
        
        try:
            result = matcher.match(db, rxn["rxn_smiles"])
            
            if result:
                print(f"   âœ?Matched")
                print(f"      Rule: {result.entry_id} - {result.entry_name}")
                print(f"      Match type: {result.match_type}")
            else:
                print(f"   âš ï¸  No match found")
        
        except Exception as e:
            print(f"   â?Error: {e}")
    
    print("\n" + "=" * 80)
    print("âœ?Test complete!")
    print("=" * 80)
    
    return results


def validate_database_structure():
    """Validate the database structure and SMARTS patterns."""
    
    db_path = Path("data/conditionDB/cn_coupling_ni.json")
    
    print("\n" + "=" * 80)
    print("Validating Database Structure")
    print("=" * 80)
    
    try:
        db = loader.load_db(str(db_path))
        
        print(f"\nâœ?Database loaded successfully")
        print(f"   Total entries: {len(db.entries)}")
        print(f"   Reaction type: {db.reaction_type}")
        print(f"   Schema version: {db.schema_version}")
        
        # Check entry details
        print(f"\nğŸ“Š Entry Details:")
        schemes = db.schemes()
        defaults = db.defaults()
        print(f"   Schemes: {len(schemes)}")
        print(f"   Defaults: {len(defaults)}")
        
        # List all scheme IDs
        print(f"\nğŸ“‹ Scheme Rules:")
        for i, entry in enumerate(schemes, 1):
            print(f"   {i}. {entry.id}: {entry.name}")
        
        if defaults:
            print(f"\nğŸ“‹ Default Conditions:")
            for i, entry in enumerate(defaults, 1):
                print(f"   {i}. {entry.id}: {entry.name}")
        
        print(f"\nâœ?Validation complete - no structural issues")
        return True
    
    except Exception as e:
        print(f"\nâ?Validation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("Ni C-N Coupling SchemeConditionDB Test Suite\n")
    
    # First validate structure
    valid = validate_database_structure()
    
    if valid:
        # Then test matching
        test_ni_cn_rules()
    else:
        print("\nDatabase has structural issues, skipping match tests")
