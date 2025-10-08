"""
Test Ni C-N Coupling Condition Database
========================================

Validates the cn_coupling_ni.json rule database with sample C-N coupling reactions.
Tests that rules match correctly and provide reasonable conditions.
"""

import sys
import json
from pathlib import Path
from typing import Dict, Any, List

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from chemtools.scdb_matcher import load_db, match
    HAS_MATCHER = True
except ImportError:
    HAS_MATCHER = False
    print("‚ö?scdb_matcher not available - skipping matcher tests")

# Sample C-N coupling reactions (from sample_reactions.py)
CN_COUPLING_REACTIONS = [
    # Primary anilines with various aryl halides
    ("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1", "Ph-Br + aniline"),
    ("Clc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1", "Ph-Cl + aniline"),
    ("Ic1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1", "Ph-I + aniline"),
    
    # Aryl chlorides with anilines (Ni strength)
    ("Clc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>c1ccccc1Nc1ccc(C(F)(F)F)cc1", "4-Cl-CF3-benzene + aniline"),
    ("Clc1ccncc1.Nc1ccccc1>>c1ccccc1Nc1ccncc1", "4-Cl-pyridine + aniline"),
    
    # Heterocyclic anilines
    ("Brc1ccccc1.Nc1ncccn1>>c1ccccc1Nc1ncccn1", "Ph-Br + 2-aminopyrimidine"),
    ("Clc1ccccc1.Nc1cccnc1>>c1ccc(Nc2cccnc2)cc1", "Ph-Cl + 3-aminopyridine"),
    
    # Secondary amines (aliphatic)
    ("Brc1ccccc1.C1CCNCC1>>c1ccccc1N1CCCCC1", "Ph-Br + piperidine"),
    ("Clc1ccccc1.C1CNCCN1>>c1ccccc1N1CCNCC1", "Ph-Cl + piperazine"),
    
    # Heteroaryl halides
    ("Brc1ccncc1.Nc1ccccc1>>c1ccccc1Nc1ccncc1", "4-Br-pyridine + aniline"),
    ("Clc1ncccn1.Nc1ccccc1>>c1ccccc1Nc1ncccn1", "2-Cl-pyrimidine + aniline"),
    
    # Functional groups present
    ("Brc1ccc(C(=O)OC)cc1.Nc1ccccc1>>COC(=O)c1ccc(Nc2ccccc2)cc1", "4-Br-methyl benzoate + aniline"),
    ("Clc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1", "4-Cl-benzonitrile + aniline"),
]


def validate_json_structure(db_path: Path) -> List[str]:
    """Validate JSON structure and identify potential issues."""
    issues = []
    
    try:
        with open(db_path, 'r', encoding='utf-8') as f:
            db = json.load(f)
    except json.JSONDecodeError as e:
        return [f"‚ù?JSON parsing error: {e}"]
    except FileNotFoundError:
        return [f"‚ù?File not found: {db_path}"]
    
    # Check top-level structure
    if "entries" not in db:
        issues.append("‚ù?Missing 'entries' key in root")
    if "reaction_type" not in db:
        issues.append("‚ö?Missing 'reaction_type' in root (optional but recommended)")
    
    entries = db.get("entries", [])
    print(f"\nüìä Database Statistics:")
    print(f"   Total entries: {len(entries)}")
    
    # Count entry types
    rule_entries = [e for e in entries if e.get("id", "").startswith("SCDB-")]
    scheme_entries = [e for e in entries if e.get("type") == "scheme" or e.get("id", "").startswith("M-")]
    default_entries = [e for e in entries if e.get("role") == "default_condition"]
    
    print(f"   Rule entries: {len(rule_entries)}")
    print(f"   Scheme entries: {len(scheme_entries)}")
    print(f"   Default conditions: {len(default_entries)}")
    
    # Validate each entry
    for i, entry in enumerate(entries):
        entry_id = entry.get("id", f"entry_{i}")
        
        # Check required fields based on type
        if entry.get("type") == "scheme" or entry.get("id", "").startswith("M-"):
            # Scheme-based entry
            if "reactant_smarts" not in entry:
                issues.append(f"‚ö?Scheme entry '{entry_id}' missing 'reactant_smarts'")
            if "conditions" not in entry:
                issues.append(f"‚ö?Scheme entry '{entry_id}' missing 'conditions'")
        
        elif entry.get("role") == "default_condition":
            # Default condition entry
            if "default_condition" not in entry:
                issues.append(f"‚ù?Default entry '{entry_id}' missing 'default_condition' block")
        
        else:
            # Rule-based entry (SCDB-*)
            if "rxn_smiles_min" not in entry:
                issues.append(f"‚ö?Rule entry '{entry_id}' missing 'rxn_smiles_min'")
            if "conditions" not in entry:
                issues.append(f"‚ù?Rule entry '{entry_id}' missing 'conditions'")
            
            # Check env.features_from_smiles
            env = entry.get("env", {})
            if "features_from_smiles" not in env:
                issues.append(f"‚ö?Rule entry '{entry_id}' missing 'env.features_from_smiles'")
        
        # Check conditions structure (if present)
        conditions = entry.get("conditions", {})
        if conditions:
            # Ni-specific checks
            if "ni_source" not in conditions and entry.get("role") != "default_condition":
                issues.append(f"‚ö?Entry '{entry_id}' missing 'ni_source' in conditions")
            if "ligands" not in conditions and "ligand" not in conditions:
                issues.append(f"‚ö?Entry '{entry_id}' missing 'ligands' in conditions")
    
    return issues


def test_with_scdb_matcher(db_path: Path, reactions: List[tuple]) -> None:
    """Test reactions with scdb_matcher if available."""
    if not HAS_MATCHER:
        return
    
    print("\n" + "=" * 70)
    print("Testing with scdb_matcher")
    print("=" * 70)
    
    try:
        db = load_db(str(db_path))
        print(f"‚ú?Loaded database: {len(db.entries)} entries")
    except Exception as e:
        print(f"‚ù?Failed to load database: {e}")
        return
    
    matches_found = 0
    no_matches = 0
    
    for rxn_smiles, description in reactions:
        print(f"\nüß™ Testing: {description}")
        print(f"   Reaction: {rxn_smiles}")
        
        try:
            result = match(db, rxn_smiles)
            
            if result.get("matched"):
                matches_found += 1
                rule_id = result.get("matched_rule_id", "unknown")
                rule_name = result.get("matched_rule_name", "")
                confidence = result.get("confidence", 0)
                
                print(f"   ‚ú?MATCHED: {rule_id}")
                if rule_name:
                    print(f"     Rule: {rule_name}")
                print(f"     Confidence: {confidence:.2f}")
                
                # Show recommended conditions
                conds = result.get("conditions", {})
                if conds:
                    print(f"     Ni source: {conds.get('ni_source', 'N/A')}")
                    print(f"     Ligands: {conds.get('ligands', 'N/A')}")
                    print(f"     Base: {conds.get('base', 'N/A')}")
                    print(f"     Solvent: {conds.get('solvent', 'N/A')}")
                    temp = conds.get('temperature_C', 'N/A')
                    time = conds.get('time_h', 'N/A')
                    print(f"     Conditions: {temp}¬∞C, {time}h")
            else:
                no_matches += 1
                print(f"   ‚ú?NO MATCH")
                reason = result.get("reason", "Unknown")
                print(f"     Reason: {reason}")
        
        except Exception as e:
            print(f"   ‚ù?ERROR: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n" + "=" * 70)
    print(f"Summary: {matches_found} matched, {no_matches} no match")
    print("=" * 70)


def manual_pattern_check(db_path: Path, reactions: List[tuple]) -> None:
    """Manually check which rules should match (without scdb_matcher)."""
    print("\n" + "=" * 70)
    print("Manual Pattern Analysis")
    print("=" * 70)
    
    with open(db_path, 'r', encoding='utf-8') as f:
        db = json.load(f)
    
    entries = db.get("entries", [])
    
    # Categorize reactions
    print("\nüìã Reaction Categories:")
    
    arcl_aniline = [r for r in reactions if "Cl" in r[0] and "Nc1ccccc1" in r[0]]
    arbr_aniline = [r for r in reactions if "Br" in r[0] and "Nc1ccccc1" in r[0]]
    ari_aniline = [r for r in reactions if "Ic1ccccc1" in r[0] and "Nc1ccccc1" in r[0]]
    aliphatic_amine = [r for r in reactions if "CNCC" in r[0] or "C1CCNCC1" in r[0]]
    
    print(f"   ArCl + aniline: {len(arcl_aniline)}")
    print(f"   ArBr + aniline: {len(arbr_aniline)}")
    print(f"   ArI + aniline: {len(ari_aniline)}")
    print(f"   Aliphatic amines: {len(aliphatic_amine)}")
    
    # Check rule coverage
    print("\nüìê Rule Coverage:")
    
    arcl_rules = [e for e in entries if "ArCl" in e.get("token_signature", []) or "ArCl" in e.get("name", "")]
    arbr_rules = [e for e in entries if "ArBr" in e.get("token_signature", []) or "ArBr" in e.get("name", "") or "ArI" in e.get("name", "")]
    otf_rules = [e for e in entries if "OTf" in e.get("token_signature", []) or "sulfonate" in e.get("name", "").lower()]
    photoredox_rules = [e for e in entries if "photoredox" in e.get("token_signature", []) or "photoredox" in e.get("name", "").lower()]
    
    print(f"   ArCl-specific rules: {len(arcl_rules)}")
    for rule in arcl_rules[:3]:  # Show first 3
        print(f"     - {rule.get('id')}: {rule.get('name', 'N/A')}")
    
    print(f"   ArBr/ArI rules: {len(arbr_rules)}")
    for rule in arbr_rules[:3]:
        print(f"     - {rule.get('id')}: {rule.get('name', 'N/A')}")
    
    print(f"   Sulfonate/triflate rules: {len(otf_rules)}")
    print(f"   Photoredox rules: {len(photoredox_rules)}")


def main():
    """Main test function."""
    print("=" * 70)
    print("Ni C-N Coupling ConditionDB Validation")
    print("=" * 70)
    
    db_path = Path("data/conditionDB/cn_coupling_ni.json")
    
    if not db_path.exists():
        print(f"\n‚ù?Database file not found: {db_path}")
        return 1
    
    print(f"\n‚ú?Found database: {db_path}")
    
    # 1. Validate JSON structure
    print("\n" + "=" * 70)
    print("JSON Structure Validation")
    print("=" * 70)
    
    issues = validate_json_structure(db_path)
    
    if not issues:
        print("\n‚ú?No structural issues found!")
    else:
        print(f"\n‚ö?Found {len(issues)} potential issues:")
        for issue in issues:
            print(f"   {issue}")
    
    # 2. Manual pattern analysis
    manual_pattern_check(db_path, CN_COUPLING_REACTIONS)
    
    # 3. Test with scdb_matcher (if available)
    if HAS_MATCHER:
        test_with_scdb_matcher(db_path, CN_COUPLING_REACTIONS)
    else:
        print("\n‚ö?scdb_matcher not available - install to run full tests")
    
    print("\n" + "=" * 70)
    print("‚ú?Validation Complete")
    print("=" * 70)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
