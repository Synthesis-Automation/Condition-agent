
import sys
from pathlib import Path
from collections import Counter

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from examples.sample_reactions import SAMPLE_REACTIONS
from chemtools.rule import RuleEngine
from chemtools.taxonomy.rule_db import resolve_rule_db_v2
from chemtools import detect_reaction

def test_all_cn_couplings():
    # Filter for C-N coupling reactions
    cn_reactions = []
    for rxn in SAMPLE_REACTIONS:
        if not isinstance(rxn, str) or ">>" not in rxn:
            continue
        
        # Look for C-N or Buchwald or Ullmann or Chan-Lam in the description
        if any(keyword in rxn.lower() for keyword in ["c-n", "buchwald", "ullmann", "chan-lam"]):
            cn_reactions.append(rxn)
    
    total = len(cn_reactions)
    print(f"Found {total} C-N coupling reactions in sample_reactions.py\n")
    
    stats = {
        "success": 0,
        "no_db": 0,
        "no_rule": 0,
        "error": 0,
        "families": Counter(),
        "rules": Counter(),
        "failures": []
    }
    
    # Cache engines to avoid reloading
    engines = {}
    rule_dir = Path("data/rule_db_v2")

    for i, rxn_str in enumerate(cn_reactions):
        smiles = rxn_str.split(" (")[0]
        desc = rxn_str.split(" (")[1].rstrip(")") if " (" in rxn_str else "No description"
        
        try:
            # Detect reaction type
            detection = detect_reaction(smiles)
            family = detection.get("family", "unknown")
            stats["families"][family] += 1
            
            # Resolve DB
            db_stem = resolve_rule_db_v2(family)
            if not db_stem:
                stats["no_db"] += 1
                stats["failures"].append(f"{desc}: No DB for family {family}")
                continue
                
            db_path = rule_dir / f"{db_stem}.json"
            if not db_path.exists():
                stats["no_db"] += 1
                stats["failures"].append(f"{desc}: DB file {db_path} not found")
                continue
            
            if db_path not in engines:
                engines[db_path] = RuleEngine.from_file(db_path)
            
            engine = engines[db_path]
            
            # Get recommendation
            recommendation = engine.recommend(smiles)
            
            if recommendation:
                stats["success"] += 1
                stats["rules"][f"{db_stem} -> {recommendation.base_rule.name}"] += 1
            else:
                stats["no_rule"] += 1
                stats["failures"].append(f"{desc}: No rule matched in {db_stem}")
                
        except Exception as e:
            stats["error"] += 1
            stats["failures"].append(f"{desc}: Error: {str(e)}")

    # Print Summary
    print("="*50)
    print("C-N COUPLING TEST SUMMARY")
    print("="*50)
    print(f"Total Reactions Tested: {total}")
    print(f"Successful Recs:      {stats['success']} ({stats['success']/total*100:.1f}%)")
    print(f"No DB Found:          {stats['no_db']}")
    print(f"No Rule Matched:      {stats['no_rule']}")
    print(f"Errors:               {stats['error']}")
    print("\nDetected Families:")
    for fam, count in stats["families"].most_common():
        print(f"  - {fam}: {count}")
    
    print("\nTop Applied Rules:")
    for rule, count in stats["rules"].most_common(10):
        print(f"  - {rule}: {count}")
    
    if stats["failures"]:
        print("\nSample Failures (first 5):")
        for fail in stats["failures"][:5]:
            print(f"  - {fail}")

if __name__ == "__main__":
    test_all_cn_couplings()
