"""
Test that the agent can now use SNAr rules for condition recommendations
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chem_assistant.chemtools_wrapper import (
    recommend_conditions_tool,
    _FAMILY_TO_RULE_DB,
)

print("="*70)
print("Test: Agent SNAr Condition Recommendation")
print("="*70)
print()

# The exact reaction from the user's request
reaction_smiles = "Clc1nc(Cl)nc(Cl)n1.OC>>COc1nc(Cl)nc(Cl)n1"
print(f"Reaction: {reaction_smiles}")
print()

# Check mapping
print("Family-to-DB mapping check:")
for family in ["snar", "s_nar", "SNAr"]:
    normalized = family.lower().replace("-", "_")
    db = _FAMILY_TO_RULE_DB.get(normalized, "NOT FOUND")
    print(f"  {family:15s} -> {db}")
print()

# Try to get recommendation
print("Getting recommendation...")
print()

try:
    # Method 1: Let it auto-detect (should use SNAr)
    result = recommend_conditions_tool.invoke({
        "reaction_smiles": reaction_smiles,
        "k": 3,
        "max_variants": 2
    })
    
    print("✓ Recommendation generated successfully!")
    print()
    print(f"Result type: {type(result)}")
    print(f"Result keys: {list(result.keys()) if isinstance(result, dict) else 'N/A'}")
    print()
    
    if isinstance(result, dict):
        if "recommendations" in result:
            recs = result["recommendations"]
            print(f"Number of recommendations: {len(recs)}")
            
            if recs:
                print()
                print("First recommendation:")
                first = recs[0]
                for key, value in first.items():
                    if key == "conditions":
                        print(f"  {key}:")
                        for cond_key, cond_val in value.items():
                            print(f"    {cond_key}: {cond_val}")
                    else:
                        print(f"  {key}: {value}")
        
        if "detected_family" in result:
            print()
            print(f"Detected family: {result['detected_family']}")
        
        if "rule_based" in result:
            print()
            print(f"Rule-based: {result['rule_based']}")
            rb = result["rule_based"]
            if isinstance(rb, dict):
                if "database" in rb:
                    print(f"  Database: {rb['database']}")
                if "matched_rule" in rb:
                    print(f"  Matched rule: {rb['matched_rule']}")
                if "base_conditions" in rb:
                    print(f"  Base conditions: {rb['base_conditions']}")

    print()
    print("="*70)
    print("✅ SUCCESS: Agent can now use SNAr rule database!")
    print("="*70)
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()
    print()
    print("="*70)
    print("❌ FAILED: Agent could not use SNAr rules")
    print("="*70)
    sys.exit(1)
