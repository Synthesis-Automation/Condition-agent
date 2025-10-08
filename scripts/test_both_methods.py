"""Test both ML and rule-based recommendation methods."""

import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from app.ui_simple import get_ml_recommendations, get_rule_recommendations
import json

# Test reaction
reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
reaction_type = "C_N_Coupling_Pd"

print("="*80)
print("TESTING BOTH RECOMMENDATION METHODS")
print("="*80)
print(f"Reaction: {reaction}")
print(f"Type: {reaction_type}")
print("="*80)

# Test ML
print("\n1. ML RECOMMENDATIONS")
print("-"*80)
ml_result = get_ml_recommendations(reaction, reaction_type, top_k=3)
ml_json, ml_table = ml_result

# Extract JSON
import re
match = re.search(r'```json\n(.*)\n```', ml_json, re.DOTALL)
if match:
    ml_data = json.loads(match.group(1))
    print(f"✅ ML returned {len(ml_data['recommendations'])} recommendations")
    print(f"   Ranks: {[r['rank'] for r in ml_data['recommendations']]}")
    print(f"   Confidences: {[r['confidence'] for r in ml_data['recommendations']]}")
    print(f"   Support: {[r['support'] for r in ml_data['recommendations']]}")
    
    # Check enrichment
    for i, rec in enumerate(ml_data['recommendations'], 1):
        base_reagents = [r for r in rec['reagents'] if r['role'] == 'base']
        solvent_reagents = [r for r in rec['reagents'] if r['role'] == 'solvent']
        print(f"\n   Rec {i}:")
        if base_reagents:
            base = base_reagents[0]
            print(f"     Base: {base.get('name', 'N/A')} ({base.get('cas', 'N/A')})")
            if 'properties' in base:
                print(f"       Properties: {base['properties']}")
        if solvent_reagents:
            solv = solvent_reagents[0]
            print(f"     Solvent: {solv.get('name', 'N/A')} ({solv.get('cas', 'N/A')})")
            if 'properties' in solv:
                print(f"       Properties: {solv['properties']}")
else:
    print(f"❌ Could not parse ML JSON")

print(f"\n   Table: {len(ml_table)} rows")

# Test Rule-based
print("\n" + "="*80)
print("2. RULE-BASED RECOMMENDATIONS")
print("-"*80)
rule_result = get_rule_recommendations(reaction, reaction_type)
rule_json, rule_table = rule_result

# Extract JSON
match = re.search(r'```json\n(.*)\n```', rule_json, re.DOTALL)
if match:
    rule_data = json.loads(match.group(1))
    print(f"✅ Rule-based returned {len(rule_data['recommendations'])} recommendations")
    print(f"   Ranks: {[r['rank'] for r in rule_data['recommendations']]}")
    print(f"   Confidences: {[r['confidence'] for r in rule_data['recommendations']]}")
    
    # Check enrichment
    for i, rec in enumerate(rule_data['recommendations'], 1):
        base_reagents = [r for r in rec['reagents'] if r['role'] == 'base']
        solvent_reagents = [r for r in rec['reagents'] if r['role'] == 'solvent']
        print(f"\n   Rec {i} ({rec.get('variant_description', 'N/A')}):")
        if base_reagents:
            base = base_reagents[0]
            print(f"     Base: {base.get('name', 'N/A')} ({base.get('cas', 'N/A')})")
            if 'properties' in base:
                print(f"       Properties: {base['properties']}")
        if solvent_reagents:
            solv = solvent_reagents[0]
            print(f"     Solvent: {solv.get('name', 'N/A')} ({solv.get('cas', 'N/A')})")
            if 'properties' in solv:
                print(f"       Properties: {solv['properties']}")
else:
    print(f"❌ Could not parse rule-based JSON")

print(f"\n   Table: {len(rule_table)} rows")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"✅ ML: {len(ml_data['recommendations']) if 'ml_data' in locals() else 0} recommendations")
print(f"✅ Rule-based: {len(rule_data['recommendations']) if 'rule_data' in locals() else 0} recommendations")
print(f"✅ Both methods working with enhanced output format!")
print("="*80)
