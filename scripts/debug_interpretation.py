"""Debug script to test automatic interpretation."""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from reaction_agent.core import analyze_reaction_smiles
import json

rxn_smiles = "COC(OC)c1ccccc1Br.C=CC(CCCC)B1OC(C)(C)C(C)(C)O1>>C=CC(CCCC)c1ccccc1C=O"

print("Testing automatic interpretation integration...")
print("=" * 80)
print()

result = analyze_reaction_smiles(rxn_smiles)

print("Result keys:", result.keys())
print()

if 'auto_interpretation' in result:
    print("✓ auto_interpretation found!")
    print()

    auto_interp = result['auto_interpretation']
    print("auto_interpretation keys:", auto_interp.keys())
    print()

    if 'error' in auto_interp:
        print(f"✗ Error: {auto_interp['error']}")
    elif 'interpretation' in auto_interp:
        print("✓ Interpretation successful!")
        interpretation = auto_interp['interpretation']
        print(f"  Complexity: {interpretation.get('reaction_complexity', 'N/A')}")
        print(f"  Patterns: {interpretation.get('patterns_detected', [])}")
        print(f"  Types: {interpretation.get('likely_reaction_types', [])}")
        print(f"  Tandem: {interpretation.get('tandem_reaction_suspected', False)}")
else:
    print("✗ auto_interpretation NOT found in result!")
    print()
    print("Available keys:", list(result.keys()))
