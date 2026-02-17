"""
Get detailed LLM validation for the tandem reaction.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from reaction_agent.scripts.llm_assisted_mapping import hybrid_mapping_workflow
import json


rxn_smiles = "COC(OC)c1ccccc1Br.C=CC(CCCC)B1OC(C)(C)C(C)(C)O1>>C=CC(CCCC)c1ccccc1C=O"

print("=" * 100)
print("LLM Validation Details for Tandem Reaction")
print("=" * 100)
print()

result = hybrid_mapping_workflow(rxn_smiles, confidence_threshold=0.7)

if result.get('llm_validation'):
    val = result['llm_validation']

    print("CONFIDENCE ADJUSTMENT:")
    print("-" * 100)
    print(f"Original (rxnmapper): {val.get('original_confidence', 'N/A'):.3f}")
    print(f"Adjusted (after LLM): {val.get('adjusted_confidence', 'N/A'):.3f}")
    print(f"Boost: +{val.get('adjusted_confidence', 0) - val.get('original_confidence', 0):.3f}")
    print()

    if val.get('validation'):
        v = val['validation']

        print("VALIDATION RESULT:")
        print("-" * 100)
        print(f"Mapping valid: {v.get('mapping_valid', 'Unknown')}")
        print(f"Reaction type: {v.get('reaction_type', 'Unknown')}")
        print(f"Recommendation: {v.get('recommendation', 'Unknown')}")
        print()

        if v.get('issues_found'):
            print(f"Issues Found: {len(v['issues_found'])}")
            for i, issue in enumerate(v['issues_found'], 1):
                print(f"{i}. {issue}")
        else:
            print("Issues Found: None (mapping is chemically reasonable)")

        print()
        print("LLM REASONING:")
        print("-" * 100)
        reasoning = v.get('reasoning', 'No reasoning provided')
        print(reasoning)
        print()

    print("=" * 100)
    print()

    # Save to file
    with open('results/tandem_rxn_llm_validation.json', 'w') as f:
        json.dump(result, f, indent=2, default=str)

    print("✓ Saved detailed results to results/tandem_rxn_llm_validation.json")

elif result.get('llm_assistance'):
    print("LLM ASSISTANCE (not validation) was triggered")
    print("This means confidence was < 0.4")

else:
    print("No LLM analysis was triggered")
    print(f"Final confidence: {result.get('final_confidence', 'N/A')}")
