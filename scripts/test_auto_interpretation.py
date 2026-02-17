"""
Enhanced hybrid mapping with automatic interpretation.

Adds instant reaction pattern analysis to help users understand results.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from chemtools._atom_mapping import analyze_bond_changes_hybrid
from chemtools.reaction_interpreter import interpret_reaction_pattern, format_interpretation_report
from typing import Dict, Any
import time


def analyze_with_interpretation(
    rxn_smiles: str,
    use_rxnmapper: bool = True,
    use_local_env: bool = True,
    use_mcs: bool = True,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Enhanced hybrid analysis with automatic interpretation.

    Performs mapping with all methods, then automatically interprets
    the results to explain what type of reaction it is and why methods
    might disagree.

    Args:
        rxn_smiles: Reaction SMILES
        use_rxnmapper: Try RXNMapper
        use_local_env: Try local environment mapping
        use_mcs: Try MCS analysis
        verbose: Print interpretation report

    Returns:
        Dictionary with mapping results + interpretation
    """

    start_time = time.time()

    # Step 1: Run hybrid mapping
    hybrid_result = analyze_bond_changes_hybrid(
        rxn_smiles,
        use_rxnmapper=use_rxnmapper,
        use_local_env=use_local_env,
        use_mcs=use_mcs,
        use_manual=True,
        auto_map=True
    )

    # Step 2: Interpret results (fast, deterministic)
    interpretation = interpret_reaction_pattern(rxn_smiles, hybrid_result)

    # Step 3: Generate report
    report = format_interpretation_report(rxn_smiles, hybrid_result, interpretation)

    elapsed = time.time() - start_time

    result = {
        'rxn_smiles': rxn_smiles,
        'hybrid_result': hybrid_result,
        'interpretation': interpretation,
        'report': report,
        'analysis_time_seconds': elapsed
    }

    if verbose:
        print(report)
        print(f"\n⚡ Analysis completed in {elapsed:.2f} seconds")
        print()

    return result


if __name__ == "__main__":
    # Test cases
    test_reactions = [
        {
            "name": "Tandem: Suzuki + Acetal Hydrolysis",
            "smiles": "COC(OC)c1ccccc1Br.C=CC(CCCC)B1OC(C)(C)C(C)(C)O1>>C=CC(CCCC)c1ccccc1C=O"
        },
        {
            "name": "Suzuki Coupling (Complex)",
            "smiles": "CC1(C)CC=C(B2OC(C)(C)C(C)(C)O2)CC1.CS(=O)(=O)N1CCN(C(=O)c2cnc3ccc(F)cc3c2Br)CC1>>CC1(C)CC=C(c2c(C(=O)N3CCN(S(C)(=O)=O)CC3)cnc3ccc(F)cc23)CC1.O=C(O)C(F)(F)F"
        },
        {
            "name": "Simple SN2",
            "smiles": "CCCBr.N>>CCCN"
        },
        {
            "name": "Esterification",
            "smiles": "CC(=O)O.CCO>>CC(=O)OCC"
        }
    ]

    print("=" * 100)
    print("TESTING AUTOMATIC REACTION INTERPRETATION")
    print("=" * 100)
    print()

    for i, test in enumerate(test_reactions, 1):
        print(f"\n{'='*100}")
        print(f"TEST {i}: {test['name']}")
        print(f"{'='*100}\n")
        print(f"SMILES: {test['smiles']}\n")

        result = analyze_with_interpretation(
            test['smiles'],
            verbose=True
        )

        print()
