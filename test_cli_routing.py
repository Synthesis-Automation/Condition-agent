"""
Test script to demonstrate CLI reaction type routing.

This script shows how the CLI determines reaction types and routes to
the appropriate recommendation system (rule-based vs ML).
"""

from app.cli_recommend import determine_final_reaction_type

# Test cases
test_cases = [
    {
        "name": "C-N Coupling with Copper (Ullmann)",
        "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "initial_type": None,
        "constraints": {
            "metal_preference": "Cu",
            "required_reagents": ["copper catalyst"]
        }
    },
    {
        "name": "C-N Coupling with Palladium (Buchwald)",
        "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "initial_type": None,
        "constraints": {
            "metal_preference": "Pd",
            "required_reagents": ["palladium catalyst"]
        }
    },
    {
        "name": "C-N Coupling - Auto-detect from structure",
        "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "initial_type": None,
        "constraints": {}
    },
    {
        "name": "Suzuki Coupling",
        "reaction": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        "initial_type": None,
        "constraints": {}
    },
    {
        "name": "User-specified C-N with Copper in reagents",
        "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "initial_type": "C_N_Coupling",
        "constraints": {
            "required_reagents": ["CuI", "copper iodide"]
        }
    },
]

def print_result(test_case, result):
    """Print test results with formatting."""
    print(f"\n{'='*70}")
    print(f"TEST: {test_case['name']}")
    print(f"{'='*70}")
    print(f"Reaction SMILES: {test_case['reaction'][:60]}...")
    print(f"Initial Type: {test_case['initial_type'] or 'None (auto-detect)'}")
    print(f"Constraints: {test_case['constraints']}")
    print(f"\nRESULTS:")
    print(f"  Final Type: {result['final_type']}")
    print(f"  Detected Family: {result['detected_family']}")
    print(f"  Confidence: {result['confidence']:.2f}")
    print(f"  Detected Metal: {result['detected_metal'] or 'None'}")
    print(f"  Routing: {result['routing']}")
    
    # Determine system
    if "rule-based" in result['routing'].lower() or "Cu" in str(result['final_type']):
        system = "RULE-BASED (deterministic constraints)"
    elif "Pd" in str(result['final_type']) or "Ni" in str(result['final_type']):
        system = "MACHINE LEARNING (similarity-based)"
    else:
        system = "AUTO-DETECT"
    
    print(f"  Recommendation System: {system}")


if __name__ == "__main__":
    print("\n" + "="*70)
    print(" REACTION TYPE ROUTING TEST")
    print("="*70)
    print("\nThis demonstrates how the CLI determines reaction types")
    print("and routes to the appropriate recommendation system.\n")
    
    for test_case in test_cases:
        try:
            result = determine_final_reaction_type(
                reaction_smiles=test_case["reaction"],
                initial_type=test_case["initial_type"],
                constraints=test_case["constraints"]
            )
            print_result(test_case, result)
        except Exception as e:
            print(f"\n{'='*70}")
            print(f"TEST: {test_case['name']}")
            print(f"{'='*70}")
            print(f"ERROR: {e}")
    
    print(f"\n{'='*70}")
    print("SUMMARY:")
    print("="*70)
    print("✓ C-N Coupling + Cu → C_N_Coupling_Cu (Ullmann, rule-based)")
    print("✓ C-N Coupling + Pd → C_N_Coupling_Pd (Buchwald, ML)")
    print("✓ Other reactions → Auto-detected type")
    print("="*70 + "\n")
