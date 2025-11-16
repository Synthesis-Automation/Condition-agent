"""
Test script to verify HTE agent can handle reaction SMILES queries
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Import the HTE tool
from chem_assistant import chemtools_wrapper

def test_reaction_smiles_with_data():
    """Test with ArBr + ArNH2 (has 112 experiments with Cu)"""
    print("=" * 80)
    print("Test 1: Reaction SMILES with available data (ArBr + ArNH2 + Cu)")
    print("=" * 80)
    
    # Call the tool invoke method
    result = chemtools_wrapper.hte_recommend_tool.invoke({
        "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        "reaction_type_filter": "C_N_Coupling",
        "catalyst_filter": "Cu",
        "top_k": 5,
        "min_experiments": 1
    })
    
    print(f"\nSuccess: {result.get('success')}")
    
    if result.get('success'):
        # _success_response merges dict directly, so fields are at top level
        print(f"Reactant A: {result.get('reactant_a_smiles')} → Type: {result.get('reactant_a_type')}")
        print(f"Reactant B: {result.get('reactant_b_smiles')} → Type: {result.get('reactant_b_type')}")
        print(f"Predicted Reaction: {result.get('predicted_reaction_type')}")
        print(f"Matching Experiments: {result.get('matching_experiments')}")
        
        recommendations = result.get('recommendations', [])
        if recommendations:
            print(f"\nTop 3 Recommendations:")
            for rec in recommendations[:3]:
                print(f"\n{rec['rank']}. Z-Score: {rec['avg_z_score']}")
                print(f"   {rec['catalyst']} / {rec['ligand']} / {rec['base']} / {rec['solvent']}")
                print(f"   Avg Yield: {rec['avg_yield']}%, Experiments: {rec['num_experiments']}")
        else:
            print("\n⚠️  No recommendations returned (check if data exists)")
    else:
        print(f"Error: {result.get('error')}")
    
    return result.get('success')


def test_reaction_smiles_no_data():
    """Test with furan-diene + ArNH2 (has 0 experiments)"""
    print("\n" + "=" * 80)
    print("Test 2: Reaction SMILES with no data (furan-diene + ArNH2 + Cu)")
    print("=" * 80)
    
    result = chemtools_wrapper.hte_recommend_tool.invoke({
        "reaction_smiles": "Brc1ccco1.Nc1ccccc1>>c1ccccc1Nc1ccco1",
        "reaction_type_filter": "C_N_Coupling",
        "catalyst_filter": "Cu",
        "top_k": 5,
        "min_experiments": 1
    })
    
    print(f"\nSuccess: {result.get('success')}")
    
    if not result.get('success'):
        print(f"Error: {result.get('error')}")
        
        # _error_response merges extra dict directly
        data_source = result
        
        print(f"\nDetected Types:")
        print(f"  Reactant A: {data_source.get('reactant_a_smiles')} → {data_source.get('reactant_a_type')}")
        print(f"  Reactant B: {data_source.get('reactant_b_smiles')} → {data_source.get('reactant_b_type')}")
        print(f"  Matching Experiments: {data_source.get('matching_experiments')}")
        
        if 'suggestion' in data_source:
            print(f"\nSuggestion:")
            print(f"  {data_source['suggestion']}")
        
        if 'filters_applied' in data_source:
            print(f"\nFilters Applied:")
            print(f"  Reaction Type: {data_source['filters_applied'].get('reaction_type')}")
            print(f"  Catalyst: {data_source['filters_applied'].get('catalyst')}")
    else:
        print("Unexpected: Found data when none expected")
    
    return not result.get('success')  # Should fail


def test_individual_reactants_legacy():
    """Test legacy format with individual reactants (backward compatibility)"""
    print("\n" + "=" * 80)
    print("Test 3: Individual reactants (legacy format)")
    print("=" * 80)
    
    result = chemtools_wrapper.hte_recommend_tool.invoke({
        "reactant_a_smiles": "Brc1ccccc1",
        "reactant_b_smiles": "Nc1ccccc1",
        "catalyst_filter": "Cu",
        "top_k": 3,
        "min_experiments": 1
    })
    
    print(f"\nSuccess: {result.get('success')}")
    
    if result.get('success'):
        # Fields merged directly at top level
        print(f"Matching Experiments: {result.get('matching_experiments')}")
        print(f"Found {len(result.get('recommendations', []))} recommendations")
    else:
        print(f"Error: {result.get('error')}")
    
    return result.get('success')


def main():
    """Run all tests"""
    print("\n🧪 Testing HTE Agent Integration\n")
    
    test1_passed = test_reaction_smiles_with_data()
    test2_passed = test_reaction_smiles_no_data()
    test3_passed = test_individual_reactants_legacy()
    
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    print(f"Test 1 (reaction SMILES with data): {'✅ PASS' if test1_passed else '❌ FAIL'}")
    print(f"Test 2 (reaction SMILES no data):   {'✅ PASS' if test2_passed else '❌ FAIL'}")
    print(f"Test 3 (individual reactants):      {'✅ PASS' if test3_passed else '❌ FAIL'}")
    
    all_passed = test1_passed and test2_passed and test3_passed
    print(f"\nOverall: {'✅ ALL TESTS PASSED' if all_passed else '❌ SOME TESTS FAILED'}")
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
