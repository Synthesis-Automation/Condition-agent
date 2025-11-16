"""
Test script for HTE screening set generation.
"""

from chemtools.HTE import HTERecommender


def test_screening_set():
    """Test screening set generation with different strategies."""
    
    recommender = HTERecommender()
    
    # Test reaction: 4-bromotoluene + aniline (C-N coupling with Cu)
    reactant_a = "Brc1ccc(C)cc1"  # 4-bromotoluene
    reactant_b = "Nc1ccccc1"       # aniline
    
    print("="*80)
    print("HTE SCREENING SET GENERATION TEST")
    print("="*80)
    print(f"Reactant A: {reactant_a} (4-bromotoluene)")
    print(f"Reactant B: {reactant_b} (aniline)")
    print(f"Reaction Type: C-N Coupling")
    print(f"Catalyst Filter: Cu")
    print()
    
    # Test 1: Balanced strategy (default) - 24 conditions
    print("-" * 80)
    print("TEST 1: Balanced Strategy (24 conditions)")
    print("-" * 80)
    
    result_balanced = recommender.generate_screening_set(
        reactant_a_smiles=reactant_a,
        reactant_b_smiles=reactant_b,
        num_conditions=24,
        reaction_type_filter="C_N_Coupling",
        catalyst_filter="Cu",
        diversity_strategy="balanced"
    )
    
    print(f"✓ Found {result_balanced.total_matching_experiments} matching experiments")
    print(f"✓ Generated {len(result_balanced.recommendations)} screening conditions")
    print(f"✓ Reactant types: {result_balanced.reactant_a_type} + {result_balanced.reactant_b_type}")
    print(f"✓ Predicted reaction: {result_balanced.predicted_reaction_type} ({result_balanced.reaction_type_confidence*100:.1f}% confidence)")
    print()
    
    # Show first 3 and last 3 conditions
    print("First 3 conditions (top performers):")
    for i, rec in enumerate(result_balanced.recommendations[:3], 1):
        print(f"  {i}. {rec.catalyst}/{rec.ligand}/{rec.base}/{rec.solvent}")
        print(f"     Z-score: {rec.avg_z_score:.2f}, Yield: {rec.avg_yield:.1f}%, Confidence: {rec.confidence_score:.1f}")
    
    print()
    print("Last 3 conditions (diverse picks):")
    for i, rec in enumerate(result_balanced.recommendations[-3:], len(result_balanced.recommendations)-2):
        print(f"  {i}. {rec.catalyst}/{rec.ligand}/{rec.base}/{rec.solvent}")
        print(f"     Z-score: {rec.avg_z_score:.2f}, Yield: {rec.avg_yield:.1f}%, Confidence: {rec.confidence_score:.1f}")
    
    print()
    
    # Test 2: Top performers strategy - 12 conditions
    print("-" * 80)
    print("TEST 2: Top Performers Strategy (12 conditions)")
    print("-" * 80)
    
    result_top = recommender.generate_screening_set(
        reactant_a_smiles=reactant_a,
        reactant_b_smiles=reactant_b,
        num_conditions=12,
        reaction_type_filter="C_N_Coupling",
        catalyst_filter="Cu",
        diversity_strategy="top_performers"
    )
    
    print(f"✓ Generated {len(result_top.recommendations)} conditions (top performers)")
    print()
    print("All conditions:")
    for i, rec in enumerate(result_top.recommendations, 1):
        print(f"  {i}. {rec.catalyst}/{rec.ligand}/{rec.base}/{rec.solvent} "
              f"(Z={rec.avg_z_score:.2f}, Y={rec.avg_yield:.1f}%)")
    
    print()
    
    # Test 3: Diverse strategy - 24 conditions
    print("-" * 80)
    print("TEST 3: Diverse Strategy (24 conditions)")
    print("-" * 80)
    
    result_diverse = recommender.generate_screening_set(
        reactant_a_smiles=reactant_a,
        reactant_b_smiles=reactant_b,
        num_conditions=24,
        reaction_type_filter="C_N_Coupling",
        catalyst_filter="Cu",
        diversity_strategy="diverse"
    )
    
    print(f"✓ Generated {len(result_diverse.recommendations)} conditions (maximizing diversity)")
    
    # Analyze reagent diversity
    catalysts = set(rec.catalyst for rec in result_diverse.recommendations)
    ligands = set(rec.ligand for rec in result_diverse.recommendations)
    bases = set(rec.base for rec in result_diverse.recommendations)
    solvents = set(rec.solvent for rec in result_diverse.recommendations)
    
    print(f"\nReagent diversity:")
    print(f"  Unique catalysts: {len(catalysts)}")
    print(f"  Unique ligands: {len(ligands)}")
    print(f"  Unique bases: {len(bases)}")
    print(f"  Unique solvents: {len(solvents)}")
    
    print()
    
    # Test 4: Agent tool integration
    print("-" * 80)
    print("TEST 4: Agent Tool Integration")
    print("-" * 80)
    
    from chem_assistant.chemtools_wrapper import hte_screening_set_tool
    
    tool_result = hte_screening_set_tool.invoke({
        "reaction_smiles": "Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1",
        "num_conditions": 24,
        "reaction_type_filter": "C_N_Coupling",
        "catalyst_filter": "Cu",
        "diversity_strategy": "balanced"
    })
    
    if tool_result['success']:
        print(f"✓ Tool returned success=True")
        print(f"✓ Reactant types: {tool_result['reactant_a_type']} + {tool_result['reactant_b_type']}")
        print(f"✓ Matching experiments: {tool_result['matching_experiments']}")
        print(f"✓ Generated conditions: {tool_result['num_conditions']}")
        print(f"✓ Plate format: {tool_result['plate_format']}")
        print(f"✓ Diversity strategy: {tool_result['diversity_strategy']}")
        
        # Show plate positions for first 6 conditions
        print(f"\nFirst 6 plate positions:")
        for cond in tool_result['screening_conditions'][:6]:
            print(f"  {cond['plate_position']}: {cond['catalyst']}/{cond['ligand']} "
                  f"(Z={cond['avg_z_score']}, Y={cond['avg_yield']}%)")
    else:
        print(f"✗ Tool returned error: {tool_result.get('error', 'Unknown error')}")
    
    print()
    print("="*80)
    print("ALL TESTS COMPLETED")
    print("="*80)


if __name__ == "__main__":
    test_screening_set()
