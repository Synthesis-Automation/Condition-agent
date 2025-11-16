"""
Test agent with HTE screening set generation.
"""

from chem_assistant.chemtools_agent import ChemToolsAgent


def test_agent_screening():
    """Test agent with screening set request."""
    
    agent = ChemToolsAgent(verbose=True)
    
    query = """
    I need to set up an HTE screening plate for C-N coupling.
    The reaction is: Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1
    (4-bromotoluene + aniline)
    
    Please generate 24 diverse conditions using copper catalyst for a 4x6 plate.
    Use balanced strategy.
    """
    
    print("="*80)
    print("AGENT TEST: HTE SCREENING SET GENERATION")
    print("="*80)
    print(f"\nQuery:\n{query}")
    print("\n" + "="*80)
    print("Agent Response:")
    print("="*80 + "\n")
    
    response = agent.run(query)
    
    print(response)
    print("\n" + "="*80)


if __name__ == "__main__":
    test_agent_screening()
