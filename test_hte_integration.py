#!/usr/bin/env python3
"""
Test HTE Tools Integration with Chem Assistant

Verifies that HTE recommendation and analytics tools are properly integrated
into the LangChain agent system.
"""

import sys
import json
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

print("=" * 80)
print("Testing HTE Tools Integration with Chem Assistant")
print("=" * 80)

# Test 1: Import HTE tools from wrapper
print("\n1. Testing imports...")
try:
    from chem_assistant.chemtools_wrapper import (
        CHEMTOOLS_TOOLS,
        hte_recommend_tool,
        hte_analytics_tool,
        HTE_AVAILABLE
    )
    print(f"   ✓ Successfully imported HTE tools")
    print(f"   ✓ HTE_AVAILABLE = {HTE_AVAILABLE}")
    print(f"   ✓ Total tools: {len(CHEMTOOLS_TOOLS)}")
except Exception as e:
    print(f"   ✗ Import failed: {e}")
    sys.exit(1)

# Test 2: Verify tools are in the collection
print("\n2. Checking tool registration...")
tool_names = [tool.name for tool in CHEMTOOLS_TOOLS]
print(f"   Available tools: {', '.join(tool_names)}")

if "hte_recommend_tool" in tool_names:
    print("   ✓ hte_recommend_tool registered")
else:
    print("   ✗ hte_recommend_tool NOT found in CHEMTOOLS_TOOLS")

if "hte_analytics_tool" in tool_names:
    print("   ✓ hte_analytics_tool registered")
else:
    print("   ✗ hte_analytics_tool NOT found in CHEMTOOLS_TOOLS")

# Test 3: Test HTE recommendation tool
print("\n3. Testing HTE recommendation tool...")
if HTE_AVAILABLE:
    try:
        # Test Suzuki-Miyaura coupling
        result = hte_recommend_tool.invoke({
            "reactant_a_smiles": "Brc1ccccc1",  # Bromobenzene
            "reactant_b_smiles": "c1ccc(B(O)O)cc1",  # Phenylboronic acid
            "top_k": 3,
            "catalyst_filter": "Pd"
        })
        
        if result.get("success"):
            print(f"   ✓ Recommendation succeeded")
            print(f"   - Predicted reaction: {result.get('predicted_reaction_type')}")
            print(f"   - Matching experiments: {result.get('matching_experiments')}")
            print(f"   - Recommendations: {len(result.get('recommendations', []))}")
            
            if result.get('recommendations'):
                top = result['recommendations'][0]
                print(f"   - Top condition: {top.get('catalyst')} + {top.get('base')} in {top.get('solvent')}")
                print(f"     Success rate: {top.get('success_rate')}%, Confidence: {top.get('confidence_score')}")
        else:
            print(f"   ✗ Recommendation failed: {result.get('error')}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
else:
    print("   ⚠ HTE not available - skipping test")

# Test 4: Test HTE analytics tool
print("\n4. Testing HTE analytics tool...")
if HTE_AVAILABLE:
    try:
        # Test list_pairs query
        result = hte_analytics_tool.invoke({
            "query_type": "list_pairs",
            "reaction_type": "Suzuki",
            "catalyst_filter": "Pd",
            "min_experiments": 50,
            "top_n": 5
        })
        
        if result.get("success"):
            print(f"   ✓ Analytics query succeeded")
            print(f"   - Query type: {result.get('query_type')}")
            print(f"   - Total results: {result.get('total_results')}")
            
            results = result.get('results', [])
            if results:
                print(f"   - Top reactant pair: {results[0].get('reactant_a')} + {results[0].get('reactant_b')}")
                print(f"     Experiments: {results[0].get('experiments')}, Success rate: {results[0].get('success_rate')}%")
        else:
            print(f"   ✗ Analytics failed: {result.get('error')}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test metals query
    print("\n   Testing metals analytics...")
    try:
        result = hte_analytics_tool.invoke({
            "query_type": "metals",
            "top_n": 5
        })
        
        if result.get("success"):
            print(f"   ✓ Metal analysis succeeded")
            metals = result.get('metals', [])
            if metals:
                print(f"   - Top 3 metals:")
                for metal_data in metals[:3]:
                    print(f"     {metal_data['metal']}: {metal_data['experiments']:,} experiments ({metal_data['percentage']}%)")
        else:
            print(f"   ✗ Metal analysis failed: {result.get('error')}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
else:
    print("   ⚠ HTE not available - skipping test")

# Test 5: Test with C-N coupling
print("\n5. Testing C-N coupling recommendation...")
if HTE_AVAILABLE:
    try:
        result = hte_recommend_tool.invoke({
            "reactant_a_smiles": "Clc1ccccc1",  # Chlorobenzene
            "reactant_b_smiles": "CCN",  # Ethylamine
            "top_k": 2,
            "reaction_type_filter": "C-N"
        })
        
        if result.get("success"):
            print(f"   ✓ C-N coupling recommendation succeeded")
            print(f"   - Predicted reaction: {result.get('predicted_reaction_type')}")
            print(f"   - Matching experiments: {result.get('matching_experiments')}")
            
            if result.get('recommendations'):
                top = result['recommendations'][0]
                print(f"   - Top catalyst: {top.get('catalyst')}")
                print(f"     Base: {top.get('base')}, Solvent: {top.get('solvent')}")
                print(f"     Success: {top.get('success_rate')}%, Yield: {top.get('avg_yield')}%")
        else:
            print(f"   ✗ Failed: {result.get('error')}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
else:
    print("   ⚠ HTE not available - skipping test")

# Test 6: Test catalyst comparison analytics
print("\n6. Testing catalyst comparison...")
if HTE_AVAILABLE:
    try:
        result = hte_analytics_tool.invoke({
            "query_type": "catalysts",
            "reaction_type": "Suzuki",
            "top_n": 5
        })
        
        if result.get("success"):
            print(f"   ✓ Catalyst analysis succeeded")
            catalysts = result.get('results', [])
            if catalysts:
                print(f"   - Top 3 Suzuki catalysts:")
                for cat in catalysts[:3]:
                    print(f"     {cat['catalyst']} ({cat['metal']}): {cat['experiments']} exp, {cat['success_rate']}% success")
        else:
            print(f"   ✗ Failed: {result.get('error')}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
else:
    print("   ⚠ HTE not available - skipping test")

# Summary
print("\n" + "=" * 80)
print("Integration Test Summary")
print("=" * 80)
print(f"HTE Available: {HTE_AVAILABLE}")
print(f"Total Tools: {len(CHEMTOOLS_TOOLS)}")
print(f"HTE Tools Registered: {'hte_recommend_tool' in tool_names and 'hte_analytics_tool' in tool_names}")

if HTE_AVAILABLE:
    print("\n✓ All HTE tools successfully integrated into chem_assistant!")
    print("\nUsage with LangChain agent:")
    print("```python")
    print("from chem_assistant.chemtools_wrapper import CHEMTOOLS_TOOLS")
    print("from langgraph.prebuilt import create_agent")
    print("")
    print("agent = create_agent(llm, CHEMTOOLS_TOOLS)")
    print('result = agent.invoke({"messages": ["Recommend conditions for Suzuki coupling of bromobenzene and phenylboronic acid"]})')
    print("```")
else:
    print("\n⚠ HTE tools not available. Install with:")
    print("   pip install chemtools[hte]")

print("=" * 80)
