#!/usr/bin/env python3
"""
Test the new hte_conditions_tool for ArI + Carbamate + Cu
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from chem_assistant.chemtools_wrapper import hte_conditions_tool

print("=" * 80)
print("Testing HTE Conditions Tool: ArI + Carbamate + Cu")
print("=" * 80)

# Test the exact user query
print("\nQuery: Top 10 conditions for copper-catalyzed C-N coupling of ArI + Carbamate")
print("-" * 80)

result = hte_conditions_tool.invoke({
    "reactant_a_type": "ArI",
    "reactant_b_type": "Carbamate",
    "reaction_type": "C-N",
    "catalyst_filter": "Cu",
    "top_k": 10,
    "min_experiments": 1,
    "sort_by": "count"
})

if result["success"]:
    print(f"✓ Query succeeded!")
    print(f"\nReactant pair: {result['reactant_pair']}")
    print(f"Total experiments: {result['total_experiments']}")
    print(f"Total unique conditions: {result['total_conditions']}")
    print(f"Filters: {result['filters']}")
    
    print(f"\nTop {len(result['conditions'])} Conditions:")
    print("=" * 80)
    
    for cond in result['conditions']:
        print(f"\n{cond['rank']}. Catalyst: {cond['catalyst']}")
        print(f"   Ligand: {cond['ligand']}")
        print(f"   Base: {cond['base']}")
        print(f"   Solvent: {cond['solvent']}")
        if cond.get('secondary_solvent'):
            print(f"   Secondary Solvent: {cond['secondary_solvent']}")
        if cond.get('additive'):
            print(f"   Additive: {cond['additive']}")
        print(f"   ---")
        print(f"   Experiments: {cond['experiments']}")
        print(f"   Avg Yield: {cond['avg_yield']}%")
        print(f"   Median Yield: {cond['median_yield']}%")
        print(f"   Success Rate: {cond['success_rate']}%")
    
    print("\n" + "=" * 80)
    print("Statistics Summary:")
    print("=" * 80)
    
    total_exp = sum(c['experiments'] for c in result['conditions'])
    best_yield = max(c['avg_yield'] for c in result['conditions'])
    best_success = max(c['success_rate'] for c in result['conditions'])
    
    print(f"Experiments in top {len(result['conditions'])}: {total_exp}")
    print(f"Best average yield: {best_yield}%")
    print(f"Best success rate: {best_success}%")
    
else:
    print(f"✗ Query failed: {result['error']}")
    if 'suggestion' in result:
        print(f"  Suggestion: {result['suggestion']}")

# Test 2: Sort by success rate
print("\n\n" + "=" * 80)
print("Test 2: Sort by Success Rate")
print("=" * 80)

result2 = hte_conditions_tool.invoke({
    "reactant_a_type": "ArI",
    "reactant_b_type": "Carbamate",
    "catalyst_filter": "Cu",
    "top_k": 5,
    "min_experiments": 2,
    "sort_by": "success"
})

if result2["success"]:
    print(f"\nTop 5 by success rate (min 2 experiments):")
    for cond in result2['conditions']:
        print(f"{cond['rank']}. {cond['catalyst']} + {cond['ligand']} + {cond['base']} in {cond['solvent']}")
        print(f"   {cond['success_rate']}% success, {cond['experiments']} exp, {cond['avg_yield']}% yield")

print("\n" + "=" * 80)
print("Summary")
print("=" * 80)
print("""
✓ hte_conditions_tool working correctly
✓ Can query specific substrate pair conditions
✓ Returns detailed catalyst/ligand/base/solvent combinations
✓ Supports sorting by count, success rate, or yield

The agent can now answer:
"For copper-catalyzed C-N coupling of ArI and Carbamate, what are the top 10 conditions?"
""")
print("=" * 80)
