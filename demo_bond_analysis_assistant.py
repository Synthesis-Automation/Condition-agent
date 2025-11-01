#!/usr/bin/env python3
"""
Demonstration of the ChemTools Assistant using the Bond Analysis Tool.

This shows how the assistant can now analyze bond breaking and formation
in chemical reactions using the newly integrated analyze_bond_changes_tool.
"""

import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chem_assistant.chemtools_wrapper import analyze_bond_changes_tool
import json

def format_result(result):
    """Format the bond analysis result for display."""
    if not result.get('success'):
        return f"❌ Error: {result.get('error', 'Unknown error')}"
    
    output = []
    output.append(f"✅ Analysis Method: {result.get('method', 'unknown').upper()}")
    output.append(f"   Confidence: {result.get('combined_confidence', 0):.2%}")
    output.append(f"   Validation: {result.get('validation', 'N/A')}")
    output.append("")
    
    # Broken bonds
    broken = result.get('broken_bonds', [])
    if broken:
        output.append(f"Bonds Broken ({len(broken)}):")
        for bond in broken:
            if isinstance(bond, list) and len(bond) == 2:
                output.append(f"  ❌ {bond[0]} — {bond[1]}")
            else:
                output.append(f"  ❌ {bond}")
    else:
        output.append("Bonds Broken: None detected")
    output.append("")
    
    # Formed bonds
    formed = result.get('formed_bonds', [])
    if formed:
        # Remove duplicates
        unique = set()
        for bond in formed:
            if isinstance(bond, list) and len(bond) == 2:
                unique.add(tuple(sorted(bond)))
        output.append(f"Bonds Formed ({len(unique)}):")
        for bond in unique:
            output.append(f"  ✅ {bond[0]} — {bond[1]}")
    else:
        output.append("Bonds Formed: None detected")
    output.append("")
    
    # Interpretation
    if 'interpretation' in result:
        output.append(f"Interpretation: {result['interpretation']}")
    
    # Method-specific info
    if result.get('rxnmapper_confidence'):
        output.append(f"RXNMapper confidence: {result['rxnmapper_confidence']:.2%}")
    if result.get('mcs_coverage'):
        output.append(f"MCS coverage: {result['mcs_coverage']:.2%}")
    
    return '\n'.join(output)


def demo_reactions():
    """Demonstrate bond analysis on various reactions."""
    
    reactions = [
        {
            "name": "Suzuki-Miyaura Cross-Coupling",
            "smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
            "description": "Classic C-C bond forming reaction with Br and B leaving"
        },
        {
            "name": "Buchwald-Hartwig C-N Coupling",
            "smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
            "description": "C-N bond formation with Br leaving group"
        },
        {
            "name": "Esterification",
            "smiles": "CC(=O)O.CCO>>CC(=O)OCC",
            "description": "Simple ester formation"
        },
        {
            "name": "Sonogashira Coupling",
            "smiles": "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",
            "description": "C-C triple bond extension"
        }
    ]
    
    print("=" * 80)
    print("CHEMTOOLS ASSISTANT - BOND ANALYSIS DEMONSTRATION")
    print("=" * 80)
    print()
    
    for i, rxn in enumerate(reactions, 1):
        print(f"\n{'=' * 80}")
        print(f"Reaction {i}: {rxn['name']}")
        print('=' * 80)
        print(f"Description: {rxn['description']}")
        print(f"SMILES: {rxn['smiles']}")
        print("-" * 80)
        
        # Analyze bonds
        result = analyze_bond_changes_tool.invoke({
            "reaction_smiles": rxn['smiles'],
            "use_hybrid": True
        })
        
        # Display results
        print(format_result(result))
        print()
    
    print("=" * 80)
    print("✅ DEMONSTRATION COMPLETE")
    print("=" * 80)
    print()
    print("The assistant can now:")
    print("  • Detect bond breaking and formation")
    print("  • Identify leaving groups (Br, I, B, etc.)")
    print("  • Validate results across multiple methods")
    print("  • Provide confidence scores")
    print("  • Interpret reaction mechanisms")
    print()


if __name__ == "__main__":
    demo_reactions()
