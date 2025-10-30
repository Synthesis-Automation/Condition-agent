"""
Demo: Reaction-Specific Reactant Analysis

Shows how to analyze reactants specific to a detected reaction type.
For example, for Suzuki reactions, analyze the aryl halide and boron reagent.
"""

from chemtools.detection import detect_reaction
from chemtools.analysis.reactants import (
    classify_reactant_smiles,
    required_reactant_categories,
    iter_reactant_matches
)

def analyze_suzuki_reactants(reaction_smiles: str):
    """Demonstrate Suzuki-specific reactant analysis"""
    print("=" * 80)
    print("SUZUKI REACTION ANALYSIS")
    print("=" * 80)
    print(f"\nReaction: {reaction_smiles}\n")
    
    # 1. Detect reaction type
    result = detect_reaction(reaction_smiles, use_ml=False)
    family = result.get('family')
    print(f"Detected Family: {family}")
    print(f"Confidence: {result.get('confidence', 0):.2f}\n")
    
    # 2. Get expected reactant categories for this reaction type
    # Map internal family name to taxonomy reaction id
    reaction_id_map = {
        'suzuki_miyaura': 'Suzuki-Miyaura',
        'buchwald_hartwig_c_n': 'Buchwald-Hartwig',
        'ullmann_cn': 'Ullmann-CN',
    }
    
    reaction_id = reaction_id_map.get(family)
    if reaction_id:
        expected = required_reactant_categories(reaction_id)
        print(f"Expected Reactant Types for {reaction_id}:")
        if expected:
            for i, group in enumerate(expected, 1):
                print(f"  Reactant {i}: {', '.join(group)}")
        print()
    
    # 3. Analyze each reactant
    details = result.get('details', {})
    reactants = details.get('reactants', [])
    
    print("Reactant Classification:")
    print("-" * 80)
    
    for i, reactant_smiles in enumerate(reactants, 1):
        print(f"\nReactant {i}: {reactant_smiles}")
        
        # Get best match
        match = classify_reactant_smiles(reactant_smiles)
        if match:
            print(f"  Category: {match.category}")
            print(f"  Type: {match.member_type}")
            print(f"  Name: {match.name}")
            print(f"  Group: {match.group}")
            print(f"  SMARTS: {match.smarts}")
            
            # Get all matches to see what else it could be
            all_matches = iter_reactant_matches(reactant_smiles)
            if len(all_matches) > 1:
                print(f"  Also matches: {', '.join(set(m.category for m in all_matches[1:4]))}")
        else:
            print("  ⚠ No classification match")
    
    # 4. Functional group detection
    fg = details.get('functional_groups', {})
    detected_groups = [k for k, v in fg.items() if v]
    print(f"\n\nDetected Functional Groups:")
    print("-" * 80)
    for group in detected_groups:
        print(f"  ✓ {group}")
    
    print("\n" + "=" * 80 + "\n")


def analyze_buchwald_hartwig_reactants(reaction_smiles: str):
    """Demonstrate Buchwald-Hartwig-specific reactant analysis"""
    print("=" * 80)
    print("BUCHWALD-HARTWIG REACTION ANALYSIS")
    print("=" * 80)
    print(f"\nReaction: {reaction_smiles}\n")
    
    # Detect reaction
    result = detect_reaction(reaction_smiles, use_ml=False)
    family = result.get('family')
    print(f"Detected Family: {family}")
    print(f"Confidence: {result.get('confidence', 0):.2f}\n")
    
    # Analyze reactants
    details = result.get('details', {})
    reactants = details.get('reactants', [])
    
    print("Reactant Classification:")
    print("-" * 80)
    
    aryl_halide = None
    amine = None
    
    for reactant_smiles in reactants:
        match = classify_reactant_smiles(reactant_smiles)
        if match:
            print(f"\n{reactant_smiles}")
            print(f"  → {match.category}: {match.name}")
            
            if 'ArBr' in match.category or 'ArCl' in match.category or 'ArI' in match.category:
                aryl_halide = reactant_smiles
                print(f"  → IDENTIFIED AS: Aryl Halide (electrophile)")
            elif 'NH' in match.category or 'amine' in match.group.lower():
                amine = reactant_smiles
                print(f"  → IDENTIFIED AS: Amine (nucleophile)")
    
    # Summary
    print("\n\nReaction Partner Summary:")
    print("-" * 80)
    if aryl_halide:
        print(f"  Aryl Halide: {aryl_halide}")
    if amine:
        print(f"  Amine: {amine}")
    
    fg = details.get('functional_groups', {})
    print(f"\n  Key functional groups:")
    if fg.get('aryl_halide'):
        print(f"    ✓ Aryl halide present")
    if fg.get('nucleophile_n'):
        print(f"    ✓ Nitrogen nucleophile present")
    
    print("\n" + "=" * 80 + "\n")


if __name__ == "__main__":
    # Test 1: Suzuki reaction
    suzuki_rxn = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    analyze_suzuki_reactants(suzuki_rxn)
    
    # Test 2: Buchwald-Hartwig reaction  
    bh_rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    analyze_buchwald_hartwig_reactants(bh_rxn)
    
    # Test 3: Different Suzuki with more complex reactants
    print("\n" + "=" * 80)
    print("COMPLEX SUZUKI EXAMPLE")
    print("=" * 80)
    complex_suzuki = "Clc1ccc(C#N)cc1.COc1ccc(B(O)O)cc1>>COc1ccc(-c2ccc(C#N)cc2)cc1"
    analyze_suzuki_reactants(complex_suzuki)
