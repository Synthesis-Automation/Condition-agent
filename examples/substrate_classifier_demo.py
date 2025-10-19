"""
Demo/Example usage of SubstrateClassifier

Shows how to use the SubstrateClassifier for various chemical substrates.
"""

from chemtools.util.substrate_classifier import (
    SubstrateClassifier,
    classify_substrate,
    get_substrate_class,
    get_substrate_family,
)


def demo_basic_usage():
    """Basic usage examples"""
    print("=" * 70)
    print("  SubstrateClassifier - Basic Usage")
    print("=" * 70)
    
    # Example 1: Primary alkyl iodide
    print("\n1. Primary Alkyl Iodide:")
    print("   SMILES: CCCCCCI")
    info = classify_substrate("CCCCCCI")
    print(f"   Family: {info.substrate_family}")
    print(f"   Class: {info.substrate_class}")
    print(f"   Functional Groups: {info.functional_groups}")
    
    # Example 2: Aryl bromide
    print("\n2. Aryl Bromide:")
    print("   SMILES: c1ccccc1Br")
    info = classify_substrate("c1ccccc1Br")
    print(f"   Family: {info.substrate_family}")
    print(f"   Class: {info.substrate_class}")
    print(f"   Has Aromatic: {info.has_aromatic}")
    
    # Example 3: Aniline
    print("\n3. Aniline:")
    print("   SMILES: c1ccccc1N")
    info = classify_substrate("c1ccccc1N")
    print(f"   Family: {info.substrate_family}")
    print(f"   Class: {info.substrate_class}")
    
    # Example 4: Carboxylic acid
    print("\n4. Carboxylic Acid:")
    print("   SMILES: CC(=O)O")
    info = classify_substrate("CC(=O)O")
    print(f"   Family: {info.substrate_family}")
    print(f"   Class: {info.substrate_class}")
    
    # Example 5: Boronic acid
    print("\n5. Boronic Acid:")
    print("   SMILES: c1ccccc1B(O)O")
    info = classify_substrate("c1ccccc1B(O)O")
    print(f"   Family: {info.substrate_family}")
    print(f"   Class: {info.substrate_class}")


def demo_special_positions():
    """Demo special position detection"""
    print("\n" + "=" * 70)
    print("  Special Position Detection")
    print("=" * 70)
    
    from chemtools.util.rdkit_helpers import rdkit_available
    
    if not rdkit_available():
        print("\n⚠️  RDKit not available - skipping special position demo")
        return
    
    # Benzylic chloride
    print("\n1. Benzylic Chloride:")
    print("   SMILES: c1ccccc1CCl")
    info = classify_substrate("c1ccccc1CCl")
    print(f"   Class: {info.substrate_class}")
    print(f"   Benzylic positions: {info.special_positions.benzylic}")
    
    # Allylic bromide
    print("\n2. Allylic Bromide:")
    print("   SMILES: C=CCBr")
    info = classify_substrate("C=CCBr")
    print(f"   Class: {info.substrate_class}")
    print(f"   Allylic positions: {info.special_positions.allylic}")
    
    # Propargylic iodide
    print("\n3. Propargylic Iodide:")
    print("   SMILES: C#CCI")
    info = classify_substrate("C#CCI")
    print(f"   Class: {info.substrate_class}")
    print(f"   Propargylic positions: {info.special_positions.propargylic}")


def demo_carbon_types():
    """Demo carbon hybridization detection"""
    print("\n" + "=" * 70)
    print("  Carbon Hybridization Analysis")
    print("=" * 70)
    
    from chemtools.util.rdkit_helpers import rdkit_available
    
    if not rdkit_available():
        print("\n⚠️  RDKit not available - skipping carbon type demo")
        return
    
    # Molecule with mixed hybridization
    print("\n1. Mixed Hybridization (C#C-C=C-C-c):")
    print("   SMILES: C#CC=CCc1ccccc1")
    info = classify_substrate("C#CC=CCc1ccccc1")
    print(f"   Carbon types: {info.carbon_types}")
    
    # Count by type
    sp_count = sum(1 for t in info.carbon_types.values() if t == 'sp')
    sp2_count = sum(1 for t in info.carbon_types.values() if t == 'sp2')
    sp3_count = sum(1 for t in info.carbon_types.values() if t == 'sp3')
    aromatic_count = sum(1 for t in info.carbon_types.values() if t == 'aromatic')
    
    print(f"   sp: {sp_count}, sp2: {sp2_count}, sp3: {sp3_count}, aromatic: {aromatic_count}")


def demo_convenience_functions():
    """Demo convenience functions"""
    print("\n" + "=" * 70)
    print("  Convenience Functions")
    print("=" * 70)
    
    substrates = [
        ("CCCCI", "Primary alkyl iodide"),
        ("c1ccccc1Br", "Aryl bromide"),
        ("c1ccccc1N", "Aniline"),
        ("CCCO", "Aliphatic alcohol"),
        ("c1ccccc1B(O)O", "Boronic acid"),
    ]
    
    print("\n| SMILES | Description | Family | Class |")
    print("|" + "-" * 68 + "|")
    
    for smiles, desc in substrates:
        family = get_substrate_family(smiles)
        cls = get_substrate_class(smiles)
        print(f"| {smiles:15} | {desc:20} | {family:10} | {cls:20} |")


def demo_real_world_examples():
    """Demo with real-world reaction substrates"""
    print("\n" + "=" * 70)
    print("  Real-World Cross-Coupling Substrates")
    print("=" * 70)
    
    examples = [
        # Alkyl halides
        ("CI", "Methyl iodide"),
        ("CCCCCCCCI", "Octyl iodide"),
        ("CC(C)Br", "Isopropyl bromide (2°)"),
        ("CC(C)(C)Cl", "tert-Butyl chloride (3°)"),
        
        # Aryl halides
        ("c1ccccc1I", "Phenyl iodide"),
        ("c1ccc(Br)cc1C", "4-Bromotoluene"),
        ("c1cnc(Cl)cc1", "2-Chloropyridine"),
        
        # Special positions
        ("c1ccccc1CI", "Benzyl iodide (benzylic)"),
        ("C=CCBr", "Allyl bromide (allylic)"),
        
        # Amines
        ("c1ccccc1N", "Aniline"),
        ("CCN", "Ethylamine"),
        
        # Boron compounds
        ("c1ccccc1B(O)O", "Phenylboronic acid"),
        ("c1ccccc1B1OC(C)(C)C(C)(C)O1", "Phenylboronic acid pinacol ester"),
    ]
    
    for smiles, desc in examples:
        info = classify_substrate(smiles)
        print(f"\n{desc}:")
        print(f"  SMILES: {smiles}")
        print(f"  Family: {info.substrate_family}")
        print(f"  Class:  {info.substrate_class}")


def demo_reusability():
    """Demo how other modules can use SubstrateClassifier"""
    print("\n" + "=" * 70)
    print("  Reusability Example: ML Feature Extraction")
    print("=" * 70)
    
    from chemtools.util.rdkit_helpers import rdkit_available
    
    def extract_substrate_features(smiles: str):
        """Example: Extract ML features from substrate"""
        info = classify_substrate(smiles)
        
        features = {
            'substrate_class': info.substrate_class,
            'substrate_family': info.substrate_family,
            'has_aromatic': info.has_aromatic,
            'has_heteroaromatic': info.has_heteroaromatic,
            'ring_count': info.ring_count,
            'aromatic_ring_count': info.aromatic_ring_count,
        }
        
        if rdkit_available():
            features.update({
                'has_benzylic': len(info.special_positions.benzylic) > 0,
                'has_allylic': len(info.special_positions.allylic) > 0,
                'sp3_carbon_count': sum(1 for t in info.carbon_types.values() if t == 'sp3'),
                'aromatic_carbon_count': sum(1 for t in info.carbon_types.values() if t == 'aromatic'),
            })
        
        return features
    
    # Extract features for a substrate
    print("\nExtracting features for: c1ccccc1CI (benzyl iodide)")
    features = extract_substrate_features("c1ccccc1CI")
    
    for key, value in features.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    demo_basic_usage()
    demo_special_positions()
    demo_carbon_types()
    demo_convenience_functions()
    demo_real_world_examples()
    demo_reusability()
    
    print("\n" + "=" * 70)
    print("  Demo Complete!")
    print("=" * 70)
    print("\nFor more examples, see:")
    print("  - tests/test_substrate_classifier.py")
    print("  - docs/SMARTS_REFACTORED_ARCHITECTURE.md")
