"""
ChemTools Basic Tools Demo
============================
Demonstrates core SMILES, router, and featurization tools.

Quick Start:
    python tests/demo_basic_tools.py
    # OR from tests directory:
    cd tests && python demo_basic_tools.py

Features Demonstrated:
    1. SMILES normalization (canonicalization)
    2. Reaction family detection (with and without catalysts)
    3. Molecular featurization (leaving groups, nucleophiles)
    4. Property lookup
    5. DRFP similarity calculations
"""

import sys
from pathlib import Path

# Add parent directory to path for chemtools import
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chemtools.smiles import normalize, normalize_reaction
from chemtools import detect_reaction
from chemtools.featurizers.molecular import featurize
from chemtools.reaction_similarity import drfp_tanimoto, drfp_available
from chemtools.reagent import find_reagent


def test_1_smiles_normalization():
    """Demonstrate SMILES canonicalization."""
    print("\n" + "="*70)
    print("  1. SMILES Normalization")
    print("="*70)
    
    # Simple molecule
    result = normalize('c1ccccc1O')
    print(f"✅ normalize('c1ccccc1O')")
    print(f"   → {result}")
    
    # Reaction normalization
    rxn = 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'
    result = normalize_reaction(rxn)
    print(f"\n✅ normalize_reaction('{rxn}')")
    print(f"   → Normalized: {result.get('normalized', 'N/A')}")
    print(f"   → Reactants: {len(result.get('reactants', []))}")
    print(f"   → Products: {len(result.get('products', []))}")
    
    # Complex reaction with agents
    rxn_with_agents = 'Brc1ccccc1.Nc1ccccc1>Pd(PPh3)4.K2CO3>c1ccccc1Nc1ccccc1'
    result = normalize_reaction(rxn_with_agents)
    print(f"\n✅ normalize_reaction (with agents)")
    print(f"   → Input: {rxn_with_agents}")
    print(f"   → Normalized: {result.get('normalized', 'N/A')}")
    print(f"   → Agents: {len(result.get('agents', []))}")


def test_2_family_detection():
    """Demonstrate reaction family detection."""
    print("\n" + "="*70)
    print("  2. Reaction Family Detection")
    print("="*70)
    
    # Method 1: From reactant list (rule-based only)
    print("Method 1: detect_reaction(reactants) - Rule-based")
    print("-" * 70)
    reactants = ['Brc1ccccc1', 'Nc1ccccc1']
    # Convert reactants to pseudo-reaction
    pseudo_rxn = ".".join(reactants) + ">>"
    result = detect_reaction(pseudo_rxn, use_ml=False)
    print(f"✅ detect_reaction(\"{pseudo_rxn}\", use_ml=False)")
    print(f"   → Family: {result.get('family', 'N/A')}")
    print(f"   → Confidence: {result.get('confidence', 0):.2f}")
    
    # Method 2: From reaction SMILES (detects catalysts!)
    print("\nMethod 2: detect_reaction() - With catalyst detection")
    print("-" * 70)
    rxn = 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'
    result = detect_reaction(rxn, use_ml=False)
    print(f"✅ detect_reaction('{rxn}', use_ml=False)")
    print(f"   → Family: {result.get('family', 'N/A')}")
    print(f"   → Confidence: {result.get('confidence', 0):.2f}")
    print(f"   💡 Automatically detects catalysts (Pd/Cu) from agents")
    
    # Example with Pd catalyst (should detect Buchwald)
    rxn_pd = 'Brc1ccccc1.Nc1ccccc1>Pd(PPh3)4>c1ccccc1Nc1ccccc1'
    result = detect_reaction(rxn_pd, use_ml=False)
    print(f"\n✅ With Pd catalyst:")
    print(f"   → Input: {rxn_pd}")
    print(f"   → Family: {result.get('family', 'N/A')}")
    print(f"   → Confidence: {result.get('confidence', 0):.2f}")
    print(f"   💡 Pd catalyst → buchwald_hartwig_c_n (higher confidence)")
    
    # Example with Cu catalyst (should detect Ullmann)
    rxn_cu = 'Brc1ccccc1.Nc1ccccc1>CuI>c1ccccc1Nc1ccccc1'
    result = detect_reaction(rxn_cu, use_ml=False)
    print(f"\n✅ With Cu catalyst:")
    print(f"   → Input: {rxn_cu}")
    print(f"   → Family: {result.get('family', 'N/A')}")
    print(f"   → Confidence: {result.get('confidence', 0):.2f}")
    print(f"   💡 Cu catalyst → ullmann_cn")
    
    # Suzuki reaction
    print("\n✅ Suzuki coupling:")
    suzuki = 'Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1'
    result = detect_reaction(suzuki, use_ml=False)
    print(f"   → Input: {suzuki}")
    print(f"   → Family: {result.get('family', 'N/A')}")
    print(f"   → Confidence: {result.get('confidence', 0):.2f}")



def test_3_molecular_featurization():
    """Demonstrate molecular featurization."""
    print("\n" + "="*70)
    print("  3. Molecular Featurization")
    print("="*70)
    
    # C-N coupling example
    lg_smiles = 'Brc1ccccc1'
    nuc_smiles = 'Nc1ccccc1'
    result = featurize(lg_smiles, nuc_smiles)
    print(f"✅ featurize('{lg_smiles}', '{nuc_smiles}')")
    print(f"   → LG: {result.get('LG', 'N/A')}")
    print(f"   → Electrophile class: {result.get('elec_class', 'N/A')}")
    print(f"   → Nucleophile class: {result.get('nuc_class', 'N/A')}")
    print(f"   → Bin: {result.get('bin', 'N/A')}")
    
    # Different leaving group
    lg_smiles = 'Clc1ccccc1'
    result = featurize(lg_smiles, nuc_smiles)
    print(f"\n✅ featurize('{lg_smiles}', '{nuc_smiles}')")
    print(f"   → LG: {result.get('LG', 'N/A')}")
    print(f"   → Electrophile class: {result.get('elec_class', 'N/A')}")
    print(f"   → Nucleophile class: {result.get('nuc_class', 'N/A')}")
    print(f"   → Bin: {result.get('bin', 'N/A')}")
    
    # Secondary amine
    nuc_smiles = 'CNc1ccccc1'
    result = featurize(lg_smiles, nuc_smiles)
    print(f"\n✅ featurize('{lg_smiles}', '{nuc_smiles}')")
    print(f"   → LG: {result.get('LG', 'N/A')}")
    print(f"   → Electrophile class: {result.get('elec_class', 'N/A')}")
    print(f"   → Nucleophile class: {result.get('nuc_class', 'N/A')}")
    print(f"   → Basicity: {result.get('n_basicity', 'N/A')}")
    print(f"   → Steric (α): {result.get('steric_alpha', 'N/A')}")
    print(f"   → Bin: {result.get('bin', 'N/A')}")


def test_4_property_lookup():
    """Demonstrate property lookup."""
    print("\n" + "="*70)
    print("  4. Property Lookup (Full Reagent Database)")
    print("="*70)
    print("  📚 Using find_reagent() - Full database (5000+ compounds)")
    print("  " + "-"*66)
    
    # Try DMF in the solvent database
    result = find_reagent('DMF', 'solvent')
    if result:
        print(f"  ✅ find_reagent('DMF', 'solvent')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
        print(f"     → SMILES: {result.get('smiles', 'N/A')}")
        props = result.get('properties', {})
        if props:
            bp = props.get('boiling_point')
            if bp:
                print(f"     → Boiling point: {bp}°C")
    
    # Try K3PO4 in base database
    result = find_reagent('K3PO4', 'base')
    if result:
        print(f"\n  ✅ find_reagent('K3PO4', 'base')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
        pka = result.get('pka')
        if pka:
            print(f"     → pKa: {pka}")
    
    # Try a ligand
    result = find_reagent('BINAP', 'ligand')
    if result:
        print(f"\n  ✅ find_reagent('BINAP', 'ligand')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
        abbr = result.get('abbreviation', [])
        if abbr:
            print(f"     → Abbreviations: {', '.join(abbr)}")
    
    # Try a base
    result = find_reagent('Cesium carbonate', 'base')
    if result:
        print(f"\n  ✅ find_reagent('Cesium carbonate', 'base')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
    
    # Try a metal precursor
    result = find_reagent('Pd(OAc)2', 'metal_precursor')
    if result:
        print(f"\n  ✅ find_reagent('Pd(OAc)2', 'metal_precursor')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
    
    print(f"\n  💡 Available databases: ligand, base, solvent, metal_precursor, additive")


def test_5_database_analytics():
    """Demonstrate database analytics."""
    from chemtools.reagent import (
        count_reagents_by_type,
        get_all_reagent_types,
        get_family_statistics,
        find_reagents_by_family,
    )
    
    print("\n" + "="*70)
    print("  5. Database Analytics")
    print("="*70)
    
    # Get all types and counts
    types = get_all_reagent_types()
    print(f"\n  ✅ Available reagent types ({len(types)}):")
    
    for reagent_type in types[:5]:  # Show first 5
        count = count_reagents_by_type(reagent_type)
        print(f"     → {reagent_type:20s}: {count:3d} reagents")
    
    if len(types) > 5:
        print(f"     ... and {len(types) - 5} more types")
    
    # Ligand family statistics
    print(f"\n  ✅ Ligand families:")
    ligand_stats = get_family_statistics('ligand')
    print(f"     → Total ligands: {ligand_stats['total_reagents']}")
    print(f"     → Total families: {ligand_stats['total_families']}")
    
    # Show top 3 families
    for family_data in ligand_stats['families'][:3]:
        name = family_data['name']
        count = family_data['count']
        print(f"     → {name}: {count} ligands")
    
    # Find reagents in a family
    print(f"\n  ✅ Finding phosphine ligands...")
    phosphines = find_reagents_by_family('ligand', 'trialkyl_triaryl_phosphines')
    print(f"     → Found {len(phosphines)} phosphine ligands")
    
    if phosphines:
        ligand = phosphines[0]
        print(f"     → Example: {ligand.get('name', 'N/A')}")
    
    print(f"\n  💡 See scripts/demo_reagent_analytics.py for more examples")


def test_6_drfp_similarity():
    """Demonstrate property lookup."""
    print("\n" + "="*70)
    print("  4. Property Lookup (Full Reagent Database)")
    print("="*70)
    print("  📚 Using find_reagent() - Full database (5000+ compounds)")
    print("  " + "-"*66)
    
    # Try DMF in the solvent database
    result = find_reagent('DMF', 'solvent')
    if result:
        print(f"  ✅ find_reagent('DMF', 'solvent')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
        print(f"     → SMILES: {result.get('smiles', 'N/A')}")
        props = result.get('properties', {})
        if props:
            bp = props.get('boiling_point')
            if bp:
                print(f"     → Boiling point: {bp}°C")
    
    # Try K3PO4 in base database
    result = find_reagent('K3PO4', 'base')
    if result:
        print(f"\n  ✅ find_reagent('K3PO4', 'base')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
        pka = result.get('pka')
        if pka:
            print(f"     → pKa: {pka}")
    
    # Try a ligand
    result = find_reagent('BINAP', 'ligand')
    if result:
        print(f"\n  ✅ find_reagent('BINAP', 'ligand')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
        abbr = result.get('abbreviation', [])
        if abbr:
            print(f"     → Abbreviations: {', '.join(abbr)}")
    
    # Try a base
    result = find_reagent('Cesium carbonate', 'base')
    if result:
        print(f"\n  ✅ find_reagent('Cesium carbonate', 'base')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
    
    # Try a metal precursor
    result = find_reagent('Pd(OAc)2', 'metal_precursor')
    if result:
        print(f"\n  ✅ find_reagent('Pd(OAc)2', 'metal_precursor')")
        print(f"     → Name: {result.get('name', 'N/A')}")
        print(f"     → CAS: {result.get('cas', 'N/A')}")
    
    print(f"\n  💡 Available databases: ligand, base, solvent, metal_precursor, additive")


def test_6_drfp_similarity():
    """Demonstrate DRFP similarity calculations."""
    print("\n" + "="*70)
    print("  6. DRFP Similarity")
    print("="*70)
    
    if not drfp_available():
        print("⚠️  DRFP not installed")
        print("   Install: pip install drfp")
        return
    
    # Compare similar reactions
    rxn1 = 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'
    rxn2 = 'Clc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'
    
    try:
        similarity = drfp_tanimoto(rxn1, rxn2)
        print(f"✅ drfp_tanimoto(rxn1, rxn2)")
        print(f"   → Rxn1: {rxn1}")
        print(f"   → Rxn2: {rxn2}")
        print(f"   → Similarity: {similarity:.3f}")
        print(f"   💡 High similarity (only LG differs: Br vs Cl)")
    except Exception as e:
        print(f"❌ Error: {e}")
    
    # Compare different reactions
    rxn3 = 'Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1'  # Suzuki
    
    try:
        similarity = drfp_tanimoto(rxn1, rxn3)
        print(f"\n✅ drfp_tanimoto(rxn1, rxn3)")
        print(f"   → Rxn1: {rxn1} (C-N)")
        print(f"   → Rxn3: {rxn3} (Suzuki)")
        print(f"   → Similarity: {similarity:.3f}")
        print(f"   💡 Lower similarity (different reaction types)")
    except Exception as e:
        print(f"❌ Error: {e}")


def show_import_patterns():
    """Display import patterns."""
    print("\n" + "="*70)
    print("  Import Patterns")
    print("="*70)
    print("""
✅ SMILES Operations:
   from chemtools.smiles import normalize, normalize_reaction

✅ Family Detection:
   from chemtools import detect_reaction  # ⭐ Unified detection API

✅ Featurization:
   from chemtools.featurizers.molecular import featurize

✅ Reagent Database:
   from chemtools.reagent import find_reagent  # Full database (5000+)
   # Databases: ligand, base, solvent, metal_precursor, additive

✅ DRFP Similarity (optional):
   from chemtools.reaction_similarity import drfp_tanimoto
""")


def show_tips():
    """Display usage tips."""
    print("\n" + "="*70)
    print("  Tips & Best Practices")
    print("="*70)
    print("""
🎯 Family Detection:
   - Use detect_reaction(rxn, use_ml=True) for ML-enhanced detection
   - Use detect_reaction(rxn, use_ml=False) for rule-based only
   - Automatically detects catalysts (Pd/Cu) from agents
   - Pd catalyst → Buchwald_CN, Cu catalyst → Ullmann_CN
   - Optional ML via use_rxn_insight=True

⚡ Performance:
   - SMILES normalization is fast (~1ms per reaction)
   - Featurization is deterministic and cached
   - DRFP calculations are slower (~10-50ms)

📚 Next Steps:
   - See demo_recommendations.py for condition recommendations
   - See CHEMTOOLS_QUICKSTART.md for detailed guide
   - See CHEMTOOLS_CLASS_GUIDE.md for advanced usage
""")


def main():
    """Run all basic tool demonstrations."""
    import sys
    import io
    # Fix Windows console encoding for emoji
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    
    print("="*70)
    print("  ChemTools Basic Tools Demo")
    print("  October 2025 - Core Utilities")
    print("="*70)
    
    tests = [
        test_1_smiles_normalization,
        test_2_family_detection,
        test_3_molecular_featurization,
        test_4_property_lookup,
        test_5_database_analytics,
        test_6_drfp_similarity,
    ]
    
    for test in tests:
        try:
            test()
        except Exception as e:
            print(f"\n❌ Error in {test.__name__}: {e}")
            import traceback
            traceback.print_exc()
    
    show_import_patterns()
    show_tips()
    
    print("\n" + "="*70)
    print("  ✅ Basic Tools Demo Complete!")
    print("="*70)
    print("\n📖 Next: python tests/demo_recommendations.py")
    print("🎨 Or try: python app/ui_gradio.py\n")


if __name__ == "__main__":
    main()
