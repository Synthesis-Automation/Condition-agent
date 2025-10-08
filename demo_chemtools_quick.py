"""
ChemTools Quick Demo - Working Examples
========================================

Simple, working demonstrations of key ChemTools features.
Tests actual API calls with the refactored modules.

Run: python demo_chemtools_quick.py
"""


def section(title):
    print("\n" + "="*70)
    print(f"  {title}")
    print("="*70)


def test_1_smiles():
    section("1. SMILES Normalization")
    from chemtools.smiles import normalize, normalize_reaction
    
    # Single molecule
    result = normalize("c1ccccc1O")
    print(f"✅ normalize('c1ccccc1O')")
    print(f"   → {result['smiles_norm']}")
    
    # Reaction
    rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = normalize_reaction(rxn)
    print(f"\n✅ normalize_reaction('{rxn}')")
    print(f"   → Normalized: {result.get('normalized', 'N/A')}")
    print(f"   → Reactants: {len(result.get('reactants', []))}")
    print(f"   → Products: {len(result.get('products', []))}")


def test_2_router():
    section("2. Reaction Family Detection")
    from chemtools.router import detect_family, detect_family_from_reaction
    
    # Method 1: From reactants
    reactants = ["Brc1ccccc1", "Nc1ccccc1"]
    result = detect_family(reactants)
    print(f"✅ detect_family({reactants})")
    print(f"   → Family: {result.get('family', 'Unknown')}")
    print(f"   → Confidence: {result.get('confidence', 0):.2f}")
    
    # Method 2: From reaction SMILES (better - detects catalysts)
    rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = detect_family_from_reaction(rxn)
    print(f"\n✅ detect_family_from_reaction('{rxn}')")
    print(f"   → Family: {result.get('family', 'Unknown')}")
    print(f"   → Confidence: {result.get('confidence', 0):.2f}")
    print(f"   💡 Also detects catalysts (Pd/Cu) from agents")


def test_3_featurize():
    section("3. Molecular Featurization")
    from chemtools.featurizers.molecular import featurize
    
    result = featurize("Brc1ccccc1", "Nc1ccccc1")
    print(f"✅ featurize('Brc1ccccc1', 'Nc1ccccc1')")
    print(f"   → LG: {result.get('LG', 'N/A')}")
    print(f"   → Nuc class: {result.get('nuc_class', 'N/A')}")
    print(f"   → Bin: {result.get('bin', 'N/A')}")


def test_4_precedent():
    section("4. Precedent Search (k-NN)")
    from chemtools.precedent import knn
    
    print("✅ knn(family='Ullmann C-N coupling', features={...}, k=5)")
    try:
        result = knn(
            family="Ullmann C-N coupling",
            features={"bin": "LG:Br|NUC:aniline"},
            k=5,
            relax={"use_drfp": False}  # Faster for demo
        )
        print(f"   → Found: {len(result.get('precedents', []))} precedents")
        if result.get('precedents'):
            p = result['precedents'][0]
            print(f"   → Top: {p.get('condition_core', 'N/A')}")
    except Exception as e:
        print(f"   ⚠️  {str(e)[:60]}...")
        print("   💡 Dataset not loaded - use selective loading")


def test_5_recommend():
    section("5. Condition Recommendation (NEW API)")
    from chemtools.recommend import recommend_from_reaction
    
    rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    print(f"✅ recommend_from_reaction('{rxn}')")
    try:
        result = recommend_from_reaction(
            reaction=rxn,
            k=5,
            relax={"use_drfp": False}  # Faster
        )
        rec = result.get('recommendation', {})
        print(f"   → Family: {result.get('family', 'N/A')}")
        print(f"   → Core: {rec.get('core', 'N/A')}")
        print(f"   → Base: {rec.get('base', 'N/A')}")
        print(f"   → Confidence: {rec.get('confidence', 0):.2f}")
    except Exception as e:
        print(f"   ⚠️  {str(e)[:60]}...")


def test_6_ml():
    section("6. ML-Enhanced Recommendations")
    from chemtools.ml.recommender import hybrid_recommend
    
    rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    print(f"✅ hybrid_recommend('{rxn}')")
    try:
        result = hybrid_recommend(
            reaction_smiles=rxn,
            k=20
        )
        print(f"   → Strategy: {result.get('strategy', 'N/A')}")
        print(f"   → Conditions: {len(result.get('recommended_conditions', []))}")
        if result.get('recommended_conditions'):
            c = result['recommended_conditions'][0]
            print(f"   → Top: {c.get('core', 'N/A')}")
    except Exception as e:
        print(f"   ⚠️  {str(e)[:60]}...")


def test_7_drfp():
    section("7. DRFP Similarity")
    
    try:
        from chemtools.reaction_similarity import drfp_tanimoto
        
        rxn1 = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        rxn2 = "Clc1ccccc1.Nc1cccnc1>>c1ccccc1Nc1cccnc1"
        
        sim = drfp_tanimoto(rxn1, rxn2)
        print(f"✅ drfp_tanimoto(rxn1, rxn2)")
        print(f"   → Similarity: {sim:.3f}")
    except ImportError:
        print("⚠️  DRFP not installed")
        print("   Install: pip install drfp")
    except Exception as e:
        print(f"   ⚠️  {str(e)[:60]}...")


def test_8_properties():
    section("8. Reagent Lookup")
    from chemtools.reagent_lookup import find_reagent
    
    try:
        result = find_reagent("K3PO4", "base")
        if result:
            print(f"✅ find_reagent('K3PO4', 'base')")
            print(f"   → Name: {result.get('name', 'N/A')}")
            print(f"   → CAS: {result.get('cas', 'N/A')}")
        
        result = find_reagent("DMF", "solvent")
        if result:
            print(f"✅ find_reagent('DMF', 'solvent')")
            print(f"   → Name: {result.get('name', 'N/A')}")
            print(f"   → CAS: {result.get('cas', 'N/A')}")
    except Exception as e:
        print(f"   ⚠️  {str(e)[:60]}...")


def show_imports():
    section("NEW Import Patterns")
    
    print("""
✅ SMILES Operations:
   from chemtools.smiles import normalize, normalize_reaction

✅ Family Detection:
   from chemtools.router import detect_family
   from chemtools.router import detect_family_from_reaction  # ⭐ Accepts reaction SMILES

✅ Featurization:
   from chemtools.featurizers.molecular import featurize

✅ Precedent Search:
   from chemtools.precedent import knn

✅ Recommendations (NEW modular package):
   from chemtools.recommend import recommend_from_reaction
   from chemtools.recommend import recommend_conditions_structured
   from chemtools.recommend import design_plate_from_reaction

✅ ML Recommendations (NEW location):
   from chemtools.ml.recommender import hybrid_recommend

✅ DRFP Similarity:
   from chemtools.reaction_similarity import drfp_tanimoto

✅ Reagent Database:
   from chemtools.reagent_lookup import find_reagent
   # Databases: ligand, base, solvent, metal_precursor, additive
""")


def show_api_signature():
    section("Key Function Signatures")
    
    print("""
recommend_from_reaction(
    reaction: str,           # NOT reaction_smiles!
    k: int = 25,
    relax: Dict | None = None,
    constraint_rules: Dict | None = None,
    family_override: str | None = None,
    max_variants: int = 3
) -> Dict

hybrid_recommend(
    reaction_smiles: str,    # Note: different param name
    model_path: str | None = None,
    k: int = 50,
    relax: Dict | None = None
) -> Dict

design_plate_from_reaction(
    reaction: str,           # NOT reaction_smiles!
    top_cores: int = 5,
    variants_per_core: int = 3,
    relax: Dict | None = None
) -> Dict
""")


def show_tips():
    section("Performance Tips")
    
    print("""
🚀 Selective Dataset Loading (50-100x faster):
   from chemtools import ChemTools
   
   fast = ChemTools(
       datasets=["C_N_Coupling_Pd"],  # Only this reaction type
       reagent_dbs=["ligand", "base"] # Only these reagents
   )
   
   precedents = fast.precedent.knn(...)

🔧 Environment Variables:
   $env:CHEMTOOLS_DISABLE_RDKIT='1'      # Faster (Windows)
   $env:CHEMTOOLS_LOAD_DATASET='0'       # Testing only
   
   export CHEMTOOLS_DISABLE_RDKIT=1      # Faster (Linux/Mac)
   export CHEMTOOLS_LOAD_DATASET=0       # Testing only

⚡ DRFP Optimization:
   relax = {
       "use_drfp": True,
       "precompute_drfp": "candidates",  # Not "all"
       "drfp_n_bits": 2048               # Smaller = faster
   }

📚 Documentation:
   - CHEMTOOLS_QUICKSTART.md
   - CHEMTOOLS_CLASS_GUIDE.md  
   - MIGRATION_GUIDE.md
   - RECOMMEND_REFACTORING_SUCCESS.md
""")


def main():
    import sys
    import io
    # Fix Windows console encoding for emoji
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    
    print("\n" + "="*70)
    print("  ChemTools Quick Demo - Working Examples")
    print("  October 2025 - Refactored Modular Architecture")
    print("="*70)
    
    # Run tests
    tests = [
        test_1_smiles,
        test_2_router,
        test_3_featurize,
        test_4_precedent,
        test_5_recommend,
        test_6_ml,
        test_7_drfp,
        test_8_properties,
    ]
    
    for test in tests:
        try:
            test()
        except Exception as e:
            print(f"\n❌ Error: {e}")
            import traceback
            traceback.print_exc()
    
    # Show reference info
    show_imports()
    show_api_signature()
    show_tips()
    
    print("\n" + "="*70)
    print("  ✅ Demo Complete!")
    print("="*70)
    print("\n📖 Next: Try the Gradio UI → python app/ui_gradio.py")
    print("🧪 Or run tests → pytest -q\n")


if __name__ == "__main__":
    main()
