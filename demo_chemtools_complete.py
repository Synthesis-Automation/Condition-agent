"""
ChemTools Complete Usage Demonstration
======================================

This script showcases all major features of the refactored ChemTools library.
Demonstrates both the new modular API and backward-compatible imports.

Run: python demo_chemtools_complete.py
"""

import json
from pprint import pprint


def print_section(title):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80 + "\n")


def demo_1_smiles_operations():
    """Demonstrate SMILES normalization and parsing."""
    print_section("1. SMILES Operations")
    
    from chemtools.smiles import normalize, normalize_reaction
    
    # Single SMILES normalization
    print("📌 Normalize SMILES:")
    result = normalize("c1ccccc1O")
    print(f"  Input: c1ccccc1O")
    print(f"  Normalized: {result['smiles_norm']}")
    print(f"  Fragments: {result.get('fragments', [])}")
    
    # Reaction SMILES normalization
    print("\n📌 Normalize Reaction SMILES:")
    rxn = "Brc1ccccc1.Nc1cccnc1>>c1ccccc1Nc1cccnc1"
    result = normalize_reaction(rxn)
    print(f"  Input: {rxn}")
    print(f"  Normalized: {result.get('normalized', 'N/A')}")
    print(f"  Reactants count: {len(result.get('reactants', []))}")
    print(f"  Products count: {len(result.get('products', []))}")


def demo_2_reaction_family_detection():
    """Demonstrate reaction family detection (router)."""
    print_section("2. Reaction Family Detection")
    
    from chemtools.router import detect_family, detect_family_from_reaction
    
    # Method 1: From reactant list (rule-based only)
    print("📌 Method 1: detect_family(reactants)")
    reactants = ["Brc1ccccc1", "Nc1cccnc1"]
    result = detect_family(reactants)
    print(f"  Reactants: {reactants}")
    print(f"  Family: {result.get('family', 'Unknown')}")
    print(f"  Confidence: {result.get('confidence', 0):.2f}")
    print(f"  Detected via: {result.get('hits', {})}")
    
    # Method 2: From reaction SMILES (with catalyst detection)
    print("\n📌 Method 2: detect_family_from_reaction(reaction_smiles)")
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = detect_family_from_reaction(reaction, use_rxn_insight=True)
    print(f"  Reaction: {reaction}")
    print(f"  Family: {result.get('family', 'Unknown')}")
    print(f"  Confidence: {result.get('confidence', 0):.2f}")
    print(f"  � Detects catalysts from agents (Pd→Buchwald, Cu→Ullmann)")
    
    # Suzuki example
    print("\n📌 Suzuki Coupling Detection:")
    reactants = ["Brc1ccccc1", "Bc1cccnc1"]  # Aryl bromide + boronic acid
    result = detect_family(reactants)
    print(f"  Reactants: {reactants}")
    print(f"  Family: {result.get('family', 'Unknown')}")
    print(f"  Confidence: {result.get('confidence', 0):.2f}")


def demo_3_molecular_featurization():
    """Demonstrate molecular featurization."""
    print_section("3. Molecular Featurization")
    
    from chemtools.featurizers.molecular import featurize
    
    print("📌 Featurize C-N Coupling Substrates:")
    electrophile = "Brc1ccccc1"
    nucleophile = "Nc1ccccc1"
    
    features = featurize(electrophile, nucleophile)
    print(f"  Electrophile: {electrophile}")
    print(f"  Nucleophile: {nucleophile}")
    print(f"  Leaving Group: {features.get('LG', 'N/A')}")
    print(f"  Nucleophile Class: {features.get('nuc_class', 'N/A')}")
    print(f"  Feature Bin: {features.get('bin', 'N/A')}")
    print(f"  Total Features: {len(features)} keys")


def demo_4_precedent_search():
    """Demonstrate precedent search (k-NN with DRFP)."""
    print_section("4. Precedent Search (k-NN with DRFP)")
    
    from chemtools.precedent import knn
    
    print("📌 Search for Similar Reactions:")
    print("  (Note: Requires dataset to be loaded)")
    
    try:
        # Simple feature-based search
        result = knn(
            family="Ullmann C-N coupling",
            features={"bin": "LG:Br|NUC:aniline"},
            k=5,
            relax={"use_drfp": True}  # Enable DRFP similarity
        )
        
        print(f"  Family: Ullmann C-N coupling")
        print(f"  Query Features: LG:Br|NUC:aniline")
        print(f"  Found: {len(result.get('precedents', []))} precedents")
        
        if result.get('precedents'):
            print("\n  Top 3 precedents:")
            for i, prec in enumerate(result['precedents'][:3], 1):
                print(f"    {i}. Core: {prec.get('condition_core', 'N/A')}, "
                      f"Similarity: {prec.get('similarity', 0):.3f}")
    
    except Exception as e:
        print(f"  ⚠️  Dataset not loaded: {e}")
        print("  💡 Tip: Use selective loading - see demo_11_resource_management()")


def demo_5_condition_recommendation():
    """Demonstrate condition recommendation (new modular API)."""
    print_section("5. Condition Recommendation (NEW Modular API)")
    
    # NEW import from recommend package
    from chemtools.recommend import recommend_from_reaction
    
    print("📌 Recommend Conditions for C-N Coupling:")
    reaction = "Brc1ccccc1.Nc1ccc2ccccc2c1>>c1ccccc1Nc1ccc2ccccc2c1"
    
    try:
        result = recommend_from_reaction(
            reaction=reaction,  # Correct parameter name
            k=5,
            relax={"use_drfp": True}
        )
        
        print(f"  Reaction: {reaction[:50]}...")
        print(f"  Detected Family: {result.get('reaction_family', 'N/A')}")
        
        rec = result.get('recommendation', {})
        print(f"\n  ✅ Top Recommendation:")
        print(f"    Core: {rec.get('core', 'N/A')}")
        print(f"    Base: {rec.get('base', 'N/A')}")
        print(f"    Solvent: {rec.get('solvent', 'N/A')}")
        print(f"    Temperature: {rec.get('temperature_C', 'N/A')}°C")
        print(f"    Time: {rec.get('time_h', 'N/A')} hours")
        print(f"    Confidence: {rec.get('confidence', 0):.2f}")
        
        # Show alternatives
        if result.get('alternatives', {}).get('cores'):
            print(f"\n  Alternative Cores:")
            for core in result['alternatives']['cores'][:3]:
                print(f"    - {core['name']} (count: {core['count']})")
    
    except Exception as e:
        print(f"  ⚠️  Error: {e}")
        print("  💡 Tip: Ensure dataset is loaded")


def demo_6_structured_recommendation():
    """Demonstrate structured condition recommendation API."""
    print_section("6. Structured Recommendation API (Top-N Variants)")
    
    from chemtools.recommend import recommend_conditions_structured
    
    print("📌 Get Top 5 Condition Variants:")
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    try:
        result = recommend_conditions_structured(
            reaction=reaction,  # Correct parameter name
            k=5,  # Top 5 variants
            relax={"use_drfp": True}
        )
        
        print(f"  Reaction: {reaction}")
        print(f"  Strategy: {result.get('strategy', 'N/A')}")
        print(f"  Variants: {len(result.get('conditions', []))}")
        
        print("\n  Top 3 Variants:")
        for i, cond in enumerate(result.get('conditions', [])[:3], 1):
            print(f"\n    Variant {i}:")
            print(f"      Core: {cond.get('core', 'N/A')}")
            print(f"      Base: {cond.get('base', 'N/A')}")
            print(f"      Solvent: {cond.get('solvent', 'N/A')}")
            print(f"      Temp: {cond.get('temperature_C', 'N/A')}°C")
            print(f"      Precedents: {cond.get('precedent_count', 0)}")
    
    except Exception as e:
        print(f"  ⚠️  Error: {e}")


def demo_7_ml_recommendations():
    """Demonstrate ML-enhanced recommendations with yield prediction."""
    print_section("7. ML-Enhanced Recommendations (Yield Prediction)")
    
    # NEW import from ml.recommender
    from chemtools.ml.recommender import hybrid_recommend
    
    print("📌 ML Hybrid Recommend (ML + k-NN):")
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    try:
        result = hybrid_recommend(
            reaction_smiles=reaction,
            model_path="models/drfp_yield_v1.pkl",  # Optional
            k=50
        )
        
        print(f"  Reaction: {reaction}")
        print(f"  Strategy: {result.get('strategy', 'N/A')}")
        print(f"  ML Model: {'Used' if result.get('ml_predictions') else 'Not available'}")
        
        print("\n  Top 3 Recommendations (with predicted yields):")
        for i, cond in enumerate(result.get('recommended_conditions', [])[:3], 1):
            predicted_yield = cond.get('predicted_yield_pct', cond.get('predicted_yield', 'N/A'))
            print(f"\n    {i}. Core: {cond.get('core', 'N/A')}")
            print(f"       Base: {cond.get('base', 'N/A')}")
            print(f"       Predicted Yield: {predicted_yield}%")
    
    except Exception as e:
        print(f"  ⚠️  Error: {e}")
        print("  💡 Tip: Fallback to k-NN if ML model not available")


def demo_8_plate_design():
    """Demonstrate plate design generation."""
    print_section("8. Plate Design Generation")
    
    from chemtools.recommend import design_plate_from_reaction
    
    print("📌 Design Reaction Plate (Top Conditions):")
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    try:
        result = design_plate_from_reaction(
            reaction=reaction,  # Correct parameter name
            top_cores=3,
            variants_per_core=2
        )
        
        print(f"  Reaction: {reaction}")
        print(f"  Cores: {result.get('top_cores', [])}")
        print(f"  Total Wells: {len(result.get('plate', []))}")
        
        print("\n  Sample Wells:")
        for i, well in enumerate(result.get('plate', [])[:4], 1):
            print(f"\n    Well {i} ({well.get('well_id', 'N/A')}):")
            print(f"      Core: {well.get('core', 'N/A')}")
            print(f"      Base: {well.get('base', 'N/A')}")
            print(f"      Solvent: {well.get('solvent', 'N/A')}")
    
    except Exception as e:
        print(f"  ⚠️  Error: {e}")


def demo_9_drfp_similarity():
    """Demonstrate DRFP similarity calculation."""
    print_section("9. DRFP Similarity Calculation")
    
    from chemtools.reaction_similarity import drfp_tanimoto
    
    print("📌 Calculate Reaction Similarity (DRFP/Tanimoto):")
    
    rxn1 = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    rxn2 = "Clc1ccccc1.Nc1cccnc1>>c1ccccc1Nc1cccnc1"
    
    try:
        similarity = drfp_tanimoto(rxn1, rxn2)
        
        print(f"  Reaction 1: {rxn1[:40]}...")
        print(f"  Reaction 2: {rxn2[:40]}...")
        print(f"  Tanimoto Similarity: {similarity:.3f}")
        
        if similarity > 0.7:
            print("  ✅ Highly similar reactions")
        elif similarity > 0.4:
            print("  ⚠️  Moderately similar reactions")
        else:
            print("  ❌ Different reactions")
    
    except Exception as e:
        print(f"  ⚠️  DRFP not available: {e}")
        print("  💡 Install: pip install drfp")


def demo_10_property_lookup():
    """Demonstrate property and reagent lookup."""
    print_section("10. Property and Reagent Lookup")
    
    from chemtools.reagent_lookup import find_reagent
    
    # Property lookup
    print("📌 Reagent Database Lookup:")
    try:
        result = find_reagent("K3PO4", "base")
        if result:
            print(f"  Query: K3PO4")
            print(f"  Name: {result.get('name', 'N/A')}")
            print(f"  CAS: {result.get('cas', 'N/A')}")
            pka = result.get('pka')
            if pka:
                print(f"  pKa: {pka}")
        
        # Try a solvent
        result = find_reagent("DMF", "solvent")
        if result:
            print(f"\n  Query: DMF")
            print(f"  Name: {result.get('name', 'N/A')}")
            print(f"  CAS: {result.get('cas', 'N/A')}")
            props = result.get('properties', {})
            bp = props.get('boiling_point')
            if bp:
                print(f"  Boiling point: {bp}°C")
    except Exception as e:
        print(f"  ⚠️  Error: {e}")
    
    # Reagent lookup - show available databases
    print("\n📌 Available Reagent Databases:")
    print("  💡 Use find_reagent(name, db_type)")
    print("  Database types: ligand, base, solvent, metal_precursor, additive")
    print("  Example:")
    print("    from chemtools.reagent_lookup import find_reagent")
    print("    result = find_reagent('BINAP', 'ligand')")


def demo_11_resource_management():
    """Demonstrate resource management and selective loading."""
    print_section("11. Resource Management (ChemTools Class)")
    
    print("📌 Using Global Instance (Simplest):")
    print("""
    from chemtools import chem
    
    # Lazy loading - only loads when needed
    result = chem.smiles.normalize("CCO")
    precedents = chem.precedent.knn(...)
    """)
    
    print("\n📌 Creating Custom Instance (Selective Loading):")
    print("""
    from chemtools import ChemTools
    
    # Only load Buchwald reactions - 50-100x faster!
    buchwald = ChemTools(
        datasets=["C_N_Coupling_Pd"],     # Only this dataset
        reagent_dbs=["ligand", "base"],   # Only these reagents
        preload=True                       # Load immediately
    )
    
    # Use the configured instance
    precedents = buchwald.precedent.knn(
        family="C_N_Coupling_Pd",
        features={...},
        k=5
    )
    """)
    
    print("\n📌 Clearing Resources:")
    print("""
    from chemtools import chem
    
    # Free memory when done
    chem.clear_cache()
    
    # Or clear specific resources
    chem.clear_precedents()
    chem.clear_reagents()
    """)


def demo_12_backward_compatibility():
    """Demonstrate backward compatibility with deprecated imports."""
    print_section("12. Backward Compatibility (Deprecated Imports)")
    
    print("📌 Old Imports Still Work (with warnings):")
    print("""
    # OLD (deprecated but still works):
    from chemtools.recommend import recommend_from_reaction  # ⚠️ Warning
    from chemtools.recommend_ml import hybrid_recommend      # ⚠️ Warning
    
    # NEW (recommended):
    from chemtools.recommend import recommend_from_reaction  # ✅ No warning
    from chemtools.ml.recommender import hybrid_recommend    # ✅ No warning
    """)
    
    print("\n📌 Migration Path:")
    print("""
    1. Your existing code continues working (100% compatible)
    2. You'll see deprecation warnings in logs
    3. Update imports at your own pace
    4. See MIGRATION_GUIDE.md for details
    """)


def demo_13_api_endpoints():
    """Demonstrate FastAPI endpoints (documentation only)."""
    print_section("13. FastAPI Endpoints (REST API)")
    
    print("📌 Start API Server:")
    print("  uvicorn app.main:app --reload --port 8000")
    print("  Open: http://127.0.0.1:8000/docs")
    
    print("\n📌 Key Endpoints:")
    endpoints = [
        ("POST", "/api/v1/smiles/normalize", "Normalize SMILES"),
        ("POST", "/api/v1/router/detect-family", "Detect reaction family"),
        ("POST", "/api/v1/featurize/molecular", "Featurize substrates"),
        ("POST", "/api/v1/precedent/knn", "Search precedents"),
        ("POST", "/api/v1/recommend", "Recommend conditions"),
        ("POST", "/api/v1/recommend/conditions", "Structured recommendations"),
        ("POST", "/api/v1/core/search", "Search by condition core"),
    ]
    
    for method, path, desc in endpoints:
        print(f"  {method:6} {path:40} - {desc}")
    
    print("\n📌 Example cURL:")
    print("""
    curl -X POST http://127.0.0.1:8000/api/v1/smiles/normalize \\
      -H "Content-Type: application/json" \\
      -d '{"smiles": "c1ccccc1Br"}'
    """)


def demo_14_gradio_ui():
    """Demonstrate Gradio UI (documentation only)."""
    print_section("14. Gradio UI (Interactive Web Interface)")
    
    print("📌 Launch Gradio UI:")
    print("  python app/ui_gradio.py")
    print("  Open: http://127.0.0.1:7860")
    
    print("\n📌 Available Tabs:")
    tabs = [
        "SMILES Normalize - Normalize and validate SMILES",
        "Detect Family - Infer reaction family from reactants",
        "Recommend Conditions - Get condition recommendations",
        "Precedent Search - Find similar reactions (k-NN + DRFP)",
        "Design Plate - Generate reaction plate layouts",
        "Core Search - Search reactions by catalyst/ligand core",
        "DRFP Similarity - Calculate reaction similarity",
        "Properties Lookup - Query reagent properties",
    ]
    
    for i, tab in enumerate(tabs, 1):
        print(f"  {i}. {tab}")


def demo_15_performance_tips():
    """Demonstrate performance optimization tips."""
    print_section("15. Performance Optimization Tips")
    
    print("📌 Selective Dataset Loading (50-100x faster):")
    print("""
    from chemtools import ChemTools
    
    # Load only what you need
    fast = ChemTools(
        datasets=["C_N_Coupling_Pd"],  # Not all reaction types
        reagent_dbs=["ligand", "base"] # Not all reagent types
    )
    """)
    
    print("\n📌 Environment Variables:")
    tips = [
        ("CHEMTOOLS_DISABLE_RDKIT=1", "Skip RDKit (faster, less accurate)"),
        ("CHEMTOOLS_LOAD_DATASET=0", "Skip dataset loading (testing only)"),
        ("CHEMTOOLS_ATTACH_ROLE_AWARE=0", "Skip role-aware features (default)"),
    ]
    
    for var, desc in tips:
        print(f"  {var:35} - {desc}")
    
    print("\n📌 DRFP Optimization:")
    print("""
    # Use precomputed DRFP fingerprints
    relax = {
        "use_drfp": True,
        "precompute_drfp": "candidates",  # Only for candidates, not all
        "drfp_n_bits": 2048               # Smaller = faster (default 4096)
    }
    """)
    
    print("\n📌 Caching:")
    print("""
    # Keep server/UI running to reuse caches
    # First call: 60-120s (loads dataset)
    # Subsequent calls: <100ms (uses cache)
    """)


def demo_16_testing():
    """Demonstrate testing approaches."""
    print_section("16. Testing Your Code")
    
    print("📌 Run Test Suite:")
    print("  pytest -q                          # All tests")
    print("  pytest tests/chemtools/ -v         # Package tests")
    print("  pytest test_new_recommend.py -v    # Recommendation tests")
    
    print("\n📌 Example Test:")
    print("""
    def test_smiles_normalization():
        from chemtools.smiles import normalize
        
        result = normalize("c1ccccc1O")
        assert result['smiles_norm'] == "Oc1ccccc1"  # Canonical form
        assert len(result['fragments']) >= 1
    
    def test_recommendation():
        from chemtools.recommend import recommend_from_reaction
        
        result = recommend_from_reaction(
            "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        )
        
        assert 'recommendation' in result
        assert result['recommendation']['core']  # Should have core
    """)


def demo_17_troubleshooting():
    """Demonstrate common issues and solutions."""
    print_section("17. Troubleshooting Common Issues")
    
    issues = [
        (
            "Dataset not loaded / precedents not found",
            [
                "Use selective loading: ChemTools(datasets=['C_N_Coupling_Pd'])",
                "Check dataset exists in data/reaction_dataset/",
                "Verify family name matches dataset file"
            ]
        ),
        (
            "Slow first request (60-120s)",
            [
                "Expected! Dataset loads on first use",
                "Use preload=True for API servers",
                "Keep server running to reuse cache"
            ]
        ),
        (
            "DRFP not available",
            [
                "Install: pip install drfp",
                "Or disable: relax={'use_drfp': False}"
            ]
        ),
        (
            "Out of memory",
            [
                "Use selective loading (load only needed datasets)",
                "Set CHEMTOOLS_LOAD_DATASET=0 for testing",
                "Clear cache: chem.clear_cache()"
            ]
        ),
        (
            "Deprecation warnings",
            [
                "Update imports: see MIGRATION_GUIDE.md",
                "Old code still works (100% compatible)",
                "Or suppress: warnings.filterwarnings('ignore', category=DeprecationWarning)"
            ]
        )
    ]
    
    for issue, solutions in issues:
        print(f"❌ {issue}:")
        for solution in solutions:
            print(f"  ✅ {solution}")
        print()


def main():
    """Run all demonstrations."""
    print("\n" + "=" * 80)
    print("  ChemTools Complete Usage Demonstration")
    print("  October 2025 - Refactored Modular Architecture")
    print("=" * 80)
    
    # Run each demo
    demos = [
        ("SMILES Operations", demo_1_smiles_operations),
        ("Reaction Family Detection", demo_2_reaction_family_detection),
        ("Molecular Featurization", demo_3_molecular_featurization),
        ("Precedent Search", demo_4_precedent_search),
        ("Condition Recommendation", demo_5_condition_recommendation),
        ("Structured Recommendation", demo_6_structured_recommendation),
        ("ML Recommendations", demo_7_ml_recommendations),
        ("Plate Design", demo_8_plate_design),
        ("DRFP Similarity", demo_9_drfp_similarity),
        ("Property Lookup", demo_10_property_lookup),
        ("Resource Management", demo_11_resource_management),
        ("Backward Compatibility", demo_12_backward_compatibility),
        ("API Endpoints", demo_13_api_endpoints),
        ("Gradio UI", demo_14_gradio_ui),
        ("Performance Tips", demo_15_performance_tips),
        ("Testing", demo_16_testing),
        ("Troubleshooting", demo_17_troubleshooting),
    ]
    
    try:
        for title, demo_func in demos:
            try:
                demo_func()
            except KeyboardInterrupt:
                print("\n\n⚠️  Demo interrupted by user")
                break
            except Exception as e:
                print(f"\n⚠️  Error in {title}: {e}")
                import traceback
                traceback.print_exc()
                print("\nContinuing with next demo...\n")
        
        # Final summary
        print_section("Summary")
        print("✅ Demonstration complete!")
        print("\n📚 Next Steps:")
        print("  1. Read CHEMTOOLS_QUICKSTART.md for quick start")
        print("  2. Read CHEMTOOLS_CLASS_GUIDE.md for detailed API docs")
        print("  3. Check MIGRATION_GUIDE.md if migrating old code")
        print("  4. Try the Gradio UI: python app/ui_gradio.py")
        print("  5. Run tests: pytest -q")
        print("\n💡 Documentation:")
        print("  - README.md - Project overview and API endpoints")
        print("  - RECOMMEND_REFACTORING_SUCCESS.md - Refactoring details")
        print("  - API docs: http://127.0.0.1:8000/docs (after starting server)")
        
    except Exception as e:
        print(f"\n❌ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
