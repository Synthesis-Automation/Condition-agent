#!/usr/bin/env python3
"""
Edge case and stress tests for refactored core.py modules.

Tests unusual inputs, boundary conditions, and potential failure modes.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))


class EdgeCaseTests:
    """Test edge cases and boundary conditions"""
    
    def test_empty_reaction_string(self):
        """Test with empty reaction string"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        print("Testing empty reaction string...")
        try:
            result = recommend_from_reaction("", k=5, rerank_strategy='none')
            print(f"  ⚠️  Empty string accepted: {type(result)}")
        except Exception as e:
            print(f"  ✅ Empty string rejected: {type(e).__name__}")
    
    def test_very_large_k_value(self):
        """Test with very large k parameter"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        print("\nTesting very large k value (k=10000)...")
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        try:
            result = recommend_from_reaction(reaction, k=10000, rerank_strategy='none')
            if isinstance(result, dict):
                prec_count = len(result.get("precedents", result.get("precedent_pack", {}).get("precedents", [])))
                print(f"  ✅ Large k handled, returned {prec_count} precedents")
        except Exception as e:
            print(f"  ⚠️  Large k failed: {type(e).__name__}: {e}")
    
    def test_zero_k_value(self):
        """Test with k=0"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        print("\nTesting k=0...")
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        try:
            result = recommend_from_reaction(reaction, k=0, rerank_strategy='none')
            print(f"  ✅ k=0 handled: {type(result)}")
        except Exception as e:
            print(f"  ⚠️  k=0 failed: {type(e).__name__}")
    
    def test_very_long_reaction_smiles(self):
        """Test with very long SMILES string"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        print("\nTesting very long SMILES...")
        # Repeat a polymer-like structure
        long_smiles = "C" * 1000 + ">>C" + "C" * 1000
        
        try:
            result = recommend_from_reaction(long_smiles, k=5, rerank_strategy='none')
            print(f"  ✅ Long SMILES handled: {type(result)}")
        except Exception as e:
            print(f"  ⚠️  Long SMILES failed: {type(e).__name__}")
    
    def test_malformed_reaction_formats(self):
        """Test various malformed reaction formats"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        print("\nTesting malformed reaction formats...")
        
        malformed_reactions = [
            "only_product",  # No >> separator
            ">>c1ccccc1",  # No reactants
            "c1ccccc1>>",  # No products
            "c1ccccc1>c1ccccc1",  # Single > instead of >>
            "c1ccccc1>>>c1ccccc1",  # Triple >>>
        ]
        
        for i, reaction in enumerate(malformed_reactions, 1):
            try:
                result = recommend_from_reaction(reaction, k=5, rerank_strategy='none')
                print(f"  ⚠️  Format {i} accepted: {reaction[:30]}...")
            except Exception as e:
                print(f"  ✅ Format {i} rejected: {type(e).__name__}")
    
    def test_special_characters_in_reaction(self):
        """Test reaction with special characters"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        print("\nTesting special characters...")
        
        special_reactions = [
            "[Na+].[Cl-].c1ccccc1>>c1ccccc1",  # Ions
            "C[C@H](O)c1ccccc1>>c1ccccc1",  # Stereochemistry
            "c1ccc([N+](=O)[O-])cc1>>c1ccccc1N",  # Charged groups
        ]
        
        for i, reaction in enumerate(special_reactions, 1):
            try:
                result = recommend_from_reaction(reaction, k=5, rerank_strategy='none')
                print(f"  ✅ Special chars {i} handled")
            except Exception as e:
                print(f"  ⚠️  Special chars {i} failed: {type(e).__name__}")
    
    def test_structured_with_zero_limit(self):
        """Test structured output with limit=0"""
        from chemtools.recommend.modules.structured import recommend_conditions_structured
        
        print("\nTesting structured with limit=0...")
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        try:
            result = recommend_conditions_structured(
                reaction=reaction,
                k=5,
                limit=0,
                rerank_strategy='none'
            )
            rec_count = len(result.get("recommendations", []))
            print(f"  ✅ limit=0 handled, returned {rec_count} recommendations")
        except Exception as e:
            print(f"  ⚠️  limit=0 failed: {type(e).__name__}")
    
    def test_structured_with_negative_limit(self):
        """Test structured output with negative limit"""
        from chemtools.recommend.modules.structured import recommend_conditions_structured
        
        print("\nTesting structured with limit=-1...")
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        try:
            result = recommend_conditions_structured(
                reaction=reaction,
                k=5,
                limit=-1,
                rerank_strategy='none'
            )
            rec_count = len(result.get("recommendations", []))
            print(f"  ✅ limit=-1 handled, returned {rec_count} recommendations")
        except Exception as e:
            print(f"  ⚠️  limit=-1 failed: {type(e).__name__}")
    
    def test_concurrent_imports(self):
        """Test importing modules concurrently"""
        import threading
        
        print("\nTesting concurrent imports...")
        errors = []
        
        def import_module(module_name):
            try:
                if module_name == "recommender":
                    from chemtools.recommend.modules.recommender import recommend_from_reaction
                elif module_name == "structured":
                    from chemtools.recommend.modules.structured import recommend_conditions_structured
                elif module_name == "output_builder":
                    from chemtools.recommend.modules.output_builder import build_formatted_output
            except Exception as e:
                errors.append((module_name, e))
        
        threads = [
            threading.Thread(target=import_module, args=("recommender",)),
            threading.Thread(target=import_module, args=("structured",)),
            threading.Thread(target=import_module, args=("output_builder",)),
        ]
        
        for t in threads:
            t.start()
        
        for t in threads:
            t.join()
        
        if errors:
            print(f"  ⚠️  Concurrent import errors: {len(errors)}")
            for mod, err in errors:
                print(f"     {mod}: {type(err).__name__}")
        else:
            print(f"  ✅ Concurrent imports successful")
    
    def test_multiple_rapid_calls(self):
        """Test making multiple rapid calls to test for state issues"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        print("\nTesting multiple rapid calls...")
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        try:
            results = []
            for i in range(10):
                result = recommend_from_reaction(reaction, k=3, rerank_strategy='none')
                results.append(result)
            
            # All results should be dict
            all_dicts = all(isinstance(r, dict) for r in results)
            print(f"  ✅ {len(results)} rapid calls successful, all dicts: {all_dicts}")
            
        except Exception as e:
            print(f"  ⚠️  Rapid calls failed: {type(e).__name__}: {e}")


class MemoryTests:
    """Test memory-related aspects"""
    
    def test_module_not_reloaded_on_import(self):
        """Test that modules are not reloaded on subsequent imports"""
        print("\nTesting module caching...")
        
        # Import once
        from chemtools.recommend.modules import recommender as mod1
        id1 = id(mod1)
        
        # Import again
        from chemtools.recommend.modules import recommender as mod2
        id2 = id(mod2)
        
        if id1 == id2:
            print(f"  ✅ Module cached properly (same id: {id1})")
        else:
            print(f"  ⚠️  Module reloaded (different ids: {id1} vs {id2})")
    
    def test_function_objects_shared(self):
        """Test that function objects are shared, not duplicated"""
        print("\nTesting function object sharing...")
        
        from chemtools.recommend.core import recommend_from_reaction as f1
        from chemtools.recommend.modules import recommend_from_reaction as f2
        from chemtools.recommend.modules.recommender import recommend_from_reaction as f3
        
        id1, id2, id3 = id(f1), id(f2), id(f3)
        
        if id1 == id2 == id3:
            print(f"  ✅ Function object shared (same id: {id1})")
        else:
            print(f"  ⚠️  Function objects not shared:")
            print(f"     core: {id1}")
            print(f"     modules: {id2}")
            print(f"     direct: {id3}")


class IntegrationTests:
    """Test integration between modules"""
    
    def test_recommender_to_structured_flow(self):
        """Test that recommender output can feed into structured"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        from chemtools.recommend.modules.structured import recommend_conditions_structured
        
        print("\nTesting recommender → structured flow...")
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        try:
            # Get basic recommendation
            basic_result = recommend_from_reaction(reaction, k=5, rerank_strategy='none')
            
            # Get structured recommendation
            struct_result = recommend_conditions_structured(
                reaction=reaction, k=5, limit=2, rerank_strategy='none'
            )
            
            # Both should succeed
            assert isinstance(basic_result, dict)
            assert isinstance(struct_result, dict)
            
            print(f"  ✅ Both recommender and structured work")
            
        except Exception as e:
            print(f"  ⚠️  Integration failed: {type(e).__name__}: {e}")
    
    def test_precedent_builder_integration(self):
        """Test precedent builder integration with recommender"""
        print("\nTesting precedent builder integration...")
        
        try:
            from chemtools.recommend.modules.recommender import recommend_from_reaction
            from chemtools.recommend.modules.precedent_builder import build_precedent_details
            
            reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
            
            # Get recommendation
            result = recommend_from_reaction(reaction, k=10, rerank_strategy='none')
            
            # Extract precedents
            precedents = result.get("precedents", result.get("precedent_pack", {}).get("precedents", []))
            
            if precedents:
                # Build details with correct signature
                details = build_precedent_details(
                    precs=precedents,
                    chosen_core=None,
                    group=precedents[:5]
                )
                assert isinstance(details, dict)
                print(f"  ✅ Precedent builder works with recommender output")
            else:
                print(f"  ⚠️  No precedents found in result")
                
        except Exception as e:
            print(f"  ⚠️  Integration failed: {type(e).__name__}: {e}")


def run_all_tests():
    """Run all edge case tests"""
    print("=" * 60)
    print("Edge Case and Stress Tests for Refactored Core.py")
    print("=" * 60)
    
    test_classes = [
        ("Edge Case Tests", EdgeCaseTests),
        ("Memory Tests", MemoryTests),
        ("Integration Tests", IntegrationTests),
    ]
    
    for test_name, test_class in test_classes:
        print(f"\n{'=' * 60}")
        print(f"{test_name}")
        print("=" * 60)
        
        instance = test_class()
        test_methods = [m for m in dir(instance) if m.startswith("test_")]
        
        for method_name in test_methods:
            try:
                method = getattr(instance, method_name)
                method()
            except Exception as e:
                print(f"  ❌ {method_name} crashed: {type(e).__name__}: {e}")
    
    print("\n" + "=" * 60)
    print("Edge case testing complete!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
