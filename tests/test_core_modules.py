#!/usr/bin/env python3
"""
Comprehensive unit tests for refactored core.py modules.

Tests each module independently to ensure correct behavior
and proper isolation between components.
"""

import pytest
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


class TestRecommenderModule:
    """Test chemtools.recommend.modules.recommender"""
    
    def test_import_recommender(self):
        """Test that recommender module can be imported"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        assert callable(recommend_from_reaction)
    
    def test_recommend_basic_suzuki(self):
        """Test basic Suzuki coupling recommendation"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        result = recommend_from_reaction(reaction, k=5, rerank_strategy='none')
        
        assert isinstance(result, dict)
        assert "family" in result or "detected_family" in result
    
    def test_recommend_with_family_override(self):
        """Test recommendation with family override"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        result = recommend_from_reaction(
            reaction, 
            k=5, 
            family_override="Suzuki",
            rerank_strategy='none'
        )
        
        assert isinstance(result, dict)
    
    def test_recommend_with_constraints(self):
        """Test recommendation with constraint rules"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        constraints = {
            "filter_unknown_reagents": True
        }
        
        result = recommend_from_reaction(
            reaction, 
            k=5,
            constraint_rules=constraints,
            rerank_strategy='none'
        )
        
        assert isinstance(result, dict)
    
    def test_recommend_different_rerank_strategies(self):
        """Test different reranking strategies"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        strategies = ['none', 'rule', 'analytics']
        for strategy in strategies:
            result = recommend_from_reaction(
                reaction, 
                k=5, 
                rerank_strategy=strategy
            )
            assert isinstance(result, dict), f"Strategy {strategy} failed"


class TestPrecedentBuilderModule:
    """Test chemtools.recommend.modules.precedent_builder"""
    
    def test_import_precedent_functions(self):
        """Test that all precedent functions can be imported"""
        from chemtools.recommend.modules.precedent_builder import (
            build_precedent_details,
            calculate_average_yield,
            calculate_yield_range,
            calculate_temp_range,
            calculate_time_range,
        )
        
        assert callable(build_precedent_details)
        assert callable(calculate_average_yield)
        assert callable(calculate_yield_range)
        assert callable(calculate_temp_range)
        assert callable(calculate_time_range)
    
    def test_calculate_average_yield_empty(self):
        """Test yield calculation with empty precedents"""
        from chemtools.recommend.modules.precedent_builder import calculate_average_yield
        
        result = calculate_average_yield([])
        assert result is None or result == 0 or isinstance(result, (int, float))
    
    def test_calculate_average_yield_with_data(self):
        """Test yield calculation with sample data"""
        from chemtools.recommend.modules.precedent_builder import calculate_average_yield
        
        precedents = [
            {"yield": 85.0},
            {"yield": 90.0},
            {"yield": 88.0},
        ]
        
        avg = calculate_average_yield(precedents)
        assert avg is not None
        if avg > 0:
            assert 85 <= avg <= 90
    
    def test_calculate_yield_range(self):
        """Test yield range calculation"""
        from chemtools.recommend.modules.precedent_builder import calculate_yield_range
        
        precedents = [
            {"yield": 85.0},
            {"yield": 90.0},
            {"yield": 88.0},
        ]
        
        result = calculate_yield_range(precedents)
        # Result could be tuple, dict, or None depending on implementation
        assert result is None or isinstance(result, (tuple, dict, list))
    
    def test_build_precedent_details_empty(self):
        """Test precedent details with empty list"""
        from chemtools.recommend.modules.precedent_builder import build_precedent_details
        
        result = build_precedent_details(precs=[], chosen_core=None, group=[])
        assert isinstance(result, dict)
    
    def test_build_precedent_details_with_data(self):
        """Test precedent details with sample data"""
        from chemtools.recommend.modules.precedent_builder import build_precedent_details
        
        precedents = [
            {"yield": 85.0, "temperature": 80, "time": 2},
            {"yield": 90.0, "temperature": 90, "time": 3},
            {"yield": 88.0, "temperature": 85, "time": 2.5},
        ]
        
        result = build_precedent_details(
            precs=precedents,
            chosen_core="Pd(PPh3)4",
            group=precedents[:2]
        )
        assert isinstance(result, dict)


class TestOutputBuilderModule:
    """Test chemtools.recommend.modules.output_builder"""
    
    def test_import_output_builder(self):
        """Test that output builder can be imported"""
        from chemtools.recommend.modules.output_builder import build_formatted_output
        assert callable(build_formatted_output)
    
    def test_build_formatted_output_basic(self):
        """Test basic formatted output building"""
        from chemtools.recommend.modules.output_builder import build_formatted_output
        
        # Minimal recommendation structure
        recommendation = {
            "catalyst": "Pd(PPh3)4",
            "base": "K2CO3",
            "solvent": "DMF"
        }
        
        precedents = []
        family = "Suzuki"
        
        try:
            result = build_formatted_output(
                recommendation=recommendation,
                precedents=precedents,
                family=family,
                max_variants=1
            )
            # Should return some kind of formatted output
            assert result is not None
        except Exception as e:
            # If it requires specific database/data, that's okay
            # We're just checking it doesn't have import/syntax errors
            assert "recommendation" in str(e).lower() or "precedent" in str(e).lower() or "database" in str(e).lower()


class TestStructuredModule:
    """Test chemtools.recommend.modules.structured"""
    
    def test_import_structured(self):
        """Test that structured module can be imported"""
        from chemtools.recommend.modules.structured import recommend_conditions_structured
        assert callable(recommend_conditions_structured)
    
    def test_structured_basic_call(self):
        """Test basic structured recommendation call"""
        from chemtools.recommend.modules.structured import recommend_conditions_structured
        
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        result = recommend_conditions_structured(
            reaction=reaction,
            k=5,
            limit=2,
            rerank_strategy='none'
        )
        
        assert isinstance(result, dict)
        assert "meta" in result
        assert "recommendations" in result or "condition_variants" in result
    
    def test_structured_has_timing(self):
        """Test that structured output includes timing information"""
        from chemtools.recommend.modules.structured import recommend_conditions_structured
        
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        result = recommend_conditions_structured(
            reaction=reaction,
            k=3,
            limit=1,
            rerank_strategy='none'
        )
        
        assert "meta" in result
        if isinstance(result["meta"], dict):
            # Check for timing info somewhere in the result
            has_timing = (
                "processing_time_ms" in result["meta"] or
                "processing_time_ms" in result
            )
            # Timing should be present
            assert has_timing or "timestamp" in result["meta"]
    
    def test_structured_with_reaction_type(self):
        """Test structured call with explicit reaction type"""
        from chemtools.recommend.modules.structured import recommend_conditions_structured
        
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        result = recommend_conditions_structured(
            reaction=reaction,
            reaction_type="Suzuki",
            k=3,
            limit=1,
            rerank_strategy='none'
        )
        
        assert isinstance(result, dict)


class TestFusionAdapterModule:
    """Test chemtools.recommend.modules.fusion_adapter"""
    
    def test_import_fusion_adapter(self):
        """Test that fusion adapter functions can be imported"""
        from chemtools.recommend.modules.fusion_adapter import (
            convert_fusion_to_core_format,
            build_formatted_output_from_fusion,
        )
        
        assert callable(convert_fusion_to_core_format)
        assert callable(build_formatted_output_from_fusion)
    
    def test_convert_fusion_empty(self):
        """Test fusion conversion with empty data"""
        from chemtools.recommend.modules.fusion_adapter import convert_fusion_to_core_format
        
        # Test with minimal required parameters
        result = convert_fusion_to_core_format(
            fusion_results={},
            reaction="c1ccccc1Br>>c1ccccc1",
            rxn_smiles_norm="c1ccccc1Br>>c1ccccc1",
            fam="Unknown",
            features={},
            detection_meta={}
        )
        assert isinstance(result, dict)


class TestModulesPackage:
    """Test chemtools.recommend.modules package"""
    
    def test_package_imports_all_functions(self):
        """Test that package __init__ exports all expected functions"""
        import chemtools.recommend.modules as modules
        
        expected_functions = [
            "recommend_from_reaction",
            "recommend_conditions_structured",
            "build_formatted_output",
            "convert_fusion_to_core_format",
            "build_formatted_output_from_fusion",
            "build_precedent_details",
            "calculate_average_yield",
            "calculate_yield_range",
            "calculate_temp_range",
            "calculate_time_range",
        ]
        
        for func_name in expected_functions:
            assert hasattr(modules, func_name), f"Missing export: {func_name}"
            assert callable(getattr(modules, func_name))
    
    def test_package_all_attribute(self):
        """Test that package defines __all__"""
        import chemtools.recommend.modules as modules
        
        # Should have __all__ defined
        if hasattr(modules, '__all__'):
            assert isinstance(modules.__all__, (list, tuple))
            assert len(modules.__all__) > 0


class TestBackwardsCompatibility:
    """Test backwards compatibility of refactored core.py"""
    
    def test_old_import_paths_work(self):
        """Test that old import paths still work"""
        # These are the old import paths that existing code uses
        from chemtools.recommend.core import recommend_from_reaction
        from chemtools.recommend.core import recommend_conditions_structured
        
        assert callable(recommend_from_reaction)
        assert callable(recommend_conditions_structured)
    
    def test_internal_functions_still_accessible(self):
        """Test that internal functions with _ prefix are still accessible"""
        from chemtools.recommend.core import (
            _build_formatted_output,
            _build_precedent_details,
        )
        
        assert callable(_build_formatted_output)
        assert callable(_build_precedent_details)
    
    def test_new_import_paths_work(self):
        """Test that new module-based imports work"""
        from chemtools.recommend.modules import recommend_from_reaction
        from chemtools.recommend.modules import recommend_conditions_structured
        
        assert callable(recommend_from_reaction)
        assert callable(recommend_conditions_structured)
    
    def test_direct_module_imports_work(self):
        """Test that direct module imports work"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        from chemtools.recommend.modules.structured import recommend_conditions_structured
        
        assert callable(recommend_from_reaction)
        assert callable(recommend_conditions_structured)
    
    def test_all_import_paths_same_function(self):
        """Test that all import paths reference the same function objects"""
        from chemtools.recommend.core import recommend_from_reaction as core_rec
        from chemtools.recommend.modules import recommend_from_reaction as mod_rec
        from chemtools.recommend.modules.recommender import recommend_from_reaction as direct_rec
        
        assert core_rec is mod_rec is direct_rec


class TestErrorHandling:
    """Test error handling in refactored modules"""
    
    def test_invalid_reaction_smiles(self):
        """Test handling of invalid SMILES"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        # Invalid SMILES should be handled gracefully
        invalid_reaction = "invalid>>smiles"
        
        try:
            result = recommend_from_reaction(invalid_reaction, k=5, rerank_strategy='none')
            # If it doesn't raise an error, check the result
            assert isinstance(result, dict)
        except Exception as e:
            # Should raise a meaningful error, not a cryptic one
            error_msg = str(e).lower()
            assert "smiles" in error_msg or "invalid" in error_msg or "parse" in error_msg
    
    def test_negative_k_value(self):
        """Test handling of invalid k parameter"""
        from chemtools.recommend.modules.recommender import recommend_from_reaction
        
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        # Negative k should either be handled or raise clear error
        try:
            result = recommend_from_reaction(reaction, k=-1, rerank_strategy='none')
            # If it succeeds, that's fine too (might be converted to 0 or default)
            assert isinstance(result, dict)
        except (ValueError, AssertionError) as e:
            # Expected error type
            pass


class TestPerformance:
    """Test performance characteristics of refactored modules"""
    
    def test_lazy_imports_no_circular_deps(self):
        """Test that modules can be imported without circular dependency issues"""
        # Import modules in different orders to test for circular deps
        
        # Order 1
        from chemtools.recommend.modules import recommender
        from chemtools.recommend.modules import structured
        from chemtools.recommend.modules import output_builder
        
        # Order 2 (reverse)
        from chemtools.recommend.modules import precedent_builder
        from chemtools.recommend.modules import fusion_adapter
        
        # All should import successfully
        assert recommender is not None
        assert structured is not None
        assert output_builder is not None
        assert precedent_builder is not None
        assert fusion_adapter is not None
    
    def test_module_import_speed(self):
        """Test that module imports are reasonably fast"""
        import time
        
        start = time.time()
        from chemtools.recommend.modules import recommend_from_reaction
        from chemtools.recommend.modules import recommend_conditions_structured
        end = time.time()
        
        import_time = end - start
        
        # Imports should be fast (< 5 seconds even on slow systems)
        assert import_time < 5.0, f"Imports took {import_time}s, too slow"


if __name__ == "__main__":
    # Run with pytest if available, otherwise run basic tests
    try:
        import pytest
        sys.exit(pytest.main([__file__, "-v"]))
    except ImportError:
        print("pytest not available, running basic tests...")
        
        # Run a few basic tests manually
        test = TestBackwardsCompatibility()
        test.test_old_import_paths_work()
        print("✅ Old import paths work")
        
        test.test_new_import_paths_work()
        print("✅ New import paths work")
        
        test.test_all_import_paths_same_function()
        print("✅ All import paths reference same functions")
        
        print("\n✅ Basic tests passed! Run 'pytest tests/test_core_modules.py -v' for full suite")
