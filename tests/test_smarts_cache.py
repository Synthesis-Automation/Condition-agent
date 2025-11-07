"""
Tests for centralized SMARTS compilation cache.

Tests cover:
- Basic compilation and caching behavior
- Batch compilation
- Cache statistics and monitoring
- Performance improvements
- Error handling and validation
"""

import pytest
import time
from chemtools.util.smarts_cache import (
    compile_smarts,
    compile_smarts_batch,
    clear_smarts_cache,
    get_cache_info,
    log_cache_stats,
)
from chemtools.util.rdkit_helpers import rdkit_available


# Skip all tests if RDKit not available
pytestmark = pytest.mark.skipif(
    not rdkit_available(),
    reason="RDKit required for SMARTS compilation tests"
)


class TestBasicCompilation:
    """Test basic SMARTS compilation functionality."""
    
    def test_compile_valid_pattern(self):
        """Test compilation of valid SMARTS pattern."""
        clear_smarts_cache()
        pattern = compile_smarts("[CX4][Cl,Br,I]")
        assert pattern is not None
    
    def test_compile_invalid_pattern_silent(self):
        """Test that invalid patterns return None by default."""
        clear_smarts_cache()
        pattern = compile_smarts("[[[invalid", validate=False)
        assert pattern is None
    
    def test_compile_invalid_pattern_validate(self):
        """Test that validation mode raises on invalid patterns."""
        clear_smarts_cache()
        with pytest.raises(ValueError, match="Invalid SMARTS"):
            compile_smarts("[[[invalid", validate=True)
    
    def test_compile_various_patterns(self):
        """Test compilation of various common SMARTS patterns."""
        clear_smarts_cache()
        
        patterns = [
            "[CX4][Cl,Br,I]",  # Alkyl halide
            "[OX2H]",           # Alcohol
            "[NX3;H1,H2]",      # Amine
            "c[Cl,Br,I]",       # Aryl halide
            "[CX3]=O",          # Carbonyl
            "C=C",              # Alkene
            "C#C",              # Alkyne
            "[BX3]",            # Boron
        ]
        
        for smarts in patterns:
            pattern = compile_smarts(smarts)
            assert pattern is not None, f"Failed to compile: {smarts}"


class TestCaching:
    """Test LRU caching behavior."""
    
    def test_cache_reuses_patterns(self):
        """Test that repeated calls return the same cached object."""
        clear_smarts_cache()
        
        pattern1 = compile_smarts("[CX4][Cl,Br,I]")
        pattern2 = compile_smarts("[CX4][Cl,Br,I]")
        
        # Should be the exact same object (cached)
        assert pattern1 is pattern2
    
    def test_cache_miss_then_hit(self):
        """Test cache statistics: miss on first call, hit on second."""
        clear_smarts_cache()
        
        # First call: cache miss
        compile_smarts("[CX4][Cl,Br,I]")
        info_after_first = get_cache_info()
        assert info_after_first["misses"] == 1
        assert info_after_first["hits"] == 0
        
        # Second call: cache hit
        compile_smarts("[CX4][Cl,Br,I]")
        info_after_second = get_cache_info()
        assert info_after_second["hits"] == 1
        assert info_after_second["misses"] == 1
    
    def test_cache_different_patterns(self):
        """Test that different patterns are cached separately."""
        clear_smarts_cache()
        
        pattern1 = compile_smarts("[Cl]")
        pattern2 = compile_smarts("[Br]")
        pattern3 = compile_smarts("[I]")
        
        # All should be different objects
        assert pattern1 is not pattern2
        assert pattern2 is not pattern3
        assert pattern1 is not pattern3
        
        # But cache size should be 3
        info = get_cache_info()
        assert info["size"] == 3
    
    def test_clear_cache(self):
        """Test that clearing cache resets statistics."""
        clear_smarts_cache()
        
        # Compile some patterns
        compile_smarts("[Cl]")
        compile_smarts("[Br]")
        compile_smarts("[I]")
        
        info_before = get_cache_info()
        assert info_before["size"] > 0
        
        # Clear and check
        clear_smarts_cache()
        info_after = get_cache_info()
        
        assert info_after["size"] == 0
        assert info_after["hits"] == 0
        assert info_after["misses"] == 0


class TestBatchCompilation:
    """Test batch compilation functionality."""
    
    def test_compile_batch_all_valid(self):
        """Test batch compilation of all valid patterns."""
        clear_smarts_cache()
        
        patterns = {
            "halide": "[Cl,Br,I]",
            "alcohol": "[OX2H]",
            "amine": "[NX3;H1,H2]",
        }
        
        compiled = compile_smarts_batch(patterns)
        
        assert len(compiled) == 3
        assert all(p is not None for p in compiled.values())
    
    def test_compile_batch_with_invalid(self):
        """Test batch compilation with some invalid patterns."""
        clear_smarts_cache()
        
        patterns = {
            "valid1": "[Cl,Br,I]",
            "invalid": "[[[bad",
            "valid2": "[OX2H]",
        }
        
        # With skip_invalid=True (default), should return None for invalid
        compiled = compile_smarts_batch(patterns, skip_invalid=True)
        
        assert compiled["valid1"] is not None
        assert compiled["invalid"] is None
        assert compiled["valid2"] is not None
    
    def test_compile_batch_validate_fail(self):
        """Test that validation mode with skip_invalid=False raises on invalid patterns."""
        clear_smarts_cache()
        
        patterns = {
            "valid": "[Cl,Br,I]",
            "invalid": "[[[bad",
        }
        
        # validate=True alone won't raise because skip_invalid=True by default
        # Need skip_invalid=False to actually raise
        with pytest.raises(ValueError):
            compile_smarts_batch(patterns, validate=True, skip_invalid=False)
    
    def test_compile_batch_no_skip(self):
        """Test batch compilation with skip_invalid=False."""
        clear_smarts_cache()
        
        patterns = {
            "valid": "[Cl,Br,I]",
            "invalid": "[[[bad",
        }
        
        with pytest.raises(ValueError):
            compile_smarts_batch(patterns, skip_invalid=False)


class TestCacheStatistics:
    """Test cache statistics and monitoring."""
    
    def test_get_cache_info_structure(self):
        """Test that get_cache_info returns expected structure."""
        clear_smarts_cache()
        
        info = get_cache_info()
        
        assert "hits" in info
        assert "misses" in info
        assert "size" in info
        assert "maxsize" in info
        assert "hit_rate" in info
        assert "total_calls" in info
        
        assert isinstance(info["hits"], int)
        assert isinstance(info["misses"], int)
        assert isinstance(info["size"], int)
        assert isinstance(info["maxsize"], int)
        assert isinstance(info["hit_rate"], float)
    
    def test_hit_rate_calculation(self):
        """Test hit rate calculation."""
        clear_smarts_cache()
        
        # 1 unique pattern called 3 times: 1 miss, 2 hits
        compile_smarts("[Cl]")  # miss
        compile_smarts("[Cl]")  # hit
        compile_smarts("[Cl]")  # hit
        
        info = get_cache_info()
        
        assert info["hits"] == 2
        assert info["misses"] == 1
        assert info["total_calls"] == 3
        assert info["hit_rate"] == pytest.approx(2/3)
    
    def test_cache_size_tracking(self):
        """Test that cache size is tracked correctly."""
        clear_smarts_cache()
        
        patterns = ["[Cl]", "[Br]", "[I]", "[F]"]
        
        for i, pattern in enumerate(patterns, 1):
            compile_smarts(pattern)
            info = get_cache_info()
            assert info["size"] == i
    
    def test_log_cache_stats(self, caplog):
        """Test logging of cache statistics."""
        clear_smarts_cache()
        
        compile_smarts("[Cl]")
        compile_smarts("[Cl]")
        
        import logging
        with caplog.at_level(logging.INFO):
            log_cache_stats(logging.INFO)
        
        assert "SMARTS cache:" in caplog.text
        assert "hit rate:" in caplog.text


class TestPerformance:
    """Test performance improvements from caching."""
    
    def test_cache_speedup(self):
        """Test that caching provides significant speedup."""
        clear_smarts_cache()
        
        smarts = "[CX4;H1,H2,H3]-[c,$(C=C-c)]"
        
        # Measure cached performance (1 miss + 999 hits)
        start = time.perf_counter()
        for _ in range(1000):
            compile_smarts(smarts)
        cached_time = time.perf_counter() - start
        
        # Measure uncached performance (1000 misses)
        start = time.perf_counter()
        for _ in range(1000):
            clear_smarts_cache()
            compile_smarts(smarts)
        uncached_time = time.perf_counter() - start
        
        speedup = uncached_time / cached_time
        
        # Expect at least 10× speedup from caching
        assert speedup > 10, f"Expected >10× speedup, got {speedup:.1f}×"
    
    def test_batch_compilation_performance(self):
        """Test that batch compilation is efficient."""
        clear_smarts_cache()
        
        patterns = {f"pattern_{i}": f"[{chr(ord('A') + i % 26)}]" 
                   for i in range(50)}
        
        # Batch compilation should be fast
        start = time.perf_counter()
        compiled = compile_smarts_batch(patterns)
        batch_time = time.perf_counter() - start
        
        assert len(compiled) == 50
        assert batch_time < 1.0  # Should complete in under 1 second


class TestRealWorldPatterns:
    """Test with real-world SMARTS patterns from the codebase."""
    
    def test_router_patterns(self):
        """Test patterns from router.py."""
        clear_smarts_cache()
        
        router_patterns = {
            "aryl_halide": "[$(c[Cl,Br,I]),$(c-[Cl,Br,I])]",
            "vinyl_halide": "C=C[Cl,Br,I]",
            "triflate": "OS(=O)(=O)C(F)(F)F",
            "boron": "[BX3;$(B(O)O),$(B(O)O),$(B(O)O)]",
            "terminal_alkyne": "C#C[H]",
            "nucleophile_n": "[NX3;H1,H2]",
            "carbonyl": "[CX3]=O",
        }
        
        compiled = compile_smarts_batch(router_patterns)
        
        # All should compile successfully
        assert all(p is not None for p in compiled.values())
    
    def test_substrate_classifier_patterns(self):
        """Test patterns from substrate_classifier.py."""
        clear_smarts_cache()
        
        special_position_patterns = {
            'benzylic': "[CX4;H1,H2,H3]-[c,$(C=C-c)]",
            'allylic': "[CX4;H1,H2,H3]-[CX3]=[CX3]",
            'propargylic': "[CX4;H1,H2,H3]-[CX2]#[CX2]",
            'alpha_to_carbonyl': "[CX4;H1,H2,H3]-[CX3]=[OX1]",
        }
        
        compiled = compile_smarts_batch(special_position_patterns)
        
        # All should compile successfully
        assert all(p is not None for p in compiled.values())
    
    def test_functional_groups_sample(self):
        """Test sample patterns from functional_groups.py."""
        clear_smarts_cache()
        
        fg_patterns = {
            "alcohol": "[OX2H]",
            "phenol": "c[OX2H]",
            "carbonyl": "[CX3]=[OX1]",
            "amine_primary": "[NX3;H2;!$(NC=O)]",
            "carboxylic_acid": "[CX3](=O)[OX2H1]",
        }
        
        compiled = compile_smarts_batch(fg_patterns)
        
        # All should compile successfully
        assert all(p is not None for p in compiled.values())


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_empty_string(self):
        """Test compilation of empty string.
        
        Note: RDKit actually returns a valid Mol object for empty string,
        representing a molecule with no atoms (universal match).
        """
        clear_smarts_cache()
        pattern = compile_smarts("", validate=False)
        # RDKit returns a valid Mol object for empty string
        assert pattern is not None
        assert hasattr(pattern, 'GetNumAtoms')
        assert pattern.GetNumAtoms() == 0
    
    def test_whitespace_only(self):
        """Test compilation of whitespace-only string."""
        clear_smarts_cache()
        pattern = compile_smarts("   ", validate=False)
        assert pattern is None
    
    def test_very_long_pattern(self):
        """Test compilation of very long pattern."""
        clear_smarts_cache()
        
        # Create a long but valid pattern
        long_pattern = "[" + ",".join([f"C{i}" for i in range(100)]) + "]"
        pattern = compile_smarts(long_pattern, validate=False)
        
        # May or may not compile depending on RDKit limits, but shouldn't crash
        assert pattern is None or pattern is not None
    
    def test_unicode_in_pattern(self):
        """Test that unicode in patterns is handled gracefully."""
        clear_smarts_cache()
        pattern = compile_smarts("[Cα]", validate=False)  # Greek alpha
        # Should return None (invalid) without crashing
        assert pattern is None


class TestThreadSafety:
    """Test thread safety of cache (functools.lru_cache is thread-safe)."""
    
    def test_concurrent_compilation(self):
        """Test that concurrent calls don't cause issues."""
        import threading
        
        clear_smarts_cache()
        results = []
        
        def compile_worker(pattern):
            result = compile_smarts(pattern)
            results.append(result)
        
        # Create multiple threads compiling same pattern
        threads = []
        for _ in range(10):
            t = threading.Thread(target=compile_worker, args=("[Cl,Br,I]",))
            threads.append(t)
            t.start()
        
        for t in threads:
            t.join()
        
        # All should succeed and return the same cached object
        assert len(results) == 10
        assert all(r is not None for r in results)
        assert all(r is results[0] for r in results)  # Same cached object


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
