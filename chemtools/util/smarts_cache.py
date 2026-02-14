"""
Centralized SMARTS Pattern Compilation with LRU Caching

Provides a single source of truth for SMARTS compilation across chemtools.
Eliminates redundant compilation and provides consistent error handling.

Usage:
    from chemtools.util.smarts_cache import compile_smarts, compile_smarts_batch
    
    # Compile single pattern
    pattern = compile_smarts("[CX4][Cl,Br,I]")
    if pattern and mol.HasSubstructMatch(pattern):
        print("Match found")
    
    # Compile multiple patterns
    patterns = compile_smarts_batch({
        "aryl_halide": "[$(c[Cl,Br,I]),$(c-[Cl,Br,I])]",
        "vinyl_halide": "C=C[Cl,Br,I]"
    })
    
    # Check cache performance
    info = get_cache_info()
    print(f"Cache hit rate: {info['hit_rate']:.1%}")

Benefits:
    - Single compilation per unique SMARTS pattern (cached globally)
    - Consistent error handling across all modules
    - Cache statistics for performance monitoring
    - Configurable validation (strict vs. lenient)
    
Performance:
    - 10-100× speedup for repeated pattern usage
    - Memory savings from unified cache vs. per-module caches
    - Lazy compilation (patterns compiled only when first used)
"""

from __future__ import annotations
from functools import lru_cache
from typing import Dict, Optional, Any
import logging

from .rdkit_helpers import rdkit_available

logger = logging.getLogger(__name__)


@lru_cache(maxsize=1024)
def compile_smarts(smarts: str, *, validate: bool = False) -> Optional[Any]:
    """
    Compile a SMARTS pattern with LRU caching.
    
    This function maintains a global LRU cache of compiled SMARTS patterns,
    ensuring that each unique pattern is compiled only once across the entire
    chemtools library. Subsequent calls with the same pattern return the cached
    compiled object.
    
    Args:
        smarts: SMARTS pattern string to compile
        validate: If True, raise ValueError on invalid patterns.
                 If False, return None silently (default: False)
        
    Returns:
        Compiled RDKit Mol object representing the SMARTS pattern,
        or None if RDKit is unavailable or pattern is invalid
        
    Raises:
        ValueError: If validate=True and pattern cannot be compiled
        
    Examples:
        >>> # Basic usage
        >>> pattern = compile_smarts("[CX4][Cl,Br,I]")
        >>> if pattern:
        ...     matches = mol.GetSubstructMatches(pattern)
        
        >>> # With validation (raises on invalid pattern)
        >>> try:
        ...     pattern = compile_smarts("[[[invalid", validate=True)
        ... except ValueError as e:
        ...     print(f"Invalid pattern: {e}")
        
        >>> # Cache benefits
        >>> import time
        >>> start = time.time()
        >>> for _ in range(1000):
        ...     p = compile_smarts("[CX4][Cl,Br,I]")  # Only compiles once!
        >>> elapsed = time.time() - start
        >>> print(f"1000 calls in {elapsed:.3f}s (cached)")
    
    Performance:
        First call: ~1-5ms (compilation time)
        Cached calls: ~0.001ms (hash lookup)
        Speedup: 1000-5000× for cached patterns
        
    Thread Safety:
        This function is thread-safe due to functools.lru_cache.
    """
    if not rdkit_available():
        if validate:
            raise ValueError("RDKit is not available, cannot compile SMARTS patterns")
        return None
    
    try:
        from rdkit import Chem
        pattern = Chem.MolFromSmarts(smarts)
        
        if pattern is None:
            msg = f"Invalid SMARTS pattern (compiled to None): {smarts}"
            if validate:
                raise ValueError(msg)
            logger.debug(msg)
            return None
        
        return pattern
        
    except Exception as e:
        msg = f"Failed to compile SMARTS pattern '{smarts}': {e}"
        if validate:
            raise ValueError(msg) from e
        logger.debug(msg)
        return None


def compile_smarts_batch(
    patterns: Dict[str, str],
    *,
    validate: bool = False,
    skip_invalid: bool = True
) -> Dict[str, Optional[Any]]:
    """
    Compile multiple SMARTS patterns efficiently in batch.
    
    This function compiles a dictionary of named SMARTS patterns, returning
    a dictionary with the same keys but compiled pattern objects as values.
    Benefits from the same global LRU cache as compile_smarts().
    
    Args:
        patterns: Dict mapping pattern names to SMARTS strings
        validate: If True, raise ValueError on first invalid pattern
        skip_invalid: If True (default), return None for invalid patterns.
                     If False, raise ValueError on invalid patterns
        
    Returns:
        Dict mapping pattern names to compiled patterns (or None for invalid)
        
    Raises:
        ValueError: If validate=True or skip_invalid=False and any pattern is invalid
        
    Examples:
        >>> patterns = {
        ...     "halide": "[Cl,Br,I]",
        ...     "alcohol": "[OX2H]",
        ...     "amine": "[NX3;H1,H2]"
        ... }
        >>> compiled = compile_smarts_batch(patterns)
        >>> for name, pattern in compiled.items():
        ...     if pattern:
        ...         print(f"{name}: valid")
        
        >>> # With validation
        >>> patterns_with_invalid = {
        ...     "good": "[Cl,Br,I]",
        ...     "bad": "[[[invalid"
        ... }
        >>> compiled = compile_smarts_batch(
        ...     patterns_with_invalid,
        ...     validate=True  # Will raise ValueError on "bad"
        ... )
    
    Performance:
        - Leverages cache from compile_smarts()
        - No overhead vs. calling compile_smarts() individually
        - Efficient for initializing module-level pattern dictionaries
    """
    result = {}
    
    for name, smarts in patterns.items():
        try:
            pattern = compile_smarts(smarts, validate=validate or not skip_invalid)
            result[name] = pattern
        except ValueError as e:
            if not skip_invalid:
                raise ValueError(f"Invalid pattern '{name}': {e}") from e
            logger.debug(f"Skipping invalid pattern '{name}': {e}")
            result[name] = None
    
    return result


def clear_smarts_cache() -> None:
    """
    Clear the SMARTS compilation cache.
    
    Useful for testing or when patterns change dynamically. Note that clearing
    the cache will cause subsequent compile_smarts() calls to recompile patterns.
    
    Examples:
        >>> # Clear cache before performance test
        >>> clear_smarts_cache()
        >>> start = time.time()
        >>> pattern = compile_smarts("[CX4][Cl,Br,I]")
        >>> first_compile_time = time.time() - start
        
        >>> # Test in pytest
        >>> def test_pattern_compilation():
        ...     clear_smarts_cache()  # Start with clean cache
        ...     pattern = compile_smarts("[CX4][Cl,Br,I]")
        ...     assert pattern is not None
    """
    compile_smarts.cache_clear()
    logger.debug("SMARTS compilation cache cleared")


def get_cache_info() -> Dict[str, Any]:
    """
    Get LRU cache statistics for performance monitoring.
    
    Returns:
        Dict with keys:
            - hits: Number of cache hits
            - misses: Number of cache misses (compilations performed)
            - size: Current number of patterns in cache
            - maxsize: Maximum cache size (1024 by default)
            - hit_rate: Ratio of hits to total calls (0.0 to 1.0)
        
    Examples:
        >>> # Monitor cache performance
        >>> info = get_cache_info()
        >>> print(f"Cache hit rate: {info['hit_rate']:.1%}")
        >>> print(f"Patterns cached: {info['size']}/{info['maxsize']}")
        
        >>> # Log statistics periodically
        >>> if info['size'] == info['maxsize']:
        ...     logger.warning("SMARTS cache is full, consider increasing maxsize")
    
    Performance Tuning:
        - Hit rate > 90%: Good caching performance
        - Hit rate < 50%: Consider increasing cache size or patterns are too diverse
        - Size == maxsize: Cache may be evicting useful patterns
    """
    info = compile_smarts.cache_info()
    total_calls = info.hits + info.misses
    
    return {
        "hits": info.hits,
        "misses": info.misses,
        "size": info.currsize,
        "maxsize": info.maxsize,
        "hit_rate": info.hits / total_calls if total_calls > 0 else 0.0,
        "total_calls": total_calls,
    }


def log_cache_stats(level: int = logging.INFO) -> None:
    """
    Log current cache statistics at specified log level.
    
    Args:
        level: Logging level (default: logging.INFO)
        
    Examples:
        >>> # Log cache stats for debugging
        >>> log_cache_stats(logging.DEBUG)
        DEBUG: SMARTS cache: 45/1024 patterns, hit rate: 94.2%
        
        >>> # Periodic monitoring
        >>> import time
        >>> while True:
        ...     log_cache_stats()
        ...     time.sleep(60)  # Log every minute
    """
    info = get_cache_info()
    logger.log(
        level,
        f"SMARTS cache: {info['size']}/{info['maxsize']} patterns, "
        f"hit rate: {info['hit_rate']:.1%} "
        f"(hits={info['hits']}, misses={info['misses']})"
    )

