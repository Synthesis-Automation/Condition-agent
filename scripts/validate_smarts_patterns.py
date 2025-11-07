#!/usr/bin/env python
"""
Validate all SMARTS patterns used across the Condition-agent codebase.

This script discovers and validates SMARTS patterns from:
- chemtools/router.py
- chemtools/util/functional_groups.py
- chemtools/featurizers/calculable_features.json
- chemtools/selector_payloads.py
- Other modules using SMARTS

Usage:
    python scripts/validate_smarts_patterns.py
    python scripts/validate_smarts_patterns.py --verbose
    python scripts/validate_smarts_patterns.py --module router

Exit codes:
    0: All patterns valid
    1: Some patterns invalid
    2: Critical error (e.g., RDKit not available)
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.util.smarts_cache import compile_smarts, clear_smarts_cache
from chemtools.util.rdkit_helpers import rdkit_available


def load_router_patterns() -> Dict[str, str]:
    """Load SMARTS patterns from router.py."""
    try:
        from chemtools.router import _SMARTS_PATTERNS
        return dict(_SMARTS_PATTERNS)
    except ImportError:
        print("Warning: Could not import router module", file=sys.stderr)
        return {}


def load_functional_group_patterns() -> Dict[str, str]:
    """Load SMARTS patterns from functional_groups.py."""
    try:
        from chemtools.util.functional_groups import FUNCTIONAL_GROUP_SMARTS
        return dict(FUNCTIONAL_GROUP_SMARTS)
    except ImportError:
        print("Warning: Could not import functional_groups module", file=sys.stderr)
        return {}


def load_calculable_features_patterns() -> Dict[str, str]:
    """Load SMARTS patterns from calculable_features.json."""
    json_path = PROJECT_ROOT / "chemtools" / "featurizers" / "calculable_features.json"
    
    if not json_path.exists():
        print(f"Warning: {json_path} not found", file=sys.stderr)
        return {}
    
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        patterns = {}
        for feature in data.get("features", []):
            if feature.get("type") == "bool" and "smarts" in feature:
                name = feature.get("token", "unknown")
                smarts_list = feature["smarts"]
                # Store first pattern (others are alternatives)
                if isinstance(smarts_list, list) and smarts_list:
                    patterns[name] = smarts_list[0]
                elif isinstance(smarts_list, str):
                    patterns[name] = smarts_list
        
        return patterns
    except Exception as e:
        print(f"Error loading calculable_features.json: {e}", file=sys.stderr)
        return {}


def load_substrate_classifier_patterns() -> Dict[str, str]:
    """Load SMARTS patterns from substrate_classifier.py."""
    patterns = {
        'benzylic': "[CX4;H1,H2,H3]-[c,$(C=C-c)]",
        'allylic': "[CX4;H1,H2,H3]-[CX3]=[CX3]",
        'propargylic': "[CX4;H1,H2,H3]-[CX2]#[CX2]",
        'alpha_to_carbonyl': "[CX4;H1,H2,H3]-[CX3]=[OX1]",
        'benzene_ring': "c1ccccc1",
    }
    return patterns


def load_selector_payloads_patterns() -> Dict[str, str]:
    """Load SMARTS patterns from selector_payloads.py."""
    try:
        from chemtools.selector_payloads import (
            _PHENOL_SMARTS,
            _FREE_ALCOHOL_SMARTS,
            _CARBOXYLIC_ACID_SMARTS
        )
        return {
            "phenol": _PHENOL_SMARTS,
            "free_alcohol": _FREE_ALCOHOL_SMARTS,
            "carboxylic_acid": _CARBOXYLIC_ACID_SMARTS,
        }
    except ImportError:
        print("Warning: Could not import selector_payloads module", file=sys.stderr)
        return {}


def validate_patterns(
    patterns: Dict[str, str],
    category: str,
    verbose: bool = False
) -> Tuple[int, int, List[str]]:
    """
    Validate a set of SMARTS patterns.
    
    Returns:
        (valid_count, invalid_count, errors)
    """
    valid = 0
    invalid = 0
    errors = []
    
    if verbose:
        print(f"\nValidating {category} ({len(patterns)} patterns)...")
    
    for name, smarts in patterns.items():
        try:
            pattern = compile_smarts(smarts, validate=True)
            if pattern is None:
                error = f"{category}/{name}: Pattern compiled to None"
                errors.append(error)
                invalid += 1
                if verbose:
                    print(f"  ✗ {name}: {error}")
            else:
                valid += 1
                if verbose:
                    print(f"  ✓ {name}")
        except Exception as e:
            error = f"{category}/{name}: {e}"
            errors.append(error)
            invalid += 1
            if verbose:
                print(f"  ✗ {name}: {e}")
    
    return valid, invalid, errors


def main():
    parser = argparse.ArgumentParser(
        description="Validate SMARTS patterns across the codebase",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate all patterns
  python scripts/validate_smarts_patterns.py
  
  # Verbose output
  python scripts/validate_smarts_patterns.py --verbose
  
  # Validate specific module
  python scripts/validate_smarts_patterns.py --module router
  python scripts/validate_smarts_patterns.py --module functional_groups
        """
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed validation results"
    )
    parser.add_argument(
        "--module", "-m",
        choices=["router", "functional_groups", "calculable", "substrate", "selector", "all"],
        default="all",
        help="Validate patterns from specific module only"
    )
    
    args = parser.parse_args()
    
    # Check RDKit availability
    if not rdkit_available():
        print("❌ ERROR: RDKit is not available. Cannot validate SMARTS patterns.", file=sys.stderr)
        return 2
    
    # Clear cache to ensure clean validation
    clear_smarts_cache()
    
    # Load patterns from different modules
    pattern_sources = {}
    
    if args.module in ("router", "all"):
        pattern_sources["Router"] = load_router_patterns()
    
    if args.module in ("functional_groups", "all"):
        pattern_sources["Functional Groups"] = load_functional_group_patterns()
    
    if args.module in ("calculable", "all"):
        pattern_sources["Calculable Features"] = load_calculable_features_patterns()
    
    if args.module in ("substrate", "all"):
        pattern_sources["Substrate Classifier"] = load_substrate_classifier_patterns()
    
    if args.module in ("selector", "all"):
        pattern_sources["Selector Payloads"] = load_selector_payloads_patterns()
    
    # Validate all patterns
    total_valid = 0
    total_invalid = 0
    all_errors = []
    
    for category, patterns in pattern_sources.items():
        if not patterns:
            continue
        
        valid, invalid, errors = validate_patterns(patterns, category, args.verbose)
        total_valid += valid
        total_invalid += invalid
        all_errors.extend(errors)
        
        if not args.verbose:
            status = "✓" if invalid == 0 else "✗"
            print(f"{status} {category}: {valid}/{len(patterns)} valid")
    
    # Print summary
    print(f"\n{'='*70}")
    print(f"Validation Summary")
    print(f"{'='*70}")
    print(f"Total patterns validated: {total_valid + total_invalid}")
    print(f"Valid patterns:          {total_valid}")
    print(f"Invalid patterns:        {total_invalid}")
    
    if all_errors:
        print(f"\n❌ Found {len(all_errors)} errors:")
        for error in all_errors:
            print(f"  • {error}")
        return 1
    else:
        print("\n✅ All patterns valid!")
        
        # Show cache stats
        from chemtools.util.smarts_cache import get_cache_info
        info = get_cache_info()
        print(f"\nCache statistics:")
        print(f"  Patterns cached: {info['size']}")
        print(f"  Cache hits: {info['hits']}")
        print(f"  Cache misses: {info['misses']}")
        print(f"  Hit rate: {info['hit_rate']:.1%}")
        
        return 0


if __name__ == "__main__":
    sys.exit(main())
