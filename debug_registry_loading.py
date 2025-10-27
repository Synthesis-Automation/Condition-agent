"""
Test registry loading with detailed error reporting.
"""

import traceback
from chemtools.taxonomy import load_registry
from chemtools.analysis._registry import get_registry, clear_registry_cache


def test_registry_loading():
    """Test registry loading with full error details."""
    
    print("=" * 80)
    print("REGISTRY LOADING DEBUG TEST")
    print("=" * 80)
    print()
    
    # Clear any cached state
    clear_registry_cache()
    
    # Try direct loading
    print("Attempting direct load_registry()...")
    try:
        registry = load_registry()
        if registry:
            print("✓ SUCCESS - Registry loaded directly")
            print(f"  Type: {type(registry)}")
            print(f"  Has reactant_types: {hasattr(registry, 'reactant_types')}")
            if hasattr(registry, 'reactant_types'):
                print(f"  Reactant types count: {len(registry.reactant_types)}")
            print(f"  Has reaction_types: {hasattr(registry, 'reaction_types')}")
            if hasattr(registry, 'reaction_types'):
                print(f"  Reaction types count: {len(registry.reaction_types)}")
        else:
            print("✗ FAILED - load_registry() returned None")
    except Exception as e:
        print(f"✗ FAILED - Exception during load_registry():")
        print(f"  Error: {e}")
        print(f"  Type: {type(e).__name__}")
        traceback.print_exc()
    
    print()
    
    # Try via get_registry
    print("Attempting get_registry()...")
    clear_registry_cache()
    try:
        registry = get_registry()
        if registry:
            print("✓ SUCCESS - Registry loaded via get_registry()")
            print(f"  Type: {type(registry)}")
        else:
            print("✗ FAILED - get_registry() returned None")
            print("  This means load_registry() threw an exception (silently caught)")
    except Exception as e:
        print(f"✗ FAILED - Exception during get_registry():")
        print(f"  Error: {e}")
        traceback.print_exc()
    
    print()
    print("=" * 80)


if __name__ == "__main__":
    test_registry_loading()
