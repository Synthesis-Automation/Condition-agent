#!/usr/bin/env python
"""Check rxn_insight installation and test detection."""

try:
    import rxn_insight
    print("✓ rxn_insight INSTALLED")
    print(f"  Version: {getattr(rxn_insight, '__version__', 'unknown')}")
except ImportError as e:
    print("✗ rxn_insight NOT INSTALLED")
    print(f"  Error: {e}")
    exit(1)

try:
    from chemtools.reaction_type_detector import is_available, detect_reaction_type
    print(f"✓ Detection available: {is_available()}")
    
    # Test with a simple esterification
    test_rxn = "CC(=O)O.CCO>>CC(=O)OCC"
    result = detect_reaction_type(test_rxn)
    print(f"\nTest reaction: {test_rxn}")
    print(f"  Success: {result.get('success')}")
    print(f"  Class: {result.get('rxn_class')}")
    print(f"  Name: {result.get('rxn_name')}")
    print(f"  Mapped family: {result.get('mapped_family')}")
    print(f"  Confidence: {result.get('confidence')}")
    
except Exception as e:
    print(f"✗ Detection error: {e}")
    import traceback
    traceback.print_exc()
