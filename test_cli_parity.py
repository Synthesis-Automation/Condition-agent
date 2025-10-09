"""
Test that web and local CLIs have feature parity.
"""

import sys
from pathlib import Path

# Add scripts directory to path
ROOT = Path(__file__).parent
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

# Test imports
try:
    import local_recommendation_cli
    import web_recommendation_cli
    
    print("Testing CLI Feature Parity")
    print("=" * 80)
    
    # Check local_ml_recommendation signatures
    import inspect
    
    local_sig = inspect.signature(local_recommendation_cli.local_ml_recommendation)
    print("\nLocal CLI - local_ml_recommendation() parameters:")
    for param_name, param in local_sig.parameters.items():
        default = f" = {param.default}" if param.default != inspect.Parameter.empty else ""
        print(f"  - {param_name}: {param.annotation.__name__ if hasattr(param.annotation, '__name__') else param.annotation}{default}")
    
    web_sig = inspect.signature(web_recommendation_cli.call_ml_recommendation)
    print("\nWeb CLI - call_ml_recommendation() parameters:")
    for param_name, param in web_sig.parameters.items():
        default = f" = {param.default}" if param.default != inspect.Parameter.empty else ""
        print(f"  - {param_name}: {param.annotation.__name__ if hasattr(param.annotation, '__name__') else param.annotation}{default}")
    
    # Check for rerank_strategy
    local_has_rerank = 'rerank_strategy' in local_sig.parameters
    web_has_rerank = 'rerank_strategy' in web_sig.parameters
    
    # Check for filter_unknown_reagents  
    local_has_filter = 'filter_unknown_reagents' in local_sig.parameters
    web_has_filter = 'filter_unknown_reagents' in web_sig.parameters
    
    print("\n" + "=" * 80)
    print("Feature Comparison:")
    print("-" * 80)
    print(f"rerank_strategy parameter:")
    print(f"  Local CLI: {'✅ YES' if local_has_rerank else '❌ NO'}")
    print(f"  Web CLI:   {'✅ YES' if web_has_rerank else '❌ NO'}")
    print()
    print(f"filter_unknown_reagents parameter:")
    print(f"  Local CLI: {'✅ YES' if local_has_filter else '❌ NO'}")
    print(f"  Web CLI:   {'✅ YES' if web_has_filter else '❌ NO'}")
    print()
    
    # Overall assessment
    if local_has_rerank and web_has_rerank and local_has_filter and web_has_filter:
        print("=" * 80)
        print("✅ FEATURE PARITY ACHIEVED!")
        print("Both CLIs now support the same parameters.")
        print("=" * 80)
    else:
        print("=" * 80)
        print("⚠️  FEATURE MISMATCH DETECTED")
        print("CLIs have different parameters.")
        print("=" * 80)
        
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()
