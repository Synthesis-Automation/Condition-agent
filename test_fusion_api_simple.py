"""
Simple test script for fusion API endpoint - no user interaction.
"""

import requests
import json

API_BASE = "http://localhost:8000"

def test_health():
    """Test server health."""
    print("=" * 80)
    print("CHECKING SERVER HEALTH")
    print("=" * 80)
    try:
        response = requests.get(f"{API_BASE}/health", timeout=5)
        response.raise_for_status()
        print("✅ Server is running\n")
        return True
    except Exception as e:
        print(f"❌ Server not accessible: {e}\n")
        return False


def test_fusion_basic():
    """Test basic fusion endpoint functionality."""
    print("=" * 80)
    print("TEST 1: BASIC FUSION ENDPOINT")
    print("=" * 80)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    print(f"Reaction: {reaction}")
    print("Type: Buchwald-Hartwig C-N coupling\n")
    
    payload = {
        "reaction": reaction,
        "k": 50,
        "max_variants": 5
    }
    
    try:
        response = requests.post(f"{API_BASE}/api/v1/recommend/fusion", json=payload, timeout=30)
        response.raise_for_status()
        result = response.json()
        
        print("✅ Response received successfully")
        print(f"Response keys: {list(result.keys())}\n")
        
        # Check fusion metadata
        if 'fusion_meta' in result:
            print("✅ Fusion metadata present")
            meta = result['fusion_meta']
            
            if 'adaptive_weights' in meta:
                weights = meta['adaptive_weights']
                print("\nAdaptive Weights:")
                for key in ['alpha', 'beta', 'gamma', 'delta', 'α', 'β', 'γ', 'δ']:
                    if key in weights:
                        print(f"  {key}: {weights[key]:.3f}")
                        
                weight_sum = sum(weights.values())
                print(f"  Sum: {weight_sum:.3f} (expected ~1.0)")
                
            if 'evidence_summary' in meta:
                evidence = meta['evidence_summary']
                print("\nEvidence Summary:")
                for key, val in evidence.items():
                    if isinstance(val, float):
                        print(f"  {key}: {val:.3f}")
                    else:
                        print(f"  {key}: {val}")
                        
            if 'reasoning' in meta:
                print(f"\nReasoning: {len(meta['reasoning'])} items")
                for i, reason in enumerate(meta['reasoning'][:3], 1):
                    print(f"  {i}. {reason}")
        else:
            print("❌ Fusion metadata missing!")
            
        # Check recommendations
        if 'formatted' in result and 'recommended_conditions' in result['formatted']:
            recs = result['formatted']['recommended_conditions']
            print(f"\n✅ Generated {len(recs)} recommendations")
            
            if recs:
                print("\nTop Recommendation:")
                top = recs[0]
                summary = top.get('summary', {})
                print(f"  Rank: {summary.get('rank', 1)}")
                print(f"  Core: {summary.get('core', 'N/A')}")
                
                base = summary.get('base', 'N/A')
                if isinstance(base, dict):
                    base = base.get('name', 'N/A')
                print(f"  Base: {base}")
                
                solvent = summary.get('solvent', 'N/A')
                if isinstance(solvent, dict):
                    solvent = solvent.get('name', 'N/A')
                print(f"  Solvent: {solvent}")
        
        print("\n" + "=" * 80)
        print("TEST 1: PASSED ✅")
        print("=" * 80 + "\n")
        return True
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        print("\n" + "=" * 80)
        print("TEST 1: FAILED ❌")
        print("=" * 80 + "\n")
        return False


def test_fusion_vs_baseline():
    """Compare fusion vs baseline endpoints."""
    print("=" * 80)
    print("TEST 2: FUSION VS BASELINE")
    print("=" * 80)
    
    reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    print(f"Reaction: {reaction}")
    print("Type: Suzuki coupling\n")
    
    # Test baseline
    print("[1] Baseline endpoint (/api/v1/recommend)")
    print("-" * 80)
    try:
        response = requests.post(
            f"{API_BASE}/api/v1/recommend",
            json={"reaction": reaction, "k": 50},
            timeout=30
        )
        response.raise_for_status()
        baseline = response.json()
        
        has_fusion_meta = 'fusion_meta' in baseline
        print(f"  Fusion metadata present: {has_fusion_meta}")
        if has_fusion_meta:
            print("  ⚠️  WARNING: Baseline should not have fusion_meta")
        else:
            print("  ✅ Correct: No fusion_meta in baseline")
            
    except Exception as e:
        print(f"  ❌ Baseline request failed: {e}")
        return False
    
    # Test fusion
    print("\n[2] Fusion endpoint (/api/v1/recommend/fusion)")
    print("-" * 80)
    try:
        response = requests.post(
            f"{API_BASE}/api/v1/recommend/fusion",
            json={"reaction": reaction, "k": 50, "max_variants": 5},
            timeout=30
        )
        response.raise_for_status()
        fusion = response.json()
        
        has_fusion_meta = 'fusion_meta' in fusion
        print(f"  Fusion metadata present: {has_fusion_meta}")
        if has_fusion_meta:
            print("  ✅ Correct: fusion_meta present")
            
            # Show weights
            if 'adaptive_weights' in fusion['fusion_meta']:
                weights = fusion['fusion_meta']['adaptive_weights']
                print("\n  Adaptive Weights:")
                for key, val in weights.items():
                    print(f"    {key}: {val:.3f}")
        else:
            print("  ❌ ERROR: Fusion metadata missing")
            return False
            
    except Exception as e:
        print(f"  ❌ Fusion request failed: {e}")
        return False
    
    print("\n" + "=" * 80)
    print("TEST 2: PASSED ✅")
    print("=" * 80 + "\n")
    return True


def main():
    """Run all tests."""
    print("\n" + "=" * 80)
    print("FUSION API ENDPOINT TEST SUITE")
    print("=" * 80)
    print(f"API Base: {API_BASE}\n")
    
    # Check server
    if not test_health():
        print("\n❌ Server not running. Start with:")
        print("  uvicorn app.main:app --reload --port 8000\n")
        return False
    
    # Run tests
    results = []
    
    results.append(("Basic Fusion Endpoint", test_fusion_basic()))
    results.append(("Fusion vs Baseline Comparison", test_fusion_vs_baseline()))
    
    # Summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    print(f"\nResults: {passed}/{total} tests passed\n")
    
    for name, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"  {status}: {name}")
    
    print("\n" + "=" * 80)
    
    return all(success for _, success in results)


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
