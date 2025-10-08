"""
Test script for FastAPI fusion recommendation endpoint.

Tests the /api/v1/recommend/fusion endpoint with various reactions.
"""

import requests
import json
from typing import Dict, Any

API_BASE = "http://localhost:8000"
FUSION_ENDPOINT = f"{API_BASE}/api/v1/recommend/fusion"

def print_section(title: str):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def test_fusion_endpoint(reaction: str, description: str, k: int = 50, max_variants: int = 5):
    """Test the fusion endpoint with a single reaction."""
    print_section(f"Testing: {description}")
    print(f"Reaction: {reaction}")
    print(f"Parameters: k={k}, max_variants={max_variants}\n")
    
    payload = {
        "reaction": reaction,
        "k": k,
        "max_variants": max_variants
    }
    
    try:
        response = requests.post(FUSION_ENDPOINT, json=payload, timeout=30)
        response.raise_for_status()
        
        result = response.json()
        
        # Check response structure
        print("✅ Response received successfully\n")
        print("Response Keys:", list(result.keys()))
        
        # Check for fusion metadata
        if 'fusion_meta' in result:
            print("\n✅ Fusion metadata present")
            meta = result['fusion_meta']
            
            # Adaptive weights
            if 'adaptive_weights' in meta:
                weights = meta['adaptive_weights']
                print(f"\nAdaptive Weights:")
                print(f"  α (precedent): {weights.get('alpha', weights.get('α', 0)):.3f}")
                print(f"  β (analytics): {weights.get('beta', weights.get('β', 0)):.3f}")
                print(f"  γ (rules):     {weights.get('gamma', weights.get('γ', 0)):.3f}")
                print(f"  δ (ML):        {weights.get('delta', weights.get('δ', 0)):.3f}")
                
                # Verify weights sum to 1
                weight_sum = sum(weights.values())
                if 0.98 <= weight_sum <= 1.02:
                    print(f"  ✅ Weights sum: {weight_sum:.3f}")
                else:
                    print(f"  ⚠️  Weights sum: {weight_sum:.3f} (expected ~1.0)")
            
            # Evidence summary
            if 'evidence_summary' in meta:
                evidence = meta['evidence_summary']
                print(f"\nEvidence Quality:")
                prec_count = evidence.get('precedent_count', evidence.get('precedents', 0))
                diversity = evidence.get('diversity_score', evidence.get('diversity', 0))
                print(f"  Precedent count: {prec_count}")
                print(f"  Diversity score: {diversity:.3f}", end="")
                if diversity < 0.3:
                    print(" (LOW - potential batch effect)")
                elif diversity < 0.5:
                    print(" (MEDIUM)")
                else:
                    print(" (OK)")
            
            # Reasoning
            if 'reasoning' in meta:
                print(f"\nAdaptive Weight Reasoning:")
                for reason in meta['reasoning'][:5]:  # Show first 5
                    print(f"  - {reason}")
        else:
            print("\n❌ Fusion metadata missing!")
        
        # Check recommendations
        if 'formatted' in result and 'recommended_conditions' in result['formatted']:
            recs = result['formatted']['recommended_conditions']
            print(f"\n✅ Recommendations: {len(recs)} generated")
            
            if recs:
                print(f"\nTop Recommendation:")
                top = recs[0]
                summary = top.get('summary', {})
                print(f"  Rank: {summary.get('rank', top.get('rank', 1))}")
                print(f"  Core: {summary.get('core', 'N/A')}")
                
                base = summary.get('base', {})
                base_name = base.get('name', 'N/A') if isinstance(base, dict) else str(base)
                print(f"  Base: {base_name}")
                
                solvent = summary.get('solvent', {})
                solvent_name = solvent.get('name', 'N/A') if isinstance(solvent, dict) else str(solvent)
                print(f"  Solvent: {solvent_name}")
                print(f"  Confidence: {summary.get('confidence', 'N/A')}")
                
                # Component scores (fusion-specific)
                if 'component_scores' in top:
                    cs = top['component_scores']
                    print(f"\n  Component Scores:")
                    print(f"    PS (Precedent):  {cs.get('PS', 0):.3f}")
                    print(f"    AS (Analytics):  {cs.get('AS', 0):.3f}")
                    print(f"    RS (Rules):      {cs.get('RS', 0):.3f}")
                    print(f"    MS (ML):         {cs.get('MS', 0):.3f}")
        
        print("\n" + "-" * 80)
        return True
        
    except requests.exceptions.RequestException as e:
        print(f"❌ Request failed: {e}")
        return False
    except json.JSONDecodeError as e:
        print(f"❌ JSON decode error: {e}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_baseline_vs_fusion():
    """Compare baseline and fusion endpoints side-by-side."""
    print_section("BASELINE VS FUSION COMPARISON")
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    print(f"Reaction: {reaction}")
    print("Description: Buchwald-Hartwig C-N coupling")
    
    # Test baseline endpoint
    print("\n[1] BASELINE ENDPOINT (/api/v1/recommend)")
    print("-" * 80)
    baseline_payload = {"reaction": reaction, "k": 50}
    
    try:
        response = requests.post(f"{API_BASE}/api/v1/recommend", json=baseline_payload, timeout=30)
        response.raise_for_status()
        baseline = response.json()
        
        print("✅ Baseline response received")
        print(f"Keys: {list(baseline.keys())}")
        
        # Check for fusion_meta absence
        if 'fusion_meta' in baseline:
            print("⚠️  WARNING: fusion_meta found in baseline (should not be present)")
        else:
            print("✅ No fusion_meta in baseline (correct)")
        
        # Check model type
        if 'formatted' in baseline and 'meta' in baseline['formatted']:
            model = baseline['formatted']['meta'].get('model', 'unknown')
            print(f"Model: {model}")
        
    except Exception as e:
        print(f"❌ Baseline request failed: {e}")
    
    # Test fusion endpoint
    print("\n\n[2] FUSION ENDPOINT (/api/v1/recommend/fusion)")
    print("-" * 80)
    fusion_payload = {"reaction": reaction, "k": 50, "max_variants": 5}
    
    try:
        response = requests.post(FUSION_ENDPOINT, json=fusion_payload, timeout=30)
        response.raise_for_status()
        fusion = response.json()
        
        print("✅ Fusion response received")
        print(f"Keys: {list(fusion.keys())}")
        
        # Check for fusion_meta presence
        if 'fusion_meta' in fusion:
            print("✅ fusion_meta present (correct)")
            
            # Show weights
            weights = fusion['fusion_meta'].get('adaptive_weights', {})
            print(f"\nAdaptive Weights:")
            for key, val in weights.items():
                print(f"  {key}: {val:.3f}")
        else:
            print("❌ fusion_meta missing (should be present)")
        
        # Check model type
        if 'formatted' in fusion and 'meta' in fusion['formatted']:
            model = fusion['formatted']['meta'].get('model', 'unknown')
            print(f"Model: {model}")
        
    except Exception as e:
        print(f"❌ Fusion request failed: {e}")


def test_various_reactions():
    """Test fusion endpoint with various reaction types."""
    print_section("TESTING VARIOUS REACTION TYPES")
    
    test_cases = [
        {
            "reaction": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
            "description": "Suzuki - Simple Ph-Ph coupling",
            "k": 50,
            "max_variants": 5
        },
        {
            "reaction": "Clc1ccc(C#N)cc1.c1ccc(B(O)O)cc1>>N#Cc1ccc(-c2ccccc2)cc1",
            "description": "Suzuki - Electron-poor ArCl",
            "k": 50,
            "max_variants": 5
        },
        {
            "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            "description": "Buchwald-Hartwig - Diphenylamine",
            "k": 50,
            "max_variants": 5
        },
        {
            "reaction": "Brc1ccccc1.NC1CCCCC1>>c1ccccc1NC1CCCCC1",
            "description": "Buchwald-Hartwig - Cyclohexylamine",
            "k": 50,
            "max_variants": 3
        },
    ]
    
    results = []
    for i, test_case in enumerate(test_cases, 1):
        print(f"\n\nTest {i}/{len(test_cases)}")
        success = test_fusion_endpoint(**test_case)
        results.append((test_case['description'], success))
    
    # Summary
    print_section("TEST SUMMARY")
    passed = sum(1 for _, success in results if success)
    total = len(results)
    print(f"\nResults: {passed}/{total} tests passed")
    
    for desc, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"  {status}: {desc}")


def check_server_health():
    """Check if the API server is running."""
    print_section("CHECKING SERVER HEALTH")
    try:
        response = requests.get(f"{API_BASE}/health", timeout=5)
        response.raise_for_status()
        print("✅ Server is running")
        print(f"Response: {response.json()}")
        return True
    except Exception as e:
        print(f"❌ Server is not accessible: {e}")
        print(f"\nPlease start the server with:")
        print(f"  uvicorn app.main:app --reload --port 8000")
        return False


def main():
    """Run all tests."""
    print_section("FUSION API ENDPOINT TEST SUITE")
    print(f"API Base: {API_BASE}")
    print(f"Fusion Endpoint: {FUSION_ENDPOINT}")
    
    # Check server
    if not check_server_health():
        return
    
    # Run tests
    print("\n")
    input("Press Enter to start tests...")
    
    # Test 1: Baseline vs Fusion comparison
    test_baseline_vs_fusion()
    
    # Test 2: Various reactions
    print("\n\n")
    input("Press Enter to test various reactions...")
    test_various_reactions()
    
    print_section("ALL TESTS COMPLETE")


if __name__ == "__main__":
    main()
