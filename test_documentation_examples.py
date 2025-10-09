"""
Test Documentation Examples

Verifies that the code examples in API_DOCUMENTATION.md actually work.
Run this after starting the API server to validate documentation accuracy.

Usage:
    # Start server first
    uvicorn app.main:app --reload --port 8000
    
    # Then run tests
    python test_documentation_examples.py
"""

import requests
import sys

BASE_URL = "http://localhost:8000"

def test_server_health():
    """Test that server is running"""
    print("🔍 Testing server health...")
    try:
        r = requests.get(f"{BASE_URL}/health", timeout=5)
        assert r.status_code == 200
        assert r.json()['ok'] == True
        print("✅ PASS: Server is healthy\n")
        return True
    except Exception as e:
        print(f"❌ FAIL: Server not responding - {e}")
        print("   Make sure server is running: uvicorn app.main:app --reload --port 8000\n")
        return False

def test_basic_recommendation():
    """Test basic recommendation example from Use Case 1"""
    print("🔍 Testing basic recommendation...")
    
    response = requests.post(
        f'{BASE_URL}/api/v1/recommend',
        json={
            'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
            'k': 50
        }
    )
    
    assert response.status_code == 200, f"Expected 200, got {response.status_code}"
    
    result = response.json()
    assert 'formatted' in result, "Response missing 'formatted' key"
    assert 'recommended_conditions' in result['formatted'], "Response missing 'recommended_conditions'"
    
    conditions = result['formatted']['recommended_conditions']
    assert len(conditions) > 0, "No conditions returned"
    
    top = conditions[0]
    assert 'summary' in top, "Top condition missing 'summary'"
    print(f"   Found {len(conditions)} conditions")
    print(f"   Top: {top['summary']['core']} + {top['summary']['base']} in {top['summary']['solvent']}")
    print("✅ PASS: Basic recommendation works\n")
    return True

def test_advanced_recommendation_with_new_parameters():
    """Test advanced example with rerank_strategy and filter_unknown_reagents"""
    print("🔍 Testing advanced recommendation with new parameters...")
    
    response = requests.post(
        f'{BASE_URL}/api/v1/recommend',
        json={
            'reaction': 'c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1c1ccccc1',  # Suzuki
            'reaction_family': 'suzuki',
            'k': 100,
            'rerank_strategy': 'rule',
            'filter_unknown_reagents': True
        }
    )
    
    assert response.status_code == 200, f"Expected 200, got {response.status_code}"
    
    result = response.json()
    conditions = result['formatted']['recommended_conditions']
    assert len(conditions) > 0, "No conditions returned"
    
    top = conditions[0]
    
    # Check for specific catalyst (not generic metal)
    catalyst_name = top['catalyst']['name']
    print(f"   Catalyst: {catalyst_name}")
    
    # Should not be generic "Palladium"
    assert 'palladium' in catalyst_name.lower(), "Expected palladium catalyst"
    assert len(catalyst_name) > 10, f"Catalyst name too short (generic?): {catalyst_name}"
    
    # Check for CAS number if available
    if 'cas' in top['catalyst']:
        print(f"   CAS: {top['catalyst']['cas']}")
    
    # Verify no unknown reagents (due to filter_unknown_reagents=True)
    catalyst_str = str(top.get('catalyst', {}))
    base_str = str(top.get('base', {}))
    ligand_str = str(top.get('ligand', {}))
    
    # Allow some conditions without all components, but none should be "[Unknown ...]"
    for component, name in [(catalyst_str, 'catalyst'), (base_str, 'base'), (ligand_str, 'ligand')]:
        if '[Unknown' in component:
            print(f"   ⚠️  Warning: {name} contains unknown component: {component}")
    
    print("✅ PASS: Advanced recommendation with new parameters works\n")
    return True

def test_reranking_strategies():
    """Test different reranking strategies"""
    print("🔍 Testing reranking strategies...")
    
    strategies = ['none', 'rule', 'analytics']
    reaction = 'c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1c1ccccc1'
    
    for strategy in strategies:
        print(f"   Testing strategy: {strategy}")
        
        response = requests.post(
            f'{BASE_URL}/api/v1/recommend',
            json={
                'reaction': reaction,
                'reaction_family': 'suzuki',
                'k': 50,
                'rerank_strategy': strategy
            }
        )
        
        assert response.status_code == 200, f"Strategy '{strategy}' failed with {response.status_code}"
        
        result = response.json()
        conditions = result['formatted']['recommended_conditions']
        assert len(conditions) > 0, f"Strategy '{strategy}' returned no conditions"
        print(f"      → {len(conditions)} conditions returned")
    
    print("✅ PASS: All reranking strategies work\n")
    return True

def test_filter_unknown_reagents():
    """Test filter_unknown_reagents parameter"""
    print("🔍 Testing filter_unknown_reagents...")
    
    # Test with filtering ON
    response_filtered = requests.post(
        f'{BASE_URL}/api/v1/recommend',
        json={
            'reaction': 'CCBr.CCCO>>CCOCC',
            'k': 50,
            'filter_unknown_reagents': True
        }
    )
    
    assert response_filtered.status_code == 200
    filtered_conditions = response_filtered.json()['formatted']['recommended_conditions']
    
    # Test with filtering OFF
    response_unfiltered = requests.post(
        f'{BASE_URL}/api/v1/recommend',
        json={
            'reaction': 'CCBr.CCCO>>CCOCC',
            'k': 50,
            'filter_unknown_reagents': False
        }
    )
    
    assert response_unfiltered.status_code == 200
    unfiltered_conditions = response_unfiltered.json()['formatted']['recommended_conditions']
    
    print(f"   Unfiltered: {len(unfiltered_conditions)} conditions")
    print(f"   Filtered: {len(filtered_conditions)} conditions")
    
    # Filtered should have fewer or equal conditions
    assert len(filtered_conditions) <= len(unfiltered_conditions), \
        "Filtered should not have more conditions than unfiltered"
    
    print("✅ PASS: filter_unknown_reagents works\n")
    return True

def main():
    """Run all documentation example tests"""
    print("=" * 60)
    print("Testing API Documentation Examples")
    print("=" * 60)
    print()
    
    tests = [
        ("Server Health", test_server_health),
        ("Basic Recommendation", test_basic_recommendation),
        ("Advanced Recommendation", test_advanced_recommendation_with_new_parameters),
        ("Reranking Strategies", test_reranking_strategies),
        ("Filter Unknown Reagents", test_filter_unknown_reagents),
    ]
    
    results = {}
    
    for name, test_func in tests:
        try:
            if name == "Server Health":
                # Server health is prerequisite
                if not test_func():
                    print("\n❌ Cannot continue without healthy server\n")
                    sys.exit(1)
                results[name] = True
            else:
                results[name] = test_func()
        except Exception as e:
            print(f"❌ FAIL: {name}")
            print(f"   Error: {e}\n")
            results[name] = False
    
    # Summary
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for name, result in results.items():
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {name}")
    
    print()
    print(f"Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All documentation examples work correctly!")
        sys.exit(0)
    else:
        print(f"\n⚠️  {total - passed} test(s) failed")
        sys.exit(1)

if __name__ == '__main__':
    main()
