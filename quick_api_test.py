"""
Quick Test: Single Suzuki Reaction API Test
============================================

A minimal script to quickly test one Suzuki reaction with all three approaches.
Perfect for quick validation or debugging.

Usage:
    python quick_api_test.py
"""

import sys
import io

# Set UTF-8 encoding for Windows console
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

import requests
import json


BASE_URL = "http://localhost:8000"
SUZUKI_DB = "data/conditionDB/Suzuki_db.json"

# Simple benchmark reaction
REACTION = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"


def check_server():
    """Quick health check"""
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=2)
        return response.status_code == 200
    except:
        return False


def test_all():
    """Test all three approaches"""
    
    print("=" * 80)
    print("QUICK API TEST - Suzuki Ph-Br + Ph-B(OH)2")
    print("=" * 80)
    print(f"Reaction: {REACTION}\n")
    
    if not check_server():
        print("❌ Server not running at", BASE_URL)
        print("\nStart server with:")
        print("  uvicorn app.main:app --reload --port 8000")
        return
    
    print("✅ Server is running\n")
    
    # Test 1: Rule-Based
    print("-" * 80)
    print("1. RULE-BASED (/match)")
    print("-" * 80)
    try:
        response = requests.post(
            f"{BASE_URL}/match",
            json={"reaction": REACTION, "db": SUZUKI_DB, "include_trace": False},
            timeout=10
        )
        if response.status_code == 200:
            result = response.json()
            conditions = result.get("conditions", {})
            catalyst = conditions.get("catalyst", {})
            
            # Handle various data structures
            def extract_name(value):
                if isinstance(value, str):
                    return value
                if isinstance(value, dict):
                    return value.get("name", value.get("core", "N/A"))
                if isinstance(value, list) and len(value) > 0:
                    item = value[0]
                    if isinstance(item, str):
                        return item
                    if isinstance(item, dict):
                        return item.get("name", item.get("core", "N/A"))
                return str(value) if value else "N/A"
            
            print(f"✅ Match: {result.get('match_type', 'N/A')}")
            
            catalyst_name = extract_name(catalyst) if catalyst else "N/A"
            if catalyst_name != "N/A":
                print(f"   Catalyst: {catalyst_name}")
            
            base_name = extract_name(conditions.get("base"))
            if base_name != "N/A":
                print(f"   Base: {base_name}")
            
            solvent_name = extract_name(conditions.get("solvent"))
            if solvent_name != "N/A":
                print(f"   Solvent: {solvent_name}")
        elif response.status_code == 503:
            print(f"⚠️  Service unavailable (rule-based matcher not loaded)")
            print(f"   This may be expected if the module is not configured")
        elif response.status_code == 404:
            print(f"⚠️  Database not found: {SUZUKI_DB}")
            print(f"   Check that the file exists")
        else:
            print(f"❌ HTTP {response.status_code}")
            try:
                error_detail = response.json()
                print(f"   Detail: {error_detail.get('detail', 'N/A')}")
            except:
                print(f"   {response.text[:200]}")
    except Exception as e:
        print(f"❌ Error: {e}")
    
    # Test 2: ML-Based
    print("\n" + "-" * 80)
    print("2. ML-BASED (/api/v1/recommend/conditions)")
    print("-" * 80)
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/recommend/conditions",
            json={"reaction": REACTION, "k": 50, "limit": 3},
            timeout=30
        )
        if response.status_code == 200:
            result = response.json()
            recs = result.get("recommendations", [])
            
            print(f"✅ Found {len(recs)} recommendations")
            
            if recs:
                top = recs[0]
                reagents = top.get("reagents", [])
                
                catalyst = next((r for r in reagents if r.get("role") == "CATALYST"), {})
                base = next((r for r in reagents if r.get("role") == "BASE"), {})
                solvent = next((r for r in reagents if r.get("role") == "SOLVENT"), {})
                
                print(f"   Confidence: {top.get('confidence', 0):.2%}")
                print(f"   Catalyst: {catalyst.get('name', 'N/A')}")
                print(f"   Base: {base.get('name', 'N/A')}")
                print(f"   Solvent: {solvent.get('name', 'N/A')}")
                print(f"   Support: {top.get('precedent_count', 0)} precedents")
        else:
            print(f"❌ HTTP {response.status_code}")
            try:
                error_detail = response.json()
                print(f"   Detail: {error_detail.get('detail', 'N/A')}")
            except:
                print(f"   Response: {response.text[:300]}")
    except Exception as e:
        print(f"❌ Error: {e}")
    
    # Test 3: Fusion
    print("\n" + "-" * 80)
    print("3. FUSION (/api/v1/recommend/fusion)")
    print("-" * 80)
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/recommend/fusion",
            json={"reaction": REACTION, "k": 50, "max_variants": 3},
            timeout=30
        )
        if response.status_code == 200:
            result = response.json()
            fusion_meta = result.get("fusion_meta", {})
            
            if fusion_meta and "error" not in fusion_meta:
                print("✅ Fusion active")
                
                weights = fusion_meta.get("adaptive_weights", {})
                if weights:
                    print(f"   Weights: α={weights.get('alpha', 0):.3f}, "
                          f"β={weights.get('beta', 0):.3f}, "
                          f"γ={weights.get('gamma', 0):.3f}, "
                          f"δ={weights.get('delta', 0):.3f}")
                
                evidence = fusion_meta.get("evidence_summary", {})
                if evidence:
                    print(f"   Precedents: {evidence.get('precedent_count', 0)}")
                    print(f"   Diversity: {evidence.get('diversity_score', 0):.3f}")
                
                # Show top recommendation
                recs = result.get("formatted", {}).get("recommended_conditions", [])
                if not recs:
                    recs = result.get("recommended_conditions", [])
                
                if recs:
                    top = recs[0]
                    summary = top.get("summary", {})
                    
                    # Extract names from potentially complex structures
                    def format_reagent(value):
                        if isinstance(value, str):
                            return value
                        if isinstance(value, dict):
                            return value.get("name", value.get("abbreviation", "N/A"))
                        return "N/A"
                    
                    print(f"\n   Top Recommendation:")
                    print(f"   Core: {format_reagent(summary.get('core'))}")
                    print(f"   Base: {format_reagent(summary.get('base'))}")
                    print(f"   Solvent: {format_reagent(summary.get('solvent'))}")
                    print(f"   Confidence: {summary.get('confidence', 'N/A')}")
            else:
                print("⚠️  Fusion metadata not available")
        else:
            print(f"❌ HTTP {response.status_code}")
    except Exception as e:
        print(f"❌ Error: {e}")
    
    print("\n" + "=" * 80)
    print("✅ Quick test complete!")
    print("=" * 80)
    print("\nFor comprehensive testing, use: python test_suzuki_api.py")


if __name__ == "__main__":
    test_all()
