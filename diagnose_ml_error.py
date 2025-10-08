"""
Diagnostic script for ML endpoint 500 error
Tests the ML recommendation endpoint and shows detailed error info
"""
import sys
import io

# UTF-8 encoding for Windows
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

import requests
import json

BASE_URL = "http://localhost:8000"

# Test reaction that's failing
REACTION = "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1"

print("=" * 80)
print("ML ENDPOINT DIAGNOSTIC")
print("=" * 80)
print(f"\nTesting reaction: {REACTION}")
print(f"Endpoint: {BASE_URL}/api/v1/recommend/conditions\n")

try:
    payload = {
        "reaction": REACTION,
        "reaction_type": None,
        "k": 50,
        "limit": 5,
        "relax": {},
        "constraints": {}
    }
    
    print("Sending request...")
    response = requests.post(
        f"{BASE_URL}/api/v1/recommend/conditions",
        json=payload,
        timeout=30
    )
    
    print(f"Status Code: {response.status_code}\n")
    
    if response.status_code == 200:
        print("✅ SUCCESS!")
        result = response.json()
        print("\nResponse structure:")
        print(json.dumps(result, indent=2)[:500])
    else:
        print(f"❌ ERROR: {response.status_code}")
        print(f"\nResponse text: {response.text}")
        
        # Try to parse as JSON for more details
        try:
            error_json = response.json()
            print("\nError details (JSON):")
            print(json.dumps(error_json, indent=2))
        except:
            pass
    
    print("\n" + "=" * 80)
    print("IMPORTANT: Check the server terminal (uvicorn) for the Python traceback!")
    print("Look for lines containing:")
    print("  - ERROR: Exception in ASGI application")
    print("  - Traceback (most recent call last)")
    print("  - The actual exception (FileNotFoundError, KeyError, etc.)")
    print("=" * 80)
    
except requests.exceptions.Timeout:
    print("❌ Request timed out (>30 seconds)")
except requests.exceptions.ConnectionError:
    print("❌ Cannot connect to server")
    print("\nIs the server running?")
    print("Start with: uvicorn app.main:app --reload --port 8000")
except Exception as e:
    print(f"❌ Unexpected error: {e}")
    import traceback
    traceback.print_exc()
