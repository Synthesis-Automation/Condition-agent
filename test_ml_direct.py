"""
Direct test for ML endpoint to diagnose the 500 error
"""
import requests
import json

BASE_URL = "http://localhost:8000"
REACTION = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"

print("Testing ML endpoint directly...")
print(f"Reaction: {REACTION}\n")

try:
    payload = {
        "reaction": REACTION,
        "reaction_type": None,
        "k": 50,
        "limit": 3,
        "relax": {},
        "constraints": {}
    }
    
    print("Sending request...")
    response = requests.post(
        f"{BASE_URL}/api/v1/recommend/conditions",
        json=payload,
        timeout=30
    )
    
    print(f"Status: {response.status_code}")
    
    if response.status_code == 200:
        result = response.json()
        print("✓ Success!")
        print(json.dumps(result, indent=2)[:500])
    else:
        print(f"✗ Error: {response.status_code}")
        print(f"Response: {response.text}")
        
except Exception as e:
    print(f"✗ Exception: {e}")
    import traceback
    traceback.print_exc()
