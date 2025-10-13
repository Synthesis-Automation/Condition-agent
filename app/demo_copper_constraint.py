"""
Demo: How "use copper catalyst" affects the API request
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from app.cli_recommend import ParsedRequest

print("=" * 70)
print("DEMO: Positive Constraints - 'use copper catalyst'")
print("=" * 70)

# Example 1: User says "use copper catalyst"
print("\n📝 User Input:")
print('   Reaction: "Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"')
print('   Requirements: "use copper catalyst"')

# Simulated LLM parsing result
request = ParsedRequest(
    reaction_smiles="Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    reaction_smiles_is_valid=True,
    reaction_type="Ullmann_C_N",
    constraints={
        "metal_preference": "Cu",
        "required_reagents": ["copper catalyst"]
    }
)

print("\n🤖 LLM Parsed Constraints:")
print(f"   - metal_preference: {request.constraints.get('metal_preference')}")
print(f"   - required_reagents: {request.constraints.get('required_reagents')}")

# Convert to API request
api_request = request.to_api_request()

print("\n📤 API Request Sent:")
print(f"   reaction: {api_request['reaction'][:50]}...")
print(f"   reaction_type: {api_request['reaction_type']}")
print(f"   constraints:")
for key, value in api_request['constraints'].items():
    print(f"      - {key}: {value}")

print("\n✅ Result:")
print("   The API will now filter/rank results to prefer copper-based catalysts")
print("   Expected recommendations: CuI, Cu(OAc)2, Cu2O, etc.")

print("\n" + "=" * 70)
print("\nComparison WITHOUT copper preference:")
print("=" * 70)

# Example 2: No metal preference
request_no_pref = ParsedRequest(
    reaction_smiles="Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    reaction_smiles_is_valid=True,
    reaction_type="Ullmann_C_N",
    constraints={}
)

api_request_no_pref = request_no_pref.to_api_request()

print("\n📤 API Request Sent:")
print(f"   reaction: {api_request_no_pref['reaction'][:50]}...")
print(f"   reaction_type: {api_request_no_pref['reaction_type']}")
print(f"   constraints: {api_request_no_pref['constraints']}")

print("\n⚠️  Result:")
print("   The API returns all suitable catalysts without metal preference")
print("   Results may include: Pd, Cu, Ni, Fe catalysts")

print("\n" + "=" * 70)
print("\nOther Metal Preferences:")
print("=" * 70)

examples = [
    ("palladium catalyst", "Pd"),
    ("nickel catalyst", "Ni"),
    ("iron catalyst", "Fe"),
]

for phrase, metal in examples:
    print(f'\n   "{phrase}" → metal_preference: "{metal}"')
    req = ParsedRequest(
        reaction_smiles="A.B>>C",
        reaction_smiles_is_valid=True,
        constraints={"metal_preference": metal}
    )
    api_req = req.to_api_request()
    print(f"   API constraint: {api_req['constraints']}")

print("\n" + "=" * 70)
