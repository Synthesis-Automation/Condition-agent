"""Test LLM workflow inline"""
import sys
from pathlib import Path

# Add paths
root = Path(__file__).parent.parent
sys.path.insert(0, str(root))
sys.path.insert(0, str(root / "data-processor"))

print("=" * 60)
print("Simple LLM Workflow Test")
print("=" * 60)

# Test basic imports
print("\n1. Importing LLMClient...")
try:
    from llmtools.clients import LLMClient
    print("   ✓ Success")
except Exception as e:
    print(f"   ❌ Failed: {e}")
    sys.exit(1)

print("\n2. Importing classifier functions...")
try:
    from llmtools.reagent_classifier import classify_role
    print("   ✓ Success")
except Exception as e:
    print(f"   ❌ Failed: {e}")
    sys.exit(1)

print("\n3. Importing reagent_taxonomy_generator...")
try:
    from reagent_taxonomy_generator import normalize_cas, resolve_identity_from_cas
    print("   ✓ Success")
except Exception as e:
    print(f"   ❌ Failed: {e}")
    sys.exit(1)

print("\n4. Testing CAS resolution (without LLM)...")
try:
    identity = resolve_identity_from_cas("121-44-8")
    if identity:
        print(f"   ✓ Resolved: {identity.get('name')}")
    else:
        print("   ❌ No data returned")
except Exception as e:
    print(f"   ❌ Failed: {e}")

print("\n5. Testing LLM client initialization...")
try:
    client = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
    print(f"   ✓ Client created: {client}")
except Exception as e:
    print(f"   ❌ Failed: {e}")
    sys.exit(1)

print("\n6. Testing role classification (LLM call)...")
try:
    result = classify_role(identity, client)
    if result.get("status") == "success":
        print(f"   ✓ Role: {result.get('role')} (confidence: {result.get('confidence'):.0%})")
    else:
        print(f"   ❌ Failed: {result.get('error')}")
except Exception as e:
    print(f"   ❌ Exception: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 60)
print("Test Complete")
print("=" * 60)
