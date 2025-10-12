"""
Complete LLM Workflow Test
Tests the full 4-step pipeline for reagent taxonomy generation
"""
import json
import sys
from pathlib import Path

# Add paths
root = Path(__file__).parent.parent
sys.path.insert(0, str(root))
sys.path.insert(0, str(root / "data-processor"))

print("\n🧪 LLM Workflow - Complete Test\n")
print("=" * 70)
print("Setup")
print("=" * 70)

# Imports
from llmtools.clients import LLMClient
from llmtools.reagent_classifier import classify_role, assign_fields, verify_entry
from reagent_taxonomy_generator import normalize_cas, resolve_identity_from_cas, dedupe_synonyms

# Initialize
client = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
registry_dir = root / "data" / "reagents"
cas = "121-44-8"  # Triethylamine

print(f"✓ LLM Client: DeepSeek V3.2 (experimental)")
print(f"✓ Registry: {registry_dir}")
print(f"✓ Test CAS: {cas}")

print("\n" + "=" * 70)
print("Step 1: Resolve Identity")
print("=" * 70)

try:
    identity = resolve_identity_from_cas(normalize_cas(cas))
    if not identity:
        print(f"❌ Failed to resolve CAS {cas}")
        sys.exit(1)
    
    print(f"✓ Name: {identity.get('name')}")
    print(f"✓ SMILES: {identity.get('smiles')}")
    print(f"✓ Formula: {identity.get('molecular_formula')}")
    print(f"✓ Synonyms: {len(identity.get('synonyms', []))} found")
    
except Exception as e:
    print(f"❌ Failed: {e}")
    sys.exit(1)

print("\n" + "=" * 70)
print("Step 2: Classify Role (LLM)")
print("=" * 70)

try:
    role_result = classify_role(identity, client)
    
    if role_result.get("status") != "success":
        print(f"❌ Failed: {role_result.get('error')}")
        sys.exit(1)
    
    role = role_result["role"]
    confidence = role_result["confidence"]
    reasoning = role_result["reasoning"]
    
    print(f"✓ Role: {role}")
    print(f"✓ Confidence: {confidence:.0%}")
    print(f"✓ Reasoning: {reasoning}")
    print(f"✓ Model: {role_result.get('model')} ({role_result.get('tokens')} tokens, {role_result.get('latency_ms')}ms)")
    
except Exception as e:
    print(f"❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 70)
print("Step 3: Assign Fields (LLM)")
print("=" * 70)

try:
    fields_result = assign_fields(
        identity=identity,
        role=role,
        registry_dir=registry_dir,
        llm_client=client,
    )
    
    if fields_result.get("status") != "success":
        print(f"❌ Failed: {fields_result.get('error')}")
        sys.exit(1)
    
    family = fields_result["family"]
    fields = fields_result["fields"]
    abbrevs = fields_result.get("abbreviations", [])
    
    print(f"✓ Family: {family}")
    print(f"✓ Fields:")
    for key, val in fields.items():
        print(f"    {key}: {val}")
    if abbrevs:
        print(f"✓ Abbreviations: {', '.join(abbrevs)}")
    print(f"✓ Confidence: {fields_result.get('confidence'):.0%}")
    print(f"✓ Model: {fields_result.get('model')} ({fields_result.get('tokens')} tokens, {fields_result.get('latency_ms')}ms)")
    
except Exception as e:
    print(f"❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 70)
print("Build Entry")
print("=" * 70)

try:
    entry = {
        "name": identity.get("name"),
        "cas": normalize_cas(cas),
        "molecular_formula": identity.get("molecular_formula"),
        "smiles": identity.get("smiles"),
        "roles": {
            role: {
                "families": [family],
                **fields,
            }
        },
    }
    
    if abbrevs:
        entry["abbreviations"] = abbrevs
    
    synonyms = identity.get("synonyms", [])
    if synonyms:
        entry["synonyms"] = synonyms[:10]  # Limit to first 10
    
    print(f"✓ Entry built successfully")
    print(f"\nEntry preview:")
    print(json.dumps(entry, indent=2, ensure_ascii=False)[:500] + "...")
    
except Exception as e:
    print(f"❌ Failed: {e}")
    sys.exit(1)

print("\n" + "=" * 70)
print("Step 4: Verify Entry (LLM)")
print("=" * 70)

try:
    verify_result = verify_entry(entry, client)
    
    if verify_result.get("status") != "success":
        print(f"❌ Failed: {verify_result.get('error')}")
        sys.exit(1)
    
    approved = verify_result["approved"]
    issues = verify_result.get("issues", [])
    suggestions = verify_result.get("suggestions", [])
    
    print(f"✓ Approved: {'YES ✓' if approved else 'NO ✗'}")
    print(f"✓ Issues: {len(issues)}")
    
    if issues:
        for i, issue in enumerate(issues, 1):
            severity = issue.get("severity", "unknown")
            field = issue.get("field", "")
            message = issue.get("message", "")
            print(f"    {i}. [{severity.upper()}] {field}: {message}")
    
    if suggestions:
        print(f"✓ Suggestions:")
        for i, sug in enumerate(suggestions, 1):
            print(f"    {i}. {sug}")
    
    print(f"✓ Model: {verify_result.get('model')} ({verify_result.get('tokens')} tokens, {verify_result.get('latency_ms')}ms)")
    
except Exception as e:
    print(f"❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 70)
print("Final Result")
print("=" * 70)

if approved and not any(iss.get("severity") == "error" for iss in issues):
    status = "ready_to_save"
    print("✅ SUCCESS - Entry ready to save!")
elif issues:
    status = "needs_review"
    error_count = sum(1 for iss in issues if iss.get("severity") == "error")
    warning_count = len(issues) - error_count
    print(f"⚠️  NEEDS REVIEW - {error_count} error(s), {warning_count} warning(s)")
else:
    status = "ready_to_save"
    print("✅ SUCCESS - Entry ready to save!")

print(f"\nFinal Status: {status}")
print(f"Entry Name: {entry['name']}")
print(f"Entry Role: {role}")
print(f"Entry Family: {family}")

print("\n" + "=" * 70)
print("✅ All workflow steps completed successfully!")
print("=" * 70)
