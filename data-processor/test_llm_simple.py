#!/usr/bin/env python3
"""
Simple test for LLM workflow - uses direct execution instead of imports
"""

import json
import sys
from pathlib import Path

# Get paths
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent

# Add to path
sys.path.insert(0, str(ROOT_DIR))
sys.path.insert(0, str(SCRIPT_DIR))

print("\n🧪 Testing LLM Workflow Implementation\n")
print("=" * 60)
print("Import Test")
print("=" * 60)

# Test imports one by one
print(f"ROOT_DIR: {ROOT_DIR}")
print(f"SCRIPT_DIR: {SCRIPT_DIR}")

try:
    print("\n1. Testing llmtools.clients import...")
    from llmtools.clients import LLMClient
    print("   ✓ LLMClient imported successfully")
except Exception as exc:
    print(f"   ❌ Failed: {exc}")
    sys.exit(1)

try:
    print("\n2. Testing llmtools.reagent_classifier import...")
    from llmtools.reagent_classifier import classify_role, assign_fields, verify_entry
    print("   ✓ Classifier functions imported successfully")
except Exception as exc:
    print(f"   ❌ Failed: {exc}")
    sys.exit(1)

# Instead of importing, load the module manually
print("\n3. Loading generate_taxonomy_entry_llm function...")
try:
    # Read and execute the function
    exec_globals = {
        '__name__': '__main__',
        'Path': Path,
        'Optional': type(None),
        'Dict': dict,
        'Any': object,
    }
    
    # Import dependencies from the new chemtools.reagent package
    from chemtools.reagent import normalize_cas, resolve_identity_from_cas, dedupe_synonyms
    exec_globals['normalize_cas'] = normalize_cas
    exec_globals['resolve_identity_from_cas'] = resolve_identity_from_cas
    exec_globals['dedupe_synonyms'] = dedupe_synonyms
    exec_globals['classify_role'] = classify_role
    exec_globals['assign_fields'] = assign_fields
    exec_globals['verify_entry'] = verify_entry
    
    # Load the function from file (moved to app/)
    qt_file = SCRIPT_DIR.parent / "app" / "reagent_taxonomy_ui.py"
    with open(qt_file, 'r', encoding='utf-8') as f:
        code = f.read()
    
    # Find and extract the function
    import re
    pattern = r'def generate_taxonomy_entry_llm\(.*?\n(?:.*?\n)*?.*?return.*?\n\n'
    match = re.search(pattern, code, re.MULTILINE | re.DOTALL)
    
    if not match:
        raise RuntimeError("Could not find generate_taxonomy_entry_llm function")
    
    func_code = match.group(0)
    exec(func_code, exec_globals)
    generate_taxonomy_entry_llm = exec_globals['generate_taxonomy_entry_llm']
    
    print("   ✓ Function loaded successfully")
    
except Exception as exc:
    print(f"   ❌ Failed: {exc}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Now run the test
print("\n" + "=" * 60)
print("Test: Triethylamine (CAS 121-44-8)")
print("=" * 60)

# Check for API key
import os
has_aliyun = os.getenv('ALIYUN_API_KEY')
has_openai = os.getenv('OPENAI_API_KEY')

if not has_aliyun and not has_openai:
    print("\n❌ No API key found!")
    print("\nPlease set environment variable:")
    print("  $env:ALIYUN_API_KEY = 'your-key'")
    print("  OR")
    print("  $env:OPENAI_API_KEY = 'your-key'")
    sys.exit(1)

# Configure client
try:
    if has_aliyun:
        client = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
        print(f"✓ Using DeepSeek V3.2 (Aliyun)")
    else:
        client = LLMClient(provider="openai", model="gpt-4o-mini")
        print(f"✓ Using GPT-4o-mini (OpenAI)")
except Exception as exc:
    print(f"❌ Failed to initialize LLM client: {exc}")
    sys.exit(1)

# Registry path
registry_dir = ROOT_DIR / "data" / "reagents"
if not registry_dir.exists():
    print(f"❌ Registry directory not found: {registry_dir}")
    sys.exit(1)
print(f"✓ Registry directory: {registry_dir}")

# Run workflow
print("\n" + "-" * 60)
print("Running LLM workflow...")
print("-" * 60)

try:
    result = generate_taxonomy_entry_llm(
        cas="121-44-8",
        registry_dir=registry_dir,
        llm_client=client,
    )
except Exception as exc:
    print(f"❌ Workflow failed: {exc}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Display result
status = result.get("status")
print(f"\nStatus: {status}")

if status == "error":
    print(f"❌ Error: {result.get('error')}")
    sys.exit(1)

# Display workflow steps
workflow = result.get("workflow", {})

print("\n" + "=" * 60)
print("Workflow Results")
print("=" * 60)

# Step 1
step1 = workflow.get("step1_identity", {})
if step1.get("status") == "success":
    identity = step1.get("identity", {})
    print(f"\n✓ Step 1: Identity Resolved")
    print(f"  Name: {identity.get('name')}")
    print(f"  CAS: {identity.get('cas')}")
    print(f"  SMILES: {identity.get('smiles')}")

# Step 2
step2 = workflow.get("step2_role", {})
if step2.get("status") == "success":
    print(f"\n✓ Step 2: Role Classification")
    print(f"  Role: {step2.get('role')}")
    print(f"  Confidence: {step2.get('confidence', 0):.0%}")
    print(f"  Reasoning: {step2.get('reasoning')}")

# Step 3
step3 = workflow.get("step3_fields", {})
if step3.get("status") == "success":
    print(f"\n✓ Step 3: Field Assignment")
    print(f"  Family: {step3.get('family')}")
    print(f"  Fields: {json.dumps(step3.get('fields', {}), indent=4)}")
    print(f"  Confidence: {step3.get('confidence', 0):.0%}")

# Step 4
step4 = workflow.get("step4_verification", {})
if step4.get("status") == "success":
    approved = step4.get("approved")
    issues = step4.get("issues", [])
    print(f"\n✓ Step 4: Verification")
    print(f"  Approved: {'✓ YES' if approved else '✗ NO'}")
    print(f"  Issues: {len(issues)}")
    
    if issues:
        for i, issue in enumerate(issues, 1):
            severity = issue.get("severity", "unknown")
            field = issue.get("field", "unknown")
            message = issue.get("message", "")
            print(f"    {i}. [{severity.upper()}] {field}: {message}")

# Final entry
entry = result.get("entry", {})
print("\n" + "=" * 60)
print("Final Entry (Summary)")
print("=" * 60)
print(f"Name: {entry.get('name')}")
print(f"CAS: {entry.get('cas')}")
print(f"Roles: {list(entry.get('roles', {}).keys())}")

# Summary
print("\n" + "=" * 60)
if status == "ready_to_save":
    print("✅ SUCCESS - Entry ready to save!")
    sys.exit(0)
elif status == "needs_review":
    print("⚠️  NEEDS REVIEW - Entry has warnings/errors")
    sys.exit(0)
else:
    print(f"❌ FAILED - Status: {status}")
    sys.exit(1)
