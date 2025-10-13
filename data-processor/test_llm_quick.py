#!/usr/bin/env python3
"""
Quick test for LLM workflow - minimal dependencies
This creates a simple Python script and runs it
"""

import json
import os
import subprocess
import sys
from pathlib import Path

# Check for API key
has_aliyun = os.getenv('ALIYUN_API_KEY')
has_openai = os.getenv('OPENAI_API_KEY')

if not has_aliyun and not has_openai:
    print("\n❌ No API key found!")
    print("\nPlease set environment variable:")
    print("  $env:ALIYUN_API_KEY = 'your-key'  (DeepSeek - recommended)")
    print("  OR")
    print("  $env:OPENAI_API_KEY = 'your-key'  (OpenAI)")
    print("\nExample:")
    print("  $env:ALIYUN_API_KEY = 'sk-...'; python test_llm_quick.py")
    sys.exit(1)

# Create test script
test_code = """
import sys
import json
from pathlib import Path

# Add paths
sys.path.insert(0, r'{root_dir}')
sys.path.insert(0, r'{data_processor_dir}')
sys.path.insert(0, r'{app_dir}')

# Import dependencies
from llmtools.clients import LLMClient
from llmtools.reagent_classifier import classify_role, assign_fields, verify_entry
from chemtools.reagent import normalize_cas, resolve_identity_from_cas, dedupe_synonyms

# Import from reagent_taxonomy_ui (now in app/)
exec(open(r'{qt_file}', encoding='utf-8').read())

# Setup
client = LLMClient(provider="{provider}", model="{model}")
registry_dir = Path(r'{registry_dir}')

# Run workflow
print("Running LLM workflow for Triethylamine (CAS 121-44-8)...")
print("=" * 60)

result = generate_taxonomy_entry_llm(
    cas="121-44-8",
    registry_dir=registry_dir,
    llm_client=client,
)

# Print result
print(json.dumps(result, indent=2, ensure_ascii=False))
"""

# Get paths
script_dir = Path(__file__).resolve().parent
root_dir = script_dir.parent
app_dir = root_dir / "app"
registry_dir = root_dir / "data" / "reagents"
qt_file = app_dir / "reagent_taxonomy_ui.py"

# Determine provider
if has_aliyun:
    provider = "aliyun"
    model = "deepseek-v3.2-exp"
    print(f"✓ Using DeepSeek V3.2 (Aliyun)")
else:
    provider = "openai"
    model = "gpt-4o-mini"
    print(f"✓ Using GPT-4o-mini (OpenAI)")

# Format code
formatted_code = test_code.format(
    root_dir=str(root_dir),
    data_processor_dir=str(script_dir),
    qt_file=str(qt_file),
    provider=provider,
    model=model,
    registry_dir=str(registry_dir),
)

# Write temp file
temp_script = script_dir / "_temp_test_llm.py"
temp_script.write_text(formatted_code, encoding='utf-8')

print(f"✓ Created temporary test script: {temp_script.name}")
print(f"✓ Registry directory: {registry_dir}")

print("\n" + "=" * 60)
print("Running Test")
print("=" * 60 + "\n")

try:
    # Run the script
    result = subprocess.run(
        [sys.executable, str(temp_script)],
        capture_output=True,
        text=True,
        timeout=60,
    )
    
    # Show output
    print(result.stdout)
    
    if result.stderr:
        print("STDERR:", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
    
    if result.returncode == 0:
        # Parse JSON output
        try:
            # Find JSON in output
            lines = result.stdout.split('\n')
            json_start = None
            for i, line in enumerate(lines):
                if line.strip().startswith('{'):
                    json_start = i
                    break
            
            if json_start is not None:
                json_text = '\n'.join(lines[json_start:])
                data = json.loads(json_text)
                
                print("\n" + "=" * 60)
                print("Test Summary")
                print("=" * 60)
                
                status = data.get("status")
                print(f"\nStatus: {status}")
                
                if status == "ready_to_save":
                    entry = data.get("entry", {})
                    workflow = data.get("workflow", {})
                    
                    print(f"\n✅ SUCCESS!")
                    print(f"   Name: {entry.get('name')}")
                    print(f"   Role: {workflow.get('step2_role', {}).get('role')}")
                    print(f"   Family: {workflow.get('step3_fields', {}).get('family')}")
                    print(f"   Approved: {workflow.get('step4_verification', {}).get('approved')}")
                    
                elif status == "needs_review":
                    print("\n⚠️  NEEDS REVIEW")
                    issues = data.get("workflow", {}).get("step4_verification", {}).get("issues", [])
                    print(f"   Issues: {len(issues)}")
                    
                else:  # error
                    print(f"\n❌ ERROR: {data.get('error')}")
                    
        except Exception as e:
            print(f"\nCould not parse JSON output: {e}")
            
    else:
        print("\n" + "=" * 60)
        print(f"❌ Test failed with exit code {result.returncode}")
        print("=" * 60)
        
finally:
    # Cleanup
    if temp_script.exists():
        temp_script.unlink()
        print(f"\n✓ Cleaned up temporary file")
