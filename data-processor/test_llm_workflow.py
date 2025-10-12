#!/usr/bin/env python3
"""
Test script for new LLM workflow
=================================

Tests the pure LLM workflow with a simple example.

Usage:
    python test_llm_workflow.py

Requirements:
    - LLM API key configured (OPENAI_API_KEY or ALIYUN_API_KEY)
    - Registry directory exists
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Add paths for imports
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent

# Add root to path for llmtools
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# Add data-processor to path for reagent_taxonomy_qt
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

try:
    from llmtools.clients import LLMClient
    from reagent_taxonomy_qt import generate_taxonomy_entry_llm
except ImportError as exc:
    print(f"❌ Import failed: {exc}")
    print(f"\nPaths checked:")
    print(f"  ROOT_DIR: {ROOT_DIR}")
    print(f"  SCRIPT_DIR: {SCRIPT_DIR}")
    print(f"\nsys.path:")
    for p in sys.path[:5]:
        print(f"  {p}")
    print("\nPlease ensure:")
    print("  1. llmtools is available")
    print("  2. LLM dependencies are installed")
    sys.exit(1)


def test_triethylamine():
    """Test with triethylamine (common base)."""
    print("=" * 60)
    print("Test: Triethylamine (CAS 121-44-8)")
    print("=" * 60)
    
    # Configure client - try DeepSeek first
    try:
        client = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
        print(f"✓ Using DeepSeek V3.2 (Aliyun)")
    except Exception:
        try:
            client = LLMClient(provider="openai", model="gpt-4o-mini")
            print(f"✓ Using GPT-4o-mini (OpenAI)")
        except Exception as exc:
            print(f"❌ Failed to initialize LLM client: {exc}")
            print("\nPlease set environment variable:")
            print("  ALIYUN_API_KEY=<your-key>  OR  OPENAI_API_KEY=<your-key>")
            return False
    
    # Registry path
    registry_dir = ROOT_DIR / "data" / "reagents"
    if not registry_dir.exists():
        print(f"❌ Registry directory not found: {registry_dir}")
        return False
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
        return False
    
    # Check result
    status = result.get("status")
    print(f"\nStatus: {status}")
    
    if status == "error":
        print(f"❌ Error: {result.get('error')}")
        return False
    
    # Display workflow steps
    workflow = result.get("workflow", {})
    
    print("\n" + "=" * 60)
    print("Workflow Results")
    print("=" * 60)
    
    # Step 1: Identity
    step1 = workflow.get("step1_identity", {})
    if step1.get("status") == "success":
        identity = step1.get("identity", {})
        print(f"\n✓ Step 1: Identity Resolved")
        print(f"  Name: {identity.get('name')}")
        print(f"  CAS: {identity.get('cas')}")
        print(f"  SMILES: {identity.get('smiles')}")
        print(f"  Formula: {identity.get('molecular_formula')}")
    
    # Step 2: Role
    step2 = workflow.get("step2_role", {})
    if step2.get("status") == "success":
        print(f"\n✓ Step 2: Role Classification")
        print(f"  Role: {step2.get('role')}")
        print(f"  Confidence: {step2.get('confidence', 0):.2f}")
        print(f"  Reasoning: {step2.get('reasoning')}")
        print(f"  Model: {step2.get('model')} ({step2.get('tokens')} tokens, {step2.get('latency_ms')}ms)")
    
    # Step 3: Fields
    step3 = workflow.get("step3_fields", {})
    if step3.get("status") == "success":
        print(f"\n✓ Step 3: Field Assignment")
        print(f"  Family: {step3.get('family')}")
        print(f"  Fields: {json.dumps(step3.get('fields', {}), indent=4)}")
        print(f"  Abbreviations: {step3.get('abbreviations')}")
        print(f"  Confidence: {step3.get('confidence', 0):.2f}")
        print(f"  Model: {step3.get('model')} ({step3.get('tokens')} tokens, {step3.get('latency_ms')}ms)")
    
    # Step 4: Verification
    step4 = workflow.get("step4_verification", {})
    if step4.get("status") == "success":
        approved = step4.get("approved")
        issues = step4.get("issues", [])
        suggestions = step4.get("suggestions", [])
        
        print(f"\n✓ Step 4: Verification")
        print(f"  Approved: {'✓ YES' if approved else '✗ NO'}")
        print(f"  Issues: {len(issues)}")
        
        if issues:
            for i, issue in enumerate(issues, 1):
                severity = issue.get("severity", "unknown")
                field = issue.get("field", "unknown")
                message = issue.get("message", "")
                print(f"    {i}. [{severity.upper()}] {field}: {message}")
        
        if suggestions:
            print(f"  Suggestions:")
            for i, suggestion in enumerate(suggestions, 1):
                print(f"    {i}. {suggestion}")
        
        print(f"  Model: {step4.get('model')} ({step4.get('tokens')} tokens, {step4.get('latency_ms')}ms)")
    
    # Final entry
    entry = result.get("entry", {})
    print("\n" + "=" * 60)
    print("Final Entry")
    print("=" * 60)
    print(json.dumps(entry, indent=2))
    
    # Summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    
    if status == "ready_to_save":
        print("✅ SUCCESS - Entry ready to save!")
        return True
    elif status == "needs_review":
        print("⚠️  NEEDS REVIEW - Entry has warnings/errors")
        print(f"   Message: {result.get('message')}")
        return True
    else:
        print(f"❌ FAILED - Status: {status}")
        return False


if __name__ == "__main__":
    print("\n🧪 Testing LLM Workflow Implementation\n")
    
    success = test_triethylamine()
    
    print("\n" + "=" * 60)
    if success:
        print("✅ All tests passed!")
    else:
        print("❌ Tests failed")
    print("=" * 60)
    
    sys.exit(0 if success else 1)
