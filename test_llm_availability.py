#!/usr/bin/env python3
"""Test script to verify LLM availability for reagent_addition_ui.py"""

import sys
from pathlib import Path

# Setup paths like reagent_addition_ui.py does
MODULE_DIR = Path(__file__).resolve().parent / "data-processor"
ROOT_DIR = MODULE_DIR.parent

if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# Test imports exactly as reagent_addition_ui.py does
print("="*60)
print("LLM AVAILABILITY CHECK")
print("="*60)

LLM_SUPPORT_ERROR = None

try:
    from llmtools.clients import AVAILABLE_MODELS as LLM_AVAILABLE_MODELS
    from llmtools.clients import RECOMMENDED_MODELS as LLM_RECOMMENDED_MODELS
    print("✓ LLM clients imported successfully")
except Exception as exc:
    LLM_SUPPORT_ERROR = f"{exc.__class__.__name__}: {exc}"
    print(f"✗ LLM clients import failed: {LLM_SUPPORT_ERROR}")

try:
    from llmtools.reagent_review import LLMReviewOptions, review_taxonomy_proposal
    from llmtools.reagent_classifier import classify_role, assign_fields, verify_entry, VALID_ROLES
    from llmtools.clients import LLMClient
    LLM_SUPPORT_AVAILABLE = True
    print("✓ LLM reagent tools imported successfully")
except Exception as exc:
    LLM_SUPPORT_AVAILABLE = False
    if not LLM_SUPPORT_ERROR:
        LLM_SUPPORT_ERROR = f"{exc.__class__.__name__}: {exc}"
    print(f"✗ LLM reagent tools import failed: {exc.__class__.__name__}: {exc}")

print("="*60)
print(f"LLM_SUPPORT_AVAILABLE: {LLM_SUPPORT_AVAILABLE}")
print(f"LLM_SUPPORT_ERROR: {LLM_SUPPORT_ERROR}")
print("="*60)

if LLM_SUPPORT_AVAILABLE:
    print("✓✓✓ LLM is AVAILABLE and ready to use!")
    print("\nAvailable LLM Providers and Models:")
    for provider, models in LLM_AVAILABLE_MODELS.items():
        print(f"  • {provider}: {len(models)} models")
        print(f"    Recommended: {LLM_RECOMMENDED_MODELS.get(provider, {})}")
    print(f"\nValid reagent roles: {len(VALID_ROLES)} roles defined")
else:
    print("✗✗✗ LLM is NOT available")
    print(f"\nError: {LLM_SUPPORT_ERROR}")

print("="*60)
