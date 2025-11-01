#!/usr/bin/env python3
"""Test script to verify default provider and model selection"""

import sys
from pathlib import Path

# Setup paths
ROOT_DIR = Path(__file__).parent
sys.path.insert(0, str(ROOT_DIR))

from llmtools.clients import AVAILABLE_MODELS, RECOMMENDED_MODELS

print("="*60)
print("DEFAULT LLM PROVIDER & MODEL SELECTION")
print("="*60)

# Simulate what the UI does
providers = sorted(AVAILABLE_MODELS.keys())
print(f"\nAvailable providers (sorted): {providers}")
print(f"  Index 0: {providers[0]}")
print(f"  Index 1: {providers[1]}")

# Find openai index
openai_index = next((i for i, p in enumerate(providers) if p == "openai"), 0)
print(f"\n✓ OpenAI provider index: {openai_index}")

default_provider = providers[openai_index]
print(f"✓ Default provider: {default_provider}")

# Get recommended model for that provider
preferences = RECOMMENDED_MODELS.get(default_provider, {})
print(f"\nRecommended models for {default_provider}:")
for key in ("balanced", "fast", "reasoning", "advanced"):
    model = preferences.get(key)
    if model:
        print(f"  {key:12} -> {model}")

# Get the actual default model (balanced is first in priority)
default_model = None
for key in ("balanced", "fast", "reasoning", "advanced"):
    candidate = preferences.get(key)
    if candidate:
        default_model = candidate
        break

print(f"\n✓ Default model (first priority): {default_model}")

print("="*60)
print("SUMMARY")
print("="*60)
print(f"When UI starts:")
print(f"  Default Provider: {default_provider}")
print(f"  Default Model: {default_model}")
print("="*60)
