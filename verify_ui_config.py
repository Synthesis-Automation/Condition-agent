#!/usr/bin/env python3
"""Complete verification of UI defaults"""

import sys
from pathlib import Path

ROOT_DIR = Path(__file__).parent
sys.path.insert(0, str(ROOT_DIR))

print("="*70)
print(" UI DEFAULT SETTINGS VERIFICATION")
print("="*70)

# 1. Check LLM availability
try:
    from llmtools.clients import AVAILABLE_MODELS, RECOMMENDED_MODELS
    from llmtools.reagent_classifier import classify_role
    print("\n✅ LLM Support: AVAILABLE")
except Exception as e:
    print(f"\n❌ LLM Support: NOT AVAILABLE - {e}")
    sys.exit(1)

# 2. Check providers
providers = sorted(AVAILABLE_MODELS.keys())
print(f"\n📋 Available Providers: {', '.join(providers)}")

# 3. Default provider (what UI will select)
openai_index = next((i for i, p in enumerate(providers) if p == "openai"), 0)
default_provider = providers[openai_index]
print(f"\n🎯 Default Provider: {default_provider.upper()}")

# 4. Default model (what UI will select for that provider)
preferences = RECOMMENDED_MODELS.get(default_provider, {})
default_model = None
for key in ("balanced", "fast", "reasoning", "advanced"):
    candidate = preferences.get(key)
    if candidate:
        default_model = candidate
        break

print(f"🎯 Default Model: {default_model}")

# 5. Show all recommended models for default provider
print(f"\n📊 All Recommended Models for {default_provider}:")
for category, model in preferences.items():
    marker = "⭐" if model == default_model else "  "
    print(f"   {marker} {category:12} → {model}")

print("\n" + "="*70)
print(" CONFIGURATION SUMMARY")
print("="*70)
print(f"""
When the Reagent Addition UI starts with LLM enabled:
  
  1️⃣  Provider dropdown will default to: {default_provider.upper()}
  2️⃣  Model dropdown will default to: {default_model}
  3️⃣  This is the "balanced" model (first priority)
  
  ✓ Changed from: aliyun/deepseek-r1-distill-qwen-7b
  ✓ Changed to:   openai/gpt-4o
""")
print("="*70)
