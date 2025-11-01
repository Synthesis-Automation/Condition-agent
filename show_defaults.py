#!/usr/bin/env python3
"""Display current LLM model defaults"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from llmtools.clients import RECOMMENDED_MODELS

print("="*60)
print("UPDATED LLM DEFAULT MODELS")
print("="*60)
print("\nOpenAI Recommended Models:")
for category, model in RECOMMENDED_MODELS['openai'].items():
    print(f"  {category:12} -> {model}")

print("\nAliyun Recommended Models:")
for category, model in RECOMMENDED_MODELS['aliyun'].items():
    print(f"  {category:12} -> {model}")
    
print("="*60)
print("✓ gpt-4o is now the default for both 'fast' and 'balanced'")
print("="*60)
