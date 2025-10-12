#!/usr/bin/env python3
"""
Test the markdown fence stripping fix for LLM responses.
"""

import json
import sys
from pathlib import Path

# Add parent to path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from llmtools.reagent_review import _strip_markdown_fences

# Your actual DeepSeek response
deepseek_response = """```json
{
  "status": "reject",
  "proposed_role": "ligand",
  "proposed_family": "phosphine",
  "confidence": 0.9,
  "justification": "The proposed family 'phosphine_pyridine_PN' is incorrect as the reagent does not contain a pyridine moiety. It should be classified under the more general 'phosphine' family.",
  "alerts": ["The current family assignment is incorrect due to the absence of a pyridine group in the reagent."],
  "suggested_synonyms": [],
  "field_suggestions": {
    "SMILES": "C=CCP(c1ccccc1)c2ccccc2",
    "molecular_formula": "C15H15P"
  }
}
```"""

print("Testing markdown fence stripping...\n")
print("Original response (first 100 chars):")
print(repr(deepseek_response[:100]))
print()

# Strip fences
cleaned = _strip_markdown_fences(deepseek_response)
print("Cleaned response (first 100 chars):")
print(repr(cleaned[:100]))
print()

# Try parsing
try:
    parsed = json.loads(cleaned)
    print("✅ SUCCESS! JSON parsed correctly")
    print()
    print("Parsed data:")
    print(f"  Status: {parsed['status']}")
    print(f"  Proposed role: {parsed['proposed_role']}")
    print(f"  Proposed family: {parsed['proposed_family']}")
    print(f"  Confidence: {parsed['confidence']}")
    print(f"  Field suggestions: {parsed['field_suggestions']}")
    print()
    print("Full parsed object:")
    print(json.dumps(parsed, indent=2))
except json.JSONDecodeError as e:
    print(f"❌ FAILED: {e}")

# Test other formats
print("\n" + "="*60)
print("Testing other markdown formats...\n")

test_cases = [
    ('```\n{"test": true}\n```', "Generic fence"),
    ('```JSON\n{"test": true}\n```', "Uppercase JSON"),
    ('{"test": true}', "No fence"),
    ('   ```json\n  {"test": true}  \n```   ', "Extra whitespace"),
]

for test_input, description in test_cases:
    cleaned = _strip_markdown_fences(test_input)
    try:
        parsed = json.loads(cleaned)
        print(f"✅ {description}: PASS")
    except Exception as e:
        print(f"❌ {description}: FAIL - {e}")
