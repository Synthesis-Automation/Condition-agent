#!/usr/bin/env python3
"""Quick test to verify the assistant can access the bond analysis tool."""

import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chem_assistant.chemtools_wrapper import CHEMTOOLS_TOOLS

print("=" * 80)
print("Checking Available Tools in ChemTools Assistant")
print("=" * 80)
print()

print(f"Total tools available: {len(CHEMTOOLS_TOOLS)}")
print()

print("Tool List:")
print("-" * 80)
for i, tool in enumerate(CHEMTOOLS_TOOLS, 1):
    tool_name = tool.name
    tool_desc = tool.description.split('\n')[0]  # First line
    print(f"{i:2d}. {tool_name:35s} - {tool_desc}")

print()
print("=" * 80)

# Check if bond analysis tool is present
bond_tool = next((t for t in CHEMTOOLS_TOOLS if 'bond' in t.name.lower()), None)
if bond_tool:
    print("✅ Bond Analysis Tool Found!")
    print(f"   Name: {bond_tool.name}")
    print(f"   Description: {bond_tool.description[:200]}...")
else:
    print("❌ Bond Analysis Tool NOT found!")

print("=" * 80)
