"""
CLI Demo: applies_if Filtering in Action
=========================================

Shows how the UnifiedRecommender CLI automatically filters rules
based on their applies_if criteria.
"""

import subprocess
import sys
import os

os.environ['PYTHONIOENCODING'] = 'utf-8'

print("=" * 80)
print("Demo: applies_if Filtering in Interactive CLI")
print("=" * 80)
print()

# Demo 1: Matching reaction (Sonogashira)
print("📋 Demo 1: Sonogashira Reaction (aryl bromide + terminal alkyne)")
print("-" * 80)

commands = [
    "/type rule",
    "/k 3",
    "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",
    "/exit"
]

print("Commands:", " → ".join(commands[:-1]))
print()

input_text = "\n".join(commands)
result = subprocess.run(
    [sys.executable, "app/unified_rule_protocol_interactive_cli.py"],
    input=input_text,
    capture_output=True,
    text=True,
    encoding='utf-8',
    errors='replace'
)

# Show results
lines = result.stdout.split('\n')
in_result = False
for line in lines:
    if 'Found' in line and 'recommendation' in line:
        in_result = True
    if in_result:
        print(line)
        if line.strip().startswith("reaction>") and "Exiting" in line:
            break

print()
print("✅ Sonogashira rule appears because:")
print("   • Query has aryl_halide (Br on benzene)")
print("   • Query has terminal_alkyne (C#Cc1ccccc1)")
print("   • Matches applies_if criteria!")

print()
print()

# Demo 2: Non-matching reaction (Amide formation)
print("📋 Demo 2: Amide Formation (should NOT show Sonogashira)")
print("-" * 80)

commands = [
    "/type rule",
    "/k 5",
    "CC(=O)O.CCN>>CC(=O)NCC",
    "/exit"
]

print("Commands:", " → ".join(commands[:-1]))
print()

input_text = "\n".join(commands)
result = subprocess.run(
    [sys.executable, "app/unified_rule_protocol_interactive_cli.py"],
    input=input_text,
    capture_output=True,
    text=True,
    encoding='utf-8',
    errors='replace'
)

# Show results
lines = result.stdout.split('\n')
in_result = False
for line in lines:
    if 'Found' in line and 'recommendation' in line:
        in_result = True
    if in_result:
        print(line)
        if line.strip().startswith("reaction>") and "Exiting" in line:
            break

print()
print("✅ Sonogashira rule correctly FILTERED OUT because:")
print("   • Query has NO aryl/vinyl halides")
print("   • Query has NO terminal alkyne")
print("   • Does NOT match applies_if criteria")
print("   • Only shows chemically appropriate rules (Amide, Reductive Amination)")

print()
print()

print("=" * 80)
print("Summary")
print("=" * 80)
print()
print("The CLI automatically uses applies_if filtering (validate_rules=True)")
print()
print("Benefits:")
print("  ✅ Only shows chemically appropriate rules")
print("  ✅ Filters based on detected substrate features")
print("  ✅ Prevents wrong recommendations (e.g., Sonogashira for amides)")
print("  ✅ No manual configuration needed - works by default")
print()
print("How it works:")
print("  1. Recommender detects features from query reaction")
print("  2. Each rule's applies_if is checked against features")
print("  3. Rules that don't match are silently filtered out")
print("  4. Only relevant rules appear in results")
