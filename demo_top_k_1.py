"""
Demo: top_k = 1 with type filtering
====================================
Shows how the new default works with /type command
"""

import subprocess
import sys
import os

os.environ['PYTHONIOENCODING'] = 'utf-8'

print("=" * 80)
print("Demo: Default top_k = 1 with Type Filtering")
print("=" * 80)
print()

# Test 1: Top 1 rule
print("Test 1: Top 1 Rule Only")
print("-" * 80)

commands = [
    "/type rule",
    "/automation on",
    "Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1",
    "/exit"
]

input_text = "\n".join(commands)
result = subprocess.run(
    [sys.executable, "app/unified_rule_protocol_interactive_cli.py"],
    input=input_text,
    capture_output=True,
    text=True,
    encoding='utf-8',
    errors='replace'
)

# Show result
lines = result.stdout.split('\n')
for line in lines:
    if 'Found' in line or line.strip().startswith('📖') or line.strip().startswith('Family:') or line.strip().startswith('Similarity:'):
        print(line)

print()
print()

# Test 2: Top 1 protocol
print("Test 2: Top 1 Protocol Only")
print("-" * 80)

commands = [
    "/type protocol",
    "/automation on",
    "Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1",
    "/exit"
]

input_text = "\n".join(commands)
result = subprocess.run(
    [sys.executable, "app/unified_rule_protocol_interactive_cli.py"],
    input=input_text,
    capture_output=True,
    text=True,
    encoding='utf-8',
    errors='replace'
)

# Show result
lines = result.stdout.split('\n')
for line in lines:
    if 'Found' in line or line.strip().startswith('📋') or line.strip().startswith('Family:') or line.strip().startswith('Similarity:'):
        print(line)

print()
print()

# Test 3: Split mode (top rule + top protocol)
print("Test 3: Split Mode (Top Rule + Top Protocol)")
print("-" * 80)

commands = [
    "/split on",
    "/automation on",
    "Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1",
    "/exit"
]

input_text = "\n".join(commands)
result = subprocess.run(
    [sys.executable, "app/unified_rule_protocol_interactive_cli.py"],
    input=input_text,
    capture_output=True,
    text=True,
    encoding='utf-8',
    errors='replace'
)

# Show result
lines = result.stdout.split('\n')
in_result = False
for line in lines:
    if 'TOP RULE' in line or 'TOP PROTOCOL' in line:
        in_result = True
    if in_result:
        print(line)
        if 'Exiting' in line:
            break

print()
print("=" * 80)
print("Summary")
print("=" * 80)
print()
print("Default behavior (top_k = 1):")
print("  • /type rule      → Shows top 1 rule")
print("  • /type protocol  → Shows top 1 protocol")
print("  • /type all       → Shows top 1 overall (rule or protocol)")
print("  • /split on       → Shows top 1 rule AND top 1 protocol")
print()
print("To see more, use /k command:")
print("  /k 3  → Show top 3 results")
print("  /k 5  → Show top 5 results")
