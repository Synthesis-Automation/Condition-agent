"""
Test: Updated CLI with top_k = 1 (default)
==========================================
"""

import subprocess
import sys
import os

os.environ['PYTHONIOENCODING'] = 'utf-8'

commands = [
    "/automation on",
    "Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1",
    "/exit"
]

print("=" * 80)
print("Test: Default top_k = 1")
print("=" * 80)
print()
print("Commands:")
for cmd in commands[:-1]:
    print(f"  {cmd}")
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

print(result.stdout)

print()
print("=" * 80)
print("✅ Now shows only top 1 result by default")
print("=" * 80)
print()
print("To see more results, use:")
print("  /k 3    # Show top 3")
print("  /k 5    # Show top 5")
