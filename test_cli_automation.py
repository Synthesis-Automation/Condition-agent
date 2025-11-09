"""
Quick test of automation format in interactive CLI.
"""

import subprocess
import sys

# Test commands
commands = [
    "/automation on",
    "/scale 2.0",
    "/k 3",
    "CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
    "/exit"
]

# Write commands to temp file
with open("_test_commands.txt", "w") as f:
    for cmd in commands:
        f.write(cmd + "\n")

print("=" * 80)
print("TESTING AUTOMATION FORMAT IN INTERACTIVE CLI")
print("=" * 80)
print()
print("Commands to execute:")
for i, cmd in enumerate(commands, 1):
    print(f"  {i}. {cmd}")
print()
print("=" * 80)
print()

# Run the CLI with input
result = subprocess.run(
    [sys.executable, "app/unified_rule_protocol_interactive_cli.py"],
    input="\n".join(commands),
    capture_output=True,
    text=True,
    encoding='utf-8'
)

print(result.stdout)

if result.stderr:
    print("STDERR:")
    print(result.stderr)

print()
print("=" * 80)
print("TEST COMPLETE")
print("=" * 80)
