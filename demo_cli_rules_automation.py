"""
Complete Interactive Example: Rule-based Recommendations with Automation Format
"""

import subprocess
import sys

commands = [
    "/show",                # Show current settings
    "/automation on",       # Enable automation format
    "/scale 5.0",          # Set 5 mmol scale
    "/type rule",          # Only rules (not protocols)
    "/k 2",                # Top 2 results
    "CCBr.c1ccccc1B(O)O>>CCc1ccccc1",  # Suzuki reaction
    "/exit"
]

print("=" * 80)
print("DEMO: Rule-Based Recommendations with Automation Format (5 mmol scale)")
print("=" * 80)
print()
print("Commands:")
for cmd in commands[:-1]:  # Don't show /exit
    print(f"  {cmd}")
print()
print("=" * 80)
print()

# Run CLI
result = subprocess.run(
    [sys.executable, "app/unified_rule_protocol_interactive_cli.py"],
    input="\n".join(commands),
    capture_output=True,
    text=True,
    encoding='utf-8'
)

print(result.stdout)

if "Addition Sequence" in result.stdout:
    print()
    print("=" * 80)
    print("✅ SUCCESS: Automation format with addition sequences displayed!")
    print("=" * 80)
