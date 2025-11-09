"""
Simple CLI Demo: Rule + Automation Format
==========================================

This shows how to get rule-based recommendations with automation format
showing the ordered addition sequence.
"""

import subprocess
import sys
import os

# Set UTF-8 encoding
os.environ['PYTHONIOENCODING'] = 'utf-8'

# Commands for interactive CLI
commands = [
    "/automation on",       # Enable automation format
    "/scale 2.0",          # Set 2 mmol scale
    "/type rule",          # Filter to rules only
    "/k 3",                # Show top 3
    "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",  # Sonogashira reaction
    "/exit"
]

print("=" * 80)
print("CLI Demo: Rules + Automation Format")
print("=" * 80)
print()
print("📝 Commands:")
for i, cmd in enumerate(commands[:-1], 1):
    print(f"  {i}. {cmd}")
print()
print("🔄 Running...")
print("=" * 80)
print()

# Run the CLI
input_text = "\n".join(commands)
result = subprocess.run(
    [sys.executable, "app/unified_rule_protocol_interactive_cli.py"],
    input=input_text,
    capture_output=True,
    text=True,
    encoding='utf-8',
    errors='replace'
)

# Print the output
print(result.stdout)

# Check for success
if "🤖 Automation Format:" in result.stdout:
    print()
    print("=" * 80)
    print("✅ SUCCESS: Automation format is showing!")
    print("=" * 80)
elif "Addition Sequence:" in result.stdout:
    print()
    print("=" * 80)
    print("✅ SUCCESS: Addition sequence displayed!")
    print("=" * 80)
else:
    print()
    print("=" * 80)
    print("ℹ️  Results shown (automation format may not be visible if no good matches)")
    print("=" * 80)
