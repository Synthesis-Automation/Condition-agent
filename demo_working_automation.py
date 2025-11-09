"""
Working Example: Protocol + Automation Format
==============================================

Since there's an index mismatch for rules (sonogashira_v2.json vs sonogashira_db.json),
let's demonstrate with PROTOCOLS which work correctly.
"""

import subprocess
import sys
import os

os.environ['PYTHONIOENCODING'] = 'utf-8'

commands = [
    "/automation on",       # Enable automation format
    "/scale 2.0",          # Set 2 mmol scale
    "/type protocol",      # Filter to protocols (these work!)
    "/k 2",                # Show top 2
    "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",  # Sonogashira reaction
    "/exit"
]

print("=" * 80)
print("✅ WORKING EXAMPLE: Protocols + Automation Format")
print("=" * 80)
print()
print("📝 Commands:")
for i, cmd in enumerate(commands[:-1], 1):
    print(f"  {i}. {cmd}")
print()
print("🔄 Running CLI...")
print("=" * 80)
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

if "🤖 Automation Format:" in result.stdout and "Addition Sequence:" in result.stdout:
    print()
    print("=" * 80)
    print("✅ SUCCESS: Automation format working perfectly!")
    print("=" * 80)
    print()
    print("💡 Note: Rules have an index mismatch issue (sonogashira_v2.json vs sonogashira_db.json)")
    print("   that needs to be fixed by rebuilding the index.")
