"""
Better Demo: Sonogashira reaction with rule filtering
"""

import subprocess
import sys

commands = [
    "/automation on",
    "/scale 2.0",
    "/type rule",
    "/k 3",
    "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",  # Sonogashira
    "/exit"
]

print("=" * 80)
print("DEMO: Sonogashira Rule with Automation Format (2 mmol scale)")
print("=" * 80)
print()
print("Reaction: Bromobenzene + Phenylacetylene → Diphenylacetylene")
print("Type: Rule-based only")
print()
print("=" * 80)
print()

result = subprocess.run(
    [sys.executable, "app/unified_rule_protocol_interactive_cli.py"],
    input="\n".join(commands),
    capture_output=True,
    text=True,
    encoding='utf-8'
)

# Extract key parts
lines = result.stdout.split('\n')
printing = False
for line in lines:
    if 'Query:' in line or printing:
        print(line)
        printing = True
    if 'Exiting' in line:
        break

if "Addition Sequence" in result.stdout:
    print()
    print("=" * 80)
    print("✅ SUCCESS: Automation format displayed for rules!")
    print("=" * 80)
else:
    print()
    print("⚠️  Note: Rule may not have matched well for this reaction")
