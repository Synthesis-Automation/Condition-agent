"""
Test if the agent can now use HTE tools correctly after system prompt update
"""

import sys
from pathlib import Path

# Add project root
project_root = Path(__file__).parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Clear any cached modules
for mod in list(sys.modules.keys()):
    if 'chemtools' in mod or 'chem_assistant' in mod:
        del sys.modules[mod]

from chem_assistant.chemtools_agent import ChemToolsAgent

print("=" * 70)
print("Testing Agent with HTE Tools")
print("=" * 70)

# Create fresh agent
agent = ChemToolsAgent(verbose=True)

# Test query
query = """use HTE system, recommend conditions for reactions: 
Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1, 
its a C-N coupling reaction, and the catalyst is copper"""

print(f"\nQuery: {query}\n")
print("=" * 70)
print("Agent Response:")
print("=" * 70)

try:
    response = agent.run(query)
    print(response)
    
    # Check if response mentions data found
    if "112" in response or "found" in response.lower():
        print("\n" + "=" * 70)
        print("✅ SUCCESS! Agent found the HTE data")
        print("=" * 70)
    elif "no data" in response.lower() or "not available" in response.lower():
        print("\n" + "=" * 70)
        print("❌ FAILED! Agent still says no data")
        print("=" * 70)
        print("\nThis might mean:")
        print("1. LLM is using cached response")
        print("2. Tool invocation has issues")
        print("3. Need to restart GUI app completely")
    
except Exception as e:
    print(f"\n❌ ERROR: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("NEXT STEPS:")
print("=" * 70)
print("1. If test shows ✅ - restart GUI app to pick up system prompt changes")
print("2. If test shows ❌ - there may be an LLM caching issue")
print("=" * 70)
