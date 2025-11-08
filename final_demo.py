"""
Final demonstration: The exact reaction from the user's request.
"""

from chem_assistant.chemtools_wrapper import rule_based_conditions_tool

# User's Sonogashira reaction
reaction = "Ic1ccccc1C(C)(C)C.C#Cc1ccccc1>>CC(C)(C)c1ccccc1C#Cc1ccccc1"

print("=" * 80)
print("USER REQUEST: Use rule to find conditions for:")
print(f"  {reaction}")
print("=" * 80)
print()

result = rule_based_conditions_tool.invoke({
    "reaction_smiles": reaction,
    "auto_detect": True,
    "include_summary": True
})

if result.get("success"):
    print("✅ AGENT RESPONSE: Successfully found conditions!\n")
    print(result["summary"])
    print()
    print(f"Database: {result.get('database_name')}")
    print(f"Time: {result.get('timing_ms')} ms")
else:
    print(f"❌ Failed: {result.get('error')}")
