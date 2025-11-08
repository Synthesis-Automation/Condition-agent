"""Test explicitly specifying RCM database."""

from chem_assistant.chemtools_wrapper import rule_based_conditions_tool

rcm_reaction = "C=CCNC(=O)C=C>>C1=CCNC(=O)C1"

print("Testing RCM with EXPLICIT database specification:")
print(f"Reaction: {rcm_reaction}\n")

# Specify the database explicitly
result = rule_based_conditions_tool.invoke({
    "reaction_smiles": rcm_reaction,
    "database": "RCM_db",  # Explicitly specify
    "auto_detect": False,  # Don't auto-detect
    "include_summary": True
})

print("=" * 80)
if result.get("success"):
    print("✅ SUCCESS - Explicit database specification works!\n")
    
    print(f"Database Used: {result.get('database_name', 'N/A')}\n")
    
    # Show summary
    if "summary" in result:
        print("Conditions Summary:")
        print("-" * 80)
        print(result["summary"])
    
    print("\n" + "=" * 80)
    print(f"⏱️  Time: {result.get('timing_ms', 0)} ms")
    
else:
    print("❌ FAILED\n")
    print(f"Error: {result.get('error', 'Unknown error')}")
