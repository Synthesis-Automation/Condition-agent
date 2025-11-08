"""Test agent detection and rule finding for Sonogashira reaction."""

from chem_assistant.chemtools_wrapper import rule_based_conditions_tool

# Test Sonogashira reaction
reaction = "Ic1ccccc1C(C)(C)C.C#Cc1ccccc1>>CC(C)(C)c1ccccc1C#Cc1ccccc1"

print("Testing Sonogashira reaction with agent tool:")
print(f"Reaction: {reaction}\n")

# Call the rule-based tool with auto-detection (it's a LangChain tool)
result = rule_based_conditions_tool.invoke({
    "reaction_smiles": reaction,
    "auto_detect": True,
    "include_summary": True
})

print("=" * 80)
if result.get("success"):
    print("✅ SUCCESS - Rule database found and conditions generated!\n")
    
    # Show detection info
    if "family_detection" in result:
        detection = result["family_detection"]
        print(f"Family Detected: {detection.get('family_name', 'N/A')}")
        print(f"Confidence: {detection.get('confidence', 0):.2f}\n")
    
    # Show database used
    print(f"Database Used: {result.get('database_name', 'N/A')}")
    print(f"Database Path: {result.get('database_path', 'N/A')}\n")
    
    # Show summary
    if "summary" in result:
        print("Conditions Summary:")
        print("-" * 80)
        print(result["summary"])
    
    print("\n" + "=" * 80)
    print(f"⏱️  Time: {result.get('timing_ms', 0)} ms")
    
else:
    print("❌ FAILED - Could not generate conditions\n")
    print(f"Error: {result.get('error', 'Unknown error')}")
    print(f"Attempted identifiers: {result.get('attempted_identifiers', [])}")
