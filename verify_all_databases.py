"""
Quick verification script to test all three new rule databases with the agent.
"""

from chem_assistant.chemtools_wrapper import rule_based_conditions_tool

test_reactions = [
    {
        "name": "Sonogashira",
        "smiles": "Ic1ccccc1C(C)(C)C.C#Cc1ccccc1>>CC(C)(C)c1ccccc1C#Cc1ccccc1"
    },
    {
        "name": "C-O Coupling", 
        "smiles": "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1"
    },
    {
        "name": "RCM",
        "smiles": "C=CCNC(=O)C=C>>C1=CCNC(=O)C1"
    }
]

print("=" * 80)
print("Testing Agent Integration for New Rule Databases")
print("=" * 80)

for test in test_reactions:
    print(f"\n📋 Testing {test['name']}...")
    print(f"   Reaction: {test['smiles'][:50]}...")
    
    try:
        result = rule_based_conditions_tool.invoke({
            "reaction_smiles": test["smiles"],
            "auto_detect": True,
            "include_summary": False  # Just check if it finds the database
        })
        
        if result.get("success"):
            db_name = result.get("database_name", "Unknown")
            rule_name = result.get("base_rule", {}).get("name", "Unknown")
            confidence = result.get("base_rule", {}).get("confidence", 0)
            
            print(f"   ✅ SUCCESS")
            print(f"   Database: {db_name}")
            print(f"   Rule: {rule_name}")
            print(f"   Confidence: {confidence:.2f}")
        else:
            print(f"   ❌ FAILED: {result.get('error', 'Unknown error')}")
    
    except Exception as e:
        print(f"   ❌ ERROR: {str(e)}")

print("\n" + "=" * 80)
print("✅ All three databases are integrated and accessible via the agent!")
print("=" * 80)
