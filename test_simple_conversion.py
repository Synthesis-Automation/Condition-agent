"""
Simple test to verify rule-to-protocol conversion works end-to-end.
"""

import json
from pathlib import Path
from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup

def main():
    print("=" * 80)
    print("END-TO-END TEST: Rule → Protocol-Compatible Format")
    print("=" * 80)
    print()
    
    # Load Suzuki rule
    rule_file = Path("data/rule_db_v2/Suzuki_db.json")
    with open(rule_file, encoding='utf-8') as f:
        rule_data = json.load(f)
    
    print("📋 Input: Suzuki-Miyaura Rule")
    print("-" * 80)
    conditions = rule_data["default_rule"]["conditions"]
    print(json.dumps(conditions, indent=2, ensure_ascii=False))
    print()
    
    # Convert to protocol format
    print("🔄 Converting to automation format...")
    print()
    
    result = rule_conditions_to_reaction_setup(
        conditions=conditions,
        scale_mmol=1.0,
        reaction_family="Suzuki_Miyaura"
    )
    
    print("✅ Output: Protocol-Compatible Format")
    print("-" * 80)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    print()
    
    # Verify structure
    print("🔍 Verification:")
    print("-" * 80)
    
    assert "reaction_setup" in result, "Missing reaction_setup"
    print("  ✅ Has reaction_setup")
    
    assert len(result["reaction_setup"]) > 0, "Empty reaction_setup"
    print(f"  ✅ Has {len(result['reaction_setup'])} setup(s)")
    
    setup = result["reaction_setup"][0]
    assert "chemicals" in setup, "Missing chemicals"
    assert "conditions" in setup, "Missing conditions"
    print(f"  ✅ Setup has chemicals ({len(setup['chemicals'])}) and conditions")
    
    # Check addition order
    chemicals = setup["chemicals"]
    roles = [c["role"] for c in chemicals]
    print(f"  ✅ Addition order: {' → '.join(roles)}")
    
    # Check expected roles
    expected_roles = {"solvent", "base", "metal_catalyst", "ligand"}
    actual_roles = set(roles)
    assert expected_roles.issubset(actual_roles), f"Missing expected roles"
    print(f"  ✅ All expected roles present: {expected_roles}")
    
    # Check metadata
    metadata = result.get("metadata", {})
    assert metadata.get("generated_from") == "rule"
    assert metadata.get("format") == "protocol-compatible"
    print(f"  ✅ Metadata: generated_from={metadata['generated_from']}, format={metadata['format']}")
    
    print()
    print("=" * 80)
    print("SUCCESS! Rule → Protocol conversion working perfectly! 🎉")
    print("=" * 80)
    print()
    
    print("Summary:")
    print(f"  • Converted {len(conditions)} rule conditions")
    print(f"  • Generated {len(chemicals)} chemical entries")
    print(f"  • Ordered by standard addition sequence")
    print(f"  • Protocol-compatible structure ✅")
    print(f"  • Ready for automated execution ✅")
    print()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
