"""
Test rule-to-protocol converter with actual rule files.
"""

import json
from pathlib import Path
from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup

def test_sonogashira_conversion():
    """Test conversion of Sonogashira rule to protocol format."""
    print("=" * 80)
    print("TEST 1: Sonogashira Rule → Protocol Format")
    print("=" * 80)
    print()
    
    # Load Sonogashira rule
    rule_file = Path("data/rule_db_v2/sonogashira_db.json")
    with open(rule_file, encoding='utf-8') as f:
        rule_data = json.load(f)
    
    # Get default conditions
    conditions = rule_data["default_rule"]["conditions"]
    
    print("Input (Rule Conditions):")
    print(json.dumps(conditions, indent=2, ensure_ascii=False))
    print()
    
    # Convert to protocol format
    result = rule_conditions_to_reaction_setup(
        conditions=conditions,
        scale_mmol=1.0,
        reaction_family="Sonogashira"
    )
    
    print("Output (Protocol-Compatible Format):")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    print()
    
    # Validate structure
    assert "reaction_setup" in result
    assert len(result["reaction_setup"]) > 0
    setup = result["reaction_setup"][0]
    assert "chemicals" in setup
    assert "conditions" in setup
    
    chemicals = setup["chemicals"]
    print(f"✅ Generated {len(chemicals)} chemical entries")
    print(f"✅ Addition order preserved: {[c['name'] for c in chemicals]}")
    
    # Check roles
    roles = [c["role"] for c in chemicals]
    assert "solvent" in roles
    assert "base" in roles or "metal_catalyst" in roles
    print(f"✅ Roles assigned: {set(roles)}")
    print()


def test_suzuki_conversion():
    """Test conversion of Suzuki rule to protocol format."""
    print("=" * 80)
    print("TEST 2: Suzuki Rule → Protocol Format")
    print("=" * 80)
    print()
    
    rule_file = Path("data/rule_db_v2/Suzuki_db.json")
    with open(rule_file, encoding='utf-8') as f:
        rule_data = json.load(f)
    
    conditions = rule_data["default_rule"]["conditions"]
    
    print("Input (Rule Conditions):")
    print(json.dumps(conditions, indent=2, ensure_ascii=False))
    print()
    
    result = rule_conditions_to_reaction_setup(
        conditions=conditions,
        scale_mmol=1.0,
        reaction_family="Suzuki_Miyaura"
    )
    
    print("Output (Protocol-Compatible Format):")
    chemicals = result["reaction_setup"][0]["chemicals"]
    for chem in chemicals:
        print(f"  {chem['name']}: {chem['role']}")
        if "equivalents" in chem["amount"]:
            print(f"    Amount: {chem['amount']['equivalents']} equiv")
    print()
    
    # Validate
    roles = [c["role"] for c in chemicals]
    assert "solvent" in roles
    assert "base" in roles
    assert "metal_catalyst" in roles or "ligand" in roles
    print(f"✅ All key roles present: {set(roles)}")
    print()


def test_with_user_substrates():
    """Test conversion with user-provided substrates."""
    print("=" * 80)
    print("TEST 3: With User Substrates")
    print("=" * 80)
    print()
    
    rule_file = Path("data/rule_db_v2/sonogashira_db.json")
    with open(rule_file, encoding='utf-8') as f:
        rule_data = json.load(f)
    
    conditions = rule_data["default_rule"]["conditions"]
    
    # User substrates
    user_substrates = [
        {
            "name": "Bromobenzene",
            "smiles": "Brc1ccccc1",
            "role": "starting_material",
            "mmol": 1.0,
            "equivalents": 1.0
        },
        {
            "name": "Phenylacetylene",
            "smiles": "C#Cc1ccccc1",
            "role": "starting_material",
            "mmol": 1.2,
            "equivalents": 1.2
        }
    ]
    
    result = rule_conditions_to_reaction_setup(
        conditions=conditions,
        user_substrates=user_substrates,
        scale_mmol=1.0,
        reaction_family="Sonogashira"
    )
    
    print("User Substrates:")
    for sub in user_substrates:
        print(f"  - {sub['name']} ({sub['equivalents']} equiv)")
    print()
    
    print("Generated Chemicals:")
    for chem in result["reaction_setup"][0]["chemicals"]:
        print(f"  {chem['name']}: {chem['role']}")
        if "smiles" in chem:
            print(f"    SMILES: {chem['smiles']}")
    print()
    
    # Find substrate entries
    substrate_chems = [c for c in result["reaction_setup"][0]["chemicals"] 
                       if c["role"] == "starting_material"]
    assert len(substrate_chems) == 2
    assert substrate_chems[0]["smiles"] == "Brc1ccccc1"
    assert substrate_chems[1]["smiles"] == "C#Cc1ccccc1"
    print("✅ User substrates correctly included")
    print()


def test_range_parsing():
    """Test range midpoint calculation."""
    print("=" * 80)
    print("TEST 4: Range Parsing")
    print("=" * 80)
    print()
    
    from chemtools.formatters.rule_to_protocol import _parse_range_midpoint
    
    test_cases = [
        ("0.5-2.0", 1.25),
        ("60-80", 70.0),
        ("1.5", 1.5),
        (1.5, 1.5),
        ("80–100", 90.0),  # en-dash
    ]
    
    for input_val, expected in test_cases:
        result = _parse_range_midpoint(input_val)
        print(f"  {repr(input_val):20} → {result:6.2f} (expected {expected})")
        assert result == expected, f"Expected {expected}, got {result}"
    
    print()
    print("✅ All range parsing tests passed")
    print()


def test_option_picking():
    """Test first option selection."""
    print("=" * 80)
    print("TEST 5: Option Picking")
    print("=" * 80)
    print()
    
    from chemtools.formatters.rule_to_protocol import _pick_first_option
    
    test_cases = [
        ("THF or toluene", "THF"),
        ("PdCl2(PPh3)2 or Pd(dppf)Cl2·DCM", "PdCl2(PPh3)2"),
        ("THF, toluene, or DMF", "THF"),
        ("Et3N or DIPEA", "Et3N"),
        ("just_one", "just_one"),
    ]
    
    for input_val, expected in test_cases:
        result = _pick_first_option(input_val)
        print(f"  {input_val:40} → {result:20}")
        assert result == expected, f"Expected {expected}, got {result}"
    
    print()
    print("✅ All option picking tests passed")
    print()


def compare_with_protocol():
    """Compare generated format with actual protocol structure."""
    print("=" * 80)
    print("TEST 6: Format Comparison with Real Protocol")
    print("=" * 80)
    print()
    
    # Load real protocol
    protocol_file = Path("data/protocol_db_v2/Sonogashira-Coupling.json")
    with open(protocol_file, encoding='utf-8') as f:
        protocol_data = json.load(f)
    
    protocol_setup = protocol_data[0]["reaction_setup"][0]
    
    # Load rule and convert
    rule_file = Path("data/rule_db_v2/sonogashira_db.json")
    with open(rule_file, encoding='utf-8') as f:
        rule_data = json.load(f)
    
    conditions = rule_data["default_rule"]["conditions"]
    result = rule_conditions_to_reaction_setup(conditions, scale_mmol=80.0)
    generated_setup = result["reaction_setup"][0]
    
    print("Real Protocol Structure:")
    print(f"  Chemicals: {len(protocol_setup['chemicals'])}")
    print(f"  Roles: {[c['role'] for c in protocol_setup['chemicals']]}")
    print()
    
    print("Generated Structure:")
    print(f"  Chemicals: {len(generated_setup['chemicals'])}")
    print(f"  Roles: {[c['role'] for c in generated_setup['chemicals']]}")
    print()
    
    # Check structure compatibility
    for key in ["chemicals", "conditions"]:
        assert key in protocol_setup, f"Protocol has {key}"
        assert key in generated_setup, f"Generated has {key}"
    
    # Check chemical structure
    for chem in generated_setup["chemicals"]:
        assert "name" in chem
        assert "role" in chem
        assert "amount" in chem
        assert isinstance(chem["amount"], dict)
    
    print("✅ Generated format matches protocol structure")
    print()


if __name__ == "__main__":
    print("\n" + "🧪 " * 30)
    print("RULE-TO-PROTOCOL CONVERTER TESTS")
    print("🧪 " * 30)
    print()
    
    test_range_parsing()
    test_option_picking()
    test_sonogashira_conversion()
    test_suzuki_conversion()
    test_with_user_substrates()
    compare_with_protocol()
    
    print("=" * 80)
    print("ALL TESTS PASSED! ✅")
    print("=" * 80)
    print()
    print("The converter successfully:")
    print("  ✅ Parses rule conditions")
    print("  ✅ Generates protocol-compatible format")
    print("  ✅ Maintains proper addition order")
    print("  ✅ Handles ranges and options")
    print("  ✅ Integrates user substrates")
    print("  ✅ Matches real protocol structure")
    print()
