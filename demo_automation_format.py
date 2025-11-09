"""
Demonstration of the complete automation format workflow.

This script showcases the end-to-end automation format feature:
1. Loading rules
2. Converting to protocol format
3. Using with UnifiedRecommender
4. Output for automated systems
"""

import json
from pathlib import Path
from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup

def demo_basic_conversion():
    """Demo 1: Basic rule-to-protocol conversion."""
    print("\n" + "=" * 80)
    print("DEMO 1: Basic Rule-to-Protocol Conversion")
    print("=" * 80)
    
    # Simple Suzuki conditions
    conditions = {
        "pd_source": "Pd(OAc)2",
        "ligand": "PPh3",
        "catalyst_loading_molpct": "5",
        "base": "K2CO3",
        "base_equiv": "3",
        "solvent": "DMF",
        "temperature_C": "80",
        "time_h": "12"
    }
    
    print("\n📝 Input (Rule Format):")
    print(json.dumps(conditions, indent=2))
    
    # Convert
    result = rule_conditions_to_reaction_setup(
        conditions=conditions,
        scale_mmol=2.0,  # 2 mmol scale
        reaction_family="Suzuki_Miyaura"
    )
    
    print("\n🔄 Output (Automation Format):")
    print(json.dumps(result, indent=2))
    
    print("\n✨ Key Features:")
    chemicals = result["reaction_setup"][0]["chemicals"]
    print(f"  • {len(chemicals)} chemicals in proper addition order")
    print(f"  • Amounts calculated for 2.0 mmol scale")
    print(f"  • Protocol-compatible structure")
    
    return result


def demo_addition_sequence():
    """Demo 2: Addition sequence extraction."""
    print("\n" + "=" * 80)
    print("DEMO 2: Addition Sequence for Lab/Robot")
    print("=" * 80)
    
    # Load Sonogashira rule
    rule_file = Path("data/rule_db_v2/sonogashira_db.json")
    with open(rule_file, encoding='utf-8') as f:
        rule_data = json.load(f)
    
    conditions = rule_data["default_rule"]["conditions"]
    
    # Convert
    result = rule_conditions_to_reaction_setup(
        conditions=conditions,
        scale_mmol=5.0,
        reaction_family="Sonogashira"
    )
    
    print("\n🤖 Automated Addition Sequence (5 mmol scale):")
    print("-" * 80)
    
    chemicals = result["reaction_setup"][0]["chemicals"]
    for i, chem in enumerate(chemicals, 1):
        name = chem["name"]
        role = chem["role"]
        amount = chem["amount"]
        
        print(f"\nStep {i}: Add {name}")
        print(f"  Role: {role}")
        
        if "equivalents" in amount:
            equiv = amount["equivalents"]
            mmol = amount.get("mmol", equiv * 5.0)
            print(f"  Amount: {equiv} equiv ({mmol:.2f} mmol)")
        elif "volume_ml" in amount:
            print(f"  Volume: As needed for target concentration")
        
        if "note" in amount:
            print(f"  Note: {amount['note']}")
    
    # Conditions
    conditions_list = result["reaction_setup"][0]["conditions"]
    if conditions_list:
        cond = conditions_list[0]
        print(f"\n⚗️  Reaction Conditions:")
        print(f"  Temperature: {cond.get('temperature_C', 'N/A')}°C")
        print(f"  Time: {cond.get('time_h', 'N/A')} hours")
        print(f"  Atmosphere: {cond.get('atmosphere', 'N/A')}")


def demo_multiple_families():
    """Demo 3: Compare different reaction families."""
    print("\n" + "=" * 80)
    print("DEMO 3: Multiple Reaction Families")
    print("=" * 80)
    
    families = [
        ("Suzuki_db.json", "Suzuki_Miyaura"),
        ("sonogashira_db.json", "Sonogashira"),
        ("heck_db.json", "Heck")
    ]
    
    print("\n📊 Addition Sequence Comparison:")
    print("-" * 80)
    
    for filename, family in families:
        rule_file = Path(f"data/rule_db_v2/{filename}")
        
        if not rule_file.exists():
            print(f"\n{family}: File not found")
            continue
        
        with open(rule_file, encoding='utf-8') as f:
            rule_data = json.load(f)
        
        conditions = rule_data["default_rule"]["conditions"]
        result = rule_conditions_to_reaction_setup(
            conditions=conditions,
            scale_mmol=1.0,
            reaction_family=family
        )
        
        chemicals = result["reaction_setup"][0]["chemicals"]
        roles = [c["role"] for c in chemicals]
        
        print(f"\n{family}:")
        print(f"  Sequence: {' → '.join(roles)}")
        print(f"  Total steps: {len(chemicals)}")


def demo_robot_integration():
    """Demo 4: Example robot integration."""
    print("\n" + "=" * 80)
    print("DEMO 4: Robot Integration Example")
    print("=" * 80)
    
    # Simulated robot controller
    class RobotController:
        def __init__(self):
            self.steps = []
        
        def add_chemical(self, name, role, mmol=None, note=None):
            step = {"name": name, "role": role}
            if mmol:
                step["mmol"] = mmol
            if note:
                step["note"] = note
            self.steps.append(step)
            print(f"  ✓ Queued: {name} ({mmol} mmol)")
        
        def set_conditions(self, temperature_C, time_h, atmosphere):
            self.conditions = {
                "temperature": temperature_C,
                "time": time_h,
                "atmosphere": atmosphere
            }
            print(f"  ✓ Conditions set: {temperature_C}°C, {time_h}h, {atmosphere}")
        
        def execute(self):
            print(f"\n  🚀 Executing {len(self.steps)} addition steps...")
            print(f"  ⚗️  Running reaction at {self.conditions['temperature']}°C for {self.conditions['time']}h")
            print(f"  ✅ Reaction complete!")
    
    # Load conditions
    rule_file = Path("data/rule_db_v2/Suzuki_db.json")
    with open(rule_file, encoding='utf-8') as f:
        rule_data = json.load(f)
    
    conditions = rule_data["default_rule"]["conditions"]
    
    # Convert to automation format
    result = rule_conditions_to_reaction_setup(
        conditions=conditions,
        scale_mmol=10.0,  # 10 mmol scale
        reaction_family="Suzuki_Miyaura"
    )
    
    print("\n🤖 Programming Robot for 10 mmol Suzuki Reaction:")
    print("-" * 80)
    
    # Send to robot
    robot = RobotController()
    
    setup = result["reaction_setup"][0]
    
    # Add chemicals in order
    for chem in setup["chemicals"]:
        mmol = chem["amount"].get("mmol")
        note = chem["amount"].get("note")
        robot.add_chemical(
            name=chem["name"],
            role=chem["role"],
            mmol=mmol,
            note=note
        )
    
    # Set conditions
    cond = setup["conditions"][0]
    robot.set_conditions(
        temperature_C=cond["temperature_C"],
        time_h=cond["time_h"],
        atmosphere=cond["atmosphere"]
    )
    
    # Execute
    robot.execute()


def main():
    """Run all demos."""
    print("\n" + "🎯 " * 30)
    print("AUTOMATION FORMAT DEMONSTRATION")
    print("🎯 " * 30)
    
    try:
        demo_basic_conversion()
        demo_addition_sequence()
        demo_multiple_families()
        demo_robot_integration()
        
        print("\n" + "=" * 80)
        print("✅ ALL DEMOS COMPLETE!")
        print("=" * 80)
        
        print("\n📚 What We Demonstrated:")
        print("  1. Basic rule-to-protocol conversion")
        print("  2. Extraction of ordered addition sequences")
        print("  3. Comparison across reaction families")
        print("  4. Robot integration workflow")
        
        print("\n🚀 This format is ready for:")
        print("  • Robotic synthesis systems")
        print("  • LLM-guided chemistry")
        print("  • Automated protocol generation")
        print("  • Lab workflow optimization")
        
        print("\n📖 Learn more:")
        print("  • docs/AUTOMATION_FORMAT.md - Full documentation")
        print("  • docs/AUTOMATION_QUICKSTART.md - Quick start guide")
        print("  • test_*.py - Example code")
        
    except Exception as e:
        print(f"\n❌ Demo failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
