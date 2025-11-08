#!/usr/bin/env python
"""
Example: Using the new Sonogashira rule database for the selected reaction.

This demonstrates how to get specific reaction conditions for:
Ic1ccccc1C(C)(C)C.C#Cc1ccccc1>>CC(C)(C)c1ccccc1C#Cc1ccccc1
(Ortho-tert-butyl iodobenzene + phenylacetylene → Sonogashira product)
"""

from chemtools.rule.engine import RuleEngine

def main():
    # The reaction from sample_reactions.py line 1169
    reaction = "Ic1ccccc1C(C)(C)C.C#Cc1ccccc1>>CC(C)(C)c1ccccc1C#Cc1ccccc1"
    
    print("="*70)
    print("Sonogashira Coupling Condition Recommendation")
    print("="*70)
    print(f"\nReaction: {reaction}")
    print("\nDescription: Ortho-tert-butyl iodobenzene + phenylacetylene")
    print("             (Sterically hindered substrate)")
    
    # Load the Sonogashira database
    engine = RuleEngine.from_file("data/rule_db/sonogashira_db.json")
    
    # Get recommendation
    rec = engine.recommend(reaction, include_reasoning=True)
    
    # Print formatted summary
    print(f"\n{rec.format_summary()}")
    
    # Access specific conditions programmatically
    print("\n" + "="*70)
    print("Key Conditions for Your Reaction:")
    print("="*70)
    
    conditions = rec.base_rule.conditions
    
    print(f"\n1. Pd Precatalyst: {conditions.get('pd_precatalyst')}")
    print(f"   Loading: {conditions.get('pd_loading_molpct')} mol%")
    
    print(f"\n2. Cu Cocatalyst: {conditions.get('cu_cocatalyst')}")
    
    print(f"\n3. Base: {conditions.get('base')}")
    
    print(f"\n4. Solvent: {conditions.get('solvent')}")
    
    print(f"\n5. Temperature: {conditions.get('temperature_C')} °C")
    
    # Show any applied modifiers
    if rec.modifiers:
        print("\n" + "="*70)
        print("Additional Recommendations (Modifiers):")
        print("="*70)
        for mod in rec.modifiers:
            print(f"\n• {mod.suggestion}")
    
    # Show detected features
    if rec.detected_features:
        print("\n" + "="*70)
        print("Detected Reaction Features:")
        print("="*70)
        for key, value in sorted(rec.detected_features.items()):
            if isinstance(value, bool) and value:
                print(f"  ✓ {key}")
    
    print("\n" + "="*70)
    print("Notes:")
    print("="*70)
    print("""
- This is an aryl iodide with ortho-tert-butyl group (sterically hindered)
- Standard Sonogashira conditions should work well
- Aryl iodides are highly reactive in Sonogashira coupling
- The ortho substitution may slightly slow the reaction
- Consider degassing thoroughly to avoid Glaser homocoupling
- Oxygen promotes diyne formation (Glaser side reaction)
    """)
    
    return rec


if __name__ == "__main__":
    rec = main()
