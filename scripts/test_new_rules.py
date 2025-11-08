#!/usr/bin/env python
"""
Test script for new rule database files.

Tests C-O coupling, RCM, and Sonogashira rule databases with sample reactions.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.rule.engine import RuleEngine


def test_c_o_coupling():
    """Test C-O coupling rule database."""
    print("\n" + "="*70)
    print("Testing C-O Coupling Database")
    print("="*70)
    
    # Sample C-O coupling: phenol + aryl bromide
    reaction = "Oc1ccccc1.Brc1ccccc1>>c1ccc(Oc2ccccc2)cc1"
    
    try:
        engine = RuleEngine.from_file("data/rule_db/C_O_coupling_db.json")
        print(f"\n✅ Loaded C_O_coupling_db.json")
        print(f"   Database: {engine.database.metadata.get('name')}")
        print(f"   Version: {engine.database.metadata.get('version')}")
        print(f"   Base rules: {len(engine.database.base_rules)}")
        print(f"   Modifiers: {len(engine.database.modifiers)}")
        
        # Test reaction
        print(f"\nTesting reaction: {reaction}")
        rec = engine.recommend(reaction, include_reasoning=True)
        
        print(f"\n✅ Recommendation generated:")
        print(f"   Rule: {rec.base_rule.name}")
        print(f"   Confidence: {rec.base_rule.confidence}")
        if rec.base_rule.conditions:
            print(f"   Conditions preview:")
            for key, value in list(rec.base_rule.conditions.items())[:5]:
                print(f"      {key}: {value}")
        
        return True
    except Exception as e:
        print(f"\n❌ Error testing C-O coupling: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_rcm():
    """Test RCM rule database."""
    print("\n" + "="*70)
    print("Testing RCM (Ring-Closing Metathesis) Database")
    print("="*70)
    
    # Sample RCM: diene to cyclic alkene
    reaction = "C=CCCC=C>>C1=CCCC1"
    
    try:
        engine = RuleEngine.from_file("data/rule_db/RCM_db.json")
        print(f"\n✅ Loaded RCM_db.json")
        print(f"   Database: {engine.database.metadata.get('name')}")
        print(f"   Version: {engine.database.metadata.get('version')}")
        print(f"   Base rules: {len(engine.database.base_rules)}")
        print(f"   Modifiers: {len(engine.database.modifiers)}")
        
        # Test reaction
        print(f"\nTesting reaction: {reaction}")
        rec = engine.recommend(reaction, include_reasoning=True)
        
        print(f"\n✅ Recommendation generated:")
        print(f"   Rule: {rec.base_rule.name}")
        print(f"   Confidence: {rec.base_rule.confidence}")
        if rec.base_rule.conditions:
            print(f"   Conditions preview:")
            for key, value in list(rec.base_rule.conditions.items())[:5]:
                print(f"      {key}: {value}")
        
        return True
    except Exception as e:
        print(f"\n❌ Error testing RCM: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_sonogashira():
    """Test Sonogashira rule database."""
    print("\n" + "="*70)
    print("Testing Sonogashira Coupling Database")
    print("="*70)
    
    # Sample Sonogashira: aryl bromide + terminal alkyne
    reaction = "Brc1ccccc1.C#C>>C#Cc1ccccc1"
    
    try:
        engine = RuleEngine.from_file("data/rule_db/sonogashira_db.json")
        print(f"\n✅ Loaded sonogashira_db.json")
        print(f"   Database: {engine.database.metadata.get('name')}")
        print(f"   Version: {engine.database.metadata.get('version')}")
        print(f"   Base rules: {len(engine.database.base_rules)}")
        print(f"   Modifiers: {len(engine.database.modifiers)}")
        
        # Test reaction
        print(f"\nTesting reaction: {reaction}")
        rec = engine.recommend(reaction, include_reasoning=True)
        
        print(f"\n✅ Recommendation generated:")
        print(f"   Rule: {rec.base_rule.name}")
        print(f"   Confidence: {rec.base_rule.confidence}")
        if rec.base_rule.conditions:
            print(f"   Conditions preview:")
            for key, value in list(rec.base_rule.conditions.items())[:5]:
                print(f"      {key}: {value}")
        
        return True
    except Exception as e:
        print(f"\n❌ Error testing Sonogashira: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("="*70)
    print("Testing New Rule Database Files")
    print("="*70)
    
    results = []
    
    # Test each database
    results.append(("C-O Coupling", test_c_o_coupling()))
    results.append(("RCM", test_rcm()))
    results.append(("Sonogashira", test_sonogashira()))
    
    # Summary
    print("\n" + "="*70)
    print("Test Summary")
    print("="*70)
    
    for name, passed in results:
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"{name:20s} {status}")
    
    all_passed = all(passed for _, passed in results)
    
    if all_passed:
        print("\n✅ All tests passed!")
        sys.exit(0)
    else:
        print("\n❌ Some tests failed. See details above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
