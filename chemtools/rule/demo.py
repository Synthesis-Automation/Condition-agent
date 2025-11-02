#!/usr/bin/env python
"""
Demonstration of Rule-Based Recommendation System
==================================================

Shows various use cases for the rule-based recommendation engine.
"""

from chemtools.rule import RuleEngine

def demo_basic_suzuki():
    """Demo 1: Basic Suzuki coupling recommendation."""
    print("\n" + "="*70)
    print("DEMO 1: Basic Suzuki Coupling (Bromobenzene + Phenylboronic Acid)")
    print("="*70)
    
    engine = RuleEngine.from_file("data/rule_db/suzuki.json")
    
    # Simple Suzuki coupling
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1"
    
    rec = engine.recommend(reaction)
    print(rec.format_summary())


def demo_activated_chloride():
    """Demo 2: Activated aryl chloride."""
    print("\n" + "="*70)
    print("DEMO 2: Activated Aryl Chloride (4-Chlorobenzoic Acid)")
    print("="*70)
    
    engine = RuleEngine.from_file("data/rule_db/suzuki.json")
    
    # Activated aryl chloride with EWG
    reaction = "Clc1ccc(C(=O)O)cc1.OB(O)c1ccccc1>>O=C(O)c1ccc(-c2ccccc2)cc1"
    
    rec = engine.recommend(reaction)
    print(rec.format_summary())


def demo_with_symptoms():
    """Demo 3: Recommendation with observed symptoms."""
    print("\n" + "="*70)
    print("DEMO 3: Recommendation with Symptoms (Low Yield)")
    print("="*70)
    
    engine = RuleEngine.from_file("data/rule_db/suzuki.json")
    
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1"
    
    # Simulate observed low yield
    rec = engine.recommend(reaction, symptoms=["low_yield"])
    print(rec.format_summary())


def demo_json_output():
    """Demo 4: JSON output for programmatic use."""
    print("\n" + "="*70)
    print("DEMO 4: JSON Output (Programmatic Access)")
    print("="*70)
    
    engine = RuleEngine.from_file("data/rule_db/suzuki.json")
    
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1"
    
    rec = engine.recommend(reaction, include_reasoning=True)
    rec_dict = rec.to_dict()
    
    import json
    print(json.dumps(rec_dict, indent=2))


def demo_batch_processing():
    """Demo 5: Batch processing multiple reactions."""
    print("\n" + "="*70)
    print("DEMO 5: Batch Processing (3 Reactions)")
    print("="*70)
    
    engine = RuleEngine.from_file("data/rule_db/suzuki.json")
    
    reactions = [
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1",  # Simple bromide
        "Clc1ccc(C(=O)O)cc1.OB(O)c1ccccc1>>O=C(O)c1ccc(-c2ccccc2)cc1",  # Activated chloride
        "Ic1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1",  # Iodide
    ]
    
    recommendations = engine.batch_recommend(reactions, include_reasoning=False)
    
    for i, rec in enumerate(recommendations, 1):
        print(f"\n--- Reaction {i} ---")
        print(f"Rule: {rec.base_rule.name}")
        print(f"Pd source: {rec.base_rule.conditions.get('pd_source', 'N/A')}")
        print(f"Modifiers: {len(rec.modifiers)}")


def demo_feature_analysis():
    """Demo 6: Feature analysis without recommendation."""
    print("\n" + "="*70)
    print("DEMO 6: Feature Analysis (No Recommendation)")
    print("="*70)
    
    from chemtools.rule.analyzer import FeatureAnalyzer
    
    analyzer = FeatureAnalyzer()
    
    # Analyze a complex substrate
    substrate = "Clc1cc(C(=O)OC)c(C(C)(C)C)cc1"
    
    features = analyzer.analyze_reactant(substrate)
    
    # Filter to show only relevant features
    relevant = analyzer.get_relevant_features(features, include_false=False)
    
    print(f"Substrate: {substrate}")
    print(f"\nDetected Features ({len(relevant)}):")
    for name, value in sorted(relevant.items()):
        if isinstance(value, bool):
            print(f"  ✓ {name}")
        else:
            print(f"  • {name} = {value}")


def demo_database_validation():
    """Demo 7: Database validation."""
    print("\n" + "="*70)
    print("DEMO 7: Database Validation")
    print("="*70)
    
    from chemtools.rule.database import RuleDatabase
    
    db = RuleDatabase.from_file("data/rule_db/suzuki.json")
    
    print(f"Database Summary:")
    print(db.get_summary())
    
    print(f"\nValidation:")
    issues = db.validate()
    if issues:
        print("❌ Issues found:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("✅ Database is valid")
    
    print(f"\nBase Rules:")
    for i, rule in enumerate(db.base_rules, 1):
        print(f"  {i}. {rule.name}")
        if rule.reactant_features:
            print(f"     Features: {rule.reactant_features}")
    
    print(f"\nModifiers:")
    for i, mod in enumerate(db.modifiers, 1):
        print(f"  {i}. {mod.suggestion[:60]}...")
        print(f"     When: {', '.join(mod.when[:3])}{'...' if len(mod.when) > 3 else ''}")


def main():
    """Run all demonstrations."""
    print("\n" + "#"*70)
    print("# RULE-BASED RECOMMENDATION SYSTEM - DEMONSTRATIONS")
    print("#"*70)
    
    demos = [
        demo_basic_suzuki,
        demo_activated_chloride,
        demo_with_symptoms,
        demo_batch_processing,
        demo_feature_analysis,
        demo_database_validation,
        # demo_json_output,  # Verbose, skip by default
    ]
    
    for demo in demos:
        try:
            demo()
        except Exception as e:
            print(f"\n❌ Error in {demo.__name__}: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "#"*70)
    print("# DEMONSTRATIONS COMPLETE")
    print("#"*70)
    print("\nFor more examples, see:")
    print("  - chemtools/rule/README.md")
    print("  - tests/test_rule_engine.py")
    print("  - python -m chemtools.rule.cli --help")
    print()


if __name__ == "__main__":
    main()
