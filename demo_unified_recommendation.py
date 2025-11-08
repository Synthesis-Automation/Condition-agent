"""
Quick demo: Unified recommendation showing protocols + rules together
"""

# Simulated unified API (proof of concept)

def unified_recommend(reaction_smiles: str):
    """
    BEFORE (2 separate calls):
        protocols = protocol_recommend(reaction_smiles)
        rules = rule_recommend(reaction_smiles)
        # User has to combine manually
    
    AFTER (1 unified call):
        result = unified_recommend(reaction_smiles)
        # Automatically returns both, ranked by similarity
    """
    
    print("=" * 80)
    print("UNIFIED RECOMMENDATION DEMO")
    print("=" * 80)
    print(f"\nQuery: {reaction_smiles[:60]}...\n")
    
    # Single DRFP computation (used for both protocols and rules)
    print("🔍 Computing DRFP fingerprint...")
    print("✅ DRFP computed (2048-bit vector)\n")
    
    # Single similarity search across ALL sources
    print("📊 Searching unified index...")
    print("   - 250 protocols")
    print("   - 9 rule databases (with reference reactions)")
    print("   - Total: 259 condition sources\n")
    
    print("🎯 Ranking by similarity...\n")
    
    # Unified results
    print("=" * 80)
    print("RESULTS")
    print("=" * 80)
    
    print("\n📖 PROTOCOLS (Detailed Experimental Procedures)")
    print("-" * 80)
    protocols = [
        {"rank": 1, "sim": 0.782, "title": "Sonogashira Coupling of Aryl Iodides (Org. Synth. 2020)"},
        {"rank": 2, "sim": 0.695, "title": "Pd-Catalyzed Alkynylation Protocol (JACS 2019)"},
        {"rank": 3, "sim": 0.621, "title": "Copper-Free Sonogashira (Tetrahedron 2021)"},
    ]
    
    for p in protocols:
        print(f"  {p['rank']}. [Similarity: {p['sim']:.3f}] {p['title']}")
        print(f"     → Full experimental details, amounts, workup, NMR data")
    
    print("\n📋 RULES (Condition Guidelines with Modifiers)")
    print("-" * 80)
    rules = [
        {"rank": 1, "sim": 0.508, "db": "sonogashira_db", "rule": "Aryl iodides/bromides (standard)"},
        {"rank": 2, "sim": 0.114, "db": "Suzuki_db", "rule": "Not applicable (different coupling)"},
        {"rank": 3, "sim": 0.028, "db": "RCM_db", "rule": "Not applicable (different reaction)"},
    ]
    
    for r in rules:
        marker = "✅" if r['rank'] == 1 else "❌"
        print(f"  {marker} {r['rank']}. [Similarity: {r['sim']:.3f}] {r['db']}")
        print(f"     Rule: {r['rule']}")
        if r['rank'] == 1:
            print(f"     → Pd: 0.5-1.5 mol%, Cu: 2-5 mol%, Base: Et3N, Temp: 40-70°C")
            print(f"     → Modifiers: tert-butyl detected → consider 60-80°C")
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("✅ Single API call returned both protocols AND rules")
    print("✅ Ranked by same similarity metric (consistent)")
    print("✅ User gets comprehensive view in one response")
    print("=" * 80)


if __name__ == '__main__':
    # Example: User's Sonogashira reaction
    reaction = "Ic1ccccc1C(C)(C)C.C#Cc1ccccc1>>CC(C)(C)c1ccccc1C#Cc1ccccc1"
    unified_recommend(reaction)
    
    print("\n\n")
    print("🎯 KEY INSIGHT:")
    print("-" * 80)
    print("Instead of maintaining 2 separate systems:")
    print("  1. Protocol recommendation (chemtools/protocol/)")
    print("  2. Rule recommendation (chemtools/rule/ + detection logic)")
    print()
    print("We have 1 unified system:")
    print("  - Single DRFP-based search")
    print("  - Protocols = detailed procedures")
    print("  - Rules = condition guidelines")
    print("  - Same ranking, same filters, same API")
    print()
    print("Result: 50% less code, easier maintenance, better UX")
    print("=" * 80)
