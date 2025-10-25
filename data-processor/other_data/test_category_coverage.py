"""
Test category-level SMARTS coverage against sample reactions.

This tests whether:
1. Category-level patterns match more reactants than member-level alone
2. Category patterns properly encompass all their members
3. Hierarchical matching improves classification
"""

import sys
import re
from collections import defaultdict

sys.path.append('data-processor/other_data')
sys.path.append('tests')

from classify_reactant import load_reactant_types, classify_reactant, get_category_matches
from sample_reactions import SAMPLE_REACTIONS


def test_category_coverage():
    """Test category-level matching against sample reactions."""
    
    # Parse reaction SMILES from SAMPLE_REACTIONS
    reactants = set()
    for reaction_smiles in SAMPLE_REACTIONS:
        if reaction_smiles == "Select a sample reaction...":
            continue
        if '>>' not in reaction_smiles:
            continue
        
        # Extract reactants from reaction SMILES
        match = re.match(r'([^>]+)>>', reaction_smiles)
        if not match:
            continue
        
        reactant_side = match.group(1)
        for smiles in reactant_side.split('.'):
            smiles = smiles.strip()
            # Filter out simple reagents
            if smiles and smiles not in ['[H][H]', '[O]', '[Zn]', 'O']:
                reactants.add(smiles)
    
    reactants = sorted(list(reactants))
    print(f"Testing {len(reactants)} unique reactants\n")
    
    reactant_types = load_reactant_types()
    
    # Track statistics
    member_matched = 0
    category_matched = 0
    category_only_matched = 0
    both_matched = 0
    neither_matched = 0
    
    category_matches_count = defaultdict(int)
    member_matches_count = defaultdict(int)
    category_only_examples = []
    
    print("=" * 80)
    print("Testing Category vs Member Coverage")
    print("=" * 80)
    
    for smiles in reactants:
        # Test member-level classification
        member_result = classify_reactant(smiles, reactant_types)
        has_member_match = member_result is not None
        
        # Test category-level matching
        categories = get_category_matches(smiles, reactant_types)
        has_category_match = len(categories) > 0
        
        # Count matches
        if has_member_match:
            member_matched += 1
            member_matches_count[member_result['category']] += 1
        
        if has_category_match:
            category_matched += 1
            for cat in categories:
                category_matches_count[cat] += 1
        
        # Classify match types
        if has_member_match and has_category_match:
            both_matched += 1
        elif has_category_match and not has_member_match:
            category_only_matched += 1
            category_only_examples.append({
                'smiles': smiles,
                'categories': categories
            })
            if len(category_only_examples) <= 5:  # Show first 5
                print(f"\n⚠️  Category-only match: {smiles}")
                print(f"   Categories: {', '.join(categories)}")
        elif has_member_match and not has_category_match:
            print(f"\n⚠️  UNEXPECTED: Member match but no category: {smiles}")
            print(f"   Member: {member_result['member_type']}")
        else:
            neither_matched += 1
    
    # Print summary
    print("\n" + "=" * 80)
    print("COVERAGE SUMMARY")
    print("=" * 80)
    print(f"Total reactants:              {len(reactants)}")
    print(f"Member-level matches:         {member_matched} ({member_matched/len(reactants)*100:.1f}%)")
    print(f"Category-level matches:       {category_matched} ({category_matched/len(reactants)*100:.1f}%)")
    print(f"Both matched:                 {both_matched}")
    print(f"Category-only (no member):    {category_only_matched}")
    print(f"Neither matched:              {neither_matched}")
    
    # Category improvement
    improvement = category_matched - member_matched
    print(f"\n📈 Category-level adds {improvement} matches ({improvement/len(reactants)*100:.1f}% improvement)")
    
    # Top categories
    print("\n" + "=" * 80)
    print("TOP 10 CATEGORIES (by coverage)")
    print("=" * 80)
    sorted_cats = sorted(category_matches_count.items(), key=lambda x: x[1], reverse=True)
    for i, (cat, count) in enumerate(sorted_cats[:10], 1):
        member_count = member_matches_count.get(cat, 0)
        print(f"{i:2}. {cat:25s} {count:3} category matches, {member_count:3} member matches")
    
    # Show category-only examples
    if category_only_examples:
        print("\n" + "=" * 80)
        print(f"CATEGORY-ONLY MATCHES ({len(category_only_examples)} total)")
        print("=" * 80)
        print("These match at category level but no specific member:")
        for ex in category_only_examples[:10]:
            print(f"  {ex['smiles']:40s} → {', '.join(ex['categories'])}")
    
    # Validation check
    print("\n" + "=" * 80)
    print("VALIDATION")
    print("=" * 80)
    if category_only_matched == 0:
        print("✅ Perfect: All member matches have corresponding category matches")
        print("✅ No category-only matches (all categories have specific members)")
    elif category_only_matched < 10:
        print(f"⚠️  {category_only_matched} category-only matches found")
        print("   This suggests some category patterns are broader than member patterns")
    else:
        print(f"❌ {category_only_matched} category-only matches")
        print("   Consider adding more specific member types or refining category patterns")
    
    return {
        'total': len(reactants),
        'member_matched': member_matched,
        'category_matched': category_matched,
        'category_only': category_only_matched,
        'improvement': improvement
    }


if __name__ == "__main__":
    results = test_category_coverage()
