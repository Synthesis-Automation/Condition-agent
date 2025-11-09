"""
Detailed problem analysis for sample compounds test
"""
import sys
from pathlib import Path
from collections import defaultdict

tests_dir = Path(__file__).parent / "tests"
sys.path.insert(0, str(tests_dir))

from sample_compounds import ALL_SAMPLE_COMPOUNDS
from chemtools.featurizers.calculable import classify_reactant_smiles, get_reactant_type_features

print("="*80)
print("PROBLEM DETECTION REPORT - Sample Compounds Analysis")
print("="*80)

# Track issues
undetected = []
unexpected_category = []
low_specificity = []
multiple_categories = []
warnings = []

for compound in ALL_SAMPLE_COMPOUNDS:
    smiles = compound['smiles']
    name = compound['name']
    role = compound.get('role', 'none')
    expected_features = compound.get('features', [])
    
    # Classify
    result = classify_reactant_smiles(smiles)
    features = get_reactant_type_features(smiles)
    
    # Issue 1: Undetected compounds with expected reactant roles
    if role in ['electrophile', 'nucleophile', 'bifunctional', 'di-electrophile', 'di-nucleophile']:
        if not result:
            undetected.append({
                'name': name,
                'smiles': smiles,
                'role': role,
                'expected': expected_features
            })
    
    # Issue 2: Check for unexpected categorization
    if result:
        detected_member = result['member_type']
        detected_category = result['category']
        
        # If we have expected features, check if detection matches
        if expected_features:
            expected_members = [f for f in expected_features if not f.endswith('_present')]
            if expected_members and detected_member not in expected_members:
                # Check if it's a more specific match
                if detected_member not in ['ArBr', 'ArCl', 'ArI', 'ArF', 'Bn-Br', 'Bn-Cl']:
                    unexpected_category.append({
                        'name': name,
                        'smiles': smiles,
                        'expected': expected_members,
                        'detected': detected_member,
                        'category': detected_category
                    })
        
        # Issue 3: Low specificity matches (general categories)
        if result['specificity'] < 5 and detected_category in ['Alkyl-C-H', 'ArH']:
            low_specificity.append({
                'name': name,
                'smiles': smiles,
                'member': detected_member,
                'category': detected_category,
                'specificity': result['specificity']
            })
        
        # Issue 4: Multiple categories detected (potential ambiguity)
        categories = features.get('categories', [])
        if len(categories) > 2:
            multiple_categories.append({
                'name': name,
                'smiles': smiles,
                'categories': categories,
                'members': features.get('member_types', [])
            })

# Print findings
print("\n🔴 CRITICAL ISSUES")
print("-" * 80)

if undetected:
    print(f"\n1. UNDETECTED COMPOUNDS ({len(undetected)} found):")
    print("   Compounds with expected reactant roles but no detection")
    for item in undetected[:5]:  # Show first 5
        print(f"   - {item['name']}")
        print(f"     Role: {item['role']}, SMILES: {item['smiles']}")
        if item['expected']:
            print(f"     Expected: {', '.join(item['expected'])}")
else:
    print("\n✅ 1. No undetected compounds with expected reactant roles")

print("\n🟡 WARNINGS")
print("-" * 80)

if unexpected_category:
    print(f"\n2. UNEXPECTED CATEGORIZATION ({len(unexpected_category)} found):")
    print("   Detection doesn't match expected features")
    for item in unexpected_category[:5]:
        print(f"   - {item['name']}")
        print(f"     Expected: {', '.join(item['expected'])}")
        print(f"     Detected: {item['detected']} ({item['category']})")
else:
    print("\n✅ 2. No unexpected categorizations")

if low_specificity:
    print(f"\n3. LOW SPECIFICITY MATCHES ({len(low_specificity)} found):")
    print("   General categories (Alkyl-C-H, ArH) detected")
    for item in low_specificity[:5]:
        print(f"   - {item['name']}: {item['member']} (specificity={item['specificity']})")
else:
    print("\n✅ 3. No low specificity issues")

if multiple_categories:
    print(f"\n4. MULTIPLE CATEGORIES ({len(multiple_categories)} found):")
    print("   Compounds matching >2 categories (potential ambiguity)")
    for item in multiple_categories[:5]:
        print(f"   - {item['name']}")
        print(f"     Categories: {', '.join(item['categories'])}")
        print(f"     Members: {', '.join(item['members'])}")
else:
    print("\n✅ 4. No multiple category issues")

# Summary statistics
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"Total compounds tested: {len(ALL_SAMPLE_COMPOUNDS)}")
print(f"Critical issues: {len(undetected)}")
print(f"Warnings: {len(unexpected_category) + len(low_specificity) + len(multiple_categories)}")

total_issues = len(undetected) + len(unexpected_category) + len(low_specificity) + len(multiple_categories)

if total_issues == 0:
    print("\n🎉 PERFECT! No issues found.")
    print("✅ System performing optimally on all sample compounds.")
elif len(undetected) == 0:
    print(f"\n✅ No critical issues!")
    print(f"⚠️  {len(unexpected_category) + len(low_specificity) + len(multiple_categories)} warnings (minor issues)")
else:
    print(f"\n⚠️  {len(undetected)} critical issue(s) requiring attention")
    print(f"⚠️  {len(unexpected_category) + len(low_specificity) + len(multiple_categories)} additional warnings")

print("="*80)
