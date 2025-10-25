"""
Test SMARTS pattern classification against real reactants from sample_reactions.py

This script extracts all reactants from sample_reactions, classifies them using
the SMARTS patterns, and generates a comprehensive report showing:
1. Classification accuracy
2. Unmatched reactants
3. Ambiguous cases (multiple matches)
4. Coverage by category
"""

import sys
import re
from collections import defaultdict, Counter
sys.path.append('data-processor/other_data')
sys.path.append('tests')

from classify_reactant import classify_reactant, get_all_matches, load_reactant_types
from sample_reactions import SAMPLE_REACTIONS

def parse_reaction_smiles(reaction_string):
    """
    Parse a reaction string to extract reactant SMILES.
    Format: "reactant1.reactant2>>product (Description)"
    """
    if reaction_string == "Select a sample reaction...":
        return []
    
    # Extract SMILES part (before >>)
    match = re.match(r'([^>]+)>>', reaction_string)
    if not match:
        return []
    
    smiles_part = match.group(1)
    
    # Split by . but be careful of disconnected structures
    # Also remove condition markers like [H][H], [O], [Zn], etc.
    reactants = []
    for part in smiles_part.split('.'):
        part = part.strip()
        # Skip simple reagents/catalysts
        if part in ['[H][H]', '[O]', '[BH4-]', '[Na+]', '[AlH4-]', '[Li+]', 
                    '[K+]', '[B-](F)(F)F', '[B-](F)(F)(F)c1ccccc1', '[Zn]', 
                    '[HCl]', '[NH3]', '[NaBH(OAc)3]', '[DIBAL]', '[CrO3]',
                    '[Dess-Martin]', '[IBX]', '[TEMPO]', '[NaOCl]', '[PCC]',
                    '[mCPBA]', '[Ir]', '[Ru]', '[MnO2]']:
            continue
        # Skip simple reagents
        if re.match(r'^\[.*\]$', part):  # Skip bracketed reagents
            continue
        reactants.append(part)
    
    return reactants

def categorize_by_functional_group(result):
    """Categorize a classification result by its group"""
    if not result:
        return 'Unmatched'
    return result['group']

def main():
    print("=" * 100)
    print("SMARTS PATTERN VALIDATION AGAINST SAMPLE REACTIONS")
    print("=" * 100)
    
    # Load reactant types
    reactant_types = load_reactant_types()
    
    # Extract all reactants
    all_reactants = []
    for reaction in SAMPLE_REACTIONS:
        reactants = parse_reaction_smiles(reaction)
        all_reactants.extend(reactants)
    
    print(f"\n📊 Total reactions: {len([r for r in SAMPLE_REACTIONS if r != 'Select a sample reaction...'])}")
    print(f"📊 Total reactants extracted: {len(all_reactants)}")
    print(f"📊 Unique reactants: {len(set(all_reactants))}")
    
    # Classify all unique reactants
    unique_reactants = list(set(all_reactants))
    results = {}
    unmatched = []
    ambiguous = []
    
    for smiles in unique_reactants:
        result = classify_reactant(smiles, reactant_types)
        results[smiles] = result
        
        if result is None:
            unmatched.append(smiles)
        else:
            # Check for ambiguous cases
            all_matches = get_all_matches(smiles, reactant_types)
            # Filter out general C-H and pi-system matches
            specific_matches = [m for m in all_matches if m['group'] not in ['C–H Donors / π-Systems']]
            if len(specific_matches) > 3:  # More than 3 specific matches
                ambiguous.append((smiles, all_matches))
    
    # Statistics
    matched = len(unique_reactants) - len(unmatched)
    match_rate = 100 * matched / len(unique_reactants)
    
    print(f"\n{'=' * 100}")
    print("CLASSIFICATION RESULTS")
    print(f"{'=' * 100}")
    print(f"✅ Matched: {matched}/{len(unique_reactants)} ({match_rate:.1f}%)")
    print(f"❌ Unmatched: {len(unmatched)}/{len(unique_reactants)} ({100-match_rate:.1f}%)")
    print(f"⚠️  Ambiguous (>3 specific matches): {len(ambiguous)}")
    
    # Group by category
    print(f"\n{'=' * 100}")
    print("CLASSIFICATION BY FUNCTIONAL GROUP")
    print(f"{'=' * 100}")
    
    by_group = defaultdict(list)
    for smiles, result in results.items():
        group = categorize_by_functional_group(result)
        by_group[group].append((smiles, result))
    
    for group in sorted(by_group.keys()):
        items = by_group[group]
        print(f"\n{group}: {len(items)} reactants")
        # Show first 5 examples
        for smiles, result in items[:5]:
            if result:
                print(f"  ✓ {smiles[:40]:40s} → {result['member_type']:20s} ({result['name']})")
            else:
                print(f"  ✗ {smiles[:40]:40s} → No match")
        if len(items) > 5:
            print(f"  ... and {len(items)-5} more")
    
    # Detailed category breakdown
    print(f"\n{'=' * 100}")
    print("CLASSIFICATION BY CATEGORY")
    print(f"{'=' * 100}")
    
    by_category = defaultdict(int)
    for smiles, result in results.items():
        if result:
            by_category[result['category']] += 1
    
    for category in sorted(by_category.keys(), key=lambda x: by_category[x], reverse=True):
        count = by_category[category]
        print(f"  {category:25s}: {count:3d} reactants")
    
    # Show unmatched reactants
    if unmatched:
        print(f"\n{'=' * 100}")
        print(f"UNMATCHED REACTANTS ({len(unmatched)})")
        print(f"{'=' * 100}")
        for smiles in unmatched[:20]:  # Show first 20
            print(f"  ❌ {smiles}")
        if len(unmatched) > 20:
            print(f"  ... and {len(unmatched)-20} more")
    
    # Show ambiguous cases
    if ambiguous:
        print(f"\n{'=' * 100}")
        print(f"AMBIGUOUS CASES ({len(ambiguous)})")
        print(f"{'=' * 100}")
        for smiles, matches in ambiguous[:10]:  # Show first 10
            print(f"\n  {smiles}")
            print(f"  Total matches: {len(matches)}")
            # Group by group
            by_group_amb = defaultdict(list)
            for m in matches:
                by_group_amb[m['group']].append(m)
            for group, group_matches in by_group_amb.items():
                print(f"    {group}:")
                for m in group_matches[:3]:
                    print(f"      - {m['member_type']:20s} ({m['name']})")
        if len(ambiguous) > 10:
            print(f"\n  ... and {len(ambiguous)-10} more ambiguous cases")
    
    # Test specific important cases
    print(f"\n{'=' * 100}")
    print("VALIDATION OF KEY SUBSTRATE TYPES")
    print(f"{'=' * 100}")
    
    test_cases = [
        ("Brc1ccccc1", "ArBr", "Simple aryl bromide"),
        ("Clc1ccncc1", "HetArCl or PyrimidineCl", "Heteroaryl chloride (pyridine)"),
        ("Ic1ccc(C(F)(F)F)cc1", "ArI", "Electron-poor aryl iodide"),
        ("Nc1ccccc1", "ArNH2", "Aniline"),
        ("NCC", "RNH2", "Primary aliphatic amine"),
        ("N1CCCC1", "pyrrolidine or cyclic amine", "Pyrrolidine"),
        ("N1CCOCC1", "morpholine", "Morpholine"),
        ("c1ccc(B(O)O)cc1", "ArB(OH)2", "Aryl boronic acid"),
        ("Oc1ccccc1", "ArOH", "Phenol"),
        ("OCC", "ROH-primary", "Primary alcohol"),
        ("CC(=O)O", "RCO2H", "Carboxylic acid"),
        ("CC(=O)Cl", "RCOCl", "Acyl chloride"),
        ("C=Cc1ccccc1", "Ar-alkene", "Styrene"),
        ("C#Cc1ccccc1", "Ar-alkyne", "Phenylacetylene"),
        ("c1ccccc1[Sn](C)(C)C", "organostannane", "Organostannane"),
        ("c1ccc([Zn]Cl)cc1", "Ar-M", "Organozinc"),
        ("c1ccc([Mg]Br)cc1", "Ar-M", "Grignard"),
    ]
    
    for smiles, expected_type, description in test_cases:
        result = classify_reactant(smiles, reactant_types)
        if result:
            status = "✅" if expected_type.lower() in result['member_type'].lower() or expected_type.lower() in result['category'].lower() else "⚠️ "
            print(f"{status} {description:30s}: {result['member_type']:20s} ({result['category']})")
        else:
            print(f"❌ {description:30s}: NO MATCH")
    
    # Summary statistics
    print(f"\n{'=' * 100}")
    print("SUMMARY")
    print(f"{'=' * 100}")
    print(f"Total unique reactants tested: {len(unique_reactants)}")
    print(f"Successfully classified: {matched} ({match_rate:.1f}%)")
    print(f"Unmatched: {len(unmatched)} ({100-match_rate:.1f}%)")
    print(f"Categories represented: {len(by_category)}")
    print(f"Functional groups represented: {len([g for g in by_group.keys() if g != 'Unmatched'])}")
    
    # Most common reactant types
    print(f"\nMost common reactant types:")
    type_counts = Counter([r['member_type'] for r in results.values() if r])
    for rtype, count in type_counts.most_common(10):
        print(f"  {rtype:20s}: {count:3d}")
    
    print(f"\n{'=' * 100}")
    print("✅ SMARTS pattern testing complete!")
    print(f"{'=' * 100}")

if __name__ == "__main__":
    main()
