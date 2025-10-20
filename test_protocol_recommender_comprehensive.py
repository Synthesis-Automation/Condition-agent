"""
Comprehensive Protocol Recommender Test
Tests the protocol recommender with all sample reactions from sample_reactions.py
"""

import sys
import json
from pathlib import Path
from collections import defaultdict
import time

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / 'tests'))

from sample_reactions import SAMPLE_REACTIONS
from chemtools.protocol.recommend import ProtocolRecommender


def parse_reaction_smiles(reaction_str):
    """Extract SMILES from reaction string (format: 'SMILES (Description)')"""
    if '>>' not in reaction_str:
        return None
    
    # The format is: "REACTANTS>>PRODUCTS (Description)"
    # But PRODUCTS can contain parentheses in SMILES!
    # Description is always at the end after a space and opening paren
    
    # Strategy: Find the last occurrence of " (" which starts the description
    desc_start = reaction_str.rfind(' (')
    if desc_start > 0 and reaction_str.endswith(')'):
        # Extract SMILES part (everything before " (Description)")
        smiles = reaction_str[:desc_start].strip()
    else:
        # No description, whole string is SMILES
        smiles = reaction_str.strip()
    
    # Validate it has >>
    if '>>' not in smiles:
        return None
    
    return smiles


def categorize_reaction(reaction_str):
    """Categorize reaction based on its description"""
    lower = reaction_str.lower()
    
    if 'suzuki' in lower:
        return 'Suzuki-Miyaura'
    elif 'buchwald' in lower or 'b-h:' in lower or 'c-n -' in lower:
        return 'Buchwald-Hartwig (C-N)'
    elif 'sonogashira' in lower:
        return 'Sonogashira'
    elif 'heck' in lower:
        return 'Heck'
    elif 'stille' in lower:
        return 'Stille'
    elif 'negishi' in lower:
        return 'Negishi'
    elif 'kumada' in lower:
        return 'Kumada'
    elif 'chan-lam' in lower:
        return 'Chan-Lam'
    elif 'ullmann' in lower:
        return 'Ullmann'
    elif 'amide' in lower and 'amidation' in lower or 'amide:' in lower:
        return 'Amide Formation'
    elif 'hydrogenation' in lower:
        return 'Hydrogenation'
    elif 'reduction' in lower or 'nabh4' in lower or 'lialh4' in lower:
        return 'Reduction'
    elif 'oxidation' in lower:
        return 'Oxidation'
    elif 'click' in lower or 'cuaac' in lower:
        return 'Click Chemistry'
    elif 'diels-alder' in lower:
        return 'Diels-Alder'
    elif 'wittig' in lower:
        return 'Wittig'
    elif 'grignard' in lower:
        return 'Grignard'
    elif 'aldol' in lower:
        return 'Aldol'
    elif 'sn2' in lower or 'sn1' in lower:
        return 'Substitution'
    elif 'e2' in lower or 'e1' in lower:
        return 'Elimination'
    else:
        return 'Other'


def test_protocol_recommender():
    """Test protocol recommender with all sample reactions"""
    
    print("=" * 80)
    print("COMPREHENSIVE PROTOCOL RECOMMENDER TEST")
    print("=" * 80)
    print()
    
    # Initialize recommender
    print("Initializing protocol recommender...")
    try:
        recommender = ProtocolRecommender()
        protocol_count = len(recommender.indexer.records)
        print(f"[OK] Loaded {protocol_count} protocols")
        print()
    except Exception as e:
        print(f"[ERROR] Failed to load recommender: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Statistics
    stats = {
        'total': 0,
        'skipped': 0,
        'parsed': 0,
        'with_matches': 0,
        'without_matches': 0,
        'errors': 0,
        'by_category': defaultdict(lambda: {'total': 0, 'matched': 0, 'no_match': 0})
    }
    
    results = []
    
    print("Testing reactions...")
    print("-" * 80)
    
    for i, reaction_str in enumerate(SAMPLE_REACTIONS[1:], 1):  # Skip first "Select..." entry
        stats['total'] += 1
        
        # Parse SMILES
        smiles = parse_reaction_smiles(reaction_str)
        if not smiles:
            stats['skipped'] += 1
            continue
        
        stats['parsed'] += 1
        category = categorize_reaction(reaction_str)
        stats['by_category'][category]['total'] += 1
        
        # Get short description
        if '(' in reaction_str and ')' in reaction_str:
            desc = reaction_str[reaction_str.find('(')+1:reaction_str.rfind(')')]
        else:
            desc = reaction_str[:50] + "..."
        
        # Test recommendation
        try:
            start_time = time.time()
            response = recommender.recommend(
                reaction_smiles=smiles,
                k=3,
                use_smarts_filter=True
            )
            elapsed_ms = (time.time() - start_time) * 1000
            
            # Standard format uses 'recommended_conditions' key
            recommendations = response.get('recommended_conditions', [])
            num_protocols = len(recommendations)
            
            if num_protocols > 0:
                stats['with_matches'] += 1
                stats['by_category'][category]['matched'] += 1
                
                # Get top match info
                top_match = recommendations[0]
                metadata = top_match.get('protocol_metadata', {})
                top_family = metadata.get('reaction_family', 'Unknown')
                top_similarity = top_match.get('confidence', 0.0)
                
                result = {
                    'index': i,
                    'category': category,
                    'description': desc,
                    'smiles': smiles,
                    'num_matches': num_protocols,
                    'top_family': top_family,
                    'top_similarity': top_similarity,
                    'time_ms': elapsed_ms,
                    'status': 'matched'
                }
                
                # Print result
                print(f"{i:3d}. [{category:25s}] {desc[:40]:40s}")
                print(f"     [MATCH] {num_protocols} match(es) | Top: {top_family} (sim={top_similarity:.3f}) | {elapsed_ms:.1f}ms")
                
            else:
                stats['without_matches'] += 1
                stats['by_category'][category]['no_match'] += 1
                
                result = {
                    'index': i,
                    'category': category,
                    'description': desc,
                    'smiles': smiles,
                    'num_matches': 0,
                    'time_ms': elapsed_ms,
                    'status': 'no_match'
                }
                
                # Print result
                print(f"{i:3d}. [{category:25s}] {desc[:40]:40s}")
                print(f"     [NO MATCH] {elapsed_ms:.1f}ms")
                
                # Check if warning was issued
                if response.get('extras', {}).get('smarts_filter_warning'):
                    print(f"     [INFO] {response['extras']['smarts_filter_warning']}")
            
            results.append(result)
            
        except Exception as e:
            stats['errors'] += 1
            stats['by_category'][category]['matched'] += 0  # Count as attempt
            
            result = {
                'index': i,
                'category': category,
                'description': desc,
                'smiles': smiles,
                'error': str(e),
                'status': 'error'
            }
            results.append(result)
            
            print(f"{i:3d}. [{category:25s}] {desc[:40]:40s}")
            print(f"     [ERROR] {str(e)[:50]}")
        
        print()
    
    # Print summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    
    print(f"Total reactions:           {stats['total']}")
    print(f"Skipped (parse error):     {stats['skipped']}")
    print(f"Successfully parsed:       {stats['parsed']}")
    print(f"With protocol matches:     {stats['with_matches']} ({stats['with_matches']/stats['parsed']*100:.1f}%)")
    print(f"Without protocol matches:  {stats['without_matches']} ({stats['without_matches']/stats['parsed']*100:.1f}%)")
    print(f"Errors:                    {stats['errors']}")
    print()
    
    # Category breakdown
    print("=" * 80)
    print("BY CATEGORY")
    print("=" * 80)
    print()
    print(f"{'Category':<30} {'Total':>8} {'Matched':>8} {'No Match':>8} {'Match %':>10}")
    print("-" * 80)
    
    for category in sorted(stats['by_category'].keys()):
        cat_stats = stats['by_category'][category]
        total = cat_stats['total']
        matched = cat_stats['matched']
        no_match = cat_stats['no_match']
        match_pct = (matched / total * 100) if total > 0 else 0
        
        print(f"{category:<30} {total:>8} {matched:>8} {no_match:>8} {match_pct:>9.1f}%")
    
    print()
    print("=" * 80)
    
    # Save detailed results to JSON
    output_file = Path(__file__).parent / 'protocol_recommender_test_results.json'
    with open(output_file, 'w') as f:
        json.dump({
            'stats': {k: dict(v) if isinstance(v, defaultdict) else v 
                     for k, v in stats.items()},
            'results': results
        }, f, indent=2)
    
    print(f"\nDetailed results saved to: {output_file}")
    print()


if __name__ == '__main__':
    test_protocol_recommender()
