#!/usr/bin/env python3
"""
Test scdb_matcher for all Suzuki coupling reactions
Generates comprehensive summary of matches, conditions, and analysis
"""

import sys
import json
from pathlib import Path
from collections import defaultdict

# Add parent directory to path to enable imports
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

# Now import after path is set
from chemtools.scdb_matcher.matcher import match
from chemtools.scdb_matcher.loader import load_db

# Read sample reactions directly from file
def load_sample_reactions():
    """Load reactions from sample_reactions.py file"""
    sample_file = parent_dir / "tests" / "sample_reactions.py"
    with open(sample_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Execute to get SAMPLE_REACTIONS
    namespace = {}
    exec(content, namespace)
    return namespace.get('SAMPLE_REACTIONS', [])

print("="*80)
print(" SUZUKI COUPLING SCDB MATCHER TEST")
print("="*80)

# Load Suzuki condition database
print("\n[1] Loading Suzuki Condition Database...")
suzuki_db_path = Path("data/conditionDB/suzuki_db.json")

if not suzuki_db_path.exists():
    print(f"ERROR: Suzuki database not found at {suzuki_db_path}")
    sys.exit(1)

suzuki_db = load_db(str(suzuki_db_path))
print(f"    Loaded database: {suzuki_db.reaction_type if hasattr(suzuki_db, 'reaction_type') else 'Unknown'}")
print(f"    Database type: {type(suzuki_db).__name__}")
if hasattr(suzuki_db, 'entries'):
    print(f"    Number of entries: {len(suzuki_db.entries)}")

# Extract Suzuki reactions
print("\n[2] Extracting Suzuki Reactions from Sample Database...")
SAMPLE_REACTIONS = load_sample_reactions()
suzuki_reactions = []

for rxn in SAMPLE_REACTIONS:
    if isinstance(rxn, str) and '>>' in rxn and 'Suzuki' in rxn:
        # Parse reaction
        if ' (' in rxn:
            label_start = rxn.rfind(' (')
            smiles = rxn[:label_start].strip()
            label = rxn[label_start+2:].rstrip(')').strip()
        else:
            smiles = rxn
            label = "Unknown"
        
        suzuki_reactions.append({
            'smiles': smiles,
            'label': label,
            'raw': rxn
        })

print(f"    Found {len(suzuki_reactions)} Suzuki reactions")

# Test each reaction
print("\n[3] Testing Reactions Against Condition Database...")
print("-"*80)

results = []
match_stats = defaultdict(int)
condition_usage = defaultdict(int)

for i, rxn in enumerate(suzuki_reactions, 1):
    print(f"\n[{i}/{len(suzuki_reactions)}] {rxn['label']}")
    print(f"    SMILES: {rxn['smiles'][:70]}...")
    
    try:
        # Match reaction using the scdb_matcher
        result = match(suzuki_db, rxn['smiles'])
        
        if result:
            print(f"    MATCH FOUND:")
            match_stats['matched'] += 1
            
            # Extract match details
            print(f"      - Match Type: {result.match_type}")
            print(f"      - Entry ID: {result.entry_id}")
            print(f"      - Entry Name: {result.entry_name}")
            print(f"      - Priority: {result.priority}")
            
            # Extract conditions
            conditions = result.conditions
            if conditions:
                print(f"      - Catalyst: {conditions.get('catalyst', 'N/A')}")
                print(f"      - Base: {conditions.get('base', 'N/A')}")
                print(f"      - Solvent: {conditions.get('solvent', 'N/A')}")
                print(f"      - Temperature: {conditions.get('temperature', 'N/A')}")
                
                # Track condition usage
                condition_usage[conditions.get('catalyst', 'unknown')] += 1
            
            # Extract electrophile/nucleophile info from trace if available
            trace = result.trace
            electrophile_type = 'matched'
            nucleophile_type = 'matched'
            
            results.append({
                'reaction': rxn['label'],
                'smiles': rxn['smiles'],
                'matched': True,
                'match_type': result.match_type,
                'entry_id': result.entry_id,
                'entry_name': result.entry_name,
                'priority': result.priority,
                'electrophile': electrophile_type,
                'nucleophile': nucleophile_type,
                'conditions': conditions
            })
        else:
            print(f"    NO MATCH: No conditions found")
            match_stats['no_match'] += 1
            
            results.append({
                'reaction': rxn['label'],
                'smiles': rxn['smiles'],
                'matched': False,
                'match_type': 'none',
                'entry_id': None,
                'entry_name': None,
                'priority': 0,
                'electrophile': 'N/A',
                'nucleophile': 'N/A',
                'conditions': {}
            })
    
    except Exception as e:
        print(f"    ERROR: {str(e)}")
        match_stats['error'] += 1
        
        results.append({
            'reaction': rxn['label'],
            'smiles': rxn['smiles'],
            'matched': False,
            'match_type': 'error',
            'entry_id': None,
            'entry_name': None,
            'priority': 0,
            'error': str(e)
        })

# Generate Summary
print("\n" + "="*80)
print(" SUMMARY STATISTICS")
print("="*80)

print(f"\n[A] MATCHING PERFORMANCE")
print("-"*80)
print(f"  Total Suzuki reactions tested:  {len(suzuki_reactions)}")
print(f"  Successfully matched:           {match_stats['matched']} ({match_stats['matched']/len(suzuki_reactions)*100:.1f}%)")
print(f"  No matches found:               {match_stats['no_match']} ({match_stats['no_match']/len(suzuki_reactions)*100:.1f}%)")
print(f"  Errors encountered:             {match_stats['error']} ({match_stats['error']/len(suzuki_reactions)*100:.1f}%)")

# Analyze by electrophile type
print(f"\n[B] ELECTROPHILE TYPE DISTRIBUTION")
print("-"*80)
electrophile_types = defaultdict(int)
for r in results:
    if r.get('matched'):
        electrophile_types[r['electrophile']] += 1

if electrophile_types:
    for etype, count in sorted(electrophile_types.items(), key=lambda x: x[1], reverse=True):
        bar = 'â–? * min(count, 40)
        print(f"  {etype:30s} {count:3d} {bar}")
else:
    print("  No electrophile data available")

# Analyze by nucleophile type
print(f"\n[C] NUCLEOPHILE TYPE DISTRIBUTION")
print("-"*80)
nucleophile_types = defaultdict(int)
for r in results:
    if r.get('matched'):
        nucleophile_types[r['nucleophile']] += 1

if nucleophile_types:
    for ntype, count in sorted(nucleophile_types.items(), key=lambda x: x[1], reverse=True):
        bar = 'â–? * min(count, 40)
        print(f"  {ntype:30s} {count:3d} {bar}")
else:
    print("  No nucleophile data available")

# Analyze catalyst usage
print(f"\n[D] CATALYST USAGE")
print("-"*80)
if condition_usage:
    for catalyst, count in sorted(condition_usage.items(), key=lambda x: x[1], reverse=True)[:10]:
        bar = 'â–? * min(count, 40)
        print(f"  {catalyst:40s} {count:3d} {bar}")
else:
    print("  No catalyst data available")

# Score distribution
print(f"\n[E] MATCH PRIORITY DISTRIBUTION")
print("-"*80)
matched_results = [r for r in results if r.get('matched')]
if matched_results:
    priorities = [r['priority'] for r in matched_results if isinstance(r.get('priority'), (int, float))]
    if priorities:
        print(f"  Average priority:  {sum(priorities)/len(priorities):.2f}")
        print(f"  Max priority:      {max(priorities)}")
        print(f"  Min priority:      {min(priorities)}")
        
        # Priority ranges
        ranges = {
            '90-100 (Highest)': 0,
            '70-89 (High)': 0,
            '50-69 (Medium)': 0,
            '0-49 (Low)': 0
        }
        
        for priority in priorities:
            if priority >= 90:
                ranges['90-100 (Highest)'] += 1
            elif priority >= 70:
                ranges['70-89 (High)'] += 1
            elif priority >= 50:
                ranges['50-69 (Medium)'] += 1
            else:
                ranges['0-49 (Low)'] += 1
        
        print(f"\n  Priority ranges:")
        for range_name, count in ranges.items():
            bar = 'â–? * min(count, 40)
            print(f"    {range_name:25s} {count:3d} {bar}")
else:
    print("  No priority data available")

# Detailed results for unmatched
print(f"\n[F] UNMATCHED REACTIONS ANALYSIS")
print("-"*80)
unmatched = [r for r in results if not r.get('matched')]
if unmatched:
    print(f"  Total unmatched: {len(unmatched)}")
    print(f"\n  Details:")
    for i, r in enumerate(unmatched[:10], 1):
        print(f"    {i}. {r['reaction'][:70]}")
        if 'error' in r:
            print(f"       Error: {r['error']}")
else:
    print("  All reactions successfully matched!")

# Export results to JSON
output_file = Path("scripts/suzuki_matcher_results.json")
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump({
        'summary': {
            'total_reactions': len(suzuki_reactions),
            'matched': match_stats['matched'],
            'no_match': match_stats['no_match'],
            'errors': match_stats['error'],
            'match_rate': f"{match_stats['matched']/len(suzuki_reactions)*100:.1f}%"
        },
        'electrophile_types': dict(electrophile_types),
        'nucleophile_types': dict(nucleophile_types),
        'catalyst_usage': dict(condition_usage),
        'detailed_results': results
    }, f, indent=2, ensure_ascii=False)

print(f"\n[G] OUTPUT FILES")
print("-"*80)
print(f"  Detailed results saved to: {output_file}")

print("\n" + "="*80)
print(" TEST COMPLETE")
print("="*80)
