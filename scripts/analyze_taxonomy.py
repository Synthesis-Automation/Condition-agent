#!/usr/bin/env python3
"""Quick analysis script for taxonomy organic groups."""

import json
from pathlib import Path

def main():
    # Load organic groups
    groups_file = Path('c:/Git-softwares/Condition-agent/chemtools/taxonomy/data/organic_groups.v1.3.json')
    with open(groups_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # Count groups by kind
    scaffolds = [g for g in data['groups'] if g.get('kind') == 'scaffold']
    substituents = [g for g in data['groups'] if g.get('kind') == 'substituent']

    print(f'Total groups: {len(data["groups"])}')
    print(f'  Scaffolds: {len(scaffolds)}')
    print(f'  Substituents: {len(substituents)}')
    print()

    # Check for missing priorities
    no_priority = [g['id'] for g in data['groups'] if 'priority' not in g]
    if no_priority:
        print(f'Groups without priority: {len(no_priority)}')
        for gid in no_priority[:10]:
            print(f'  - {gid}')
        if len(no_priority) > 10:
            print(f'  ... and {len(no_priority) - 10} more')
    else:
        print('All groups have priorities: OK')
    print()

    # Check for duplicate IDs
    ids = [g['id'] for g in data['groups']]
    duplicates = [id for id in set(ids) if ids.count(id) > 1]
    if duplicates:
        print(f'Duplicate IDs found: {duplicates}')
    else:
        print('No duplicate IDs: OK')
    print()

    # Check attachpoint mappings
    print('Attachpoint mapping validation:')
    issues = []
    for g in data['groups']:
        smarts = g.get('smarts', '')
        kind = g.get('kind', '')
        gid = g.get('id', '')
        
        has_1 = ':1]' in smarts or ':1>' in smarts
        has_2 = ':2]' in smarts or ':2>' in smarts
        
        if kind == 'scaffold' and not has_1:
            issues.append(f"  {gid:20s} (scaffold): Missing :1 mapping in SMARTS")
        elif kind == 'substituent' and not has_2:
            issues.append(f"  {gid:20s} (substituent): Missing :2 mapping in SMARTS")
    
    if issues:
        print(f'  Found {len(issues)} attachpoint issues:')
        for issue in issues[:20]:
            print(issue)
        if len(issues) > 20:
            print(f'  ... and {len(issues) - 20} more')
    else:
        print('  All attachpoint mappings correct: OK')
    print()

    # Check SMARTS syntax basics
    print('SMARTS syntax checks:')
    syntax_issues = []
    for g in data['groups']:
        smarts = g.get('smarts', '')
        gid = g.get('id', '')
        
        # Basic bracket matching
        if smarts.count('[') != smarts.count(']'):
            syntax_issues.append(f"  {gid}: Unmatched brackets")
        if smarts.count('(') != smarts.count(')'):
            syntax_issues.append(f"  {gid}: Unmatched parentheses")
    
    if syntax_issues:
        print(f'  Found {len(syntax_issues)} syntax issues:')
        for issue in syntax_issues:
            print(issue)
    else:
        print('  Basic syntax looks good: OK')
    print()

    # Priority distribution
    print('Priority distribution:')
    priority_counts = {}
    for g in data['groups']:
        p = g.get('priority', 1)
        priority_counts[p] = priority_counts.get(p, 0) + 1
    
    for p in sorted(priority_counts.keys()):
        print(f'  Priority {p}: {priority_counts[p]} groups')
    print()

    # Scaffold groups list
    print('Scaffold groups:')
    for g in scaffolds:
        print(f"  {g['id']:20s} priority={g.get('priority', '?'):2} {g.get('description', '')[:50]}")

if __name__ == '__main__':
    main()
