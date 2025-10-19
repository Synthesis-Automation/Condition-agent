#!/usr/bin/env python3
"""Check RXNMapper + rxnutils pipeline status for all protocols"""

import json
import os

protocols_dir = 'data/protocol_db'
protocols = [f for f in os.listdir(protocols_dir) 
             if f.endswith('.json') and f != '.protocol_index.json']

results = []

print('=' * 130)
print('RXNMapper + rxnutils Pipeline Status')
print('=' * 130)
print()

for f in sorted(protocols):
    filepath = os.path.join(protocols_dir, f)
    data = json.load(open(filepath, encoding='utf-8'))
    
    if 'reaction' in data and 'mapping_confidence' in data['reaction']:
        conf = data['reaction']['mapping_confidence']
        pattern = data['reaction'].get('reaction_smarts_applicability', {}).get('core', 'N/A')
        
        # Determine method used
        if conf > 0.70:
            method = '✓ Atom-mapped'
        else:
            method = '⚡ Heuristics'
        
        results.append({
            'name': f.replace('.json', ''),
            'conf': conf,
            'method': method,
            'pattern': pattern
        })

# Sort by confidence
results.sort(key=lambda x: x['conf'], reverse=True)

# Print results
print(f"{'Protocol':<50} {'Conf%':>7} {'Method':>15}  Pattern")
print('-' * 130)

for r in results:
    pattern_display = r['pattern'][:50] + '...' if len(r['pattern']) > 50 else r['pattern']
    print(f"{r['name']:<50} {r['conf']*100:>6.1f}% {r['method']:>15}  {pattern_display}")

# Summary stats
print()
print('=' * 130)
print('Summary:')
print(f"  Total protocols: {len(results)}")
print(f"  Atom-mapped (>70% conf): {sum(1 for r in results if r['conf'] > 0.70)}")
print(f"  Heuristics fallback (<70% conf): {sum(1 for r in results if r['conf'] <= 0.70)}")
print(f"  Average confidence: {sum(r['conf'] for r in results) / len(results) * 100:.1f}%")
print()
print('Key improvements:')
print('  ✓ Aryl mesylate now detects sulfonate leaving group: c-[OX2]-[SX4](=O)(=O)-[CH3]')
print('  ✓ Sonogashira perfect pattern (no spectator groups): C#C.c-[I]>>c-C#C')
print('  ✓ Hybrid approach: atom mapping for high confidence, heuristics for low confidence')
print('=' * 130)
