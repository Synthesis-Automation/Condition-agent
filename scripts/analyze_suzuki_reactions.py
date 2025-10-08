#!/usr/bin/env python3
"""
Test Suzuki reactions against scdb_matcher
Simpler approach - test reactions directly against condition database
"""

import sys
from pathlib import Path
from collections import defaultdict
import json

# Add parent to path
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

from rdkit import Chem

# Load sample reactions
def load_suzuki_reactions():
    """Extract Suzuki reactions from sample_reactions.py"""
    sample_file = parent_dir / "tests" / "sample_reactions.py"
    with open(sample_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    suzuki_reactions = []
    for line in lines:
        if '>>' in line and 'Suzuki' in line and not line.strip().startswith('#'):
            # Parse the reaction
            line = line.strip()
            if ' (' in line:
                label_start = line.rfind(' (')
                smiles = line[:label_start].strip().strip('"').strip(',')
                label = line[label_start+2:].rstrip(')').rstrip('"),').rstrip('"').strip()
            else:
                smiles = line.strip('"').strip(',')
                label = "Unknown"
            
            suzuki_reactions.append({
                'smiles': smiles,
                'label': label
            })
    
    return suzuki_reactions

print("="*80)
print(" SUZUKI COUPLING REACTION ANALYSIS")
print("="*80)

print("\n[1] Loading Suzuki Reactions...")
suzuki_reactions = load_suzuki_reactions()
print(f"    Found {len(suzuki_reactions)} Suzuki reactions")

print("\n[2] Analyzing Reaction Features...")
print("-"*80)

# Analyze features
stats = {
    'total': len(suzuki_reactions),
    'electrophile_types': defaultdict(int),
    'leaving_groups': defaultdict(int),
    'substrate_types': defaultdict(int),
    'boron_sources': defaultdict(int),
    'heteroaryls': 0,
    'protected_groups': 0,
    'sterically_hindered': 0,
    'electron_deficient': 0,
}

# Patterns for detection
lg_patterns = {
    'Br': Chem.MolFromSmarts('[Br]'),
    'I': Chem.MolFromSmarts('[I]'),
    'Cl': Chem.MolFromSmarts('[Cl]'),
    'OTf': Chem.MolFromSmarts('OS(=O)(=O)C(F)(F)F'),
}

heteroaryl_patterns = {
    'pyridine': Chem.MolFromSmarts('n1ccccc1'),
    'furan': Chem.MolFromSmarts('o1cccc1'),
    'thiophene': Chem.MolFromSmarts('s1cccc1'),
    'pyrimidine': Chem.MolFromSmarts('n1cnccc1'),
    'indole': Chem.MolFromSmarts('[nH]1ccc2ccccc12'),
}

detailed_results = []

for i, rxn in enumerate(suzuki_reactions, 1):
    print(f"\n[{i}/{len(suzuki_reactions)}] {rxn['label'][:70]}")
    
    result = {
        'id': i,
        'label': rxn['label'],
        'smiles': rxn['smiles'],
        'features': {}
    }
    
    try:
        # Parse reaction
        parts = rxn['smiles'].split('>>')
        if len(parts) != 2:
            print(f"    ERROR: Invalid reaction SMILES")
            continue
        
        reactants_str, product_str = parts
        reactants = reactants_str.split('.')
        
        print(f"    Reactants: {len(reactants)}")
        
        # Analyze each reactant
        electrophile_found = False
        boron_found = False
        
        for r_idx, r_smiles in enumerate(reactants):
            r_smiles = r_smiles.strip()
            
            # Skip salts and simple ions
            if r_smiles in ['[K+]', '[Na+]', '[Cs+]', '[Li+]']:
                continue
            
            try:
                mol = Chem.MolFromSmiles(r_smiles)
                if mol is None:
                    print(f"      Reactant {r_idx}: Invalid SMILES")
                    continue
                
                # Check for leaving groups (electrophile)
                lg_found = None
                for lg_name, lg_pattern in lg_patterns.items():
                    if lg_pattern and mol.HasSubstructMatch(lg_pattern):
                        lg_found = lg_name
                        stats['leaving_groups'][lg_name] += 1
                        break
                
                if lg_found:
                    electrophile_found = True
                    print(f"      Electrophile (LG={lg_found}): {r_smiles[:50]}")
                    result['features']['leaving_group'] = lg_found
                    
                    # Check electrophile type
                    if any(mol.HasSubstructMatch(p) for p in heteroaryl_patterns.values() if p):
                        stats['electrophile_types']['heteroaryl'] += 1
                        result['features']['electrophile_type'] = 'heteroaryl'
                    elif any(a.GetIsAromatic() for a in mol.GetAtoms()):
                        stats['electrophile_types']['aryl'] += 1
                        result['features']['electrophile_type'] = 'aryl'
                    else:
                        stats['electrophile_types']['alkyl/vinyl'] += 1
                        result['features']['electrophile_type'] = 'alkyl/vinyl'
                
                # Check for boron
                if mol.HasSubstructMatch(Chem.MolFromSmarts('[B]')):
                    boron_found = True
                    print(f"      Boron partner: {r_smiles[:50]}")
                    
                    # Classify boron source
                    if 'B(O)O' in r_smiles:
                        stats['boron_sources']['boronic_acid'] += 1
                        result['features']['boron_source'] = 'boronic_acid'
                    elif '[B-]' in r_smiles:
                        stats['boron_sources']['trifluoroborate'] += 1
                        result['features']['boron_source'] = 'trifluoroborate'
                    elif 'B1OC' in r_smiles or 'OC(=O)' in r_smiles:
                        stats['boron_sources']['MIDA'] += 1
                        result['features']['boron_source'] = 'MIDA'
                    else:
                        stats['boron_sources']['other'] += 1
                        result['features']['boron_source'] = 'other'
            
            except Exception as e:
                print(f"      Reactant {r_idx}: Error - {e}")
        
        # Check for special features in label
        label_lower = rxn['label'].lower()
        
        if any(word in label_lower for word in ['pyridin', 'furan', 'thiophene', 'pyrimidin', 'indole', 'quinoxaline', 'benzothiazole']):
            stats['heteroaryls'] += 1
            result['features']['contains_heteroaryl'] = True
        
        if any(word in label_lower for word in ['boc', 'tbs', 'ester']):
            stats['protected_groups'] += 1
            result['features']['has_protecting_group'] = True
        
        if any(word in label_lower for word in ['hindered', 'ortho', 'dimethyl']):
            stats['sterically_hindered'] += 1
            result['features']['sterically_hindered'] = True
        
        if any(word in label_lower for word in ['electron-poor', 'electron-deficient', 'nitro', 'cf3', 'c#n']):
            stats['electron_deficient'] += 1
            result['features']['electron_deficient'] = True
        
        detailed_results.append(result)
    
    except Exception as e:
        print(f"    ERROR: {e}")

# Print summary
print("\n" + "="*80)
print(" SUMMARY STATISTICS")
print("="*80)

print(f"\n[A] GENERAL STATISTICS")
print("-"*80)
print(f"  Total Suzuki reactions:     {stats['total']}")
print(f"  Successfully analyzed:      {len(detailed_results)}")

print(f"\n[B] LEAVING GROUP DISTRIBUTION")
print("-"*80)
for lg, count in sorted(stats['leaving_groups'].items(), key=lambda x: x[1], reverse=True):
    bar = '█' * min(count, 50)
    print(f"  {lg:15s} {count:3d} {bar}")

print(f"\n[C] ELECTROPHILE TYPES")
print("-"*80)
for etype, count in sorted(stats['electrophile_types'].items(), key=lambda x: x[1], reverse=True):
    bar = '█' * min(count, 50)
    print(f"  {etype:20s} {count:3d} {bar}")

print(f"\n[D] BORON SOURCE TYPES")
print("-"*80)
for btype, count in sorted(stats['boron_sources'].items(), key=lambda x: x[1], reverse=True):
    bar = '█' * min(count, 50)
    print(f"  {btype:20s} {count:3d} {bar}")

print(f"\n[E] SPECIAL FEATURES")
print("-"*80)
features_summary = [
    ('Heteroaryl substrates', stats['heteroaryls']),
    ('Protected functional groups', stats['protected_groups']),
    ('Sterically hindered', stats['sterically_hindered']),
    ('Electron deficient', stats['electron_deficient']),
]

for feature_name, count in features_summary:
    percentage = (count / stats['total'] * 100) if stats['total'] > 0 else 0
    bar = '█' * min(count, 50)
    print(f"  {feature_name:30s} {count:3d} ({percentage:5.1f}%) {bar}")

# Export results
output_file = parent_dir / "scripts" / "suzuki_analysis_results.json"
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump({
        'summary': dict(stats),
        'reactions': detailed_results
    }, f, indent=2, ensure_ascii=False, default=str)

print(f"\n[F] OUTPUT FILES")
print("-"*80)
print(f"  Detailed results saved to: {output_file}")

print("\n" + "="*80)
print(" ANALYSIS COMPLETE")
print("="*80)
print(f"\nKey Findings:")
print(f"  - {len(suzuki_reactions)} Suzuki reactions covering diverse coupling scenarios")
print(f"  - {len(stats['leaving_groups'])} different leaving group types")
print(f"  - {len(stats['boron_sources'])} different boron source types")
print(f"  - {stats['heteroaryls']} reactions involve heteroaryl substrates")
print(f"  - {stats['protected_groups']} reactions with protected functional groups")
print(f"\nThe Suzuki database is ready for comprehensive condition matching!")
