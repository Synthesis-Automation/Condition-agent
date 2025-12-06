"""
Streaming refactor of calculable_features.json:
1. Build mapping of _reactant → _present features
2. Stream through file, merging metadata
3. Skip specific compound entries
"""

import json
import shutil
from pathlib import Path
from collections import defaultdict

# Specific compounds to remove
REMOVE_TOKENS = {
    "NaBH4_reactant", "LiBH4_reactant", "LiAlH4_reactant", "NaBH(OAc)3_reactant", 
    "DIBAL_reactant", "BH3_reactant", "B2H6_reactant", "9-BBN_reactant", 
    "catecholborane_reactant", "pinacolborane_reactant",
    "RMgBr_reactant", "RMgCl_reactant", "RMgI_reactant", 
    "RZnBr_reactant", "RZnCl_reactant", "R2Zn_reactant",
    "RLi_reactant", "ArLi_reactant", "nBuLi_reactant", "tBuLi_reactant",
    "mCPBA_reactant", "PCC_reactant", "DMP_reactant", "IBX_reactant", 
    "TEMPO_reactant", "CrO3_reactant", "KMnO4_reactant", "DDQ_reactant", "NaOCl_reactant",
    "H2O2_reactant", "tBuOOH_reactant", "H2_reactant", "formic_acid_reactant", "isopropanol_reactant",
    "NaCN_reactant", "KCN_reactant", "CN-_reactant",
    "NaN3_reactant", "N3-_reactant",
    "NaI_reactant", "KI_reactant", "I-_reactant",
    "NaOH_reactant", "KOH_reactant", "OH-_reactant",
    "NaOMe_reactant", "NaOEt_reactant", "KOtBu_reactant", "RO-_reactant",
    "malonate_reactant", "acetoacetate_reactant",
}

def main():
    json_path = Path('chemtools/taxonomy/data/calculable_features.json')
    
    # Create backup
    backup_path = json_path.with_suffix('.json.v4-backup')
    if backup_path.exists():
        print(f"Backup already exists: {backup_path}")
    else:
        shutil.copy2(json_path, backup_path)
        print(f"✓ Created backup: {backup_path}")
    
    # Load data
    print("Loading calculable_features.json...")
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    original_count = len(data['features'])
    print(f"Original features: {original_count}")
    
    # Build index: base_token → {present_feature, reactant_feature}
    token_map = defaultdict(lambda: {'present': None, 'reactant': []})
    
    for feat in data['features']:
        token = feat['token']
        if token.endswith('_reactant'):
            base = token[:-9]
            token_map[base]['reactant'].append(feat)
        else:
            # Try to extract base token
            if token.endswith('_present'):
                base = token[:-8]
                token_map[base]['present'] = feat
    
    # Process merges
    merged_count = 0
    skipped_count = 0
    removed_count = 0
    
    new_features = []
    
    for feat in data['features']:
        token = feat['token']
        
        # Skip reactant features (will be merged or removed)
        if token.endswith('_reactant'):
            if token in REMOVE_TOKENS:
                removed_count += 1
                print(f"  ✗ Removed: {token}")
            else:
                skipped_count += 1
            continue
        
        # Keep non-reactant feature, but check if we need to add metadata
        base_token = token[:-8] if token.endswith('_present') else None
        
        if base_token and base_token in token_map and token_map[base_token]['reactant']:
            # Merge metadata from first reactant feature
            reactant_feat = token_map[base_token]['reactant'][0]
            metadata = reactant_feat.get('metadata', {})
            
            if metadata:
                feat['reactant_metadata'] = {
                    'is_reactant': True,
                    'compound_type': metadata.get('reactant_name', ''),
                    'reactant_category': metadata.get('reactant_category', ''),
                    'reactant_member': metadata.get('reactant_member', ''),
                    'coupling_role': metadata.get('coupling_role', ''),
                }
                # Clean empty fields
                feat['reactant_metadata'] = {k: v for k, v in feat['reactant_metadata'].items() if v}
                merged_count += 1
                print(f"  ✓ Merged metadata: {token}")
        
        new_features.append(feat)
    
    # Update data
    data['features'] = new_features
    data['version'] = 'v5.0-unified-metadata'
    
    # Save
    print("\nSaving refactored version...")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print("\n" + "="*70)
    print("REFACTORING COMPLETE")
    print("="*70)
    print(f"Original features:           {original_count}")
    print(f"Merged to _present:          {merged_count}")
    print(f"Removed (specific cpds):     {removed_count}")
    print(f"Removed (_reactant suffix):  {skipped_count}")
    print(f"Final feature count:         {len(new_features)}")
    print(f"Reduction:                   {original_count - len(new_features)} ({(original_count - len(new_features))/original_count*100:.1f}%)")
    print(f"\nBackup: {backup_path}")

if __name__ == '__main__':
    main()
