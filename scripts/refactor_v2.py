"""
REVISED Strategy:
1. Reactant features (ArBr_reactant, etc.) are MORE SPECIFIC than present features (aryl_halide_present)
2. KEEP reactant features but REMOVE _reactant suffix → becomes ArBr_present (or just drop _present too?)
3. REMOVE general _present features that are superseded by specific ones
4. ADD reactant_metadata to kept features
5. REMOVE specific compound reagents (H2, mCPBA, OH-, etc.)
"""

import json
import shutil
from pathlib import Path

# Specific compounds to remove (reagents, not functional groups)
REMOVE_SPECIFIC_COMPOUNDS = {
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
    json_path = Path('chemtools/featurizers/calculable_features.json')
    backup_path = json_path.with_suffix('.json.v4-backup')
    
    # Restore from backup
    print(f"Restoring from: {backup_path}")
    shutil.copy2(backup_path, json_path)
    
    # Load
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    original_count = len(data['features'])
    print(f"Original features: {original_count}")
    
    # Process features
    new_features = []
    stats = {
        'kept_reactant': 0,
        'removed_specific': 0,
        'kept_present': 0,
        'kept_other': 0,
    }
    
    for feat in data['features']:
        token = feat['token']
        
        # Handle _reactant features
        if token.endswith('_reactant'):
            if token in REMOVE_SPECIFIC_COMPOUNDS:
                stats['removed_specific'] += 1
                print(f"  ✗ Removed specific: {token}")
                continue
            
            # Keep but rename by removing _reactant suffix and adding metadata
            new_token = token[:-9] + '_present'  # ArBr_reactant → ArBr_present
            feat['token'] = new_token
            
            # Move metadata to reactant_metadata
            if 'metadata' in feat:
                metadata = feat.pop('metadata')
                feat['reactant_metadata'] = {
                    'is_reactant': True,
                    'compound_type': metadata.get('reactant_name', ''),
                    'reactant_category': metadata.get('reactant_category', ''),
                    'reactant_member': metadata.get('reactant_member', ''),
                    'coupling_role': metadata.get('coupling_role', ''),
                }
                # Clean empty fields
                feat['reactant_metadata'] = {k: v for k, v in feat['reactant_metadata'].items() if v}
            
            new_features.append(feat)
            stats['kept_reactant'] += 1
            if stats['kept_reactant'] <= 10:
                print(f"  ✓ Kept & renamed: {token} → {new_token}")
        
        # Handle _present features
        elif token.endswith('_present'):
            stats['kept_present'] += 1
            new_features.append(feat)
        
        # Handle other features
        else:
            stats['kept_other'] += 1
            new_features.append(feat)
    
    # Update
    data['features'] = new_features
    data['version'] = 'v5.0-unified-metadata'
    
    # Save
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print("\n" + "="*70)
    print("REFACTORING COMPLETE")
    print("="*70)
    print(f"Original features:              {original_count}")
    print(f"  Reactant (kept & renamed):    {stats['kept_reactant']}")
    print(f"  Reactant (removed specific):  {stats['removed_specific']}")
    print(f"  Present (kept):               {stats['kept_present']}")
    print(f"  Other (kept):                 {stats['kept_other']}")
    print(f"Final feature count:            {len(new_features)}")
    print(f"Reduction:                      {original_count - len(new_features)} ({(original_count - len(new_features))/original_count*100:.1f}%)")

if __name__ == '__main__':
    main()
