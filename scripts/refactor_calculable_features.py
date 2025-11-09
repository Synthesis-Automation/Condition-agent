"""
Refactor calculable_features.json to eliminate redundancy:
1. Remove _reactant suffix duplicates by merging metadata into _present features
2. Remove specific compound entries (reagents like H2, mCPBA, OH-, etc.)
"""

import json
import shutil
from pathlib import Path
from typing import Dict, List, Set

# Specific compounds to remove (reagents, not structural features)
SPECIFIC_COMPOUNDS_TO_REMOVE = {
    # Reducing agents
    "NaBH4_reactant", "LiBH4_reactant", "LiAlH4_reactant", "NaBH(OAc)3_reactant", 
    "DIBAL_reactant", "BH3_reactant", "B2H6_reactant", "9-BBN_reactant", 
    "catecholborane_reactant", "pinacolborane_reactant",
    
    # Organometallics (specific compounds, not functional groups)
    "RMgBr_reactant", "RMgCl_reactant", "RMgI_reactant", 
    "RZnBr_reactant", "RZnCl_reactant", "R2Zn_reactant",
    "RLi_reactant", "ArLi_reactant", "nBuLi_reactant", "tBuLi_reactant",
    
    # Oxidizing agents
    "mCPBA_reactant", "PCC_reactant", "DMP_reactant", "IBX_reactant", 
    "TEMPO_reactant", "CrO3_reactant", "KMnO4_reactant", "DDQ_reactant", "NaOCl_reactant",
    "H2O2_reactant", "tBuOOH_reactant",
    
    # Specific inorganic reagents
    "H2_reactant", "formic_acid_reactant", "isopropanol_reactant",
    
    # Salts and nucleophiles
    "NaCN_reactant", "KCN_reactant", "CN-_reactant",
    "NaN3_reactant", "N3-_reactant",
    "NaI_reactant", "KI_reactant", "I-_reactant",
    "NaOH_reactant", "KOH_reactant", "OH-_reactant",
    "NaOMe_reactant", "NaOEt_reactant", "KOtBu_reactant", "RO-_reactant",
    
    # Other specific reagents
    "malonate_reactant", "acetoacetate_reactant",
}

def load_json(file_path: Path) -> dict:
    """Load JSON file"""
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def save_json(data: dict, file_path: Path):
    """Save JSON file with formatting"""
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

def create_backup(file_path: Path):
    """Create backup of file"""
    backup_path = file_path.with_suffix('.json.v4-backup')
    shutil.copy2(file_path, backup_path)
    print(f"✓ Created backup: {backup_path}")
    return backup_path

def extract_base_token(token: str) -> str:
    """Extract base token from _reactant suffix"""
    if token.endswith('_reactant'):
        return token[:-9]  # Remove '_reactant'
    return token

def find_present_equivalent(base_token: str, features: List[dict]) -> dict | None:
    """Find _present feature that matches the base token"""
    # Direct match
    present_token = f"{base_token}_present"
    for feat in features:
        if feat['token'] == present_token:
            return feat
    
    # Try without underscores/hyphens
    normalized = base_token.replace('-', '_').replace('_', '').lower()
    for feat in features:
        feat_normalized = feat['token'].replace('-', '_').replace('_', '').lower()
        if feat_normalized == normalized and feat['token'].endswith('_present'):
            return feat
    
    return None

def merge_reactant_metadata(present_feat: dict, reactant_feat: dict):
    """Add reactant_metadata to present feature from reactant feature"""
    metadata = reactant_feat.get('metadata', {})
    
    if not metadata:
        return  # No metadata to merge
    
    # Create reactant_metadata structure
    present_feat['reactant_metadata'] = {
        'is_reactant': True,
        'compound_type': metadata.get('reactant_name', ''),
        'reactant_category': metadata.get('reactant_category', ''),
        'reactant_member': metadata.get('reactant_member', ''),
        'coupling_role': metadata.get('coupling_role', ''),
        'reaction_families': [],  # Can be populated later
        'specificity': len(reactant_feat.get('detect', {}).get('smarts_any', [''])[0]) if 'detect' in reactant_feat else 0
    }
    
    # Clean up empty fields
    present_feat['reactant_metadata'] = {k: v for k, v in present_feat['reactant_metadata'].items() if v}

def refactor_features(json_path: Path) -> Dict[str, int]:
    """Main refactoring logic"""
    stats = {
        'total_features': 0,
        'reactant_features_found': 0,
        'merged_to_present': 0,
        'removed_specific_compounds': 0,
        'removed_unmatched_reactants': 0,
        'final_feature_count': 0,
    }
    
    # Load data
    data = load_json(json_path)
    features = data['features']
    stats['total_features'] = len(features)
    
    # Separate reactant and non-reactant features
    reactant_features = [f for f in features if f['token'].endswith('_reactant')]
    other_features = [f for f in features if not f['token'].endswith('_reactant')]
    
    stats['reactant_features_found'] = len(reactant_features)
    
    print(f"\nProcessing {len(reactant_features)} reactant features...")
    
    # Process each reactant feature
    for reactant_feat in reactant_features:
        token = reactant_feat['token']
        
        # Check if this is a specific compound to remove
        if token in SPECIFIC_COMPOUNDS_TO_REMOVE:
            stats['removed_specific_compounds'] += 1
            print(f"  ✗ Removing specific compound: {token}")
            continue
        
        # Try to find corresponding _present feature
        base_token = extract_base_token(token)
        present_feat = find_present_equivalent(base_token, other_features)
        
        if present_feat:
            # Merge metadata
            merge_reactant_metadata(present_feat, reactant_feat)
            stats['merged_to_present'] += 1
            print(f"  ✓ Merged {token} → {present_feat['token']}")
        else:
            # No _present equivalent found - keep as structural feature but remove reactant suffix
            # But skip if it's overly specific
            if any(specific in token for specific in ['Wittig', 'HWE', 'ylide', 'enolate', 'silyl_enol_ether', 'beta_dicarbonyl']):
                stats['removed_unmatched_reactants'] += 1
                print(f"  ✗ Removing unmatched specific: {token}")
            else:
                stats['removed_unmatched_reactants'] += 1
                print(f"  ⚠ No _present match for: {token} (removing)")
    
    # Remove derived expressions referencing removed _reactant tokens
    if 'derived' in data:
        derived = data['derived']
        updated_derived = []
        removed_derived = 0
        
        for deriv in derived:
            expr = deriv.get('derived', '')
            # Check if expression references removed tokens
            has_removed_token = any(token in expr for token in SPECIFIC_COMPOUNDS_TO_REMOVE)
            if has_removed_token:
                removed_derived += 1
                print(f"  ✗ Removing derived expression: {deriv['token']}")
            else:
                updated_derived.append(deriv)
        
        data['derived'] = updated_derived
        print(f"\nRemoved {removed_derived} derived expressions with specific compound references")
    
    # Update features list
    data['features'] = other_features
    stats['final_feature_count'] = len(other_features)
    
    # Update version
    data['version'] = 'v5.0-unified-metadata'
    
    return data, stats

def main():
    json_path = Path(__file__).parent.parent / 'chemtools' / 'featurizers' / 'calculable_features.json'
    
    print(f"Refactoring: {json_path}")
    
    # Create backup
    backup_path = create_backup(json_path)
    
    # Refactor
    refactored_data, stats = refactor_features(json_path)
    
    # Save
    save_json(refactored_data, json_path)
    
    # Print summary
    print("\n" + "="*70)
    print("REFACTORING COMPLETE")
    print("="*70)
    print(f"Original features:            {stats['total_features']}")
    print(f"Reactant features found:      {stats['reactant_features_found']}")
    print(f"  - Merged to _present:       {stats['merged_to_present']}")
    print(f"  - Removed (specific cpds):  {stats['removed_specific_compounds']}")
    print(f"  - Removed (unmatched):      {stats['removed_unmatched_reactants']}")
    print(f"Final feature count:          {stats['final_feature_count']}")
    print(f"Reduction:                    {stats['total_features'] - stats['final_feature_count']} features ({(stats['total_features'] - stats['final_feature_count'])/stats['total_features']*100:.1f}%)")
    print(f"\nBackup saved to: {backup_path}")
    print(f"Updated file: {json_path}")

if __name__ == '__main__':
    main()
