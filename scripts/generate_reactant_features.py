#!/usr/bin/env python3
"""
Legacy: reactant_types.json is retired.

Reactant type SMARTS are now derived from organic_compounds.v1.3.json via
chemtools.taxonomy.calculable_spec. This script is kept for reference only.
"""

import json
from pathlib import Path
from typing import Dict, List

# Paths
REPO_ROOT = Path(__file__).parent.parent
REACTANT_TYPES_FILE = REPO_ROOT / "chemtools" / "taxonomy" / "data" / "reactant_types.json"
OUTPUT_FILE = REPO_ROOT / "scripts" / "reactant_features_generated.json"

# Role mapping for coupling chemistry
ROLE_MAPPING = {
    # Electrophiles
    "ArX*": "electrophile",
    "HetAr-X": "electrophile",
    "VinylX*": "electrophile",
    "Alkyl-X": "electrophile",
    "Allylic-X": "electrophile",
    "Benzylic-X": "electrophile",
    "RSO2Cl": "electrophile",
    "Acyl-electrophile": "electrophile",
    
    # Nucleophiles / Coupling partners
    "ArB*": "nucleophile",
    "RB*": "nucleophile",
    "R-M": "nucleophile",
    "RMgX": "nucleophile",
    "RZnX": "nucleophile",
    "RLi": "nucleophile",
    "RNH2/R2NH": "nucleophile",
    "ArNH2/Ar2NH": "nucleophile",
    "ROH": "nucleophile",
    "ArOH": "nucleophile",
    "RSH": "nucleophile",
    
    # Both roles
    "Alkyne": "both",
    "Alkene": "both",
    "Aldehyde": "electrophile",
    "Ketone": "electrophile",
    "Ester": "electrophile",
    "RCO2H": "nucleophile",
    "Amide": "nucleophile",
    
    # Special
    "Metal-H": "reductant",
    "H2": "reductant",
    "Oxidant": "oxidant",
    "Diene": "dienophile",
    "Dienophile": "dienophile",
}


def load_reactant_types() -> List[dict]:
    """Load reactant types from JSON."""
    if not REACTANT_TYPES_FILE.exists():
        raise SystemExit(
            "reactant_types.json has been removed; use organic_compounds.v1.3.json instead."
        )
    with open(REACTANT_TYPES_FILE) as f:
        return json.load(f)


def generate_member_features(reactant_data: List[dict]) -> List[dict]:
    """Generate individual member features."""
    features = []
    
    for category in reactant_data:
        cat_id = category['id']
        cat_name = category['name']
        cat_desc = category.get('description', '')
        role = ROLE_MAPPING.get(cat_id, 'other')
        
        for member in category.get('members', []):
            member_id = member['id']
            member_name = member['name']
            smarts = member.get('smarts', '')
            
            if not smarts:
                continue
            
            # Generate feature token
            token = f"{member_id}_reactant"
            
            # Create feature definition
            feature = {
                "token": token,
                "type": "bool",
                "scope": "global",
                "category": "reactant_types",
                "detect": {
                    "smarts_any": [smarts]
                },
                "why": f"{member_name} - {cat_name}",
                "metadata": {
                    "reactant_category": cat_id,
                    "reactant_member": member_id,
                    "reactant_name": member_name,
                    "category_name": cat_name,
                    "coupling_role": role,
                    "legacy_taxonomy_id": member_id,
                    "category_description": cat_desc
                }
            }
            
            features.append(feature)
    
    return features


def generate_category_features(reactant_data: List[dict]) -> List[dict]:
    """Generate category-level derived features."""
    derived_features = []
    
    for category in reactant_data:
        cat_id = category['id']
        cat_name = category['name']
        members = category.get('members', [])
        
        if not members:
            continue
        
        # Generate category token
        # Replace special characters in category ID
        safe_cat_id = cat_id.replace('*', '').replace('-', '_').replace('/', '_')
        token = f"{safe_cat_id}_reactant"
        
        # Build OR expression from all members
        member_tokens = [f"{m['id']}_reactant" for m in members if m.get('smarts')]
        
        if not member_tokens:
            continue
        
        derive_expr = " OR ".join(member_tokens)
        
        # Create derived feature
        derived = {
            "token": token,
            "derive": derive_expr,
            "why": f"Category-level: {cat_name}",
            "metadata": {
                "reactant_category": cat_id,
                "category_name": cat_name,
                "is_category_level": True,
                "member_count": len(member_tokens)
            }
        }
        
        derived_features.append(derived)
    
    return derived_features


def main():
    print("=" * 80)
    print("Generating Reactant Type Features for calculable_features.json")
    print("=" * 80)
    
    # Load reactant types
    print(f"\nLoading reactant types from: {REACTANT_TYPES_FILE}")
    reactant_data = load_reactant_types()
    print(f"Loaded {len(reactant_data)} reactant categories")
    
    # Generate member features
    print("\nGenerating member-level features...")
    member_features = generate_member_features(reactant_data)
    print(f"Generated {len(member_features)} member features")
    
    # Generate category features
    print("\nGenerating category-level derived features...")
    category_features = generate_category_features(reactant_data)
    print(f"Generated {len(category_features)} category features")
    
    # Combine results
    output = {
        "features": member_features,
        "derived_shortcuts": category_features,
        "summary": {
            "total_categories": len(reactant_data),
            "total_member_features": len(member_features),
            "total_category_features": len(category_features),
            "total_features": len(member_features) + len(category_features)
        }
    }
    
    # Save to file
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n✓ Generated features saved to: {OUTPUT_FILE}")
    print("\nSummary:")
    print(f"  - Member features: {len(member_features)}")
    print(f"  - Category features: {len(category_features)}")
    print(f"  - Total features: {len(member_features) + len(category_features)}")
    
    print("\nSample member features:")
    for feat in member_features[:5]:
        print(f"  • {feat['token']}: {feat['metadata']['reactant_name']}")
    
    print("\nSample category features:")
    for feat in category_features[:5]:
        print(f"  • {feat['token']}: {feat['metadata']['category_name']}")
    
    print("\n" + "=" * 80)
    print("Next steps:")
    print("1. Review generated features in: scripts/reactant_features_generated.json")
    print("2. Merge into chemtools/taxonomy/data/calculable_features.json")
    print("3. Update calculable.py to handle reactant type queries")
    print("=" * 80)


if __name__ == "__main__":
    main()
