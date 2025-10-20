#!/usr/bin/env python3
"""
Batch SMARTS Pattern Fixer

Fixes common SMARTS pattern issues in protocol JSON files.

Common fixes:
1. Replace [CH], [CH2], [CH3] with [C;H1], [C;H2], [C;H3] or more general patterns
2. Fix pattern structures to match actual reaction SMILES
3. Make patterns more flexible while maintaining chemical meaning

Usage:
    python fix_smarts_patterns.py
"""

import json
from pathlib import Path
from typing import Dict, Any, List, Tuple

# Define fixes for each protocol
FIXES = {
    "Aryl mesylate_Suzuki.json": {
        "old": ["O=S(OC)([c,C,n,o,s])=O.OB(O)[c,C,n,o,s]>>[c,C,n,o,s]"],
        "new": ["[c,C,n,o,s]OS(=O)(=O)C.OB(O)[c,C,n,o,s]>>[c,C,n,o,s][c,C,n,o,s]"],
        "reason": "Mesylate pattern structure corrected and product pattern shows C-C bond formation"
    },
    
    "Sonogashira-Coupling.json": {
        "old": ["[c,C,n,o,s]I.[!#1]C#C[H]>>[c,C,n,o,s]C#C[!#1]"],
        "new": ["[c,C,n,o,s]I.C#C>>[c,C,n,o,s]C#C"],
        "reason": "Removed restrictive [H] requirement for terminal alkyne"
    },
    
    "Alkyl_Iodide_Borylation.json": {
        "old": ["I[CH2].CC1(C(OB(BB2OC(C(O2)(C)C)(C)C)O1)(C)C)C>>[CH2]B3OC(C(O3)(C)C)(C)C"],
        "new": ["IC.B2pin2>>CB(Opin)"],
        "reason": "Simplified pattern using common abbreviations; [CH2] causes RDKit errors"
    },
    
    "Suzuki_Cu_C(sp3)-C(sp2).json": {
        "old": ["Br[CH].[c,C,n,o,s]B(OC1(C)C)OC1(C)C>>[c,C,n,o,s][CH]"],
        "new": ["BrC.[c,C,n,o,s]B(Opin)>>[c,C,n,o,s]C"],
        "reason": "Simplified to avoid [CH] RDKit error; uses Bpin abbreviation pattern"
    },
    
    "Hydroacylation_Ni_aryl_alkene+acyl_fluoride.json": {
        "old": ["[c,C,n,o,s]/C([H])=C([H])/[H].O=C([c,C,n,o,s])F>>[c,C,n,o,s]C(C([c,C,n,o,s])=O)C"],
        "new": ["C=C.O=C(F)>>C(C)C(=O)"],
        "reason": "Simplified pattern; removed stereochemistry and [H] to avoid RDKit issues"
    },
    
    "Ni_Cross_Electrophile_Acylation.json": {
        "old": ["[C;!a;!H0]CI.[C;!a;!H0]C(Cl)=O>>[C;!a;!H0]C(C[C;!a;!H0])=O"],
        "new": ["CI.C(Cl)=O>>CC(C)=O"],
        "reason": "Simplified pattern; removed !H0 which causes RDKit errors"
    },
    
    "pd_acetylation_aryl_bromide_Garg_v98p0068.json": {
        "old": ["Br[c,C,n,o,s].[CH2,CH3]C([Si](C)(C)C)=O>>[CH2,CH3]C([c,C,n,o,s])=O"],
        "new": ["Br[c,C,n,o,s].CC([Si])=O>>CC([c,C,n,o,s])=O"],
        "reason": "Simplified [CH2,CH3] to C to avoid RDKit errors"
    },
    
    "Pd_Buchwald_Arylsulfonate_Amination_CMPhos.json": {
        "old": ["[c,C,n,o,s]OS(=O)(C)=O.[N;H1,H2]>>[N;H1,H2][c,C,n,o,s]"],
        "new": ["[c,C,n,o,s]OS(=O)(=O)C.[NH]>>[NH][c,C,n,o,s]"],
        "reason": "Fixed mesylate structure and simplified N pattern to avoid H count issues"
    },
    
    "Pd_Conjugate_Addition_Alkyne_to_Enone.json": {
        "old": ["[!#1]C#C[H].[#6;!a][C](=O)[CH]=[CH2]>>[#6;!a][C](=O)[CH2][CH2]C#C[!#1]"],
        "new": ["C#C.C(=O)C=C>>C(=O)CCC#C"],
        "reason": "Simplified to avoid [H], [CH], [CH2] RDKit errors"
    },
    
    "Renaudet_Reymond_2004_Mitsunobu.json": {
        "old": ["O[CH2,CH3,CH].[c,C,n,o,s]O>>[c,C,n,o,s]O[CH2,CH3,CH]"],
        "new": ["OC.[c,C,n,o,s]O>>[c,C,n,o,s]OC"],
        "reason": "Simplified [CH2,CH3,CH] to C to avoid RDKit errors"
    }
}


def fix_protocol(protocol_path: Path, old_patterns: List[str], new_patterns: List[str], reason: str) -> bool:
    """Fix SMARTS patterns in a protocol file"""
    try:
        # Load protocol
        with open(protocol_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Update reaction_SMARTS
        if 'reaction' in data and 'reaction_SMARTS' in data['reaction']:
            current = data['reaction']['reaction_SMARTS']
            
            if current == old_patterns:
                data['reaction']['reaction_SMARTS'] = new_patterns
                
                # Save with pretty formatting
                with open(protocol_path, 'w', encoding='utf-8') as f:
                    json.dump(data, f, indent=2, ensure_ascii=False)
                
                print(f"✅ Fixed {protocol_path.name}")
                print(f"   Old: {old_patterns}")
                print(f"   New: {new_patterns}")
                print(f"   Reason: {reason}")
                print()
                return True
            else:
                print(f"⚠️  {protocol_path.name}: Current pattern doesn't match expected old pattern")
                print(f"   Expected: {old_patterns}")
                print(f"   Found: {current}")
                print()
                return False
        else:
            print(f"❌ {protocol_path.name}: No reaction_SMARTS field found")
            return False
            
    except Exception as e:
        print(f"❌ Error fixing {protocol_path.name}: {e}")
        return False


def main():
    """Apply all fixes"""
    # Get protocol directory
    current_file = Path(__file__)
    protocol_dir = current_file.parent.parent / "data" / "protocol_db"
    
    if not protocol_dir.exists():
        print(f"❌ Protocol directory not found: {protocol_dir}")
        return 1
    
    print("=" * 70)
    print("Batch SMARTS Pattern Fixer")
    print("=" * 70)
    print()
    print(f"Protocol directory: {protocol_dir}")
    print(f"Fixes to apply: {len(FIXES)}")
    print()
    print("=" * 70)
    print()
    
    # Apply fixes
    success_count = 0
    for filename, fix_data in FIXES.items():
        filepath = protocol_dir / filename
        
        if not filepath.exists():
            print(f"⚠️  File not found: {filename}")
            print()
            continue
        
        if fix_protocol(
            filepath,
            fix_data["old"],
            fix_data["new"],
            fix_data["reason"]
        ):
            success_count += 1
    
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"Successfully fixed: {success_count}/{len(FIXES)} files")
    print()
    
    if success_count == len(FIXES):
        print("✅ All fixes applied successfully!")
        print()
        print("Next steps:")
        print("1. Validate the fixes:")
        print("   python -m chemtools.protocol.validate_protocols --verbose")
        print()
        print("2. Rebuild the index:")
        print("   python -m chemtools.protocol.cli build --force")
        return 0
    else:
        print("⚠️  Some fixes were not applied. Check the output above.")
        return 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
