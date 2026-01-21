#!/usr/bin/env python3
"""Chemistry validation for taxonomy organic groups using RDKit."""

import json
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("RDKit not available - skipping SMARTS validation")

def validate_smarts(smarts, group_id):
    """Validate a SMARTS pattern with RDKit."""
    if not HAS_RDKIT:
        return True, "RDKit not available"
    
    try:
        mol = Chem.MolFromSmarts(smarts)
        if mol is None:
            return False, "Invalid SMARTS - RDKit returned None"
        return True, "Valid"
    except Exception as e:
        return False, f"Exception: {str(e)}"

def check_overlap_issues(groups_data):
    """Check for potential overlap issues between groups."""
    if not HAS_RDKIT:
        return []
    
    issues = []
    groups = groups_data['groups']
    
    # Check scaffold overlaps (Ar vs HeteroAr, Alkyl vs specific alkyls, etc.)
    scaffolds = {g['id']: g for g in groups if g.get('kind') == 'scaffold'}
    
    # Specific checks
    if 'Ar' in scaffolds and 'HeteroAr' in scaffolds:
        # HeteroAr is more specific - should have higher or equal priority
        ar_pri = scaffolds['Ar'].get('priority', 1)
        het_pri = scaffolds['HeteroAr'].get('priority', 1)
        if het_pri < ar_pri:
            issues.append(f"Priority issue: HeteroAr ({het_pri}) should have priority >= Ar ({ar_pri}) since it's more specific")
    
    if 'Alkyl' in scaffolds and 'RCH2' in scaffolds:
        # RCH2/R2CH/R3C are more specific than generic Alkyl
        alkyl_pri = scaffolds['Alkyl'].get('priority', 1)
        for specific in ['RCH2', 'R2CH', 'R3C']:
            if specific in scaffolds:
                spec_pri = scaffolds[specific].get('priority', 1)
                if spec_pri < alkyl_pri:
                    issues.append(f"Priority issue: {specific} ({spec_pri}) should have priority >= Alkyl ({alkyl_pri})")
    
    # Check for Alkynyl vs Alkynyl_terminal
    if 'Alkynyl' in scaffolds and 'Alkynyl_terminal' in scaffolds:
        alk_pri = scaffolds['Alkynyl'].get('priority', 1)
        term_pri = scaffolds['Alkynyl_terminal'].get('priority', 1)
        if term_pri < alk_pri:
            issues.append(f"Priority issue: Alkynyl_terminal ({term_pri}) should have priority >= Alkynyl ({alk_pri}) since it's more specific")
    
    return issues

def main():
    groups_file = Path('c:/Git-softwares/Condition-agent/chemtools/taxonomy/data/organic_groups.v1.3.json')
    with open(groups_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    print("=" * 70)
    print("CHEMISTRY VALIDATION FOR ORGANIC GROUPS TAXONOMY")
    print("=" * 70)
    print()

    if not HAS_RDKIT:
        print("⚠️  RDKit not available - limited validation possible")
        print()
    
    # Validate all SMARTS patterns
    print("SMARTS Pattern Validation:")
    invalid_patterns = []
    
    for g in data['groups']:
        gid = g.get('id', '?')
        smarts = g.get('smarts', '')
        valid, msg = validate_smarts(smarts, gid)
        
        if not valid:
            invalid_patterns.append((gid, smarts, msg))
    
    if invalid_patterns:
        print(f"  ❌ Found {len(invalid_patterns)} invalid SMARTS patterns:")
        for gid, smarts, msg in invalid_patterns:
            print(f"     {gid:20s}: {msg}")
            print(f"       Pattern: {smarts}")
    else:
        print("  ✅ All SMARTS patterns are valid")
    print()

    # Check overlap/priority issues
    print("Logical Priority Checks:")
    overlap_issues = check_overlap_issues(data)
    
    if overlap_issues:
        print(f"  ⚠️  Found {len(overlap_issues)} potential priority issues:")
        for issue in overlap_issues:
            print(f"     {issue}")
    else:
        print("  ✅ No obvious priority conflicts detected")
    print()

    # Check chemistry-specific patterns
    print("Chemistry-Specific Checks:")
    chem_issues = []
    
    for g in data['groups']:
        gid = g.get('id', '')
        smarts = g.get('smarts', '')
        desc = g.get('description', '')
        
        # Check halogen patterns
        if gid in ['-Cl', '-Br', '-I', '-F']:
            if not (smarts == f'[{gid[1:]}:2]'):
                chem_issues.append(f"{gid}: Halogen SMARTS should be simpler")
        
        # Check -CF3 pattern
        if gid == '-CF3':
            if 'F)(F)F' not in smarts:
                chem_issues.append(f"{gid}: CF3 should have three fluorines")
        
        # Check boronic acid vs ester distinction
        if gid == '-B(OH)2':
            if 'O[H]' not in smarts and 'OH' not in smarts:
                chem_issues.append(f"{gid}: Boronic acid should have OH groups")
        
        # Check sulfonyl vs sulfone
        if 'SO2' in gid and gid != '-SO2R':
            if '=O)(=O)' not in smarts:
                chem_issues.append(f"{gid}: Sulfonyl group should have two oxygens with =O")
    
    if chem_issues:
        print(f"  ⚠️  Found {len(chem_issues)} chemistry issues:")
        for issue in chem_issues:
            print(f"     {issue}")
    else:
        print("  ✅ Chemistry patterns look consistent")
    print()

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    total_issues = len(invalid_patterns) + len(overlap_issues) + len(chem_issues)
    
    if total_issues == 0:
        print("✅ Taxonomy system looks chemically sound!")
    else:
        print(f"⚠️  Found {total_issues} issues that should be reviewed")
    
    print()
    print(f"Total groups: {len(data['groups'])}")
    print(f"  Scaffolds: {len([g for g in data['groups'] if g.get('kind') == 'scaffold'])}")
    print(f"  Substituents: {len([g for g in data['groups'] if g.get('kind') == 'substituent'])}")

if __name__ == '__main__':
    main()
