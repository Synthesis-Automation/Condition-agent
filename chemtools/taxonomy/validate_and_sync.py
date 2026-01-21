#!/usr/bin/env python3
"""
Validate and synchronize organic_compounds.v1.3.json with organic_groups.v1.3.json

This script ensures that:
1. All group references (A, B) in compounds exist in organic_groups
2. Compound IDs follow the pattern A-B (using actual group IDs)
3. Compound names match their IDs (simplified naming system)
4. Reports missing groups, invalid references, and inconsistencies
5. Can optionally auto-fix naming mismatches

Usage:
    python validate_and_sync.py                  # Validate only
    python validate_and_sync.py --fix            # Fix naming issues
    python validate_and_sync.py --check-only     # Exit with code 1 if issues found
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any
from collections import defaultdict


class TaxonomyValidator:
    def __init__(self, taxonomy_dir: Path):
        self.taxonomy_dir = taxonomy_dir
        self.groups_file = taxonomy_dir / "organic_groups.v1.3.json"
        self.compounds_file = taxonomy_dir / "organic_compounds.v1.3.json"
        
        self.groups: Dict[str, Any] = {}
        self.compounds: Dict[str, Any] = {}
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.fixes_applied: List[str] = []
        
    def load_data(self) -> bool:
        """Load both JSON files"""
        try:
            with open(self.groups_file, 'r', encoding='utf-8') as f:
                groups_data = json.load(f)
                self.groups = {g['id']: g for g in groups_data.get('groups', [])}
                print(f"✓ Loaded {len(self.groups)} organic groups from {self.groups_file.name}")
        except Exception as e:
            self.errors.append(f"Failed to load groups: {e}")
            return False
            
        try:
            with open(self.compounds_file, 'r', encoding='utf-8') as f:
                compounds_data = json.load(f)
                self.compounds_data = compounds_data
                self.compounds_list = compounds_data.get('compounds', [])
                print(f"✓ Loaded {len(self.compounds_list)} organic compounds from {self.compounds_file.name}")
        except Exception as e:
            self.errors.append(f"Failed to load compounds: {e}")
            return False
            
        return True
    
    def validate_group_references(self) -> Tuple[Set[str], Set[str]]:
        """Validate that all A and B references exist in groups"""
        missing_refs = set()
        valid_refs = set()
        
        # Special handling for substituent variants (e.g., Ar_Subst, Alkyl_Subst)
        # These reference the same group but indicate it's used as a substituent
        def normalize_group_ref(ref: str) -> str:
            """Remove _Subst suffix for validation"""
            if ref.endswith('_Subst'):
                return ref.replace('_Subst', '')
            return ref
        
        for compound in self.compounds_list:
            comp_id = compound.get('id', 'unknown')
            a_ref = compound.get('A')
            b_ref = compound.get('B')
            
            # Check A reference
            if a_ref:
                normalized_a = normalize_group_ref(a_ref)
                if normalized_a in self.groups:
                    valid_refs.add(a_ref)
                else:
                    missing_refs.add(a_ref)
                    self.errors.append(f"Compound '{comp_id}': Group A='{a_ref}' not found in organic_groups")
            
            # Check B reference
            if b_ref:
                normalized_b = normalize_group_ref(b_ref)
                if normalized_b in self.groups:
                    valid_refs.add(b_ref)
                else:
                    missing_refs.add(b_ref)
                    self.errors.append(f"Compound '{comp_id}': Group B='{b_ref}' not found in organic_groups")
        
        return valid_refs, missing_refs
    
    def validate_naming_consistency(self) -> List[Dict[str, Any]]:
        """Check if compound names match their IDs"""
        mismatches = []
        
        for compound in self.compounds_list:
            comp_id = compound.get('id', '')
            comp_name = compound.get('name', '')
            
            if comp_id != comp_name:
                mismatches.append({
                    'id': comp_id,
                    'current_name': comp_name,
                    'expected_name': comp_id,
                    'compound': compound
                })
                self.warnings.append(
                    f"Compound '{comp_id}': name='{comp_name}' should be '{comp_id}' "
                    f"(use --fix to auto-correct)"
                )
        
        return mismatches
    
    def validate_id_format(self) -> None:
        """Warn about compound IDs that don't follow A-B pattern"""
        for compound in self.compounds_list:
            comp_id = compound.get('id', '')
            a_ref = compound.get('A', '')
            b_ref = compound.get('B', '')
            
            # Expected format: A-B
            expected_id = f"{a_ref}-{b_ref}" if a_ref and b_ref else None
            
            if expected_id and comp_id != expected_id:
                self.warnings.append(
                    f"Compound '{comp_id}': ID doesn't match A-B pattern "
                    f"(expected '{expected_id}' from A='{a_ref}', B='{b_ref}')"
                )
    
    def validate_dependencies(self) -> None:
        """Check version dependencies"""
        depends_on = self.compounds_data.get('depends_on', {})
        groups_version = depends_on.get('organic_groups')
        
        if not groups_version:
            self.warnings.append("No 'depends_on.organic_groups' version specified in compounds file")
        elif groups_version != 'v1.3':
            self.warnings.append(
                f"Compounds depend on organic_groups '{groups_version}' but using 'v1.3'"
            )
    
    def generate_usage_report(self, valid_refs: Set[str]) -> None:
        """Show which groups are used and which are unused"""
        # Normalize references (remove _Subst suffix)
        normalized_refs = set()
        for ref in valid_refs:
            if ref.endswith('_Subst'):
                normalized_refs.add(ref.replace('_Subst', ''))
            else:
                normalized_refs.add(ref)
        
        unused_groups = set(self.groups.keys()) - normalized_refs
        
        print(f"\n{'='*70}")
        print(f"Group Usage Statistics:")
        print(f"  Total groups defined: {len(self.groups)}")
        print(f"  Groups used in compounds: {len(normalized_refs)}")
        print(f"  Unused groups: {len(unused_groups)}")
        
        if unused_groups:
            print(f"\n  Unused groups:")
            for group_id in sorted(unused_groups):
                group = self.groups[group_id]
                print(f"    - {group_id:20s} ({group.get('name', 'N/A')})")
    
    def fix_naming_mismatches(self, mismatches: List[Dict[str, Any]]) -> None:
        """Auto-fix compounds where name != id"""
        if not mismatches:
            return
        
        print(f"\n{'='*70}")
        print(f"Fixing {len(mismatches)} naming mismatches...")
        
        for mismatch in mismatches:
            compound = mismatch['compound']
            old_name = compound['name']
            new_name = mismatch['expected_name']
            
            # Update name to match id
            compound['name'] = new_name
            
            # Preserve old name in aliases if not already there
            aliases = compound.get('aliases', [])
            if old_name not in aliases and old_name != new_name:
                aliases.append(old_name)
                compound['aliases'] = aliases
            
            self.fixes_applied.append(
                f"  ✓ {mismatch['id']}: '{old_name}' → '{new_name}' (old name preserved in aliases)"
            )
        
        # Save the updated compounds file
        try:
            with open(self.compounds_file, 'w', encoding='utf-8') as f:
                json.dump(self.compounds_data, f, indent=2, ensure_ascii=False)
            print(f"\n✓ Updated {self.compounds_file.name} with {len(mismatches)} fixes")
        except Exception as e:
            self.errors.append(f"Failed to save fixes: {e}")
    
    def validate_all(self, fix_mode: bool = False) -> bool:
        """Run all validations"""
        print(f"\n{'='*70}")
        print("Validating Organic Taxonomy Consistency")
        print(f"{'='*70}\n")
        
        if not self.load_data():
            return False
        
        # 1. Check group references
        print("1. Validating group references (A, B)...")
        valid_refs, missing_refs = self.validate_group_references()
        if missing_refs:
            print(f"  ✗ Found {len(missing_refs)} missing group references")
        else:
            print(f"  ✓ All group references are valid")
        
        # 2. Check naming consistency
        print("\n2. Validating naming consistency (name == id)...")
        mismatches = self.validate_naming_consistency()
        if mismatches:
            print(f"  ⚠ Found {len(mismatches)} naming mismatches")
            if fix_mode:
                self.fix_naming_mismatches(mismatches)
        else:
            print(f"  ✓ All compound names match their IDs")
        
        # 3. Check ID format
        print("\n3. Validating compound ID formats...")
        self.validate_id_format()
        
        # 4. Check dependencies
        print("\n4. Validating version dependencies...")
        self.validate_dependencies()
        
        # 5. Generate usage report
        self.generate_usage_report(valid_refs)
        
        return True
    
    def print_summary(self) -> int:
        """Print final summary and return exit code"""
        print(f"\n{'='*70}")
        print("Validation Summary")
        print(f"{'='*70}")
        
        if self.fixes_applied:
            print(f"\n✓ Fixes Applied ({len(self.fixes_applied)}):")
            for fix in self.fixes_applied:
                print(fix)
        
        if self.warnings:
            print(f"\n⚠ Warnings ({len(self.warnings)}):")
            for warning in self.warnings:
                print(f"  - {warning}")
        
        if self.errors:
            print(f"\n✗ Errors ({len(self.errors)}):")
            for error in self.errors:
                print(f"  - {error}")
            print(f"\nValidation FAILED with {len(self.errors)} error(s).")
            return 1
        elif self.warnings and not self.fixes_applied:
            print(f"\nValidation completed with {len(self.warnings)} warning(s).")
            return 0
        else:
            print(f"\n✓ Validation PASSED - Taxonomy is consistent!")
            return 0


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Validate and synchronize organic taxonomy files'
    )
    parser.add_argument(
        '--fix',
        action='store_true',
        help='Automatically fix naming mismatches (name != id)'
    )
    parser.add_argument(
        '--check-only',
        action='store_true',
        help='Exit with code 1 if any issues found (useful for CI/CD)'
    )
    parser.add_argument(
        '--taxonomy-dir',
        type=Path,
        help='Path to taxonomy data directory (default: auto-detect)'
    )
    
    args = parser.parse_args()
    
    # Auto-detect taxonomy directory
    if args.taxonomy_dir:
        taxonomy_dir = args.taxonomy_dir
    else:
        script_dir = Path(__file__).parent
        taxonomy_dir = script_dir / "data"
        if not taxonomy_dir.exists():
            print(f"Error: Cannot find taxonomy directory at {taxonomy_dir}")
            print("Please specify --taxonomy-dir")
            sys.exit(1)
    
    validator = TaxonomyValidator(taxonomy_dir)
    
    success = validator.validate_all(fix_mode=args.fix)
    if not success:
        sys.exit(1)
    
    exit_code = validator.print_summary()
    
    if args.check_only and (validator.errors or validator.warnings):
        sys.exit(1)
    
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
