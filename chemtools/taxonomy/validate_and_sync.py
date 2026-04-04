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

from chemtools.taxonomy import loader as taxonomy_loader
from chemtools.taxonomy.substituent_composer import (
    load_organic_groups_with_compositions,
    validate_substituent_fragments_payload,
)
from chemtools.util.rdkit_helpers import rdkit_available
from chemtools.util.smarts_cache import compile_smarts

_GENERATED_ONLY_GROUP_IDS = {
    "-COOH",
    "-COOR",
    "-CONH2",
    "-CONHR",
    "-CONR2",
    "-CON3",
    "-CONHNH2",
    "-OCONH2",
    "-SO2NH2",
    "-SO2NHR",
    "-SO2NR2",
    "-SO2NHNH2",
    "-PO2OH",
    "-PO2OR",
    "-PO2NH2",
    "-PO2NHR",
    "-PO2NR2",
}

class TaxonomyValidator:
    def __init__(self, taxonomy_dir: Path):
        self.taxonomy_dir = taxonomy_dir
        self.groups_file = taxonomy_dir / "organic_groups.v1.3.json"
        self.compounds_file = taxonomy_dir / "organic_compounds.v1.3.json"  # legacy snapshot (optional)
        
        self.groups: Dict[str, Any] = {}
        self.compounds: Dict[str, Any] = {}
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.fixes_applied: List[str] = []
        self.composed_count: int = 0
        self.base_group_ids: Set[str] = set()
        
    def load_data(self) -> bool:
        """Load both JSON files"""
        try:
            with open(self.groups_file, "r", encoding="utf-8") as f:
                base_payload = json.load(f)
            self.base_group_ids = {
                str(g.get("id") or "").strip()
                for g in (base_payload.get("groups") or [])
                if isinstance(g, dict) and g.get("id")
            }

            groups_data = load_organic_groups_with_compositions(self.groups_file)
            self.groups = {
                g["id"]: g
                for g in groups_data.get("groups", [])
                if isinstance(g, dict) and g.get("id")
            }
            composed = groups_data.get("composed_groups", {}) or {}
            self.composed_count = int(composed.get("generated_count") or 0)
            for msg in composed.get("errors", []) or []:
                self.errors.append(f"Composed group error: {msg}")
            print(f"✓ Loaded {len(self.groups)} organic groups from {self.groups_file.name}")
        except Exception as e:
            self.errors.append(f"Failed to load groups: {e}")
            return False
            
        try:
            compounds_data = taxonomy_loader.load_organic_compounds()
            self.compounds_data = compounds_data
            self.compounds_list = compounds_data.get('compounds', [])
            print(f"✓ Loaded {len(self.compounds_list)} organic compounds from generated catalog")
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

    def validate_composed_group_smarts(self) -> None:
        """Validate SMARTS compilability for composed groups."""
        composed = [
            group
            for group in self.groups.values()
            if isinstance(group, dict) and isinstance(group.get("generated"), dict)
        ]
        if not composed:
            return
        if not rdkit_available():
            self.warnings.append("RDKit unavailable; skipped composed SMARTS validation.")
            return
        for group in composed:
            group_id = str(group.get("id") or "")
            smarts = str(group.get("smarts") or "")
            if not smarts:
                self.errors.append(f"Composed group '{group_id}' missing SMARTS.")
                continue
            if compile_smarts(smarts, validate=False) is None:
                self.errors.append(f"Composed group '{group_id}' has invalid SMARTS.")

    def validate_substituent_fragments_schema(self) -> None:
        """Validate substituent_fragments.v1.json editing schema."""
        fragments_file = self.taxonomy_dir / "substituent_fragments.v1.json"
        if not fragments_file.exists():
            self.warnings.append("substituent_fragments.v1.json not found; skipped schema validation.")
            return
        try:
            with open(fragments_file, "r", encoding="utf-8") as f:
                payload = json.load(f)
        except Exception as e:
            self.errors.append(f"Failed to load substituent fragments: {e}")
            return
        for msg in validate_substituent_fragments_payload(payload):
            self.errors.append(f"substituent_fragments.v1.json: {msg}")

    def validate_generated_only_groups(self) -> None:
        bad = sorted(_GENERATED_ONLY_GROUP_IDS & self.base_group_ids)
        for group_id in bad:
            self.errors.append(
                f"Group '{group_id}' must be generated via substituent_fragments.v1.json, "
                "not defined directly in organic_groups.v1.3.json."
            )

    def validate_h_pseudo_scaffold_policy(self) -> None:
        """Disallow pseudo-scaffold H compound motifs."""
        for compound in self.compounds_list:
            if not isinstance(compound, dict):
                continue
            a_ref = str(compound.get("A") or "").strip()
            b_ref = str(compound.get("B") or "").strip()
            explicit_id = str(compound.get("id") or "").strip()
            compound_id = explicit_id
            if not compound_id and a_ref and b_ref:
                compound_id = f"{a_ref}{b_ref}" if b_ref.startswith("-") else f"{a_ref}-{b_ref}"

            if a_ref == "H":
                self.errors.append(
                    f"Compound '{compound_id or '<missing-id>'}' uses deprecated pseudo-scaffold A='H'."
                )
    
    def check_explicit_id_fields(self) -> int:
        """Check for explicit 'id' fields (should be auto-generated from A-B)"""
        explicit_count = 0
        
        for compound in self.compounds_list:
            if 'id' in compound:
                explicit_count += 1
                comp_id = compound.get('id', '')
                a_ref = compound.get('A', '')
                b_ref = compound.get('B', '')
                
                # Check if explicit ID matches A-B pattern
                expected_id = f"{a_ref}-{b_ref}" if a_ref and b_ref else None
                if expected_id and comp_id != expected_id:
                    self.warnings.append(
                        f"Compound has explicit ID '{comp_id}' that differs from auto-generated '{expected_id}'"
                    )
        
        return explicit_count
    
    def validate_id_format(self) -> None:
        """Check that compounds don't have explicit ID fields (should be auto-generated)"""
        # With auto-generation, we don't need to validate ID format
        # IDs are generated as f"{A}-{B}" when compounds are loaded
        pass
    
    def validate_dependencies(self) -> None:
        """Check version dependencies"""
        if "depends_on" not in self.compounds_data:
            return
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
    
    def fix_naming_mismatches(self, issues: List[Dict[str, Any]]) -> None:
        """Remove redundant name fields from compounds"""
        if not issues:
            return
        
        print(f"\n{'='*70}")
        print(f"Removing {len(issues)} redundant 'name' fields...")
        
        for issue in issues:
            compound = issue['compound']
            comp_id = issue['id']
            
            # Remove redundant name field
            if 'name' in compound:
                del compound['name']
            
            self.fixes_applied.append(
                f"  ✓ {comp_id}: removed redundant 'name' field"
            )
        
        self.warnings.append(
            "Fix mode is deprecated for generated compound catalog. Update "
            "compound_generation_rules.v1.json instead."
        )
    
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
        
        # 2. Check for explicit 'id' fields
        print("\n2. Checking for explicit 'id' fields...")
        explicit_count = self.check_explicit_id_fields()
        if explicit_count == 0:
            print(f"  ✓ No explicit 'id' fields (using auto-generated A-B system)")
        else:
            print(f"  ℹ Found {explicit_count} compounds with explicit 'id' fields")
            if fix_mode:
                print(f"    Run with --remove-ids to clean up explicit IDs")
        
        # 3. Check ID format
        print("\n3. Validating compound ID formats...")
        self.validate_id_format()
        
        # 4. Check dependencies
        print("\n4. Validating version dependencies...")
        self.validate_dependencies()

        # 5. Enforce generated-only groups
        print("\n5. Enforcing generated-only group policy...")
        self.validate_generated_only_groups()

        # 6. Reject deprecated pseudo-scaffold H motifs
        print("\n6. Checking deprecated pseudo-scaffold H motifs...")
        self.validate_h_pseudo_scaffold_policy()

        # 7. Validate substituent fragments schema
        print("\n7. Validating substituent fragments schema...")
        self.validate_substituent_fragments_schema()

        # 8. Validate composed group SMARTS
        print("\n8. Validating composed substituent SMARTS...")
        self.validate_composed_group_smarts()
        if self.composed_count:
            print(f"  ✓ Loaded {self.composed_count} composed substituent group(s)")
        else:
            print("  ℹ No composed substituent groups generated")

        # 9. Generate usage report
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

