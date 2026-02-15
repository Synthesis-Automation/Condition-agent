#!/usr/bin/env python3
"""
Auto-generate compound suggestions when organic groups are added or modified.

This script helps maintain the organic_compounds file by:
1. Detecting new groups that aren't used in any compounds
2. Suggesting potential compound combinations
3. Identifying groups that should typically pair together

Usage:
    python suggest_compounds.py                           # Show suggestions
    python suggest_compounds.py --generate               # Generate JSON for new compounds
    python suggest_compounds.py --scaffold Ar            # Suggest compounds for specific scaffold
    python suggest_compounds.py --substituent Cl         # Suggest compounds for specific substituent
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any

from chemtools.taxonomy.substituent_composer import load_organic_groups_with_compositions


class CompoundSuggester:
    def __init__(self, taxonomy_dir: Path):
        self.taxonomy_dir = taxonomy_dir
        self.groups_file = taxonomy_dir / "organic_groups.v1.3.json"
        self.compounds_file = taxonomy_dir / "organic_compounds.v1.3.json"
        
        self.groups: Dict[str, Any] = {}
        self.scaffolds: Dict[str, Any] = {}
        self.substituents: Dict[str, Any] = {}
        self.compounds_list: List[Any] = []
        self.existing_pairs: Set[Tuple[str, str]] = set()
        
    def load_data(self) -> bool:
        """Load and categorize groups"""
        try:
            groups_data = load_organic_groups_with_compositions(self.groups_file)
            for group in groups_data.get('groups', []):
                group_id = group['id']
                self.groups[group_id] = group
                
                # Categorize by kind
                kind = group.get('kind', 'unknown')
                if kind == 'scaffold':
                    self.scaffolds[group_id] = group
                elif kind == 'substituent':
                    self.substituents[group_id] = group
                        
            print(f"✓ Loaded {len(self.groups)} groups ({len(self.scaffolds)} scaffolds, {len(self.substituents)} substituents)")
        except Exception as e:
            print(f"✗ Failed to load groups: {e}")
            return False
            
        try:
            with open(self.compounds_file, 'r', encoding='utf-8') as f:
                compounds_data = json.load(f)
                self.compounds_list = compounds_data.get('compounds', [])
                
                # Track existing A-B pairs
                for compound in self.compounds_list:
                    a = compound.get('A')
                    b = compound.get('B')
                    if a and b:
                        self.existing_pairs.add((a, b))
                        
            print(f"✓ Loaded {len(self.compounds_list)} existing compounds")
        except Exception as e:
            print(f"✗ Failed to load compounds: {e}")
            return False
            
        return True
    
    def normalize_group_ref(self, ref: str) -> str:
        """Remove _Subst suffix"""
        return ref.replace('_Subst', '') if ref.endswith('_Subst') else ref
    
    def get_used_groups(self) -> Tuple[Set[str], Set[str]]:
        """Get sets of used scaffold and substituent groups"""
        used_scaffolds = set()
        used_substituents = set()
        
        for compound in self.compounds_list:
            a = compound.get('A')
            b = compound.get('B')
            
            if a:
                normalized_a = self.normalize_group_ref(a)
                if normalized_a in self.scaffolds:
                    used_scaffolds.add(normalized_a)
            
            if b:
                normalized_b = self.normalize_group_ref(b)
                if normalized_b in self.substituents:
                    used_substituents.add(normalized_b)
        
        return used_scaffolds, used_substituents
    
    def find_unused_groups(self) -> Tuple[List[str], List[str]]:
        """Find groups that aren't used in any compounds"""
        used_scaffolds, used_substituents = self.get_used_groups()
        
        unused_scaffolds = [gid for gid in self.scaffolds.keys() if gid not in used_scaffolds]
        unused_substituents = [gid for gid in self.substituents.keys() if gid not in used_substituents]
        
        return unused_scaffolds, unused_substituents
    
    def suggest_common_pairs(self, scaffold_id: str = None, substituent_id: str = None) -> List[Dict[str, Any]]:
        """Suggest compound pairs that should typically exist"""
        suggestions = []
        
        # Common substituents that should pair with most scaffolds
        common_substituents = ['-Cl', '-Br', '-I', '-F', '-H', '-OH', '-NH2', '-B(OH)2']
        
        # Common scaffolds
        common_scaffolds = ['Ar', 'Alkyl', 'Bn', 'Allyl', 'Alkenyl', 'Alkynyl']
        
        # If specific scaffold provided, suggest with common substituents
        if scaffold_id:
            scaffolds_to_check = [scaffold_id] if scaffold_id in self.scaffolds else []
        else:
            scaffolds_to_check = common_scaffolds if not substituent_id else list(self.scaffolds.keys())
        
        # If specific substituent provided, suggest with common scaffolds
        if substituent_id:
            substituents_to_check = [substituent_id] if substituent_id in self.substituents else []
        else:
            substituents_to_check = common_substituents if not scaffold_id else list(self.substituents.keys())
        
        # Generate suggestions for pairs that don't exist
        for scaffold in scaffolds_to_check:
            if scaffold not in self.scaffolds:
                continue
                
            scaffold_info = self.scaffolds[scaffold]
            
            for substituent in substituents_to_check:
                if substituent not in self.substituents:
                    continue
                
                # Check if pair already exists (with or without _Subst suffix)
                if (scaffold, substituent) in self.existing_pairs:
                    continue
                if (scaffold, substituent + '_Subst') in self.existing_pairs:
                    continue
                
                substituent_info = self.substituents[substituent]
                
                # Determine template based on substituent type
                template = "single_bond"
                if substituent in ['OTf', 'OTs', 'OMs', 'OSO2R']:
                    template = "via_oxygen"
                
                # Generate compound suggestion
                # Avoid double dash if substituent starts with dash
                if substituent.startswith("-"):
                    compound_id = f"{scaffold}{substituent}"
                else:
                    compound_id = f"{scaffold}-{substituent}"
                compound_name = compound_id  # Simplified naming
                
                # Generate description
                scaffold_name = scaffold_info.get('name', scaffold)
                substituent_name = substituent_info.get('name', substituent)
                description = f"{scaffold_name} paired with {substituent_name} substituent."
                
                suggestion = {
                    "id": compound_id,
                    "template": template,
                    "A": scaffold,
                    "B": substituent,
                    "description": description,
                    "scaffold_kind": scaffold_info.get('kind'),
                    "substituent_kind": substituent_info.get('kind'),
                    "scaffold_tags": scaffold_info.get('tags', []),
                    "substituent_tags": substituent_info.get('tags', [])
                }
                
                suggestions.append(suggestion)
        
        return suggestions
    
    def print_suggestions(self, suggestions: List[Dict[str, Any]], title: str = "Suggested Compounds"):
        """Print suggestions in a readable format"""
        if not suggestions:
            print(f"\n✓ No new {title.lower()} needed - all common combinations exist!")
            return
        
        print(f"\n{'='*70}")
        print(f"{title} ({len(suggestions)} total)")
        print(f"{'='*70}\n")
        
        for i, suggestion in enumerate(suggestions, 1):
            print(f"{i}. {suggestion['id']}:")
            print(f"   Template: {suggestion['template']}")
            print(f"   A: {suggestion['A']} (scaffold, tags: {', '.join(suggestion['scaffold_tags'][:3])})")
            print(f"   B: {suggestion['B']} (substituent, tags: {', '.join(suggestion['substituent_tags'][:3])})")
            print(f"   Description: {suggestion['description']}")
            print()
    
    def generate_compound_json(self, suggestions: List[Dict[str, Any]]) -> str:
        """Generate JSON array for suggested compounds (ID auto-generated from A-B)"""
        compounds_json = []
        
        for suggestion in suggestions:
            compound = {
                "template": suggestion["template"],
                "A": suggestion["A"],
                "B": suggestion["B"],
                "description": suggestion["description"]
            }
            compounds_json.append(compound)
        
        return json.dumps(compounds_json, indent=2, ensure_ascii=False)
    
    def report_unused_groups(self):
        """Report groups that aren't used in any compounds"""
        unused_scaffolds, unused_substituents = self.find_unused_groups()
        
        if unused_scaffolds:
            print(f"\n{'='*70}")
            print(f"Unused Scaffolds ({len(unused_scaffolds)}):")
            print(f"{'='*70}")
            for scaffold_id in sorted(unused_scaffolds):
                scaffold = self.scaffolds[scaffold_id]
                print(f"  - {scaffold_id:20s} ({scaffold.get('name', 'N/A')})")
                print(f"    Description: {scaffold.get('description', 'N/A')[:60]}")
        
        if unused_substituents:
            print(f"\n{'='*70}")
            print(f"Unused Substituents ({len(unused_substituents)}):")
            print(f"{'='*70}")
            for substituent_id in sorted(unused_substituents):
                substituent = self.substituents[substituent_id]
                print(f"  - {substituent_id:20s} ({substituent.get('name', 'N/A')})")
                print(f"    Description: {substituent.get('description', 'N/A')[:60]}")
        
        return unused_scaffolds, unused_substituents
    
    def suggest_for_new_groups(self):
        """Main suggestion workflow"""
        print("\n" + "="*70)
        print("Organic Compound Suggestion Tool")
        print("="*70)
        
        if not self.load_data():
            return
        
        # Report unused groups
        unused_scaffolds, unused_substituents = self.report_unused_groups()
        
        # Suggest compounds for unused groups
        if unused_scaffolds or unused_substituents:
            print(f"\n{'='*70}")
            print("Suggested Compounds for Unused Groups:")
            print(f"{'='*70}")
            
            all_suggestions = []
            
            # Suggest for each unused scaffold
            for scaffold_id in unused_scaffolds:
                suggestions = self.suggest_common_pairs(scaffold_id=scaffold_id)
                all_suggestions.extend(suggestions)
            
            # Suggest for each unused substituent
            for substituent_id in unused_substituents:
                suggestions = self.suggest_common_pairs(substituent_id=substituent_id)
                all_suggestions.extend(suggestions)
            
            # Remove duplicates
            seen = set()
            unique_suggestions = []
            for s in all_suggestions:
                key = (s['A'], s['B'])
                if key not in seen:
                    seen.add(key)
                    unique_suggestions.append(s)
            
            self.print_suggestions(unique_suggestions, "Compounds for Unused Groups")
            
            return unique_suggestions
        else:
            print("\n✓ All groups are used in at least one compound!")
            return []


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Suggest new organic compounds based on available groups'
    )
    parser.add_argument(
        '--generate',
        action='store_true',
        help='Generate JSON output for suggested compounds'
    )
    parser.add_argument(
        '--scaffold',
        type=str,
        help='Generate suggestions for specific scaffold group ID'
    )
    parser.add_argument(
        '--substituent',
        type=str,
        help='Generate suggestions for specific substituent group ID'
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
    
    suggester = CompoundSuggester(taxonomy_dir)
    
    if args.scaffold or args.substituent:
        # Specific suggestions
        if not suggester.load_data():
            sys.exit(1)
        
        suggestions = suggester.suggest_common_pairs(
            scaffold_id=args.scaffold,
            substituent_id=args.substituent
        )
        
        if args.generate:
            print(suggester.generate_compound_json(suggestions))
        else:
            title = f"Suggestions for "
            if args.scaffold:
                title += f"scaffold '{args.scaffold}'"
            if args.substituent:
                if args.scaffold:
                    title += " with "
                title += f"substituent '{args.substituent}'"
            suggester.print_suggestions(suggestions, title)
    else:
        # Full workflow for unused groups
        suggestions = suggester.suggest_for_new_groups()
        
        if args.generate and suggestions:
            print(f"\n{'='*70}")
            print("JSON Output (copy to compounds file):")
            print(f"{'='*70}\n")
            print(suggester.generate_compound_json(suggestions))


if __name__ == "__main__":
    main()
