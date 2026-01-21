#!/usr/bin/env python3
"""
Further simplify taxonomy by removing 'aliases' and 'tags' fields.
Move important information to 'description' if needed.
"""

import json
from pathlib import Path
from typing import List, Dict, Any


def simplify_groups(groups_file: Path):
    """Remove tags from groups (keep description as-is)"""
    with open(groups_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    groups = data.get('groups', [])
    removed_tags = 0
    
    for group in groups:
        if 'tags' in group:
            del group['tags']
            removed_tags += 1
    
    with open(groups_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write('\n')
    
    print(f"✓ Updated {groups_file.name}")
    print(f"  - Removed 'tags' from {removed_tags} groups")


def simplify_compounds(compounds_file: Path):
    """Remove aliases from compounds (keep description as-is)"""
    with open(compounds_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    compounds = data.get('compounds', [])
    removed_aliases = 0
    
    for compound in compounds:
        if 'aliases' in compound:
            del compound['aliases']
            removed_aliases += 1
    
    with open(compounds_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write('\n')
    
    print(f"✓ Updated {compounds_file.name}")
    print(f"  - Removed 'aliases' from {removed_aliases} compounds")


def main():
    print("=" * 70)
    print("Removing aliases and tags for maximum simplicity")
    print("=" * 70)
    print()
    
    taxonomy_dir = Path("chemtools/taxonomy/data")
    groups_file = taxonomy_dir / "organic_groups.v1.3.json"
    compounds_file = taxonomy_dir / "organic_compounds.v1.3.json"
    
    simplify_groups(groups_file)
    print()
    simplify_compounds(compounds_file)
    
    print("\n" + "=" * 70)
    print("✓ Maximum Simplicity Achieved!")
    print("=" * 70)
    print("\nRemaining fields per group:")
    print("  • id: Group identifier")
    print("  • kind: scaffold | substituent")
    print("  • priority: Matching priority (optional)")
    print("  • smarts: SMARTS pattern")
    print("  • description: Human-readable description")
    print("\nRemaining fields per compound:")
    print("  • A: Scaffold group ID")
    print("  • B: Substituent group ID")
    print("  • template: Bonding pattern")
    print("  • description: Human-readable description")


if __name__ == "__main__":
    main()
