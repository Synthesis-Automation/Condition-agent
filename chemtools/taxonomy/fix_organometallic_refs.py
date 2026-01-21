#!/usr/bin/env python3
"""Quick fix for organometallic group references"""
import json
from pathlib import Path

# Get the absolute path
taxonomy_dir = Path(__file__).parent / "data"
compounds_file = taxonomy_dir / "organic_compounds.v1.3.json"

# Load
with open(compounds_file, 'r', encoding='utf-8') as f:
    data = json.load(f)

# Fix references
count = 0
for compound in data['compounds']:
    b_value = compound.get('B')
    if b_value in ['Sn', 'Zn', 'Mg', 'Si']:
        compound['B'] = b_value + '*'
        count += 1
        print(f"  Fixed: {compound['id']} - B: '{b_value}' → '{b_value}*'")

# Save
with open(compounds_file, 'w', encoding='utf-8') as f:
    json.dump(data, f, indent=2, ensure_ascii=False)

print(f"\n✓ Fixed {count} organometallic group references in {compounds_file.name}")
