#!/usr/bin/env python3
"""Fix the corrupted Suzuki.jsonl file."""

# Read corrupted file
with open('data/reaction_dataset/Suzuki.jsonl', 'rb') as f:
    content = f.read()

print(f"Original file size: {len(content)} bytes")
print(f"First 50 bytes: {content[:50]}")

# Remove "updat" prefix
fixed = content.lstrip(b'updat')

print(f"\nFixed file size: {len(fixed)} bytes")
print(f"First 50 bytes: {fixed[:50]}")
print(f"Removed: {len(content) - len(fixed)} bytes")

# Write fixed content
with open('data/reaction_dataset/Suzuki.jsonl', 'wb') as f:
    f.write(fixed)

print(f"\n✓ File fixed! Wrote {len(fixed)} bytes")

# Verify it's valid JSON now
import json
with open('data/reaction_dataset/Suzuki.jsonl', 'r', encoding='utf-8') as f:
    first_line = f.readline()
    try:
        rxn = json.loads(first_line)
        print(f"✓ First line is valid JSON")
        print(f"  reaction_id: {rxn.get('reaction_id')}")
        print(f"  reaction_type: {rxn.get('reaction_type')}")
    except json.JSONDecodeError as e:
        print(f"✗ Still corrupted: {e}")
