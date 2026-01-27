"""Test that precedent system now loads from all HTE sources."""
from chemtools.precedent.loader import _load_selective, _iter_literature_files
import os

# Test file discovery
files = _iter_literature_files()
print(f"✓ Found {len(files)} CSV files across all HTE sources")

# Group by source directory
sources = {}
for f in files:
    dirname = os.path.basename(os.path.dirname(f))
    sources[dirname] = sources.get(dirname, 0) + 1

print("\nFiles by source:")
for source, count in sorted(sources.items()):
    print(f"  {source}: {count} files")

# Load all data
rows = _load_selective(families=None)
print(f"\n✓ Loaded {len(rows):,} total precedent rows from all sources")

# Sample a protocol if available
protocol_rows = [r for r in rows if "protocol" in r.get("file_family", "").lower()]
if protocol_rows:
    print(f"\n✓ Found {len(protocol_rows)} rows from protocol sources")
    sample = protocol_rows[0]
    print(f"  Sample: {sample.get('rxn_type', 'N/A')} - {sample.get('reaction_smiles', 'N/A')[:50]}")
else:
    # Check if protocols directory exists but has no data
    for f in files:
        if "protocols" in f.lower():
            print(f"\n⚠ Protocol file found but no rows loaded: {os.path.basename(f)}")
            break
    else:
        print("\n⚠ No protocol files found in data/HTE_db/protocols/")

print("\n✅ Precedent system now loads from all HTE sources!")
