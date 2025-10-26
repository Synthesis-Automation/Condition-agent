#!/usr/bin/env python3
"""Check which reactions are being filtered out."""

import sys
sys.path.insert(0, 'tests')
from sample_reactions import SAMPLE_REACTIONS, BUCHWALD_HARTWIG_REACTIONS

def extract_smiles_from_entry(entry: str) -> str:
    """Extract SMILES from a sample reaction entry - same logic as test script."""
    if ">>>" in entry or entry.startswith("Select"):
        return None
    
    # Find the last opening parenthesis (description starts there)
    # This avoids splitting SMILES like B(O)O
    last_paren = entry.rfind(" (")
    if last_paren == -1:
        # No description in parentheses - check if it's still a valid reaction
        return entry.strip() if ">>" in entry else None
    
    # Extract SMILES before the description
    smiles = entry[:last_paren].strip()
    return smiles if ">>" in smiles else None

# Count SAMPLE_REACTIONS
valid_sample = 0
filtered_sample = 0
for entry in SAMPLE_REACTIONS:
    if extract_smiles_from_entry(entry) is not None:
        valid_sample += 1
    elif ">>" in entry:
        filtered_sample += 1
        if filtered_sample <= 5:
            print(f"SAMPLE filtered: {entry[:100]}")

# Count BUCHWALD_HARTWIG_REACTIONS
valid_bh = 0
filtered_bh = 0
for entry in BUCHWALD_HARTWIG_REACTIONS:
    if extract_smiles_from_entry(entry) is not None:
        valid_bh += 1
    elif ">>" in entry:
        filtered_bh += 1
        if filtered_bh <= 5:
            print(f"B-H filtered: {entry[:100]}")

print(f"\nSAMPLE_REACTIONS:")
print(f"  Total entries: {len(SAMPLE_REACTIONS)}")
print(f"  Valid (tested): {valid_sample}")
print(f"  Filtered out: {filtered_sample}")

print(f"\nBUCHWALD_HARTWIG_REACTIONS:")
print(f"  Total entries: {len(BUCHWALD_HARTWIG_REACTIONS)}")
print(f"  Valid (tested): {valid_bh}")
print(f"  Filtered out: {filtered_bh}")

print(f"\nTOTAL TESTED: {valid_sample + valid_bh}")
print(f"TOTAL FILTERED: {filtered_sample + filtered_bh}")
print(f"GRAND TOTAL: {len(SAMPLE_REACTIONS) + len(BUCHWALD_HARTWIG_REACTIONS)}")
