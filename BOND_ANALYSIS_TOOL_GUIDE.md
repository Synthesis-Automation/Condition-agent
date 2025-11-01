# Bond Analysis Tool Integration - Quick Guide

## Overview

The `analyze_bond_changes_tool` has been successfully integrated into the ChemTools Assistant, enabling AI-powered analysis of bond breaking and formation in chemical reactions.

## Features

### Multi-Method Hybrid Approach
- **Manual Mapping**: Uses existing atom mapping if present (100% confidence)
- **RXNMapper**: ML-based automatic atom mapping (70-90% confidence)  
- **MCS**: Maximum Common Substructure fallback (50-80% confidence)

### Capabilities
✅ Detects broken bonds (including leaving groups like Br, I, B)
✅ Identifies formed bonds (new C-C, C-N, etc.)
✅ Tracks leaving groups and joining groups accurately
✅ Provides confidence scores and method agreement
✅ Validates results across multiple methods

## Usage with Assistant

### Example Questions You Can Ask:

1. **Basic Bond Analysis**
   ```
   "What bonds break and form in this Suzuki coupling: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
   ```

2. **Detailed Analysis**
   ```
   "Analyze the bond changes in CC(=O)O.CCO>>CC(=O)OCC and tell me which method has higher confidence"
   ```

3. **Reaction Mechanism Understanding**
   ```
   "For the reaction Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1, identify all leaving groups and which bonds form"
   ```

4. **Comparative Analysis**
   ```
   "Compare the bond breaking patterns between Suzuki coupling and Buchwald-Hartwig amination"
   ```

## Tool Parameters

```python
analyze_bond_changes_tool(
    reaction_smiles: str,      # Reaction SMILES (reactants>>products)
    use_hybrid: bool = True    # Use hybrid approach (recommended)
)
```

## Output Format

```json
{
  "success": true,
  "method": "hybrid",
  "combined_confidence": 0.75,
  "broken_bonds": [
    [5, "Br (leaving group)"],
    [4, "B (leaving group)"]
  ],
  "formed_bonds": [
    [4, 5]
  ],
  "leaving_groups": [
    [5, "Br", "SINGLE"],
    [4, "B", "SINGLE"]
  ],
  "agreement": {
    "rxnmapper_vs_mcs": true
  },
  "validation": "Both methods agree",
  "interpretation": "Substitution/coupling: 2 bond(s) break, 1 bond(s) form"
}
```

## Key Improvements

### Before
- Only RXNMapper
- Missed leaving groups (Br, I not detected)
- No validation
- Lower confidence

### After  
- Hybrid approach (Manual + RXNMapper + MCS)
- ✅ Detects leaving groups accurately
- ✅ Cross-validation between methods
- ✅ Higher combined confidence when methods agree
- ✅ Warnings when methods disagree

## Testing

Run the test script to verify integration:

```bash
python test_bond_tool.py
```

Or use the assistant CLI:

```bash
python -m chem_assistant.chemtools_cli
```

Then ask: "Analyze bonds in Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"

## Technical Details

- **Location**: `chem_assistant/chemtools_wrapper.py`
- **Tool Name**: `analyze_bond_changes_tool`
- **Backend**: `chemtools._atom_mapping` module
- **Dependencies**: RXNMapper (optional but recommended), RDKit

## Installation

If RXNMapper is not installed:

```bash
pip install rxnmapper
```

## Notes

- The hybrid approach is **recommended** for best accuracy
- Manual atom-mapped SMILES get 100% confidence (ground truth)
- When methods disagree, the tool flags it for review
- Leaving groups are now correctly identified (fixes previous bug)
