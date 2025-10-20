# Protocol Validation Tool Summary

## Overview

Created a comprehensive CLI tool to validate protocol JSON files in the database. The tool ensures that each protocol's `reaction_smiles` matches its declared `reaction_SMARTS` patterns using RDKit substructure matching.

## Tool Location

- **File**: `chemtools/protocol/validate_protocols.py`
- **Module**: `chemtools.protocol.validate_protocols`

## Features

### 1. Comprehensive Validation
- ✅ Validates JSON structure
- ✅ Checks for required fields (`reaction`, `reaction_smiles`, `reaction_SMARTS`)
- ✅ Uses RDKit to verify SMARTS pattern matching
- ✅ Detects both errors and warnings

### 2. Flexible Usage
```bash
# Validate all protocols
python -m chemtools.protocol.validate_protocols

# Validate specific file
python -m chemtools.protocol.validate_protocols --file "Protocol.json"

# Show only errors (hide valid protocols)
python -m chemtools.protocol.validate_protocols --errors-only

# Verbose output with detailed information
python -m chemtools.protocol.validate_protocols --verbose

# Export validation report to JSON
python -m chemtools.protocol.validate_protocols --output report.json

# Exit with code 1 if any errors (for CI/CD)
python -m chemtools.protocol.validate_protocols --fail-on-error
```

### 3. Detailed Reporting

**Console Output:**
```
======================================================================
Protocol Validation Summary
======================================================================

Total protocols: 17
✅ Valid: 4
❌ Invalid: 13

======================================================================
Invalid Protocols:
======================================================================

❌ Aryl mesylate_Suzuki.json
   Reaction: CC(C)(C)c1ccc(OS(=O)(=O)C)cc1.B(O)(O)c2ccccc2>>CC(C)(C)c1ccc(-c2ccccc2)cc1
   SMARTS patterns: 1
   ERROR: Reaction SMILES does not match any of the 1 SMARTS pattern(s)
   Patterns:
     - O=S(OC)([c,C,n,o,s])=O.OB(O)[c,C,n,o,s]>>[c,C,n,o,s]
```

**JSON Report Format:**
```json
{
  "total": 17,
  "valid": 4,
  "invalid": 13,
  "results": [
    {
      "filename": "Aryl mesylate_Suzuki.json",
      "valid": false,
      "reaction_smiles": "...",
      "reaction_smarts": ["..."],
      "errors": ["Reaction SMILES does not match..."],
      "warnings": []
    }
  ]
}
```

## Validation Logic

The tool performs the following checks:

1. **JSON Integrity**: Validates that the file is well-formed JSON
2. **Structure Check**: Ensures `reaction` field exists
3. **Required Fields**: Checks for `reaction_smiles` and `reaction_SMARTS`
4. **SMARTS Matching**: Uses RDKit to verify:
   - Parse reaction SMILES
   - Parse each SMARTS pattern
   - Check if reaction matches at least one pattern
   - Validate both reactants and products

## Current Validation Results

From your database (17 protocols):
- ✅ **Valid**: 4 protocols
- ❌ **Invalid**: 13 protocols

### Common Issues Found

1. **Invalid JSON files** (2):
   - `Evano_2016_Cu_cyanation_alkenyl_iodides_stepA.json`
   - `Grubbs_RCM_Ferguson_2003.json`

2. **SMARTS pattern mismatch** (2):
   - `Aryl mesylate_Suzuki.json`
   - `Sonogashira-Coupling.json`

3. **RDKit SMARTS parsing errors** (9):
   - Issues with `[CH]`, `[CH2]`, `[CH3]` in SMARTS
   - RDKit pre-condition violation: `getNumImplicitHs()`

4. **Missing reaction field** (1):
   - `.protocol_index.json` (not a protocol file)

## Recommended Fixes

### Issue: SMARTS Pattern with `[CH]` Causes RDKit Error

**Problem**: SMARTS like `Br[CH]` cause RDKit to fail with:
```
Pre-condition Violation: getNumImplicitHs() called without preceding call to calcImplicitValence()
```

**Solution**: Use more general patterns:
- Replace `[CH]` with `[C;H1]` or `[c,C,n,o,s]`
- Replace `[CH2]` with `[C;H2]` or `[C]`
- Replace `[CH3]` with `[C;H3]` or `[C]`

**Example Fix for `Aryl mesylate_Suzuki.json`:**
```json
// Current (doesn't match):
"reaction_SMARTS": ["O=S(OC)([c,C,n,o,s])=O.OB(O)[c,C,n,o,s]>>[c,C,n,o,s]"]

// Fixed (matches reaction):
"reaction_SMARTS": ["[c,C,n,o,s]OS(=O)(=O)C.OB(O)[c,C,n,o,s]>>[c,C,n,o,s][c,C,n,o,s]"]
```

### Issue: Sonogashira Pattern Not Matching

**Problem**: Pattern expects terminal alkyne `[H]` but product has internal alkyne
```json
"reaction_SMARTS": ["[c,C,n,o,s]I.[!#1]C#C[H]>>[c,C,n,o,s]C#C[!#1]"]
```

**Solution**: Make pattern more flexible:
```json
"reaction_SMARTS": ["[c,C,n,o,s]I.C#C>>[c,C,n,o,s]C#C"]
```

## Integration with Workflow

### Before Building Index
```bash
# 1. Validate all protocols
python -m chemtools.protocol.validate_protocols --errors-only

# 2. Fix any errors in protocol JSON files

# 3. Validate again to confirm
python -m chemtools.protocol.validate_protocols

# 4. Build index
python -m chemtools.protocol.cli build --force
```

### CI/CD Integration
```yaml
# In your CI/CD pipeline
- name: Validate Protocols
  run: python -m chemtools.protocol.validate_protocols --fail-on-error
  
- name: Build Index
  run: python -m chemtools.protocol.cli build --force
```

## Code Structure

```python
# Main validation function
def validate_protocol(filepath: Path) -> ValidationResult:
    """
    Validates:
    1. JSON can be loaded
    2. Required fields exist
    3. Reaction SMILES matches SMARTS patterns
    """

# SMARTS matching using RDKit
def match_reaction_smarts(
    reaction_smiles: str, 
    smarts_patterns: List[str]
) -> Tuple[bool, List[str]]:
    """
    Uses RDKit to:
    1. Parse reaction SMILES
    2. Parse SMARTS patterns
    3. Check substructure matching
    4. Return match status and errors
    """
```

## Benefits

1. **Data Quality**: Ensures protocols are valid before indexing
2. **Early Detection**: Catches SMARTS/SMILES mismatches early
3. **CI/CD Ready**: Can be integrated into automated workflows
4. **Detailed Reports**: Helps identify and fix specific issues
5. **Flexible**: Can validate all or specific protocols

## Next Steps

1. Fix the 13 invalid protocols in the database
2. Add validation to the index building process
3. Create a pre-commit hook to validate new protocols
4. Consider adding more validation rules (e.g., required chemical roles, source metadata)

---

**Created**: October 20, 2025
**Tool**: `chemtools/protocol/validate_protocols.py`
**Documentation**: Updated in `chemtools/protocol/readme.md`
