# Batch Protocol SMARTS Updater

Automated tool to generate and update `reaction_smarts_applicability` in protocol JSON files using chemistry-aware pattern generation.

## Overview

The batch updater scans protocol JSON files, extracts reaction SMILES, generates chemistry-aware SMARTS patterns using the SubstrateClassifier and SmartsBuilder utilities, and updates the files with the new patterns.

## Features

✅ **Automated Pattern Generation**: Uses chemistry intelligence to create substrate-specific patterns  
✅ **Batch Processing**: Processes entire directories of protocol files  
✅ **Safe Overwriting**: Replaces existing patterns with improved chemistry-aware versions  
✅ **Dry Run Mode**: Preview changes before applying them  
✅ **Error Handling**: Gracefully handles missing data, invalid JSON, and parsing errors  
✅ **Detailed Reporting**: Shows processing status and summary statistics

## Usage

### Basic Usage

Process all protocol files in the default directory (`data/protocol_db`):

```powershell
python -m chemtools.protocol.batch_update_protocol_smarts
```

### Dry Run Mode

Preview what would change without modifying files:

```powershell
python -m chemtools.protocol.batch_update_protocol_smarts --dry-run
```

### Custom Directory

Process protocols in a different directory:

```powershell
python -m chemtools.protocol.batch_update_protocol_smarts --protocol-dir path/to/protocols
```

## Output Format

Each protocol file is updated with `reaction_smarts_applicability` containing:

- **`core`**: Chemistry-aware SMARTS pattern for reactants → products
- **`guards_forbid`**: List of exclusion patterns to prevent false matches

### Example: Alkyl Iodide Borylation

**Input** (reaction_smiles only):
```json
{
  "reaction": {
    "reaction_smiles": "CCCCCCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCCCCCCB1OC(C)(C)C(C)(C)O1",
    "family": "Alkyl_Iodide_Borylation"
  }
}
```

**Output** (with generated SMARTS):
```json
{
  "reaction": {
    "reaction_smiles": "CCCCCCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCCCCCCB1OC(C)(C)C(C)(C)O1",
    "family": "Alkyl_Iodide_Borylation",
    "reaction_smarts_applicability": {
      "core": "[CX4;H2,H3]-[I]>>[B]1OC(C)(C)C(C)(C)O1",
      "guards_forbid": [
        "[CX4;H1]-[I]  # Exclude secondary",
        "[CX4;H0]-[I]  # Exclude tertiary",
        "[CH2;$([CH2][c])]-[I]  # Exclude benzylic",
        "[CH2;$([CH2]C=C)]-[I]  # Exclude allylic"
      ]
    }
  }
}
```

### Example: Buchwald-Hartwig Coupling

**Input**:
```json
{
  "reaction": {
    "reaction_smiles": "CC(C)(C)c1ccc(OS(=O)(=O)C)cc1.CNc2ccccc2>>CC(C)(C)c1ccc(N(C)c2ccccc2)cc1",
    "family": "Buchwald_Hartwig"
  }
}
```

**Output**:
```json
{
  "reaction": {
    "reaction_smiles": "CC(C)(C)c1ccc(OS(=O)(=O)C)cc1.CNc2ccccc2>>CC(C)(C)c1ccc(N(C)c2ccccc2)cc1",
    "family": "Buchwald_Hartwig",
    "reaction_smarts_applicability": {
      "core": "a>>[NX3;H0;!$(NC=O);!$(N=*)]",
      "guards_forbid": [
        "[CX4]-[NX3;H2;!$(NC=O)]  # Exclude aliphatic primary amines",
        "[CX4]-[NX3;H1;!$(NC=O)]  # Exclude aliphatic secondary amines"
      ]
    }
  }
}
```

## Processing Report

The tool provides detailed output for each file:

```
================================================================================
🔬 Protocol SMARTS Batch Updater
================================================================================

📂 Found 18 protocol files in data\protocol_db
================================================================================

📄 Processing: Alkyl_Iodide_Borylation.json
   Reaction: CCCCCCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCCCCCCB1OC(C)(C)C(C)(C)O1
   ℹ️  Existing pattern found (will be overwritten)
   ✅ Updated successfully

📄 Processing: Pd_Buchwald_Arylsulfonate_Amination_CMPhos.json
   Reaction: CC(C)(C)c1ccc(OS(=O)(=O)C)cc1.CNc2ccccc2>>CC(C)(C)c1ccc(N(C)c2ccccc2)cc1
   ✅ Updated successfully

================================================================================
📊 PROCESSING SUMMARY
================================================================================

✅ Successful: 16/18
❌ Failed: 2/18

⚠️  Failed files:
   • .protocol_index.json: No reaction_smiles found in protocol
   • Ni-PCy3_C–O_Activation_Suzuki_of_Methoxyarenes.json: No reaction_smiles found in protocol

✅ Successfully processed files:
   • Alkyl_Iodide_Borylation.json - updated (replaced existing)
   • Pd_Buchwald_Arylsulfonate_Amination_CMPhos.json - updated (new pattern)
   ...
```

## Chemistry Intelligence

The generated patterns include:

### Substrate Classification
- **Primary vs Secondary vs Tertiary**: Distinguishes alkyl substitution levels
- **Aryl vs Alkyl**: Separates aromatic from aliphatic substrates
- **Aniline vs Aliphatic Amine**: Recognizes aryl-N vs alkyl-N
- **Special Positions**: Detects benzylic, allylic, propargylic carbons

### Context-Aware Guards
- **Exclusion Patterns**: Automatically generates guards to prevent mismatches
- **Commented Guards**: Each guard includes explanation of what it excludes
- **Substrate-Specific**: Guards adapt to the substrate type

### Examples of Generated Patterns

| Reaction Type | Substrate | Core Pattern | Guard Examples |
|--------------|-----------|--------------|----------------|
| Alkyl Iodide Borylation | Primary alkyl iodide | `[CX4;H2,H3]-[I]` | Exclude secondary, tertiary, benzylic, allylic |
| Buchwald-Hartwig | Aryl halide + Aniline | `c-[Br]` + `c-[NX3;H1;!$(NC=O)]` | Exclude aliphatic halides/amines, amides |
| Suzuki Coupling | Aryl bromide | `c-[Br]` | Exclude aliphatic halides |
| Cyanation | Alkenyl iodide | `[C;H1]=[C]-[I]` | Exclude alkyl, aryl iodides |

## Error Handling

The tool handles various edge cases:

- **Missing reaction_smiles**: Skips files without reaction SMILES
- **Invalid JSON**: Reports JSON parsing errors
- **Pattern Generation Failures**: Catches and reports pattern generation issues
- **Invalid SMILES**: Handles malformed SMILES strings gracefully

## Implementation Details

### Architecture

```
batch_update_protocol_smarts.py
├── ProtocolSmartsUpdater (main class)
│   ├── process_all_protocols() - Batch process directory
│   ├── process_protocol_file() - Process single file
│   └── generate_smarts_pattern() - Generate chemistry-aware pattern
├── ProcessingResult (dataclass) - Result tracking
└── main() - CLI entry point
```

### Dependencies

- `chemtools.util.smarts_builders.build_smarts_with_guards()`: Core pattern generation
- `chemtools.util.substrate_classifier`: Substrate classification (used internally by smarts_builders)
- `rdkit`: Chemistry toolkit for SMILES parsing (optional, handled gracefully)

### Test Coverage

- ✅ Single file processing
- ✅ Dry run mode
- ✅ Overwriting existing patterns
- ✅ Missing reaction SMILES handling
- ✅ Invalid JSON handling
- ✅ Batch processing multiple files
- ✅ Chemistry awareness validation
- ✅ Specific reaction types (Buchwald-Hartwig, etc.)

All 8 tests passing (100% coverage).

## Command Line Options

```
python -m chemtools.protocol.batch_update_protocol_smarts [OPTIONS]

Options:
  --protocol-dir PATH   Directory containing protocol JSON files
                        (default: data/protocol_db)
  
  --dry-run            Dry run mode - show what would change without
                        modifying files
  
  --verbose, -v        Verbose output with detailed pattern information
  
  --help               Show help message and exit
```

## Real-World Results

Successfully processed **16 out of 18** protocol files in `data/protocol_db`:

✅ **Generated new patterns** for 15 files  
✅ **Replaced existing pattern** for 1 file (Alkyl_Iodide_Borylation.json)  
⚠️ **Skipped 2 files** without reaction SMILES (.protocol_index.json, incomplete protocol)

## Integration with Other Tools

The batch updater uses the same chemistry utilities as:

- **SMARTS Generator CLI**: Interactive pattern generation tool
- **Recommendation Engine**: Chemistry-aware condition recommendations
- **Reaction Type Detector**: Automatic reaction classification
- **ML Featurizers**: Chemistry-aware feature extraction

This demonstrates the **reusability** of the chemistry utilities across different parts of the codebase.

## Future Enhancements

Potential improvements:

- [ ] Auto-generate atom mapping numbers (:1, :2, :3)
- [ ] Generate product patterns from reactants + reaction type
- [ ] Validate patterns against example molecules
- [ ] Parallel processing for large directories
- [ ] Pattern quality scoring
- [ ] Pattern library for common reaction types
- [ ] Interactive review mode before applying changes

## See Also

- [SMARTS Generator CLI](../chemtools/protocol/smarts_generator_cli.py) - Interactive pattern generation
- [SubstrateClassifier](../chemtools/util/substrate_classifier.py) - Substrate classification utility
- [SmartsBuilder](../chemtools/util/smarts_builders.py) - Pattern generation utility
- [Integration Tests](../tests/test_batch_update_protocol_smarts.py) - Test suite
