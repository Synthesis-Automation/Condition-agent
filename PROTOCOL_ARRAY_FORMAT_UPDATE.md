# Protocol Array Format Update

## Summary

Updated the protocol-based recommendation system to support the new JSON format where each protocol file can contain multiple protocols as an array.

## Changes Made

### 1. **Core Indexer** (`chemtools/protocol/indexer.py`)
- Modified `_process_protocol_file()` to return lists of records and fingerprints
- Added `_process_single_protocol()` to handle individual protocol processing
- Updated protocol identification scheme:
  - Single protocol: `filename.json` → ID: `filename`
  - Multiple protocols: `filename.json` → IDs: `filename[0]`, `filename[1]`, etc.
- Updated `build_index()` to handle multiple records per file

### 2. **Protocol Recommender** (`chemtools/protocol/recommend.py`)
- Updated `get_protocol_details()` to parse protocol identifiers with array indices
- Handles both formats:
  - Legacy: `Suzuki_Cu_C(sp3)-C(sp2).json` (single protocol)
  - New: `Suzuki_protocols[0]` (protocol at index 0 in array)

### 3. **Protocol Matcher** (`chemtools/protocol/matcher.py`)
- Updated `_extract_metadata()` to return list of metadata objects
- Updated `get_protocol()` to handle array-indexed filenames
- Maintains backward compatibility with legacy single-protocol format

### 4. **Batch Update SMARTS** (`chemtools/protocol/batch_update_protocol_smarts.py`)
- Modified `process_protocol_file()` to iterate over protocol arrays
- Preserves original file format (array or single object)
- Reports aggregate results for multi-protocol files

### 5. **Atom Mapping Tool** (`chemtools/protocol/add_atom_mapping.py`)
- Updated `process_protocol_file()` to handle protocol arrays
- Generates mappings for each protocol in array
- Preserves original file format when writing updates

### 6. **Protocol Validator** (`chemtools/protocol/validate_protocols.py`)
- Modified `load_protocol()` to always return list of protocols
- Updated `validate_protocol()` to return list of validation results
- Reports validation status for each protocol in multi-protocol files

## File Format Support

### Legacy Format (Single Protocol)
```json
{
  "source": { ... },
  "reaction": { ... },
  "reaction_setup": [ ... ]
}
```

### New Format (Protocol Array)
```json
[
  {
    "source": { ... },
    "reaction": { ... },
    "reaction_setup": [ ... ]
  },
  {
    "source": { ... },
    "reaction": { ... },
    "reaction_setup": [ ... ]
  }
]
```

Both formats are fully supported. Single protocols are automatically wrapped in an array internally for consistent processing.

## Protocol Identification

- **Single protocol file**: Uses base filename as ID
  - Example: `Ni_Suzuki_ArylHalide+BoronicAcid_tAmOH.json` → ID: `Ni_Suzuki_ArylHalide+BoronicAcid_tAmOH`

- **Multi-protocol file**: Uses filename with array index
  - Example: `Suzuki_protocols.json` with 4 protocols → IDs: `Suzuki_protocols[0]`, `Suzuki_protocols[1]`, `Suzuki_protocols[2]`, `Suzuki_protocols[3]`

## Next Steps

### 1. Rebuild Protocol Index
After adding new protocol files, rebuild the index:

```powershell
python -m chemtools.protocol.cli build --compute-drfp
```

This will:
- Scan all JSON files in `data/protocol_db/`
- Create index entries for each protocol (including array elements)
- Compute DRFP fingerprints for similarity search
- Save index to `.protocol_index.json` and `.protocol_drfp.npz`

### 2. Validate Protocols
Check for SMARTS pattern matching issues:

```powershell
python -m chemtools.protocol.validate_protocols
```

### 3. Update SMARTS Patterns (if needed)
For protocols with validation errors, regenerate SMARTS patterns:

```powershell
python -m chemtools.protocol.batch_update_protocol_smarts
```

## Current Validation Status

Running validation shows:
- **Total protocols**: 20 (across 17 files)
- **Valid**: 17
- **Invalid**: 3 (SMARTS pattern mismatches in `Suzuki_protocols.json[1]`, `[2]`, `[3]`)

The validation errors indicate that some SMARTS patterns need to be updated to match the actual reaction SMILES. This is expected when using manually-defined patterns.

## Backward Compatibility

All changes maintain full backward compatibility:
- Existing code using single-protocol files continues to work
- Index loading handles both old and new formats
- Protocol lookup supports both filename styles
- All tools auto-detect and adapt to the format

## Testing

To test the updated system:

```powershell
# Rebuild index with new multi-protocol files
python -m chemtools.protocol.cli build --compute-drfp

# Test recommendation with a Suzuki reaction
python -m chemtools.protocol.cli recommend "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --k 5

# Validate all protocols
python -m chemtools.protocol.validate_protocols
```
