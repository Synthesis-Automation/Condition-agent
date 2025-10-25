# Rule-Based System Updates for Standardized Rule Format

## Overview
Updated the rule-based recommendation system and output formatter to support the new standardized rule file format while maintaining backward compatibility with legacy formats.

## Changes Made

### 1. **Updated MatchResult Type** (`chemtools/rule_scdb_matcher/types.py`)
- Added `reagents` field to `MatchResult` dataclass
- Updated `to_json_dict()` to include reagents when present
- Maintains backward compatibility - reagents field is optional

```python
@dataclass(slots=True)
class MatchResult:
    ...
    reagents: Any = None  # New field for standardized reagents array
```

### 2. **Updated Matcher** (`chemtools/rule_scdb_matcher/matcher.py`)
- Modified `_match_scheme()` to extract reagents from entry's raw data
- Passes reagents to MatchResult constructor
- Reagents are extracted from `selected_entry.raw.get("reagents")`

### 3. **Updated Output Formatter** (`chemtools/output_formatter.py`)
- Completely rewrote `_convert_rule_match_to_recommendations()` function
- **New format support**: Extracts reagents from `reagents` array
- **Legacy format support**: Falls back to old keys (`pd_source`, `ligand`, `base`, `solvent`)
- **Role mapping**: Maps standardized roles to output roles
  - `metal_source` â†?`metal_precursor`
  - `ligand` â†?`ligand`
  - `base` â†?`base`
  - `solvent` â†?`solvent`
  - `additive` â†?`additive`
  - `boron_partner` â†?`partner`

### 4. **Fixed Rule Database** (`data/rule_db/Suzuki_db.json`)
- Added missing SPhos ligand to SCDB-SUZ-ARBRI-GENERAL-SPhos entry
- Ensured all reagents mentioned in `token_signature` are in `reagents` array

### 5. **Created Utility Script** (`scripts/fix_rule_reagents.py`)
- Automated script to fix missing ligands in rule databases
- Scans `token_signature` and adds missing reagents to `reagents` array
- Supports Suzuki and C-N coupling databases
- Predefined ligand mappings with appropriate amounts

## New Rule File Format

### Structure
```json
{
  "entries": [
    {
      "id": "SCDB-SUZ-ARBRI-GENERAL-SPhos",
      "reaction_type": "Suzuki_Miyaura",
      "name": "Aryl iodides/bromides + aryl boron (SPhos set)",
      "rxn_smiles_min": "[c:1]-[I,Br:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
      "token_signature": [...],
      "conditions": {
        "temperature_C": [45.0, 60.0],
        "time_h": [6.0, 10.0]
      },
      "reagents": [
        {
          "name": "Pd2(dba)3",
          "role": "metal_source",
          "amount": "0.5â€?.5%"
        },
        {
          "name": "SPhos",
          "role": "ligand",
          "amount": "1.0â€?.0%"
        },
        {
          "name": "K2CO3",
          "role": "base",
          "amount": "200%"
        },
        {
          "name": "THF/H2O (4:1)",
          "role": "solvent",
          "amount": ""
        }
      ],
      "env": {...},
      "evidence": {...},
      "notes": [...],
      "priority": 46
    }
  ]
}
```

### Key Changes from Legacy Format
1. **`reagents` array**: New top-level field with structured reagent information
2. **Role-based organization**: Each reagent has explicit `role`, `name`, and `amount`
3. **Consolidated conditions**: `temperature_C` and `time_h` in `conditions` object
4. **Removed old keys**: No longer using `pd_source`, `ligand`, `base`, `solvent` as direct keys

## Output Format

### Standardized Recommendation Output
```json
{
  "recommended_conditions": [
    {
      "rank": 1,
      "chemicals": [
        {
          "smiles": "Ic1ccncc1",
          "role": "electrophile",
          "name": null,
          ...
        },
        {
          "smiles": "c1ccc(B(O)O)cc1",
          "role": "nucleophile",
          ...
        },
        {
          "name": "Pd2(dba)3",
          "role": "metal_precursor",
          ...
        },
        {
          "name": "SPhos",
          "role": "ligand",
          ...
        },
        {
          "name": "K2CO3",
          "role": "base",
          ...
        },
        {
          "name": "THF/H2O (4:1)",
          "role": "solvent",
          ...
        }
      ],
      "conditions": {
        "temperature": [45.0, 60.0],
        "time": [6.0, 10.0],
        "atmosphere": null
      },
      "source": {
        "match_type": "scheme",
        "entry_id": "SCDB-SUZ-ARBRI-GENERAL-SPhos",
        "entry_name": "Aryl iodides/bromides + aryl boron (SPhos set)",
        "priority": 46
      }
    }
  ]
}
```

## Backward Compatibility

The system maintains full backward compatibility:
- **Tries new format first**: Looks for `reagents` array
- **Falls back to legacy**: Uses old keys if `reagents` not found
- **No breaking changes**: Existing rule files without `reagents` array continue to work

## Testing

### Test Results
âœ?Rule-based matching: Working with new format
âœ?Reagent extraction: All reagents (metal, ligand, base, solvent) extracted correctly
âœ?Output formatting: Proper role mapping and chemical structure
âœ?Backward compatibility: Legacy format still supported
âœ?CLI integration: Local recommendation CLI works end-to-end

### Test Command
```powershell
echo "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1`n2" | python scripts/local_recommendation_cli.py
```

### Expected Output
```
Rule-based family: Suzuki_Miyaura
Top match (rank 1):
  Entry: Aryl iodides/bromides + aryl boron (SPhos set)
  Match type: scheme
  Catalyst: Pd2(dba)3
  Temperature: [45.0, 60.0]Â°C
```

## Next Steps

### Recommended Actions
1. **Update remaining rule databases**: Apply same reagent format to all condition DBs
2. **Add validation**: Create schema validator for new format
3. **Update documentation**: Document new format in API guides
4. **Add tests**: Create unit tests for new format parsing

### Optional Enhancements
1. **Amount parsing**: Extract numeric amounts from amount strings
2. **Reagent enrichment**: Look up CAS numbers and SMILES for reagents
3. **Equivalents normalization**: Convert percentages to equivalents
4. **Validation warnings**: Warn if token_signature doesn't match reagents array

## Files Modified

1. `chemtools/rule_scdb_matcher/types.py` - Added reagents field to MatchResult
2. `chemtools/rule_scdb_matcher/matcher.py` - Extract and pass reagents
3. `chemtools/output_formatter.py` - Support new and legacy formats
4. `data/rule_db/Suzuki_db.json` - Added missing SPhos ligand
5. `scripts/fix_rule_reagents.py` - New utility script (created)
6. `scripts/recommendation_cli_utils.py` - Already updated for single output key

## Summary

The rule-based recommendation system has been successfully updated to support the new standardized rule file format with `reagents` arrays. All reagents (metal sources, ligands, bases, solvents) are now properly extracted and formatted in the output. The system maintains full backward compatibility with legacy formats, ensuring a smooth transition.
