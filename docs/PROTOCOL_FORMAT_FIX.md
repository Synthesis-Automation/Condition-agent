# Protocol Recommendation Standard Format Update

## Issue Fixed

The protocol-based recommendation system was not following the standard JSON output format used by other recommendation types (ML-based, Rule-based). 

### Previous Format (Incorrect)
```json
{
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 0.3578,
      "protocol": {  // ❌ Metadata was at top level
        "filename": "...",
        "title": "...",
        ...
      },
      "similarity": 0.3578,
      "source": "protocol_database",
      "conditions": { ... }  // ⚠️ Only simplified conditions
    }
  ]
}
```

### New Format (Correct) ✅
```json
{
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 1.0,
      "chemicals": [  // ✅ Full chemicals list from protocol
        {
          "name": "Cyclopentanone",
          "cas": "120-92-3",
          "smiles": "O=C1CCCC1",
          "amount": {
            "weight_g": 1.68,
            "mmol": 20.0,
            "equivalents": 1.0
          },
          "role": "starting_material"
        },
        ...
      ],
      "conditions": {  // ✅ Extracted conditions
        "catalyst": "Pd(OAc)2",
        "ligand": "P(o-tol)3",
        "base": "NaOAc",
        "solvent": "1,4-Dioxane",
        "temperature_C": 130,
        "time_h": 24.0,
        "atmosphere": "Ar; sealed heavy-walled pressure vessel"
      },
      "source": "protocol_database",
      "similarity": 1.0,
      "protocol_metadata": {  // ✅ Metadata in nested object
        "filename": "alpha_arylation_Pd_enamine_Dong_v100p0099.json",
        "title": "α-Arylation of Cyclopentanones...",
        "journal": "Organic Syntheses",
        "year": 2023,
        "doi": "10.15227/orgsyn.100.0099",
        "url": "https://www.orgsyn.org/demo.aspx?prep=v100p0099",
        "reaction_smiles": "O=C1CCCC1.Brc2ccc(C(C)=O)cc2>>...",
        "reaction_smarts": ["O=C1CCCC1.Br[a]>>[a]C2C(CCC2)=O"],
        "reaction_family": "Pd/Enamine_α-Arylation_C(sp3)-C(sp2)",
        "tags": ["palladium", "enamine", ...],
        "notes": "Step A only..."
      }
    }
  ]
}
```

## Changes Made

### File: `chemtools/protocol/recommend.py`

1. **Updated `_format_standard_output()` method**:
   - Now loads the full protocol JSON to extract the chemicals list
   - Extracts chemicals from `reaction_setup[0].chemicals`
   - Extracts conditions using the existing `extract_conditions()` method
   - Moves protocol metadata to a nested `protocol_metadata` object
   - Includes `chemicals` and `conditions` at the top level (standard format)

2. **Simplified `recommend_with_details()` method**:
   - No longer needs to add conditions separately
   - Conditions are already included in standard format
   - Now just calls `recommend()` directly
   - Added documentation note about this

## Benefits

### ✅ **Standard Compliance**
- Now follows the same format as ML and Rule-based recommendations
- Compatible with existing formatters and UI components
- Can be processed by `normalize_recommendation_entry()`

### ✅ **Complete Information**
- Includes full chemicals list with amounts, CAS numbers, SMILES
- Includes detailed conditions (temperature, time, atmosphere)
- Retains all protocol metadata for reference

### ✅ **Better User Experience**
- Users get actionable information (what chemicals to use, how much)
- Conditions are extracted and formatted consistently
- Protocol source information available in `protocol_metadata`

## Example Output

```json
{
  "meta": {
    "generated_at": "2025-10-19T15:55:23.498845Z",
    "model": "Protocol-DRFP",
    "model_version": "1.0.0",
    "status": "success",
    "version": "2.0",
    "processing_time_ms": 291.18
  },
  "input": {
    "reaction_smiles": "O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>O=C(C)c1ccc(C2C(CCC2)=O)cc1",
    "options": {"k": 1}
  },
  "detection": {
    "family": "Pd/Enamine_α-Arylation_C(sp3)-C(sp2)",
    "detected_reaction_type": "Pd/Enamine_α-Arylation_C(sp3)-C(sp2)",
    "method": "protocol-similarity",
    "confidence": 1.0
  },
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 1.0,
      "chemicals": [
        {
          "name": "Cyclopentanone",
          "cas": "120-92-3",
          "smiles": "O=C1CCCC1",
          "amount": {"weight_g": 1.68, "mmol": 20.0, "equivalents": 1.0},
          "role": "starting_material"
        },
        {
          "name": "Palladium(II) acetate",
          "abbreviation": "Pd(OAc)2",
          "cas": "3375-31-3",
          "amount": {"weight_g": 0.045, "mmol": 0.2, "equivalents": 0.01},
          "role": "metal_precursor"
        }
        // ... 6 more chemicals
      ],
      "conditions": {
        "catalyst": "Pd(OAc)2",
        "ligand": "P(o-tol)3",
        "base": "NaOAc",
        "solvent": "1,4-Dioxane",
        "additives": ["Pyrrolidine", "tert-octylamine"],
        "temperature_C": 130,
        "time_h": 24.0,
        "atmosphere": "Ar; sealed heavy-walled pressure vessel; oil bath"
      },
      "source": "protocol_database",
      "similarity": 1.0,
      "protocol_metadata": {
        "filename": "alpha_arylation_Pd_enamine_Dong_v100p0099.json",
        "title": "α-Arylation of Cyclopentanones by Palladium/Enamine Cooperative Catalysis",
        "journal": "Organic Syntheses",
        "year": 2023,
        "doi": "10.15227/orgsyn.100.0099",
        "reaction_smarts": ["O=C1CCCC1.Br[a]>>[a]C2C(CCC2)=O"],
        "tags": ["palladium", "enamine", "α-arylation", ...]
      }
    }
  ],
  "extras": {
    "num_candidates": 1,
    "num_total_protocols": 4,
    "num_matches": 1,
    "smarts_filtering_enabled": true
  }
}
```

## Testing

Run the test script to verify:
```bash
python show_protocol_format.py
```

This will display the complete formatted output for inspection.

## Migration Guide

If you have code that accesses the old format:

### Old Code ❌
```python
title = rec['protocol']['title']
family = rec['protocol']['reaction_family']
```

### New Code ✅
```python
title = rec['protocol_metadata']['title']
family = rec['protocol_metadata']['reaction_family']
chemicals = rec['chemicals']  # Now available!
conditions = rec['conditions']  # Now properly formatted!
```

## Summary

The protocol-based recommendation system now:
- ✅ Follows the standard output format
- ✅ Includes complete chemicals list with amounts
- ✅ Includes properly formatted conditions
- ✅ Maintains all protocol metadata in `protocol_metadata`
- ✅ Compatible with existing formatters and normalization functions
- ✅ Provides actionable information for users

This makes protocol recommendations consistent with other recommendation types and ready for production use!
