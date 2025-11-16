# README Update Summary

## Changes Made

Updated the HTE README (`chemtools/HTE/README.md`) to include comprehensive CLI usage examples with reaction type and catalyst filters.

## Key Additions

### 1. Expanded CLI Examples Section

Reorganized and expanded the CLI section with clear subsections:

- **Basic Queries**: Simple reactant pair queries
- **Filtering by Reaction Type**: Examples for Suzuki, C_N_Coupling, amide_formation
- **Filtering by Catalyst Type**: Examples for Cu, Pd, Ni (symbols and full names)
- **Combined Filters**: Examples combining reaction type and catalyst filters
- **Output Formats**: Compact and JSON examples with filters
- **Batch Processing**: Batch processing with filters
- **Full Example**: Comprehensive example using all options

### 2. Updated Python API Examples

Added examples showing:
- Basic queries
- Reaction type filtering
- Catalyst metal filtering
- Combined filters (reaction type + catalyst)
- Accessing z-score results in output

### 3. Updated Features List

- Added "Z-Score ranking" as primary feature
- Clarified catalyst filtering supports multiple formats (Cu, Pd, Ni, copper, palladium, etc.)

### 4. Updated Ranking & Scoring Section

- Changed from "Confidence Scoring" to "Ranking & Scoring"
- Documented z-score as primary ranking metric
- Added z-score interpretation guide
- Updated confidence score formula (60% z-score, 25% sample size, 15% yield)

### 5. Updated API Reference

- Added `catalyst_filter` parameter to `HTERecommender.recommend()`
- Added `avg_z_score` field to `ConditionRecommendation` dataclass
- Documented that avg_z_score is the PRIMARY ranking metric

## Example Commands Added

```bash
# Copper-catalyzed C-N coupling
python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "c1ccc(N)cc1" \
    --reaction C_N_Coupling --catalyst Cu -k 10 --min-exp 1

# Palladium-catalyzed Suzuki coupling
python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "c1ccc(B(O)O)cc1" \
    --reaction Suzuki --catalyst Pd -k 10

# With JSON output
python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "c1ccc(N)cc1" \
    --catalyst Cu --json -o cu_cn_results.json

# Full example with all options
python -m chemtools.HTE.cli \
    -a "c1ccc(Br)cc1" \
    -b "c1ccc(N)cc1" \
    --reaction C_N_Coupling \
    --catalyst Cu \
    -k 20 \
    --min-exp 1 \
    --json \
    -o cu_cn_coupling_conditions.json
```

## Documentation Improvements

1. **Clearer organization**: Separate sections for different filter types
2. **More examples**: 15+ new CLI examples covering various use cases
3. **Better explanations**: Comments explaining what each filter does
4. **Z-score prominence**: Z-score is now clearly documented as the primary ranking metric
5. **Catalyst flexibility**: Documented that catalyst names can be symbols (Cu, Pd) or full names (copper, palladium)

## Testing

All examples tested and verified working with the updated z-score ranking system:
- ✅ Catalyst filters work (Cu, Pd, Ni, copper, palladium, nickel)
- ✅ Reaction type filters work (C_N_Coupling, Suzuki, etc.)
- ✅ Combined filters work correctly
- ✅ Output shows z-score prominently
- ✅ Compact format includes z-score

## Files Modified

- `chemtools/HTE/README.md` - Comprehensive CLI documentation update
