# z-Score CSV Standardization Summary

## Overview
Successfully processed the z-Score Peaks with FG.csv dataset (66,308 reactions) to use standardized reactant and reaction type nomenclature from our `reactant_types.json` and `reaction_types.json` taxonomy systems.

## Files Created

### 1. `process_zscore_csv.py`
Main processing script that:
- Loads reactant and reaction type definitions from JSON files
- Creates mapping dictionaries from CSV values to standardized IDs
- Processes all 66,308 reactions
- Adds 5 new standardized columns
- Saves output to `z-Score Peaks with FG_STANDARDIZED.csv`

### 2. `z-Score Peaks with FG_STANDARDIZED.csv`
Output file with **28 columns** (original 23 + 5 new standardized columns):

**New Columns Added:**
1. `Reactant_Type_Electrophile` - Standardized reactant type ID for electrophile (e.g., "ArBr", "ArCl")
2. `Reactant_Category_Electrophile` - Category key for electrophile (e.g., "ArX*", "Alkyl-X")
3. `Reactant_Type_Nucleophile` - Standardized reactant type ID for nucleophile (e.g., "ArB(OH)2", "RNH2")
4. `Reactant_Category_Nucleophile` - Category key for nucleophile (e.g., "ArB*", "Aliphatic-amine")
5. `Reaction_Type_Standardized` - Standardized reaction type ID (e.g., "Suzuki-Miyaura", "Buchwald-Hartwig")

### 3. `analyze_standardized_types.py`
Analysis script showing usage statistics of standardized types.

## Processing Statistics

### Coverage
- **Total reactions processed:** 66,308
- **Electrophiles mapped:** 54,077 / 54,077 (100% of non-null values)
- **Nucleophiles mapped:** 51,745 / 51,745 (100% of non-null values)
- **Reactions mapped:** 66,308 / 66,308 (100%)

### Reactant Types Used

**Electrophile Types (12 unique):**
- ArBr (32,826) - Most common
- ArCl (15,729)
- ArI (3,680)
- Alkyl-Br (1,368)
- ArOH (1,344)
- Alkyl-Cl (510)
- Alkyl-OSO2R (480)
- Alkyl-I (432)
- ArOSO2R (432)
- ArF (384)
- alkene-Br (192)
- alkene-I (96)

**Nucleophile Types (31 unique):**
- RNH2 (9,466) - Most common
- R2NH (6,682)
- RNH2-a-branch (6,051)
- ArNH2 (5,860)
- ArB(OR)2 (5,308)
- Alkyl-H (4,776)
- Alkyl-H-acidic (4,464)
- RCO2H (3,984)
- R2NH-a-branch (3,648)
- arom-NH (3,424)
- ArB(OH)2 (3,188)
- alkene (2,448)
- Lactam (1,536)
- alkeneB(OR)2 (1,440)
- ROH-a-branch (1,368)
- ArOH (1,344)
- ArNHR (1,032)
- Carbamate (912)
- alkyne (849)
- RCONH2 (840)
- Urea (606)
- Alkyl-B(OH)2 (528)
- Alkyl-M (456)
- Ar2NH (384)
- ArBF3K (336)
- Alkyl-B(OR)2 (210)
- ROH-primary (192)
- Ar-H (168)
- Alkyl-BF3K (96)
- enol-ether (96)
- RSH (24)

### Reaction Types Distribution (Top 15)

1. **Buchwald-Hartwig** - 20,286 (30.6%)
2. **Suzuki-Miyaura** - 11,588 (17.5%)
3. **Arylation-acidic-C-H** - 4,152 (6.3%)
4. **Amide-coupling** - 3,960 (6.0%)
5. **CN-Coupling** - 3,726 (5.6%)
6. **CO-Coupling** - 3,123 (4.7%)
7. **Condensation** - 2,220 (3.3%)
8. **CH-Activation** - 1,952 (2.9%)
9. **Negishi-in-situ** - 1,752 (2.6%)
10. **Cyclization** - 1,656 (2.5%)
11. **Borylation-Miyaura** - 1,402 (2.1%)
12. **Suzuki-Miyaura-in-situ** - 1,296 (2.0%)
13. **Alkylation** - 984 (1.5%)
14. **Deprotection** - 936 (1.4%)
15. **Negishi** - 840 (1.3%)

## Mapping Logic

### Reactant Mapping Strategy
1. **Direct ID matching** - Exact matches to member IDs in `reactant_types.json`
2. **Category grouping** - Multiple specific types grouped under category keys (e.g., ArBr, ArCl, ArI → ArX*)
3. **Separation of concerns** - Electrophiles vs. Nucleophiles tracked separately
4. **Mixed value handling** - Comma-separated values split and mapped individually

### Reaction Mapping Strategy
1. **Name and alias matching** - Uses both primary name and all aliases from `reaction_types.json`
2. **Variant handling** - Distinguishes between base reactions and variants (e.g., "Suzuki-Miyaura" vs. "Suzuki-Miyaura-in-situ")
3. **100% coverage** - All 42 unique reaction types from CSV successfully mapped

## Benefits of Standardization

### Data Quality
✅ **Consistent nomenclature** - All reactant and reaction types use standardized IDs
✅ **No free-text descriptions** - Structured, machine-readable format
✅ **Hierarchical organization** - Both specific types and broader categories available
✅ **100% coverage** - No unmapped values for non-null data

### Analysis Capabilities
✅ **Category-level analysis** - Can analyze by broad categories (ArX*, ArB*, etc.)
✅ **Type-level analysis** - Can drill down to specific types (ArBr, ArCl, etc.)
✅ **Cross-referencing** - Standardized IDs enable linking to other databases
✅ **Statistical aggregation** - Easy to count, group, and analyze reaction patterns

### Integration
✅ **Compatible with `reaction_types.json`** - Can now link reactions to their definitions
✅ **Compatible with `reactant_types.json`** - Can look up detailed reactant properties
✅ **API-ready** - Structured format suitable for programmatic access
✅ **Machine learning ready** - Standardized features for ML models

## Example Transformations

### Before (Original CSV)
```
Aryl Halide: "ArBr"
N-Nucleophile/Boronate Type: "ArB(OR)2"
Reaction Type: "Suzuki-Miyaura"
```

### After (Standardized CSV)
```
Aryl Halide: "ArBr"
Reactant_Type_Electrophile: "ArBr"
Reactant_Category_Electrophile: "ArX*"

N-Nucleophile/Boronate Type: "ArB(OR)2"
Reactant_Type_Nucleophile: "ArB(OR)2"
Reactant_Category_Nucleophile: "ArB*"

Reaction Type: "Suzuki-Miyaura"
Reaction_Type_Standardized: "Suzuki-Miyaura"
```

## Usage Examples

### Load and Filter by Reaction Type
```python
import pandas as pd

df = pd.read_csv("z-Score Peaks with FG_STANDARDIZED.csv")

# Get all Buchwald-Hartwig reactions
bh_reactions = df[df['Reaction_Type_Standardized'] == 'Buchwald-Hartwig']
print(f"Found {len(bh_reactions)} Buchwald-Hartwig reactions")
```

### Analyze by Reactant Category
```python
# Get all reactions using aryl halides
aryl_halide_reactions = df[df['Reactant_Category_Electrophile'] == 'ArX*']

# Count by specific type
aryl_bromide = df[df['Reactant_Type_Electrophile'] == 'ArBr']
aryl_chloride = df[df['Reactant_Type_Electrophile'] == 'ArCl']
```

### Cross-Reference with Reaction Types Definition
```python
import json

# Load reaction type definition
with open('reaction_types.json', 'r') as f:
    reaction_types = json.load(f)

# For each reaction in dataset, look up its definition
for idx, row in df.iterrows():
    reaction_id = row['Reaction_Type_Standardized']
    # Find definition in reaction_types.json...
```

## Validation

### Perfect Mapping Achievement
- ✅ Every non-null electrophile value mapped (54,077/54,077)
- ✅ Every non-null nucleophile value mapped (51,745/51,745)
- ✅ Every reaction type mapped (66,308/66,308)
- ✅ No unmapped values or errors

### Data Integrity
- ✅ All standardized IDs exist in `reactant_types.json` or `reaction_types.json`
- ✅ Categories correctly reference parent category keys
- ✅ Multiple values in same field properly split and mapped
- ✅ Original data preserved alongside standardized versions

## Files Location
```
data-processor/other_data/
  ├── reactant_types.json                      # Source taxonomy
  ├── reaction_types.json                      # Source taxonomy
  ├── z-Score Peaks with FG.csv               # Original dataset
  ├── z-Score Peaks with FG_STANDARDIZED.csv  # Processed output ⭐
  ├── process_zscore_csv.py                    # Processing script
  └── analyze_standardized_types.py            # Analysis script
```

## Next Steps

The standardized CSV can now be used for:
1. **Training ML models** with consistent feature representation
2. **Statistical analysis** of reaction success patterns
3. **Recommendation systems** linking similar reactions
4. **API integration** with structured, validated data
5. **Cross-dataset comparisons** using standardized nomenclature
