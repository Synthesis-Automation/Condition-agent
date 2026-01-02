Legacy note: reactant_types.json has been removed; use organic_compounds.v1.3.json for reactant type definitions.

# Reactant-Reaction Type Alignment Summary

## Overview

Successfully aligned the `reactants` field in `reaction_types.json` with the standardized reactant types defined in `reactant_types.json`.

## Changes Made

### Before

- **Issues found**: 76 mismatches
- **Problems**:
  - Free-text descriptions ("ArX", "ArB(OH)2 or ArB(OR)2")
  - Non-standardized names ("R2NH or RNH2", "Carbonyl compound")
  - Generic terms ("Nucleophile", "Substrate", "Various")
  - Reagents mixed with reactants ("Oxidant", "H2O", "Protecting reagent")

### After

- **Issues found**: 0 ✅
- **Reactant type usage**: 18/89 (20.2%)
- **All reactions now reference valid reactant types**

## Mapping Logic

### Electrophiles

- "ArX" → "ArX\*" (category key for all aryl halides)
- "R-X" → "Alkyl-X" (alkyl halides category)
- "Ar-X (activated)" → "ArX\*"

### Nucleophiles

- "R2NH or RNH2" → "Aliphatic-amine"
- "ArNH2" → "Aniline-type"
- "N-nucleophile" → ["Aliphatic-amine", "Aniline-type"]
- "ROH or ArOH" → ["ROH", "ArOH"]

### Boron Reagents

- "ArB(OH)2 or ArB(OR)2" → "ArB\*"
- "B2pin2 or similar" → "RB\*"

### Organometallics

- "RZnX" → "Organometallic"
- "RMgX (Grignard reagent)" → "Organometallic"
- "R-SnR3 (organostannane)" → "Organometallic"

### Carbonyls

- "Carbonyl compound" → ["Aldehyde", "Ketone"]
- "R2C=O" → ["Aldehyde", "Ketone"]

### C-H Donors

- "C-H substrate" → ["Alkyl-C-H", "ArH"]
- "Acidic C-H substrate" → "Alkyl-C-H"

### Unsaturated

- "alkene" → "Alkene"
- "terminal alkyne" → "Alkyne"
- "Unsaturated substrate" → ["Alkene", "Alkyne"]

## Examples of Fixed Reactions

### Suzuki-Miyaura

**Before**: `["ArX", "ArB(OH)2 or ArB(OR)2"]`  
**After**: `["ArX*", "ArB*"]`

### Buchwald-Hartwig

**Before**: `["ArX", "R2NH or RNH2"]`  
**After**: `["ArX*", "Aliphatic-amine", "Aniline-type"]`

### C-N Coupling

**Before**: `["ArX or R-X", "N-nucleophile"]`  
**After**: `["ArX*", "Alkyl-X", "Aliphatic-amine", "Aniline-type"]`

### SNAr

**Before**: `["Ar-X (activated)", "Nucleophile (NH2, OH, SR, etc.)"]`  
**After**: `["ArX*", "Aliphatic-amine", "Aniline-type", "ROH", "ArOH", "RSH"]`

### Alkylation

**Before**: `["Nucleophile", "R-X"]`  
**After**: `["Aliphatic-amine", "Aniline-type", "ROH", "ArOH", "RSH", "Alkyl-X"]`

## Non-Reactant Items Removed

The following were removed as they are reagents/conditions, not reactants:

- "Oxidant"
- "H2O"
- "Zn" (metal)
- "CN source (NaCN, KCN, Zn(CN)2, etc.)"
- "Chlorinating agent"
- "Electrophilic F source"
- "Fluorinating reagent"
- "Protecting reagent"
- "Deprotecting reagent"
- "Activating reagent"
- "Ph3P=CR2" (Wittig reagent)
- "NaNO2/HX"
- "CuX"
- "Base"
- "Acid"

## Reactant Type Categories Used

Currently using 18 reactant type categories:

1. **ArX\*** - Aryl halides (most common)
2. **Alkyl-X** - Alkyl halides
3. **ArB\*** - Aryl boron reagents
4. **RB\*** - Alkyl/alkenyl boron reagents
5. **Organometallic** - Grignard, organozinc, organostannanes
6. **Aliphatic-amine** - Aliphatic amines
7. **Aniline-type** - Aromatic amines
8. **ROH** - Aliphatic alcohols
9. **ArOH** - Phenols
10. **RSH** - Thiols
11. **Alkyl-C-H** - Aliphatic C-H donors
12. **ArH** - Aryl C-H donors
13. **Aldehyde** - Aldehydes
14. **Ketone** - Ketones
15. **Alkene** - Alkenes
16. **Alkyne** - Alkynes
17. **Acyl-source** - Carboxylic acids/esters
18. **Amide-type** - Amides, lactams, carbamates

## Validation

### Tools Created

1. **analyze_reactant_mapping.py** - Analyzes alignment between reaction_types.json and reactant_types.json
2. **fix_reactants_auto.py** - Automatically fixes reactants field using mapping logic

### Results

```
✅ 0 issues found (down from 76)
✅ All 48 reactions validated
✅ 100% compliance with reactant_types.json
```

## Benefits

1. **Consistency**: All reactions now use standardized reactant type references
2. **Maintainability**: Easy to update reactant types in one place (reactant_types.json)
3. **Extensibility**: New reactions can easily reference existing reactant types
4. **Clarity**: Clear separation between reactants (substrates) and reagents (catalysts, oxidants, etc.)
5. **Validation**: Automated tools ensure ongoing compliance

## Files

- **reactant_types.json**: Defines all valid reactant types (21 categories, 71 members)
- **reaction_types.json**: Defines all reactions with standardized reactants field (48 reactions)
- **reaction_types_BACKUP.json**: Backup of original file before update
- **analyze_reactant_mapping.py**: Analysis and validation tool
- **fix_reactants_auto.py**: Automated fixing tool
- **REACTANT_REACTION_ALIGNMENT.md**: This documentation

## Usage

### Validate alignment:

```bash
python data-processor/other_data/analyze_reactant_mapping.py
```

### Check reaction type for specific reactants:

Look up reactant types in `reactant_types.json`, then find reactions using those types in `reaction_types.json`.

### Add new reaction:

Use existing reactant type IDs from `reactant_types.json` in the `reactants` field.

Example:

```json
{
  "id": "New-Reaction",
  "name": "New Reaction",
  "aliases": ["Alternative name"],
  "description": "Description",
  "reactants": ["ArX*", "Aliphatic-amine"],
  "typical_catalysts": ["Catalyst"],
  "typical_conditions": "Conditions"
}
```

## Future Enhancements

1. ~~Add more specific reactant types as needed (e.g., "Heterocyclic-halide", "Enamine")~~ ✅ **COMPLETED** - See `REACTANT_TYPES_ENHANCEMENT.md`
   - Added 7 new categories: Heterocyclic-halide, Enamine, Imines, Allylic-halide, Benzylic-halide, Azide, Nitrile
   - Total categories: 21 → 28
   - Total member types: 71 → 98
2. Create reverse mapping (reactant type → compatible reactions)
3. ~~Add SMARTS patterns to reactant types for automatic detection~~ ✅ **COMPLETED** - See `SMARTS_AUTO_CLASSIFICATION.md`
   - Added SMARTS patterns to all 98 member types
   - Created `classify_reactant.py` for automatic classification
   - Structure-based detection using RDKit
   - Smart prioritization (specific > general functional groups)
   - **Tested on 308 real reactants from sample_reactions.py: 93.5% accuracy** ✅
   - See `SMARTS_TESTING_REPORT.md` for detailed test results
4. Link to specific examples in sample_reactions.py
