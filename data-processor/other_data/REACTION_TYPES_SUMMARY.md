# Reaction Types JSON - Summary

## Overview
A comprehensive, logically organized JSON database of reaction types covering all reactions in the z-Score Peaks dataset (66,308 reactions) and compatibility with the main reaction dataset.

## Statistics

### Coverage
- **Total reaction types defined**: 42
- **Total aliases**: 66
- **z-Score dataset coverage**: 42/42 (100.0%)
- **Reaction dataset compatibility**: 5/5 (100.0%)

### Organization
Reactions are organized into **10 logical categories**:

1. **Cross-Coupling Reactions** (8 reactions)
   - Suzuki-Miyaura, Negishi, Heck, Sonogashira, Stannylation, etc.
   
2. **C-N Coupling Reactions** (3 reactions)
   - Buchwald-Hartwig, CN-Coupling (general), SNAr
   
3. **C-O Coupling Reactions** (3 reactions)
   - CO-Coupling, Mitsunobu, Glycosidation
   
4. **C-S Coupling Reactions** (1 reaction)
   - CS-Coupling
   
5. **Amide and Carbonyl Chemistry** (4 reactions)
   - Amide coupling, Condensation, Wittig, Stetter
   
6. **C-H Activation and Functionalization** (4 reactions)
   - CH-Activation, Arylation (acidic C-H), Borylation (C-H and Miyaura)
   
7. **Functional Group Transformations** (9 reactions)
   - Alkylation, Oxidation, Reduction, Halogenation, Cyanation, etc.
   
8. **Protection and Deprotection** (2 reactions)
   - Protection, Deprotection
   
9. **Hydrolysis and Water Chemistry** (3 reactions)
   - Hydrolysis, Dehydration, Hydration
   
10. **Miscellaneous Reactions** (5 reactions)
    - Cyclization, Addition, Dimerization, Salt formation, Activation

## Reaction Entry Structure

Each reaction entry contains:
- **id**: Unique identifier (kebab-case)
- **name**: Full descriptive name
- **aliases**: Alternative names (including z-Score dataset naming)
- **description**: Brief explanation of the reaction mechanism/purpose
- **reactants**: Typical reactant types required
- **typical_catalysts**: Common catalysts used
- **typical_conditions**: Standard reaction conditions

## z-Score Dataset Mapping

All 42 unique reaction types from the z-Score Peaks dataset are successfully mapped:

### High-frequency reactions (>3000 occurrences):
- Buchwald-Hartwig: 20,286 reactions
- Suzuki-Miyaura: 11,588 reactions
- Arylation, acidic C-H: 4,152 reactions
- Amide coupling: 3,960 reactions
- CN-Coupling: 3,726 reactions
- CO-Coupling: 3,123 reactions

### Medium-frequency reactions (500-3000):
- Condensation: 2,220 reactions
- CH-Activation: 1,952 reactions
- Negishi, in-situ: 1,752 reactions
- Cyclization: 1,656 reactions
- Borylation, Miyaura: 1,402 reactions
- Suzuki-Miyaura, in situ: 1,296 reactions
- Alkylation: 984 reactions
- Deprotection: 936 reactions
- Negishi: 840 reactions
- Heck: 696 reactions
- CC-Coupling: 680 reactions
- SNAr: 648 reactions
- Hydrolysis: 624 reactions

### Low-frequency reactions (<500):
- All other 24 reaction types

## Reaction Dataset Compatibility

The JSON is fully compatible with the existing reaction dataset families:

| Dataset Family | JSON Reaction ID | Status |
|----------------|------------------|--------|
| Amide formation | Amide-coupling | ✓ |
| Suzuki | Suzuki-Miyaura | ✓ |
| C_N_Coupling | CN-Coupling | ✓ |
| C_O_Coupling | CO-Coupling | ✓ |
| C_S_Coupling | CS-Coupling | ✓ |

## Usage Examples

### Lookup by z-Score dataset name:
```json
"Buchwald-Hartwig" → reaction_types["C-N_bond_formation"]["reactions"][0]
```

### Lookup by reaction dataset family:
```json
"Suzuki" → reaction_types["C-C_bond_formation"]["reactions"][0]
```

### Get all cross-coupling reactions:
```json
reaction_types["C-C_bond_formation"]["reactions"] → 8 reactions
```

## Validation

Use `validate_reaction_types.py` to:
- Check JSON validity
- Verify coverage of z-Score dataset types
- Confirm reaction dataset compatibility
- Generate statistics and mapping reports

## Files

- **reaction_types.json**: Main reaction type database
- **validate_reaction_types.py**: Validation and analysis script
- **extract_reaction_types.py**: Data extraction utility
- **REACTION_TYPES_SUMMARY.md**: This documentation

## Design Principles

1. **Complete Coverage**: Every reaction in both datasets is mapped
2. **Logical Organization**: Grouped by bond formation type and chemistry theme
3. **Flexible Naming**: Aliases support different naming conventions
4. **Rich Metadata**: Includes reactants, catalysts, and conditions for each reaction
5. **Extensible**: Easy to add new reactions or categories
6. **Validated**: 100% verified against actual datasets

## Future Enhancements

Potential additions:
- Reaction conditions details (temperature ranges, typical yields)
- SMARTS patterns for automatic reaction detection
- Literature references for standard procedures
- Compatibility notes (functional group tolerance)
- Reaction mechanism classifications
