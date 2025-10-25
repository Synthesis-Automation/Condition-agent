# Dataset Field Analysis: z-Score Peaks with FG_STANDARDIZED.csv

## Summary
- **Total Rows**: 66,308
- **Total Columns**: 28

## Field Analysis & Recommendations

### ✅ ESSENTIAL FIELDS (Keep)

#### Core Reaction Data
1. **Reaction_Type_Standardized** ✅
   - 42 unique values
   - 0 nulls
   - **Keep**: Primary classification for reaction matching
   - Note: Replaces "Reaction Type" (has same info but standardized format)

2. **AREA_TOTAL_REDUCED** ✅
   - 9,456 unique values (continuous)
   - 1 null only
   - **Keep**: Reaction yield/performance metric

3. **z-Score** ✅
   - 815 unique values (continuous)
   - 1 null only
   - **Keep**: Statistical measure of reaction quality

#### Reactant Classification (Standardized)
4. **Reactant_Type_Electrophile** ✅
   - 18 unique types
   - 12,231 nulls (18.4%)
   - **Keep**: Specific electrophile type (e.g., ArBr, ArCl)

5. **Reactant_Type_Nucleophile** ✅
   - 42 unique types
   - 14,563 nulls (22.0%)
   - **Keep**: Specific nucleophile type (e.g., R2NH, ArB(OR)2)

6. **Reactant_Category_Electrophile** ✅
   - 7 unique categories
   - 12,231 nulls (18.4%)
   - **Keep**: High-level category (e.g., ArX*)

7. **Reactant_Category_Nucleophile** ✅
   - 21 unique categories
   - 14,563 nulls (22.0%)
   - **Keep**: High-level category (e.g., ArB*, Aliphatic-amine)

#### Reaction Conditions
8. **Catalyst** ✅
   - 229 unique catalysts
   - 13,796 nulls (20.8%)
   - **Keep**: Critical for condition recommendation

9. **Ligand** ✅
   - 153 unique ligands
   - 0 nulls
   - **Keep**: Important for catalyst performance

10. **Base** ✅
    - 132 unique bases
    - 9,408 nulls (14.2%)
    - **Keep**: Common reaction component

11. **Solvent** ✅
    - 67 unique solvents
    - 595 nulls (0.9%)
    - **Keep**: Primary solvent

12. **Secondary Solvent** ⚠️
    - 45 unique solvents
    - 46,078 nulls (69.5%)
    - **Consider**: Useful but sparse; could combine with Solvent

13. **Additive** ⚠️
    - 172 unique additives
    - 49,446 nulls (74.6%)
    - **Consider**: Useful but very sparse

14. **Coupling Reagent** ⚠️
    - 36 unique reagents
    - 59,876 nulls (90.3%)
    - **Consider**: Only needed for specific reactions (amide coupling)

---

### ❌ REDUNDANT/LOW-VALUE FIELDS (Remove)

#### Duplicate Information
15. **Reaction Type** ❌
    - **Remove**: Duplicates "Reaction_Type_Standardized" but inconsistent format
    - Same 42 unique values, just different naming convention

16. **FG A** ❌
    - **Remove**: Redundant with Reactant_Type fields
    - Example: "ArBr" duplicates Reactant_Type_Electrophile

17. **FG B** ❌
    - **Remove**: Redundant with Reactant_Type fields
    - Example: "ArB(OR)2" duplicates Reactant_Type_Nucleophile

18. **FG_sorted** ❌
    - **Remove**: Just concatenation of FG A + FG B
    - 139 unique combinations already captured by individual reactant fields

19. **Aryl Halide** ❌
    - **Remove**: Redundant with Reactant_Type_Electrophile
    - Only 18 unique values, same as Reactant_Type_Electrophile

20. **N-Nucleophile/Boronate Type** ❌
    - **Remove**: Redundant with Reactant_Type_Nucleophile
    - 42 unique values, same as Reactant_Type_Nucleophile

#### Lab Metadata (Not Needed for Condition Recommendation)
21. **ELN_ID** ❌
    - **Remove**: Lab notebook reference, not useful for ML/matching
    - 535 unique IDs

22. **PLATENUMBER** ❌
    - **Remove**: Plate number (1-24), lab artifact

23. **Coordinate** ❌
    - **Remove**: Well position (A1-H12), lab artifact

24. **PlateRow** ❌
    - **Remove**: Row letter (A-H), derived from Coordinate

25. **PlateColumn** ❌
    - **Remove**: Column number (1-12), derived from Coordinate

26. **RowLabel** ❌
    - **Remove**: 4,832 unique verbose row descriptions
    - Example: "A: DMAc (10vol)/ water (10vol)/"
    - Information already in individual condition fields

27. **ColumnLabel** ❌
    - **Remove**: 7,152 unique verbose column descriptions
    - Example: "1: XantPhos Pd(allyl)Cl (2mol%)/ K2CO3 (3 eq)/"
    - Information already in individual condition fields

28. **Tertiary Solvent** ❌
    - **Remove**: Only 2 unique values, 65,363 nulls (98.6%)
    - Almost never used

---

## Recommended Simplified Dataset

### Keep These 14 Core Fields:

1. **Reaction_Type_Standardized** - Primary reaction classification
2. **AREA_TOTAL_REDUCED** - Yield/performance
3. **z-Score** - Statistical quality measure
4. **Reactant_Type_Electrophile** - Specific electrophile
5. **Reactant_Type_Nucleophile** - Specific nucleophile
6. **Reactant_Category_Electrophile** - Electrophile category
7. **Reactant_Category_Nucleophile** - Nucleophile category
8. **Catalyst** - Metal catalyst
9. **Ligand** - Ligand
10. **Base** - Base
11. **Solvent** - Primary solvent
12. **Secondary Solvent** - Optional secondary solvent
13. **Additive** - Optional additive
14. **Coupling Reagent** - Optional coupling reagent (for amide coupling)

### Remove These 14 Fields:

- Reaction Type (use standardized version)
- FG A, FG B, FG_sorted (redundant reactant info)
- Aryl Halide, N-Nucleophile/Boronate Type (redundant)
- ELN_ID, PLATENUMBER, Coordinate (lab metadata)
- PlateRow, PlateColumn (lab metadata)
- RowLabel, ColumnLabel (verbose duplicates)
- Tertiary Solvent (98.6% null)

---

## Size Reduction
- **Original**: 28 columns
- **Simplified**: 14 columns
- **Reduction**: 50% fewer columns

## Data Quality Notes
- Nulls are expected for condition fields (not all reactions need all components)
- Reactant nulls (18-22%) indicate reactions without paired electrophile/nucleophile (e.g., Wittig, Deprotection)
- This is normal and expected based on reaction type diversity

## Alignment with reaction_types.json
The simplified dataset aligns well with the updated `reaction_types.json` structure:
- `Reactant_Category_*` maps to category names (ArX*, ArB*, etc.)
- `Reactant_Type_*` maps to specific member types
- `Reaction_Type_Standardized` maps to reaction IDs
