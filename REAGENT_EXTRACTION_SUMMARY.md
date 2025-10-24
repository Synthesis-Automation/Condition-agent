# Reagent Registry Extraction - Summary

## Overview

Successfully extracted **552 unique reagents** from the z-Score Peaks dataset to enrich the reagent registry.

## Source Data

- **File:** `z-Score Peaks with FG.csv`
- **Total Reactions:** 66,308
- **Reaction Types:** 42 different types
- **Date Extracted:** October 24, 2025

## Extracted Reagent Categories

| Role                 | Count   | Description                                                    |
| -------------------- | ------- | -------------------------------------------------------------- |
| **Catalyst**         | 141     | Palladium complexes, copper salts, nickel catalysts, etc.      |
| **Ligand**           | 146     | Phosphine ligands, bidentate ligands, NHC ligands, etc.        |
| **Base**             | 63      | Inorganic and organic bases (K₃PO₄, Cs₂CO₃, NaOtBu, DBU, etc.) |
| **Additive**         | 97      | Phase transfer catalysts, acids, protecting groups, etc.       |
| **Solvent**          | 69      | Organic solvents (THF, DMF, toluene, etc.)                     |
| **Coupling Reagent** | 36      | Amide coupling reagents (EDC, PyBrOP, HATU, etc.)              |
| **TOTAL**            | **552** |                                                                |

## Top 10 Most Common Reagents by Category

### Bases (Most Used)

1. **K₃PO₄** - 10,296 occurrences
2. **Cs₂CO₃** - 8,801 occurrences
3. **K₂CO₃** - 4,977 occurrences
4. **NaOtBu** - 4,951 occurrences
5. **Na₂CO₃** - 3,103 occurrences
6. **DBU** - 3,058 occurrences
7. **KOH** - 2,985 occurrences
8. **NEt₃** - 1,888 occurrences
9. **DIPEA** - 1,780 occurrences
10. **NaH** - 1,682 occurrences

### Catalysts (Most Used)

1. **CuI** - 3,786 occurrences
2. **tBuBrettPhos Pd(allyl)OTf** - 2,413 occurrences
3. **dtbpfPdCl₂** - 2,190 occurrences
4. **XantPhos Pd(allyl)Cl** - 2,048 occurrences
5. **P(tBu)₃ Pd(crotyl)Cl** - 1,656 occurrences
6. **SPhos Pd(crotyl)Cl** - 1,532 occurrences
7. **XPhos Pd(crotyl)Cl** - 1,496 occurrences
8. **QPhos Pd(crotyl)Cl** - 1,304 occurrences
9. **AmPhos Pd(crotyl)Cl** - 1,295 occurrences
10. **dppfPdCl₂** - 1,276 occurrences

### Solvents (Most Used)

1. **DMF** - 6,348 occurrences
2. **THF** - 5,848 occurrences
3. **Dioxane** - 5,640 occurrences
4. **PhMe (Toluene)** - 4,760 occurrences
5. **MeCN** - 3,300 occurrences
6. **DMAc** - 2,464 occurrences
7. **iPrOAc** - 2,320 occurrences
8. **MeTHF** - 2,280 occurrences
9. **TBME** - 1,512 occurrences
10. **NMP** - 1,088 occurrences

### Additives (Most Used)

1. **water** - 2,696 occurrences
2. **Potassium 2-ethylhexanoate** - 1,722 occurrences
3. **TBAB** - 1,540 occurrences
4. **Brij 35** - 1,072 occurrences
5. **TBAC** - 1,006 occurrences
6. **HOPO** - 896 occurrences
7. **TMEDA** - 820 occurrences
8. **Octanoic acid** - 728 occurrences
9. **KOPiv** - 544 occurrences

## Reaction Type Coverage

The reagents are used across 42 different reaction types:

**Major Reaction Types:**

- Buchwald-Hartwig (20,286 reactions)
- Suzuki-Miyaura (11,588 reactions)
- Arylation, acidic C-H (4,152 reactions)
- Amide coupling (3,960 reactions)
- C-N Coupling (3,726 reactions)
- C-O Coupling (3,123 reactions)
- And 36 more...

## Output Files

### Main Output: `extracted_reagents_registry.md`

A comprehensive markdown file containing:

- Complete list of all 552 reagents
- Role classification for each reagent
- Usage frequency (number of occurrences)
- Reaction types where each reagent is used
- Cross-reference data for enriching the existing reagent registry

## Data Exclusions

The following were **excluded** as they are reactants/substrates, not reagents:

- Aryl Halide (ArBr, ArCl, ArI)
- N-Nucleophile/Boronate Type (ArB(OR)₂, ArNH₂, R₂NH, etc.)
- Functional Groups A & B (FG A, FG B)

Also excluded:

- "Missing" entries
- "No Ligand" entries
- Empty/null values

## Use Cases for Reagent Registry Enrichment

This extracted data can be used to:

1. **Validate existing registry entries** - Cross-check reagent names and aliases
2. **Add usage statistics** - Show which reagents are most commonly used
3. **Add reaction type associations** - Link reagents to the reactions they're used in
4. **Identify popular combinations** - Find which catalysts pair with which ligands
5. **Coverage analysis** - Identify gaps in the current registry
6. **Real-world validation** - These are reagents used in actual high-throughput experiments

## Next Steps for Integration

1. **Compare with existing registry** - Identify new reagents not in current database
2. **Add CAS numbers** - Look up chemical identifiers for each reagent
3. **Standardize names** - Reconcile naming variations (e.g., "PhMe" vs "Toluene")
4. **Add SMILES structures** - For computational work
5. **Create usage rankings** - Prioritize important reagents based on frequency
6. **Extract reagent combinations** - Analyze which catalyst-ligand-base combinations work best

## Files Generated

- ✅ **extracted_reagents_registry.md** (2,853 lines) - Complete reagent database
- ✅ **analyze_csv.py** - Analysis script
- ✅ **extract_reagents.py** - Extraction script
- ✅ **REAGENT_EXTRACTION_SUMMARY.md** (this file)

---

**Data Quality Notes:**

- Some entries may contain multiple reagents in a single field (comma-separated)
- Reagent names are as-written in the original dataset (may need standardization)
- Usage counts reflect experimental observations from real high-throughput screening data
- This represents validated, experimentally-tested reagent combinations
