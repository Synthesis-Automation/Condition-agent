Legacy note: reactant_types.json has been removed; use organic_compounds.v1.3.json for reactant type definitions.

# Reactant Types Enhancement Summary

## Overview

Enhanced `reactant_types.json` with **7 new specific reactant type categories** to improve condition recommendation precision and coverage of common substrate classes.

## Previous Stats

- **Categories**: 21
- **Total member types**: 71

## Updated Stats

- **Categories**: 28 (+7 new)
- **Total member types**: 98 (+27 new)

## New Reactant Type Categories Added

### 1. Heterocyclic-halide (8 members)

**Group**: Electrophiles

**Description**: Heteroaryl halides (pyridines, pyrimidines, thiazoles, indoles, etc.) - often more challenging than simple ArX due to electronic effects and coordination.

**Members**:

- `HetArBr`: heteroaryl bromide
- `HetArCl`: heteroaryl chloride
- `HetArI`: heteroaryl iodide
- `HetArOTf`: heteroaryl triflate
- `PyridineBr`: pyridine bromide
- `PyrimidineCl`: pyrimidine chloride
- `ThiazoleBr`: thiazole bromide
- `IndoleBr`: indole bromide

**Rationale**: Heteroaryl halides are extremely common in medicinal chemistry but often require different conditions than simple aryl halides due to:

- N-coordination to metal catalysts
- Different electronic properties (electron-deficient vs. electron-rich)
- Regioselectivity challenges

### 2. Enamine (3 members)

**Group**: Nucleophiles

**Description**: Enamines - versatile carbon nucleophiles formed from aldehydes/ketones and secondary amines; used in organocatalysis and C-C bond formation.

**Members**:

- `enamine`: enamine (R2C=C-NR2)
- `morpholine-enamine`: morpholine-derived enamine
- `pyrrolidine-enamine`: pyrrolidine-derived enamine

**Rationale**: Enamines are key intermediates in:

- Organocatalytic reactions (Stork enamine chemistry)
- Mannich reactions
- Michael additions
  Different enamine sources (morpholine vs. pyrrolidine) give different reactivity and selectivity.

### 3. Imines (3 members)

**Group**: Nucleophiles

**Description**: Imines and Schiff bases - electrophilic and nucleophilic depending on activation; common in reductive amination and cycloadditions.

**Members**:

- `imine`: imine (R2C=NR)
- `Ar-imine`: aromatic imine
- `N-protected-imine`: N-protected imine (e.g., N-Ts, N-Boc)

**Rationale**: Imines are crucial for:

- Reductive amination
- Cycloadditions
- Condensation reactions
  Protection groups (Ts, Boc) dramatically affect reactivity.

### 4. Allylic-halide (5 members)

**Group**: Electrophiles

**Description**: Allylic halides and leaving groups - undergo SN2', allylic substitution, and cross-coupling with retention or inversion.

**Members**:

- `Allyl-Br`: allylic bromide
- `Allyl-Cl`: allylic chloride
- `Allyl-I`: allylic iodide
- `Allyl-OAc`: allylic acetate
- `Allyl-OCO2R`: allylic carbonate

**Rationale**: Allylic systems have unique reactivity:

- Allylic substitution (SN2' mechanism)
- π-allyl complexes in metal catalysis
- Regioselectivity (α vs. γ attack)
- Different conditions needed vs. simple alkyl halides

### 5. Benzylic-halide (3 members)

**Group**: Electrophiles

**Description**: Benzylic halides - more reactive than simple alkyl halides due to resonance stabilization of carbocation/radical intermediates.

**Members**:

- `Bn-Br`: benzylic bromide
- `Bn-Cl`: benzylic chloride
- `Bn-I`: benzylic iodide

**Rationale**: Benzylic halides are:

- More reactive than alkyl halides
- Prone to SN1 and radical pathways
- Common protecting group intermediates
- Require milder conditions than Alkyl-X

### 6. Azide (3 members)

**Group**: Nucleophiles

**Description**: Azides - versatile nucleophiles for click chemistry, reduction to amines, and cycloadditions.

**Members**:

- `R-N3`: alkyl azide
- `Ar-N3`: aryl azide
- `NaN3`: sodium azide

**Rationale**: Azides are essential for:

- Click chemistry (CuAAC, SPAAC)
- Aza-transfer reactions
- Reduction to primary amines
- Staudinger ligation
  Different from typical nitrogen nucleophiles (amines, anilines).

### 7. Nitrile (2 members)

**Group**: Neutral / Directing

**Description**: Nitriles - directing groups for C-H activation and precursors to amides, amines, carboxylic acids.

**Members**:

- `R-CN`: alkyl nitrile
- `Ar-CN`: aryl nitrile

**Rationale**: Nitriles are important as:

- Directing groups in C-H activation
- Electrophiles in nucleophilic addition
- Precursors to carboxylic acids, amides, amines
- Different reactivity vs. carbonyls

## Impact on Existing Categories

The new categories complement the existing taxonomy:

### Electrophiles (now 8 categories, 37 total members)

- `ArX*` (8) - Simple aryl halides
- **`Heterocyclic-halide` (8)** - **NEW**: Heteroaryl halides
- `VinylX*` (4) - Vinyl halides
- `Alkyl-X` (4) - Simple alkyl halides
- **`Allylic-halide` (5)** - **NEW**: Allylic systems
- **`Benzylic-halide` (3)** - **NEW**: Benzylic systems
- `RSO2Cl` (2) - Sulfonyl chlorides
- `Acyl-source-electrophile` (3) - Acyl chlorides/anhydrides

### Nucleophiles (now 11 categories, 38 total members)

- `ArB*` (3) - Aryl boron
- `RB*` (5) - Alkyl/vinyl boron
- `Organometallic` (2) - Grignards, organozincs
- `Aliphatic-amine` (5) - Aliphatic amines
- `Aniline-type` (4) - Aromatic amines
- **`Enamine` (3)** - **NEW**: Enamine nucleophiles
- **`Imines` (3)** - **NEW**: Imine/Schiff bases
- `ROH` (4) - Aliphatic alcohols
- `ArOH` (1) - Phenols
- `RSH` (2) - Thiols
- **`Azide` (3)** - **NEW**: Azide nucleophiles

### Neutral / Directing (now 7 categories, 21 total members)

- `Amide-type` (6) - Amides, ureas, carbamates
- `Acyl-source` (3) - Carboxylic acids/esters
- `Alkene` (4) - Alkenes
- `Alkyne` (3) - Alkynes
- `Aldehyde` (2) - Aldehydes
- `Ketone` (2) - Ketones
- **`Nitrile` (2)** - **NEW**: Nitrile groups

## Benefits for Condition Recommendation

### 1. **Improved Specificity**

- Heteroaryl halides vs. simple aryl halides → different catalyst/ligand preferences
- Allylic vs. alkyl halides → different mechanisms and conditions
- Enamines vs. simple amines → different nucleophilicity

### 2. **Better Coverage**

- Now covers ~98 distinct substrate types across common reaction classes
- Captures fine-grained reactivity differences that affect condition selection

### 3. **Enhanced Precision**

- Recommendations can now distinguish:
  - Pyridine bromide vs. phenyl bromide
  - Benzylic bromide vs. primary alkyl bromide
  - Morpholine enamine vs. primary amine
  - Alkyl azide vs. aniline

### 4. **Mechanistic Insight**

- Each new category represents substrates with distinct:
  - Electronic properties
  - Steric profiles
  - Preferred reaction pathways
  - Optimal conditions

## Validation & Testing

### Next Steps

1. **Re-run standardization** on z-Score dataset to incorporate new types:

   ```bash
   python data-processor/other_data/process_zscore_csv.py
   ```

2. **Update reaction_types.json** if any reactions specifically use these new types:

   ```bash
   python data-processor/other_data/analyze_reactant_mapping.py
   ```

3. **Test recommendations** with new substrate types:

   ```python
   # Example: Test heteroaryl halide recommendation
   recommender.recommend(
       reaction_type="Buchwald-Hartwig",
       electrophile_type="HetArBr",  # NEW
       nucleophile_type="RNH2"
   )
   ```

4. **Analyze coverage** in dataset:
   ```python
   python data-processor/other_data/analyze_standardized_types.py
   ```

## Files Modified

- **reactant_types.json**: Added 7 new categories with 27 new member types
- **verify_new_types.py**: Created to validate the additions

## Backward Compatibility

✅ **Fully backward compatible**

- All existing 21 categories retained
- All existing 71 member types retained
- Only additions made, no modifications to existing types
- Existing code using old types will continue to work

## Example Use Cases

### Heterocyclic-halide

```json
{
  "reaction": "Buchwald-Hartwig",
  "electrophile": "2-bromopyridine",
  "electrophile_type": "PyridineBr",
  "nucleophile": "morpholine",
  "conditions": "...different from ArBr due to N-coordination"
}
```

### Enamine

```json
{
  "reaction": "Organocatalytic-Michael",
  "nucleophile": "pyrrolidine-cyclohexanone enamine",
  "nucleophile_type": "pyrrolidine-enamine",
  "electrophile": "methyl vinyl ketone"
}
```

### Allylic-halide

```json
{
  "reaction": "Allylic-substitution",
  "electrophile": "cinnamyl chloride",
  "electrophile_type": "Allyl-Cl",
  "nucleophile": "malonate",
  "conditions": "...SN2' vs SN2 selectivity"
}
```

## Summary

Successfully enhanced the reactant type taxonomy from **21 → 28 categories** and **71 → 98 member types**, adding critical substrate classes that are common in synthesis but have distinct reactivity profiles requiring different reaction conditions. This will significantly improve the precision and coverage of condition recommendations.
