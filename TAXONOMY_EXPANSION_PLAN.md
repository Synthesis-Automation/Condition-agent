# Taxonomy Expansion Plan

**Date:** 2025-10-26  
**Current Status:** Analysis of gaps in reaction and reactant taxonomy  
**Goal:** Expand coverage from current ~60% to 90%+ of sample reactions

---

## Current Coverage Analysis

### ✅ Well-Covered Reaction Families (Currently Supported)

1. **C-C Coupling**
   - Suzuki-Miyaura ✓
   - Stille ✓
   - Sonogashira ✓
   - Heck ✓
   - Negishi ✓
   - Kumada ✓

2. **C-N Coupling**
   - Buchwald-Hartwig ✓
   - Ullmann C-N ✓
   - Chan-Lam ✓

3. **C-O Coupling**
   - Ullmann Ether ✓
   - Some SNAr examples ✓

4. **C-S Coupling**
   - Thioether formation ✓

### ❌ Missing Reaction Families (40 reactions marked UNKNOWN)

From `test_sample_reactions.py` output, these reaction types are NOT recognized:

1. **Reduction Reactions (~15 reactions)**
   - Hydrogenation (H2, Pd/C, Pt/C)
   - Metal hydride reductions (NaBH4, LiAlH4, DIBAL)
   - Transfer hydrogenation (Ir, Ru catalysts)
   - Dissolving metal reductions (Zn/HCl, Birch)
   - Selective reductions (NaBH(OAc)3)

2. **Oxidation Reactions (~10 reactions)**
   - Alcohol oxidation (PCC, DMP, Swern, Jones)
   - Epoxidation (mCPBA)
   - Baeyer-Villiger oxidation
   - Sulfide oxidation
   - TEMPO oxidation

3. **Cycloaddition Reactions (~5 reactions)**
   - Diels-Alder [4+2]
   - 1,3-Dipolar cycloaddition (Click chemistry, CuAAC)
   - [2+2] cycloaddition

4. **Substitution Reactions (~8 reactions)**
   - SN2 (improved coverage needed)
   - SN1
   - Nucleophilic aromatic substitution (SNAr)
   - Finkelstein reaction
   - Appel reaction
   - Mitsunobu reaction

5. **Elimination Reactions (~3 reactions)**
   - E2 elimination
   - E1 elimination

6. **Carbonyl Chemistry (~12 reactions)**
   - Aldol reaction
   - Wittig reaction
   - Mannich reaction
   - Michael addition
   - Claisen condensation
   - Knoevenagel condensation
   - Henry (nitroaldol) reaction
   - Reformatsky reaction
   - Horner-Wadsworth-Emmons

7. **Organometallic Reactions (~3 reactions)**
   - Grignard addition

8. **Metathesis (~2 reactions)**
   - Ring-closing metathesis (RCM)
   - Cross-metathesis (CM)
   - ROMP

9. **Heterocycle Synthesis (~10 reactions)**
   - Paal-Knorr (pyrrole synthesis)
   - Hantzsch (dihydropyridine synthesis)
   - Biginelli (dihydropyrimidinone synthesis)
   - Fischer indole synthesis
   - Pictet-Spengler
   - Benzimidazole synthesis
   - Benzoxazole synthesis
   - Benzothiazole synthesis

10. **Amide Formation (~15 reactions)**
    - Carboxylic acid + amine (EDC, DCC, HATU, PyBOP, T3P)
    - Acyl chloride + amine
    - Anhydride + amine
    - Ester aminolysis

11. **Esterification (~5 reactions)**
    - Fischer esterification
    - Acyl chloride + alcohol
    - Mitsunobu esterification

12. **Protecting Group Chemistry (~10 reactions)**
    - Boc protection/deprotection
    - Cbz protection/deprotection
    - FMOC protection/deprotection
    - Silyl ether (TBS, TBDMS) protection/deprotection
    - PMB protection/deprotection
    - Acetonide protection/deprotection
    - Benzyl ether protection/deprotection
    - Trityl protection/deprotection

---

## Missing Reactant Types

### Current Reactant Coverage
- ArX halides ✓
- Amines (ArNH2, RNH2, R2NH) ✓
- Alcohols (ArOH, ROH) ✓
- Thiols (RSH) ✓
- Alkynes ✓

### Need to Add:

1. **Carbonyl Compounds**
   - `RCHO` - Aldehyde ⚠️ (partially detected)
   - `RCOR` - Ketone ⚠️ (partially detected)
   - `RCO2H` - Carboxylic acid ❌
   - `RCOCl` - Acyl chloride ❌
   - `RCOOR` - Ester ❌
   - `(RCO)2O` - Anhydride ❌
   - `RCONH2` - Amide ❌

2. **Organometallic Reagents**
   - `RMgX` - Grignard reagent ❌
   - `RLi` - Organolithium ❌
   - `R2Zn` - Organozinc ❌
   - `R3B` or `R-B(OH)2` - Organoboron (partially covered)
   - `RSnR3` - Organostannane ❌

3. **Hydrides**
   - `NaBH4` - Sodium borohydride ❌
   - `LiAlH4` - Lithium aluminum hydride ❌
   - `DIBAL` - Diisobutylaluminum hydride ❌
   - `Red-Al` - Sodium bis(2-methoxyethoxy)aluminum hydride ❌
   - `H2` - Hydrogen gas ❌

4. **Oxidants**
   - `mCPBA` - meta-Chloroperoxybenzoic acid ❌
   - `PCC` - Pyridinium chlorochromate ❌
   - `DMP` - Dess-Martin periodinane ❌
   - `IBX` - 2-Iodoxybenzoic acid ❌
   - `TEMPO` - (2,2,6,6-Tetramethylpiperidin-1-yl)oxyl ❌
   - `CrO3` - Chromium trioxide ❌
   - `KMnO4` - Potassium permanganate ❌
   - `DDQ` - 2,3-Dichloro-5,6-dicyano-1,4-benzoquinone ❌

5. **Alkenes & Dienes**
   - `terminal-alkene` ✓
   - `internal-alkene` ❌
   - `conjugated-diene` ❌
   - `dienophile` ❌

6. **Nitriles & Isocyanates**
   - `RCN` - Nitrile ❌
   - `RNCO` - Isocyanate ❌

7. **Epoxides & Aziridines**
   - `epoxide` ❌
   - `aziridine` ❌

8. **Protecting Groups (as detectible entities)**
   - `Boc` - tert-Butoxycarbonyl ❌
   - `Cbz` - Benzyloxycarbonyl ❌
   - `FMOC` - Fluorenylmethoxycarbonyl ❌
   - `TBS`/`TBDMS` - tert-Butyldimethylsilyl ❌
   - `PMB` - para-Methoxybenzyl ❌
   - `Bn` - Benzyl ❌
   - `Tr` - Trityl ❌

9. **Ylides & Enolates**
   - `Wittig-ylide` - Phosphonium ylide ❌
   - `HWE-phosphonate` - Phosphonate ester ❌
   - `enolate` - Enolate ion ❌

10. **Azides & Diazo**
    - `RN3` - Azide ❌
    - `RN2+` - Diazonium ❌
    - `R2CN2` - Diazo compound ❌

---

## Proposed Taxonomy Expansion

### Phase 1: High-Priority Reaction Families (Weeks 1-2)

**Add to `reaction_types.json`:**

1. **Reduction Family**
   ```json
   {
     "id": "hydrogenation",
     "name": "Hydrogenation",
     "description": "Catalytic addition of H2 to unsaturated bonds",
     "category": "reduction",
     "required_reactants": ["unsaturated"],
     "conditions": {"reagents": ["H2", "Pd/C", "Pt/C", "Raney Ni"]}
   }
   ```

2. **Metal Hydride Reduction Family**
   ```json
   {
     "id": "metal_hydride_reduction",
     "name": "Metal Hydride Reduction",
     "category": "reduction",
     "required_reactants": ["carbonyl", "hydride"],
     "conditions": {"reagents": ["NaBH4", "LiAlH4", "DIBAL"]}
   }
   ```

3. **Oxidation Family**
   ```json
   {
     "id": "alcohol_oxidation",
     "name": "Alcohol Oxidation",
     "category": "oxidation",
     "required_reactants": ["alcohol"],
     "conditions": {"reagents": ["PCC", "DMP", "Swern", "Jones"]}
   }
   ```

4. **Amide Formation Family**
   ```json
   {
     "id": "amide_coupling",
     "name": "Amide Coupling",
     "category": "condensation",
     "required_reactants": ["acid_or_derivative", "amine"],
     "conditions": {"coupling_agents": ["EDC", "DCC", "HATU", "PyBOP", "T3P"]}
   }
   ```

5. **Cycloaddition Family**
   ```json
   {
     "id": "diels_alder",
     "name": "Diels-Alder Reaction",
     "category": "cycloaddition",
     "required_reactants": ["diene", "dienophile"]
   }
   ```

### Phase 2: Medium-Priority (Weeks 3-4)

6. **Carbonyl Condensation Family**
   - Aldol
   - Wittig
   - Mannich
   - Michael
   - Claisen

7. **Substitution Family Enhancement**
   - SN1/SN2 improved detection
   - SNAr
   - Mitsunobu
   - Appel

8. **Heterocycle Synthesis**
   - Paal-Knorr
   - Hantzsch
   - Biginelli
   - Fischer Indole
   - Pictet-Spengler

### Phase 3: Lower-Priority (Weeks 5-6)

9. **Protecting Group Operations**
   - Protection reactions
   - Deprotection reactions

10. **Metathesis**
    - RCM
    - CM
    - ROMP

11. **Rearrangements**
    - Various named rearrangements

---

## Proposed Reactant Type Expansion

### Phase 1: Essential Reactant Types

**Add to `reactant_types.json`:**

```json
{
  "id": "carboxylic_acid_group",
  "description": "Carboxylic acids and derivatives",
  "members": [
    {
      "id": "RCO2H",
      "smarts": "[CX3](=O)[OX2H]",
      "name": "carboxylic acid"
    },
    {
      "id": "RCOCl",
      "smarts": "[CX3](=O)[Cl]",
      "name": "acyl chloride"
    },
    {
      "id": "RCOOR",
      "smarts": "[CX3](=O)[OX2][#6]",
      "name": "ester"
    },
    {
      "id": "anhydride",
      "smarts": "[CX3](=O)[OX2][CX3](=O)",
      "name": "carboxylic anhydride"
    }
  ]
}
```

```json
{
  "id": "carbonyl_group",
  "description": "Aldehydes and ketones",
  "members": [
    {
      "id": "RCHO",
      "smarts": "[CX3H1](=O)[#6]",
      "name": "aldehyde"
    },
    {
      "id": "RCOR",
      "smarts": "[#6][CX3](=O)[#6]",
      "name": "ketone"
    }
  ]
}
```

```json
{
  "id": "hydride_reagents",
  "description": "Metal hydride reducing agents",
  "members": [
    {
      "id": "NaBH4",
      "smarts": "[Na+].[BH4-]",
      "name": "sodium borohydride"
    },
    {
      "id": "LiAlH4",
      "smarts": "[Li+].[AlH4-]",
      "name": "lithium aluminum hydride"
    }
  ]
}
```

```json
{
  "id": "grignard_organozinc",
  "description": "Organometallic nucleophiles",
  "members": [
    {
      "id": "RMgX",
      "smarts": "[#6][Mg][Br,Cl,I]",
      "name": "Grignard reagent"
    },
    {
      "id": "RZnX",
      "smarts": "[#6][Zn][Br,Cl]",
      "name": "organozinc halide"
    }
  ]
}
```

---

## Implementation Strategy

### Step 1: Extend Reactant Types First
**Rationale:** Reaction detection depends on reactant identification

**Priority Order:**
1. Carbonyl compounds (RCO2H, RCOCl, RCOOR, RCHO, RCOR)
2. Hydride reagents (NaBH4, LiAlH4, DIBAL)
3. Organometallics (RMgX, RZnX, RLi)
4. Oxidants (mCPBA, PCC, DMP)
5. Alkenes/dienes (diene, dienophile)

### Step 2: Add Reaction Families
**Rationale:** Use new reactant types to classify reactions

**Priority Order:**
1. Amide coupling (most common in pharma)
2. Reductions (metal hydrides, hydrogenation)
3. Oxidations (alcohol oxidation, epoxidation)
4. Cycloadditions (Diels-Alder, Click)
5. Carbonyl chemistry (Aldol, Wittig, etc.)

### Step 3: Create Detection Rules
**Location:** `chemtools/reaction_type_detector.py`

Add detection logic for each new family based on:
- Required reactant types
- Product structure analysis
- Reagent/catalyst presence
- Functional group transformations

### Step 4: Validate with Sample Reactions
**Goal:** Reduce "UNKNOWN" reactions from 40 to <10

Run `test_sample_reactions.py` after each addition to verify:
- Correct family assignment
- No false positives
- Proper reactant type detection

---

## Expected Impact

### Before Expansion:
- **Reactions classified:** 62/102 (60.8%)
- **Reactions UNKNOWN:** 40/102 (39.2%)

### After Phase 1:
- **Reactions classified:** ~85/102 (83.3%)
- **Reactions UNKNOWN:** ~17/102 (16.7%)

### After Full Expansion:
- **Reactions classified:** 95+/102 (93%+)
- **Reactions UNKNOWN:** <7/102 (<7%)

---

## File Structure Changes

### New Files to Create:

1. **`chemtools/taxonomy/data/carbonyl_reactants.json`**
   - Aldehydes, ketones, acids, esters, acyl chlorides, anhydrides

2. **`chemtools/taxonomy/data/organometallic_reactants.json`**
   - Grignard, organozinc, organolithium, organoboron

3. **`chemtools/taxonomy/data/reagent_library.json`**
   - Hydrides (NaBH4, LiAlH4)
   - Oxidants (PCC, DMP, mCPBA)
   - Coupling agents (EDC, HATU, DCC)

4. **`chemtools/taxonomy/data/reduction_reactions.json`**
   - Hydrogenation, metal hydride reduction, transfer hydrogenation

5. **`chemtools/taxonomy/data/oxidation_reactions.json`**
   - Alcohol oxidation, epoxidation, Baeyer-Villiger

6. **`chemtools/taxonomy/data/cycloaddition_reactions.json`**
   - Diels-Alder, Click chemistry, [2+2] cycloaddition

7. **`chemtools/taxonomy/data/condensation_reactions.json`**
   - Amide coupling, esterification, aldol, Claisen

### Files to Extend:

1. **`chemtools/taxonomy/data/reactant_types.json`**
   - Add carbonyl, organometallic, reagent categories

2. **`chemtools/taxonomy/data/reaction_types.json`**
   - Add new reaction families

3. **`chemtools/taxonomy/data/manifest.json`**
   - Register new data files

4. **`chemtools/reaction_type_detector.py`**
   - Add detection logic for new families

---

## Testing Plan

### Unit Tests to Create:

1. **`tests/test_carbonyl_reactants.py`**
   - Test RCHO, RCOR, RCO2H, RCOCl detection

2. **`tests/test_reduction_reactions.py`**
   - Test hydrogenation, NaBH4, LiAlH4 reactions

3. **`tests/test_oxidation_reactions.py`**
   - Test alcohol oxidation, epoxidation

4. **`tests/test_cycloaddition_reactions.py`**
   - Test Diels-Alder, Click chemistry

5. **`tests/test_amide_coupling.py`**
   - Test acid+amine, EDC coupling, etc.

### Integration Tests:

Update `test_sample_reactions.py` to verify:
- All 102 reactions get correct family assignment
- Reactant types detected for all substrates
- No regression in existing classifications

---

## Success Criteria

1. ✅ **Coverage:** 90%+ of sample reactions classified
2. ✅ **Accuracy:** <5% false positive rate
3. ✅ **Completeness:** All major organic reaction types covered
4. ✅ **Usability:** Clear human-readable names and descriptions
5. ✅ **Performance:** Classification <1 second per reaction

---

## Timeline

- **Week 1:** Phase 1 reactant types (carbonyl, hydrides)
- **Week 2:** Phase 1 reaction families (reduction, oxidation, amide)
- **Week 3:** Phase 2 reactant types (organometallics, reagents)
- **Week 4:** Phase 2 reaction families (condensations, cycloadditions)
- **Week 5:** Phase 3 (protecting groups, metathesis)
- **Week 6:** Testing, validation, documentation

---

## Next Steps

**Immediate Actions:**

1. Create JSON files for Phase 1 reactant types
2. Add SMARTS patterns for carbonyl compounds
3. Implement detection logic for amide coupling
4. Create unit tests for new types
5. Run validation on sample reactions

**Would you like me to:**
- ✅ Start implementing Phase 1 reactant types?
- ✅ Create the JSON files for carbonyl compounds?
- ✅ Add detection logic for amide coupling?
- ✅ All of the above?

---

**Document Version:** 1.0  
**Author:** GitHub Copilot  
**Date:** 2025-10-26
