# Calculable Features JSON - Expansion Plan

**Date**: November 2, 2025  
**Current Version**: 2025-11-02.v1.1  
**Proposed Version**: 2025-11-02.v2.0

## Executive Summary

Based on validation testing with 119 sample compounds, we've identified several areas where `calculable_features.json` can be expanded to improve feature detection accuracy and coverage. This document outlines specific enhancements needed.

---

## 🔴 Critical Issues (Fix Required)

### 1. Heteroaryl Halide Detection Pattern Too Restrictive

**Current Issue:**
```json
"smarts_any": [
  "[n,o,s]1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[Cl,Br,I,F]",
  "[n,o,s]OS(=O)(=O)"
]
```

**Problem:**
- ❌ Misses 4-bromopyridine (`Brc1ccncc1`) - halogen not adjacent to N
- ❌ Misses 3-bromopyridine (`Brc1cccnc1`)  
- ❌ Misses 2-bromothiophene (`Brc1cccs1`) - 5-membered ring
- ❌ Misses 3-bromofuran (`Brc1ccoc1`) - halogen position
- ✅ Matches 2-bromopyridine (`Brc1ccccn1`) - halogen adjacent to N

**Root Cause:** Pattern requires halogen to be directly bonded to the heteroatom in a 6-membered ring.

**Proposed Fix:**
```json
"smarts_any": [
  "c1:[n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[Cl,Br,I,F]",  
  "c1:[c,n,o,s]:[n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[Cl,Br,I,F]",
  "c1:[c,n,o,s]:[c,n,o,s]:[n,o,s]:[c,n,o,s]:[c,n,o,s]:1[Cl,Br,I,F]",
  "c1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[n,o,s]:[c,n,o,s]:1[Cl,Br,I,F]",
  "c1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[n,o,s]:1[Cl,Br,I,F]",
  "c1:[n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[Cl,Br,I,F]",
  "c1:[c,n,o,s]:[n,o,s]:[c,n,o,s]:[c,n,o,s]:1[Cl,Br,I,F]",
  "c1:[c,n,o,s]:[c,n,o,s]:[n,o,s]:[c,n,o,s]:1[Cl,Br,I,F]",
  "c1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[n,o,s]:1[Cl,Br,I,F]",
  "[n,o,s]OS(=O)(=O)"
]
```

**Simpler Alternative:**
```json
"smarts_any": [
  "[a;r5,r6][Cl,Br,I,F]",  // Any aromatic halogen in 5/6-membered ring
  "[n,o,s]OS(=O)(=O)"      // Heteroaryl sulfonate
]
```
Then use derived feature:
```json
"token": "heteroaryl_halide_present",
"derive": "(aryl_halide_present OR sp2_pseudohalide_present) AND heteroaryl_present"
```

**Impact:** Will correctly detect 14+ heteroaryl halides in sample compounds.

---

## 🟡 High Priority Additions

### 2. Distinguish Boronic Acids from Boronic Esters

**Current State:**
- Only `sp2_boron_present` exists (generic)
- `boron_bpin_present` exists but is for pinacol ester specifically

**Proposed Addition:**
```json
{
  "token": "boronic_acid_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[#6]B(O)O",
      "[#6]B([OH])[OH]"
    ]
  },
  "why": "Free boronic acid (not ester); hygroscopic, air-sensitive"
},
{
  "token": "boronic_ester_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "heuristic": "sp2_boron_present AND NOT boronic_acid_present"
  },
  "why": "Protected boronic acid; more stable than free acid"
}
```

**Impact:** Enables differentiation of storage/handling requirements.

---

### 3. Distinguish Primary, Secondary, Tertiary Amines

**Current State:**
- Only `aliphatic_amine_present` (all N-H amines)
- Only `aniline_present` (aryl-NH2/NH-R)

**Proposed Addition:**
```json
{
  "token": "primary_amine_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[NX3H2;!$([N][C,c]=O);!$([N]S(=O)(=O))]"
    ]
  },
  "why": "Primary amine (NH2); not amide/sulfonamide; highest nucleophilicity"
},
{
  "token": "secondary_amine_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[NX3H1;!$([N][C,c]=O);!$([N]S(=O)(=O))]"
    ]
  },
  "why": "Secondary amine (NHR); not amide/sulfonamide"
},
{
  "token": "tertiary_amine_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[NX3H0;!$([N][C,c]=O);!$([N]S(=O)(=O));!$(n)]"
    ]
  },
  "why": "Tertiary amine (NR3); not amide/sulfonamide; potential base/ligand"
}
```

**Impact:** Critical for Buchwald-Hartwig coupling selectivity and base selection.

---

### 4. Add Amide Detection

**Current State:** Not present

**Proposed Addition:**
```json
{
  "token": "amide_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[NX3][CX3](=[OX1])[#6]",
      "[NX3][CX3](=[OX1])[#1]"
    ]
  },
  "why": "Amide group; resonance-stabilized, low nucleophilicity"
},
{
  "token": "primary_amide_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[NX3H2][CX3](=[OX1])"
    ]
  },
  "why": "Primary amide (CONH2); potential hydrogen bonding"
},
{
  "token": "secondary_amide_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[NX3H1][CX3](=[OX1])"
    ]
  },
  "why": "Secondary amide (CONHR); reduced nucleophilicity"
}
```

**Impact:** Important for substrate compatibility assessment.

---

### 5. Add Ester Detection

**Current State:** Not present

**Proposed Addition:**
```json
{
  "token": "ester_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[#6][CX3](=[OX1])[OX2][#6]"
    ]
  },
  "why": "Ester group; potential saponification/transesterification"
},
{
  "token": "acetate_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[OX2][CX3](=[OX1])[CH3]"
    ]
  },
  "why": "Acetate ester; common protecting group"
}
```

**Impact:** Identifies potential side reactions under basic conditions.

---

### 6. Add Ketone and Aldehyde Detection

**Current State:** Not present

**Proposed Addition:**
```json
{
  "token": "carbonyl_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[CX3]=[OX1]"
    ]
  },
  "why": "Any carbonyl (ketone/aldehyde/ester/amide/acid)"
},
{
  "token": "ketone_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[#6][CX3](=[OX1])[#6]"
    ]
  },
  "why": "Ketone; potential condensation/reduction"
},
{
  "token": "aldehyde_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[CX3H1](=[OX1])[#6]",
      "[CX3H1](=[OX1])[H]"
    ]
  },
  "why": "Aldehyde; highly reactive, oxidation-prone"
}
```

**Impact:** Flags potential side reactions and stability issues.

---

### 7. Add Nitrile Detection

**Current State:** Not present

**Proposed Addition:**
```json
{
  "token": "nitrile_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[NX1]#[CX2]"
    ]
  },
  "why": "Nitrile group; potential metal coordination"
},
{
  "token": "aryl_nitrile_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "c[CX2]#[NX1]"
    ]
  },
  "why": "Aryl nitrile; electron-withdrawing substituent"
}
```

**Impact:** Important for electronic effects and coordination chemistry.

---

### 8. Add Nitro Group Detection

**Current State:** Not present

**Proposed Addition:**
```json
{
  "token": "nitro_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[NX3](=[OX1])=[OX1]"
    ]
  },
  "why": "Nitro group; strong electron-withdrawing, potential reduction"
},
{
  "token": "aryl_nitro_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "c[NX3](=[OX1])=[OX1]"
    ]
  },
  "why": "Aryl nitro; activates SNAr, can be reduced"
}
```

**Impact:** Identifies strong EWG and potential incompatibilities.

---

### 9. Add Phenol Detection

**Current State:** Only generic `alcohol_present`

**Proposed Addition:**
```json
{
  "token": "phenol_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "c[OX2H]"
    ]
  },
  "why": "Phenolic OH; more acidic than aliphatic OH"
},
{
  "token": "aliphatic_alcohol_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[CX4][OX2H]"
    ]
  },
  "why": "Aliphatic alcohol; less acidic than phenol"
}
```

**Impact:** Important for base compatibility and nucleophilicity.

---

### 10. Add Halogen Counting Features

**Current State:** Only binary (present/absent)

**Proposed Addition:**
```json
{
  "token": "halogen_count",
  "type": "int",
  "scope": "global",
  "detect": {
    "count_substructure": "[F,Cl,Br,I]"
  },
  "why": "Total number of halogens; polyhalogenation indicator"
},
{
  "token": "polyhalogenated",
  "type": "bool",
  "scope": "global",
  "detect": {
    "heuristic": "halogen_count >= 2"
  },
  "why": "Multiple halogens; potential for sequential coupling"
}
```

**Impact:** Identifies substrates for sequential/iterative coupling.

---

### 11. Add Steric Hindrance Indicators

**Current State:** Not present

**Proposed Addition:**
```json
{
  "token": "tert_butyl_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[CH0](C)(C)C"
    ]
  },
  "why": "tert-Butyl group; sterically demanding"
},
{
  "token": "isopropyl_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[CH](C)C"
    ]
  },
  "why": "Isopropyl group; moderate steric bulk"
},
{
  "token": "ortho_substitution_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "heuristic": "exists aryl ring with substituents at ortho positions"
  },
  "why": "Ortho-substituted aryl; increased steric hindrance"
}
```

**Impact:** Helps predict reaction rates and selectivity.

---

### 12. Add Protecting Group Detection

**Current State:** Only `carbamate_present`

**Proposed Addition:**
```json
{
  "token": "boc_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "N[CX3](=[OX1])[OX2]C(C)(C)C"
    ]
  },
  "why": "Boc protecting group; acid-labile"
},
{
  "token": "cbz_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "N[CX3](=[OX1])[OX2]Cc1ccccc1"
    ]
  },
  "why": "Cbz protecting group; removed by hydrogenolysis"
},
{
  "token": "fmoc_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "N[CX3](=[OX1])OCC1c2ccccc2-c2ccccc12"
    ]
  },
  "why": "Fmoc protecting group; base-labile"
},
{
  "token": "silyl_ether_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "O[Si]([#6])([#6])[#6]"
    ]
  },
  "why": "Silyl ether protecting group; fluoride-labile"
}
```

**Impact:** Identifies deprotection requirements and orthogonal strategies.

---

### 13. Add Electron-Withdrawing/Donating Group Detection

**Current State:** Not systematically covered

**Proposed Addition:**
```json
{
  "token": "strong_ewg_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "c[NX3](=[OX1])=[OX1]",  // Nitro
      "c[CX2]#[NX1]",          // Nitrile
      "c[CX3](=[OX1])[F,Cl]",  // Acyl halide
      "cS(=O)(=O)[#6]",        // Sulfonyl
      "c[CX3](=O)O"            // Carboxylic acid
    ]
  },
  "why": "Strong electron-withdrawing group on aryl ring"
},
{
  "token": "strong_edg_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "c[OX2H]",               // Phenol
      "c[OX2][CH3]",           // Methoxy
      "c[NX3H2]",              // Aniline
      "c[NX3H1]",              // N-substituted aniline
      "c[NX3]([CH3])[CH3]"     // N,N-dimethylaniline
    ]
  },
  "why": "Strong electron-donating group on aryl ring"
}
```

**Impact:** Predicts reactivity and regioselectivity.

---

### 14. Add Chelating Group Detection

**Current State:** Only `pyridine_poison_risk`

**Proposed Addition:**
```json
{
  "token": "bidentate_chelator_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "c1ccccc1-c1ccccn1",     // 2-Phenylpyridine
      "c1ccc(N)cc1N",           // o-Phenylenediamine
      "c(N)C(=O)O",             // Amino acid
      "[#6]~[#7]~[#6]~[#7]",   // Diamine backbone
      "[#6]~[#8]~[#6]~[#8]",   // Diol backbone
      "[#6]~[#7]~[#6]~[#8]"    // Amino alcohol
    ]
  },
  "why": "Potential bidentate ligand; can chelate metals"
},
{
  "token": "phosphine_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[PX3]([#6])([#6])[#6]"
    ]
  },
  "why": "Phosphine; common ligand, air-sensitive"
}
```

**Impact:** Flags potential catalyst deactivation or ligand effects.

---

## 🟢 Nice-to-Have Additions

### 15. Molecular Weight Categories

```json
{
  "token": "low_molecular_weight",
  "type": "bool",
  "scope": "global",
  "detect": {
    "heuristic": "molecular_weight < 200"
  },
  "why": "Low MW; volatile, potential loss during workup"
},
{
  "token": "high_molecular_weight",
  "type": "bool",
  "scope": "global",
  "detect": {
    "heuristic": "molecular_weight > 500"
  },
  "why": "High MW; potential solubility issues"
}
```

### 16. Ring System Complexity

```json
{
  "token": "fused_ring_system",
  "type": "bool",
  "scope": "global",
  "detect": {
    "heuristic": "has fused ring system (naphthalene, quinoline, etc.)"
  },
  "why": "Fused rings; increased rigidity, π-system extension"
},
{
  "token": "spirocyclic_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "smarts_any": [
      "[C@@]12([C,c][C,c]1)[C,c][C,c]2"
    ]
  },
  "why": "Spirocyclic center; 3D character, increased complexity"
}
```

### 17. Chiral Center Detection

```json
{
  "token": "chiral_center_present",
  "type": "bool",
  "scope": "global",
  "detect": {
    "heuristic": "has stereogenic center"
  },
  "why": "Chirality; potential racemization or stereoselective requirements"
},
{
  "token": "chiral_center_count",
  "type": "int",
  "scope": "global",
  "detect": {
    "heuristic": "count stereogenic centers"
  },
  "why": "Number of chiral centers; diastereomer complexity"
}
```

---

## 📊 Summary Statistics

### Current State (v1.1)
- **Total features defined**: 71
  - Boolean (SMARTS): 55
  - Integer (count): 2
  - Heuristic: 5
  - Derived: 14

### Proposed State (v2.0)
- **Total features**: 110+ (55% increase)
  - Boolean (SMARTS): 85+
  - Integer (count): 4+
  - Heuristic: 12+
  - Derived: 18+

### Coverage Improvements
- ✅ Heteroaryl halide detection: 14+ compounds
- ✅ Amine classification: 20+ compounds
- ✅ Carbonyl detection: 30+ compounds
- ✅ Protecting group identification: 15+ compounds
- ✅ Electronic effect prediction: 40+ compounds

---

## 🚀 Implementation Priority

### Phase 1: Critical Fixes (Week 1)
1. Fix heteroaryl_halide_present SMARTS pattern
2. Add boronic_acid_present vs boronic_ester_present
3. Add primary/secondary/tertiary amine detection

### Phase 2: High Priority (Week 2)
4. Add amide detection (primary/secondary)
5. Add ester detection
6. Add ketone/aldehyde detection
7. Add nitrile and nitro detection

### Phase 3: Medium Priority (Week 3)
8. Add phenol vs aliphatic alcohol
9. Add halogen counting
10. Add steric hindrance indicators
11. Add protecting group detection

### Phase 4: Enhancement (Week 4)
12. Add EWG/EDG detection
13. Add chelating group detection
14. Add molecular weight categories
15. Add ring system complexity
16. Add chiral center detection

---

## 📝 Testing Strategy

For each new feature:
1. ✅ Add SMARTS pattern or heuristic to JSON
2. ✅ Add unit test in `tests/test_calculable_features.py`
3. ✅ Add example compound to `tests/sample_compounds.py`
4. ✅ Run `pytest tests/test_calculable_features.py -v`
5. ✅ Run `python scripts/validate_features_simple.py`
6. ✅ Update `FEATURE_DETECTION_REPORT.md`

---

**Next Steps:**
1. Review and approve this expansion plan
2. Implement Phase 1 critical fixes
3. Validate with sample compounds
4. Proceed to subsequent phases

