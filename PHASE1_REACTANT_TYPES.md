# Phase 1 Implementation: Critical Reactant Types

This document provides the exact JSON additions needed for Phase 1 of the taxonomy expansion.

---

## 1. Carbonyl Compounds (Add to `reactant_types.json`)

### Carboxylic Acids and Derivatives

```json
{
  "id": "acyl_source",
  "description": "Carboxylic acids and activated acyl derivatives for amide/ester formation",
  "smarts": "[CX3](=O)[OX2H,OX2,Cl,$([OX2][CX3](=O))]",
  "members": [
    {
      "id": "RCO2H",
      "smarts": "[#6][CX3](=O)[OX2H]",
      "name": "carboxylic acid"
    },
    {
      "id": "RCO2M",
      "smarts": "[#6][CX3](=O)[O-].[Na+,K+,Li+]",
      "name": "carboxylate salt"
    },
    {
      "id": "RCOCl",
      "smarts": "[#6][CX3](=O)[Cl]",
      "name": "acyl chloride"
    },
    {
      "id": "RCOBr",
      "smarts": "[#6][CX3](=O)[Br]",
      "name": "acyl bromide"
    },
    {
      "id": "Anhydride",
      "smarts": "[#6][CX3](=O)[OX2][CX3](=O)[#6]",
      "name": "carboxylic anhydride"
    },
    {
      "id": "RCOOR",
      "smarts": "[#6][CX3](=O)[OX2][#6;!$(C=O)]",
      "name": "ester"
    },
    {
      "id": "activated_ester",
      "smarts": "[#6][CX3](=O)[OX2]c1ccc([N+](=O)[O-])cc1",
      "name": "activated ester (p-nitrophenyl)"
    }
  ]
}
```

### Aldehydes and Ketones

```json
{
  "id": "carbonyl",
  "description": "Aldehydes and ketones for nucleophilic addition reactions",
  "smarts": "[#6][CX3](=O)[#6,H]",
  "members": [
    {
      "id": "RCHO",
      "smarts": "[CX3H1](=O)[#6]",
      "name": "aldehyde"
    },
    {
      "id": "ArCHO",
      "smarts": "[CX3H1](=O)c",
      "name": "aromatic aldehyde"
    },
    {
      "id": "RCOR",
      "smarts": "[#6][CX3](=O)[#6;!c]",
      "name": "aliphatic ketone"
    },
    {
      "id": "ArCOR",
      "smarts": "c[CX3](=O)[#6]",
      "name": "aryl ketone"
    },
    {
      "id": "alpha_beta_unsaturated_carbonyl",
      "smarts": "[#6]=[#6][CX3](=O)[#6,H]",
      "name": "α,β-unsaturated carbonyl"
    }
  ]
}
```

### Nitriles and Amides

```json
{
  "id": "nitrogen_carbonyl",
  "description": "Nitriles, amides, and related functional groups",
  "smarts": "[CX3](=[OX1])[NX3,NX4+],C#[NX1]",
  "members": [
    {
      "id": "RCN",
      "smarts": "[#6][CX2]#[NX1]",
      "name": "nitrile"
    },
    {
      "id": "RCONH2",
      "smarts": "[#6][CX3](=O)[NX3H2]",
      "name": "primary amide"
    },
    {
      "id": "RCONHR",
      "smarts": "[#6][CX3](=O)[NX3H1][#6]",
      "name": "secondary amide"
    },
    {
      "id": "RCONR2",
      "smarts": "[#6][CX3](=O)[NX3]([#6])[#6]",
      "name": "tertiary amide"
    },
    {
      "id": "isocyanate",
      "smarts": "[NX2]=[CX2]=[OX1]",
      "name": "isocyanate"
    }
  ]
}
```

---

## 2. Organometallic Reagents (Add to `reactant_types.json`)

```json
{
  "id": "grignard_organometallic",
  "description": "Grignard and organozinc reagents for nucleophilic addition",
  "smarts": "[#6][Mg,Zn][Br,Cl,I]",
  "members": [
    {
      "id": "RMgBr",
      "smarts": "[#6][Mg][Br]",
      "name": "alkyl/aryl magnesium bromide"
    },
    {
      "id": "RMgCl",
      "smarts": "[#6][Mg][Cl]",
      "name": "alkyl/aryl magnesium chloride"
    },
    {
      "id": "RZnBr",
      "smarts": "[#6][Zn][Br]",
      "name": "organozinc bromide"
    },
    {
      "id": "RZnCl",
      "smarts": "[#6][Zn][Cl]",
      "name": "organozinc chloride"
    },
    {
      "id": "R2Zn",
      "smarts": "[#6][Zn][#6]",
      "name": "dialkylzinc"
    }
  ]
}
```

```json
{
  "id": "organolithium",
  "description": "Organolithium reagents (strong nucleophiles/bases)",
  "smarts": "[#6][Li]",
  "members": [
    {
      "id": "RLi",
      "smarts": "[#6;!c][Li]",
      "name": "alkyllithium"
    },
    {
      "id": "ArLi",
      "smarts": "c[Li]",
      "name": "aryllithium"
    },
    {
      "id": "LDA",
      "smarts": "[Li].[N-](C(C)C)C(C)C",
      "name": "lithium diisopropylamide"
    },
    {
      "id": "nBuLi",
      "smarts": "[Li]CCCC",
      "name": "n-butyllithium"
    }
  ]
}
```

---

## 3. Reducing Agents (Add to `reactant_types.json`)

```json
{
  "id": "metal_hydride",
  "description": "Metal hydride reducing agents",
  "smarts": "[BH4-,AlH4-,BH3,AlH3]",
  "members": [
    {
      "id": "NaBH4",
      "smarts": "[Na+].[BH4-]",
      "name": "sodium borohydride"
    },
    {
      "id": "LiBH4",
      "smarts": "[Li+].[BH4-]",
      "name": "lithium borohydride"
    },
    {
      "id": "LiAlH4",
      "smarts": "[Li+].[AlH4-]",
      "name": "lithium aluminum hydride"
    },
    {
      "id": "NaBH(OAc)3",
      "smarts": "[Na+].[BH]([O-]C(=O)C)([O-]C(=O)C)([O-]C(=O)C)",
      "name": "sodium triacetoxyborohydride"
    },
    {
      "id": "DIBAL",
      "smarts": "[Al](CC(C)C)(CC(C)C)[H]",
      "name": "diisobutylaluminum hydride"
    },
    {
      "id": "BH3",
      "smarts": "[BH3]",
      "name": "borane"
    },
    {
      "id": "9-BBN",
      "smarts": "B1CCCCCCCC1",
      "name": "9-borabicyclo[3.3.1]nonane"
    }
  ]
}
```

```json
{
  "id": "hydrogen_source",
  "description": "Molecular hydrogen and hydrogen donors",
  "smarts": "[H][H]",
  "members": [
    {
      "id": "H2",
      "smarts": "[H][H]",
      "name": "hydrogen gas"
    },
    {
      "id": "formic_acid_donor",
      "smarts": "[H]C(=O)[OH]",
      "name": "formic acid (H-donor)"
    },
    {
      "id": "cyclohexadiene",
      "smarts": "C1=CCC=CC1",
      "name": "1,4-cyclohexadiene (H-donor)"
    }
  ]
}
```

---

## 4. Oxidizing Agents (Add to `reactant_types.json`)

```json
{
  "id": "oxidant",
  "description": "Oxidizing reagents for alcohol/alkene oxidation",
  "smarts": "[Cr,Mn,I,N](=O)",
  "members": [
    {
      "id": "mCPBA",
      "smarts": "c1cc(Cl)cc(c1)C(=O)OO",
      "name": "meta-chloroperoxybenzoic acid"
    },
    {
      "id": "PCC",
      "smarts": "[Cr](=O)(=O)(Cl).c1ccncc1",
      "name": "pyridinium chlorochromate"
    },
    {
      "id": "DMP",
      "smarts": "c1ccc2c(c1)I(OC(=O)C)(OC(=O)C)OC2(OC(=O)C)OC(=O)C",
      "name": "Dess-Martin periodinane"
    },
    {
      "id": "IBX",
      "smarts": "c1ccc2c(c1)I(=O)(=O)OC2(O)O",
      "name": "2-iodoxybenzoic acid"
    },
    {
      "id": "TEMPO",
      "smarts": "CC1(C)CCCC(C)(C)N1[O]",
      "name": "TEMPO radical"
    },
    {
      "id": "CrO3",
      "smarts": "[Cr](=O)(=O)=O",
      "name": "chromium trioxide"
    },
    {
      "id": "KMnO4",
      "smarts": "[K+].[Mn](=O)(=O)(=O)[O-]",
      "name": "potassium permanganate"
    },
    {
      "id": "DDQ",
      "smarts": "ClC1=C(Cl)C(=O)C(C#N)=C(C#N)C1=O",
      "name": "DDQ"
    },
    {
      "id": "NaOCl",
      "smarts": "[Na+].[O-]Cl",
      "name": "sodium hypochlorite (bleach)"
    }
  ]
}
```

---

## 5. Dienes and Dienophiles (Add to `reactant_types.json`)

```json
{
  "id": "diene_dienophile",
  "description": "Dienes and dienophiles for Diels-Alder reactions",
  "smarts": "[#6]=[#6][#6]=[#6]",
  "members": [
    {
      "id": "conjugated-diene",
      "smarts": "[#6]=[#6][#6]=[#6]",
      "name": "1,3-conjugated diene"
    },
    {
      "id": "cyclopentadiene",
      "smarts": "C1=CC=CC1",
      "name": "cyclopentadiene"
    },
    {
      "id": "furan-diene",
      "smarts": "o1cccc1",
      "name": "furan (aromatic diene)"
    },
    {
      "id": "dienophile-alkene",
      "smarts": "[#6]=[#6][#6]=O",
      "name": "α,β-unsaturated carbonyl (dienophile)"
    },
    {
      "id": "dienophile-simple",
      "smarts": "[#6]=[#6]",
      "name": "simple alkene (dienophile)"
    }
  ]
}
```

---

## 6. Azides and Alkynes (for Click Chemistry)

```json
{
  "id": "click_chemistry_partners",
  "description": "Azides and alkynes for 1,3-dipolar cycloaddition",
  "smarts": "[NX1]=[NX2+]=[NX1-]",
  "members": [
    {
      "id": "RN3",
      "smarts": "[#6][NX2]=[NX2+]=[NX1-]",
      "name": "alkyl azide"
    },
    {
      "id": "ArN3",
      "smarts": "c[NX2]=[NX2+]=[NX1-]",
      "name": "aryl azide"
    },
    {
      "id": "BnN3",
      "smarts": "[#6]c[NX2]=[NX2+]=[NX1-]",
      "name": "benzyl azide"
    }
  ]
}
```

---

## 7. Ylides and Phosphonates (for Wittig/HWE)

```json
{
  "id": "ylides_phosphonates",
  "description": "Phosphorus ylides and phosphonate esters",
  "smarts": "[P+][#6-],[P](=O)([O,OC])[#6]",
  "members": [
    {
      "id": "Wittig-ylide",
      "smarts": "[P+](c1ccccc1)(c2ccccc2)(c3ccccc3)[C-]",
      "name": "phosphonium ylide (Wittig reagent)"
    },
    {
      "id": "stabilized-ylide",
      "smarts": "[P+](c1ccccc1)(c2ccccc2)(c3ccccc3)[C-]C(=O)[O,OC]",
      "name": "stabilized ylide (ester)"
    },
    {
      "id": "HWE-phosphonate",
      "smarts": "[P](=O)(OCC)(OCC)[#6]",
      "name": "phosphonate ester (HWE reagent)"
    }
  ]
}
```

---

## 8. Coupling Reagents (Special category)

```json
{
  "id": "coupling_reagent",
  "description": "Amide/ester coupling activation reagents",
  "smarts": "[#6]=[NX2+]=[NX2]",
  "members": [
    {
      "id": "EDC",
      "smarts": "CCN=C=NCCCN(C)C",
      "name": "EDC (carbodiimide)"
    },
    {
      "id": "DCC",
      "smarts": "C1CCCCC1N=C=NC2CCCCC2",
      "name": "DCC (dicyclohexylcarbodiimide)"
    },
    {
      "id": "HOBt",
      "smarts": "n1nc([OH])c2ccccc12",
      "name": "HOBt"
    },
    {
      "id": "HATU",
      "smarts": "CN(C)c1ncc[n+]([O-])c1N([P](N(C)C)(N(C)C)=N)N(C)C",
      "name": "HATU"
    }
  ]
}
```

---

## Implementation Notes

### Priority Order:

1. **Carbonyl compounds** (acyl_source, carbonyl) - Enables amide/ester detection
2. **Metal hydrides** - Enables reduction detection
3. **Organometallics** - Enables Grignard detection
4. **Oxidants** - Enables oxidation detection
5. **Dienes/dienophiles** - Enables Diels-Alder detection
6. **Azides** - Enables Click chemistry detection
7. **Ylides** - Enables Wittig detection
8. **Coupling reagents** - Helps with amide coupling confidence

### Integration Steps:

1. Add these JSON objects to `chemtools/taxonomy/data/reactant_types.json`
2. Update `manifest.json` to include new types
3. Run validation: `python -m chemtools.taxonomy.validate`
4. Test with sample reactions: `python test_sample_reactions.py`
5. Verify SMARTS patterns match expected molecules

### Expected Impact:

- **Before:** 40 reactions marked UNKNOWN
- **After Phase 1:** ~20 reactions marked UNKNOWN
- **Coverage increase:** +20% (from 60% to 80%)

---

## Next: Reaction Family Detection

After adding these reactant types, we need to add reaction family detection logic. See `PHASE1_REACTIONS.md` for:

- Amide coupling detection
- Reduction reaction detection
- Oxidation reaction detection
- Cycloaddition detection

---

**Document:** Phase 1 Reactant Types Implementation  
**Date:** 2025-10-26  
**Status:** Ready for implementation
