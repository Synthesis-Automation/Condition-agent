# ChemTools Reaction Type Detection Methods

## Overview

ChemTools provides **multiple methods** for detecting and classifying reaction types, ranging from simple rule-based SMARTS pattern matching to advanced ML-based classification. This document catalogs all available detection methods.

---

## 🔍 Detection Methods Summary

| Method | Module | Type | Accuracy | Speed | Dependencies |
|--------|--------|------|----------|-------|--------------|
| **Rule-Based SMARTS** | `router.detect_family()` | Deterministic | High | Fast | RDKit (optional) |
| **Full Reaction Analysis** | `router.detect_family_from_reaction()` | Hybrid | Very High | Fast | RDKit (optional) |
| **ML-Based (rxn-insight)** | `reaction_type_detector.detect_reaction_type()` | ML | Very High | Medium | rxn-insight |
| **Reactant Classification** | `featurizers.analysis.reactants` | Taxonomy | High | Fast | RDKit |
| **Context-Aware Classification** | `featurizers.analysis.reaction_context` | Advanced | Very High | Medium | RDKit |

---

## 1. Rule-Based SMARTS Detection

### `chemtools.router.detect_family(reactants: List[str])`

**Primary deterministic method** using SMARTS patterns to detect reaction families.

#### Features
- ✅ Fast and deterministic
- ✅ No ML dependencies
- ✅ Works offline
- ✅ 30+ reaction families supported
- ✅ Priority-based family assignment

#### Supported Reactions

**Priority 1: Organometallic Addition to Carbonyl**
- Grignard addition
- Organolithium addition
- Organozinc addition

**Priority 2: C-C Couplings**
- Suzuki-Miyaura coupling
- Sonogashira coupling
- Kumada coupling
- Negishi coupling
- Heck coupling

**Priority 3: Heteroatom Couplings**
- C-N coupling (unified Buchwald/Ullmann)
- C-O coupling (Ullmann-type)
- C-S coupling (Ullmann-type)

**Priority 4: Amide/Ester Formation**
- Amide coupling
- Esterification

**Priority 5: SN2 & Substitutions**
- Hydroboration
- Nitrile formation
- Finkelstein reaction
- Williamson ether synthesis

**Phase 2 Additions:**
- Hydrogenation
- Carbonyl reduction
- Alcohol oxidation
- Epoxidation
- E2 elimination

#### SMARTS Patterns Used

```python
{
    "aryl_halide": "[$(c[Cl,Br,I]),$(c-[Cl,Br,I])]",
    "vinyl_halide": "C=C[Cl,Br,I]",
    "triflate": "OS(=O)(=O)C(F)(F)F",
    "boron": "[BX3;$(B(O)O),$(B(O)O),$(B(O)O)]",
    "terminal_alkyne": "C#C[H]",
    "acid": "C(=O)[OH]",
    "nucleophile_n": "[NX3;H1,H2]",
    "nucleophile_o": "[OX2H]",
    "nucleophile_s": "[SX2H]",
    "carbonyl": "[CX3]=O",
    "aldehyde": "[CX3H](=O)",
    "ketone": "[CX3](=O)[C]",
    "ester": "[CX3](=O)[OX2][C,H]",
    "alcohol": "[CX4][OX2H]",
    "grignard": "[C,c][Mg][Br,Cl,I]",
    "organozinc": "[C,c][Zn][Br,Cl,I]",
    "organolithium": "[C,c][Li]",
    "cyanide": "[C-]#N",
    "alkyl_halide": "[CX4][Cl,Br,I]",
    "alkene": "C=C",
    "borane": "[BH3,BH2,BH]",
    # ... and more
}
```

#### Usage Example

```python
from chemtools.router import detect_family

reactants = ["Brc1ccccc1", "c1ccc(B(O)O)cc1"]
result = detect_family(reactants)

# Result:
{
    "family": "suzuki_miyaura",
    "confidence": 0.9,
    "hits": {
        "aryl_halide": True,
        "boron": True,
        "nucleophile_n": False,
        ...
    }
}
```

---

## 2. Full Reaction Analysis (Recommended)

### `chemtools.router.detect_family_from_reaction(reaction_smiles: str, use_rxn_insight: bool = True)`

**Main entry point** - combines rule-based + ML detection with catalyst analysis.

#### Features
- ✅ Accepts full reaction SMILES
- ✅ Normalizes reaction first
- ✅ Detects catalyst metals from agents
- ✅ Applies catalyst-based overrides (Pd→Buchwald, Cu→Ullmann)
- ✅ Falls back to rule-based if ML unavailable
- ✅ Returns comprehensive detection metadata

#### Detection Flow

```
Input: reaction SMILES
    ↓
1. Normalize reaction (chemtools.smiles)
    ↓
2. Extract reactants and agents
    ↓
3. Rule-based detection (detect_family)
    ↓
4. Detect catalyst metals (Pd, Cu, Ni, Co)
    ↓
5. Apply catalyst overrides
    ↓
6. Optional: ML detection (rxn-insight)
    ↓
7. Merge and rank results
    ↓
Output: Final family + confidence
```

#### Catalyst Override Logic

For **C-N coupling** reactions:
- `Pd` catalyst → `buchwald_cn` (conf: 0.95)
- `Cu` catalyst → `ullmann_cn` (conf: 0.90)
- No catalyst or Ni/Co → `cn_coupling` (generic)

#### Usage Example

```python
from chemtools.router import detect_family_from_reaction

reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

# With ML (default)
result = detect_family_from_reaction(reaction, use_rxn_insight=True)

# Rule-based only
result = detect_family_from_reaction(reaction, use_rxn_insight=False)

# Result:
{
    "family": "buchwald_cn",  # or "ullmann_cn" based on catalyst
    "confidence": 0.95,
    "hits": {
        "aryl_halide": True,
        "nucleophile_n": True,
        "catalyst_pd": True,
        ...
    },
    "auto": {  # ML results if available
        "rxn_class": "C-N Coupling",
        "rxn_name": "Buchwald-Hartwig amination",
        "mapped_family": "buchwald_cn",
        "confidence": 0.92
    }
}
```

---

## 3. ML-Based Detection (rxn-insight)

### `chemtools.reaction_type_detector.detect_reaction_type(reaction_smiles: str)`

**Advanced ML method** using the optional `rxn-insight` package.

#### Features
- ✅ Deep learning-based classification
- ✅ Trained on large reaction database
- ✅ Provides reaction class and specific name
- ✅ Automatic mapping to ChemTools taxonomy
- ✅ Graceful fallback if unavailable

#### Requirements
```bash
pip install rxn-insight
```

#### Detection Capabilities

The ML model can detect:
- Reaction class (broad category)
- Reaction name (specific type)
- Confidence score
- Mapped to canonical taxonomy

#### Usage Example

```python
from chemtools.reaction_type_detector import detect_reaction_type, is_available

if is_available():
    result = detect_reaction_type("Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1")
    
    # Result:
    {
        "available": True,
        "success": True,
        "rxn_class": "C-C Coupling",
        "rxn_name": "Suzuki coupling with boronic acids",
        "mapped_family": "suzuki_miyaura",
        "confidence": 0.92,
        "raw": {...},  # Raw rxn-insight output
        "catalysts": ["Pd"]  # Detected from agents
    }
```

#### Family Mapping

The ML detector maps rxn-insight classes to ChemTools families:

| rxn-insight Class | ChemTools Family |
|-------------------|------------------|
| "C-C Coupling" + boron | `suzuki_miyaura` |
| "C-N Coupling" + Pd | `buchwald_cn` |
| "C-N Coupling" + Cu | `ullmann_cn` |
| "Heteroatom Alkylation" | `ullmann_cn` |
| "Acylation" | `amide_coupling` |
| ... | ... |

---

## 4. Reactant Classification

### `chemtools.featurizers.analysis.reactants.classify_reactant_smiles(smiles: str)`

**Taxonomy-based reactant classification** for identifying reactant types.

#### Features
- ✅ Classifies individual reactants
- ✅ Multi-category matching
- ✅ Hierarchical taxonomy (category → group → specific)
- ✅ No RDKit required (uses taxonomy definitions)

#### Reactant Categories

- `aryl_halide`
- `amine`
- `boronic_acid`
- `alkyne`
- `carboxylic_acid`
- `alcohol`
- `grignard`
- `organozinc`
- `aldehyde`
- `ketone`
- ... and more

#### Usage Example

```python
from chemtools.reagent import classify_reactant_smiles

result = classify_reactant_smiles("Brc1ccccc1")

# Result:
{
    "category": "aryl_halide",
    "group": "electrophile",
    "matches": [
        {
            "category": "aryl_halide",
            "confidence": 0.95,
            "method": "smarts"
        }
    ]
}
```

---

## 5. Context-Aware Classification

### `chemtools.featurizers.analysis.reaction_context.classify_reactants_with_context(reaction_smiles: str)`

**Advanced method** combining SMILES normalization, reactant classification, and reaction family detection.

#### Features
- ✅ Full reaction context analysis
- ✅ Normalized SMILES for each component
- ✅ Reactant role assignment
- ✅ Family detection with confidence
- ✅ Comprehensive metadata

#### Usage Example

```python
from chemtools.featurizers.analysis.reaction_context import classify_reactants_with_context

result = classify_reactants_with_context("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")

# Result:
{
    "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "normalized": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "reactants": [
        {
            "smiles": "Brc1ccccc1",
            "normalized": "Brc1ccccc1",
            "category": "aryl_halide",
            "role": "electrophile"
        },
        {
            "smiles": "Nc1ccccc1",
            "normalized": "Nc1ccccc1",
            "category": "amine",
            "role": "nucleophile"
        }
    ],
    "products": [...],
    "family": "buchwald_cn",
    "confidence": 0.95,
    "detection_method": "catalyst_override"
}
```

---

## 🎯 Which Method Should I Use?

### Quick Decision Guide

**For most use cases:**
```python
from chemtools.router import detect_family_from_reaction

result = detect_family_from_reaction(reaction_smiles)
```

**For reactant analysis only:**
```python
from chemtools.router import detect_family

result = detect_family(["Brc1ccccc1", "Nc1ccccc1"])
```

**For ML-based detection:**
```python
from chemtools.reaction_type_detector import detect_reaction_type

result = detect_reaction_type(reaction_smiles)
```

**For complete context analysis:**
```python
from chemtools.featurizers.analysis.reaction_context import classify_reactants_with_context

result = classify_reactants_with_context(reaction_smiles)
```

**For individual reactants:**
```python
from chemtools.reagent import classify_reactant_smiles

result = classify_reactant_smiles(smiles)
```

---

## 📊 Method Comparison

| Use Case | Recommended Method | Why |
|----------|-------------------|-----|
| **Production API** | `detect_family_from_reaction()` | Best accuracy, hybrid approach |
| **Offline/No network** | `detect_family_from_reaction(use_rxn_insight=False)` | Rule-based, no ML needed |
| **Reactant lists only** | `detect_family()` | Fast, no normalization needed |
| **Highest accuracy** | `detect_family_from_reaction(use_rxn_insight=True)` | Combines all methods |
| **Individual reactants** | `classify_reactant_smiles()` | Taxonomy-based classification |
| **Full analysis** | `classify_reactants_with_context()` | Comprehensive metadata |

---

## 🔧 Integration Examples

### FastAPI Endpoint

```python
from chemtools.router import detect_family_from_reaction

@app.post("/api/v1/router/detect-family")
def api_router_detect(reaction: str):
    return detect_family_from_reaction(reaction)
```

### LangChain Tool

```python
from lang_chain.chemtools_wrapper import detect_reaction_family_tool

result = detect_reaction_family_tool.invoke({
    "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
})
```

### CLI

```python
from chemtools import chem

result = chem.router.detect_family(reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")
```

---

## 🎓 Advanced Features

### Catalyst Metal Detection

All methods detect catalyst metals from reaction agents:

```python
result = detect_family_from_reaction("Brc1ccccc1.Nc1ccccc1>Pd(OAc)2.BINAP>c1ccccc1Nc1ccccc1")

# Detects Pd from agents, overrides to buchwald_cn
result["hits"]["catalyst_pd"]  # True
result["family"]  # "buchwald_cn"
```

### Priority-Based Family Assignment

The rule-based detector uses priority levels to handle ambiguous cases:

1. **Highest**: Organometallic addition to carbonyl
2. **High**: C-C couplings
3. **Medium**: Heteroatom couplings
4. **Low**: Amide/ester formation
5. **Lowest**: SN2 reactions

### Confidence Scoring

All methods return confidence scores:

- `0.90-1.00`: Very high confidence
- `0.80-0.89`: High confidence
- `0.70-0.79`: Medium confidence
- `0.50-0.69`: Low confidence
- `<0.50`: Very low confidence

---

## 📚 Related Documentation

- [Router Module Documentation](./CHEMTOOLS_INTERFACE_OVERVIEW.md#routing-utilities)
- [Reaction Type Detector](../chemtools/reaction_type_detector.py)
- [Analysis Module](../chemtools/analysis/)
- [API Documentation](./API_DOCUMENTATION.md)

---

## 🛠️ Troubleshooting

### rxn-insight not available

If ML detection fails:
```python
result = detect_family_from_reaction(reaction, use_rxn_insight=False)
```

### RDKit not available

SMARTS patterns gracefully degrade to text-based matching.

### Low confidence results

1. Check SMILES validity
2. Ensure reactants are properly separated
3. Include agents in reaction SMILES
4. Try ML method if available

---

**Last Updated**: October 29, 2025
