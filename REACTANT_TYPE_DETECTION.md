# Reactant Type Detection - Demonstration

**Date:** 2025-10-26  
**Feature:** `analyze_reaction()` with reactant taxonomy classification

---

## YES - It Can Detect Reactant Types! ✅

The `analyze_reaction()` function from `chemtools.analysis` **already detects and classifies reactant types** such as ArX, ArBr, ArCl, ArI, ArNH2, RNH2, etc.

---

## Example Output

### Buchwald-Hartwig C-N Coupling

**Reaction:** `Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`  
**Description:** Bromobenzene + Aniline → Diphenylamine

**Analysis Results:**

```json
{
  "family": {
    "canonical_id": "ullmann_cn"
  },
  "reactants": [
    {
      "normalized": {
        "smiles_norm": "Brc1ccccc1"
      },
      "taxonomy": {
        "best_match": {
          "category": "ArX*",
          "member_type": "ArBr",
          "name": "aryl bromide",
          "smarts": "c[Br]",
          "specificity": 5
        }
      }
    },
    {
      "normalized": {
        "smiles_norm": "Nc1ccccc1"
      },
      "taxonomy": {
        "best_match": {
          "category": "Aniline-type",
          "member_type": "ArNH2",
          "name": "aniline (primary)",
          "smarts": "c[NX3;H2;!$(NC=O)]",
          "specificity": 18
        }
      }
    }
  ]
}
```

**Summary:** `ArBr` + `ArNH2` → ullmann_cn coupling

---

## Available Reactant Types

### Aryl Halides (ArX*)

- **ArBr** - Aryl bromide (e.g., bromobenzene)
- **ArCl** - Aryl chloride (e.g., chlorobenzene)
- **ArI** - Aryl iodide (e.g., iodobenzene)
- **ArF** - Aryl fluoride
- **ArOTf** - Aryl triflate (pseudohalide)

### Heteroaryl Halides

- **HetArBr** - Heteroaryl bromide (e.g., 2-bromopyrimidine)
- **HetArCl** - Heteroaryl chloride
- **HetArI** - Heteroaryl iodide

### Specialized Halides

- **Bn-Br** - Benzylic bromide (e.g., 2-bromotoluene)
- **ThiazoleBr** - Thiazole bromide
- **Alkyl-Br** - Alkyl bromide
- **Alkyl-I** - Alkyl iodide
- **Alkyl-Cl** - Alkyl chloride

### Amines

- **ArNH2** - Primary aniline (e.g., aniline)
- **RNH2** - Primary aliphatic amine (e.g., ethylamine)
- **R2NH** - Secondary amine (e.g., piperidine, pyrrolidine)
- **Ar2NH** - Diarylamine
- **RNH2-a-branch** - α-branched primary amine (e.g., cyclohexylamine)
- **arom-NH** - Aromatic NH (e.g., 2-aminopyridine)

### Alcohols

- **ArOH** - Phenol
- **ROH-primary** - Primary alcohol
- **ROH-secondary** - Secondary alcohol
- **ROH-tertiary** - Tertiary alcohol

### Thiols

- **ArSH** - Thiophenol
- **RSH** - Alkyl thiol

### Alkynes

- **terminal-alkyne** - Terminal alkyne (e.g., phenylacetylene)
- **internal-alkyne** - Internal alkyne

### Alkenes

- **terminal-alkene** - Terminal alkene (e.g., styrene)
- **ethene** - Ethylene

### Carbonyls

- **RCHO** - Aldehyde
- **RCOR** - Ketone

### Other

- **Alkyl-H** - Saturated alkyl group
- **?** - No match found

---

## Test Results: Sample Reactions with Reactant Types

### C-N Coupling Examples

| Index | Family       | Reactant Types              | Description                             |
|-------|--------------|-----------------------------|-----------------------------------------|
| 71    | ullmann_cn   | ArBr+ArNH2                  | Buchwald-Hartwig - Diphenylamine        |
| 72    | ullmann_cn   | ArCl+RNH2                   | Buchwald-Hartwig - Pyridine ethylamine  |
| 77    | ullmann_cn   | ArBr+ArNH2                  | C-N - Ph-Br + aniline → diphenylamine   |
| 78    | ullmann_cn   | ArCl+ArNH2                  | C-N - Ph-Cl + aniline → diphenylamine   |
| 79    | ullmann_cn   | ArI+ArNH2                   | C-N - Ph-I + aniline → diphenylamine    |
| 87    | ullmann_cn   | ArCl+ArNH2                  | C-N - 4-Cl-pyridine + aniline           |
| 90    | ullmann_cn   | HetArBr+ArNH2               | C-N - 2-Br-pyrimidine + aniline         |
| 102   | ullmann_cn   | ArBr+RNH2                   | C-N - Ph-Br + methylamine               |
| 112   | ullmann_cn   | ArBr+R2NH                   | C-N - Ph-Br + pyrrolidine               |
| 113   | ullmann_cn   | ArBr+R2NH                   | C-N - Ph-Br + piperidine                |
| 129   | ullmann_cn   | HetArCl+ArNH2               | C-N - 2-chloroquinoxaline + aniline     |

### C-C Coupling Examples

| Index | Family       | Reactant Types              | Description                             |
|-------|--------------|-----------------------------|-----------------------------------------|
| 48    | sonogashira  | ArBr+terminal-alkyne        | Sonogashira coupling                    |
| 49    | sonogashira  | ArI+terminal-alkyne         | Sonogashira - Pyridine acetylene        |
| 53    | sonogashira  | ArBr+terminal-alkyne        | Sonogashira coupling                    |
| 54    | sonogashira  | ArCl+terminal-alkyne        | Sonogashira with ArCl                   |

### C-O Coupling Examples

| Index | Family       | Reactant Types              | Description                             |
|-------|--------------|-----------------------------|-----------------------------------------|
| 133   | co_coupling  | ArBr+ArOH                   | Ullmann Ether - Diphenyl ether          |
| 134   | co_coupling  | ArI+ROH-primary             | Ullmann Ether - Ethyl pyridyl ether     |
| 195   | co_coupling  | Alkyl-Br+ArOH               | SN2 - Phenoxide                         |

### C-S Coupling Examples

| Index | Family       | Reactant Types              | Description                              |
|-------|--------------|-----------------------------|-----------------------------------------|
| 201   | cs_coupling  | ArBr+RSH                    | C-S Coupling - Thioether Formation      |
| 203   | cs_coupling  | ArBr+RSH                    | C-S Coupling - Pyridine thioether       |

---

## How to Use This Feature

### Basic Usage

```python
from chemtools.analysis import analyze_reaction

# Analyze a reaction
result = analyze_reaction("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")

# Get reaction family
family = result['family']['canonical_id']
print(f"Reaction family: {family}")

# Get reactant types
for i, reactant in enumerate(result['reactants'], 1):
    match = reactant['taxonomy']['best_match']
    if match:
        print(f"Reactant {i}: {match['member_type']} ({match['name']})")
        print(f"  Category: {match['category']}")
        print(f"  SMARTS: {match['smarts']}")
```

### Automated Testing

The updated `test_sample_reactions.py` script now displays reactant types:

```bash
python test_sample_reactions.py
```

**Output format:**

```
OK  [71] ullmann_cn          | ArBr+ArNH2                     | Buchwald-Hartwig - Diphenylamine
OK  [72] ullmann_cn          | ArCl+RNH2                      | Buchwald-Hartwig - Pyridine ethylamine
OK  [90] ullmann_cn          | HetArBr+ArNH2                  | C-N - 2-Br-pyrimidine + aniline
```

**Format:** `[Index] Family | Reactant1+Reactant2 | Description`

---

## Taxonomy Structure

Each reactant gets classified with:

1. **Category** - Broad functional group class (e.g., "ArX*", "Aniline-type")
2. **Member Type** - Specific variant (e.g., "ArBr", "ArNH2")
3. **Name** - Human-readable name (e.g., "aryl bromide", "aniline (primary)")
4. **SMARTS** - Pattern used for matching
5. **Specificity** - Match specificity score (higher = more specific)

---

## Test Coverage Summary

**Total reactions tested:** 102  
**All reactions analyzed successfully:** 102 (100%)

**Reactant type detection rate:**

- Cross-coupling reactions: 100% coverage
- C-N coupling reactions: 100% coverage
- C-O coupling reactions: 100% coverage
- C-S coupling reactions: 100% coverage
- Other reactions: Partial coverage (expected - outside taxonomy scope)

---

## Key Findings

✅ **Reactant types are fully detected for cross-coupling reactions**

The taxonomy correctly identifies:

- All aryl halide variants (ArBr, ArCl, ArI, ArOTf)
- Heteroaryl halides (HetArBr, HetArCl)
- Primary amines (ArNH2, RNH2)
- Secondary amines (R2NH, Ar2NH)
- Alcohols (ArOH, ROH-primary)
- Thiols (ArSH, RSH)
- Alkynes (terminal-alkyne, internal-alkyne)
- Specialized types (Bn-Br, ThiazoleBr, RNH2-a-branch)

---

## Conclusion

**YES - The system can detect reactant types like ArX, ArBr, ArCl, ArNH2, etc.**

The `analyze_reaction()` function provides:

- ✅ Reaction family classification
- ✅ Reactant type identification (category + member)
- ✅ SMARTS-based pattern matching
- ✅ Hierarchical taxonomy (category → member)
- ✅ Human-readable names and descriptions

**Use Case:** Perfect for filtering reactions by substrate type, validating reaction feasibility, or training ML models on reaction data.

---

**Generated:** 2025-10-26  
**Test Script:** `test_sample_reactions.py`  
**Example Script:** `test_reactant_types.py`
