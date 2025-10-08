# Ullmann C-N Coupling - DRFP Model Test Report

## Summary

- **Total tests**: 24
- **Successful predictions**: 24
- **Failed predictions**: 0
- **Average predicted yield**: 74.4%
- **Yield range**: 63.6% - 83.0%

## Test Conditions

### Standard Cu/DMF

- **Catalyst**: Cu
- **Base**: 7778-53-2
- **Solvent**: 68-12-2
- **Temperature**: 110.0°C
- **Time**: 12.0h
- **Note**: Simple copper salt in DMF (most common)

**Average predicted yield**: 76.6%

### Cu/Phen/DMSO

- **Catalyst**: Cu/phen
- **Base**: 584-08-7
- **Solvent**: 67-68-5
- **Temperature**: 120.0°C
- **Time**: 18.0h
- **Note**: Phenanthroline ligand for more challenging substrates

**Average predicted yield**: 74.7%

### Cu/L-Proline/Dioxane

- **Catalyst**: Cu/L-Proline
- **Base**: 534-17-8
- **Solvent**: 123-91-1
- **Temperature**: 90.0°C
- **Time**: 24.0h
- **Note**: Amino acid ligand (Buchwald-Goldberg conditions)

**Average predicted yield**: 71.8%

## Detailed Results

### Pyrrolidine + Aryl Chloride (Simple)

**Description**: Simple Cu-catalyzed C-N coupling of pyrrolidine with aryl chloride

**Reaction SMILES**:
```
C1CCNC1.Clc1ccc(C(F)(F)F)cc1 >> FC(F)(F)c1ccc(N2CCCC2)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Cu/DMF | Cu | 110.0°C | 12.0h | 77.7% |
| Cu/Phen/DMSO | Cu/phen | 120.0°C | 18.0h | 75.4% |
| Cu/L-Proline/Dioxane | Cu/L-Proline | 90.0°C | 24.0h | 70.4% |

**Best conditions**: Standard Cu/DMF (77.7%)

---

### Piperidine + Aryl Bromide

**Description**: Aryl bromide coupling with piperidine

**Reaction SMILES**:
```
C1CCNCC1.Brc1ccccc1 >> c1ccc(N2CCCCC2)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Cu/DMF | Cu | 110.0°C | 12.0h | 83.0% |
| Cu/Phen/DMSO | Cu/phen | 120.0°C | 18.0h | 80.7% |
| Cu/L-Proline/Dioxane | Cu/L-Proline | 90.0°C | 24.0h | 78.0% |

**Best conditions**: Standard Cu/DMF (83.0%)

---

### Aniline + Aryl Iodide

**Description**: Diaryl amine formation

**Reaction SMILES**:
```
Nc1ccccc1.Ic1ccc(C#N)cc1 >> N#Cc1ccc(Nc2ccccc2)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Cu/DMF | Cu | 110.0°C | 12.0h | 68.9% |
| Cu/Phen/DMSO | Cu/phen | 120.0°C | 18.0h | 66.1% |
| Cu/L-Proline/Dioxane | Cu/L-Proline | 90.0°C | 24.0h | 63.6% |

**Best conditions**: Standard Cu/DMF (68.9%)

---

### Morpholine + Heteroaryl Chloride

**Description**: Heteroaryl coupling

**Reaction SMILES**:
```
C1COCCN1.Clc1ncccc1 >> c1cc(N2CCOCC2)ncc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Cu/DMF | Cu | 110.0°C | 12.0h | 78.4% |
| Cu/Phen/DMSO | Cu/phen | 120.0°C | 18.0h | 76.8% |
| Cu/L-Proline/Dioxane | Cu/L-Proline | 90.0°C | 24.0h | 72.6% |

**Best conditions**: Standard Cu/DMF (78.4%)

---

### Indole + Aryl Iodide

**Description**: N-arylation of indole

**Reaction SMILES**:
```
c1ccc2[nH]ccc2c1.Ic1ccccc1 >> c1ccc(n2ccc3ccccc32)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Cu/DMF | Cu | 110.0°C | 12.0h | 81.6% |
| Cu/Phen/DMSO | Cu/phen | 120.0°C | 18.0h | 81.0% |
| Cu/L-Proline/Dioxane | Cu/L-Proline | 90.0°C | 24.0h | 76.3% |

**Best conditions**: Standard Cu/DMF (81.6%)

---

### Carbazole + Aryl Bromide

**Description**: N-arylation of carbazole

**Reaction SMILES**:
```
c1ccc2c(c1)[nH]c1ccccc12.Brc1ccc(F)cc1 >> Fc1ccc(n2c3ccccc3c3ccccc32)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Cu/DMF | Cu | 110.0°C | 12.0h | 74.9% |
| Cu/Phen/DMSO | Cu/phen | 120.0°C | 18.0h | 72.8% |
| Cu/L-Proline/Dioxane | Cu/L-Proline | 90.0°C | 24.0h | 72.6% |

**Best conditions**: Standard Cu/DMF (74.9%)

---

### Imidazole + Aryl Chloride

**Description**: Heteroaryl amine coupling

**Reaction SMILES**:
```
c1cnc[nH]1.Clc1ccc(OC)cc1 >> COc1ccc(n2ccnc2)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Cu/DMF | Cu | 110.0°C | 12.0h | 75.3% |
| Cu/Phen/DMSO | Cu/phen | 120.0°C | 18.0h | 74.2% |
| Cu/L-Proline/Dioxane | Cu/L-Proline | 90.0°C | 24.0h | 70.2% |

**Best conditions**: Standard Cu/DMF (75.3%)

---

### Benzylamine + Aryl Bromide

**Description**: Primary amine coupling

**Reaction SMILES**:
```
NCc1ccccc1.Brc1ccc(C(=O)OC)cc1 >> COC(=O)c1ccc(NCc2ccccc2)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Cu/DMF | Cu | 110.0°C | 12.0h | 72.9% |
| Cu/Phen/DMSO | Cu/phen | 120.0°C | 18.0h | 70.9% |
| Cu/L-Proline/Dioxane | Cu/L-Proline | 90.0°C | 24.0h | 70.5% |

**Best conditions**: Standard Cu/DMF (72.9%)

---

## Condition Set Comparison

| Condition Set | Avg Yield | Min Yield | Max Yield | Success Rate |
|--------------|-----------|-----------|-----------|--------------|
| Standard Cu/DMF | 76.6% | 68.9% | 83.0% | 100% |
| Cu/Phen/DMSO | 74.7% | 66.1% | 81.0% | 100% |
| Cu/L-Proline/Dioxane | 71.8% | 63.6% | 78.0% | 100% |

## Key Insights

1. **Best overall conditions**: Standard Cu/DMF (76.6% average yield)

2. **Substrate scope**:
   - Aryl halide electrophiles work well (Cl, Br, I all compatible)
   - Heterocyclic amines (indole, carbazole, imidazole) are viable nucleophiles
   - Simple aliphatic amines (pyrrolidine, piperidine, morpholine) couple efficiently

3. **Cu catalyst types**:
   - **Simple Cu salts**: Most common, suitable for activated substrates
   - **Cu/phenanthroline**: For more challenging couplings (electron-rich aromatics)
   - **Cu/L-Proline**: Buchwald-Goldberg conditions, milder temperatures

4. **Comparison to Buchwald (Pd)**:
   - Ullmann (Cu) is cheaper and simpler than Pd catalysts
   - Similar predicted yields (~70-80% range)
   - Higher temperatures typically needed (90-120°C vs 80-100°C for Pd)
   - No expensive ligands like XPhos/RuPhos required

## Model Information

- **Model**: Ullmann DRFP Yield Predictor
- **Training dataset**: 4,367 Cu-catalyzed C-N coupling reactions
- **Test MAE**: 9.61% (validation MAE: 10.01%)
- **DRFP fingerprints**: 2048-bit, radius=3
- **Algorithm**: LightGBM gradient boosting

---

*Report generated by scripts/test_ullmann_reactions.py*
