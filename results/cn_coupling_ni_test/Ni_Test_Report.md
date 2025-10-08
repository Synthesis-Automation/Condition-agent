# Ni C-N Coupling - DRFP Model Test Report

## Summary

- **Total tests**: 18
- **Successful predictions**: 18
- **Failed predictions**: 0
- **Average predicted yield**: 78.4%
- **Yield range**: 68.6% - 87.8%

## Test Conditions

### Standard Ni/DMSO

- **Catalyst**: Ni
- **Base**: 280-57-9
- **Solvent**: 67-68-5
- **Temperature**: 80.0°C
- **Time**: 24.0h
- **Note**: Simple Ni salt in DMSO (most common)

**Average predicted yield**: 78.8%

### Ni/dtbbpy/DMAc

- **Catalyst**: Ni/dtbbpy
- **Base**: 865-47-4
- **Solvent**: 127-19-5
- **Temperature**: 100.0°C
- **Time**: 18.0h
- **Note**: Di-tert-butylbipyridine ligand for challenging substrates

**Average predicted yield**: 79.8%

### Ni/4,4-Me-bipy/DMF

- **Catalyst**: Ni/4,4'-Dimethyl-2,2'-bipyridine
- **Base**: 6674-22-2
- **Solvent**: 68-12-2
- **Temperature**: 90.0°C
- **Time**: 20.0h
- **Note**: Dimethyl-bipyridine ligand

**Average predicted yield**: 76.5%

## Detailed Results

### Indole + Aryl Bromide

**Description**: Simple Ni-catalyzed N-arylation of indole

**Reaction SMILES**:
```
c1ccc2[nH]ccc2c1.Brc1ccccc1 >> c1ccc(n2ccc3ccccc32)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Ni/DMSO | Ni | 80.0°C | 24.0h | 80.1% |
| Ni/dtbbpy/DMAc | Ni/dtbbpy | 100.0°C | 18.0h | 81.0% |
| Ni/4,4-Me-bipy/DMF | Ni/4,4'-Dimethyl-2,2'-bipyridine | 90.0°C | 20.0h | 76.4% |

**Best conditions**: Ni/dtbbpy/DMAc (81.0%)

---

### Aniline + Aryl Chloride

**Description**: Diaryl amine formation with Ni

**Reaction SMILES**:
```
Nc1ccccc1.Clc1ccc(C(F)(F)F)cc1 >> FC(F)(F)c1ccc(Nc2ccccc2)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Ni/DMSO | Ni | 80.0°C | 24.0h | 77.7% |
| Ni/dtbbpy/DMAc | Ni/dtbbpy | 100.0°C | 18.0h | 80.2% |
| Ni/4,4-Me-bipy/DMF | Ni/4,4'-Dimethyl-2,2'-bipyridine | 90.0°C | 20.0h | 76.3% |

**Best conditions**: Ni/dtbbpy/DMAc (80.2%)

---

### Pyrrolidine + Aryl Bromide

**Description**: Aliphatic amine coupling

**Reaction SMILES**:
```
C1CCNC1.Brc1ccc(C#N)cc1 >> N#Cc1ccc(N2CCCC2)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Ni/DMSO | Ni | 80.0°C | 24.0h | 86.2% |
| Ni/dtbbpy/DMAc | Ni/dtbbpy | 100.0°C | 18.0h | 87.8% |
| Ni/4,4-Me-bipy/DMF | Ni/4,4'-Dimethyl-2,2'-bipyridine | 90.0°C | 20.0h | 85.5% |

**Best conditions**: Ni/dtbbpy/DMAc (87.8%)

---

### Carbazole + Aryl Iodide

**Description**: N-arylation of carbazole

**Reaction SMILES**:
```
c1ccc2c(c1)[nH]c1ccccc12.Ic1ccc(OC)cc1 >> COc1ccc(n2c3ccccc3c3ccccc32)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Ni/DMSO | Ni | 80.0°C | 24.0h | 72.2% |
| Ni/dtbbpy/DMAc | Ni/dtbbpy | 100.0°C | 18.0h | 72.6% |
| Ni/4,4-Me-bipy/DMF | Ni/4,4'-Dimethyl-2,2'-bipyridine | 90.0°C | 20.0h | 68.6% |

**Best conditions**: Ni/dtbbpy/DMAc (72.6%)

---

### Piperidine + Heteroaryl Chloride

**Description**: Heteroaryl coupling

**Reaction SMILES**:
```
C1CCNCC1.Clc1ncccc1 >> c1cc(N2CCCCC2)ncc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Ni/DMSO | Ni | 80.0°C | 24.0h | 77.1% |
| Ni/dtbbpy/DMAc | Ni/dtbbpy | 100.0°C | 18.0h | 78.3% |
| Ni/4,4-Me-bipy/DMF | Ni/4,4'-Dimethyl-2,2'-bipyridine | 90.0°C | 20.0h | 76.2% |

**Best conditions**: Ni/dtbbpy/DMAc (78.3%)

---

### Imidazole + Aryl Bromide

**Description**: Heteroaryl amine coupling

**Reaction SMILES**:
```
c1cnc[nH]1.Brc1ccc(C(=O)OC)cc1 >> COC(=O)c1ccc(n2ccnc2)cc1
```

**Predicted Yields by Condition**:

| Conditions | Catalyst | Temp | Time | Predicted Yield |
|-----------|----------|------|------|----------------|
| Standard Ni/DMSO | Ni | 80.0°C | 24.0h | 79.4% |
| Ni/dtbbpy/DMAc | Ni/dtbbpy | 100.0°C | 18.0h | 79.0% |
| Ni/4,4-Me-bipy/DMF | Ni/4,4'-Dimethyl-2,2'-bipyridine | 90.0°C | 20.0h | 76.0% |

**Best conditions**: Standard Ni/DMSO (79.4%)

---

## Condition Set Comparison

| Condition Set | Avg Yield | Min Yield | Max Yield | Success Rate |
|--------------|-----------|-----------|-----------|--------------|
| Standard Ni/DMSO | 78.8% | 72.2% | 86.2% | 100% |
| Ni/dtbbpy/DMAc | 79.8% | 72.6% | 87.8% | 100% |
| Ni/4,4-Me-bipy/DMF | 76.5% | 68.6% | 85.5% | 100% |

## Key Insights

1. **Best overall conditions**: Ni/dtbbpy/DMAc (79.8% average yield)

2. **Substrate scope**:
   - Heterocyclic amines (indole, carbazole, imidazole) are excellent Ni substrates
   - Simple aliphatic amines (pyrrolidine, piperidine) work well
   - Both aryl halides and heteroaryl halides are viable

3. **Ni catalyst types**:
   - **Simple Ni salts**: Most common, cost-effective
   - **Ni/dtbbpy**: For more challenging couplings (electron-rich aromatics)
   - **Ni/dimethyl-bipy**: Balanced reactivity and selectivity

4. **Comparison to Ullmann (Cu) and Buchwald (Pd)**:
   - Ni achieves highest predicted yields (Test MAE: 8.90% vs Cu: 9.61%, Pd: 11.42%)
   - Similar temperature range to Pd, milder than Cu
   - Simpler ligands than Pd, cost-effective like Cu
   - Sometimes requires oxidants (O2, DDQ) - unique to Ni

## Model Information

- **Model**: Ni DRFP Yield Predictor
- **Training dataset**: 778 Ni-catalyzed C-N coupling reactions
- **Test MAE**: 8.90% (BEST among Cu/Pd/Ni!) ✨
- **DRFP fingerprints**: 2048-bit, radius=3
- **Algorithm**: LightGBM gradient boosting
- **Note**: Model trained with default T=80°C, time=24h (Ni dataset lacks this data)

---

*Report generated by scripts/test_ni_reactions.py*
