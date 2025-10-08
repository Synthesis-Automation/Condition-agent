# C-N Coupling Reaction Test Report

**Generated:** test_cn_coupling_reactions.py

**Model:** Buchwald DRFP v1 (LightGBM + DRFP 2048-bit fingerprints)

---

## Executive Summary

- **Total C-N Reactions Tested:** 62
- **Successful Predictions:** 62 (100.0%)
- **Average Predicted Yield:** 74.9%

### Breakdown by Reaction Type

| Reaction Type | Count | Avg Predicted Yield | Success Rate |
|--------------|-------|---------------------|-------------|
| Buchwald-Hartwig | 6 | 72.3% | 100.0% |
| C-N Coupling | 53 | 75.2% | 100.0% |
| Chan-Lam | 3 | 73.4% | 100.0% |

## Top 10 Reactions by Predicted Yield

| Rank | Predicted Yield | Description | Best Conditions |
|------|----------------|-------------|----------------|
| 1 | 90.0% | Nc2ccc3cc(Nc4ccccc4)ccc3n2)cc1 (C-N - 4,7-dichloro... | Standard Buchwald |
| 2 | 82.3% | Nc2ccccc2)cn1 (C-N - 2-Br-4-methylpyrimidine + ani... | Standard Buchwald |
| 3 | 81.8% | Nc2ccc3c(c2)OCO3)cc1 (C-N - 5-bromobenzo[d][1,3]di... | Standard Buchwald |
| 4 | 81.3% | C-N - 3-Br-furan + aniline | SPhos Conditions |
| 5 | 81.0% | C-N - Ph-Br + pyrrolidine | SPhos Conditions |
| 6 | 80.7% | c2ccccc2)CC1 (C-N - Ph-Br + N-methylpiperazine | SPhos Conditions |
| 7 | 80.1% | C-N - Ph-Br + morpholine | SPhos Conditions |
| 8 | 80.1% | =O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-NO2-Ph-Br + anilin... | Standard Buchwald |
| 9 | 79.8% | C-N - Ph-Br + tetrahydroisoquinoline | SPhos Conditions |
| 10 | 78.1% | N2CCCCC2)cc1 (C-N - 4-Br-anisole + piperidine | SPhos Conditions |

## Bottom 10 Reactions by Predicted Yield

| Rank | Predicted Yield | Description | Best Conditions |
|------|----------------|-------------|----------------|
| 1 | 65.3% | Nc2ccccc2)cc1 (C-N - 4-CN-Ph-Br + aniline | Standard Buchwald |
| 2 | 67.9% | =O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-acetyl-Ph-Br + ani... | Standard Buchwald |
| 3 | 68.7% | =O)c1ccccc1Nc1ccccc1 (C-N - 2-bromoacetophenone + ... | Standard Buchwald |
| 4 | 69.1% | C-N - 2-chloroquinoxaline + aniline | SPhos Conditions |
| 5 | 69.2% | Nc2ccccc2)cc1 (B-H - Aldehyde substrate | SPhos Conditions |
| 6 | 70.2% | C-N - 2-bromoanisole + aniline | Standard Buchwald |
| 7 | 70.3% | C-N - 3-Cl-pyridine + aniline | Standard Buchwald |
| 8 | 70.3% | C)Nc1ccccc1 (C-N - Ph-Br + isopropylamine | Standard Buchwald |
| 9 | 70.5% | Nc2ccccc2)cc1 (Chan-Lam - Methoxy boronic acid | SPhos Conditions |
| 10 | 70.7% | Nc2ccccc2)c1 (C-N - Ph-Br + 3-methylaniline | Standard Buchwald |

---

## Detailed Results

### All Tested Reactions

#### Buchwald-Hartwig

##### 1. Buchwald-Hartwig - Diphenylamine

**SMILES:** `Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.5%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 74.0% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 71.1% |

---

##### 2. Buchwald-Hartwig - Pyridine ethylamine

**SMILES:** `Clc1ccncc1.NCC>>CCNc1ccncc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **71.3%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 70.6% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 68.6% |

---

##### 3. F)(F)c1ccc(NC2CCCCC2)cc1 (B-H - Cyclohexylamine

**SMILES:** `Brc1ccc(C(F)(F)F)cc1.NC1CCCCC1>>FC`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **71.0%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 70.7% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 67.9% |

---

##### 4. Nc2ccccc2)cc1 (B-H - Aldehyde substrate

**SMILES:** `Ic1ccc(C=O)cc1.Nc1ccccc1>>O=Cc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **69.2%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 68.5% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 66.4% |

---

##### 5. B-H - Naphthylamine

**SMILES:** `Brc1ccc2ccccc2c1.NCC>>CCNc1ccc2ccccc2c1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **75.8%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 75.3% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 72.1% |

---

##### 6. B-H - Benzimidazole

**SMILES:** `Clc1nc2ccccc2[nH]1.Nc1ccccc1>>c1ccccc1Nc1nc2ccccc2[nH]1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **72.2%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 71.9% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 69.4% |

---

#### C-N Coupling

##### 1. C-N - Ph-Br + aniline → diphenylamine

**SMILES:** `Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.5%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 74.0% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 71.1% |

---

##### 2. C-N - Ph-Cl + aniline → diphenylamine

**SMILES:** `Clc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **77.0%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 76.8% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 75.4% |

---

##### 3. C-N - Ph-I + aniline → diphenylamine

**SMILES:** `Ic1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **73.7%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 72.5% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 70.6% |

---

##### 4. Nc2ccccc2)cc1 (C-N - 4-MeO-Ph-Br + aniline

**SMILES:** `Brc1ccc(OC)cc1.Nc1ccccc1>>COc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **76.5%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 75.9% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 73.5% |

---

##### 5. Nc2ccccc2)cc1 (C-N - 4-Me-Ph-Br + aniline

**SMILES:** `Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **71.6%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 71.2% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 68.0% |

---

##### 6. F)(F)c1ccc(Nc2ccccc2)cc1 (C-N - 4-CF3-Ph-Br + aniline

**SMILES:** `Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **72.6%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 72.3% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 69.1% |

---

##### 7. Nc2ccccc2)cc1 (C-N - 4-CN-Ph-Br + aniline

**SMILES:** `Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **65.3%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 64.7% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 61.7% |

---

##### 8. =O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-NO2-Ph-Br + aniline

**SMILES:** `Brc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>[O-][N+]`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **80.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 79.4% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 76.6% |

---

##### 9. =O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-acetyl-Ph-Br + aniline

**SMILES:** `Brc1ccc(C(=O)C)cc1.Nc1ccccc1>>CC`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **67.9%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 67.2% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 64.3% |

---

##### 10. Nc2ccccc2)cc1 (C-N - 4-F-Ph-Br + aniline

**SMILES:** `Brc1ccc(F)cc1.Nc1ccccc1>>Fc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **76.2%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 76.1% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 72.4% |

---

##### 11. C-N - 4-Cl-pyridine + aniline

**SMILES:** `Clc1ccncc1.Nc1ccccc1>>c1ccccc1Nc1ccncc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **77.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 76.9% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 75.4% |

---

##### 12. C-N - 4-Br-pyridine + aniline

**SMILES:** `Brc1ccncc1.Nc1ccccc1>>c1ccccc1Nc1ccncc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **77.6%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 77.4% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 75.5% |

---

##### 13. C-N - 3-Cl-pyridine + aniline

**SMILES:** `Clc1cccnc1.Nc1ccccc1>>c1ccccc1Nc1cccnc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **70.3%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 70.0% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 68.4% |

---

##### 14. C-N - 2-Br-pyrimidine + aniline

**SMILES:** `Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **73.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 72.0% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 71.2% |

---

##### 15. C-N - 4-Cl-quinoline + aniline

**SMILES:** `Clc1cccc2ncccc12.Nc1ccccc1>>c1ccccc1Nc1cccc2ncccc12`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **77.0%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 76.3% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 75.0% |

---

##### 16. C-N - 5-Br-indole + aniline

**SMILES:** `Brc1ccc2[nH]ccc2c1.Nc1ccccc1>>c1ccccc1Nc1ccc2[nH]ccc2c1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **75.6%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 74.9% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 72.4% |

---

##### 17. C-N - 3-Br-furan + aniline

**SMILES:** `Brc1ccco1.Nc1ccccc1>>c1ccccc1Nc1ccco1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **81.3%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 80.8% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 78.9% |

---

##### 18. C-N - 2-Cl-thiophene + aniline

**SMILES:** `Clc1cccs1.Nc1ccccc1>>c1ccccc1Nc1cccs1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **73.8%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 73.4% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 71.5% |

---

##### 19. Nc2ccccc2)cn1 (C-N - 2-Br-4-methylpyrimidine + aniline

**SMILES:** `Brc1cnc(C)cn1.Nc1ccccc1>>Cc1cnc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **82.3%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 81.9% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 80.1% |

---

##### 20. Nc2ccccc2)cc1 (C-N - Ph-Br + 4-methylaniline

**SMILES:** `Brc1ccccc1.Nc1ccc(C)cc1>>Cc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **75.5%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 74.6% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 71.9% |

---

##### 21. Nc2ccccc2)cc1 (C-N - Ph-Br + 4-methoxyaniline

**SMILES:** `Brc1ccccc1.Nc1ccc(OC)cc1>>COc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **73.4%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 72.9% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 69.8% |

---

##### 22. Nc2ccccc2)cc1 (C-N - Ph-Br + 4-fluoroaniline

**SMILES:** `Brc1ccccc1.Nc1ccc(F)cc1>>Fc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.7%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 74.5% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 70.9% |

---

##### 23. F)(F)c1ccc(Nc2ccccc2)cc1 (C-N - Ph-Br + 4-CF3-aniline

**SMILES:** `Brc1ccccc1.Nc1ccc(C(F)(F)F)cc1>>FC`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **71.3%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 70.9% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 67.7% |

---

##### 24. Nc2ccccc2)c1 (C-N - Ph-Br + 3-methylaniline

**SMILES:** `Brc1ccccc1.Nc1cccc(C)c1>>Cc1cccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **70.7%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 70.3% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 67.2% |

---

##### 25. C)cc(Nc2ccccc2)c1 (C-N - Ph-Br + 3,5-dimethylaniline

**SMILES:** `Brc1ccccc1.Nc1cc(C)cc(C)c1>>Cc1cc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **71.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 70.8% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 67.7% |

---

##### 26. C-N - Ph-Br + methylamine

**SMILES:** `Brc1ccccc1.CN>>CNc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.5%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 74.0% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 70.9% |

---

##### 27. C-N - Ph-Br + ethylamine

**SMILES:** `Brc1ccccc1.CCN>>CCNc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.9%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 74.5% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 71.4% |

---

##### 28. C-N - Ph-Br + propylamine

**SMILES:** `Brc1ccccc1.CCCN>>CCCNc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **77.5%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 77.0% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 73.9% |

---

##### 29. C)Nc1ccccc1 (C-N - Ph-Br + isopropylamine

**SMILES:** `Brc1ccccc1.CC(C)N>>CC`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **70.3%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 69.8% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 67.1% |

---

##### 30. C)(C)Nc1ccccc1 (C-N - Ph-Br + tert-butylamine

**SMILES:** `Brc1ccccc1.CC(C)(C)N>>CC`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 73.4% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 70.5% |

---

##### 31. C-N - Ph-Br + benzylamine

**SMILES:** `Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.7%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 74.2% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 71.0% |

---

##### 32. C-N - Ph-Br + 2-methoxyethylamine

**SMILES:** `Brc1ccccc1.NCCOC>>COCCNc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **75.9%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 75.5% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 72.4% |

---

##### 33. C-N - 4-Br-pyridine + propylamine

**SMILES:** `Brc1ccncc1.NCCC>>CCCNc1ccncc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 73.7% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 71.2% |

---

##### 34. C)Nc1ccc(C#N)cc1 (C-N - 4-Cl-benzonitrile + isopropylamine

**SMILES:** `Clc1ccc(C#N)cc1.CC(C)N>>CC`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **72.2%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 71.8% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 70.6% |

---

##### 35. CC)c1ccccc1 (C-N - Ph-Br + diethylamine

**SMILES:** `Brc1ccccc1.CCN(CC)CC>>CCN`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 73.4% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 70.5% |

---

##### 36. C-N - Ph-Br + pyrrolidine

**SMILES:** `Brc1ccccc1.N1CCCC1>>c1ccccc1N1CCCC1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **81.0%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 78.2% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 74.8% |

---

##### 37. C-N - Ph-Br + piperidine

**SMILES:** `Brc1ccccc1.N1CCCCC1>>c1ccccc1N1CCCCC1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **77.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 76.3% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 73.3% |

---

##### 38. C-N - Ph-Br + morpholine

**SMILES:** `Brc1ccccc1.N1CCOCC1>>c1ccccc1N1CCOCC1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **80.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 77.0% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 74.3% |

---

##### 39. c2ccccc2)CC1 (C-N - Ph-Br + N-methylpiperazine

**SMILES:** `Brc1ccccc1.N1CCN(C)CC1>>CN1CCN`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **80.7%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 78.8% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 75.3% |

---

##### 40. c2ccccc2)CC1 (C-N - Ph-Br + 4-hydroxypiperidine

**SMILES:** `Brc1ccccc1.N1CCC(O)CC1>>OC1CCN`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.2%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 73.5% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 69.9% |

---

##### 41. C-N - Ph-Br + tetrahydroisoquinoline

**SMILES:** `Brc1ccccc1.N1Cc2ccccc2C1>>c1ccccc1N1Cc2ccccc2C1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **79.8%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 78.3% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 75.0% |

---

##### 42. C-N - 4-Br-pyridine + pyrrolidine

**SMILES:** `Brc1ccncc1.N1CCCC1>>c1ccncc1N1CCCC1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **77.0%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 74.7% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 72.3% |

---

##### 43. F)(F)c1ccc(N2CCOCC2)cc1 (C-N - 4-Cl-benzotrifluoride + morpholine

**SMILES:** `Clc1ccc(C(F)(F)F)cc1.N1CCOCC1>>FC`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **78.0%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 77.9% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 76.7% |

---

##### 44. N2CCCCC2)cc1 (C-N - 4-Br-anisole + piperidine

**SMILES:** `Brc1ccc(OC)cc1.N1CCCCC1>>COc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **78.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 77.8% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 74.8% |

---

##### 45. Nc2ccc(Nc3ccccc3)cc2)cc1 (C-N - 1,4-dibromobenzene + aniline

**SMILES:** `Brc1ccc(Br)cc1.Nc1ccccc1>>c1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 73.7% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 70.5% |

---

##### 46. C-N - 2-methylbromobenzene + aniline

**SMILES:** `Brc1ccccc1C.Nc1ccccc1>>Cc1ccccc1Nc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **77.6%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 76.7% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 74.6% |

---

##### 47. C-N - 2-bromoanisole + aniline

**SMILES:** `Brc1ccccc1OC.Nc1ccccc1>>COc1ccccc1Nc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **70.2%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 69.6% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 67.2% |

---

##### 48. =O)c1ccccc1Nc1ccccc1 (C-N - 2-bromoacetophenone + aniline

**SMILES:** `Brc1ccccc1C(=O)C.Nc1ccccc1>>CC`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **68.7%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 67.6% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 65.3% |

---

##### 49. Nc2ccc3cc(Nc4ccccc4)ccc3n2)cc1 (C-N - 4,7-dichloroquinoline + aniline

**SMILES:** `Clc1ccc2nc(Cl)ccc2c1.Nc1ccccc1>>c1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **90.0%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 89.7% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 88.5% |

---

##### 50. Nc2ccc3c(c2)OCO3)cc1 (C-N - 5-bromobenzo[d][1,3]dioxole + aniline

**SMILES:** `Brc1ccc2c(c1)OCO2.Nc1ccccc1>>c1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **81.8%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 80.4% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 78.2% |

---

##### 51. C-N - 2-chlorobenzothiazole + aniline

**SMILES:** `Clc1nc2ccccc2s1.Nc1ccccc1>>c1ccccc1Nc1nc2ccccc2s1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **74.7%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 74.3% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 72.3% |

---

##### 52. C-N - 2-bromobenzoxazole + aniline

**SMILES:** `Brc1nc2ccccc2o1.Nc1ccccc1>>c1ccccc1Nc1nc2ccccc2o1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **71.9%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 71.5% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 69.7% |

---

##### 53. C-N - 2-chloroquinoxaline + aniline

**SMILES:** `Clc1cnc2ccccc2n1.Nc1ccccc1>>c1ccccc1Nc1cnc2ccccc2n1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **69.1%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 69.0% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 66.9% |

---

#### Chan-Lam

##### 1. Chan-Lam - Oxidative

**SMILES:** `c1ccccc1B(O)O.Nc1ccccc1.[O]>>c1ccccc1Nc1ccccc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **71.9%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 71.9% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 70.0% |

---

##### 2. Nc2ccccc2)cc1 (Chan-Lam - Methoxy boronic acid

**SMILES:** `c1ccc(OC)cc1B(O)O.Nc1ccccc1.[O]>>COc1ccc`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/SPhos
- **Base:** 1907-33-1
- **Solvent:** 123-91-1
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **70.5%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 70.4% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 68.4% |

---

##### 3. Chan-Lam - Pyridine boronic acid

**SMILES:** `c1ccncc1B(O)O.Nc1ccccc1.[O]>>c1ccccc1Nc1ccncc1`

**Recommended Conditions:**

- **Catalyst/Ligand:** Pd/XPhos
- **Base:** 1907-33-1
- **Solvent:** 108-88-3
- **Temperature:** 100.0°C
- **Time:** 12.0 hours
- **Predicted Yield:** **77.9%**

**Alternative Conditions Tested:**

| Conditions | Predicted Yield |
|-----------|----------------|
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 77.2% |
| RuPhos Conditions (Pd/RuPhos, K3PO4, Toluene) | 75.9% |

---

## Analysis and Insights

### Condition Preferences

| Condition Set | Times Recommended |
|--------------|------------------|
| Standard Buchwald (Pd/XPhos, NaOtBu, Toluene) | 44 |
| SPhos Conditions (Pd/SPhos, NaOtBu, Dioxane) | 18 |

### Yield Distribution

- **Min:** 65.3%
- **Max:** 90.0%
- **Mean:** 74.9%
- **Median:** 74.5%

**Yield Brackets:**

- **High (≥80%):** 8 reactions (12.9%)
- **Medium (60-79%):** 54 reactions (87.1%)
- **Low (<60%):** 0 reactions (0.0%)

## Conclusions

1. **Model Performance:** The Buchwald DRFP model successfully predicted yields for 62/62 C-N coupling reactions.

2. **Average Yield:** Predicted average yield of 74.9% indicates good prospects

3. **Condition Diversity:** Model tested multiple condition sets (XPhos, RuPhos, SPhos) to find optimal conditions for each reaction.

4. **Next Steps:**
   - Validate top predictions experimentally
   - Compare predicted vs actual yields
   - Expand model training data with additional C-N reaction types
   - Integrate with automated synthesis workflows

