# ML vs Rule-Based Verification Report

**Purpose:** Validate ML predictions against precedent-based recommendations

---

## Executive Summary

- **Total Reactions Verified:** 59
- **Fully Verified:** 0 (0.0%)
- **Minor Discrepancies:** 21 (35.6%)
- **Major Discrepancies:** 36 (61.0%)
- **Unverifiable:** 2 (3.4%)

## Common Discrepancy Types

| Discrepancy Type | Count | % of Total |
|-----------------|-------|------------|
| Base Mismatch | 57 | 96.6% |
| Core Mismatch | 50 | 84.7% |
| Yield Large Diff | 26 | 44.1% |
| Solvent Mismatch | 21 | 35.6% |
| No Precedents | 2 | 3.4% |

## ⚠⚠ Major Discrepancies (Requires Review)

### Buchwald-Hartwig - Diphenylamine

**SMILES:** `Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 74.5%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### Buchwald-Hartwig - Pyridine ethylamine

**SMILES:** `Clc1ccncc1.NCC>>CCNc1ccncc1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 71.3%

**Top Precedent:**
- Core: Pd/XPhos
- Actual Yield: 91.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### Nc2ccccc2)cc1 (B-H - Aldehyde substrate

**SMILES:** `Ic1ccc(C=O)cc1.Nc1ccccc1>>O=Cc1ccc`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 69.2%

**Top Precedent:**
- Core: Pd/XantPhos
- Actual Yield: 90.0%
- Precedent Count: 1

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

### B-H - Benzimidazole

**SMILES:** `Clc1nc2ccccc2[nH]1.Nc1ccccc1>>c1ccccc1Nc1nc2ccccc2[nH]1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 72.2%

**Top Precedent:**
- Core: Pd/XPhos
- Actual Yield: 91.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### C-N - Ph-Br + aniline → diphenylamine

**SMILES:** `Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 74.5%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### C-N - Ph-Cl + aniline → diphenylamine

**SMILES:** `Clc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 77.0%

**Top Precedent:**
- Core: Pd/CyclohexylJohnPhos
- Actual Yield: 99.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

### C-N - Ph-I + aniline → diphenylamine

**SMILES:** `Ic1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 73.7%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

### Nc2ccccc2)cc1 (C-N - 4-CN-Ph-Br + aniline

**SMILES:** `Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 65.3%

**Top Precedent:**
- Core: Pd/XantPhos
- Actual Yield: 90.0%
- Precedent Count: 1

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### =O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-acetyl-Ph-Br + aniline

**SMILES:** `Brc1ccc(C(=O)C)cc1.Nc1ccccc1>>CC`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 67.9%

**Top Precedent:**
- Core: Pd/XantPhos
- Actual Yield: 90.0%
- Precedent Count: 1

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### C-N - 3-Br-furan + aniline

**SMILES:** `Brc1ccco1.Nc1ccccc1>>c1ccccc1Nc1ccco1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 81.3%

**Top Precedent:**
- Core: Pd/XantPhos
- Actual Yield: 90.0%
- Precedent Count: 1

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### C-N - 2-Cl-thiophene + aniline

**SMILES:** `Clc1cccs1.Nc1ccccc1>>c1ccccc1Nc1cccs1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 73.8%

**Top Precedent:**
- Core: Pd/XPhos
- Actual Yield: 91.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### Nc2ccccc2)cc1 (C-N - Ph-Br + 4-methylaniline

**SMILES:** `Brc1ccccc1.Nc1ccc(C)cc1>>Cc1ccc`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 75.5%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### Nc2ccccc2)cc1 (C-N - Ph-Br + 4-methoxyaniline

**SMILES:** `Brc1ccccc1.Nc1ccc(OC)cc1>>COc1ccc`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 73.4%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### Nc2ccccc2)cc1 (C-N - Ph-Br + 4-fluoroaniline

**SMILES:** `Brc1ccccc1.Nc1ccc(F)cc1>>Fc1ccc`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 74.7%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### F)(F)c1ccc(Nc2ccccc2)cc1 (C-N - Ph-Br + 4-CF3-aniline

**SMILES:** `Brc1ccccc1.Nc1ccc(C(F)(F)F)cc1>>FC`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 71.3%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### Nc2ccccc2)c1 (C-N - Ph-Br + 3-methylaniline

**SMILES:** `Brc1ccccc1.Nc1cccc(C)c1>>Cc1cccc`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 70.7%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### C)cc(Nc2ccccc2)c1 (C-N - Ph-Br + 3,5-dimethylaniline

**SMILES:** `Brc1ccccc1.Nc1cc(C)cc(C)c1>>Cc1cc`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 71.1%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### C-N - Ph-Br + methylamine

**SMILES:** `Brc1ccccc1.CN>>CNc1ccccc1`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 74.5%

**Top Precedent:**
- Core: Pd/CyclohexylJohnPhos
- Actual Yield: 99.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

### C-N - Ph-Br + ethylamine

**SMILES:** `Brc1ccccc1.CCN>>CCNc1ccccc1`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 74.9%

**Top Precedent:**
- Core: Pd/CyclohexylJohnPhos
- Actual Yield: 99.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

### C-N - Ph-Br + propylamine

**SMILES:** `Brc1ccccc1.CCCN>>CCCNc1ccccc1`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 77.5%

**Top Precedent:**
- Core: Pd/CyclohexylJohnPhos
- Actual Yield: 99.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

### C)Nc1ccccc1 (C-N - Ph-Br + isopropylamine

**SMILES:** `Brc1ccccc1.CC(C)N>>CC`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 70.3%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### C)(C)Nc1ccccc1 (C-N - Ph-Br + tert-butylamine

**SMILES:** `Brc1ccccc1.CC(C)(C)N>>CC`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 74.1%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### C-N - Ph-Br + benzylamine

**SMILES:** `Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 74.7%

**Top Precedent:**
- Core: Pd/CyclohexylJohnPhos
- Actual Yield: 99.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

### C-N - Ph-Br + 2-methoxyethylamine

**SMILES:** `Brc1ccccc1.NCCOC>>COCCNc1ccccc1`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 75.9%

**Top Precedent:**
- Core: Pd/CyclohexylJohnPhos
- Actual Yield: 99.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

### C-N - 4-Br-pyridine + propylamine

**SMILES:** `Brc1ccncc1.NCCC>>CCCNc1ccncc1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 74.1%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 81.0%
- Precedent Count: 1

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### CC)c1ccccc1 (C-N - Ph-Br + diethylamine

**SMILES:** `Brc1ccccc1.CCN(CC)CC>>CCN`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 74.1%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### C-N - Ph-Br + pyrrolidine

**SMILES:** `Brc1ccccc1.N1CCCC1>>c1ccccc1N1CCCC1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 81.0%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### C-N - Ph-Br + piperidine

**SMILES:** `Brc1ccccc1.N1CCCCC1>>c1ccccc1N1CCCCC1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 77.1%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### C-N - Ph-Br + morpholine

**SMILES:** `Brc1ccccc1.N1CCOCC1>>c1ccccc1N1CCOCC1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 80.1%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### c2ccccc2)CC1 (C-N - Ph-Br + N-methylpiperazine

**SMILES:** `Brc1ccccc1.N1CCN(C)CC1>>CN1CCN`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 80.7%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### c2ccccc2)CC1 (C-N - Ph-Br + 4-hydroxypiperidine

**SMILES:** `Brc1ccccc1.N1CCC(O)CC1>>OC1CCN`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 74.2%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

### C-N - Ph-Br + tetrahydroisoquinoline

**SMILES:** `Brc1ccccc1.N1Cc2ccccc2C1>>c1ccccc1N1Cc2ccccc2C1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 79.8%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### C-N - 4-Br-pyridine + pyrrolidine

**SMILES:** `Brc1ccncc1.N1CCCC1>>c1ccncc1N1CCCC1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 77.0%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 81.0%
- Precedent Count: 1

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch

---

### C-N - 2-bromoanisole + aniline

**SMILES:** `Brc1ccccc1OC.Nc1ccccc1>>COc1ccccc1Nc1ccccc1`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 70.2%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### =O)c1ccccc1Nc1ccccc1 (C-N - 2-bromoacetophenone + aniline

**SMILES:** `Brc1ccccc1C(=O)C.Nc1ccccc1>>CC`

**ML Prediction:**
- Core: Pd/XPhos
- Predicted Yield: 68.7%

**Top Precedent:**
- Core: Pd/Tri-tert-butylphosphinetetrafluoroborate
- Actual Yield: 96.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Yield Large Diff

---

### C-N - 2-chloroquinoxaline + aniline

**SMILES:** `Clc1cnc2ccccc2n1.Nc1ccccc1>>c1ccccc1Nc1cnc2ccccc2n1`

**ML Prediction:**
- Core: Pd/SPhos
- Predicted Yield: 69.1%

**Top Precedent:**
- Core: Pd/XPhos
- Actual Yield: 91.0%
- Precedent Count: 5

**Flags:**
- 🚩 Core Mismatch
- 🚩 Base Mismatch
- 🚩 Solvent Mismatch
- 🚩 Yield Large Diff

---

## ⚠ Minor Discrepancies

| Reaction | ML Yield | Precedent Yield | Flags |
|----------|----------|----------------|-------|
| B-H - Naphthylamine | 75.8% | 89.0% | core mismatch, base mismatch |
| Nc2ccccc2)cc1 (C-N - 4-MeO-Ph-Br + anili... | 76.5% | 90.0% | core mismatch, base mismatch |
| Nc2ccccc2)cc1 (C-N - 4-Me-Ph-Br + anilin... | 71.6% | 90.0% | core mismatch, base mismatch |
| F)(F)c1ccc(Nc2ccccc2)cc1 (C-N - 4-CF3-Ph... | 72.6% | 90.0% | core mismatch, base mismatch |
| =O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-NO2-Ph-B... | 80.1% | 90.0% | core mismatch, base mismatch |
| Nc2ccccc2)cc1 (C-N - 4-F-Ph-Br + aniline | 76.2% | 90.0% | core mismatch, base mismatch |
| C-N - 4-Cl-pyridine + aniline | 77.1% | 91.0% | base mismatch |
| C-N - 4-Br-pyridine + aniline | 77.6% | 90.0% | core mismatch, base mismatch |
| C-N - 3-Cl-pyridine + aniline | 70.3% | 91.0% | base mismatch, yield large diff |
| C-N - 2-Br-pyrimidine + aniline | 73.1% | 90.0% | core mismatch, base mismatch |
| C-N - 4-Cl-quinoline + aniline | 77.0% | 91.0% | base mismatch |
| C-N - 5-Br-indole + aniline | 75.6% | 90.0% | core mismatch, base mismatch |
| Nc2ccccc2)cn1 (C-N - 2-Br-4-methylpyrimi... | 82.3% | 90.0% | core mismatch, base mismatch |
| C)Nc1ccc(C#N)cc1 (C-N - 4-Cl-benzonitril... | 72.2% | 91.0% | base mismatch |
| F)(F)c1ccc(N2CCOCC2)cc1 (C-N - 4-Cl-benz... | 78.0% | 91.0% | base mismatch |
| Nc2ccc(Nc3ccccc3)cc2)cc1 (C-N - 1,4-dibr... | 74.1% | 90.0% | core mismatch, base mismatch |
| C-N - 2-methylbromobenzene + aniline | 77.6% | 96.0% | core mismatch, base mismatch |
| Nc2ccc3cc(Nc4ccccc4)ccc3n2)cc1 (C-N - 4,... | 90.0% | 91.0% | base mismatch |
| Nc2ccc3c(c2)OCO3)cc1 (C-N - 5-bromobenzo... | 81.8% | 90.0% | core mismatch, base mismatch |
| C-N - 2-chlorobenzothiazole + aniline | 74.7% | 91.0% | base mismatch |
| C-N - 2-bromobenzoxazole + aniline | 71.9% | 90.0% | core mismatch, base mismatch |

## Recommendations

1. **Review 36 major discrepancies** - These predictions may be unreliable
2. **2 predictions lack precedents** - Consider running these experimentally to validate

