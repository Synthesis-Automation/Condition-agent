# Suzuki Database Validation Report

**Generated:** validate_all_suzuki_rules.py

## Summary Statistics

- **Total Rules:** 23
- **Passed:** 10 (43.5%)
- **Failed:** 13
- **Errors:** 0

## Detailed Results

### ✓ PASS (10)

#### SCDB-SUZ-ARBRI-GENERAL-SPhos

- **Description:** 4-Iodotoluene + PhB(OH)2 (ArI, SPhos set)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `SCDB-SUZ-ARBRI-GENERAL-SPhos`
- **Priority:** 46
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `SPhos`
  - Temperature: 45°C
  - Solvent: `THF/H2O (4:1)`

#### SCDB-SUZ-BULKY-NUC-XPHOS

- **Description:** PhBr + 2,6-dimethylphenylboronic acid (hindered nucleophile)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `SCDB-SUZ-BULKY-NUC-XPHOS`
- **Priority:** 48
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `XPhos`
  - Temperature: 80°C
  - Solvent: `toluene`

#### SCDB-SUZ-ARI-META-TBUXPHOS

- **Description:** 3,5-Dimethyl iodobenzene + PhB(OH)2 (meta-hindered ArI)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `SCDB-SUZ-ARI-META-TBUXPHOS`
- **Priority:** 55
- **Conditions:**
  - Pd source: `Pd(OAc)2`
  - Ligand: `tBuXPhos`
  - Temperature: 70°C
  - Solvent: `THF`

#### SCDB-SUZ-ARCL-EPoor-XPhos

- **Description:** 4-Chlorobenzonitrile + PhB(OH)2 (electron-poor ArCl)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `SCDB-SUZ-ARCL-EPoor-XPhos`
- **Priority:** 60
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `XPhos`
  - Temperature: 90°C
  - Solvent: `toluene/H2O (5:1)`

#### SCDB-SUZ-OTf-DPPF

- **Description:** Phenyl triflate + PhB(OH)2 (aryl triflate)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `SCDB-SUZ-OTf-DPPF`
- **Priority:** 70
- **Conditions:**
  - Pd source: `PdCl2(dppf)`
  - Ligand: `dppf`
  - Temperature: 60°C
  - Solvent: `THF`

#### SCDB-SUZ-ALKYL-PRIMARYI-9BBN

- **Description:** Propyl iodide + 9-BBN (primary alkyl sp3)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `SCDB-SUZ-ALKYL-PRIMARYI-9BBN`
- **Priority:** 62
- **Conditions:**
  - Pd source: `Pd(PPh3)4`
  - Ligand: `PPh3 (intrinsic)`
  - Temperature: 55°C
  - Solvent: `dioxane`

#### SCDB-SUZ-HET-THIOPHENE-FURAN-TRIBORONATE

- **Description:** PhBr + 2-thiophenylboronic acid (protodeboronation-prone)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `SCDB-SUZ-HET-THIOPHENE-FURAN-TRIBORONATE`
- **Priority:** 66
- **Conditions:**
  - Pd source: `Pd(PPh3)4`
  - Ligand: `PPh3 (intrinsic)`
  - Temperature: 45°C
  - Solvent: `THF/H2O (4:1)`

#### M-SUZ-VINYL-RT

- **Description:** Vinyl iodide + styrylboronic acid (scheme vinyl)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `M-SUZ-VINYL-RT`
- **Priority:** 75
- **Conditions:**
  - Pd source: `Pd(PPh3)4`
  - Ligand: `PPh3 (intrinsic)`
  - Temperature: 20°C
  - Solvent: `THF/H2O (9:1)`

#### SCDB-SUZ-HET-AZINE-BORON

- **Description:** PhBr + pyrimidine-5-boronic acid (azine boron partner)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `SCDB-SUZ-HET-AZINE-BORON`
- **Priority:** 64
- **Conditions:**
  - Pd source: `Pd(PPh3)4`
  - Ligand: `PPh3 (intrinsic)`
  - Temperature: 70°C
  - Solvent: `dioxane/H2O (4:1)`

#### SCDB-SUZ-DEFAULT-ArCl

- **Description:** Chlorobenzene + PhB(OH)2 (general ArCl default, no special activation)
- **Reason:** Correctly matched expected rule
- **Matched ID:** `SCDB-SUZ-DEFAULT-ArCl`
- **Priority:** 0
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `N/A`
  - Temperature: 100°C
  - Solvent: `toluene/H2O (5:1)`

### ✗ FAIL (13)

#### SCDB-SUZ-ARBRI-GENERAL-PPh3

- **Description:** Simple PhBr + PhB(OH)2 → biphenyl (ArBr, boronic acid)
- **Reason:** Matched different rule: SCDB-SUZ-ARBRI-GENERAL-SPhos (priority 46)
- **Matched ID:** `SCDB-SUZ-ARBRI-GENERAL-SPhos`
- **Priority:** 46
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `SPhos`
  - Temperature: 45°C
  - Solvent: `THF/H2O (4:1)`

#### SCDB-SUZ-2HETARYL-SPHOS

- **Description:** 4-Bromopyridine + PhB(OH)2 (2-hetaryl halide)
- **Reason:** Matched different rule: SCDB-SUZ-ARBRI-GENERAL-SPhos (priority 46)
- **Matched ID:** `SCDB-SUZ-ARBRI-GENERAL-SPhos`
- **Priority:** 46
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `SPhos`
  - Temperature: 45°C
  - Solvent: `THF/H2O (4:1)`

#### SCDB-SUZ-ARBR-ORTHO-XPhos

- **Description:** 2-Bromotoluene + PhB(OH)2 (ortho-hindered ArBr)
- **Reason:** Matched different rule: SCDB-SUZ-ARBRI-GENERAL-SPhos (priority 46)
- **Matched ID:** `SCDB-SUZ-ARBRI-GENERAL-SPhos`
- **Priority:** 46
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `SPhos`
  - Temperature: 45°C
  - Solvent: `THF/H2O (4:1)`

#### SCDB-SUZ-ARCL-ERich-L95

- **Description:** 4-Chloroanisole + PhB(OH)2 (electron-rich ArCl)
- **Reason:** Matched different rule: SCDB-SUZ-DEFAULT-ArCl (priority 0)
- **Matched ID:** `SCDB-SUZ-DEFAULT-ArCl`
- **Priority:** 0
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `N/A`
  - Temperature: 100°C
  - Solvent: `toluene/H2O (5:1)`

#### SCDB-SUZ-HET-2PYRIDYL-SLOWRELEASE

- **Description:** PhBr + 4-pyridyl-BF3K (slow-release 2-pyridyl partner)
- **Reason:** Matched different rule: SCDB-SUZ-DEFAULT-ArI_ArBr (priority 0)
- **Matched ID:** `SCDB-SUZ-DEFAULT-ArI_ArBr`
- **Priority:** 0
- **Conditions:**
  - Pd source: `Pd(PPh3)4`
  - Ligand: `N/A`
  - Temperature: 60°C
  - Solvent: `1,4-dioxane/H2O (4:1)`

#### SCDB-SUZ-VINYL-DPPF-RT

- **Description:** Vinyl bromide + vinylboronic acid (vinyl-vinyl coupling at RT)
- **Reason:** Matched different rule: M-SUZ-VINYL-RT (priority 75)
- **Matched ID:** `M-SUZ-VINYL-RT`
- **Priority:** 75
- **Conditions:**
  - Pd source: `Pd(PPh3)4`
  - Ligand: `PPh3 (intrinsic)`
  - Temperature: 20°C
  - Solvent: `THF/H2O (9:1)`

#### SCDB-SUZ-MIYAURA-BORYLATION

- **Description:** PhBr + B2pin2 → PhBpin (Miyaura borylation precursor - using correct B2pin2 structure)
- **Reason:** Matched different rule: SCDB-SUZ-DEFAULT-ArI_ArBr (priority 0)
- **Matched ID:** `SCDB-SUZ-DEFAULT-ArI_ArBr`
- **Priority:** 0
- **Conditions:**
  - Pd source: `Pd(PPh3)4`
  - Ligand: `N/A`
  - Temperature: 60°C
  - Solvent: `1,4-dioxane/H2O (4:1)`

#### M-SUZ-OTf-DPPF

- **Description:** 4-Methyl phenyl triflate + PhB(OH)2 (scheme-based OTf)
- **Reason:** Matched different rule: SCDB-SUZ-OTf-DPPF (priority 70)
- **Matched ID:** `SCDB-SUZ-OTf-DPPF`
- **Priority:** 70
- **Conditions:**
  - Pd source: `PdCl2(dppf)`
  - Ligand: `dppf`
  - Temperature: 60°C
  - Solvent: `THF`

#### M-SUZ-BF3K-GENERAL

- **Description:** 4-Bromoanisole + PhBF3K (simplified - use boronic acid as proxy)
- **Reason:** Matched different rule: SCDB-SUZ-ARBRI-GENERAL-SPhos (priority 46)
- **Matched ID:** `SCDB-SUZ-ARBRI-GENERAL-SPhos`
- **Priority:** 46
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `SPhos`
  - Temperature: 45°C
  - Solvent: `THF/H2O (4:1)`

#### SCDB-SUZ-DEFAULT-ArI_ArBr

- **Description:** 4-Fluorobromobenzene + PhB(OH)2 (general ArBr default)
- **Reason:** Matched different rule: SCDB-SUZ-ARBRI-GENERAL-SPhos (priority 46)
- **Matched ID:** `SCDB-SUZ-ARBRI-GENERAL-SPhos`
- **Priority:** 46
- **Conditions:**
  - Pd source: `Pd2(dba)3`
  - Ligand: `SPhos`
  - Temperature: 45°C
  - Solvent: `THF/H2O (4:1)`

#### SCDB-SUZ-DEFAULT-2HET_prone

- **Description:** PhBr + 2-furylboronic acid (protodeboronation-prone default)
- **Reason:** Matched different rule: SCDB-SUZ-HET-THIOPHENE-FURAN-TRIBORONATE (priority 66)
- **Matched ID:** `SCDB-SUZ-HET-THIOPHENE-FURAN-TRIBORONATE`
- **Priority:** 66
- **Conditions:**
  - Pd source: `Pd(PPh3)4`
  - Ligand: `PPh3 (intrinsic)`
  - Temperature: 45°C
  - Solvent: `THF/H2O (4:1)`

#### SCDB-SUZ-DEFAULT-vinyl

- **Description:** β-Bromostyrene + vinylboronic acid (vinyl default)
- **Reason:** Matched different rule: M-SUZ-VINYL-RT (priority 75)
- **Matched ID:** `M-SUZ-VINYL-RT`
- **Priority:** 75
- **Conditions:**
  - Pd source: `Pd(PPh3)4`
  - Ligand: `PPh3 (intrinsic)`
  - Temperature: 20°C
  - Solvent: `THF/H2O (9:1)`

#### SCDB-SUZ-DEFAULT-OTf

- **Description:** 4-CF3 phenyl triflate + PhB(OH)2 (general OTf default)
- **Reason:** Matched different rule: SCDB-SUZ-OTf-DPPF (priority 70)
- **Matched ID:** `SCDB-SUZ-OTf-DPPF`
- **Priority:** 70
- **Conditions:**
  - Pd source: `PdCl2(dppf)`
  - Ligand: `dppf`
  - Temperature: 60°C
  - Solvent: `THF`


## Recommendations

### Failed Rules Analysis

- **SCDB-SUZ-ARBRI-GENERAL-PPh3**: Matched `SCDB-SUZ-ARBRI-GENERAL-SPhos` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 46 beat expected rule

- **SCDB-SUZ-2HETARYL-SPHOS**: Matched `SCDB-SUZ-ARBRI-GENERAL-SPhos` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 46 beat expected rule

- **SCDB-SUZ-ARBR-ORTHO-XPhos**: Matched `SCDB-SUZ-ARBRI-GENERAL-SPhos` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 46 beat expected rule

- **SCDB-SUZ-ARCL-ERich-L95**: Matched `SCDB-SUZ-DEFAULT-ArCl` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 0 beat expected rule

- **SCDB-SUZ-HET-2PYRIDYL-SLOWRELEASE**: Matched `SCDB-SUZ-DEFAULT-ArI_ArBr` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 0 beat expected rule

- **SCDB-SUZ-VINYL-DPPF-RT**: Matched `M-SUZ-VINYL-RT` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 75 beat expected rule

- **SCDB-SUZ-MIYAURA-BORYLATION**: Matched `SCDB-SUZ-DEFAULT-ArI_ArBr` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 0 beat expected rule

- **M-SUZ-OTf-DPPF**: Matched `SCDB-SUZ-OTf-DPPF` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 70 beat expected rule

- **M-SUZ-BF3K-GENERAL**: Matched `SCDB-SUZ-ARBRI-GENERAL-SPhos` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 46 beat expected rule

- **SCDB-SUZ-DEFAULT-ArI_ArBr**: Matched `SCDB-SUZ-ARBRI-GENERAL-SPhos` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 46 beat expected rule

- **SCDB-SUZ-DEFAULT-2HET_prone**: Matched `SCDB-SUZ-HET-THIOPHENE-FURAN-TRIBORONATE` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 66 beat expected rule

- **SCDB-SUZ-DEFAULT-vinyl**: Matched `M-SUZ-VINYL-RT` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 75 beat expected rule

- **SCDB-SUZ-DEFAULT-OTf**: Matched `SCDB-SUZ-OTf-DPPF` instead
  - **Action:** Check SMARTS patterns, feature requirements, and priority
  - **Conflict:** Priority 70 beat expected rule

