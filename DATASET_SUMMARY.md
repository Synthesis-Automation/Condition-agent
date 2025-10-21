# Reaction Dataset Summary

**Last Updated:** October 21, 2025

## Overview

The Condition-agent system contains **140,427 total reactions** across 5 reaction families, with comprehensive condition and yield data.

## Dataset Breakdown

### 1. **Suzuki Coupling** (50,214 reactions)

- **File:** `data/reaction_dataset/Suzuki.jsonl` (81.5 MB)
- **DRFP Fingerprints:** `Suzuki_drfp.npz` (4.04 MB)
- **Yield Coverage:** 99.9% (50,158 reactions with yield data)
- **Yield Stats:** Mean 80.3%, Median 82.0%, Range 1-100%
- **Unique Components:**
  - 184 condition cores
  - 1,316 unique catalysts
  - 41 unique bases
  - 137 unique solvents

**Top Catalytic Systems:**

- Tetrakis(triphenylphosphine)palladium(0): 16,183× (80.1% avg yield)
- Bis(diphenylphosphino)ferrocene palladium(II) dichloride: 5,487× (76.6% avg yield)
- Bis(triphenylphosphine)palladium(II) dichloride: 3,658× (79.3% avg yield)

**Top Bases:** K₂CO₃ (18,560×), Na₂CO₃ (10,195×), K₃PO₄ (7,165×)

**Top Solvents:** Water (35,432×), 1,4-Dioxane (14,870×), Toluene (10,493×)

---

### 2. **Amide Formation** (41,427 reactions)

- **File:** `data/reaction_dataset/Amide_formation.jsonl` (75.4 MB)
- **DRFP Fingerprints:** `Amide_formation_drfp.npz` (2.89 MB)
- **Yield Coverage:** 99.8% (41,357 reactions with yield data)
- **Yield Stats:** Mean 78.6%, Median 81.0%, Range 1-100%
- **Unique Components:**
  - 101 condition cores
  - 267 unique catalysts
  - 38 unique bases
  - 98 unique solvents

**Top Catalytic Systems:**

- 4-(Dimethylamino)pyridine (DMAP): 3,126× (80.4% avg yield)
- N,N-Dimethylformamide (DMF): 2,086× (80.9% avg yield)
- Palladium: 1,120× (70.4% avg yield)

**Top Bases:** DIPEA (16,108×), Et₃N (9,330×), N-Methylmorpholine (2,019×)

**Top Solvents:** Dichloromethane (18,176×), DMF (16,305×), Water (7,285×)

---

### 3. **C-N Coupling** (24,571 reactions)

- **File:** `data/reaction_dataset/C_N_Coupling.jsonl` (38.0 MB)
- **DRFP Fingerprints:** `C_N_Coupling_drfp.npz` (1.79 MB)
- **Yield Coverage:** 99.9% (24,535 reactions with yield data)
- **Yield Stats:** Mean 79.9%, Median 82.0%, Range 1-100%
- **Unique Components:**
  - 176 condition cores
  - 676 unique catalysts
  - 44 unique bases
  - 130 unique solvents

**Top Catalytic Systems:**

- Copper(I) iodide: 1,393× (73.7% avg yield)
- Pd₂(dba)₃ + Xantphos: 455× (78.9% avg yield)
- Copper(II) acetate: 325× (72.3% avg yield)

**Top Bases:** K₂CO₃ (3,720×), Cs₂CO₃ (2,877×), DIPEA (2,464×)

**Top Solvents:** DMF (3,758×), Toluene (3,499×), Water (3,482×)

---

### 4. **C-O Coupling** (12,795 reactions)

- **File:** `data/reaction_dataset/C_O_Coupling.jsonl` (19.0 MB)
- **DRFP Fingerprints:** `C_O_Coupling_drfp.npz` (0.84 MB)
- **Yield Coverage:** 99.7% (12,754 reactions with yield data)
- **Yield Stats:** Mean 83.2%, Median 84.0%, Range 1-100%
- **Unique Components:**
  - 112 condition cores
  - 480 unique catalysts
  - 34 unique bases
  - 119 unique solvents

**Top Catalytic Systems:**

- Copper(I) iodide: 198× (78.1% avg yield)
- CuI + 98-98-6: 143× (82.1% avg yield)
- 18-Crown-6: 93× (82.4% avg yield)

**Top Bases:** K₂CO₃ (4,401×), Cs₂CO₃ (1,914×), NaH (953×)

**Top Solvents:** DMF (4,088×), Water (2,202×), DMSO (1,804×)

---

### 5. **C-S Coupling** (11,420 reactions)

- **File:** `data/reaction_dataset/C_S_Coupling.jsonl` (16.5 MB)
- **DRFP Fingerprints:** `C_S_Coupling_drfp.npz` (0.77 MB)
- **Yield Coverage:** 99.7% (11,382 reactions with yield data)
- **Yield Stats:** Mean 85.0%, Median 87.0%, Range 1-100%
- **Unique Components:**
  - 103 condition cores
  - 402 unique catalysts
  - 32 unique bases
  - 79 unique solvents

**Top Catalytic Systems:**

- Copper(I) iodide: 381× (84.9% avg yield)
- Pd₂(dba)₃ + Xantphos: 359× (88.2% avg yield)
- Copper: 177× (87.1% avg yield)

**Top Bases:** K₂CO₃ (2,520×), Cs₂CO₃ (1,030×), Et₃N (697×)

**Top Solvents:** DMF (3,193×), Water (2,100×), DMSO (1,035×)

---

## Additional Database Resources

### Condition Databases (`data/conditionDB/`)

- **amide_formation_db.json** - Curated amide formation conditions
- **C_N_Coupling_Cu_db.json** - Copper-catalyzed C-N coupling conditions
- **C_N_Coupling_Ni_db.json** - Nickel-catalyzed C-N coupling conditions
- **C_N_Coupling_Pd_db.json** - Palladium-catalyzed C-N coupling conditions
- **suzuki_db.json** - Curated Suzuki coupling conditions

### Protocol Database (`data/protocol_db/`)

17 curated protocol files including:

- Alkyl Iodide Borylation
- Aryl Iodide Cyanation
- Various Suzuki protocols
- Sonogashira coupling
- Mitsunobu reactions
- Specialized protocols (RCM, hydroacylation, etc.)

### Reagent Taxonomy (`data/reagents/`)

14 reagent classification files:

- **acid.json** - Acid reagents
- **additive.json** - Additives
- **base.json** - Base reagents
- **condensation_agent.json** - Condensation/coupling agents
- **enzyme.json** - Enzymatic catalysts
- **ligand.json** - Ligands
- **metal_precursor.json** - Metal precursors
- **organo_catalyst.json** - Organocatalysts
- **oxidant.json** - Oxidizing agents
- **preformed_metal_catalyst.json** - Pre-formed catalysts
- **reductant.json** - Reducing agents
- **solvent.json** - Solvents
- **other_reagent.json** - Other reagents
- **not_determined_reagents.json** - Unclassified reagents

---

## Dataset Features

Each reaction entry includes:

- **reaction_id**: Unique identifier
- **reaction_type**: Reaction family classification
- **condition_core**: Core catalytic/reagent condition
- **catalytic_system**: List of catalysts and ligands with CAS numbers
- **reagents**: List of reagents with roles (BASE, ACID, OXIDANT, etc.) and CAS numbers
- **solvents**: List of solvents with CAS numbers
- **conditions**: Temperature (°C), time (h), yield (%)
- **smiles**: Reactants and products as SMILES strings
- **reference**: Literature citation
- **precomputed**: Additional computed features including:
  - Reaction SMILES
  - Normalized reactants
  - Detected reaction family
  - Substrate features (for C-N coupling)

---

## Analytics Tools

Use the built-in analytics module for data exploration:

```python
from chemtools.dataset_analytics import (
    get_dataset_stats,
    get_common_catalysts,
    get_common_ligands,
    get_common_bases,
    get_common_solvents,
    get_common_catalytic_systems,
    get_condition_cores,
    print_analytics_summary
)

# Get statistics for a family
stats = get_dataset_stats("Suzuki")

# Get top 10 catalysts with yields
catalysts = get_common_catalysts("Suzuki", top_n=10)

# Print comprehensive summary
print_analytics_summary("C_N_Coupling", top_n=10)
```

### Quick Summary Generation

```bash
# Generate full dataset summary
python generate_dataset_summary.py
```

---

## Key Statistics Summary

| Reaction Family | Total Reactions | Avg Yield | Top Catalyst | Top Base | Top Solvent |
|----------------|----------------|-----------|--------------|----------|-------------|
| **Suzuki** | 50,214 | 80.3% | Pd(PPh₃)₄ (16,183×) | K₂CO₃ (18,560×) | Water (35,432×) |
| **Amide Formation** | 41,427 | 78.6% | DMAP (3,530×) | DIPEA (16,108×) | DCM (18,176×) |
| **C-N Coupling** | 24,571 | 79.9% | CuI (3,849×) | K₂CO₃ (3,720×) | DMF (3,758×) |
| **C-O Coupling** | 12,795 | 83.2% | CuI (1,022×) | K₂CO₃ (4,401×) | DMF (4,088×) |
| **C-S Coupling** | 11,420 | 85.0% | CuI (1,040×) | K₂CO₃ (2,520×) | DMF (3,193×) |

---

## Data Quality Notes

- **High yield coverage:** 99.7-99.9% of reactions have yield data
- **Limited temperature/time data:** Most reactions lack explicit temperature (0-0.1% coverage) and time data (0.1-0.6% coverage)
- **Comprehensive chemical coverage:** Thousands of unique catalysts, ligands, bases, and solvents
- **Literature-backed:** All reactions include literature references
- **DRFP fingerprints:** Pre-computed for similarity-based searching

---

## Files Location

- **Main datasets:** `data/reaction_dataset/*.jsonl`
- **Fingerprints:** `data/reaction_dataset/*_drfp.npz`
- **Condition DBs:** `data/conditionDB/*.json`
- **Protocols:** `data/protocol_db/*.json`
- **Reagent taxonomy:** `data/reagents/*.json`

---

*For detailed analytics and queries, see `chemtools/dataset_analytics.py` or run `python generate_dataset_summary.py`*
