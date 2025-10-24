# Reagent Registry Role Mapping Summary

## Overview

Successfully mapped 552 reagents from the z-Score Peaks dataset to the reagent registry system roles.

## Files Generated

- **Input:** `extracted_reagents_registry.md` (original extraction with dataset roles)
- **Output:** `reagents_mapped_to_registry_roles.md` (mapped to registry system roles)
- **Script:** `map_reagents_to_registry_roles.py` (mapping logic)

## Role Mapping Table

| Dataset Role | Registry System Role | Count | Description |
|--------------|---------------------|-------|-------------|
| Additive | `additive` | 97 | Phase-transfer agents, halide scavengers, fluoride sources, and related modifiers |
| Base | `base` | 63 | Bronsted/Lewis bases spanning amides, alkoxides, carbonates, superbases |
| **Coupling Reagent** | **`condensation_agent`** | **36** | **Carbodiimides, uronium/phosphonium activators, and similar condensers** |
| Ligand | `ligand` | 146 | Ligands including phosphines, NHCs, diimines, and ancillary donor sets |
| Solvent | `solvent` | 69 | Reaction media categorized by polarity, coordinating ability, and safety profile |
| Catalyst (simple salts) | `metal_precursor` | 114 | Metal salts or complexes that generate the catalytically active species |
| Catalyst (pre-coordinated) | `preformed_metal_catalyst` | 27 | Precatalysts supplied with ligands; typically used as-is |

**Total Reagents:** 552

## Key Mapping Details

### 1. Coupling Reagent → condensation_agent ✅

As requested, all coupling reagents (36 total) have been mapped to `condensation_agent` role.

**Top Coupling Reagents:**
- DIC (1,718 occurrences) - Diisopropylcarbodiimide
- T3P (862 occurrences) - Propylphosphonic anhydride
- BOP (604 occurrences) - Benzotriazol-1-yl-oxytris(dimethylamino)phosphonium hexafluorophosphate
- PyBrOP (524 occurrences)
- EDC (336 occurrences) - 1-Ethyl-3-(3-dimethylaminopropyl)carbodiimide

### 2. Catalyst Split into Two Categories

The original "Catalyst" role (141 reagents) was intelligently split based on structure:

**Metal Precursor (114 reagents):**
- Simple metal salts: PdCl2, CuI, Pd(OAc)2, etc.
- Top: CuI (3,786 occurrences), dtbpfPdCl2 (2,190), XantPhos Pd(allyl)Cl (2,048)

**Preformed Metal Catalyst (27 reagents):**
- Pre-coordinated complexes: Pd(PPh3)4, PEPPSI variants, XPhos Pd G3, etc.
- Top: tBuBrettPhos Pd(allyl)OTf (2,413), SPhos Pd(crotyl)Cl (1,532), XPhos Pd(crotyl)Cl (1,496)

### 3. All Other Roles Mapped 1:1

- Additive → additive
- Base → base  
- Ligand → ligand
- Solvent → solvent

## Usage Statistics

### Most Common Reagents by Role

| Role | Top Reagent | Occurrences |
|------|-------------|-------------|
| `additive` | water | 2,696 |
| `base` | K₃PO₄ | 10,296 |
| `condensation_agent` | DIC | 1,718 |
| `ligand` | N (generic ligand) | 20,740 |
| `metal_precursor` | CuI | 3,786 |
| `preformed_metal_catalyst` | tBuBrettPhos Pd(allyl)OTf | 2,413 |
| `solvent` | water | 18,343 |

## Output File Structure

The `reagents_mapped_to_registry_roles.md` file contains:

1. **Overview** - Dataset statistics and role counts
2. **Role Mapping Table** - Shows dataset → registry role conversions
3. **Detailed Sections** - One section per registry role with:
   - Role description
   - All reagents sorted by occurrence count (descending)
   - For each reagent:
     - Registry role (standardized)
     - Original role (from dataset)
     - Occurrence count
     - Reaction types where used
4. **Summary Statistics** - Top 3 reagents per role
5. **Notes** - Methodology and next steps

## Next Steps

1. **Cross-reference** with existing reagent registry files in `data/reagents/`
2. **Identify new entries** not already in the system
3. **Add chemical identifiers:**
   - CAS registry numbers
   - SMILES structures
   - InChI keys
4. **Standardize naming:**
   - Resolve abbreviations (e.g., "DIC" vs "Diisopropylcarbodiimide")
   - Handle synonyms
5. **Assign family classifications** within each role (e.g., carbodiimides vs uronium reagents for condensation_agent)
6. **Integrate** with reagent registry UI for curation

## Validation

✅ All 552 reagents successfully mapped  
✅ Coupling Reagent → condensation_agent mapping confirmed  
✅ Catalyst intelligently split into metal_precursor vs preformed_metal_catalyst  
✅ Output file generated (5,080 lines)  
✅ Summary statistics included  
✅ Ready for registry integration
