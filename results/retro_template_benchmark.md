# Retro Template Benchmark Results

**Date**: 2026-03-30
**Script**: `scripts/benchmark_retro_templates.py`

## Configuration

| Parameter | Value |
|---|---|
| Templates used | 500 (top by quality score, from 10,156 extracted) |
| Template source | `results/extracted_templates.json` |
| Holdout reactions | 4,416 (20% deterministic split of 20,831 total) |
| Reaction types evaluated | 102 |
| Quality threshold | 0.0 (no filtering) |
| Data | All 119 HTE CSVs, 200 rows/file max |

## Overall Metrics

| Metric | Value | Description |
|---|---|---|
| **Coverage** | **80.9%** | Products with ≥1 valid template disconnection |
| **Exact match (SMILES)** | **7.7%** | Predicted reactants exactly match actual (canonical SMILES) |
| **InChIKey match** | **20.5%** | Any predicted reactant matches actual by connectivity (stereo-insensitive) |
| **Partial match (Tan > 0.7)** | **12.1%** | Best-match Tanimoto similarity > 0.7 |
| Mean similarity | 0.3456 | Mean best-match Tanimoto across all holdout reactions |
| Median similarity | 0.2941 | Median best-match Tanimoto |
| Mean component hit fraction | 0.1658 | Fraction of actual reactant components recovered |
| Templates per hit (mean) | 10.5 | Average number of templates that fire per product |
| Elapsed time | 588.0s | ~7.5 rxn/s with substructure pre-screening |

## Critical Fix: Leaving-Group Atom Maps

During benchmarking, a critical bug was discovered: **rdchiral returned 0 outcomes for every template** because extracted templates contained atom maps on leaving-group atoms (e.g., `[Br;H0;D1;+0:1]`, `[B;H0;D3;+0:6]`) that have no corresponding atom in the product pattern.

**Fix**: Added `_clean_template_atom_maps()` which strips atom map numbers from reactant-side atoms that don't appear on the product side. This restored rdchiral functionality and moved results from 0% → 7.7% exact match.

## Per Reaction Type Results

### Top Performers (n ≥ 10, sorted by InChIKey match)

| Reaction Type | n | Coverage% | Exact% | IKey% |
|---|---|---|---|---|
| LGDisp+C-O | 10 | 100.0 | 80.0 | 100.0 |
| Chan-Lam C-N Coupling | 47 | 100.0 | 0.0 | 97.9 |
| C-O Coupling | 56 | 100.0 | 51.8 | 94.6 |
| Suzuki-Miyaura | 59 | 100.0 | 13.6 | 93.2 |
| Minisci acylation | 21 | 95.2 | 90.5 | 90.5 |
| Hiyama | 53 | 98.1 | 0.0 | 88.7 |
| Kumada | 17 | 88.2 | 17.6 | 82.4 |
| Sonogashira | 28 | 100.0 | 46.4 | 82.1 |
| Negishi | 16 | 100.0 | 12.5 | 75.0 |
| LGDisp+C-N | 31 | 96.8 | 22.6 | 67.7 |
| C-N Coupling | 47 | 93.6 | 38.3 | 66.0 |
| Stille | 54 | 96.3 | 0.0 | 57.4 |
| C-S Coupling | 47 | 80.9 | 23.4 | 53.2 |
| Amide formation | 98 | 100.0 | 18.4 | 53.1 |
| Reductive amination | 92 | 95.7 | 15.2 | 46.7 |
| Ir-catalyzed C-H borylation | 13 | 61.5 | 23.1 | 46.2 |
| Alkyl nucleophilic substitution | 131 | 90.1 | 13.0 | 44.3 |
| Alkylation Friedel-Crafts | 50 | 92.0 | 44.0 | 44.0 |
| C-N+C-C | 30 | 100.0 | 23.3 | 43.3 |
| Esterification | 81 | 77.8 | 12.3 | 27.2 |

### High Coverage but Zero Precision (n ≥ 20)

| Reaction Type | n | Coverage% | Exact% | IKey% |
|---|---|---|---|---|
| Reduction carbonyl→alcohol | 218 | 56.0 | 0.0 | 0.0 |
| Oxidation 1° alcohol→aldehyde | 164 | 78.0 | 0.0 | 0.0 |
| Alcohol→Alkyl halide | 133 | 62.4 | 0.0 | 0.0 |
| Oxidation 2° alcohol→ketone | 99 | 77.8 | 0.0 | 0.0 |
| Baeyer-Villiger oxidation | 65 | 90.8 | 0.0 | 0.0 |
| Benzylic C-H oxidation→carbonyl | 64 | 71.9 | 0.0 | 0.0 |
| Acyl halide formation | 57 | 68.4 | 0.0 | 0.0 |
| Click (CuAAC) | 55 | 94.5 | 0.0 | 0.0 |
| Trifluoromethylation | 26 | 76.9 | 0.0 | 0.0 |
| Halogenation (unsaturated) | 28 | 35.7 | 0.0 | 0.0 |
| Decarboxylative coupling | 20 | 80.0 | 0.0 | 0.0 |
| Minisci alkylation | 20 | 85.0 | 0.0 | 0.0 |

### All Reaction Types (sorted by count)

| Reaction Type | n | Cov | Cov% | Exact | Exact% | IKey | IKey% |
|---|---|---|---|---|---|---|---|
| Unknown | 421 | 376 | 89.3 | 74 | 17.6 | 101 | 24.0 |
| Event:C-N | 218 | 177 | 81.2 | 11 | 5.0 | 29 | 13.3 |
| Reduction_carbonyl_to_alcohol | 218 | 122 | 56.0 | 0 | 0.0 | 0 | 0.0 |
| Alkene_Hydrogenation | 202 | 151 | 74.8 | 0 | 0.0 | 6 | 3.0 |
| Multi-Event:LGDisp+C-C | 185 | 142 | 76.8 | 2 | 1.1 | 28 | 15.1 |
| Event:C-C | 180 | 145 | 80.6 | 17 | 9.4 | 36 | 20.0 |
| Wittig_reaction | 172 | 138 | 80.2 | 9 | 5.2 | 28 | 16.3 |
| Oxidation_primary_alcohol_to_aldehyde | 164 | 128 | 78.0 | 0 | 0.0 | 0 | 0.0 |
| Alcohol_to_Alkyl_Halide | 133 | 83 | 62.4 | 0 | 0.0 | 0 | 0.0 |
| Alkyl_Nucleophilic_Substitution | 131 | 118 | 90.1 | 17 | 13.0 | 58 | 44.3 |
| Alkene_oxidative_cleavage_to_carbonyl | 108 | 88 | 81.5 | 0 | 0.0 | 2 | 1.9 |
| Oxidation_secondary_alcohol_to_ketone | 99 | 77 | 77.8 | 0 | 0.0 | 0 | 0.0 |
| Amide_formation | 98 | 98 | 100.0 | 18 | 18.4 | 52 | 53.1 |
| Multi-Event:Ann+C-C | 97 | 63 | 64.9 | 1 | 1.0 | 1 | 1.0 |
| Reductive_amination | 92 | 88 | 95.7 | 14 | 15.2 | 43 | 46.7 |
| Esterification | 81 | 63 | 77.8 | 10 | 12.3 | 22 | 27.2 |
| Event:C-O | 74 | 51 | 68.9 | 5 | 6.8 | 9 | 12.2 |
| Multi-Event:Ann+C-N | 70 | 65 | 92.9 | 0 | 0.0 | 2 | 2.9 |
| Baeyer_villiger_oxidation | 65 | 59 | 90.8 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Ann+C-N+C-C | 64 | 58 | 90.6 | 1 | 1.6 | 3 | 4.7 |
| Benzylic_CH_oxidation_to_carbonyl | 64 | 46 | 71.9 | 0 | 0.0 | 0 | 0.0 |
| Event:Ann | 63 | 49 | 77.8 | 1 | 1.6 | 1 | 1.6 |
| Suzuki_miyaura | 59 | 59 | 100.0 | 8 | 13.6 | 55 | 93.2 |
| Acyl_Halides_formation | 57 | 39 | 68.4 | 0 | 0.0 | 0 | 0.0 |
| C_O_Coupling | 56 | 56 | 100.0 | 29 | 51.8 | 53 | 94.6 |
| Sandmeyer | 55 | 26 | 47.3 | 0 | 0.0 | 2 | 3.6 |
| Click_cuaac | 55 | 52 | 94.5 | 0 | 0.0 | 0 | 0.0 |
| Stille | 54 | 52 | 96.3 | 0 | 0.0 | 31 | 57.4 |
| Hiyama | 53 | 52 | 98.1 | 0 | 0.0 | 47 | 88.7 |
| Aromatic_Halogen_Exchange | 52 | 32 | 61.5 | 0 | 0.0 | 0 | 0.0 |
| Alkylation_friedel_crafts | 50 | 46 | 92.0 | 22 | 44.0 | 22 | 44.0 |
| C_N_Coupling | 47 | 44 | 93.6 | 18 | 38.3 | 31 | 66.0 |
| C_S_Coupling | 47 | 38 | 80.9 | 11 | 23.4 | 25 | 53.2 |
| Chan_Lam_C_N_Coupling | 47 | 47 | 100.0 | 0 | 0.0 | 46 | 97.9 |
| Multi-Event:Ann+C-O | 46 | 34 | 73.9 | 0 | 0.0 | 7 | 15.2 |
| Reduction_ester_to_alcohol | 37 | 30 | 81.1 | 1 | 2.7 | 3 | 8.1 |
| Multi-Event:Amid+C-N | 36 | 33 | 91.7 | 3 | 8.3 | 5 | 13.9 |
| Event:C-S | 35 | 30 | 85.7 | 0 | 0.0 | 2 | 5.7 |
| Multi-Event:LGDisp+C-N | 31 | 30 | 96.8 | 7 | 22.6 | 21 | 67.7 |
| Multi-Event:C-N+C-C | 30 | 30 | 100.0 | 7 | 23.3 | 13 | 43.3 |
| Azide_coupling | 30 | 16 | 53.3 | 0 | 0.0 | 0 | 0.0 |
| Halogenation_unsaturated | 28 | 10 | 35.7 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:C-N+C-O | 28 | 28 | 100.0 | 0 | 0.0 | 2 | 7.1 |
| Sonogashira | 28 | 28 | 100.0 | 13 | 46.4 | 23 | 82.1 |
| Trifluoromethylation | 26 | 20 | 76.9 | 0 | 0.0 | 0 | 0.0 |
| Halogenation_aromatic | 21 | 15 | 71.4 | 0 | 0.0 | 1 | 4.8 |
| Minisci_acylation | 21 | 20 | 95.2 | 19 | 90.5 | 19 | 90.5 |
| Decarboxylative_Coupling | 20 | 16 | 80.0 | 0 | 0.0 | 0 | 0.0 |
| Minisci_alkylation | 20 | 17 | 85.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Ann+C-O+C-C | 19 | 14 | 73.7 | 1 | 5.3 | 2 | 10.5 |
| Alkene_oxidative_cleavage_to_carboxylic_acid | 17 | 17 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Ann+C-S+C-C | 17 | 12 | 70.6 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Ann+C-S | 17 | 16 | 94.1 | 0 | 0.0 | 0 | 0.0 |
| Kumada | 17 | 15 | 88.2 | 3 | 17.6 | 14 | 82.4 |
| Aliphatic_Halide_Exchange | 16 | 7 | 43.8 | 0 | 0.0 | 0 | 0.0 |
| Acylation_amide | 16 | 16 | 100.0 | 0 | 0.0 | 3 | 18.8 |
| Negishi | 16 | 16 | 100.0 | 2 | 12.5 | 12 | 75.0 |
| Benzylic_CH_oxidation_to_carboxylic_acid | 14 | 4 | 28.6 | 0 | 0.0 | 0 | 0.0 |
| Michael_addition | 13 | 11 | 84.6 | 0 | 0.0 | 0 | 0.0 |
| Ircatalyzed_CH_borylation | 13 | 8 | 61.5 | 3 | 23.1 | 6 | 46.2 |
| Unresolved:NoMotifEvidence | 13 | 10 | 76.9 | 0 | 0.0 | 3 | 23.1 |
| Reduction_nitrile_to_amine | 12 | 12 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Event:LGDisp | 10 | 6 | 60.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:LGDisp+C-O | 10 | 10 | 100.0 | 8 | 80.0 | 10 | 100.0 |
| Arylation_acidic_C_H | 9 | 8 | 88.9 | 0 | 0.0 | 1 | 11.1 |
| Multi-Event:C-O+C-C | 8 | 6 | 75.0 | 0 | 0.0 | 2 | 25.0 |
| Quinazolinone_annulation | 8 | 8 | 100.0 | 0 | 0.0 | 1 | 12.5 |
| Multi-Event:Ann+C-N+C-O+C-C | 8 | 7 | 87.5 | 1 | 12.5 | 1 | 12.5 |
| Multi-Event:Ann+C-N+C-O | 8 | 8 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:LGDisp+C-S | 6 | 3 | 50.0 | 0 | 0.0 | 1 | 16.7 |
| Multi-Event:LGDisp+C-N+C-O | 6 | 6 | 100.0 | 0 | 0.0 | 6 | 100.0 |
| Reduction_imine_to_amine | 5 | 4 | 80.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Amid+Ann+C-N | 5 | 5 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:EsterHyd+C-N | 4 | 4 | 100.0 | 0 | 0.0 | 2 | 50.0 |
| Multi-Event:Ann+LGDisp+C-N+C-S | 4 | 4 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Rosenmund_reduction | 4 | 4 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Ann+LGDisp+C-S+C-C | 4 | 3 | 75.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:LGDisp+C-N+C-C | 3 | 3 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Ann+LGDisp+C-N+C-O | 3 | 2 | 66.7 | 0 | 0.0 | 1 | 33.3 |
| Multi-Event:Amid+Ann+C-N+C-C | 3 | 3 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:C-N+C-S | 3 | 3 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Amid+C-N+C-C | 3 | 3 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Arylation_Ar_H | 3 | 3 | 100.0 | 0 | 0.0 | 3 | 100.0 |
| Multi-Event:Ann+LGDisp+C-N+C-C | 3 | 3 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:C-S+C-C | 3 | 2 | 66.7 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Ann+C-N+C-S | 2 | 2 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Ann+LGDisp+C-C | 2 | 2 | 100.0 | 0 | 0.0 | 1 | 50.0 |
| Multi-Event:Ann+LGDisp+C-O+C-C | 2 | 2 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Ann+LGDisp+C-N | 2 | 1 | 50.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:BzOAlk+LGDisp+C-O | 2 | 2 | 100.0 | 1 | 50.0 | 2 | 100.0 |
| Multi-Event:Ann+C-N+C-S+C-C | 2 | 2 | 100.0 | 1 | 50.0 | 1 | 50.0 |
| Cyanation_coupling | 2 | 2 | 100.0 | 0 | 0.0 | 1 | 50.0 |
| Multi-Event:BzOAlk+EsterHyd+LGDisp | 2 | 2 | 100.0 | 0 | 0.0 | 2 | 100.0 |
| Multi-Event:C-N+C-O+C-C | 1 | 1 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Miyaura_borylation | 1 | 0 | 0.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Amid+C-N+C-S | 1 | 1 | 100.0 | 0 | 0.0 | 1 | 100.0 |
| Event:EsterHyd | 1 | 1 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Amid+Ann+C-N+C-O | 1 | 1 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Amid+Ann+C-N+C-S | 1 | 1 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:Amid+LGDisp+C-N | 1 | 1 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Multi-Event:C-O+C-S+C-C | 1 | 1 | 100.0 | 0 | 0.0 | 0 | 0.0 |
| Reduction_nitro | 1 | 1 | 100.0 | 0 | 0.0 | 0 | 0.0 |

## Summary Statistics

- **Total reaction types**: 102
- **Types with any exact match**: 32 (31.4%)
- **Types with any InChIKey match**: 57 (55.9%)
- **Types with zero coverage**: 1 (Miyaura borylation, n=1)

## Key Observations

1. **Cross-coupling reactions are well-served**: Suzuki (93% IKey), C-O coupling (95%), Hiyama (89%), Sonogashira (82%), Kumada (82%), Chan-Lam (98%), Negishi (75%), Stille (57%) all show strong connectivity-level matching.

2. **Bond-forming reactions work well**: Amide formation (53%), reductive amination (47%), alkyl nucleophilic substitution (44%), Friedel-Crafts alkylation (44%), esterification (27%) show meaningful precision.

3. **Redox reactions fail**: Carbonyl reductions, alcohol oxidations, Baeyer-Villiger, and C-H oxidations show high coverage (templates fire) but 0% precision (wrong precursors). The templates match the product substructure but disconnect at wrong bonds.

4. **InChIKey >> Exact match**: The gap between IKey% and Exact% (especially for Suzuki: 93% vs 14%) indicates templates correctly identify the disconnection site but predict different leaving groups (e.g., ArBr instead of ArI, ArB(OH)₂ instead of ArB(pin)).

5. **Coverage is high across the board**: 80.9% overall, with 101/102 types having >0% coverage. The 500 top-quality templates provide broad applicability.

## Metric Definitions

- **Coverage**: Product has ≥1 valid template disconnection (rdchiral returns valid SMILES)
- **Exact match**: Canonical SMILES of predicted reactants exactly match actual (allowing subset/superset)
- **InChIKey match**: At least one predicted reactant matches an actual reactant by InChIKey first block (connectivity layer, stereo-insensitive)
- **Partial match**: Best-match average Tanimoto similarity between predicted and actual reactant fingerprints > 0.7
- **Component hit fraction**: Fraction of actual reactant components found in predicted set (by InChIKey connectivity)
