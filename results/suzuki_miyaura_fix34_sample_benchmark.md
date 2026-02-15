# Suzuki_miyaura Fix 3/4 Sample Investigation

- Generated: 2026-02-15 17:12:50
- Source: `data\HTE_db\literature\suzuki_miyaura_canonical.csv`
- Sample: 800 rows (seed=20260215)
- Sample input: `results\suzuki_miyaura_fix34_sample_input.csv`
- Sample output: `results\suzuki_miyaura_fix34_sample_output.csv`
- Audit CSV: `results\suzuki_miyaura_fix34_sample_audit.csv`

## Problematic Metrics (Before vs After)

| Metric | Before | After |
|---|---:|---:|
| `total_rows` | 800 | 789 (-11 vs before) |
| `unknown_detected_reaction_type` | 4 | 0 (-4 vs before) |
| `unresolved_detected_reaction_type` | 0 | 0 (0 vs before) |
| `missing_reaction_key` | 0 | 0 (0 vs before) |
| `missing_reaction_events` | 0 | 0 (0 vs before) |
| `reaction_key_empty_or_degenerate` | 0 | 0 (0 vs before) |
| `reaction_smiles_no_arrow` | 0 | 0 (0 vs before) |
| `reaction_smiles_missing_side` | 0 | 0 (0 vs before) |

## Converter Audit

- `routing_excluded`: 11
- `cleanup_applied`: 16
- coordination components removed: 0
- counterion components removed: 23

## Top Detected Reaction Types (Sample)

### Before
| type | count |
|---|---:|
| `Suzuki_miyaura` | 692 |
| `Multi-Event:LGDisp+C-C` | 39 |
| `Miyaura_borylation` | 16 |
| `Event:Ann` | 12 |
| `Alkylation_friedel_crafts` | 9 |
| `Unknown` | 4 |
| `Multi-Event:LGDisp+C-O+C-C` | 3 |
| `Kumada` | 3 |
| `Halogenation_aromatic` | 2 |
| `Sonogashira` | 2 |

### After
| type | count |
|---|---:|
| `Suzuki_miyaura` | 692 |
| `Multi-Event:LGDisp+C-C` | 40 |
| `Miyaura_borylation` | 16 |
| `Alkylation_friedel_crafts` | 9 |
| `Kumada` | 3 |
| `Multi-Event:LGDisp+C-O+C-C` | 3 |
| `Multi-Event:LGDisp+C-N+C-O+C-C` | 2 |
| `Ircatalyzed_CH_borylation` | 2 |
| `Halogenation_aromatic` | 2 |
| `C_O_Coupling` | 2 |

## Potential Improvement Buckets (Suzuki-specific)

| bucket | count |
|---|---:|
| `borylation_like_non_suzuki_cc` | 19 |
| `non_cc_without_borylation_label` | 6 |

## Manual Review Examples (borylation-like in Suzuki dataset)

| sample_row | detected_reaction_type | formed_motifs | formed_bond_classes | reaction_smiles (preview) |
|---:|---|---|---|---|
| 37 | `Miyaura_borylation` | `HeteroAr-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.COc1cc(-n2cnn3cc(Br)cc3c2=O)ccc1OCC(C)(C)O>>COc1cc(-n2cnn3cc(B4OC(C)(C)C(C)(C)O4...` |
| 79 | `Event:LGDisp` | `Alkyl-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.CC(C)(C)OC(=O)N1CC2(CC(I)C2)C1>>CC(C)(C)OC(=O)N1CC2(CC(B3OC(C)(C)C(C)(C)O3)C2)C1` |
| 116 | `Miyaura_borylation` | `HeteroAr-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.CNc1cnc2ccc(Br)cc2n1>>CNc1cnc2ccc(B3OC(C)(C)C(C)(C)O3)cc2n1` |
| 137 | `Miyaura_borylation` | `HeteroAr-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.COc1cc2cc(Br)ccc2nc1-c1ccccc1>>COc1cc2cc(B3OC(C)(C)C(C)(C)O3)ccc2nc1-c1ccccc1` |
| 158 | `Miyaura_borylation` | `HeteroAr-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.Cc1ccc(S(=O)(=O)N2CC(c3ccccc3)n3c2cc2cc(Br)ccc23)cc1>>Cc1ccc(S(=O)(=O)N2CC(c3ccc...` |
| 161 | `Miyaura_borylation` | `Ar-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.C[C@@H]1C(=O)N(C)c2ccc(Br)cc2N1C1CCC1>>C[C@@H]1C(=O)N(C)c2ccc(B3OC(C)(C)C(C)(C)O...` |
| 247 | `Miyaura_borylation` | `Ar-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.FC(F)(F)c1ccc2c(c1)CN1CN2Cc2cc(Br)ccc21>>CC1(C)OB(c2ccc3c(c2)CN2CN3Cc3cc(C(F)(F)...` |
| 307 | `Miyaura_borylation` | `HeteroAr-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.Clc1cc(I)cc(N2C[C@@H]3C[C@H]2CO3)n1>>CC1(C)OB(c2cc(Cl)nc(N3C[C@@H]4C[C@H]3CO4)c2...` |
| 336 | `Miyaura_borylation` | `HeteroAr-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.C[C@@H]1CN(c2ccc(Br)cn2)CCN1C(=O)OC(C)(C)C>>C[C@@H]1CN(c2ccc(B3OC(C)(C)C(C)(C)O3...` |
| 403 | `Miyaura_borylation` | `Ar-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.O=C1COc2cc(Br)ccc2N1>>CC1(C)OB(c2ccc3c(c2)OCC(=O)N3)OC1(C)C` |

## Recommended Next Improvements

1. Add an optional Suzuki scope gate: keep `C-C` coupling rows, route `B-C`/`*-Bpin` outcomes to `Miyaura_borylation` or exclude for pure Suzuki training sets.
2. Add explicit bucket export for `borylation_like_non_suzuki_cc` rows during conversion so curation is traceable.
3. Add 2-3 tests for Suzuki scope behavior (true Suzuki C-C keep, borylation route out, mixed multi-event flag).
