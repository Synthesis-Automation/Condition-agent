# Suzuki Scope Policy Update (Taxonomy-Driven)

- Generated: 2026-02-15 17:35:36
- Policy source: `chemtools/taxonomy/data/reaction_types.v4.0.json` (`Suzuki_miyaura` constraints/products)
- Scope logic: require Suzuki-consistent evidence from taxonomy (C-C coupling and/or allowed Suzuki product motifs)

## Output Size Comparison (same 800-row sample)

- Previous (cleanup+routing only): 789 rows
- Updated (cleanup+routing+taxonomy scope): 770 rows
- Additional scope exclusions: 19 rows

## Exclusion Summary

- Total scope exclusions: 19
- Reason counts:
  - `formed_motif_conflict`: 19

## Example Excluded Rows

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
| 509 | `Miyaura_borylation` | `Ar-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.O=C1NCc2ccc(Br)cc21>>CC1(C)OB(c2ccc3c(c2)CNC3=O)OC1(C)C` |
| 576 | `Miyaura_borylation` | `Ar-Bpin` | `B-C` | `CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C.CN1CCc2ccc(Br)cc2C1=O>>CN1CCc2ccc(B3OC(C)(C)C(C)(C)O3)cc2C1=O` |
