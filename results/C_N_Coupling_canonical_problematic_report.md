# C_N_Coupling_canonical.csv Problematic/Unsolved Check

- Dataset: `data/HTE_db/literature/C_N_Coupling_canonical.csv`
- Generated: 2026-02-15 12:22:05
- Total rows: 25445
- Unique problematic rows (union of rules below): 220

## Rules Used

- `unresolved_detected_reaction_type`: `detected_reaction_type` starts with `Unresolved`
- `unknown_detected_reaction_type`: `detected_reaction_type` is `Unknown`
- `missing_reaction_key`: empty `Reaction_Key`
- `missing_reaction_events`: empty `Reaction_Events`
- `reaction_key_empty_or_degenerate`: key contains `|[] -> []` or `| ->`
- `reaction_smiles_no_arrow`: missing `>>` in `reaction_smiles`
- `reaction_smiles_missing_side`: `reaction_smiles` has empty reactant or product side

## Counts

| Bucket | Count |
|---|---:|
| `unresolved_detected_reaction_type` | 4 |
| `unknown_detected_reaction_type` | 215 |
| `missing_reaction_key` | 7 |
| `missing_reaction_events` | 7 |
| `reaction_key_empty_or_degenerate` | 2 |
| `reaction_smiles_no_arrow` | 0 |
| `reaction_smiles_missing_side` | 0 |

## Examples (up to 10 per bucket)

### unresolved_detected_reaction_type (4)

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | Reaction_Key / Reaction_Events (preview) |
|---:|---|---|---|---|
| 165 | `C_N_Coupling` | `Unresolved:NoMotifEvidence` | `NC=O.Cc1cc(Br)c2c(c1)C(C)(C)c1cc(C)cc(Br)c1O2>>Cc1cc(NC=O)c2c(c1)C(C)(C)c1cc(C)cc(NC=O)c1O2` | `|Ar-Br|Unclassified-Reactant -> [] | spectators: Ar-Alkyl|Ar-OR` |
| 5183 | `C_N_Coupling` | `Unresolved:NoMotifEvidence` | `COc1ccc2nc(Cl)c(N3CCN(C(=O)OC(C)(C)C)CC3)nc2c1.O=C(Nc1ncc([N+](=O)[O-])s1)N1CCNCC1>>COc1ccc2nc(N3CC…` | `|RCH2-NHR -> [] | spectators: Ar-NO2|Ar-NR2|Ar-Urea|HeteroAr-Cl|HeteroAr-H|HeteroAr-OR|Pyrimidine|Thiazole` |
| 16541 | `C_N_Coupling` | `Unresolved:NoMotifEvidence` | `COc1ccc(Nc2ccc(OC)cc2)cc1.O=c1c(O)c(O)c1=O>>COC1=CC=C(N(C2=CC=C(OC)C=C2)C2=C(O)C(=[NH](C3=CC=C(OC)C…` | `|[] -> []` |
| 24077 | `C_N_Coupling` | `Unresolved:NoMotifEvidence` | `CN1CCOCC1.COc1nc(Cl)nc(OC)n1>>COC1=NC(OC)=NC([N+]2(C)CCOCC2)=N1.[F-][B+3]([F-])([F-])[F-]` | `|HeteroAr-Cl|RCH2-NR2 -> [] | spectators: HeteroAr-OR|RCH2-OR` |

### unknown_detected_reaction_type (215)

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | Reaction_Key / Reaction_Events (preview) |
|---:|---|---|---|---|
| 344 | `C_N_Coupling` | `Unknown` | `Fc1cc(N2CCOCC2)cc(F)n1>>NN=c1cc(N2CCOCC2)cc(F)[nH]1` | `|HeteroAr-F -> AromN-H | spectators: Ar-NR2|HeteroAr-H|Pyridine|RCH2-OR` |
| 421 | `C_N_Coupling` | `Unknown` | `COc1cc2c(cc1OC)-c1c(cnc(Cl)c1C#N)C2>>COc1cc2c(cc1OC)-c1c(c[nH]c(=NN)c1C#N)C2` | `|HeteroAr-Cl -> AromN-H | spectators: Ar-Alkyl|Ar-Ar|Ar-OR|HeteroAr-CN|HeteroAr-H|Pyridine` |
| 480 | `C_N_Coupling` | `Unknown` | `Clc1ncnc2sc3c(c12)CCCC3>>NN=c1nc[nH]c2sc3c(c12)CCCC3` | `|HeteroAr-Cl -> AromN-H | spectators: Ar-Alkyl|HeteroAr-H|Pyrimidine|Thiophene` |
| 481 | `C_N_Coupling` | `Unknown` | `Clc1ncnc2sc3c(c12)CCC3>>NN=c1nc[nH]c2sc3c(c12)CCC3` | `|HeteroAr-Cl -> AromN-H | spectators: Ar-Alkyl|HeteroAr-H|Pyrimidine|Thiophene` |
| 806 | `C_N_Coupling` | `Unknown` | `Cn1cc(-c2cnc(F)c(F)c2)cn1>>Cn1cc(-c2c[nH]c(=NN)c(F)c2)cn1` | `|HeteroAr-F -> AromN-H | spectators: Alkyl-AromN|Ar-Ar|HeteroAr-H|Pyrazole|Pyridine` |
| 1066 | `C_N_Coupling` | `Unknown` | `CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C.CC(C)(C)[C@H](NC(=O)c1ccccc1)[C@@H](O)c1cnccc1Cl>>CN(C)c1ccncc1[…` | `` |
| 1067 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1.Cl` | `` |
| 1070 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.CC(C)(C)OC(=O)NC1CCN(->[Zn+2]([Cl-])([Cl-])<-N2CCC(NC(=O)OC(C)(C)C)C2)C1>>CC(C)(C)OC(…` | `` |
| 1073 | `C_N_Coupling` | `Unknown` | `CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C.OC(c1cnccc1Cl)C(c1ccccc1)(c1ccccc1)c1ccccc1>>CN(C)c1ccncc1C(O)C(…` | `` |
| 1074 | `C_N_Coupling` | `Unknown` | `Clc1ccncc1-c1c(-c2ccccc2)ccc2ccccc12.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1-c1c(-c2ccccc…` | `` |

### missing_reaction_key (7)

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | Reaction_Key / Reaction_Events (preview) |
|---:|---|---|---|---|
| 1066 | `C_N_Coupling` | `Unknown` | `CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C.CC(C)(C)[C@H](NC(=O)c1ccccc1)[C@@H](O)c1cnccc1Cl>>CN(C)c1ccncc1[…` | `` |
| 1067 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1.Cl` | `` |
| 1070 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.CC(C)(C)OC(=O)NC1CCN(->[Zn+2]([Cl-])([Cl-])<-N2CCC(NC(=O)OC(C)(C)C)C2)C1>>CC(C)(C)OC(…` | `` |
| 1073 | `C_N_Coupling` | `Unknown` | `CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C.OC(c1cnccc1Cl)C(c1ccccc1)(c1ccccc1)c1ccccc1>>CN(C)c1ccncc1C(O)C(…` | `` |
| 1074 | `C_N_Coupling` | `Unknown` | `Clc1ccncc1-c1c(-c2ccccc2)ccc2ccccc12.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1-c1c(-c2ccccc…` | `` |
| 1075 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.C[C@H](N->[Zn+2]([Cl-])([Cl-])<-N[C@@H](C)C1=CC=CC=C1)C1=CC=CC=C1>>C[C@H](Nc1ccncc1)c…` | `` |
| 23043 | `C_N_Coupling` | `Unknown` | `O=Cc1c(Cl)c2ccccc2oc1=O.Nc1ccncc1.CCO[PH](=O)OCC>>CCOP(OCC)(=O->[Zn+2]([Cl-])([Cl-])[Cl-])C(O)C1=C(…` | `` |

### missing_reaction_events (7)

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | Reaction_Key / Reaction_Events (preview) |
|---:|---|---|---|---|
| 1066 | `C_N_Coupling` | `Unknown` | `CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C.CC(C)(C)[C@H](NC(=O)c1ccccc1)[C@@H](O)c1cnccc1Cl>>CN(C)c1ccncc1[…` | `` |
| 1067 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1.Cl` | `` |
| 1070 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.CC(C)(C)OC(=O)NC1CCN(->[Zn+2]([Cl-])([Cl-])<-N2CCC(NC(=O)OC(C)(C)C)C2)C1>>CC(C)(C)OC(…` | `` |
| 1073 | `C_N_Coupling` | `Unknown` | `CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C.OC(c1cnccc1Cl)C(c1ccccc1)(c1ccccc1)c1ccccc1>>CN(C)c1ccncc1C(O)C(…` | `` |
| 1074 | `C_N_Coupling` | `Unknown` | `Clc1ccncc1-c1c(-c2ccccc2)ccc2ccccc12.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1-c1c(-c2ccccc…` | `` |
| 1075 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.C[C@H](N->[Zn+2]([Cl-])([Cl-])<-N[C@@H](C)C1=CC=CC=C1)C1=CC=CC=C1>>C[C@H](Nc1ccncc1)c…` | `` |
| 23043 | `C_N_Coupling` | `Unknown` | `O=Cc1c(Cl)c2ccccc2oc1=O.Nc1ccncc1.CCO[PH](=O)OCC>>CCOP(OCC)(=O->[Zn+2]([Cl-])([Cl-])[Cl-])C(O)C1=C(…` | `` |

### reaction_key_empty_or_degenerate (2)

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | Reaction_Key / Reaction_Events (preview) |
|---:|---|---|---|---|
| 8117 | `C_N_Coupling` | `Multi-Event:Ann+C-N` | `COc1ccc(NP(=O)(OCC/C=C/CCOC(=O)OC(C)(C)C)Oc2ccc(Cl)c3cccnc23)cc1>>COc1ccc(N2[C@@H](/C=C/COC(=O)OC(C…` | `|[] -> [] | spectators: Inorganic` |
| 16541 | `C_N_Coupling` | `Unresolved:NoMotifEvidence` | `COc1ccc(Nc2ccc(OC)cc2)cc1.O=c1c(O)c(O)c1=O>>COC1=CC=C(N(C2=CC=C(OC)C=C2)C2=C(O)C(=[NH](C3=CC=C(OC)C…` | `|[] -> []` |

## Notes

- Most problematic entries are `Unknown` detection labels; unresolved explicit labels are rare.
- Missing `Reaction_Key` and missing `Reaction_Events` occur on the same 7 rows in this dataset.
- No malformed `reaction_smiles` rows were found by the `>>` checks.
