# C_N_Coupling_canonical.csv Problematic Rows by Suggested Fix Action

- Dataset: `data/HTE_db/literature/C_N_Coupling_canonical.csv`
- Generated: 2026-02-15 12:28:55
- Total problematic rows classified: 220

## Action Summary

| Action Bucket | Count | With explicit salt/counterion patterns |
|---|---:|---:|
| `taxonomy_expand_hydrazone_diazene_like` | 99 | 0 |
| `taxonomy_general_unknown` | 66 | 0 |
| `taxonomy_expand_unclassified_reactant` | 43 | 0 |
| `smiles_coordination_notation_cleanup` | 7 | 3 |
| `manual_curation_low_evidence` | 5 | 0 |

## Recommended Actions

### Taxonomy expansion (hydrazone/diazene-like N-N outcomes) (`taxonomy_expand_hydrazone_diazene_like`) - 99 rows

- Recommendation: Add reaction-family coverage and motifs for N-N forming pathways currently mapped as `AromN-H` style unknowns.

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | reaction_key/reaction_events (preview) |
|---:|---|---|---|---|
| 344 | `C_N_Coupling` | `Unknown` | `Fc1cc(N2CCOCC2)cc(F)n1>>NN=c1cc(N2CCOCC2)cc(F)[nH]1` | `|HeteroAr-F -> AromN-H | spectators: Ar-NR2|HeteroAr-H|Pyridine|RCH2-OR` |
| 421 | `C_N_Coupling` | `Unknown` | `COc1cc2c(cc1OC)-c1c(cnc(Cl)c1C#N)C2>>COc1cc2c(cc1OC)-c1c(c[nH]c(=NN)c1C#N)C2` | `|HeteroAr-Cl -> AromN-H | spectators: Ar-Alkyl|Ar-Ar|Ar-OR|HeteroAr-CN|HeteroAr-H|Pyridine` |
| 480 | `C_N_Coupling` | `Unknown` | `Clc1ncnc2sc3c(c12)CCCC3>>NN=c1nc[nH]c2sc3c(c12)CCCC3` | `|HeteroAr-Cl -> AromN-H | spectators: Ar-Alkyl|HeteroAr-H|Pyrimidine|Thiophene` |
| 481 | `C_N_Coupling` | `Unknown` | `Clc1ncnc2sc3c(c12)CCC3>>NN=c1nc[nH]c2sc3c(c12)CCC3` | `|HeteroAr-Cl -> AromN-H | spectators: Ar-Alkyl|HeteroAr-H|Pyrimidine|Thiophene` |
| 806 | `C_N_Coupling` | `Unknown` | `Cn1cc(-c2cnc(F)c(F)c2)cn1>>Cn1cc(-c2c[nH]c(=NN)c(F)c2)cn1` | `|HeteroAr-F -> AromN-H | spectators: Alkyl-AromN|Ar-Ar|HeteroAr-H|Pyrazole|Pyridine` |
| 1906 | `C_N_Coupling` | `Unknown` | `Nc1c(Cl)nnc(Cl)c1[N+](=O)[O-]>>NN=c1[nH]nc(Cl)c(N)c1[N+](=O)[O-]` | `|HeteroAr-Cl -> AromN-H | spectators: Ar-NO2|HeteroAr-NH2` |
| 1917 | `C_N_Coupling` | `Unknown` | `Nc1nc(Cl)nc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O>>NN=c1nc(N)c2ncn([C@@H]3O[C@H](CO)[C@@H](O)[C@H…` | `|HeteroAr-Cl -> AromN-H | spectators: Alkyl-AromN|HeteroAr-H|HeteroAr-NH2|Pyrimidine|R2CH-OH|R2CH-OR|RCH2-OH` |
| 2515 | `C_N_Coupling` | `Unknown` | `NS(=O)(=O)c1ccc(-n2cc(-c3ccc(Cl)cc3)c3c(Cl)ncnc32)cc1>>NN=c1nc[nH]c2c1c(-c1ccc(Cl)cc1)cn2-c1ccc(S(N…` | `|HeteroAr-Cl -> AromN-H | spectators: Ar-Ar|Ar-AromN|Ar-Cl|Ar-SO2NH2|HeteroAr-H|Pyrimidine` |
| 2898 | `C_N_Coupling` | `Unknown` | `Clc1cnc2ccccc2n1>>NN=c1cnc2ccccc2[nH]1` | `|HeteroAr-Cl -> AromN-H | spectators: HeteroAr-H|Pyrimidine` |
| 3508 | `C_N_Coupling` | `Unknown` | `Clc1cncc(Cl)n1>>NN=c1cncc(Cl)[nH]1` | `|HeteroAr-Cl -> AromN-H | spectators: HeteroAr-H|Pyrimidine` |

### General taxonomy unknowns (`taxonomy_general_unknown`) - 66 rows

- Recommendation: Queue remaining unknowns for taxonomy triage and incremental motif/rule additions.

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | reaction_key/reaction_events (preview) |
|---:|---|---|---|---|
| 2389 | `C_N_Coupling` | `Unknown` | `C1COCCOCCNCCOCCOCCN1.Clp1oc2ccc3ccccc3c2c2c(ccc3ccccc32)o1>>c1ccc2c(c1)ccc1op(N3CCOCCOCCN(p4oc5ccc6…` | `|RCH2-NHR -> RCH2-NR2 | spectators: HeteroAr-H|RCH2-OR` |
| 2390 | `C_N_Coupling` | `Unknown` | `C1COCCNCCOCCOCCOCCN1.Clp1oc2ccc3ccccc3c2c2c(ccc3ccccc32)o1>>c1ccc2c(c1)ccc1op(N3CCOCCOCCOCCN(p4oc5c…` | `|RCH2-NHR -> RCH2-NR2 | spectators: HeteroAr-H|RCH2-OR` |
| 2391 | `C_N_Coupling` | `Unknown` | `C1CNCCOCCOCCOCCOCCN1.Clp1oc2ccc3ccccc3c2c2c(ccc3ccccc32)o1>>c1ccc2c(c1)ccc1op(N3CCOCCOCCOCCOCCN(p4o…` | `|RCH2-NHR -> RCH2-NR2 | spectators: HeteroAr-H|RCH2-OR` |
| 2759 | `C_N_Coupling` | `Unknown` | `Clc1nc2ccccc2s1.Cn1ccnc1>>CN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1.[F-][P+5]([F-])([F-])([F-])([F-])[F-]` | `|HeteroAr-Cl -> Ar-AromN | spectators: Alkyl-AromN|HeteroAr-H|Thiazole` |
| 2760 | `C_N_Coupling` | `Unknown` | `Clc1nc2ccccc2s1.CC(C)n1ccnc1>>CC(C)N1C=C[N+](C2=NC3=CC=CC=C3S2)=C1.[F-][P+5]([F-])([F-])([F-])([F-]…` | `|HeteroAr-Cl -> Ar-AromN | spectators: Alkyl-AromN|HeteroAr-H|Thiazole` |
| 2761 | `C_N_Coupling` | `Unknown` | `Clc1nc2ccccc2s1.CCCCn1ccnc1>>CCCCN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1.[F-][P+5]([F-])([F-])([F-])([F-])[…` | `|HeteroAr-Cl -> Ar-AromN | spectators: Alkyl-AromN|HeteroAr-H|Thiazole` |
| 2762 | `C_N_Coupling` | `Unknown` | `Clc1nc2ccccc2s1.c1ccc(Cn2ccnc2)cc1>>C1=CC=C(CN2C=C[N+](C3=NC4=CC=CC=C4S3)=C2)C=C1.[F-][P+5]([F-])([…` | `|HeteroAr-Cl -> Ar-AromN | spectators: Alkyl-AromN|Ar-Alkyl|HeteroAr-H|Thiazole` |
| 2763 | `C_N_Coupling` | `Unknown` | `Clc1nc2ccccc2s1.Cc1c(C)c(C)c(Cn2ccnc2)c(C)c1C>>CC1=C(C)C(C)=C(CN2C=C[N+](C3=NC4=CC=CC=C4S3)=C2)C(C)…` | `|HeteroAr-Cl -> Ar-AromN | spectators: Alkyl-AromN|Ar-Alkyl|HeteroAr-H|Thiazole` |
| 4121 | `C_N_Coupling` | `Unknown` | `Clc1ccccn1.[N-]=[N+]=NCc1ccc(C2=NC2)cc1>>[N-]=[N+]=NCc1ccc(-c2cnc3ccccn23)cc1` | `|Ar-C=N|HeteroAr-Cl -> Ar-Ar|HeteroAr-H | spectators: Bn-N3|HeteroAr-H|Pyridine` |
| 4122 | `C_N_Coupling` | `Unknown` | `Brc1ccc(C2=NC2)cc1.Clc1ncccc1Br>>Brc1ccc(-c2cnc3c(Br)cccn23)cc1` | `|Ar-C=N|HeteroAr-Cl -> Ar-Ar|HeteroAr-H | spectators: Ar-Br|HeteroAr-Br|HeteroAr-H|Pyridine` |

### Taxonomy expansion (Unclassified-Reactant) (`taxonomy_expand_unclassified_reactant`) - 43 rows

- Recommendation: Add missing reactant motifs / improve reactant typing for currently unclassified partners.

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | reaction_key/reaction_events (preview) |
|---:|---|---|---|---|
| 2324 | `C_N_Coupling` | `Unknown` | `Cc1cccc2c1NCC2.COC1=CC(C(N)=O)=C([I+]C2=C(C)C=C(C)C=C2C)C(OC)=C1.[F-][B+3]([F-])([F-])[F-]>>COc1cc(…` | `|Ar-Alkyl|Ar-Iodonium|Ar-NHR|Unclassified-Reactant -> Ar-NR2 | spectators: Ar-CONH2|Ar-OR` |
| 2325 | `C_N_Coupling` | `Unknown` | `Cc1cccc2c1NCC21CCCCC1.COC1=CC(C(N)=O)=C([I+]C2=C(C)C=C(C)C=C2C)C(OC)=C1.[F-][B+3]([F-])([F-])[F-]>>…` | `|Ar-Alkyl|Ar-Iodonium|Ar-NHR|Unclassified-Reactant -> Ar-NR2 | spectators: Ar-CONH2|Ar-OR` |
| 2326 | `C_N_Coupling` | `Unknown` | `Cc1cccc2c1NCC2.CC1=CC(C)=C([I+]C2=C(OC(C)C)C=C(OC(C)C)C=C2C(N)=O)C(C)=C1.[F-][B+3]([F-])([F-])[F-]>…` | `|Ar-Alkyl|Ar-Iodonium|Ar-NHR|Unclassified-Reactant -> Ar-NR2 | spectators: Ar-CONH2|Ar-OR` |
| 2483 | `C_N_Coupling` | `Unknown` | `O=c1[nH][nH]c(=O)c2ccccc12.CC1=CC=C2C(=C1)[I+]C1=CC=CC=C12.[F-][B+3]([F-])([F-])[F-]>>Cc1ccc2c3cccc…` | `|Ar-Ar|Ar-Iodonium|AromN-H|Unclassified-Reactant -> HeteroAr-H | spectators: Ar-Alkyl|HeteroAr-H` |
| 3070 | `C_N_Coupling` | `Unknown` | `C1=CC=C([I+]C2=CC=CC=C2)C=C1.[F-][B+3]([F-])([F-])[F-].O=c1[nH]ccc2ccnc(Cl)c12>>O=c1c2c(Cl)nccc2ccn…` | `|Ar-Iodonium|AromN-H|Unclassified-Reactant -> Ar-AromN | spectators: HeteroAr-Cl|HeteroAr-H|Pyridine` |
| 3071 | `C_N_Coupling` | `Unknown` | `O=c1[nH]ccc2ccnc(Cl)c12.CCC1=CC=C([I+]C2=CC=C(CC)C=C2)C=C1.[F-][B+3]([F-])([F-])[F-]>>CCc1ccc(-n2cc…` | `|Ar-Alkyl|Ar-Iodonium|AromN-H|Unclassified-Reactant -> Ar-AromN | spectators: HeteroAr-Cl|HeteroAr-H|Pyridine` |
| 3072 | `C_N_Coupling` | `Unknown` | `CC(C)C1=CC=C([I+]C2=CC=C(C(C)C)C=C2)C=C1.[F-][B+3]([F-])([F-])[F-].O=c1[nH]ccc2ccnc(Cl)c12>>CC(C)c1…` | `|Ar-Alkyl|Ar-Iodonium|AromN-H|Unclassified-Reactant -> Ar-AromN | spectators: HeteroAr-Cl|HeteroAr-H|Pyridine` |
| 3073 | `C_N_Coupling` | `Unknown` | `CC(C)(C)C1=CC=C([I+]C2=CC=C(C(C)(C)C)C=C2)C=C1.[F-][B+3]([F-])([F-])[F-].O=c1[nH]ccc2ccnc(Cl)c12>>C…` | `|Ar-Alkyl|Ar-Iodonium|AromN-H|Unclassified-Reactant -> Ar-AromN | spectators: HeteroAr-Cl|HeteroAr-H|Pyridine` |
| 3074 | `C_N_Coupling` | `Unknown` | `CC1=CC=C([I+]C2=CC=C(C)C=C2)C=C1.[F-][B+3]([F-])([F-])[F-].O=c1[nH]ccc2ccnc(Cl)c12>>Cc1ccc(-n2ccc3c…` | `|Ar-Alkyl|Ar-Iodonium|AromN-H|Unclassified-Reactant -> Ar-AromN | spectators: HeteroAr-Cl|HeteroAr-H|Pyridine` |
| 3075 | `C_N_Coupling` | `Unknown` | `COC1=CC=C([I+]C2=CC=C(OC)C=C2)C=C1.[F-][B+3]([F-])([F-])[F-].O=c1[nH]ccc2ccnc(Cl)c12>>COc1ccc(-n2cc…` | `|Ar-Iodonium|Ar-OR|AromN-H|Unclassified-Reactant -> Ar-AromN | spectators: HeteroAr-Cl|HeteroAr-H|Pyridine` |

### SMILES normalization cleanup (coordination/dative notation) (`smiles_coordination_notation_cleanup`) - 7 rows

- Recommendation: Normalize/remove dative bond syntax (`->`, `<-`) and re-standardize metal complexes before featurization.

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | reaction_key/reaction_events (preview) |
|---:|---|---|---|---|
| 1066 | `C_N_Coupling` | `Unknown` | `CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C.CC(C)(C)[C@H](NC(=O)c1ccccc1)[C@@H](O)c1cnccc1Cl>>CN(C)c1ccncc1[…` | `` |
| 1067 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1.Cl` | `` |
| 1070 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.CC(C)(C)OC(=O)NC1CCN(->[Zn+2]([Cl-])([Cl-])<-N2CCC(NC(=O)OC(C)(C)C)C2)C1>>CC(C)(C)OC(…` | `` |
| 1073 | `C_N_Coupling` | `Unknown` | `CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C.OC(c1cnccc1Cl)C(c1ccccc1)(c1ccccc1)c1ccccc1>>CN(C)c1ccncc1C(O)C(…` | `` |
| 1074 | `C_N_Coupling` | `Unknown` | `Clc1ccncc1-c1c(-c2ccccc2)ccc2ccccc12.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1-c1c(-c2ccccc…` | `` |
| 1075 | `C_N_Coupling` | `Unknown` | `Cl.Clc1ccncc1.C[C@H](N->[Zn+2]([Cl-])([Cl-])<-N[C@@H](C)C1=CC=CC=C1)C1=CC=CC=C1>>C[C@H](Nc1ccncc1)c…` | `` |
| 23043 | `C_N_Coupling` | `Unknown` | `O=Cc1c(Cl)c2ccccc2oc1=O.Nc1ccncc1.CCO[PH](=O)OCC>>CCOP(OCC)(=O->[Zn+2]([Cl-])([Cl-])[Cl-])C(O)C1=C(…` | `` |

### Manual curation (low-evidence / degenerate key) (`manual_curation_low_evidence`) - 5 rows

- Recommendation: Manually inspect rows with `-> []` or empty formed motifs; many need data correction or explicit exclusions.

| Row | reaction_id | detected_reaction_type | reaction_smiles (preview) | reaction_key/reaction_events (preview) |
|---:|---|---|---|---|
| 165 | `C_N_Coupling` | `Unresolved:NoMotifEvidence` | `NC=O.Cc1cc(Br)c2c(c1)C(C)(C)c1cc(C)cc(Br)c1O2>>Cc1cc(NC=O)c2c(c1)C(C)(C)c1cc(C)cc(NC=O)c1O2` | `|Ar-Br|Unclassified-Reactant -> [] | spectators: Ar-Alkyl|Ar-OR` |
| 5183 | `C_N_Coupling` | `Unresolved:NoMotifEvidence` | `COc1ccc2nc(Cl)c(N3CCN(C(=O)OC(C)(C)C)CC3)nc2c1.O=C(Nc1ncc([N+](=O)[O-])s1)N1CCNCC1>>COc1ccc2nc(N3CC…` | `|RCH2-NHR -> [] | spectators: Ar-NO2|Ar-NR2|Ar-Urea|HeteroAr-Cl|HeteroAr-H|HeteroAr-OR|Pyrimidine|Thiazole` |
| 8117 | `C_N_Coupling` | `Multi-Event:Ann+C-N` | `COc1ccc(NP(=O)(OCC/C=C/CCOC(=O)OC(C)(C)C)Oc2ccc(Cl)c3cccnc23)cc1>>COc1ccc(N2[C@@H](/C=C/COC(=O)OC(C…` | `|[] -> [] | spectators: Inorganic` |
| 16541 | `C_N_Coupling` | `Unresolved:NoMotifEvidence` | `COc1ccc(Nc2ccc(OC)cc2)cc1.O=c1c(O)c(O)c1=O>>COC1=CC=C(N(C2=CC=C(OC)C=C2)C2=C(O)C(=[NH](C3=CC=C(OC)C…` | `|[] -> []` |
| 24077 | `C_N_Coupling` | `Unresolved:NoMotifEvidence` | `CN1CCOCC1.COc1nc(Cl)nc(OC)n1>>COC1=NC(OC)=NC([N+]2(C)CCOCC2)=N1.[F-][B+3]([F-])([F-])[F-]` | `|HeteroAr-Cl|RCH2-NR2 -> [] | spectators: HeteroAr-OR|RCH2-OR` |

## Detection Reason Breakdown

| Reason token | Count |
|---|---:|
| `nn_like_unknown` | 99 |
| `unknown_generic` | 66 |
| `unclassified_reactant` | 43 |
| `coordination_notation` | 7 |
| `degenerate_or_unresolved` | 5 |

## Notes

- This is a triage report: each problematic row is assigned a single primary action bucket by priority.
- Priorities used: coordination-notation cleanup > missing outputs recompute > unresolved/degenerate manual curation > taxonomy expansion buckets.
