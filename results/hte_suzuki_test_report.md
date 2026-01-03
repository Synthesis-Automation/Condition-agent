# HTE Recommendation Test Report: Suzuki Reactions
Source: examples\sample_reactions.csv
Database: data\HTE_db\HTE_canonical.csv

| ID | Example | Reactants | Predicted Type | Matches | Top Condition (Catalyst/Ligand/Base/Solvent) | Avg Z-Score |
|---|---|---|---|---|---|---|
| 0 | Simple Ph-Ph | `Brc1ccccc1.c1ccc(B(O)O)cc1` | suzuki_miyaura | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 1 | Electron-poor ArCl | `Clc1ccc(C#N)cc1.c1ccc(B(O)O)cc1` | None | 960 | P(tBu)2Ph)2PdCl2 / P(tBu)2Ph / K2CO3 / Dioxane | 2.45 |
| 2 | Electron-rich ArBr | `Brc1ccc(OC)cc1.c1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 3 | Heteroaryl pyridine | `Ic1ccncc1.c1ccc(B(O)O)cc1` | None | 0 | N/A | N/A |
| 4 | CF3 substrate | `Brc1ccc(C(F)(F)F)cc1.c1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 5 | Naphthyl chloride | `Clc1ccc2ccccc2c1.c1ccc(B(O)O)cc1` | suzuki_miyaura | 960 | P(tBu)2Ph)2PdCl2 / P(tBu)2Ph / K2CO3 / Dioxane | 2.45 |
| 6 | Sterically hindered | `Brc1cc(C)ccc1C.c1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 7 | Aryl triflate electrophile | `Fc1ccc(OS(=O)(=O)C(F)(F)F)cc1.c1ccc(B(O)O)cc1` | None | 0 | N/A | N/A |
| 8 | Vinyl bromide + vinyl boronic acid | `C/C=C/Br.C=CB(O)O` | None | 0 | N/A | N/A |
| 9 | Potassium aryltrifluoroborate | `Brc1ccc(OC)cc1.[K+].[B-](F)(F)c2ccccc2` | None | 4621 | dtbpfPdCl2 / dtbpf / K3PO4 / Dioxane, water | 5.24 |
| 10 | 2-pyridyl MIDA slow-release | `Clc1ccncc1.c1ccncc1B1OC(=O)CN(CC(=O)O)C(=O)CN(CC(=O)O)C(=O)O1` | None | 2656 | cataCXiumA Pd G6 Br)2 / CataCXiumA / NaHCO3 / tAmOH, water | 2.11 |
| 11 | Dichloropyridine | `Brc1ccc(Cl)nc1.c2ccc(B(O)O)cc2` | None | 960 | P(tBu)2Ph)2PdCl2 / P(tBu)2Ph / K2CO3 / Dioxane | 2.45 |
| 12 | Mixed halide pyridine | `Brc1cc(Cl)nc(Cl)c1.c2ccc(B(O)O)nc2` | None | 960 | P(tBu)2Ph)2PdCl2 / P(tBu)2Ph / K2CO3 / Dioxane | 2.45 |
| 13 | Indole heterobiaryl | `Brc1cccc2[nH]ccc12.c3ccc(B(O)O)cc3` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 14 | Trifluoromethyl aryl chloride + pyrimidine boronate | `FC(F)(F)c1ccc(Cl)cc1.Nc1ccnc(B(O)O)n1` | None | 960 | P(tBu)2Ph)2PdCl2 / P(tBu)2Ph / K2CO3 / Dioxane | 2.45 |
| 15 | Chloropyridyl boronic acid with aryl bromide | `Clc1cncc(B(O)O)c1.Brc2ccc(F)cc2` | None | 1032 |  / No Ligand / Na2CO3 / THF | 2.76 |
| 16 | Bis-coupling to biphenyl | `Brc1ccc(Br)cc1.c1ccc(B(O)O)cc1` | suzuki_miyaura | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 17 | Meta-dibromobenzene | `Brc1cccc(Br)c1.c1ccc(B(O)O)cc1` | suzuki_miyaura | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 18 | Pyrimidine-5-boronic acid | `Brc1ccccc1.c1ncc(B(O)O)cn1` | suzuki_miyaura | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 19 | Furan-2-boronic acid | `Clc1ccc(C#N)cc1.c1coc(B(O)O)c1` | None | 960 | P(tBu)2Ph)2PdCl2 / P(tBu)2Ph / K2CO3 / Dioxane | 2.45 |
| 20 | Thiophene-2-boronic acid | `Brc1ccncc1.c1csc(B(O)O)c1` | suzuki_miyaura | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 21 | Pyrrole-3-boronic acid | `Ic1ccc(C(=O)OC)cc1.c1c[nH]c(B(O)O)c1` | None | 408 | Cu(MeCN)4BF4 / trans-NN-Dimethylcyclohexane-12-diamine / Cs2CO3 / CPME | 1.39 |
| 22 | 3-Bromoquinoxaline | `Brc1cnc2ccccc2n1.c1ccc(B(O)O)cc1` | suzuki_miyaura | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 23 | 4-Chloroindole + pyridine boronic acid | `Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1` | None | 960 | P(tBu)2Ph)2PdCl2 / P(tBu)2Ph / K2CO3 / Dioxane | 2.45 |
| 24 | 2-Bromobenzothiazole | `Brc1nc2ccccc2s1.Cc1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 25 | 2-Iodobenzoxazole | `Ic1nc2ccccc2o1.COc1ccc(B(O)O)cc1` | None | 0 | N/A | N/A |
| 26 | Benzothiadiazole | `Brc1cccc2nsnc12.c1ccc(B(O)O)cc1` | suzuki_miyaura | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 27 | 2-Isopropyl-bromobenzene | `Brc1ccccc1C(C)C.c1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 28 | 2,6-Dimethyl-chlorobenzene | `Clc1c(C)cccc1C.COc1ccc(B(O)O)cc1` | None | 960 | P(tBu)2Ph)2PdCl2 / P(tBu)2Ph / K2CO3 / Dioxane | 2.45 |
| 29 | Ortho-ethoxy + ortho-methyl | `Brc1ccccc1OCC.c1ccc(B(O)O)cc1C` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 30 | Pentafluorobromobenzene | `Fc1c(F)c(F)c(Br)c(F)c1F.c1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 31 | 3,5-Dichloro-bromobenzene | `Clc1cc(Cl)cc(Br)c1.c1ccc(B(O)O)cc1` | None | 960 | P(tBu)2Ph)2PdCl2 / P(tBu)2Ph / K2CO3 / Dioxane | 2.45 |
| 32 | 2,5-Dinitro-bromobenzene | `Brc1ccc([N+](=O)[O-])cc1[N+](=O)[O-].c1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 33 | Ethyl ester protected | `Brc1ccc(C(=O)OCC)cc1.c1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 34 | Boc-protected amine | `Ic1ccc(NC(=O)OC(C)(C)C)cc1.c1ccc(B(O)O)cc1` | None | 0 | N/A | N/A |
| 35 | TBS-protected phenol | `Brc1ccc(O[Si](C)(C)C(C)(C)C)cc1.c1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 36 | Vinylboronic acid to styrene | `Brc1ccccc1.C=CB(O)O` | None | 0 | N/A | N/A |
| 37 | (E)-Propenylboronic acid | `Ic1ccc(C=O)cc1.C/C=C/B(O)O` | None | 0 | N/A | N/A |
| 38 | Isopropenylboronic acid | `Brc1ccncc1.C=C(C)B(O)O` | None | 0 | N/A | N/A |
| 39 | Ethynylboronic acid | `Brc1ccc(OC)cc1.C#CB(O)O` | None | 0 | N/A | N/A |
| 40 | Macrocyclization precursor | `Brc1ccc(Br)cc1CCCCCCCC(=O)O.c1ccc(B(O)O)cc1` | None | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 41 | Potassium phenyltrifluoroborate | `Brc1ccccc1.[K+].[B-](F)(F)(F)c1ccccc1` | miyaura_borylation | 4621 | dtbpfPdCl2 / dtbpf / K3PO4 / Dioxane, water | 5.24 |
| 42 | BF3K + ArCl | `Clc1ccc(C(F)(F)F)cc1.[K+].[B-](F)(F)(F)c1ccc(OC)cc1` | None | 1512 | RuPhos Pd(crotyl)Cl / RuPhos / K3PO4 / PhMe, water | 4.92 |
| 43 | Pyridine N-oxide | `[O-][n+]1ccccc1Br.c1ccc(B(O)O)cc1` | suzuki_miyaura | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |
| 44 | Cyclopropyl bromide | `BrC1(c2ccccc2)CC1.c1ccc(B(O)O)cc1` | None | 0 | N/A | N/A |
| 437 | Bromobenzene + phenylboronic acid | `Brc1ccccc1.c1ccc(B(O)O)cc1` | suzuki_miyaura | 2024 | dtbpfPdCl2 / dtbpf / NaHCO3 / THF, water | 2.75 |