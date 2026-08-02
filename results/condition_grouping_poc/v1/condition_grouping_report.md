# Weak-label condition grouping POC

This experimental model groups role-aware material systems only. Reaction
labels, yields, temperature, and time were excluded from clustering.

## Dataset

- Rows: 59,277
- Exact recipe IDs: 44,977
- Exact normalized material systems: 25,047
- Input SHA-256: `a74f864d7d2a8f980405be977befc48bfcd430e5381dd03e5e6edf621ee326c2`

## Model

- Requested groups: 384
- Populated groups: 384
- TF-IDF features: 1,284
- Latent dimensions: 48
- Seed: 17

## Diagnostics

- Sampled silhouette score: 0.311056
- Mean centroid similarity: 0.927297
- Median assignment margin: 0.078871
- Median material systems per group: 49.0
- Median observations per group: 111.0
- Singleton groups: 1
- Supported assignments: 14,833
- Provisional assignments: 8,135
- Review assignments: 2,079

## Caveats

- `CGPOC1` identifiers are bound to this dataset snapshot and model settings.
- Groups are statistical proposals, not claims of reagent interchangeability.
- Reaction-type purity is an audit metric only; reaction labels were not features.
- Production identities require chemistry review and frozen definitions.

## Largest groups

| Group | Rows | Material systems | Prototype | Mean similarity |
| --- | ---: | ---: | --- | ---: |
| `CGPOC1:13080393c6a5e2957f82b3c1ac3fe42c3c1c32cf2ca7c8886a6807bad7786044` | 1,098 | 437 | CsOH.H2O [base]; NBP [solvent]; No ligand | 0.9520 |
| `CGPOC1:d5204693db07c70118c0d170cdff2ba90e2cb16bd53eba8ee257da25a2f958a0` | 967 | 541 | HOPO [additive]; NMI [base]; CDT [coupling_reagent]; acetone [solvent]; No ligand | 0.9391 |
| `CGPOC1:5fe5e2ff8917bcf819356d850eb741feb2f7f11c4d4b66d7b361716c7c6b54c4` | 867 | 165 | NaBHT [base]; QPhos Pd(crotyl)Cl [catalyst]; QPhos [ligand]; TBME [solvent] | 0.9456 |
| `CGPOC1:a334b814944448a5a95c0249ca72d3f4555d0781b171b8cf13b02d5850feafb1` | 811 | 175 | CsF [base]; tBuBrettPhos Pd(allyl)OTf [catalyst]; tBuBrettPhos [ligand]; PhCN [solvent] | 0.9438 |
| `CGPOC1:ced0d9c9aa264a48e4dbf22f5c515d59a69574f9a83083a06c636e80fb4b49d6` | 726 | 205 | CsF [base]; (R)-BINAP Pd(allyl)Cl [catalyst]; BINAP [ligand]; PhCN [solvent] | 0.9224 |
| `CGPOC1:504b1359f9d0e17c0514c964cce7c945938957c0bbd6967e3e545990490eed27` | 608 | 200 | N-Methyldicyclohexylamin [base]; tBuXPhos Pd(allyl)OTf [catalyst]; Me4tBuXPhos [ligand]; ethylene glycol [solvent] | 0.9434 |
| `CGPOC1:38d10d8ae0114a70ac3984c91278c47208ff5b9d365a7df6a91a65331c7554c5` | 604 | 211 | LiOH anh [base]; cBridP Pd(allyl)Cl [catalyst]; cBridP [ligand]; CPME [solvent] | 0.9229 |
| `CGPOC1:790925c3f0ac14d309541f4422607c04817106804bdd6a7d731e01775c06a3a6` | 598 | 164 | KOtPent [base]; BrettPhos Pd(crotyl)OTf [catalyst]; BrettPhos [ligand]; Ph2O [solvent] | 0.9313 |
| `CGPOC1:dfdf0e268de656417e9e0d0cbda8c92340b70d8009fc4b6c427aa1d66627f200` | 598 | 180 | CsF [base]; XantPhos Pd(allyl)Cl [catalyst]; XantPhos [ligand]; PhCN [solvent] | 0.9393 |
| `CGPOC1:69d567895a18e9ccf6564c6b6e9ba09ae46bc7e56c2c3b0399b79f004ed7dfe0` | 590 | 156 | CsF [base]; JosiphosSL-J009-1 Pd G3 [catalyst]; JosiphosSL-J009-1 [ligand]; PhCN [solvent] | 0.9494 |
| `CGPOC1:82ac41fd619448bdb5d6f7d89ceabcdb8586f4d1d87cf3a00df6daf606693fdd` | 527 | 183 | N-Methyldicyclohexylamin [base]; AlphosPd)2cod [catalyst]; AlPhos [ligand]; ethylene glycol [solvent] | 0.9309 |
| `CGPOC1:1ddfcf4a4fda9e1422429971b6ad0f21d85c33c6fadf050a2390363fcca26d35` | 494 | 174 | KOtPent [base]; (DiMeIHept Cl) Pd(cinnamyl)Cl [catalyst]; DiMeIHept Cl [ligand]; Ph2O [solvent] | 0.9452 |
| `CGPOC1:26dbd3f1556df3b1b420605138035726c6a2b64a6e574640d24f3229c2f16b92` | 492 | 137 | KOtPent [base]; BippyPhos Pd(allyl)OTf [catalyst]; BippyPhos [ligand]; Ph2O [solvent] | 0.9500 |
| `CGPOC1:706e29d9879b98739ccd1b0551da103f95d3580356bfcab56a47761467ac2a48` | 488 | 199 | BSTFA [additive]; HOSu [additive]; NMM [base]; DIC [coupling_reagent]; MeCN [solvent]; No ligand | 0.9684 |
| `CGPOC1:ae683d0f2f1e05a90c5a7ea1f6f4685b5b39a93e2a5caeac9760894a788a0c5b` | 488 | 332 | Li2CO3 [base]; T3P [coupling_reagent]; acetone [solvent]; No ligand | 0.9608 |
| `CGPOC1:3f05406bfa64bbf62ce14915861bf5be803fed66538900ddb6226050ad4d118b` | 472 | 123 | NMI, NMM [base]; T3P [coupling_reagent]; MeCN [solvent]; No ligand | 0.9797 |
| `CGPOC1:6710d1e5ff9c50fac49c17d9999345264b7a2e0c4d62a814181df81cf5c4f7f9` | 447 | 147 | NaOtPent [base]; pyBOP [coupling_reagent]; THF [solvent]; No ligand | 0.9797 |
| `CGPOC1:e5adf913b0f1c9a75358706dce540cabdfc8b71eeec6dca7159242cc4bcfca71` | 447 | 83 | NaOH [base]; XantPhos Pd(allyl)Cl [catalyst]; XantPhos [ligand]; iPrOH [solvent]; water [solvent] | 0.9139 |
| `CGPOC1:c1f795b902e6aba4050a4b470efaecca4fe8c43ff1e6511d5603953e9d9fe183` | 443 | 124 | K2HPO4 [base]; AdBrettPhosPd G6 Br [catalyst]; Ad-BrettPhos [ligand]; TPGS-750M [solvent] | 0.8940 |
| `CGPOC1:943fca0451c24dc75c214d7c822812357ae9dd21dbeeb542e5c40343e12f2d5d` | 442 | 93 | NaHCO3 [base]; SPhos Pd(crotyl)Cl [catalyst]; SPhos [ligand]; TBME [solvent]; water [solvent] | 0.9140 |
