# Chemist-oriented weak-label condition grouping POC

The learned identity covers decisive condition cores. Solvents, additives,
temperature, and time remain attached as cross-referenced protocol context.
Reaction labels and outcomes were excluded from clustering.

## Dataset

- Rows: 59,277
- Exact recipe IDs: 44,977
- Complete material systems: 25,047
- Exact decisive cores: 3,801
- Material systems without a decisive core: 220
- Input SHA-256: `a74f864d7d2a8f980405be977befc48bfcd430e5381dd03e5e6edf621ee326c2`

## Model

- Requested groups: 256
- Populated groups: 256
- TF-IDF features: 993
- Latent dimensions: 48
- Seed: 17

## Diagnostics

- Sampled silhouette score: 0.56231
- Mean centroid similarity: 0.973684
- Median assignment margin: 0.078458
- Median exact cores per group: 5.0
- Median material variants per group: 46.0
- Median observations per group: 122.5
- Supported core assignments: 2,692
- Provisional core assignments: 704
- Review core assignments: 405

## Caveats

- `CGPOC2` identifiers are bound to this dataset snapshot and model settings.
- Groups are statistical proposals, not reagent-interchangeability claims.
- Additives are contextual in this POC even though some are mechanistically critical.
- Production identities require chemistry review and frozen definitions.

## Largest groups

| Group | Rows | Cores | Material variants | Prototype core |
| --- | ---: | ---: | ---: | --- |
| `CGPOC2:ce467fcd32aeef0c28fc89ecf7a6ca369c1501ca5ac21a64c6bed96eb59dc312` | 3,086 | 43 | 1,306 | Cs2CO3 [base] |
| `CGPOC2:8246b1f2de5d1a8fa561bb0c59f5b27c9119e5e86980f24229ad3bac6d685e17` | 1,688 | 42 | 460 | K3PO4 [base]; tBuBrettPhos Pd(allyl)OTf [catalyst]; tBuBrettPhos [ligand] |
| `CGPOC2:9635021a0ac1fb2366d1eb1987b3c6ba1287e66c0c032f457cf8cf1d4f4fcee9` | 1,286 | 46 | 520 | K3PO4 [base]; dppfPdCl2 [catalyst]; dppf [ligand] |
| `CGPOC2:f635431829dc4e78e98f6e38819be685db8e122946602936cbfb40e991bf4c2d` | 1,243 | 129 | 567 | NMI [base]; DIC [coupling_reagent] |
| `CGPOC2:f8eaa3867ad073d9f2b1e7e9767343e6d8a52a6367822bfe9785c00927981066` | 1,144 | 39 | 445 | Na2CO3 [base]; dtbpfPdCl2 [catalyst]; dtbpf [ligand] |
| `CGPOC2:a055495df2bc6449ab288cfd1597e765727922724ff49dad604d1dbfa3da6dc2` | 1,005 | 50 | 345 | Cs2CO3 [base]; CuI [catalyst]; trans-NN-Dimethylcyclohexane-12-diamine [ligand] |
| `CGPOC2:1ec74c098ac713d9287824e3c9e2e5e86a58a7ff5805076ad2803199d9f377a4` | 999 | 42 | 246 | Cs2CO3 [base]; QPhos Pd(crotyl)Cl [catalyst]; QPhos [ligand] |
| `CGPOC2:84d94085b40f177c6a4e85485f2fa48e85acdaf9902763e603f976d8c73dac52` | 962 | 31 | 367 | K3PO4 [base]; P(tBu)3 Pd(crotyl)Cl [catalyst]; P(tBu)3 [ligand] |
| `CGPOC2:b70c9f6b0bcb3a4b68950e1dd6741db5b8cc15a8808dd08b4f3ddfed921e74d5` | 949 | 36 | 378 | K3PO4 [base]; (DiMeIHept Cl) Pd(cinnamyl)Cl [catalyst]; DiMeIHept Cl [ligand] |
| `CGPOC2:b24ec5ae0196c1c3d892e95aedb39af784b69ed5f910933312e8d8212381e471` | 908 | 48 | 314 | Cs2CO3 [base]; tBuXPhos Pd(allyl)OTf [catalyst]; Me4tBuXPhos [ligand] |
| `CGPOC2:453ed8394a86f21ef59ecd576761b19ad0f0347f0086d013942c99d7b212c731` | 891 | 44 | 356 | VPhos Pd G4 [catalyst]; VPhos [ligand] |
| `CGPOC2:04c858ed0e12dafc25ce9f74b60c8afc9124ed2ef16c812ee09ed28a29adaa2d` | 843 | 40 | 308 | Cs2CO3 [base]; JosiphosSL-J009-1 Pd G3 [catalyst]; JosiphosSL-J009-1 [ligand] |
| `CGPOC2:46da7957c67de8958a6e51e310e5b61016945d9aa41d29edf52f74889bfbd3ab` | 792 | 29 | 447 | K2CO3 [base]; meCgPPh Pd G3 [catalyst]; meCgPPh [ligand] |
| `CGPOC2:986f20cb7efae48872337388fea5567228834a3c860006dd869727f69ccca6bc` | 792 | 22 | 316 | DIC [coupling_reagent] |
| `CGPOC2:2e34d0b79b1b21091394cb3c8a4368ec5d7862e4b75afd35bcab586ddbc95e53` | 785 | 40 | 229 | Cs2CO3 [base]; BrettPhos Pd(crotyl)OTf [catalyst]; BrettPhos [ligand] |
| `CGPOC2:896899974226b9d60c74658f7766721fe56cfc056c552459715f1e26dd8c4e2b` | 765 | 41 | 239 | Cs2CO3 [base]; BippyPhos Pd(allyl)OTf [catalyst]; BippyPhos [ligand] |
| `CGPOC2:558a3288cd6a1a4d855647238ed801e8f966ccb98b412effd2e32e070dc57981` | 764 | 34 | 259 | K3PO4 [base]; cBridP Pd(allyl)Cl [catalyst]; cBridP [ligand] |
| `CGPOC2:4c95ef950455ea61fc7987b7a24c83346dacdcfadfd90e602373780ba31c08c5` | 756 | 31 | 242 | Cs2CO3 [base]; SPhos Pd(crotyl)Cl [catalyst]; SPhos [ligand] |
| `CGPOC2:bf914a5022d53b58905323554be7b21238b0a442b46fb70a1be41887d75d92e9` | 733 | 58 | 364 | (Ru(p-cymene)Cl2)2 [catalyst]; No ligand |
| `CGPOC2:4b31d62d21f65c3313b33bdf5a18424f9cf4166e8c37ec1dff6fb0045e6c8f73` | 719 | 38 | 441 | NMM [base]; DIC [coupling_reagent] |
