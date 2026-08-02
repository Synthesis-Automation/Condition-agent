# Weak-label condition grouping POC

This experimental model groups role-aware material systems only. Reaction
labels, yields, temperature, and time were excluded from clustering.

## Dataset

- Rows: 59,277
- Exact recipe IDs: 44,977
- Exact normalized material systems: 25,047
- Input SHA-256: `a74f864d7d2a8f980405be977befc48bfcd430e5381dd03e5e6edf621ee326c2`

## Model

- Requested groups: 768
- Populated groups: 768
- TF-IDF features: 1,284
- Latent dimensions: 64
- Seed: 17

## Diagnostics

- Sampled silhouette score: 0.280545
- Mean centroid similarity: 0.92785
- Median assignment margin: 0.076227
- Median material systems per group: 25.0
- Median observations per group: 53.5
- Singleton groups: 6

## Caveats

- `CGPOC1` identifiers are bound to this dataset snapshot and model settings.
- Groups are statistical proposals, not claims of reagent interchangeability.
- Reaction-type purity is an audit metric only; reaction labels were not features.
- Production identities require chemistry review and frozen definitions.

## Largest groups

| Group | Rows | Material systems | Prototype | Mean similarity |
| --- | ---: | ---: | --- | ---: |
| `CGPOC1:261dcfecf4a61ea9be3085d13b648e7f812928a6c6a714f00b98448f78585e35` | 611 | 120 | CsF [base]; tBuBrettPhos Pd(allyl)OTf [catalyst]; tBuBrettPhos [ligand]; PhCN [solvent] | 0.9318 |
| `CGPOC1:cedfdf5c5a1e7ca1fbbd9c916674310f45677ccce6bd92eb36033b2004c209b8` | 554 | 168 | TMG [base]; RuPhos Pd(crotyl)Cl [catalyst]; RuPhos [ligand]; TPGS-750M [solvent] | 0.9154 |
| `CGPOC1:9b92c985080bc0d01c90ca4acae2bb127fee0be8ef87b2b047418b3138093f62` | 551 | 153 | CsF [base]; XantPhos Pd(allyl)Cl [catalyst]; XantPhos [ligand]; PhCN [solvent] | 0.9257 |
| `CGPOC1:d33b46595507036d27231e2e62c3577de7bb8eb625e7ade899c847d81ef075ca` | 524 | 258 | Ag2O [base]; EtOAc [solvent]; No ligand | 0.9640 |
| `CGPOC1:8d81fb6917030599f7aadd1c6e2fa25685a6ebfc01771494f59a5d1929e12059` | 405 | 137 | CsF [base]; tBuXPhos Pd(allyl)OTf [catalyst]; Me4tBuXPhos [ligand]; PhCN [solvent] | 0.9407 |
| `CGPOC1:fd77e3952443f00d3b2b29b346ec0cc3301b9d2f2ff91c5ed0ee6bab89cf02cd` | 382 | 122 | KOtPent [base]; tBuXPhos Pd G3 [catalyst]; tBuXPhos [ligand]; Ph2O [solvent] | 0.9443 |
| `CGPOC1:fa772a04d1feebc70aad627aa5795104a63d4b6df5ff5999f0c1fb371a5dc474` | 364 | 110 | CsF [base]; JosiphosSL-J009-1 Pd G3 [catalyst]; JosiphosSL-J009-1 [ligand]; PhCN [solvent] | 0.9360 |
| `CGPOC1:fd4b0e5c8dea22b65ff4cd318f44d5e7b8a0f301252b8f0ce0ec61b6934d7ebd` | 363 | 98 | NaBHT [base]; QPhos Pd(crotyl)Cl [catalyst]; QPhos [ligand]; xylenes [solvent] | 0.9376 |
| `CGPOC1:7a9c25650726b0e56139b5b2d244c9c050b6a50f1ac850d193e177291e962543` | 353 | 32 | NaOtBu [base]; tBuBrettPhos Pd(allyl)OTf [catalyst]; tBuBrettPhos [ligand]; tAmOH [solvent] | 0.9398 |
| `CGPOC1:a007fb69a18d90cbbd4b8b594c07c01ff3e880db6fbdb6074db5ac7326229de6` | 337 | 130 | KOPiv [base]; (DiMeIHept Cl) Pd(cinnamyl)Cl [catalyst]; DiMeIHept Cl [ligand]; nBuOAc [solvent] | 0.9371 |
| `CGPOC1:2e4b528a7a06183f2d3f03f8a96a86c22974989772e2c57ce095a59235cdcfcd` | 332 | 97 | Imidazole [base]; T3P [coupling_reagent]; MeCN [solvent]; No ligand | 0.9852 |
| `CGPOC1:8343d1acb7e678660e7ea096a8bd6f690b5696b3a06f3cee8a14e9af2f35064c` | 330 | 106 | N-Methyldicyclohexylamin [base]; MorDalphos Pd G3 [catalyst]; MorDalphos [ligand]; BnOH [solvent] | 0.9011 |
| `CGPOC1:c37bf5d55ee8e245d85897707dc8a0bcf81453bbc40a26bed4d292e608039185` | 329 | 119 | KF [base]; XPhos Pd(crotyl)Cl [catalyst]; XPhos [ligand]; TBME [solvent] | 0.9154 |
| `CGPOC1:448542eed63104b78077ea239a052e3bb9b11c9f7dc913811357c2c5f731c88e` | 323 | 189 | NMM [base]; TCFH [coupling_reagent]; DCM [solvent]; No ligand | 0.9694 |
| `CGPOC1:1b5f7177485b3dc7cc7385ff8b287c0bee4947490963ef37b5da3205f026b7f2` | 321 | 63 | K3PO4 [base]; tBuBrettPhos Pd(allyl)OTf [catalyst]; tBuBrettPhos [ligand]; PrOH [solvent]; water [solvent] | 0.9120 |
| `CGPOC1:59d1b291fb2f684315fa4c100bb9c98cd313a897042f7ca86cd13e59bdace15e` | 319 | 116 | KOtPent [base]; Pd-PEPPSI-IPent Cl o-picoline [catalyst]; IPENT Cl [ligand]; Ph2O [solvent] | 0.9235 |
| `CGPOC1:0a0e78ce3784d51dc476a214d3d3743d4ad8a672a3b2cd9e2b8dfc4fd781b1e2` | 310 | 120 | DMAP, NEt3 [base]; CDT [coupling_reagent]; THF [solvent]; No ligand | 0.9830 |
| `CGPOC1:b5b9e8795e6afb250b5d57c7c1cbf52ba5fb75878bd89151f8e7b8912c5161db` | 309 | 66 | Brij 35, Octanoic acid [additive]; P(tBu)3 Pd(crotyl)Cl [catalyst]; P(tBu)3 [ligand]; ethylene glycol [solvent]; THF [solvent] | 0.9464 |
| `CGPOC1:46e57d2ae4523be5c5f9d37d58b4bbb34eff526361a7c786197c165703f208bb` | 299 | 137 | Sodium phenoxide [base]; triisobutylphosphatrane Pd G6 [catalyst]; triisobutylphosphatrane [ligand]; Propionitrile [solvent] | 0.9241 |
| `CGPOC1:8bb03f92779b179a8d862f2687edb60b80175d3c2a87796044bb8cce2300fb66` | 299 | 67 | HOSu [additive]; lutidine [base]; DIC [coupling_reagent]; MeCN [solvent]; No ligand | 0.9826 |
