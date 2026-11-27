# Registry and dataset review, 25 November 2026

The registry corrections apply to new identifier lookups and conversions. The published condition and fragment indexes continue to carry their original build identities. A fresh canonical conversion is required to propagate revised ingredient identities to stored recipes.

| Area | Reproduced defect | Verification |
| --- | --- | --- |
| Material identity | Salts, isotopologues, formulations and metal complexes shared incomplete or incorrect graphs. The sodium methoxide entry lacked the alias `MeONa`. | Corrected 144 stable registry entries. All 27,180 supplied graphs parse. Eight unsupported assertions were withdrawn. |
| Condition recommendations | A failed experiment could supply a positive recipe vote or stop retrieval from searching further. | On HiTEA query `RKB-03582582`, 35 of 50 suggestions had only 0% UV-area yields. The updated engine returned 50 exact matches with zero failed-only recipes. |
| Fragment search | Splitting reaction SMILES on every `>` discarded molecular coordinate bonds. | The old index omitted 1,295 literature records. A complete replay recovered searchable product components for 1,141. |
| Patent routes | Comparing map numbers alone allowed a product atom to acquire another element. | RDKit verified 274 affected steps. The reader rechecked 278 previously qualified route trees and returned them as unresolved with their source memberships. |
| Forward and retro checks | Removing stereo before comparing products could erase a real inversion or an isotope change. | On stored patent and CbyD controls, the corrected check retained the reported product and rejected explicitly inverted and carbon-13-substituted targets. |

The chemistry identity ledger records the primary references and correction reasons. A remaining structural flag requires review of the specific source material. An isolated ion cannot establish the counterion, and component-only representations of metal complexes may omit coordination geometry. [Registry identity evidence](condition_registry_identity_evidence_20261125.jsonl)

The two common structure representations of ordinary DMSO, `CS(C)=O` and `C[S+](C)[O-]`, now resolve to CAS 67-68-5. The registry retains an independent identity for the deuterated solvent, CAS 2206-27-1. A supplied deuterated graph paired with the name of the ordinary solvent still reports a conflict. [PubChem DMSO](https://pubchem.ncbi.nlm.nih.gov/compound/679)

The usage-ranked screen retained 714 possible name/structure discrepancies. Of these, 122 appear as resolved ingredients in 598 of 227,101 default-release recipe representatives. The other residual issue is spreadsheet date conversion in the older registry export. The trisodium citrate dihydrate record now has a verified CAS and structure. The other 431 date-form CAS cells still require identity corroboration.

The default release contains 1,726,556 indexed condition observations. For a recipe variant with supporting precedents, the historical outcome summary can still include its failed experiments. An unreported yield remains unknown. Explicitly ambiguous raw ingredients appear in the recommendation's caution.

The frozen CbyD development partition yielded 1,322 additional review observations across 926 routes. The common ingredients with corrected identities include potassium carbonate, sodium borohydride and lithium aluminium hydride. Together with the 1,684 previously qualified steps, 3,006 observations are available for independent review. In 641 of them, a map-free trial met machine admission criteria. Zero CbyD recipes have been approved for condition transfer. The local dataset folder contains the new records, checksum manifest and the link to the original queue. [Review dataset](../../datasets/processed_datasets/curation/cbyd-v2026.11-condition-review-20261125-r2/README.md)

Older route and operator files can be read with runtime checks on their retained source evidence. The historic operator support totals were aggregated before stricter replay. Recompile from the original canonical observations to update those totals. Both known positive controls also demonstrate a blind-search coverage limitation under the default operator budget. A bounded miss is inconclusive even when targeted replay can reproduce the known transformation.

Of the historically stereo-relaxed CbyD records, 1,821 of 6,207 retro templates and 881 of 1,943 forward operators failed current source replay. Under a paired replay on identical generated outcomes, the older rule accepted explicit stereo conflicts for 1,794 of the retro templates and all 881 forward operators. A single source reaction can contribute several records at different levels of abstraction.
