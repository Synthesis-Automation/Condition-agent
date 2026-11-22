# Standalone precedent support for a route step

Use `condition_recommender.step_plausibility_cli` to assess a proposed precursor-to-product transformation against a local reaction catalog. The command returns source reaction records, structural and functional-group comparisons, reported outcomes, and an ordinal evidence grade. It runs locally, through the standalone chemistry and recommendation packages.

## Run one step

From the repository root, give the reaction in the forward direction, even when you are planning the route retrosynthetically.

```powershell
python -m condition_recommender.step_plausibility_cli 'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1' --match-limit 3
```

The default uses the active processed release. The text output includes literature citations when its evidence catalog supplies them. To pin a release, pass `--release PATH` with a dataset root or release manifest. To supply your own indexes, use the explicitly bound pair. Add `--json` for the versioned Python assessment, including individual source observations, reference IDs, source reaction structures and available recipes.

```powershell
python -m condition_recommender.step_plausibility_cli '<reactants>>product>' --index path/to/generic_index.sqlite --shared-index path/to/generic_index.shared_core.sqlite --json
```

The public Python entry point is `condition_recommender.assess_step_plausibility(reaction_smiles, index, shared_index, *, candidate_limit=128, match_limit=20)`. It returns an immutable assessment with a `to_dict()` method. Use `candidate_limit` to bound the search, and `match_limit` to select the number of evidence rows displayed. All graph-qualified candidates visited contribute to the statistics before the display limit is applied. The result records both truncation flags.

## Read the result

| Rank | Grade | Interpretation |
|---:|---|---|
| — | `unresolved` | The supplied structures do not support the graph comparison needed for a grade. |
| 0 | `no_qualified_precedent` | The bounded lookup found zero qualifying transformations. |
| 1 | `structural_analogue` | A source matches the transformation at a qualified graph level. |
| 2 | `context_supported` | A sufficiently close reaction-core match retains the queried functional-group context on the corresponding reacting partner or partners. |
| 3 | `measured_support` | A qualifying source group reports an outcome above 10% on an exact or sufficiently similar, context-covered analogue. |
| 4 | `repeated_measured_support` | At least three qualifying source groups supply such measurements. |

The levels describe the strength of evidence that the intended transformation has worked in related chemistry. A score of 4 leaves open the yield of a proposed recipe, selectivity on new substrates, and experimental feasibility on a particular scale. The >10% yield threshold and three-group rule are explicit development choices. A source above 10% can still be a poor preparative result.

Read the grade together with `measured_evidence_status` and the counts of above-threshold, low, zero, unreported and unusable yields among comparable observations. Repeated positive reports can coexist with many low-yield experiments. The grade captures precedent for achieving the transformation under some source recipes. The distribution and the cited source conditions reveal where that evidence is mixed.

The lookup uses observed edits and a verified reaction core, the same qualification used by the condition recommender. Product-only hits involving different edits cannot raise the grade. A detailed local analogue needs a minimum reacting-partner Morgan fingerprint similarity of 0.50 to count toward measured support. This cutoff applies individually to each aligned reacting partner. The source substrate must retain the required context on the corresponding participant. A free acid and a phenol on unrelated source components cannot establish tolerance of an ortho-hydroxybenzoic acid on the aryl substrate. For queries with specified stereochemistry, v1 requires an exact whole-reaction identity for measured support. The functional-group vocabulary and explicit positional relationships limit what v1 can assess. Inspect reported differences, unreported outcomes and measured low yields alongside positive examples.

Measured support from a nonexact analogue also requires coverage in the reverse direction. Extra source features can make an analogy misleading, such as an electron-withdrawing cyano group on the source aryl chloride that is absent from the proposed substrate. The text and JSON list extra source features and unresolved differences. Rank 2 requires only that the source covers the query's tracked context.

Multiple observations from one reference contribute at most one positive source group. The grouping also merges detected copies of the same experiment across references. Separate groups are evidence of depth within the catalog. Independent laboratories, unreported failed attempts, publication bias, equivalent experiments encoded differently, and transfer of an exact operating procedure require additional review.

Read positive source groups alongside the distinct positive source reaction count. Reports of one successful substrate pair can establish repeat support under reported conditions while leaving broad substrate tolerance uncertain.

## Development examples

The real-data examples use a pinned processed release. The Suzuki coupling to biphenyl illustrates a transformation with exact source records. The salicylate borylation asks how far the many core matches transfer to a free acid adjacent to a phenol. An analogous benzonitrile borylation changes the aryl context. A constructed ethoxy homologue of a reported methoxypyridyl sulfonamide tests transfer from close analogues. Omitting the salicylate query's boron reagent illustrates a missing atom source. These examples exercise the implementation. Independent chemist review and held-out measured experiments are still needed to evaluate predictive value.

Release `05fd3b97e436c3ba79147e26` produced these results with a 128-candidate budget. Each structurally resolved query reached the search cap.

| Proposed reaction | Rank | Evidence among visited records |
|---|---:|---|
| Bromobenzene + phenylboronic acid → biphenyl | 4 | 73 exact records, 69 with usable yield above 10%, in 65 positive source groups. |
| Free 4-bromosalicylic acid → aryl Bpin | 1 | 128 detailed local matches from 113 reference IDs. None covered the free salicylic-acid context on the reacting aryl partner. |
| 4-Bromobenzonitrile → aryl Bpin | 2 | One of 128 detailed local matches covered the nitrile and its tracked positional context. Its weakest-partner similarity was 0.275. |
| Ethoxy-substituted bromopyridyl sulfonamide → aryl Bpin | 4 | Zero exact hits. Four observations all document one closely related methoxy substrate pair. The source yields were 60%, 75%, 75% and 77%, contributing three positive source groups and one distinct positive reaction. |
| Same salicylate target, diboron reagent omitted | — | The product graph could not be accounted for from the supplied precursor. |

Replay the borylation inputs against a local release using PowerShell.

```powershell
$diboron = 'CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C'
$salicylateProduct = 'O=C(O)c1ccc(B2OC(C)(C)C(C)(C)O2)cc1O'
python -m condition_recommender.step_plausibility_cli "O=C(O)c1ccc(Br)cc1O.$diboron>>$salicylateProduct"
python -m condition_recommender.step_plausibility_cli "N#Cc1ccc(Br)cc1.$diboron>>N#Cc1ccc(B2OC(C)(C)C(C)(C)O2)cc1"
python -m condition_recommender.step_plausibility_cli "CCOc1ncc(Br)cc1NS(=O)(=O)c1ccc(F)cc1.$diboron>>CCOc1ncc(B2OC(C)(C)C(C)(C)O2)cc1NS(=O)(=O)c1ccc(F)cc1"
python -m condition_recommender.step_plausibility_cli "O=C(O)c1ccc(Br)cc1O>>$salicylateProduct"
```
