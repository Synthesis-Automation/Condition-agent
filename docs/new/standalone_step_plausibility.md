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

## Compare a proposed recipe with measured controls

Add `--recipe` to see what happened on the same substrates or close analogues under the specified conditions. For a reproducible example, build the small experiment catalog from the six published observations in [Liu et al., Supporting Information, Table S1](https://doi.org/10.1021/jacs.7b07390.s001). Run these commands from the repository root with an active processed release.

```powershell
$rxn = 'CC(=O)c1ccc(Br)cc1.Cc1ccc(S)cc1>>CC(=O)c1ccc(Sc2ccc(C)cc2)cc1'
$catalog = 'results/ai_native/recipe_step_evidence_20261123/liu.sqlite'
python -m condition_recommender.experiment_catalog examples/common/liu_thioether_controls.jsonl $catalog
python -m condition_recommender.step_plausibility_cli $rxn --recipe examples/common/liu_thioether_standard_recipe.json --experiment-catalog $catalog --match-limit 6
python -m condition_recommender.step_plausibility_cli $rxn --recipe examples/common/liu_thioether_no_base_recipe.json --experiment-catalog $catalog --match-limit 6
```

The source reports 97% isolated yield under the standard conditions, 0% for the no-base experiment, and trace product in the dark. Each refers to a distinct experiment for the same target product. The analysis preserves the original measurement and its source. It also keeps the route-step grade so you can distinguish precedent for making a bond from evidence for making it under a particular recipe.

The [standard recipe](../../examples/common/liu_thioether_standard_recipe.json) specifies component identities, base equivalents, solvent, illumination, atmosphere, temperature and time. The [no-base recipe](../../examples/common/liu_thioether_no_base_recipe.json) explicitly declares the base absent. Leaving a role blank leaves its presence open. A result of `matches_stated_fields` means all the conditions you specified agree. Read the additional conditions reported by the source before transferring a procedure. An observed difference, an incompletely reported procedure and a procedure with no comparable fields each have their own comparison status.

Where the registry identifies an ingredient but cannot establish its role, the comparison can still confirm the ingredient and its amount. A query requiring that role can remain a partial match. The JSON field `reported_conditions` retains the original structured source information when available. Check it for details outside the comparison, such as pressure and reactor setup.

The output shows each experiment's source structure, weakest reacting-partner similarity, functional-group differences, recipe agreement and observed outcome. For an analogue, the yield belongs to the displayed source substrates. Counts are grouped by recipe agreement and include all structurally comparable experiments retrieved, even when `--match-limit` shortens the display. The step grade draws on the bound recommendation index. Supplemental measurements appear in the recipe view. Estimating a success probability requires an independently evaluated set of attempted reactions.

Use `--json` for the full evidence and source locators. The active release supplies its admitted precedents and evidence catalog. The older v2 evidence catalog can additionally recover exact reactions through its release's fragment index. Newly built v3 evidence catalogs also index local-core analogues, including verified observations excluded from condition recommendation. Use `--experiment-catalog` to add local studies independently of a release rebuild. Its input is a JSON array or JSONL of canonical `source_observation.v1` records, such as the [six transcribed controls](../../examples/common/liu_thioether_controls.jsonl).

In Python, call `condition_recommender.assess_step_recipe_evidence(reaction_smiles, recipe, index, shared_index, *, catalog=None, fragment_index=None, experiment_catalog=None, candidate_limit=128, match_limit=20)`. The returned `StepRecipeEvidence` includes the original `step_assessment`, the normalized query recipe, experiment-level evidence, counts and the scope of each search. Each index receives its own candidate budget.

## Read the result

| Rank | Grade | Interpretation |
|---:|---|---|
| — | `unresolved` | The supplied structures do not support the graph comparison needed for a grade. |
| 0 | `no_qualified_precedent` | The bounded lookup found zero qualifying transformations. |
| 1 | `structural_analogue` | A source matches the transformation at a qualified graph level. |
| 2 | `context_supported` | A sufficiently close reaction-core match retains the queried functional-group context on the corresponding reacting partner or partners. |
| 3 | `measured_support` | A qualifying source group reports an outcome above 10% on an exact or sufficiently similar, context-covered analogue. |
| 4 | `repeated_measured_support` | At least three qualifying source groups supply such measurements. |

For an unresolved query, the CLI reports any product atom counts that exceed the
submitted reactant inventory. JSON exposes the element and count pairs in
`query_context.product_element_excess`. Check the structures, donor assignment
and stoichiometry when a deficit is reported. A balanced inventory can also
remain unresolved if the atom correspondence is ambiguous.

The levels describe the strength of evidence that the intended transformation has worked in related chemistry. A score of 4 leaves open the yield of a proposed recipe, selectivity on new substrates, and experimental feasibility on a particular scale. The >10% yield threshold and three-group rule are explicit development choices. A source above 10% can still be a poor preparative result.

Read the grade together with `measured_evidence_status` and the counts of above-threshold, low, zero, unreported and unusable yields among comparable observations. Repeated positive reports can coexist with many low-yield experiments. The grade captures precedent for achieving the transformation under some source recipes. The distribution and the cited source conditions reveal where that evidence is mixed.

The lookup uses observed edits and a verified reaction core, the same qualification used by the condition recommender. Product-only hits involving different edits cannot raise the grade. A detailed local analogue needs a minimum reacting-partner Morgan fingerprint similarity of 0.50 to count toward measured support. This cutoff applies individually to each aligned reacting partner. The source substrate must retain the required context on the corresponding participant. A free acid and a phenol on unrelated source components cannot establish tolerance of an ortho-hydroxybenzoic acid on the aryl substrate. For queries with specified stereochemistry, v1 requires an exact whole-reaction identity for measured support. The functional-group vocabulary and explicit positional relationships limit what v1 can assess. Inspect reported differences, unreported outcomes and measured low yields alongside positive examples.

Measured support from a nonexact analogue also requires coverage in the reverse direction. Extra source features can make an analogy misleading, such as an electron-withdrawing cyano group on the source aryl chloride that is absent from the proposed substrate. The text and JSON list extra source features and unresolved differences. Rank 2 requires only that the source covers the query's tracked context.

Multiple observations from one reference contribute at most one positive source group. The grouping also merges detected copies of the same experiment across references. Separate groups are evidence of depth within the catalog. Independent laboratories, unreported failed attempts, publication bias, equivalent experiments encoded differently, and transfer of an exact operating procedure require additional review.

Read positive source groups alongside the distinct positive source reaction count. Reports of one successful substrate pair can establish repeat support under reported conditions while leaving broad substrate tolerance uncertain.

## Development examples

The [matched Organic Syntheses controls](../../examples/common/org_syn_1_matched_controls.json)
contain seven literature-backed substrate–product pairs, four donor omissions
and seven alternative products. The source citation and row number identify each
parent. A different product is a challenge to selectivity under the published
conditions. Its experimental yield is unknown. To inspect the α- versus
β-arylation pair with the active local release, run from the repository root.

```powershell
$controls = Get-Content examples/common/org_syn_1_matched_controls.json -Raw -Encoding utf8 | ConvertFrom-Json
$controls.cases | Where-Object { $_.case_id -in 'ketone_alpha_arylation', 'beta_arylation' } | ForEach-Object {
    Write-Output $_.case_id
    python -m condition_recommender.step_plausibility_cli $_.reaction_smiles --match-limit 1
}
```

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
