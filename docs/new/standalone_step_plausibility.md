# Standalone precedent support for a route step

Use `condition_recommender.step_plausibility_cli` to assess a proposed precursor-to-product transformation against a local reaction catalog. The command returns source reaction records, structural and functional-group comparisons, reported outcomes, and an ordinal evidence grade. It runs locally, through the standalone chemistry and recommendation packages.

## Run one step

From the repository root, give the reaction in the forward direction, even when you are planning the route retrosynthetically.

```powershell
python -m condition_recommender.step_plausibility_cli 'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1' --match-limit 3
```

For a short terminal report, add `--concise`. The first line gives the step score out of four. The next line gives the evidence counts, followed by one line per displayed precedent with the graph match, source yield and citation. Concise mode shows three precedents by default. Set `--match-limit` to display a different number. With `--recipe`, the short report shows which measurements match the proposed conditions.

```powershell
$rxn = 'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1'
python -m condition_recommender.step_plausibility_cli $rxn --concise
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
| 3 | `transferable_precedent` | One or two source groups document an exact reaction or a structurally transferable analogue that qualifies under the yield rule below. |
| 4 | `repeated_transferable_precedent` | At least three source groups document qualifying exact reactions or transferable analogues. |

**Yield affects the score only when a comparable source has a usable numeric target yield below 10%.** A cited exact precedent with an unreported yield, 10% yield or 90% yield contributes equally. A 0% or 9.9% observation is adverse evidence for its recorded recipe. It remains visible and cannot supply rank 3 or 4 support. A separate qualifying experiment can still support the step. A step supported only by low-yield observations can reach rank 2 on its structural evidence. Qualitative or unusable measurements carry no assumed percentage.

The standalone grade reads the yield value and its usability flag from the bound condition index. Older index records can lack assay details or a link to the measured product. Use the source record or the recipe evidence view to check attribution when a low yield affects an important decision.

For an unresolved query, the CLI reports any product atom counts that exceed the
submitted reactant inventory. JSON exposes the element and count pairs in
`query_context.product_element_excess`. Check the structures, donor assignment
and stoichiometry when a deficit is reported. A balanced inventory can also
remain unresolved if the atom correspondence is ambiguous.

The grade describes the structural support for a transformation in the records visited. Check `measured_evidence_status` for the availability of numeric outcomes. A step can earn 4/4 even if every supporting source omits yield. Its expected yield, selectivity and performance under a particular recipe remain separate experimental questions.

For mixed evidence, inspect the counts of low nonzero, zero, unreported, unusable and at-least-10% yields, together with the cited recipes. Each sub-10% observation is excluded from the support count. The grade can still reflect other experiments on the same chemistry.

The lookup qualifies the observed bond edits and reaction core. For an analogue to count toward rank 3 or 4, each aligned reacting partner must have a Morgan fingerprint similarity of at least 0.50 and cover the tracked functional groups on the corresponding partner. For example, evidence for an ortho-hydroxybenzoic acid requires the acid and phenol on the same aryl substrate in the appropriate positions. A query with specified stereochemistry requires an exact whole-reaction match for ranks 3 and 4.

A transferable analogue also requires functional-group coverage in the reverse direction. For example, a cyano group present only on the source aryl chloride can make the precedent a poor model for the proposed substrate. The text and JSON show extra source features and unresolved differences. Rank 2 requires only that the source covers the query's tracked context.

Multiple observations from a reference contribute at most one supporting source group. Grouping also merges probable copies across references using source observation identity or matching reaction graphs and a known recipe ID. The magnitude or availability of a yield does not distinguish groups. The count is a proxy for depth of documentation. Claims about independent replication require inspection of the underlying publications.

Read `support_counts.supporting_source_group_count` alongside `support_counts.supporting_distinct_reaction_count`. Multiple groups can document one substrate pair while leaving tolerance of other substrates uncertain. The JSON assessment identifies this scoring contract as `step_plausibility.v2`.

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

Policy 2.0 produced these results on release `05fd3b97e436c3ba79147e26` with a 128-candidate budget. Each structurally resolved query reached the search cap.

| Proposed reaction | Rank | Evidence among visited records |
|---|---:|---|
| Saved route step T01-R1-S01, preparation of a Weinreb amide | 3 | The exact precursor and product pair is indexed to US 9,029,536 B2. Its yield is unreported. One source group supplies transferable precedent. |
| Bromobenzene + phenylboronic acid → biphenyl | 4 | 73 exact records across 67 supporting groups. Of these records, 69 report a yield of at least 10% and four omit yield. |
| Free 4-bromosalicylic acid → aryl Bpin | 1 | 128 detailed local matches from 113 reference IDs. None covered the free salicylic-acid context on the reacting aryl partner. |
| 4-Bromobenzonitrile → aryl Bpin | 2 | One of 128 detailed local matches covered the nitrile and its tracked positional context. Its weakest-partner similarity was 0.275. |
| Ethoxy-substituted bromopyridyl sulfonamide → aryl Bpin | 3 | Four observations document the same closely related methoxy substrate pair across two supporting groups. The different reported yields do not establish additional groups. |
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
