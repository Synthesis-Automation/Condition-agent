# Core retrosynthesis

`core_retrosynthesis` is the canonical chemistry-first retrosynthesis package.
Its primary path compiles type-agnostic graph operators from structure-derived
reaction observations, admits them by source round trip, and forward-validates
every proposed single-step disconnection before ranking.

Shared chemistry and atom-mapping utilities live at the package root. Historical
RDChiral C-X behavior is isolated under `baselines/cx_rdchiral` and is used only
for controlled evaluation; it is not a second application or routing path.

The reaction core remains observed evidence rather than an executable prediction.
This package composes `reactive_taxonomy`, `condition_registry`, and optional
`condition_recommender` evidence without changing their ownership boundaries.

For the current architecture and implementation matrix, see
[`../docs/new/type_agnostic_retrosynthesis_design_and_status.md`](../docs/new/type_agnostic_retrosynthesis_design_and_status.md).

## Curate an external route corpus

The released higher-level retrosynthesis route JSONL can be converted into a
small, deterministic POC corpus without making its unordered reaction arrays
part of our public route contract:

```powershell
python -m core_retrosynthesis curate-route-corpus `
  datasets/external/higher_level_retrosynthesis/figshare_v2/extracted/datasets/routes/uspto.higher-level.routes.jsonl.gz `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.v1.jsonl.gz `
  --testset datasets/external/higher_level_retrosynthesis/figshare_v2/extracted/datasets/route_testset/uspto190_canon_reactions.smi `
  --maximum-routes 5000
```

The curator reconstructs retrosynthetic step order from canonical molecular
identity, validates complete atom mapping and route connectivity, retains the
full reagent field, excludes routes whose targets occur as products in the
released USPTO-190 reaction list, and assigns splits by patent. The default POC
uses 3-6-step linear routes with at least one algorithmic higher-level step
reduction and selects at most one route per patent. These are sampling
constraints for a clean first experiment, not production rejection rules for
branched chemistry.

See
[`../docs/new/higher_level_route_dataset_audit_and_poc.md`](../docs/new/higher_level_route_dataset_audit_and_poc.md)
for the source audit, limitations, and required follow-up cleanup.

Render a reproducible random sample for chemistry review:

```powershell
python -m core_retrosynthesis render-route-review `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.v1.jsonl.gz `
  results/core_retrosynthesis/route_reviews/higher_level_routes_random50.html `
  --sample-size 50 --seed 20260814
```

The self-contained report opens with a chemist-oriented synthetic sequence:
starting materials, intermediates, added reactants and conditions, then final
product. A toolbar control reverses it into retrosynthetic reading order. The
original and higher-level mapped reaction drawings remain below for auditability.
Review status and notes persist in browser storage and can be exported as JSON.

For a smaller report containing only the sequence and review controls, add
`--sequence-only`.

Convert the curated observations to the shared evidence-neutral route-tree
schema:

```powershell
python -m core_retrosynthesis convert-route-trees `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.v1.jsonl.gz `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.tree.v2.jsonl.gz
```

`ReactionRouteTree` alternates molecule occurrences and reaction occurrences,
supports convergent branches, and distinguishes observed, inferred, and
predicted evidence. Existing planned routes now serialize through the same v2
contract while retaining temporary compatibility fields. See
[`../docs/new/route_tree_contract_and_conversion.md`](../docs/new/route_tree_contract_and_conversion.md).

Build context-preserving route-core projections from those trees:

```powershell
python -m core_retrosynthesis build-route-cores `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.tree.v2.jsonl.gz `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.core.v1.jsonl.gz `
  --workers 8
```

Render the minimized chemistry and cross-step atom lineage:

```powershell
python -m core_retrosynthesis render-route-core-review `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.random50.core.v1.jsonl.gz `
  results/core_retrosynthesis/route_reviews/higher_level_routes_random50_route_core.html `
  --sample-size 50 --seed 20260814
```

Route cores contain complete generic step signatures and cores, display-only
minimum reactions, ambiguity-preserving cross-step atom lineage, route-core
identities, and explicit two-/three-event motif definitions. They are projection
and learning evidence, not executable multistep templates. See
[`../docs/new/route_core_projection_design_and_status.md`](../docs/new/route_core_projection_design_and_status.md).

Replay observed route actions through the validated single-step operator ladder:

```powershell
python -m core_retrosynthesis evaluate-route-actions `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.random50.tree.v2.jsonl.gz `
  results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz `
  results/core_retrosynthesis/route_action_evaluation/routes.poc.random50.action_labels.v3.jsonl.gz `
  --labels-only --workers 8 --overwrite
```

The benchmark retains every step and reports independent site, retained-edit,
synthon, exact-precursor, strategy, realization, and executable-operator evidence.
Candidate replay is optional; returned alternatives are validated but remain
unchosen alternatives rather than asserted negatives. See
[`../docs/new/route_action_replay_benchmark_design_and_status.md`](../docs/new/route_action_replay_benchmark_design_and_status.md).

For long label or replay jobs, partition routes by stable identity and resume
only manifest-validated shards:

```powershell
python -m core_retrosynthesis evaluate-route-actions SOURCE LIBRARY PART `
  --labels-only --workers 8 --shard-count 16 --shard-index 0 --resume

python -m core_retrosynthesis merge-route-action-shards MERGED PARTS...
```

The strict merger requires a complete, compatible shard set, verifies each
checksum and route assignment, rejects duplicates, and writes deterministic
route order. The complete 5,000-route labels audit contains 18,647 steps:
18,280 retained-edit labels, 18,078 synthon labels, 18,307 exact-precursor
labels, 13,099 STRAT1 labels, and 1,201 independently executable operators.

Render a stratified review of labels promoted beyond strict operator admission:

```powershell
python -m core_retrosynthesis render-route-action-review `
  results/core_retrosynthesis/route_action_evaluation/routes.poc.full5000.action_labels.v3.jsonl.gz `
  results/core_retrosynthesis/route_reviews/higher_level_routes_promoted_action_labels_v3.html `
  --sample-size 120
```

Review decisions and notes persist in browser storage and export as JSON. A
chemist accepting a route label does not admit an executable template.

## Representation

Two initial abstraction levels are generated:

- `L1`: edit atoms plus connected precursor-only handle subgraphs;
- `L2`: L1 chemistry plus the first molecular shell.

The SMARTS contain structural requirements only. Reaction-core substituent
profiles, unchanged spectator groups, steric accessibility, and electronic
activation are retained as separate precedent context. Context can be used for
ranking or disabled to isolate candidate-generation performance.

Only single-heavy-atom-event, single-product C-N, C-O, and C-S bond formations
are enabled in this first experiment. A directly connected N-H, O-H, or S-H
loss may be normalized into that event. The compiler writes explicit mapped
hydrogen-count queries into the SMARTS so reverse application restores a
valence-correct X-H precursor, including adjacent aromatic tautomer carriers.
Every admitted SMARTS must reproduce its recorded contributing precursors when
applied to its source product.

## Full and Compact input libraries

Commands that consume converted reactions accept the recommendation-library
root and a mode selector:

```powershell
--library-mode full     # complete converted corpus; the default
--library-mode compact  # deterministic sampled corpus for development
```

For example, `datasets/literature --library-mode compact` resolves to
`datasets/literature/compact`. Row-oriented commands read that mode's
deduplicated `combined_records.jsonl.gz`. Shard-oriented benchmarks and the
resumable full-scale builder read only the completed batch shards named by
`combined_batch_manifest.json`; incomplete or saved-but-not-combined batches
are excluded. An explicit JSONL file, legacy shard directory, or already
selected `full/` or `compact/` directory remains valid and does not require a
path migration.

Use separate output directories for Full and Compact builds. Compact is the
recommended iteration mode; use Full for final operator-library construction
and release evaluation.

## Build a balanced library

Use a separate limit for every included cohort so C-N rows do not consume a
global limit before C-O and C-S are read:

```powershell
python -m core_retrosynthesis build-library `
  datasets/literature `
  results/core_retrosynthesis/compact/core_templates.json.gz `
  --library-mode compact `
  --include "C_N_Coupling*.jsonl.gz" `
  --include "C_O_Coupling*.jsonl.gz" `
  --include "C_S_Coupling*.jsonl.gz" `
  --max-rows-per-include 1000
```

Use `--level L1` or `--level L2` to build only one abstraction level.

## Disconnect a target

```powershell
python -m core_retrosynthesis disconnect `
  results/core_retrosynthesis/core_templates.json.gz `
  "CCNc1ccc2ccccc2c1" `
  --concise --top-k 5
```

This searches all supported bonds and both abstraction levels. Add `--bond C-N`
or `--level L2` to restrict the search. Add `--no-context` to use the same
chemistry/similarity score without the explicit context contribution.
Candidates are similarity-ranked before proposal-level featurization; by
default, only the best 50 template/precursor pairs are forward-validated. Use
`--max-candidates-to-validate` to change that bounded shortlist.

## Recommended ensemble

The completed comparison supports retaining the previous RDChiral ranking and
using contextual L1 core rules only to fill unused ranks:

```powershell
python -m core_retrosynthesis disconnect-ensemble `
  results/core_retrosynthesis_comparison/balanced_500_final/baseline_templates.json.gz `
  results/core_retrosynthesis_comparison/balanced_500_final/core_templates.json.gz `
  "CCNc1ccc2ccccc2c1" `
  --concise --top-k 10
```

Full JSON output records whether each proposal came from `rdchiral_baseline` or
`core_l1_context`.

## Paired held-out comparison

The comparison command creates a deterministic reference-group split, builds
both libraries from training references only, and evaluates baseline,
core-neutral, and core-context ranking on held-out products:

```powershell
python -m core_retrosynthesis compare `
  datasets/literature `
  results/core_retrosynthesis_comparison/compact `
  --library-mode compact `
  --include "C_N_Coupling*.jsonl.gz" `
  --include "C_O_Coupling*.jsonl.gz" `
  --include "C_S_Coupling*.jsonl.gz" `
  --max-rows-per-include 1000 `
  --test-fraction 0.2 --top-k 10 `
  --max-test-targets 100 `
  --max-candidates-to-validate 20
```

Omit `--max-test-targets` for the complete held-out evaluation. A bound is
useful during development because forward validation and contextual scoring
featurize generated reactions.

Outputs are:

- `baseline_templates.json.gz`;
- `core_templates.json.gz`; and
- `comparison.json` with extraction, coverage, validity, exact-precursor, and
  reaction-center recovery metrics.

Exact precursor recovery is intentionally reported separately from correct
reaction-center recovery because an alternative precursor can be chemically
reasonable.

## Completed comparison

The final deterministic experiment used 500 rows from each C-N, C-O, and C-S
cohort: 1,500 source rows, 1,400 training rows, and all 86 eligible held-out
targets from 62 held-out reference groups. The held-out set contained 17 C-N,
44 C-O, and 25 C-S targets. The evaluator forward-validated a maximum of 20
template/precursor pairs per method and rejected invalid or structurally
unresolved proposals before ranking.

| Method | Top-1 exact | Top-5 exact | Top-10 exact | Valid | Candidates/target |
| --- | ---: | ---: | ---: | ---: | ---: |
| RDChiral baseline | 0.372 | 0.744 | 0.884 | 1.000 | 7.35 |
| Core L1 + context | 0.291 | 0.721 | 0.942 | 1.000 | 8.56 |
| Core L2 + context | 0.337 | 0.709 | 0.802 | 1.000 | 5.74 |
| Core L1/L2 + context | 0.337 | 0.756 | 0.884 | 1.000 | 6.66 |
| Baseline-first L1 ensemble | **0.372** | **0.802** | **0.965** | **1.000** | 8.74 |

Both extraction methods admitted the same 400 training observations after the
core compiler was corrected to retain complete connected precursor-only handle
subgraphs. The baseline produced 76 templates. Core compilation produced 49 L1
operators and 84 L2 operators; the executable library contains multiple SMARTS
variants per operator because distinct precursor handles remain separate.

The result does not support replacing the baseline. L1 is broader and has lower
top-1 precision, but it recovers C-N disconnections absent from the baseline.
The baseline-first ensemble preserved every baseline result and rescued seven
additional held-out C-N reactions: five within the top five and two at rank six.
It introduced no top-10 losses. This is therefore the recommended POC workflow.

Hydrogen-event normalization also removed the zero-candidate failure for
`N#Cc1ccccc1-n1cccn1`. All methods now propose the independently supported,
valence-correct alternative `N#Cc1ccccc1Br.c1cn[nH]c1`. It is not an exact
match to the recorded protected precursor, so it remains uncredited by the
strict exact-precursor and full-transition-center metrics.

The canonical local report is
`results/core_retrosynthesis_comparison/balanced_500_hydrogen_normalized/comparison.json`.
It includes aggregate and per-bond metrics plus every held-out target, expected
precursor, candidate list, exact-precursor rank, and center rank for all method
ablations. The complete offline comparison took about 16.5 minutes; normal
ensemble target search executes only baseline and contextual L1 retrieval.

## Render an HTML chemistry review

Generate a self-contained report with inline SVG drawings for target molecules,
recorded reactions, and generated reactions:

```powershell
python -m core_retrosynthesis render-report `
  results/core_retrosynthesis_comparison/balanced_500_final_filtered/comparison.json `
  results/core_retrosynthesis_comparison/balanced_500_final_filtered/review.html `
  --top-k 5
```

The default report compares the baseline, contextual L1 core method, and the
recommended ensemble. Repeat `--method` to select other ablations, use
`--max-targets` for a smaller report, or increase `--top-k` to render more
proposals. The HTML has no external assets and can be opened directly in a
browser.

## Structurally diverse benchmark

The generic extension routes by normalized graph edits rather than reaction
names. Named dataset cohorts are used only to construct a balanced sample. The
supported structural archetypes are acyl substitution, C-C coupling, carbonyl
oxidation and reduction, conjugate addition, carbonyl condensation/reductive
amination, and azide-alkyne ring formation.

Every RDChiral or core-derived template must regenerate its source precursors.
Generated proposals are remapped, featurized, and rejected unless their forward
graph transformation is classified as the same archetype as the template.
Stereochemistry-relaxed source round-trip is retained explicitly when a reverse
edit creates a stereocenter absent from the product.

Run the balanced comparison with round-robin sampling across all shards:

```powershell
python -m core_retrosynthesis compare-diverse `
  datasets/literature `
  results/diverse_retrosynthesis_poc/compact/balanced_round_robin_100 `
  --library-mode compact `
  --max-rows-per-cohort 100 `
  --max-targets-per-transformation 5 `
  --top-k 10 --max-candidates-to-validate 20
```

The completed POC used 700 source rows, 536 training rows, and 27 eligible
held-out targets from 60 held-out reference groups. The small target count,
especially one carbonyl-condensation and one ring-formation target, makes this
an architecture and coverage test rather than a production accuracy estimate.

| Method | Coverage | Top-1 exact | Top-10 exact | Top-10 center | Valid |
| --- | ---: | ---: | ---: | ---: | ---: |
| Generic RDChiral | 0.519 | 0.407 | 0.444 | 0.444 | 1.000 |
| Generic core L1/L2 + context | **1.000** | **0.741** | 0.778 | 0.778 | 1.000 |
| Baseline-first generic-core ensemble | **1.000** | **0.741** | **0.815** | **0.815** | 1.000 |

Core exact top-10 recovery by archetype was 0.80 for acyl substitution, 0.40
for C-C coupling, 1.00 for carbonyl condensation, 1.00 for oxidation, 1.00 for
reduction, 0.60 for conjugate addition, and 1.00 for ring formation. C-C
coupling and conjugate addition are the useful stress cases: they generate
multiple structurally valid precursor alternatives and expose ranking and
chemoselectivity limitations hidden by C-X tests.

Render the visual review with:

```powershell
python -m core_retrosynthesis render-report `
  results/diverse_retrosynthesis_poc/compact/balanced_round_robin_100/comparison.json `
  results/diverse_retrosynthesis_poc/compact/balanced_round_robin_100/review.html `
  --method baseline --method core_context `
  --method ensemble_baseline_core_context --top-k 5
```

Apply the generic core library directly with, for example:

```powershell
python -m core_retrosynthesis disconnect-generic `
  results/diverse_retrosynthesis_poc/compact/balanced_round_robin_100/core_templates.json.gz `
  "OCc1ccc(F)cc1" `
  --transformation carbonyl_reduction --concise --top-k 5
```

## C-C coupling and conjugate-addition stress test

The focused stress benchmark separates exact recorded precursor recovery from
product-side disconnection-site recovery. The `SITE1` identity uses canonical
product symmetry classes and the formed/order-changed bonds, so changing a
Suzuki halide or boronate form does not turn the same product cut into a false
center miss.

Run the larger reference-group-held-out comparison with:

```powershell
python -m core_retrosynthesis compare-stress `
  datasets/literature `
  results/diverse_retrosynthesis_poc/compact/stress_cc_michael_250 `
  --library-mode compact `
  --max-rows-per-cohort 250 --test-fraction 0.30 `
  --max-targets-per-transformation 30 `
  --top-k 10 --max-candidates-to-validate 75
```

The completed run used 500 source rows, 355 training rows, 58 held-out
reference groups, and 60 eligible targets split equally between Suzuki C-C
coupling and conjugate addition. Forty-eight targets had more than one
structurally valid candidate product site.

| Method | Coverage | Exact @1 | Exact @10 | Site @1 | Site @5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Generic RDChiral | 0.433 | 0.117 | 0.317 | 0.400 | 0.417 |
| Generic core + context | **1.000** | 0.567 | 0.883 | 0.917 | **1.000** |
| Baseline-first core ensemble | **1.000** | **0.600** | **0.933** | **0.950** | 0.983 |

For C-C coupling alone, core exact @1 was 0.233 while site @1 was 0.833 and
site @5 was 1.000. Much of the apparent exact failure is therefore ranking of
halide/boronate form rather than selection of the wrong bond. For conjugate
addition, core exact @1 was 0.900 and site @1 was 1.000.

The experimental `--diversify-sites` search option interleaves precursor forms
across product sites. It is retained as an ablation, but it is not the default:
on this benchmark it did not improve site recall and reduced exact @10 from
0.883 to 0.867. The ordinary core ranker already recovered every expected site
by rank 5.

Render the focused chemistry review with:

```powershell
python -m core_retrosynthesis render-report `
  results/diverse_retrosynthesis_poc/compact/stress_cc_michael_250/comparison.json `
  results/diverse_retrosynthesis_poc/compact/stress_cc_michael_250/review.html `
  --method baseline --method core_context `
  --method core_site_diverse `
  --method ensemble_baseline_core_context --top-k 5 `
  --title "C-C coupling and conjugate-addition stress review"
```

## Coverage-first data-derived graph operators

The unrestricted operator path removes the seven-archetype admission gate.
Named transformations remain optional annotations; operator admission and
identity use only mapped graph edits, product topology, precursor handles,
hydrogen events, local environments, and source round-trip evidence.

The hierarchy separates:

- `OP1`: handle-independent normalized graph edits;
- `REAL1`: one executable precursor-handle realization of an operator;
- `SYN1`: mapped precursor skeletons after precursor-only handles are removed;
- `SITE1`: edited atoms and bonds on a particular target product.

Every admitted realization must regenerate its recorded precursors. Every
generated candidate is remapped, forward-featurized, and rejected unless its
observed edit signature agrees with the source operator. L0 may be absent for
an individual realization when a sufficiently general template cannot pass
the source round-trip check.

Run the all-shard census and held-out comparison with:

```powershell
python -m core_retrosynthesis compare-operators `
  datasets/literature `
  results/operator_retrosynthesis_poc/compact/all_shards_400 `
  --library-mode compact `
  --max-rows 400 --test-fraction 0.25 `
  --max-targets 24 --max-targets-per-operator 5 `
  --top-k 25 --max-templates 500 `
  --max-candidates-to-validate 100
```

The completed census sampled every one of the 399 shards before taking a
second row. Of 298 training rows, unrestricted compilation admitted 178
(59.7%) into 58 operators and 177 realizations. The seven-archetype view
admitted 78 rows (26.2%). One hundred admitted observations had no supported
named annotation.

The broad 54-target audit increased candidate coverage from 0.704 to 0.963 and
correct-site @25 from 0.315 to 0.630. A balanced 24-target policy audit was then
used to compare the specificity-preserving fallback:

| Method | Coverage | Exact @25 | Synthon @25 | Operator @25 | Site @25 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Seven-archetype L1/L2 | 0.708 | 0.292 | 0.417 | 0.417 | 0.458 |
| Data-derived L2 | 0.917 | 0.625 | 0.708 | 0.792 | 0.708 |
| Data-derived L2 then L1 | **1.000** | **0.625** | **0.875** | **0.875** | **0.875** |
| Data-derived L2/L1/L0 | **1.000** | **0.625** | **0.875** | **0.875** | **0.875** |

For the 12 unannotated targets in that audit, site/operator/synthon @25 reached
0.833 versus 0.167 site recall for the seven-archetype view. L0 produced no
additional recovery in this sample, so it remains a final fill tier rather
than competing with L2 or L1 candidates.

Render the equivalence-aware review with:

```powershell
python -m core_retrosynthesis render-report `
  results/operator_retrosynthesis_poc/compact/all_shards_400/specificity_comparison.json `
  results/operator_retrosynthesis_poc/compact/all_shards_400/review.html `
  --method supported_l1_l2 --method data_l2 `
  --method data_l2_l1 --method data_ladder --top-k 5 `
  --title "Data-derived graph-operator coverage review"
```

Apply the compiled library to a new target with the measured fallback policy:

```powershell
python -m core_retrosynthesis disconnect-operators `
  results/operator_retrosynthesis_poc/compact/all_shards_400/operator_library.json.gz `
  "Cc1ccnc(-c2ccccc2)c1" --top-k 10 --concise
```

### Strategy-grouped single-step results

The flat validated-candidate API remains available for benchmarks and detailed
precursor review. For chemist-facing single-step results, candidates can also
be grouped by the handle-independent strategy identity:

```python
from core_retrosynthesis import disconnect_strategies

strategies = disconnect_strategies(
    target_smiles,
    operator_library,
    top_k_strategies=10,
    max_realizations_per_strategy=3,
)
```

`STRAT1` is the deterministic combination of `OP1`, `SITE1`, and `SYN1`.
Different halides, boronates, protecting groups, templates, or other concrete
`REAL1` handle choices therefore consume one strategy slot and are returned as
nested realizations. Grouping happens only after the existing hard forward-graph
and operator-signature validation; incomplete or core-only candidates are not
promoted into strategy proposals. Strategy support uses the maximum independent
template support rather than a sum, avoiding double-counting overlapping
template evidence.

Exact precursor recovery remains intentionally separate from site, operator,
and synthon recovery. A different source-supported handle realization is not
called the recorded reaction, and condition compatibility is reported as not
evaluated when the query supplies no conditions.

## Full-scale operator library v3

The v3 builder is a resumable, shard-oriented path for datasets that are too
large to collect in memory. Each source shard writes three deterministic
artifacts under `OUTPUT/shards/`:

- a partial executable operator library;
- an admission JSONL ledger containing the exact rejection stage or accepted
  template/operator/completion identities;
- a manifest keyed by the source size, modification time, and build config.

Unchanged shard artifacts are reused. The merge uses `support.sqlite3` to
deduplicate observation and reference support exactly across shards, then
writes `operator_library_v3.json.gz` and `build_report.json`.

Run a progressive Compact build first:

```powershell
python -m core_retrosynthesis build-operators-full `
  datasets/literature `
  results/operator_retrosynthesis_poc/full_scale_v3/compact `
  --library-mode compact `
  --max-shards 50 --max-rows-per-shard 10 --workers 4
```

Then remove both limits and select Full for the complete corpus. Repeating the
same command resumes unchanged shards. Use `--force` only to intentionally
recompile them.

```powershell
python -m core_retrosynthesis build-operators-full `
  datasets/literature `
  results/operator_retrosynthesis_poc/full_scale_v3/full `
  --library-mode full `
  --workers 6
```

## Test a Compact operator build

Set the generated library path once for the following PowerShell examples:

```powershell
$library = "results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz"
$buildReport = "results/operator_retrosynthesis_poc/full_scale_v3/compact/build_report.json"
```

First confirm that both final artifacts exist and inspect the build census:

```powershell
Test-Path $library
Get-Content $buildReport |
  ConvertFrom-Json |
  Select-Object source_shards, source_rows, accepted_observations, `
    operator_count, completion_group_count, realization_count, template_count
```

For a quick functional smoke test, pass a product SMILES—not a reaction—to the
specificity-preserving L2-to-L1-to-L0 operator ladder. By default, each level
retrieves a pool four times the requested result count and then diversifies it
by graph operator, disconnection site, and synthon signature:

```powershell
python -m core_retrosynthesis disconnect-operators `
  $library `
  "Cc1ccnc(-c2ccccc2)c1" `
  --top-k 5 --concise
```

Every printed line is a proposed `precursors>>product` reaction. An empty result
is valid for an unsupported target, so test several products representative of
the chemistry present in the Compact batches. Remove `--concise` and request
diagnostics when investigating retrieval or validation:

```powershell
python -m core_retrosynthesis disconnect-operators `
  $library `
  "Cc1ccnc(-c2ccccc2)c1" `
  --top-k 5 `
  --max-templates 500 `
  --max-candidates-to-validate 100 `
  --diagnostics
```

Diversification is conservative: it round-robins chemistry-distinct groups only
within 0.05-wide structural-score bands and never moves an L1 or L0 proposal
ahead of an available L2 proposal. Use `--no-diversity` only for an ablation or
before/after comparison. The versioned defaults live in
`definitions/retrosynthesis_ranking.v3.json`; output candidates record the
policy ID, original structural rank, diversity rank, group key, and score band.
The same policy reserves up to two already-generated strategic alternatives
for outputs of at least five candidates, with a maximum two-band displacement.
The reserve never generates chemistry or bypasses forward validation.

The default ranking now adds a second conservative hierarchy after mandatory
forward validation and structural band assignment. It first ranks distinct
product disconnection sites (`SITE1`), then synthon skeletons (`SYN1`), and
finally precursor-handle realizations (`REAL1`). Completion evidence is derived
only from the operator library's independent-reference counts. A smoothed prior
backs off from `operator + SYN1` to `operator` and then to the global completion
distribution; unavailable evidence remains explicit and is not treated as a
failure. Sites and synthons are round-robin interleaved inside the same
abstraction-level and effective score-band partition, so the hierarchy cannot
promote L1/L0 over L2 or cross a structural/realism band. Every candidate
serializes its completion group, probability, backoff level, support
denominator, stage scores/ranks, and
`hierarchical_retrosynthesis_ranking.v2` definition ID. SITE1 includes the
deterministic strategic-complexity reduction score inside its existing
abstraction-level and structural-score-band partition. `--no-diversity`
disables both this hierarchy and the earlier diversity pass for a controlled
ablation. Use `--no-hierarchical-ranking` to retain the legacy diversity pass
and isolate the hierarchy's contribution. The policy is declared in
`definitions/hierarchical_ranking.v2.json`.

Strategic complexity is derived from mapped product bonds and product-derived
precursor atoms by
`../reactive_taxonomy/definitions/strategic_complexity.v1.json`.
The trace reports core fragmentation, ring-topology reduction, graph-complexity
change, stereochemical simplification, convergency, and tactical penalties.
Functional-group interconversions and protection changes remain valid but rank
as tactical; missing scaffold-level proposals are reported separately as an
operator-coverage warning.

Precursor compatibility is assessed before hierarchical ranking. Chemical
pair knowledge is declared in
`../reactive_taxonomy/definitions/reactive_pair_interactions.v1.json`; the
generic evaluator matches existing typed reactive sites and calculates graph
distance or possible closure size. Only pairs within one dot-separated
molecular component are intrinsic precursor warnings. The planner actions are
declared separately in `definitions/precursor_compatibility_policy.v1.json`.
The initial policy emits a strong warning and a structural-band demotion for
high or critical conflicts; the same groups on separate precursor components
receive no intrinsic molecular penalty.

The JSON output exposes product-index retrieval, SMARTS applicability,
generated precursors, validation attempts, operator mismatches, and valid
candidates. Retained candidates should have
`forward_validation_status: "verified_signature"`. Use `--skip-l0` to test only
the more specific L2/L1 tiers and `--no-context` as a ranking ablation.

That field is the inexpensive reaction-signature sanity check, not the complete
forward audit. The web runtime can additionally replay each retained operator
from its proposed precursors, run a target-blind product search, verify signature
and edit agreement plus reverse precursor recovery, and disclose competing
products. This independent audit is advisory and is serialized separately as
`forward_assessment`.

Next run a quantitative audit against a JSONL dataset excluded from every batch
used to build the library. For example, if the Dess–Martin observations were
not part of the Compact build:

```powershell
python -m core_retrosynthesis audit-operator-coverage `
  $library `
  "datasets/intermediate/DessMartin_periodinane_DMP_Alcohols__aldehydesketones.csv.observations.jsonl.gz" `
  "results/operator_retrosynthesis_poc/full_scale_v3/compact/heldout_audit" `
  --max-rows 100 `
  --top-k 25 `
  --max-templates 500 `
  --max-candidates-to-validate 100
```

The audit writes `coverage_audit.json` with aggregate exact, synthon, operator,
and site recall, plus `coverage_cases.jsonl.gz` with each target's successful
stage or precise failure stage. A file is genuinely held out only when its
observation/reference identities do not occur in any training batch; merely
choosing a different filename is not sufficient evidence of independence.

Finally run the focused deterministic regression suite:

```powershell
pytest -q tests/core_retrosynthesis_tests
```

## Rank verified disconnections with condition evidence

Condition recommendation is an optional application-layer step after the
operator ladder. Supply a condition index to recommend recipes for every
candidate whose forward validation status is `verified_signature`:

```powershell
$operatorLibrary = "results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz"
$conditionIndex = "datasets/literature/compact/generic_index.sqlite"

python -m core_retrosynthesis disconnect-operators `
  $operatorLibrary `
  "Cc1ccnc(-c2ccccc2)c1" `
  --top-k 5 `
  --condition-index $conditionIndex `
  --condition-top-k 3
```

Each candidate retains `retrosynthesis_rank` and gains
`condition_informed_rank` plus a `condition_evidence` object containing the
recommended canonical recipes, retrieval level, compatible and independent
support counts, precedent reference IDs, cautions, warnings, and errors.
Condition support is labeled as:

- `recommended_direct` for verified-signature condition evidence;
- `recommended_fallback` when type-agnostic, core, or topology fallback was
  required;
- `insufficient_evidence` when no compatible condition precedent was found.

Reranking is deterministic and bounded. Condition evidence can reorder
candidates only within the same abstraction level and 0.05-wide structural
score band, so a well-supported recipe cannot promote a substantially weaker
or less-specific disconnection. Inside that scope, direct evidence precedes
fallback evidence, which precedes insufficient evidence. Independent compatible
support, best-recipe reference support, and recipe score break ties before the
original retrosynthesis rank. The output records the policy ID, structural score
band, and reranking scope. Use `--keep-retrosynthesis-order` to attach the same
evidence without changing candidate order.

Graph validation remains authoritative. Candidates without a verified forward
signature are not sent to the condition recommender, and condition evidence
cannot rescue them. Conversely, missing condition evidence does not invalidate
a structurally verified disconnection; it records uncertainty caused by sparse
precedent coverage. Use `--condition-use-rxnmapper` only when the external
mapper is installed, and use `--condition-unrestricted-fallback` only when the
condition artifacts explicitly permit review-core access or trusted-index
reuse.

The command prints a progress heartbeat to stderr every 30 seconds while
keeping the final JSON report on stdout. Compilation messages include completed
shards, processed rows, accepted observations, reused checkpoints, throughput,
and elapsed time. Merge messages include partial libraries, rows, and template
counts. For example:

```text
[compile] 84/399 shards (21.1%), 69,412 completed rows, 41,506 accepted, 12 reused, completed-output average 20.4 rows/s, 6 active, 309 queued, elapsed 00:56:43
[merge] 120/399 shards (30.1%), 98,630 rows, 42,118 templates, elapsed 03:48:10
```

Rows from active workers are not counted until their entire shard is
checkpointed, so `completed-output average` is a conservative, bursty measure,
not live per-row throughput. A heartbeat with no newly completed shard reports
the number of active and queued shards without showing a declining rate.

Use `--progress-interval 60` to report once per minute, or
`--quiet-progress` to suppress progress messages.

Admission no longer requires an already serialized core and observation. The
compiler recomputes mapping, observation, core, and completeness from the
source reaction, but admission still requires a verified recomputed core,
verified product completeness, and source round-trip. Missing serialized data
is recovered only when those checks pass. Rejections are separated into source,
mapping, observation, core, completeness, canonicalization, operator,
template, and round-trip stages.

Templates are organized at four data-derived levels:

- graph operator;
- handle-completion group (`operator + direct precursor-handle signature`);
- executable realization identified by its normalized precursor-handle
  subgraph;
- L0/L1/L2 product-context SMARTS.

Precedents are retained deterministically across distinct context bins instead
of keeping the first observations encountered. A necessary-feature product
index uses product-observable edited atom/bond tokens before SMARTS matching;
it cannot exclude a chemically applicable template, and final candidates still
pass mapped graph-edit validation.

Audit a genuinely held-out shard or dataset with:

```powershell
python -m core_retrosynthesis audit-operator-coverage `
  results/operator_retrosynthesis_poc/full_scale_v3/full/operator_library_v3.json.gz `
  path/to/heldout.jsonl.gz `
  results/operator_retrosynthesis_poc/full_scale_v3/full/heldout_audit `
  --max-rows 1000 --top-k 25
```

The audit attributes misses to source compilation, operator absence, product
index retrieval, product SMARTS applicability, precursor generation,
structural validation, or global ranking. It also writes exact, synthon,
operator, and site recall plus per-target `coverage_cases.jsonl.gz`. Do not use
a training shard as the audit source.

The completed 399-shard, one-row-per-shard smoke build admitted 239 unique
observations (59.9%) after removing five cross-shard duplicates. It produced 72
operators, 105 handle-completion groups, 156 realizations, and 463 templates;
an unchanged resumability rerun completed in under one second. A denser
50-shard x 10-row progressive validation admitted 349/500 rows
(69.8%) into 37 operators, 56 handle-completion groups, 116 realizations, and
382 templates. Normalized local map identity reduced that library from the
initial 321 realizations and 869 templates. The detailed rejections were 111 unavailable mappings, 31
unverified recomputed cores, and 9 source round-trip failures. On 25 reactions
from an unused shard, candidate/operator/site coverage was 1.0 and synthon
recall was 0.92; exact precursor recall was zero because the library proposed a
different source-supported alkylating handle. This small, single-shard audit is
a pipeline regression, not a diverse accuracy estimate.

Candidate validation now consumes the atom correspondence returned by
RDChiral. This prevents unchanged duplicate atoms, such as an aryl bromine and
a bromide leaving group in the same proposal, from being swapped by a second
global remapping pass and incorrectly creating extra graph edits.

## Condition-aware selectivity POC

`selectivity_poc.py` tests a dataset-driven alternative to hand-authored
functional-group compatibility rules. For an observed single-connection
replacement, it holds the partners and edit topology fixed, enumerates other
connection endpoints on the selected partner component, applies each graph
edit, and retains only products that pass RDKit sanitization. Named reaction
families do not participate.

Every observed reaction becomes a listwise choice set. The mapped endpoint is
the selected outcome; other available endpoints are weak `not observed as
major` alternatives rather than asserted failures. A small deterministic
hashed softmax model learns endpoint and endpoint-by-condition preferences:

```python
from core_retrosynthesis import (
    ConditionalEditChoiceModel,
    build_reaction_choice_set,
)

acidic_s_alkylation = build_reaction_choice_set(
    "Cc1[nH]cnc1CCl.NCCS>>NCCSCc1c(C)[nH]cn1",
    {"medium": "acidic", "salt_state": "hydrochloride"},
    reference_id="example-s",
    label_strength=1.0,
)
neutral_n_alkylation = build_reaction_choice_set(
    "Cc1[nH]cnc1CCl.NCCS>>SCCNCc1c(C)[nH]cn1",
    {"medium": "neutral", "salt_state": "free_base"},
    reference_id="example-n",
    label_strength=1.0,
)

model = ConditionalEditChoiceModel()
model.fit((acidic_s_alkylation, neutral_n_alkylation))
assessment = model.assess(acidic_s_alkylation)
```

Canonical converted dataset rows can be consumed without source-column logic
using `build_reaction_choice_set_from_record(row)`. Its resolved-recipe
projection retains canonical component identities, roles, temperature, time,
concentration, atmosphere, and stages while excluding recipe IDs, source
fields, warnings, and provenance that could leak document identity.

The assessment reports the desired and best-competitor probabilities,
probability margin, normalized entropy, exact-condition independent reference
support, and every ranked counterfactual product. This is an architectural POC,
not a calibrated production predictor. Its current graph scope is one broken
and one formed heavy-atom bond sharing an electrophilic endpoint. It does not
yet model no-reaction outcomes, stoichiometric overreaction, multiple events,
or unexpected transformations outside the enumerated choice set.

The structural enumerator is integrated into generic retrosynthesis as a
review-only audit. A validated candidate receives a
`POSSIBLE_FUNCTIONAL_GROUP_COMPETITION` warning when another endpoint on the
same partner produces a distinct sanitized product. The warning serializes the
selected endpoint and each alternative product into CLI/API output, and the web
review displays those structures. It is explicitly condition-unaware and is
not used by candidate admission, scoring, deduplication, or diversity ranking.
Unsupported topologies fail open without changing the retrosynthesis result.

## Bounded multistep route search

The initial multistep planner recursively composes the existing validated
L2-to-L1-to-L0 operator ladder without changing single-step chemistry. Search
depth is restricted to two or three reactions along the longest branch. A
route is solved only after at least one disconnection and when every leaf:

- has RDKit molecular weight at or below the configured threshold (150 by
  default); or
- exactly matches the configured local stock portfolio with physical or
  supplier-in-stock evidence; or, as a legacy fallback, matches the canonical
  literature-molecule index with reactant-like source-role provenance.

Build the reusable SQLite index as documented in `cas_tools/README.md`, then
run:

```powershell
python -m core_retrosynthesis plan-routes `
  results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz `
  cas_tools/data/stock_portfolio.sqlite `
  "TARGET_SMILES" `
  --max-depth 3 `
  --top-k-routes 5
```

Catalyst, solvent, reagent, product, and untyped legacy-index matches are not
strong literature terminals. Rebuild the index with `source_role` provenance;
the planner exposes an explicit compatibility switch for legacy untyped
catalogs but does not enable it by default. Low-molecular-weight leaves remain
labelled heuristic terminals and add a declarative route-cost penalty.

The target itself is never accepted as a zero-step stock, literature, or
molecular-weight terminal. The deterministic best-first beam expands the
largest unresolved leaf first, caches repeated molecule expansions, reserves
specificity-first L2 candidates and widens to L1/L0 only when a narrower tier
cannot fill the candidate quota, rejects ancestor cycles, and reports
depth-, candidate-, and search-limited partial routes separately from solved
routes. It searches the configured budget before selecting the lowest-cost
top-k solutions and retains bounded alternative paths to the same molecular
state. Step costs, fallback and selectivity penalties, heuristic-terminal
penalties, precursor-compatibility and realism bands, condition-availability
penalties, strategic-progress deficits, and path breadth come from
`multistep_ranking.v3.json`. Literature
presence is disclosed as a corpus observation rather than
commercial-availability evidence.

Multi-step expansion uses staged action validation: indexed templates and
RDChiral precursor proposals are ranked cheaply, then hard forward-graph and
operator-signature validation proceeds only until the planner's per-step
action quota is filled. The exhaustive one-step API remains the default for
interactive disconnection review. Search diagnostics distinguish proposed,
validated, and accepted actions; record the first-solution expansion, dead
ends, beam pruning, cache use, and abstraction levels attempted; and the web
API adds elapsed time plus the configured search budget.

Every route also serializes a `route_tree` alternating molecule occurrences
and validated reaction nodes. Occurrence IDs prevent identical molecules on
different branches from collapsing. `route_postprocessing` exposes stable tree
IDs and a deterministic multiset-Jaccard distance matrix, providing a canonical
contract for route clustering and diversity selection without changing the
chemistry-first cost ranking in this release.

### Deterministic route refinement

The assistance layer can turn typed route issues into one bounded refinement
request. The model may select a closed objective and ask for either an alternate
disconnection or an alternate realization, but it cannot provide reaction
SMILES, atom edits, conditions, or replacement precursors. The deterministic
planner derives the excluded candidate identity from the stored route, reruns
validated search and condition recommendation, and reports issue-count changes
with lineage back to the original route. Both the original and refined results
remain available for comparison.

Reaction-regime compatibility is also derived from molecular evidence. The
initial rule recognizes a linked C-halogen cleavage, C-C formation, and
carbonyl C=O-to-C-O change, then checks unchanged spectator groups for a free
O-H, N-H, or S-H handle. It reports a strong protic-quench warning with the
matched edit and spectator provenance. This is an inferred strongly basic
carbon-transfer regime, not a named-reaction assignment or proof of specific
conditions; alternative realizations remain eligible for deterministic search.
The assessment is attached to each validated candidate before ranking. Its
versioned policy applies a structural-band demotion in one-step ranking and an
explicit route-cost penalty in multistep search, so raw similarity cannot hide
a known strong conflict.

Each typed step issue also exposes deterministic repair proposals. Alternate
realization and alternate disconnection proposals are actionable because they
reuse the validated search engine. A temporary-masking proposal is reported as
unavailable until registered protection and deprotection sequence operators can
construct and forward-validate every added structure. This explicit unavailable
state prevents the model from presenting an invented protected route as a
system result.

The assistance-facing route tools reuse these contracts rather than exposing
general RDKit execution. Step-precedent lookup returns only source-round-tripped
observations attached to the selected admitted template. Whole-route
verification independently checks route-tree integrity, target and step graph
identity, forward-validation status, terminal leaves, deterministic chemistry
issues, and condition-evidence coverage. Neither tool generates candidates or
changes route ranking.

### Optional precedent-route action policy

Candidate replay from observed route trees can train a deterministic listwise
policy over validated single-step actions. It is a residual over the existing
candidate order: it cannot generate candidates, bypass graph validation, or
change terminal-material rules. Exact operator/template/precursor identities
are excluded from the default features, and an artifact has zero planner
influence until its held-out validation population satisfies the versioned
activation gate.

```powershell
python -m core_retrosynthesis train-route-action-policy REPLAY MODEL `
  --report REPORT --overwrite

python -m core_retrosynthesis plan-routes LIBRARY STOCK_INDEX TARGET `
  --route-action-policy MODEL
```

Planner JSON includes the model and definition IDs, residual scale, activation
status, per-step policy diagnostics, and the number of reordered expansions.
The initial 50-route artifact is an inactive safety POC. A deterministic
500-route artifact passes the validation activation gate and improves held-out
next-action ranking, but its fixed complex-target A/B does not improve solved
routes, expansion count, or observed-strategy recovery. It therefore remains
optional and experimental rather than becoming the default planner policy. See
`docs/new/multistep_route_action_policy_poc.md` for the split metrics,
whole-route result, artifacts, limitations, and next evidence gate.

Before application use, run `calibrate-route-action-policy` on a fixed
validation-route panel. The whole-route gate compares every versioned residual
scale by solved targets, observed-strategy recovery, terminal progress, and
search effort, preferring zero influence on ties. The current 500-route model
fails that gate and is frozen at residual scale zero; the untouched test panel
then reproduces baseline routes exactly. Raw trained models are research
artifacts, while route-calibrated models are the deployment-facing contract.

The same status document links a 12-target familiar-chemistry comparison panel.
The reusable `multistep_panel_review` module renders baseline and policy routes
side by side with inline reaction drawings, search diagnostics, filters,
persistent accept/question/reject decisions, notes, and review JSON export.

### Two-step precedent-route chemical-space expansion

The route-expansion POC asks a different question from action ranking: given an
observed two-step route, which products are structurally reachable when declared
building blocks are substituted and each step is executed forward? It reports a
cumulative ladder from exact replay (`R0`), through exact-context L2 expansion
(`R1`), to explicit local-context relaxation with L1 operators (`R2`).

```powershell
python -m core_retrosynthesis expand-precedent-routes `
  core_retrosynthesis/definitions/two_step_observed_route_expansion_poc.v1.json `
  results/core_retrosynthesis/precedent_route_expansion/two_step_observed_poc.v1.json `
  --stock-index cas_tools/data/stock_portfolio.sqlite `
  --route-core-source datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.core.v1.jsonl.gz
```

Every retained pathway passes the existing forward product, reverse recovery,
and operator-edit agreement checks. The primary panel contains three
source-verified patent-route segments and exact supplier-stock-verified
substitutes. It yields 3 exact R0, 16 cumulative R1, and 33 cumulative R2
products. These are structural possibilities rather than experimental
feasibility or yield predictions. The former N/O/S graph panel remains a fast
regression fixture. See
`docs/new/two_step_precedent_route_expansion_poc.md` for the design, result, and
next evidence gate.
