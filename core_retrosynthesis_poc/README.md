# Reaction-core-derived retrosynthesis POC

This standalone POC compiles executable retrosynthetic SMARTS from the
structure-derived `ReactionCoreProjection`, normalized edits, source molecular
graphs, and attachment ports. It is designed for a controlled comparison with
the existing `retrosynthesis_poc`, whose SMARTS are selected by the RDChiral
template extractor.

The core remains an observation rather than an executable prediction. This
package is a separate compiler and application layer; it does not change
`reactive_taxonomy`, `condition_registry`, or `condition_recommender`.

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

## Build a balanced library

Use a separate limit for every included cohort so C-N rows do not consume a
global limit before C-O and C-S are read:

```powershell
python -m core_retrosynthesis_poc build-library `
  datasets/literature/shards `
  results/core_retrosynthesis_poc/core_templates.json.gz `
  --include "C_N_Coupling*.jsonl.gz" `
  --include "C_O_Coupling*.jsonl.gz" `
  --include "C_S_Coupling*.jsonl.gz" `
  --max-rows-per-include 1000
```

Use `--level L1` or `--level L2` to build only one abstraction level.

## Disconnect a target

```powershell
python -m core_retrosynthesis_poc disconnect `
  results/core_retrosynthesis_poc/core_templates.json.gz `
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
python -m core_retrosynthesis_poc disconnect-ensemble `
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
python -m core_retrosynthesis_poc compare `
  datasets/literature/shards `
  results/core_retrosynthesis_comparison `
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
python -m core_retrosynthesis_poc render-report `
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
python -m core_retrosynthesis_poc compare-diverse `
  datasets/literature/shards `
  results/diverse_retrosynthesis_poc/balanced_round_robin_100 `
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
python -m core_retrosynthesis_poc render-report `
  results/diverse_retrosynthesis_poc/balanced_round_robin_100/comparison.json `
  results/diverse_retrosynthesis_poc/balanced_round_robin_100/review.html `
  --method baseline --method core_context `
  --method ensemble_baseline_core_context --top-k 5
```

Apply the generic core library directly with, for example:

```powershell
python -m core_retrosynthesis_poc disconnect-generic `
  results/diverse_retrosynthesis_poc/balanced_round_robin_100/core_templates.json.gz `
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
python -m core_retrosynthesis_poc compare-stress `
  datasets/literature/shards `
  results/diverse_retrosynthesis_poc/stress_cc_michael_250 `
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
python -m core_retrosynthesis_poc render-report `
  results/diverse_retrosynthesis_poc/stress_cc_michael_250/comparison.json `
  results/diverse_retrosynthesis_poc/stress_cc_michael_250/review.html `
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
python -m core_retrosynthesis_poc compare-operators `
  datasets/literature/shards `
  results/operator_retrosynthesis_poc/all_shards_400 `
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
python -m core_retrosynthesis_poc render-report `
  results/operator_retrosynthesis_poc/all_shards_400/specificity_comparison.json `
  results/operator_retrosynthesis_poc/all_shards_400/review.html `
  --method supported_l1_l2 --method data_l2 `
  --method data_l2_l1 --method data_ladder --top-k 5 `
  --title "Data-derived graph-operator coverage review"
```

Apply the compiled library to a new target with the measured fallback policy:

```powershell
python -m core_retrosynthesis_poc disconnect-operators `
  results/operator_retrosynthesis_poc/all_shards_400/operator_library.json.gz `
  "Cc1ccnc(-c2ccccc2)c1" --top-k 10 --concise
```

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

Run a progressive build first:

```powershell
python -m core_retrosynthesis_poc build-operators-full `
  datasets/literature/shards `
  results/operator_retrosynthesis_poc/full_scale_v3 `
  --max-shards 50 --max-rows-per-shard 10 --workers 4
```

Remove both limits for the full corpus. Repeating the same command resumes
unchanged shards. Use `--force` only to intentionally recompile them.

```powershell
python -m core_retrosynthesis_poc build-operators-full `
  datasets/literature/shards `
  results/operator_retrosynthesis_poc/full_scale_v3 `
  --workers 4
```

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
python -m core_retrosynthesis_poc audit-operator-coverage `
  results/operator_retrosynthesis_poc/full_scale_v3/operator_library_v3.json.gz `
  path/to/heldout.jsonl.gz `
  results/operator_retrosynthesis_poc/full_scale_v3/heldout_audit `
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
