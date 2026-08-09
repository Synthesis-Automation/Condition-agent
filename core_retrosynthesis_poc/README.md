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
