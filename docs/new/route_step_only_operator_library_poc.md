# Route-Step-Only Operator Library POC

## Outcome

The canonical multistep-route corpus can now be projected into the same
`source_observation.v1` intermediate contract used by the standard condition
dataset converter. The extractor emits one physical reaction observation per
unique source reaction, retains every route occurrence as provenance, and
preserves the patent-level train, validation, and test split.

The first strict operator build is intentionally conservative. It proves the
end-to-end path but accepts only 966 of 14,933 training observations because the
current compiler requires a recomputed reaction core with status `pass`. Most
route cores are structurally supported but have status `review`, commonly when
not every edit can be graph-checked. This POC does not weaken that gate.

## Extraction contract

`core_retrosynthesis.route_step_observations` consumes paired canonical
`ReactionRouteTree v2` and `RouteCoreProjection v1` artifacts. It:

- requires routes with at least two physical reactions;
- validates tree/core identity, reaction counts, patent IDs, and splits;
- rejects a source reaction observed across data splits;
- writes separate patent-preserved train, validation, and test files;
- converts three-field reaction SMILES to the two-field mapped transformation
  expected by the intermediate converter;
- retains reaction-middle-field structures as condition claims without
  assigning unsupported reagent, catalyst, or solvent roles;
- records route, tree, core, reaction node, molecular occurrence, depth,
  internal-precursor, and terminal-precursor provenance; and
- produces deterministic gzip bytes and a hash manifest.

Only the training artifact is eligible to build operators. Validation and test
observations are reserved for later coverage and retrosynthesis evaluation.

## POC corpus

Inputs:

```text
datasets/external/higher_level_retrosynthesis/figshare_v2/curated/
  routes.poc.tree.v2.jsonl.gz
  routes.poc.core.v1.jsonl.gz
```

Output manifest:

```text
datasets/intermediate/route_steps_poc_v1/
  route_step_observation_manifest.json
```

The 5,000 multistep routes contain 18,647 unique physical step occurrences:

| Split | Steps | Patents |
|---|---:|---:|
| train | 14,933 | 4,003 |
| validation | 1,987 | 525 |
| test | 1,727 | 472 |

Route-core quality before standard conversion is 1,233 `pass`, 17,085
`review`, 326 `blocked`, and 3 `unavailable`.

## Standard conversion

The training artifact was processed through the unmodified sharded generic
converter. Its final report is:

```text
results/core_retrosynthesis/route_step_operator_library/v1/
  train_conversion_final/conversion_report.json
```

All 14,933 input rows were emitted and the integrity check found no duplicate
observations. Chemistry is `verified` for 14,930 rows and `rejected` for 3.
There are 14,930 reaction signatures, 4,003 patent references, and 14,128 rows
eligible for the recommendation index. Condition evidence remains the principal
data limitation: 13,543 rows retain unresolved conditions, 588 are completely
resolved, and 802 are unusable. Route reaction middle fields do not justify
inventing condition roles.

## Strict executable-operator build

The provenance-correct training shards produced:

```text
results/core_retrosynthesis/route_step_operator_library/v1/operators_final/
  build_report.json
  operator_library_v3.json.gz
  support.sqlite3
```

| Metric | Result |
|---|---:|
| source observations | 14,933 |
| accepted observations | 966 (6.47%) |
| templates | 477 |
| operators | 104 |
| realizations | 149 |
| completion groups | 104 |

The compiler rejected 13,937 observations as
`materialized_core_not_verified`, 27 after `source_round_trip_failed`, and 3 as
`materialized_core_missing`. The low acceptance is an executable-template
admission result, not an extraction failure: the standard converted reaction
records remain available for analysis and future policy work.

The motivating two-step aryl sequence is present as source reactions `21308_0`
(nitration) and `21309_0` (nitro reduction) under patent `US04011323`. Both have
verified stored completeness and supplied-map confidence 1.0, but both are
currently rejected as `materialized_core_not_verified` because the recomputed
core status is `review`. Consequently, this strict route-only library does not
yet expand retrosynthesis coverage for that example.

## Reproduction

```powershell
python -m core_retrosynthesis extract-route-step-observations `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.tree.v2.jsonl.gz `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.core.v1.jsonl.gz `
  datasets/intermediate/route_steps_poc_v1

python -m condition_recommender.sharded_conversion_cli `
  datasets/intermediate/route_steps_poc_v1/route_steps.train.observations.jsonl.gz `
  results/core_retrosynthesis/route_step_operator_library/v1/train_conversion_final `
  --shard-size 1000 --mode full --workers 4 --checkpoint-interval 5

python -m core_retrosynthesis build-operators-full `
  results/core_retrosynthesis/route_step_operator_library/v1/train_conversion_final/shards `
  results/core_retrosynthesis/route_step_operator_library/v1/operators_final `
  --workers 4 --quiet-progress
```

## Next evidence gate

The next iteration should define and test an explicit route-core review policy,
instead of silently treating all review records as pass. A candidate policy
should require validated supplied mapping, verified reaction completeness, no
blocking reasons, an enumerated non-blocking review reason, and successful
forward round-trip reconstruction. It should be measured on validation routes
for accepted-operator coverage, exact observed-precursor recovery, invalid
product rate, and the two-step strategy examples before changing the production
compiler contract.
