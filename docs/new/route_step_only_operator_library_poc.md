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
rejected from the strict baseline as `materialized_core_not_verified` because
the recomputed core status is `review`.

An opt-in `validated_departures` compiler policy now distinguishes this narrow
case from general review evidence. It requires complete active-atom mapping,
verified product completeness, no blocking reason, all retained edits to be
graph-checked, and every unchecked edit to be a mapped broken bond from a
retained atom to an atom omitted from the reported product. The default remains
`pass_only`.

The policy exposed and fixed a template-materialization defect: product-absent
mapped atoms must become unmapped precursor handles before RDChiral execution.
After that fix, nitration exactly recovers the observed substrate and nitric
acid at L0, L1, and L2. Nitro reduction exactly recovers the observed nitro
precursor at L0; its L1 and L2 templates produce no outcome and are discarded.
Both reactions therefore pass the mandatory source round trip under the opt-in
policy, but they are not part of the unchanged strict baseline library.

The per-level, provenance-bearing result is:

```text
results/core_retrosynthesis/route_step_operator_library/v1/
  validated_departure_round_trip_audit.v1.json
```

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

python -m core_retrosynthesis audit-operator-round-trips `
  results/core_retrosynthesis/route_step_operator_library/v1/train_conversion_final/shards `
  results/core_retrosynthesis/route_step_operator_library/v1/validated_departure_round_trip_audit.v1.json `
  --reaction-id US04011323:21308_0 `
  --reaction-id US04011323:21309_0 `
  --core-admission-policy validated_departures
```

## Next evidence gate

The separate experimental build is complete:

```text
results/core_retrosynthesis/route_step_operator_library/v1/
  operators_validated_departures/
```

| Metric | Strict `pass_only` | `validated_departures` |
|---|---:|---:|
| accepted observations | 966 (6.47%) | 14,218 (95.21%) |
| templates | 477 | 7,298 |
| operators | 104 | 644 |
| realizations | 149 | 1,971 |
| completion groups | 104 | 985 |

The experimental rejections are 276 `materialized_core_not_verified`, 411
`source_round_trip_failed`, 24 `missing_generic_operator_signature`, 3
`materialized_core_missing`, and 1 `generic_l0_compilation_failed`. The large
admission change confirms that product-omitted departing fragments dominate
route-core review status; it also makes held-out validation mandatory before
promotion.

### Library-independent held-out comparison

Panel selection is now separated from library coverage. The frozen panel has a
content-derived ID and contains 12 distinct validation/test patents, six
`handle_progression` and six `same_site_coupled` cases. No exact target occurs
in training; two targets have a training Murcko-scaffold match. The panel forces
in the recurrent nitration-reduction strategy so that the original capability
gap cannot disappear through case filtering.

```text
results/core_retrosynthesis/coupled_strategy_evaluation/
  route_only_fixed_panel.v1.json
  route_only_strict_fixed_panel.v1.json
  route_only_strict_fixed_panel.v1.html
  route_only_validated_departures_fixed_panel.v1.json
  route_only_validated_departures_fixed_panel.v1.html
```

| Fixed-panel metric | Strict route-only | Validated departures |
|---|---:|---:|
| selected cases | 12 | 12 |
| operator-pair coverage gaps | 12 | 0 |
| ordinary depth-two recovery | 0/12 | 1/12 |
| promoted v1 pair recovery | 0/12 | 9/12 |
| promoted validation attempts | 0 | 399 |

The held-out nitration-reduction case is from validation patent `US08598374B2`.
Its exact target and Murcko scaffold are absent from strategy training. The
experimental pair-guided search recovers the observed nitro intermediate at
rank 1; strict search has neither constituent operator.

Three experimental cases remain unrecovered. Two return alternate validated
macro actions but not the observed intermediate within top 3, while one returns
no macro action under the fixed budget. These are search/selectivity or
realization misses rather than operator-coverage gaps.

The evidence supports retaining `validated_departures` as an experimental
route-library policy and using it for the next unified-library comparison. It
does not yet justify replacing the production `pass_only` default: prospective
invalid-product, selectivity, and broader route-quality measurements remain
required.
