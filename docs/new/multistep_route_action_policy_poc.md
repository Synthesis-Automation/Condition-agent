# Data-Driven Multistep Route-Action Policy POC

**Status date:** 2026-08-14

## Outcome

The multistep planner can now use precedent routes to rank its already
validated single-step disconnections. The integration is optional and does not
change template generation, graph validation, terminal-material checks, or the
existing planner when no model is supplied.

The first 50-route training artifact exercises the complete path from route
replay to multistep search, but it is deliberately inactive. Only five
validation choice sets were available, below the versioned activation minimum
of 20. Its learned residual therefore has zero planner influence. This is a
successful plumbing and safety POC, not evidence that route quality improved.

## Design

```text
observed route trees
        |
        v
route-action labels + bounded candidate replay
        |
        v
listwise choice sets (chosen precedent action vs validated alternatives)
        |
        v
held-out-gated route-action policy
        |
        v
rerank only graph-validated one-step candidates
        |
        v
existing best-first multistep search and terminal checks
```

The chemistry boundary is strict:

- deterministic single-step search creates and validates every candidate;
- the policy cannot create a reaction, bypass validation, or make a terminal;
- an unchosen candidate is a relative alternative, not a chemical negative;
- exact observed precursors are strong labels and strategy-equivalent actions
  are confidence-weighted labels;
- train, validation, and test route splits are retained; and
- validation selects the learned residual strength, while test is audit-only.

The policy uses coarse candidate properties, strategic-complexity evidence,
search depth, and product Morgan bits. Exact operator, template, precursor, and
synthon identities are excluded from the default feature definition to reduce
memorization. Feature groups, optimization values, planner cost weight, and
the activation threshold are versioned in
`core_retrosynthesis/definitions/route_action_policy.v1.json`.

## Implemented contracts

`core_retrosynthesis.route_action_policy` owns:

- `RoutePolicyExample`, a deterministic listwise choice set;
- `RouteActionPolicyModel`, a sparse serializable hashed-softmax residual model;
- replay conversion, training, held-out evaluation, and deterministic gzip
  model I/O; and
- baseline top-1 and mean reciprocal rank metrics for every retained split.

`plan_multistep_routes(..., route_action_policy=model)` records:

- model and definition identities, residual scale, and activation status;
- per-step probability, learned rank, and original single-step rank;
- the optional probability-deficit cost component; and
- scored-action and reordered-expansion diagnostics.

Schema `multistep_route.v1.5` owns these fields. A model with residual scale
zero has a zero planner cost weight and cannot alter route ordering.

## First bounded result

Source replay:

```text
results/core_retrosynthesis/route_action_evaluation/
  routes.poc.random50.action_replay.bounded.v3.jsonl.gz
```

The replay contains 50 routes and 188 reported steps. It yielded 81 usable
choice sets: 57 exact-precursor labels and 24 strategy-equivalent labels. The
split was 65 train, 5 validation, and 11 untouched test examples.

The learned training objective decreased from 1.50176266 to 1.21677395. The
five-example validation set was too small for activation, so the selected
residual scale is `0.0`. Consequently held-out scores equal the existing-rank
baseline:

| Split | Examples | Baseline top-1 | Baseline MRR |
| --- | ---: | ---: | ---: |
| Validation | 5 | 0.6000 | 0.8000 |
| Test | 11 | 0.4545 | 0.6848 |

Model and report:

```text
results/core_retrosynthesis/route_policy/
  random50.route_action_policy.v1.json.gz
  random50.route_action_policy.v1.report.json
```

The model ID is
`RAPM1:2481c8795b1c61bbe80215d9235ea221ca0c33739b48bfaae0fa89107acb2640`.
Its artifact SHA-256 is
`f277d3b55d09d9ab53c33fe954cfab1e04cb59bafb3e6016f6c67efa189c9e5b`.

A bounded depth-two planner exercise loaded this model and scored 12 validated
actions across four molecule expansions. It reordered zero expansions, as
required for an inactive artifact. The small search budget produced two
partial routes and no solved route in both baseline and policy modes; this was
an integration check, not a route-quality benchmark.

## Commands

Train a deterministic model from a candidate-replay artifact:

```powershell
python -m core_retrosynthesis train-route-action-policy `
  results/core_retrosynthesis/route_action_evaluation/routes.poc.random50.action_replay.bounded.v3.jsonl.gz `
  results/core_retrosynthesis/route_policy/random50.route_action_policy.v1.json.gz `
  --report results/core_retrosynthesis/route_policy/random50.route_action_policy.v1.report.json `
  --overwrite
```

Load it into multistep planning:

```powershell
python -m core_retrosynthesis plan-routes `
  results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz `
  cas_tools/data/stock_portfolio.sqlite `
  "TARGET_SMILES" `
  --max-depth 3 `
  --top-k-routes 5 `
  --route-action-policy results/core_retrosynthesis/route_policy/random50.route_action_policy.v1.json.gz
```

## Next evidence gate

The next useful experiment is larger leakage-controlled replay, not a more
complex model. It should provide at least 20 validation choice sets and a
meaningful untouched test set, preferably from patent-group and scaffold-aware
splits. Activate a learned residual only if validation improves over the
existing ranker, then compare whole-route recovery, solved-route rate, search
efficiency, route diversity, and chemist preference on test routes.

The 5,000-route labels remain useful positive supervision, but training against
planner choices requires replayed alternatives. Scaling replay and excluding
source-patent or close-analogue precedents are therefore the next data tasks.
