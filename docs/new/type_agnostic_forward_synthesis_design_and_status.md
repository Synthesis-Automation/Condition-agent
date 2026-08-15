# Type-Agnostic Forward Synthesis: Design and Status

**Status date:** 2026-08-15  
**Primary implementation:** [`forward_synthesis`](../../forward_synthesis/)  
**Release posture:** deterministic structural predictor; route impact advisory

## Purpose

Forward synthesis fills the structural gap between condition recommendation and
retrosynthesis. It accepts starting materials, optionally accepts a resolved
condition recipe, and returns possible main products. Its principal route use is
to challenge a retrosynthetic step with products that can compete with the
intended target.

The target is never required by blind generation. When an intended target is
available, targeted replay and blind competition are executed separately so the
answer cannot be manufactured from the proposed target.

## Package boundary

```text
reactive_taxonomy       condition_registry
        |                       |
        +----> forward_synthesis <---- condition_recommender
                        ^
                        |
              core_retrosynthesis
```

`reactive_taxonomy.reaction_operators` owns the direction-neutral operator
contract and deterministic multi-component application. `forward_synthesis`
owns precursor retrieval, forward admission, generation, validation, ranking,
competition grouping, evaluation, and explanations. `condition_recommender`
owns supplied-recipe compatibility. `core_retrosynthesis` consumes an advisory
assessment and does not own forward chemistry.

No new-system package imports legacy `chemtools`.

## Implemented data flow

1. Project each generic graph template into explicit precursor and product
   sides.
2. Admit it only when forward application reconstructs every retained source
   product and reverse application recovers its contributing precursors.
3. Build a conservative precursor index from component counts and necessary
   fixed elements.
4. Enumerate all injective matches between precursor templates and supplied
   components. Unmatched input components remain explicit nonparticipants.
5. Apply the graph transformation with RDKit and reject unsanitizable products.
6. Re-featurize each generated reaction through `reactive_taxonomy`.
7. Require a reaction signature, verified reaction core, generating edit-token
   agreement, and reverse precursor recovery.
8. Apply hard recipe conflicts before scoring when a resolved recipe is
   supplied.
9. Rank with versioned structural evidence, group duplicate products, and
   preserve alternative operator/template pathway identities.

The current ranking definition is
[`forward_ranking.v1.json`](../../forward_synthesis/definitions/forward_ranking.v1.json).
Its score is explicitly an uncalibrated deterministic priority.

## Public contracts

- `BidirectionalReactionOperator` (`OP1` chemistry plus `FOP1` directional
  execution identity)
- `ForwardOperatorLibrary`
- `ForwardProductCandidate`
- `ForwardCompetitionGroup`
- `ForwardPredictionResult`
- `RouteStepForwardAssessment`
- `RouteForwardAssessmentReport`
- `ForwardReplayEvaluationReport`

`predict_products` is blind. `assess_proposed_step` adds an intended product and
optional operator hint only to the separate targeted replay. Its blind result is
unchanged by those values.

## Retrosynthesis integration

For every planned route reaction, the route adapter obtains the child precursor
occurrences and independently evaluates the parent molecule. Results distinguish:

- `clear`: intended product leads the enumerated set by the guarded score margin;
- `competitive`: intended product is present but a close or higher product exists;
- `unsupported`: targeted chemistry may be plausible but blind retrieval did not
  recover the intended product;
- `structurally_inconsistent`: the proposed operator cannot reproduce the target;
- `out_of_scope`: the query or current operator library cannot be assessed.

These labels are review evidence only. They do not currently reject or rerank a
route.

## Evidence limits

- The initial scope is a single reaction event and one main organic product.
- Named reaction families never generate or admit a candidate.
- Missing conditions mean that returned products are possibilities, not a major
  product prediction.
- Recipe compatibility is neither a yield estimate nor proof of product
  formation.
- Patent and literature rows usually report one isolated product. Generated
  alternatives are unobserved counterfactuals, not labelled negative outcomes.
- Stereo-relaxed matches are disclosed separately from exact product matches.
- Cascades, rearrangements, unusual stoichiometry, and unsupported atom-origin
  cases may abstain.

## Evaluation gate

`evaluate_product_hidden_replay` removes the product from each held-out query and
reports exact product top-1/top-5/top-10 recovery and reciprocal rank. Reference
overlap with operator precedents fails by default. Production selectivity claims
additionally require scaffold- and time-disjoint panels with experimentally
observed competing outcomes or negative/high-throughput data.

Route-ranking influence remains disabled until evaluation establishes product
recall, competition recall, false-warning rate, abstention behavior, and runtime
on an untouched route panel.

## Web integration

The local workbench exposes Forward synthesis as a first-class analysis mode at
`POST /api/v1/forward-synthesis`. Starting materials are required; intended
product, retrosynthesis operator hint, resolved recipe, broad L0 fallback, result
count, and Full/Compact operator-library mode are explicit request fields.

Without an intended product the response contains blind predictions. With an
intended product it also contains a separate route-step assessment while keeping
the target out of the blind prediction pass. The UI presents product ranks,
valid pathway and operator counts, alternative pathway identities, competition
groups, condition compatibility, graph/operator traces, and the audit
disposition. Results remain locally exportable as versioned JSON.
