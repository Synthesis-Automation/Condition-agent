# Condition recommendation corrections — 2026-09-12

The review reproduced five deterministic defects and one uncertainty-reporting
gap. These changes correct those cases; ranking weights remain provisional.

| Problem | Implemented behavior |
| --- | --- |
| Filtering after top-k and ingredient-core aggregation removed valid cold variants | Registry constraints apply to individual precedents before retrieval stopping, bounded neighbor scoring, aggregation, and ranking. The final response retains a defensive constraint check. |
| A hot stage passed a cooler summary-temperature limit | Every reported reaction stage is checked. An unknown stage temperature cannot establish compliance with a hard maximum. Atmosphere constraints use gas tokens and preserve unknown/mixed reports. |
| One common recipe exhausted the global 50-reference budget | All eligible recipe cores are grouped before taking up to 50 independent scoring neighbors per core. Total support remains separate from that sampling limit. |
| Similarity ties selected the highest yield from a publication | Representative evidence no longer uses yield as a tie-break. Historical outcomes are averaged within references and then equally across references, restricted to the selected full recipe variant. |
| Duplicate amounts changed with input order | Registry quantity rules normalize compatible units and preserve sorted source observations. Conflicting amounts remain unresolved; only explicitly identified distinct additions to the same stage may be summed. |
| Missing conditions received a perfect compatibility score | Assessments distinguish `unknown`, `no_known_conflict`, and `conflict`. Unknown coverage has no affirmative compatibility contribution. A narrow edit-based reduction capability check records missing hydrogen-source evidence without assuming a mechanism or forcing a named family. |

## Contract and API impact

- Resolved components and recipes: `2.2`; quantity definitions:
  `quantity_normalization.v1@1.0`. Recipe identifiers change because the quantity
  definition is part of canonical identity.
- Recommendation records/converter: `10.3` / `generic_conversion.v10.3`.
- Persisted generic index: `6.5`; compatibility definition: `1.3`;
  ranking definition: `1.2`.
- Generic result: `4.0`; weak-label result: `2.0`.
- `historical_yield_pct` replaces `expected_yield_pct` in public results and all
  application consumers. Generic results also expose `historical_yield_summary`
  with observed range, observation count, independent evidence count, method,
  selected recipe variant, and `is_prediction: false`.
- No route names change. API clients consuming the old yield field must migrate.

The observed range is not a confidence interval. Published yields are selected
historical evidence, and the summary does not estimate success probability for
the submitted query. Existing weights are explicitly marked
`prior_weights_pending_current_pipeline_validation`.

## Canonical artifact refresh

`python -m condition_recommender.refresh_artifacts_cli SOURCE_LIBRARY OUTPUT
--workers 8` refreshes inventoried canonical shards in a separate directory.
It verifies source-shard hashes and structural contracts, restores typed graph
observations, reruns the shared taxonomy annotation path, and invokes the normal
converter. It does not invent new atom correspondence or maintain a second
conversion implementation. Fresh-conversion parity regressions cover mapped
reactions, multiple bond types, ambiguity, and invalid maps.

The refresh rejects incompatible structural definitions and externally enriched
mapping provenance. Those inputs require ordinary conversion with their original
mapping provider. Observation IDs and existing signature IDs must remain stable.
Checkpoints bind source hashes, definitions, the refresh algorithm version, and
output hashes. Conversion validation checks all rows, source coverage, duplicate
identities, definitions, and checksums before index building.

Local validation evidence is under
`results/condition_recommendation_improvements/`. Large generated datasets and
indexes are local artifacts and must not be committed.

## Release limits

The initial deterministic baseline was 1,323 passing tests. The six reproduced
failures are development regressions, not a held-out accuracy benchmark.
The local rebuild makes the current application artifacts usable with current
contracts. Independent chemist adjudication and the untouched evaluation remain
required before declaring a production release or calibrated recommendation
accuracy. No untouched labels are used to tune these fixes.

The source-balanced development diagnostic contains 1,427 records from 119
sources, of which 711 enter the trusted development index. Its connected
reference/reaction split contains 593 training and 118 test observations, with
zero reference or canonical-reaction overlap. Hybrid and generic-only retrieval
both cover 58/118 queries; top-five recovery is 5/8 for recipes represented in
training. Family-only retrieval abstains on all queries. All three runs report
zero rule-detected hard incompatibilities. The small seen-recipe denominator,
mean top-result support of about one independent reference, and incomplete
condition evidence preclude a calibrated accuracy claim. Scaffold and source
overlap remain in this grouped diagnostic; it is not a disjoint release test.

An 80-case development blind packet contains 162 candidates, including 80
controls, with drawings and a separate answer key. Its blank review form is for
an independent chemist; no approval or adjudication has been filled in. This
packet does not replace the roadmap's complete stratified release review.
