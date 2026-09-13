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

Profiling the index rebuild also found redundant reactivity annotation in
`departing_fragment_tokens`. It now requests structural parsing explicitly.
The extractor consumes graph connectivity and input component SMILES only;
annotation-parity regressions preserve fragment identities while preventing
reactivity annotation from re-entering this indexing and retrieval path. This
performance correction does not change schemas or definition identities.
On 64 development records, both paths produced identical fragment tokens; the
structural-only path took 0.253 seconds versus 1.876 seconds (7.42 times faster
for this extraction step, not an end-to-end recommendation benchmark).

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
identities, definitions, and checksums before activation. The local full-corpus
run uses the same validation, catalog, and index builders concurrently after all
immutable shards exist; no runtime artifact is activated until every check passes.

Local validation evidence is under
`results/condition_recommendation_improvements/`. Large generated datasets and
indexes are local artifacts and must not be committed.

The Full rebuild contains 660,190 canonical records across 725 shards and 119
sources, with 571,157 trusted indexed precedents. All source observations and
existing signature identities are preserved, with no changes to the frozen
Full admission tiers, chemistry classifications, condition status, or index
eligibility. Canonical validation finds no duplicate observations or integrity
issues.

The refreshed Compact library preserves all 117,232 observation identities and
uses their frozen Full structural observations. Its trusted index contains
94,643 precedents, compared with 94,638 previously: 11 additions and six
removals. The 11 additions already had verified global atom correspondence in
the frozen Full artifact, while the older Compact artifact was unresolved.
Five removals have unresolved transformation evidence in frozen Full; the sixth
has unaccounted product atoms and a suspected missing reactant. All 17 refreshed
observations exactly match their frozen Full observations. This consolidates
pre-existing artifact disagreements; it does not establish that either older
atom-correspondence result was chemically correct. The removed observations
remain in canonical records, and the old Compact artifact is retained for
review. Per-record evidence is in `compact_change_provenance.json`.

Both rebuilt libraries passed deep SQLite record and lookup validation with no
issues and are installed at `datasets/literature/full` and
`datasets/literature/compact`. Previous artifacts are preserved under
`datasets/literature/_backup_condition_refresh_20260912`. Installation verifies
the staged report and index hashes, backs up prior files, and verifies copied
bytes before replacing each destination atomically. `activation_report.json`
records installed paths, counts, index identities, and file hashes.

API integration checks passed against both staged and default installed paths:
health and frontend routes return HTTP 200, and normal and unrestricted
recommendation requests return nonempty schema `4.0` results with historical
yield fields. These positive integration checks establish artifact and API
compatibility; they are not an accuracy or latency benchmark.

## Release limits

The initial deterministic baseline was 1,323 passing tests. The six reproduced
failures are development regressions, not a held-out accuracy benchmark.
The corrected code passes all 1,357 tests, including 34 added regressions. The
frontend build, Python static checks, and validation of 27,432 registry records
also pass. Code hashes bind these checks to `code_freeze.json`.
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
