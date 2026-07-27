# Generic Condition Recommendation Implementation Plan

**Status:** In progress
**Primary entry point:** `condition_recommender.recommend_generic_conditions()`  
**Source corpus:** `data-processor/reaction_dataset/`  
**Scope:** sample-first conversion, reference-aware precedent modeling, generic
retrieval, recipe ranking, evaluation, and eventual full-corpus conversion

## 1. Purpose

Implement a chemistry-first, type-agnostic condition recommendation workflow
over the structure-rich reaction datasets in `data-processor/reaction_dataset/`.

The system must recommend conditions from observed molecular transformations,
not from source reaction names or dataset filenames. A named reaction family is
an optional interpretation with confidence. It may improve retrieval when
strongly supported, but it is never a mandatory routing key.

The implementation extends the existing generic pilot rather than creating a
second recommendation path. The intended workflow is:

```text
source CSVs
    |
    v
deterministic sampling and source audit
    |
    v
canonical reaction and condition conversion
    |
    v
reference-aware condition-series modeling
    |
    v
versioned generic precedent index
    |
    v
chemistry compatibility filtering
    |
    v
structural retrieval and recipe aggregation
    |
    v
chemist-readable recommendation or explicit abstention
```

This plan follows:

- `docs/type_agnostic_reaction_recommendation_implementation.md`;
- `docs/phased_implementation_and_evaluation_plan.md`;
- `condition_recommender/README.md`;
- the package boundaries and current priorities in `AGENTS.md`.

## 2. Non-Negotiable Design Rules

1. The molecular graph is the source of truth.
2. Source reaction names, filenames, and references are provenance and evidence,
   not substitutes for structural chemistry.
3. Hard chemistry and recipe-compatibility checks run before similarity.
4. `named_family=None` is valid when mapped edits or exact reconstruction provide
   sufficient evidence.
5. Conflicting mapping and reconstruction evidence is retained for review and
   never resolved silently.
6. Missing condition values are reported as missing. They are not imputed as
   observed facts.
7. Canonical nested JSON or Parquet is the machine artifact. CSV is a review
   view.
8. A publication with many examples provides substrate-scope evidence, not many
   independent literature confirmations.
9. Development, validation, and test records must be reference-disjoint.
10. The full corpus is converted only after the sampled workflow and schemas
    pass machine and chemist review gates.

## 3. Current Starting Point

The repository already contains a functional pilot:

- generic record conversion in `condition_recommender/conversion/`;
- canonical `ResolvedConditionRecipe` support in `condition_registry`;
- generic signature indexing in `condition_recommender/generic_indexing.py`;
- hierarchical retrieval in `condition_recommender/generic_retrieval.py`;
- compatibility checks in `condition_recommender/compatibility.py`;
- orchestration and aggregation in `condition_recommender/generic_api.py`;
- grouped held-out evaluation in `condition_recommender/evaluation.py`.

The pilot is a useful baseline, but it is not yet the final production
implementation. The following gaps must be addressed before full conversion:

- finish the reaction-signature foundation and its human chemistry gate;
- include schema-level hydrogen changes in generic edit identity;
- separate chemistry admission from condition-data quality;
- retain unresolved condition identities with explicit uncertainty;
- prevent multi-stage conditions from being assigned to the wrong reaction
  event;
- normalize references and identify multiple condition series within one
  publication;
- count independent reactions and references rather than raw rows;
- make the environment retrieval tier reachable and meaningful;
- add reference- and scaffold-safe evaluation;
- preload or cache the runtime index rather than rebuilding it per request;
- calibrate retrieval and ranking by transformation class.

### 3.1 Implemented slices

The initial implementation slice now provides:

- reaction-signature schema 1.4, with schema-level hydrogen changes included in
  event, L2, and L3 identity plus a typed product-atom completeness assessment;
- deterministic DOI, patent, and normalized bibliographic reference identity;
- reference provenance in canonical converted records and persisted indices;
- deterministic source-balanced smoke, development, validation, and
  untouched-test sampling;
- connected partition groups for shared references and exact normalized
  reaction text;
- original source provenance in sampled CSV conversion;
- reference- and canonical-reaction-connected held-out evaluation groups.

The second implementation slice now also provides:

- independent chemistry, condition, outcome, and index-eligibility statuses,
  with the old admission tier retained as a compatibility summary;
- recipe-core identity that is stable across amounts and operating-parameter
  variants;
- deterministic reference-local condition-series identity using compatible
  generic chemistry, recipe core, stage/step structure, temperature, and time;
- quarantine of unstructured multi-stage records and solvent-only recipes;
- index admission based on explicit eligibility rather than the old tier;
- retention of chemistry/condition precedents with missing outcomes while
  excluding those observations from yield modeling;
- aggregation by recipe core with observed variants, canonical-reaction,
  condition-series, reference, dataset, and observation support reported
  separately;
- publication-level deduplication for similarity and outcome aggregation.

The third implementation slice now also provides:

- a definition-ordered retrieval ladder with validated configuration;
- local-environment neighbor narrowing inside the hard-compatible bond-edit
  pool, before the broad bond-edit fallback;
- a persisted inverted environment-feature map in generic index schema 1.4, so
  neighbor selection does not scan unrelated bond-edit precedents;
- minimum-support decisions based on independent references, or canonical
  reactions when a reference is unavailable;
- a typed per-tier retrieval trace containing row, independent-support,
  compatibility, exclusion, and selection counts;
- a reusable `GenericConditionRecommender` that loads and validates a persisted
  index once for repeated runtime queries;
- evaluation-case serialization of independent counts and retrieval traces.

The fourth implementation slice now also provides:

- separate, validated `generic_similarity.v1` and `generic_ranking.v1`
  definitions, leaving retrieval configuration responsible only for pool
  selection;
- focused similarity and recipe-ranking modules with no aggregation logic in
  the application-facing orchestration layer;
- a complete score trace containing structural feature values, weighted
  contributions, ranking components, applied weights, definition versions,
  independent evidence counts, and usable outcome counts;
- reference-aware evidence deduplication plus saturating support, reaction
  breadth, and dataset-diversity terms, so repeated examples from one
  publication do not act as independent confirmations;
- explicit omission and weight renormalization when a recipe has no usable
  yield evidence;
- a centralized registry-owned recipe-component vocabulary used by conversion
  and compatibility validation;
- compatibility definition schema 1.1, validated against registry families,
  component buckets, and taxonomy spectator tags;
- an added acid/hydrolysis caution while preserving hard chemistry exclusions
  before similarity.

The fifth implementation slice now also provides:

- generic index schema 1.5 with publication year and molecularly derived
  reactive-partner scaffold tokens;
- grouped-random, strict scaffold-disjoint, source-corpus-disjoint, and
  forward-time evaluation modes;
- family-only, generic-only, hybrid, transformation-prior, and legacy-pilot
  baselines behind the same chemistry and compatibility gates;
- metrics for abstention, variants, outcomes, uncertainty, explanation
  completeness, independent support, transformation class, and evidence;
- development-only calibration with a separate validation promotion gate;
- the promoted `support_heavy_min2` ranking profile and minimum independent
  support of two;
- blind review packets with a chemist-readable HTML view, highlighted
  structures, randomized candidates, unmarked unsuitable controls, blank
  review fields, and a separate answer key;
- post-review adjudication that refuses premature unblinding, checks exact
  form/key candidate coverage, reports recommendation/control acceptance, and
  binds independent reviewer sign-off to artifact hashes.

The Phase 6 implementation foundation now also provides restartable parallel
canonical JSONL-gzip shards, periodic manifest checkpoints, deterministic
recipe/reference/series catalogs, source and artifact checksums, definition
contracts, and conversion/index/release integrity validators.

Structured multi-stage assignment and reference condition inference remain
future conversion-quality improvements. Untouched-test and full-corpus
execution remain gated by independent review of the generated chemist packet.

The first local pilot artifacts use the documented default sizes: 500 smoke,
5,000 development, 2,000 validation, and 2,000 untouched-test rows. Reference
and exact-reaction-text leakage are both zero. The balanced 500-row smoke
conversion under the new dimensional policy yields 22 legacy-tier verified,
280 review, and 198 rejected records. It identifies 76 chemistry-verified rows,
145 rows with signatures, and 24 explicitly index-eligible precedents. The
balanced sample deliberately contains 355 multi-stage rows, all now quarantined
from automatic indexing; 11 missing and one invalid outcome are tracked
separately from chemistry and condition usability. These are diagnostic
coverage measurements, not production performance claims. The first smoke run
also exposed and fixed a mixed-edit sorting defect in the conserved-scaffold
fallback, now covered by a deterministic regression.

The Phase 4 grouped smoke diagnostic used the 24 index-eligible precedents with
19 train rows and 5 test rows across disjoint connected groups. It produced
zero reference/reaction-group leakage, covered all five queries, and returned
zero hard-incompatible recommendations. None of the five held-out recipes
occurred in training, so the observed zero top-k recipe recovery is not a
ranking comparison; the 43.34 percentage-point yield MAE is likewise a
small-sample diagnostic. At that slice the ranking definition was marked
`uncalibrated_pilot`; it was superseded only by the later development/validation
calibration described below.

The frozen sharded sample artifacts contain:

- development: 5,000 converted rows, 534 index-eligible precedents, 347 recipe
  cores, 390 references, and zero failed shards;
- validation: 2,000 converted rows, 251 index-eligible precedents, 186 recipe
  cores, 177 references, and zero failed shards.

The selected profile improved mean development seen-recipe top-1 recovery from
approximately 11.0% to 12.2%. On validation, selected and prior profiles both
achieved 25% seen-recipe top-1 and 75% top-5 recovery, full coverage, and zero
hard-incompatible recommendations. Strict scaffold-, source-, and time-disjoint
validation also returned zero hard-incompatible recommendations. These remain
sample-level measurements rather than broad production-accuracy claims.

The blind packet at `results/generic_evaluation/v2/chemist_review/` contains
111 cases, 428 randomized candidates, and one unsuitable control per case.
Review `review_packet.html` and record decisions in `review_form.csv` without
opening `answer_key.jsonl`. Once every decision is recorded, the adjudication
command unblinds the packet and emits a provenance-bound `review_summary.json`.
Untouched-test evaluation and full conversion remain gated until an independent
chemist signs that summary with no unresolved systematic defect.

## 4. Source Corpus Profile

The current corpus contains 118 CSV files and 238,214 rows.

| Field | Non-empty rows | Coverage |
| --- | ---: | ---: |
| Reaction SMILES | 238,213 | approximately 100% |
| Yield | 235,077 | 98.68% |
| Reagent CAS | 199,295 | 83.66% |
| Catalyst CAS | 145,699 | 61.16% |
| Solvent CAS | 223,618 | 93.87% |
| Experimental procedure | 19,193 | 8.06% |
| Time | 1,102 | 0.46% |
| Temperature | 212 | 0.09% |

There are 83,032 records with more than one reported stage, or 34.86% of the
corpus. These records are not automatically safe precedents because the source
condition lists may combine multiple operations.

The first useful release should therefore recommend a condition-component
regime:

- catalysts and ligands;
- bases, acids, oxidants, reductants, activators, and additives;
- solvents;
- any reliably reported operating parameters.

It must not present temperature, time, concentration, atmosphere, catalyst
loading, or equivalents as known when the source did not report them.

## 5. Sample-First Development Strategy

### 5.1 Why sampling comes first

Full conversion requires molecular parsing, reaction analysis, registry
resolution, serialization, and indexing for more than 238,000 rows. Repeating
that work while contracts are changing would be slow and would create stale
artifacts.

Create deterministic pilot partitions before changing recommendation logic:

| Partition | Default target | Intended use |
| --- | ---: | --- |
| Smoke | 500 rows | Fast conversion and index debugging |
| Development | 5,000 rows | Implement chemistry, conversion, and retrieval |
| Validation | 2,000 rows | Tune thresholds and weights |
| Untouched test | 2,000 rows | Final machine evaluation |
| Chemist review | 300 rows | Blind review of outputs and abstentions |

These are configurable defaults, not schema constants.

### 5.2 Sampling unit and partition ownership

The normalized publication reference is the primary partition owner. All
selected rows from the same reference must remain in one partition, including
when that reference appears in multiple source CSVs.

Canonical duplicate reactions must also remain in one partition. Where
reference and canonical-reaction grouping conflict, construct connected groups
so neither identity crosses partitions.

Scaffold overlap is audited after chemistry conversion. The untouched test set
must report cross-partition substrate-scaffold overlap and should use a stricter
scaffold-disjoint variant for final calibration.

### 5.3 Sampling strata

The initial metadata sampler should balance:

- all 118 source files;
- large and small source categories;
- singleton and repeated references;
- references with one raw recipe and multiple raw recipes;
- references appearing in multiple files;
- single-stage and multi-stage records;
- complete, partial, and missing condition source fields;
- yield bins and missing yield;
- mapped and unmapped reaction SMILES;
- records with and without procedure text.

After the first chemistry conversion, create a second deterministic
stratification pass covering:

- exact product reconstruction;
- exact multi-event reconstruction;
- valid supplied mapping;
- unknown-family verified signatures;
- grammar-only candidates;
- conflicting edit evidence;
- invalid maps and invalid products;
- single- and multi-event signatures;
- inter-, intra-, mixed-, and unimolecular topology;
- named-family confidence ranges;
- condition identity and role-confidence ranges.

Balanced samples are used to find failures. A separate prevalence-weighted
evaluation sample is required to estimate expected corpus behavior.

### 5.4 Sampling artifacts

Add:

```text
condition_recommender/conversion/sampling.py
condition_recommender/sample_cli.py
```

Generate:

```text
datasets/generic_pilot/v1/
  sample_manifest.v1.json
  smoke.csv
  development.csv
  validation.csv
  untouched_test.csv
  chemist_review.csv
  sampling_report.json
```

Each manifest entry retains:

- source path and source dataset;
- original row number and reaction ID;
- raw and normalized reference identity;
- canonical reaction ID when available;
- assigned partition;
- sampling strata;
- deterministic selection hash;
- sampler schema, version, and seed.

If source records cannot be committed, commit only stable IDs, hashes, and
sampling metadata. Keep extracted source rows local.

## 6. Reference-Aware Precedent Modeling

### 6.1 Corpus behavior

Reference grouping is important in this corpus:

- 234,322 rows contain a reference;
- 42,336 normalized references are represented;
- 25,789 references contain multiple rows;
- 14,388 repeated references currently show one raw condition signature;
- 11,401 repeated references show multiple raw condition signatures;
- 5,566 references occur in more than one source CSV.

Consequently, a reference is not equivalent to a recipe. One publication may
contain an optimization table, several substrate-scope regimes, control
experiments, or multiple reaction classes.

### 6.2 Reference identity

Add an immutable contract:

```python
@dataclass(frozen=True)
class ReferenceIdentity:
    reference_id: str
    raw_reference: str
    normalized_citation: str
    doi: str | None
    patent_number: str | None
    publication_year: int | None
    resolution_status: str
    warnings: tuple[str, ...]
    schema_version: str
```

Resolve reference identity in this order:

1. normalized DOI;
2. normalized patent publication number;
3. normalized journal, year, volume, and page or article number;
4. deterministic normalized-citation hash;
5. deterministic raw-reference hash with an uncertainty warning.

Never discard or overwrite the original source citation.

### 6.3 Condition series inside a reference

Define:

```text
reference_condition_series_id =
    reference_id
    + canonical recipe core
    + compatible stage/procedure identity
```

Rows belong to the same series only when:

- resolved recipe components and contextual roles agree;
- stage structures are compatible;
- reported operating parameters do not conflict;
- the reactions represent compatible chemistry;
- no procedure evidence distinguishes separate regimes.

Preserve legitimate alternative recipes as separate series.

Add:

```python
@dataclass(frozen=True)
class ReferenceConditionSeries:
    series_id: str
    reference_id: str
    recipe_core_id: str
    recipe_variant_ids: tuple[str, ...]
    observation_ids: tuple[str, ...]
    canonical_reaction_ids: tuple[str, ...]
    transformation_classes: tuple[str, ...]
    condition_consistency: str
    evidence: tuple[str, ...]
    warnings: tuple[str, ...]
    schema_version: str
```

### 6.4 Missing conditions and reference inference

Reference-related conditions must never overwrite missing source values.

When a row has missing conditions, the converter may record separate,
non-observed candidates only if:

- the reference contains one dominant condition series;
- the row has compatible structural chemistry;
- no competing series applies;
- stage assignment is unambiguous.

The derived field must retain:

- candidate recipe ID;
- confidence;
- supporting observation IDs;
- inference method;
- warnings such as `CONDITIONS_INFERRED_FROM_REFERENCE_SERIES`.

Reference-inferred recipes are lower-quality evidence and are not equal-weight
observed precedents. Multi-series or multi-stage ambiguity prevents automatic
inference.

## 7. Canonical Conversion Contract

### 7.1 Separate quality dimensions

Replace the single overloaded admission tier with explicit dimensions:

```text
chemistry_status:
  verified | review | rejected

condition_status:
  resolved_complete
  resolved_partial
  unresolved_retained
  multistage_ambiguous
  unusable

outcome_status:
  usable | missing | invalid

index_eligibility:
  eligible | review_only | ineligible
```

Important behavior:

- missing `named_family` does not reduce verified chemistry;
- exact reconstruction and valid supplied mapping can verify chemistry;
- conflicting edit evidence remains review-only;
- unresolved identities remain searchable with uncertainty when a raw recipe is
  usable;
- missing yield does not erase useful chemistry and condition observations, but
  excludes the row from yield modeling;
- multi-stage ambiguity is review-only until conditions are assigned to the
  reaction-forming stage;
- solvent-only records are not complete reaction recipes;
- source reaction labels remain provenance.

### 7.2 Canonical record

Evolve `RecommendationRecord` to include:

- raw source identity and row number;
- `ReferenceIdentity`;
- reference condition-series identity;
- raw and canonical reaction identity;
- complete versioned `ReactionSignature`;
- observation, interpretation, and conflict evidence;
- independent chemistry, condition, outcome, and index statuses;
- optional named-family candidates with confidence and evidence;
- canonical resolved recipe;
- recipe completeness and uncertainty;
- stage and procedure provenance;
- source yield and other reported outcomes;
- taxonomy, registry, converter, reference, and schema versions.

### 7.3 Recipe identity

Use two related identities:

```text
recipe_core_id
    resolved substances + contextual roles

recipe_variant_id / RCR1
    recipe core + known amounts + known operating parameters + stages
```

This prevents sparse temperature and time fields from fragmenting otherwise
equivalent condition regimes while retaining exact procedure variants when
reported.

Every component retains:

- raw identifier;
- resolved substance ID and canonical name;
- source field;
- contextual roles and confidence;
- amount and unit when available;
- provenance;
- identity and role warnings.

### 7.4 Canonical artifacts

For pilot development, write:

```text
results/generic_conversion/pilot_v1/
  records.jsonl
  recipe_catalog.jsonl
  reference_catalog.jsonl
  reference_condition_series.jsonl
  verified.csv
  review.csv
  rejected.csv
  multistage_review.csv
  unresolved_conditions.csv
  conversion_report.json
  conversion_report.md
```

For the full corpus, use partitioned Parquet as the primary record format and
retain JSONL or CSV only for manageable review subsets.

## 8. Generic Index

The persisted index must contain or reference:

- exact, handle, transformation, edit, and environment signature keys;
- complete event and topology information;
- named-family candidates and confidence;
- spectator groups and compatibility tags;
- recipe core and variant IDs;
- observation, canonical reaction, condition-series, and reference IDs;
- source-corpus identity;
- outcomes and evidence quality;
- all schema and definition versions.

Index admission uses `index_eligibility`, not the old combined tier.

Minimum support is measured using independent evidence units:

- unique canonical reactions;
- unique substrate or scaffold combinations;
- unique reference condition series;
- unique references;
- unique source corpora.

Raw row count is reported but does not determine adequate support.

## 9. Retrieval and Compatibility Logic

### 9.1 Retrieval ladder

Select the narrowest chemistry-compatible tier with adequate independent
support:

1. exact L0 reaction signature;
2. L1 handles, contexts, events, and topology;
3. high-confidence named family intersected with compatible edits;
4. L2 generic transformation and topology;
5. L3 compatible edit set ranked or narrowed by L4 environment similarity;
6. a validated broad condition-regime prior with a strong caution;
7. typed abstention.

The environment tier must not follow an already broader adequate bond-edit
pool. Environment similarity is used to select neighbors inside the compatible
edit pool.

Family evidence never crosses an incompatible bond-edit gate. Low-confidence
family evidence contributes only to scoring.

### 9.2 Hard compatibility before scoring

Hard exclusions include:

- incompatible complete net edits;
- incompatible multi-event sets;
- incompatible valence or reactive-handle class;
- essential topology conflicts;
- oxidant/reductant conflicts with unchanged sensitive groups;
- acid, moisture, or atmosphere conflicts when evidence is sufficient;
- certainly missing mandatory regime components;
- conflicting or invalid transformation evidence.

Soft penalties and cautions include:

- catalyst poisoning or strong coordination;
- acid/base interactions;
- hydrolysis risk;
- known temperature risk;
- incomplete recipes;
- identities or roles too uncertain for a hard decision;
- topology relaxation at a documented fallback.

Unknown identity means uncertain presence or role. It must not be interpreted as
certain absence.

### 9.3 Similarity

Score only after hard filtering. Use interpretable features:

- complete edit and event similarity;
- reaction topology and ring-formation similarity;
- reactive handles and attachment contexts;
- local steric and electronic environments;
- nearby functional groups;
- unchanged spectator groups;
- condition-regime compatibility;
- confidence-weighted family agreement;
- precedent chemistry and condition-data quality.

Missing features never count as matches. Their treatment must be explicit in
the score trace.

Weights remain in a versioned definition and are calibrated on the development
and validation partitions. Source reaction labels, filenames, and references
do not participate in structural similarity.

## 10. Reference-Aware Recipe Aggregation

Group retrieved precedents first by `recipe_core_id`, then expose reported
variants.

For every recipe report:

- raw observation count;
- unique canonical reaction count;
- unique scaffold-pair count;
- reference condition-series count;
- independent reference count;
- independent source-corpus count;
- yield count, mean, spread, and uncertainty;
- recipe completeness distribution;
- identity and role uncertainty.

Multiple rows from one publication receive diminishing support. A large
substrate scope demonstrates breadth within that publication but does not count
as independent replication.

Ranking combines:

- structural similarity;
- hard-filtered compatibility score;
- independent reference support;
- within-reference substrate breadth with saturation;
- dataset and evidence quality;
- recipe completeness;
- cautiously weighted, shrinkage-adjusted outcome evidence;
- penalties for uncertainty and mismatches.

Yield is secondary because the corpus is dominated by published successful
reactions and may contain strong selection bias.

## 11. Public API and Explanation Contract

Keep `recommend_generic_conditions()` as the public orchestration entry point.
Move similarity, aggregation, and trace construction into focused modules, for
example:

```text
condition_recommender/
  generic_api.py
  generic_retrieval.py
  compatibility.py
  similarity.py
  recipe_ranking.py
  retrieval_trace.py
```

The API requires a complete product-specified reaction SMILES. It must:

1. validate and featurize the reaction;
2. reject or abstain on conflicting or insufficient query evidence;
3. validate index and definition versions;
4. retrieve through the explicit ladder;
5. apply hard recipe compatibility;
6. score precedents;
7. aggregate canonical recipe cores and variants;
8. return recommendations or a typed abstention.

Every result reports:

- query signature and transformation;
- optional family and confidence;
- retrieval level;
- candidate counts at every attempted tier;
- exclusion reasons;
- matching edits, events, handles, topology, and environments;
- important mismatches;
- compatibility evidence;
- recipe components and missing parameters;
- raw, reaction, series, reference, and corpus support;
- representative precedent reaction and reference IDs;
- uncertainty and cautions;
- artifact and definition versions.

If support is concentrated in one publication, include a caution such as:

```text
Most supporting examples come from one publication and one experimental
condition series; independent reproducibility is not established.
```

Runtime applications should preload a validated persisted index. They should
not rebuild the JSONL index on every query.

## 12. Evaluation

### 12.1 Leakage control

Evaluation partitions must be disjoint by:

- normalized reference ID;
- canonical reaction ID;
- source-corpus identity.

Also report substrate-scaffold leakage. The final benchmark should include a
strict scaffold-disjoint evaluation and, when publication year is available, a
forward time split.

### 12.2 Baselines

Compare:

1. family-only retrieval;
2. generic-only retrieval;
3. hybrid family-to-generic retrieval;
4. a simple transformation-level recipe prior;
5. the existing generic pilot.

### 12.3 Metrics

Report globally and by transformation and evidence quality:

- query and recommendation coverage;
- typed abstention rate and correctness;
- top-1 and top-k recipe-core recovery;
- variant recovery when the exact variant exists in training;
- recovery conditional on recipe availability;
- fallback-level distribution;
- independent-reference support distribution;
- yield MAE and calibration where yield is usable;
- compatibility exclusions;
- hard-incompatible recommendation count;
- identity and recipe-completeness uncertainty;
- explanation-field completeness.

### 12.4 Chemist review

Generate a blind review packet containing:

- query structures and highlighted reaction centers;
- concise observed transformations;
- randomized candidate recipe order;
- representative precedents and citations;
- matching and mismatching chemistry;
- compatibility warnings;
- deliberately unsuitable negative-control recipes;
- blank reviewer decision and comment fields.

Review accepted, rejected, ambiguous, low-support, and broad-fallback cases.
Confirmed disagreements become regression tests or definition changes.

## 13. Full-Corpus Conversion

Begin the full conversion only after:

- signature and converter schemas are frozen;
- sample conversion snapshots are reviewed;
- reference normalization and series grouping are validated;
- retrieval and compatibility tests pass;
- the untouched evaluation is complete;
- the chemist review has no unresolved systematic defect.

The full conversion must be deterministic, restartable, and sharded:

```text
results/generic_conversion/v1/
  shards/
    suzuki_miyaura.part-000.parquet
    C_N_Coupling.part-000.parquet
    ...
  recipe_catalog.parquet
  reference_catalog.parquet
  reference_condition_series.parquet
  shard_manifest.json
  conversion_report.json
  conversion_report.md
```

Each shard records:

- source path and checksum;
- source row range;
- conversion status;
- input and output row counts;
- output checksum;
- converter and record schemas;
- taxonomy and registry definition versions;
- sampling or full-corpus mode;
- warnings and failure counts.

Reuse only shards whose source checksum and every participating schema and
definition version still match.

Cache deterministic work by:

- canonical reaction identity for reaction analysis;
- raw condition signature plus reaction context for role resolution;
- normalized citation text for reference resolution.

Parallel processing may operate on independent shards. Merge ordering and
artifact IDs must remain deterministic.

## 14. Implementation Phases

### Phase 0: Foundation gate

- complete reaction-edit human review;
- correct hydrogen-change signature identity;
- verify L0-L4 determinism and partner-order invariance;
- preserve mapped unknown-family signatures;
- pass the complete test suite.

**Gate:** Phase A signature requirements and parity tests pass.

### Phase 1: Sampler and reference identity

- implement deterministic reference-safe sampling;
- add `ReferenceIdentity`;
- generate smoke, development, validation, untouched-test, and review manifests;
- add sampling and reference normalization tests.

**Gate:** repeated runs produce identical partitions with zero reference
leakage.

### Phase 2: Canonical conversion

- separate chemistry, condition, outcome, and index statuses;
- add recipe core identity and completeness;
- add reference condition-series modeling;
- quarantine stage-ambiguous records;
- write canonical catalogs and review views.

**Gate:** sampled conversion is deterministic and chemist review explains every
admission class.

### Phase 3: Generic index and retrieval

- extend the index with reference and series identities;
- correct the fallback ladder;
- count independent support;
- make environment-neighbor retrieval reachable;
- add a preloaded-index runtime path.

**Gate:** all retrieval tiers, support thresholds, and stale-artifact failures
have deterministic tests.

### Phase 4: Compatibility, similarity, and aggregation

- expand declarative compatibility rules;
- extract similarity and recipe aggregation modules;
- add reference-aware diminishing support;
- return recipe cores plus reported variants;
- implement complete retrieval traces.

**Gate:** no hard-incompatible recipe is returned in the curated and validation
benchmarks.

### Phase 5: Evaluation and calibration

- run reference- and canonical-reaction-disjoint evaluation;
- add scaffold and time-split reports;
- compare family-only, generic-only, and hybrid baselines;
- calibrate thresholds and weights;
- conduct blind chemist review.

**Gate:** hybrid retrieval expands justified coverage without increasing hard
chemistry violations or materially degrading established-family performance.

### Phase 6: Full conversion and release

- freeze versions;
- convert the corpus in restartable shards;
- build and validate the production index;
- run full-corpus coverage and integrity reports;
- document supported and abstaining chemistry scope.

**Gate:** index integrity, artifact provenance, performance, and final regression
suite pass.

## 15. Initial Chemistry Scope

The first calibrated pilot should prioritize the required parity chemistry:

- Suzuki C-C coupling;
- C-N coupling;
- C-O coupling;
- C-S coupling;
- mapped unknown-family reactions.

The next wave may include already supported graph operators:

- amide and other acyl N/O/S formation;
- reductive amination;
- alkyl C-N/O/S substitution;
- hydrogenation and other validated redox edits;
- exactly reconstructed Heck reactions.

Wittig, aldol, Minisci, Grignard, cycloaddition, protection, cascade, and
rearrangement datasets remain review-only until their partner-site, operator,
and multi-edit contracts are implemented and tested. A matching source filename
does not establish support.

## 16. Testing Requirements

Add deterministic tests for:

- sampling stability and reference partition ownership;
- DOI, patent, journal, and unresolved reference normalization;
- one-reference/one-series and one-reference/multiple-series cases;
- cross-file references;
- missing-condition reference inference and prohibited inference;
- multi-stage quarantine;
- independent chemistry and condition statuses;
- unresolved condition retention;
- recipe core and variant identity;
- hydrogen changes in edit identity;
- reachable environment retrieval;
- support counted by reaction, series, and reference;
- duplicate-row resistance;
- multi-event retrieval;
- hard compatibility exclusions before similarity;
- missing features not treated as matches;
- reference-aware recipe aggregation;
- typed abstentions and stale-index errors;
- reference-, reaction-, and scaffold-leakage reports.

Required chemistry regressions remain:

- Suzuki;
- C-N, C-O, and C-S;
- positive, negative, ambiguous, and conflicting evidence;
- partner-order invariance;
- mapped unknown-family reactions;
- invalid maps;
- deterministic signature and recipe IDs.

Run:

```powershell
pytest -q
pytest -q tests/reactive_taxonomy
pytest -q tests/condition_registry
pytest -q tests/condition_recommender
```

Dataset snapshot changes require chemistry review. Do not accept a snapshot
solely because implementation output changed.

## 17. Recommended Change Sequence

Implement as small, reviewable changes:

1. `fix(reactive-taxonomy): complete generic signature identity`
2. `feat(condition-recommender): add reference-safe pilot sampling`
3. `feat(condition-recommender): normalize reference identities`
4. `refactor(condition-recommender): separate admission dimensions`
5. `feat(condition-registry): add recipe core identity`
6. `feat(condition-recommender): model reference condition series`
7. `fix(condition-recommender): correct generic retrieval ladder`
8. `refactor(condition-recommender): extract similarity and ranking`
9. `test(condition-recommender): add leakage-safe evaluation`
10. `feat(condition-recommender): add sharded full conversion`

Do not tune final ranking weights before the conversion schema, support units,
and leakage-safe evaluation partitions are stable.

## 18. Definition of Done

The generic recommendation workflow is ready for full-corpus use when:

1. sampled conversion is deterministic and versioned;
2. references and within-reference condition series are modeled explicitly;
3. no reference or canonical reaction crosses evaluation partitions;
4. verified unknown-family reactions are indexable;
5. multi-stage ambiguity cannot silently contaminate recipes;
6. unresolved identities remain traceable with uncertainty;
7. retrieval uses the narrowest adequate chemistry-compatible tier;
8. hard compatibility runs before similarity;
9. support distinguishes observations, reactions, series, and references;
10. recommendations report missing fields rather than guessing them;
11. every result cites representative reactions and references;
12. unsupported chemistry produces a typed abstention;
13. held-out evaluation and blind chemist review pass;
14. the full conversion is restartable, sharded, and integrity checked;
15. the complete test suite passes.
