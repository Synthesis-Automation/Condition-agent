# Condition Recommender

## Source-data preprocessing

Heterogeneous raw CSV files can first be normalized into the versioned,
chemistry-free `source_observation.v1` contract. Each selected source file
produces one independently reusable `<source>.observations.jsonl.gz` artifact
without a separate JSON log. Embedded source provenance, adapter version, and
intermediate schema metadata control cache reuse; the returned in-memory summary
still includes the output checksum and row statistics.

Launch the Qt application from the repository root:

```powershell
python -m app.data_preprocessor_gui
```

Or preprocess explicitly selected files from the command line:

```powershell
python -m condition_recommender.preprocess_cli `
  raw_dataset/HiTEA/8_SEPT_APPROVED_full_dataset.csv `
  --output-dir datasets/intermediate
```

For the cleaned USPTO dataset (automatic detection selects the USPTO adapter):

```powershell
python -m condition_recommender.preprocess_cli `
  raw_dataset/USPTO/USPTO_condition_reactions_cleaned.csv `
  --output-dir datasets/intermediate
```

Available adapters cover the literature CSV contract, the HiTEA approved CSV,
the cleaned USPTO condition-reaction CSV, and the original weak-label v2.1 CSV.
The USPTO adapter prefers `remapped_rxn` when it contains atom-map syntax and
complete reaction sides, otherwise retaining `canonical_rxn` as a clearly
flagged fallback. Automatic selection requires one exact, unambiguous header
match. Conditions are stored as source-faithful component claims, identifier
evidence, quantities, and ordered process stages. This stage does not parse
molecules, validate atom mapping, assign reaction families, resolve registry
substances, or make admission decisions.

The generic conversion engine accepts either raw CSV files or the generated
`*.observations.jsonl.gz` artifacts. Structure-backed observations continue
through normal chemistry and condition-registry conversion. Label-only records
are retained and condition-normalized but remain structurally ineligible for
the generic precedent index.

`condition_recommender` exposes one canonical structure-backed recommendation
path. It consumes reaction observations from `reactive_taxonomy` and canonical
condition identities and recipes from `condition_registry`. Former expert-rule
and weak-label recommenders have been removed; label-only source observations
may be preserved for audit but cannot become structural condition precedents.

## Recommendation paths and status

| Path | Public entry point | Current status | Intended use |
| --- | --- | --- | --- |
| Generic structure-backed retrieval | `recommend_generic_conditions()` | Functional pilot; coverage and calibration depend on converted structure-rich records | Type-agnostic precedent retrieval and canonical recipe aggregation |
| Exploratory reaction discovery | `ReactionDiscoveryExplorer.discover()` | Functional MVP over the validated SQLite precedent index | Multi-tier structural analogue search with observed-condition evidence |

Supporting audit, conversion, indexing, and evaluation commands prepare data for
this path; they do not recommend conditions by themselves.

## Reaction discovery mode

The desktop application includes **Reaction discovery** as an option in the
existing recommender window. It deliberately uses a separate result contract
and backend method from condition recommendation:

```python
from condition_recommender import ReactionDiscoveryExplorer

explorer = ReactionDiscoveryExplorer.from_path(
    "datasets/literature/generic_index.sqlite"
)
result = explorer.discover(
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
    top_k=20,
    view="diverse_strategies",
)
```

Discovery unions exact-signature, bond-edit, reaction-core,
local-environment, and anonymous edit-graph candidates instead of stopping at
the first sufficiently supported recommendation tier. Ranking is declared in
`definitions/discovery_retrieval.v1.json`. Bond-edit similarity has the largest
weight, followed by reaction-center transitions and local environment; partner
category, unchanged spectator groups, reaction topology, and reactive scaffold
provide additional evidence. Missing factors remain `null`, and the available
weights are explicitly renormalized in each `DiscoveryScoreTrace`.

Each hit is an individual reaction observation with its yield, evidence tier,
reference, structural matches and mismatches, and observed recipe. Low-yield
and unreported outcomes are included by default because they can identify
failure boundaries. They can be hidden with `include_low_yield=False` or
`include_unreported_outcomes=False`. For an ambiguous query, edit hypotheses
are searched and labeled separately rather than collapsed into a false
consensus.

Discovery conditions are evidence for expert ideation, **not recommended
conditions**. The MVP searches the same admitted SQLite precedent corpus as the
recommender; a future broader discovery-only corpus may add structurally
verified failed or incomplete records without granting them recommendation
authority. A complete `reactants>>product` query is currently required.

## Which path should I use?

- Use generic structure-backed retrieval when a compatible converted record set
  or persisted generic index is available.
- Abstain when structural evidence, compatible precedents, or independent
  support is insufficient. Review-qualified core, edit-hypothesis, external-map,
  partial-transformation, and structure-fallback modes remain inside the generic
  API and explicitly preserve their evidence level and cautions.

## Generic system: quick start

Run commands from the repository root with Python 3.10 or newer and RDKit
available. The generic system is composed of three standalone packages:

```text
reactive_taxonomy + condition_registry -> condition_recommender
```

`reactive_taxonomy` determines what changed in the molecular graph,
`condition_registry` resolves condition identities and canonical recipes, and
`condition_recommender` retrieves compatible precedents and ranks recipes.

For verified queries with a trustworthy reaction core, the recommender first
compares generic graph-derived reaction facets: normalized edits, active-atom
states, reaction topology, retained attachment classes, and unambiguous active
X-H site classes. It fills the requested `top_k` progressively with distinct
canonical recipe cores: exact structural-class recipes first, then recipes
whose retained attachment classes match through the versioned parent
hierarchy, and only then the broader signature/edit ladder. Each recommendation
retains its actual retrieval tier and relaxation caution.

### 1. Create a small local index

First create reference-safe samples, then convert only one 100-row shard for a
mechanical smoke test. The sampler prevents a smoke test from accidentally
consisting only of the first source file. This bounded index is deliberately
small and is not an accuracy benchmark:

```powershell
python -m condition_recommender.sample_cli `
  data-processor/reaction_dataset `
  results/quickstart/samples

python -m condition_recommender.sharded_conversion_cli `
  results/quickstart/samples/smoke.csv `
  results/quickstart/conversion `
  --mode smoke --shard-size 100 --max-shards 1

python -m condition_recommender.generic_index_cli `
  results/quickstart/conversion/records.jsonl.gz `
  results/quickstart/generic_index.sqlite

python -m condition_recommender.generic_index_integrity_cli `
  results/quickstart/generic_index.sqlite `
  --output-path results/quickstart/index_integrity.json
```

The integrity command exits with a nonzero status when the index is stale,
internally inconsistent, or contains duplicate observations. Because
`--max-shards 1` intentionally covers only part of `smoke.csv`, use this
artifact only to exercise the software path.

### 2. Request conditions

The CLI accepts reaction SMILES and returns JSON:

```powershell
python -m condition_recommender.generic_recommend_cli `
  "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1" `
  --records results/quickstart/generic_index.sqlite `
  --top-k 5
```

Chemists can choose a transparent ranking preset without changing chemistry
admission or compatibility gates. For example:

```powershell
python -m condition_recommender.generic_recommend_cli `
  "Brc1ccccc1.CCNCC>>CCN(CC)c1ccccc1" `
  --records results/quickstart/generic_index.sqlite `
  --ranking-profile reactant_category `
  --ranking-weight functional_group_tolerance=0.30
```

Available presets are `default`, `reactant_category`,
`functional_group_tolerance`, `evidence`, `yield`, and
`procedure_completeness`. Repeated `--ranking-weight COMPONENT=VALUE`
arguments override the chosen preset; weights are validated and normalized.

For repeated recommendations, load the validated index once:

```python
from condition_recommender import (
    ChemistRankingPreferences,
    GenericConditionRecommender,
)

recommender = GenericConditionRecommender.from_path(
    "results/quickstart/generic_index.sqlite"
)
result = recommender.recommend(
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
    top_k=5,
    ranking_preferences=ChemistRankingPreferences(
        profile_id="functional_group_tolerance"
    ),
)

if not result.valid:
    print(result.error)
else:
    print("retrieval level:", result.retrieval_level)
    for recommendation in result.recommendations:
        print(
            recommendation.rank,
            recommendation.recipe_core_id,
            recommendation.score,
            recommendation.resolved_recipe,
        )
```

Use `result.to_dict()` when a JSON-serializable audit payload is needed.
Applications should inspect:

- `valid`, `error`, and `warnings` before using any recommendation;
- `recommendation_mode`, which distinguishes verified-signature retrieval from
  the explicitly unverified structure fallback;
- `retrieval_level` and `retrieval_trace` to see which fallback was required;
- `resolved_recipe` for canonical condition components and operating variants;
- independent `reference_support` and `precedent_reference_ids`, not only raw
  observation counts;
- `compatibility_evidence`, `cautions`, `explanation`, and `score_trace`;
- `ranking_preferences`, `default_rank`, `rank_change`, and `factor_evidence`
  to audit how a chemist-selected priority changed the ordering;
- `expected_yield_pct`, which is `None` when no usable outcome evidence exists.

An empty recommendation list or typed error is a valid abstention. It must not
be replaced by a reaction-name guess.

### Transparent chemist-controlled ranking

Preferences rerank only recipes that have already passed structure-based
admission, the retrieval support threshold, reaction compatibility, and hard
condition-conflict filtering. A weight cannot admit a rejected precedent,
disable a chemical incompatibility, or force a named family.

The nine visible factors are:

- reaction similarity;
- reacting-partner category similarity, including graph-derived distinctions
  such as secondary aliphatic amine versus primary aryl amine;
- functional-group tolerance, based only on directly observed unchanged
  spectator-group matches and weighted by distance from the reaction center;
- similarity-weighted observed yield;
- independent evidence support;
- distinct reaction breadth;
- dataset diversity;
- post-filter recipe compatibility; and
- condition identity/stage certainty.

Missing functional-group evidence means unknown tolerance, not demonstrated
failure or success. Unavailable factors receive zero applied weight and the
remaining priorities are renormalized. Each recommendation retains its default
rank and score alongside the custom rank, effective weights, contributions,
category labels, matched/missing spectator groups, cautions, and precedent IDs.
The Qt application exposes the same presets and an editable priority dialog.

### Canonical reaction layers

Converted records serialize `reaction_observation` and the optional
`reaction_interpretation` separately. The nested `reaction_signature` is built
from the observation and remains generic; interpretation roles and named-family
evidence stay in the interpretation. Review CSVs expose one `reaction_label`
from the shared terminal renderer, plus status, basis, confidence, and warning
metadata. The reaction core has no display field. Label text does not
participate in identity or retrieval.

### Reaction-core retrieval

Mapped edit observations carry a template-free reaction-core projection beside
the ordinary reaction signature. Converted records and persisted indices retain
its exact, typed, robust-shape, and diagnostic center keys.

For a signed query, reaction-core retrieval is a conservative fallback after
exact bond-edit retrieval and before anonymous edit-graph neighbors. For an
unsigned query with a usable mapped core, a separate query-only branch runs
before the ordinary structure fallback. Both use a narrow-to-broad ladder:
exact `RCX2`, typed-local `RCT2`, then mapping-robust contextual `RSH2`.
Every tier uses verified indexed precedents, matching event count, hard
condition compatibility, and independent support.

Independent support prefers publication identity. When publication identity
is absent, `RME1` mapping equivalence prevents alternate atom-map assignments
from being counted as independent evidence; every source observation remains
available for audit. Passing and review-qualified cores may retrieve, while a
blocked core abstains.

The center-only `RCS2` key is never used for retrieval. External-mapper and
core-only queries remain review-qualified and never become verified indexed
precedents. Successful unsigned-core results use
`recommendation_mode="reaction_core_review"` and require expert review.

### Unresolved-query structure fallback

When parsing succeeds but no verified `ReactionSignature` can be constructed,
the generic API may use a separate `ReactionFallbackDescriptor`. Conversion
creates this descriptor from canonical reactant and product inventories,
taxonomy reactive sites, functional groups, local contexts, candidate interpretation annotations
and edit hypotheses, and net bond/element inventory changes. It is deliberately
not a `ReactionSignature`: candidate edits are hypotheses and are never
serialized as observed bond edits.

The fallback normally searches structure-verified, condition-resolved
precedents in the generic index. A separate review-qualified path handles
source-incomplete attachment replacements such as
`R-C(=O)-OH -> R-C(=O)-F`, `Ar-Br -> Ar-C#N`, or analogous bounded fragment
changes. Taxonomy emits:

- a deterministic `PTS1` partial-transformation key describing the center,
  removed branch, installed fragment graph, and attachment bond;
- one or more typed `FSR1` fragment-source requirements; and
- an explicit product-origin gap. This remains an observation, not a verified
  atom correspondence or `ReactionSignature`.

During conversion,
[`fragment_source_capabilities.v1.json`](definitions/fragment_source_capabilities.v1.json)
matches each requirement against curated resolved condition substances,
families, or narrowly allowlisted raw identities. A precedent is eligible for
the partial index only when every required fragment has a reported
reagent/catalyst source capability. Elemental presence alone is insufficient:
for example, a fluorinated solvent is not inferred to be a fluoride donor.
The capability match supports the recipe as a precedent but does not prove
which atoms in the product came from that component.

At query time, the recommender first searches the exact `PTS1` pool and applies
a hard requirement-ID/source-support filter. Thus an acid-to-acyl-fluoride
query can retrieve acid-to-acyl-fluoride precedents with curated fluorinating
conditions, while an amidation cannot enter that pool. Only after this isolated
path is unavailable does the ordinary unsigned-query policy decide whether a
different fallback is permitted.

Candidate generation is event-first when partial or observed edit tokens are
available: precedents must share the normalized edit before structural ranking.
Similarity then compares deterministic reactant-side reaction-center features
at radius zero through radius three. Radius zero records center element,
hybridization, charge, hydrogen count, ring state, and substitution class;
successive shells compare atom, bond, branching, and topology features without
requiring identical whole-molecule SMILES. Whole-component identity remains a
low-weight tie-breaker.

The ordinary structure-neighbor route requires multiple shared high-signal
feature groups, a configured similarity threshold, conservative compatibility
screening against all observed query functional-group tags, and normally two
independent evidence units. A single precedent is accepted only at the stricter
limited-support threshold. Exact source-supported partial transformations use
the same independent-evidence accounting, but do not admit structurally similar
rows with a different partial-transformation identity.
Results use
`recommendation_mode="unverified_structure_fallback"` and include
`QUERY_TRANSFORMATION_NOT_VERIFIED`,
`UNVERIFIED_REACTION_FALLBACK_USED`, and
`FALLBACK_RECOMMENDATIONS_REQUIRE_EXPERT_REVIEW`. Partial-correspondence
retrieval additionally reports
`QUERY_PRODUCT_ATOM_SOURCE_UNVERIFIED:<element>` and
`SOURCE_SUPPORTED_PARTIAL_TRANSFORMATION_USED`.

The source-supported route still abstains for invalid or conflicting atom maps,
ambiguous scaffold correspondence, transformations outside the bounded
attachment-replacement definition, insufficient product coverage, unsupported
fragment requirements, suspected reactant multiplicity problems, and reactions
without a discriminating structural change. Its current curated capability
seed covers selected single-atom F/Cl/Br/I and multi-atom cyanide/azide sources;
coverage expands by validating the general capability registry, not by adding a
reaction-name branch. An index created before the descriptor, capability,
converter, or index definition changed must be rebuilt.

The desktop recommender and Python API also expose an opt-in
`unrestricted_fallback` expert mode. It bypasses fallback eligibility,
similarity, independent-support, and condition-compatibility gates while
retaining feature-index candidate generation and the finite candidate limit.
It also loads the separate review-core precedent index. The mode is disabled
by default and reports distinct unrestricted-fallback and review-core-index
warnings and cautions; its output is exploratory and may be chemically
incompatible with the query.

Core use is assigned once during canonical conversion by the versioned
`core_eligibility.v1` policy:

- `trusted_core`: may support the default trusted index;
- `review_core`: internally corroborated, full-product-coverage external
  mapping may support only the expert review-core index;
- `query_only`: usable to describe or retrieve for a query, never as a
  condition precedent;
- `blocked` or `unavailable`: excluded from core lookup maps.

The canonical record separately persists `precedent_tier`. A trusted
precedent may remain indexed through a verified signature even when an unsafe
core view is removed. Review-only records enter the expert index only when
both `precedent_tier=review_core` and `core_eligibility=review_core` are
present. Index building does not reinterpret old records. Schema-9.1 records
and schema-5.0 indexes must be regenerated.

## Removed expert rule-based recommendation (historical)

> Historical design record only. The expert-rule modules, definitions, CLIs,
> and public APIs described in this section have been removed. None of the
> following examples is executable or part of the current recommendation
> contract. The material is retained temporarily to support migration review
> and should move to an archived design document once useful rationale has been
> reconciled.

### Design

Expert rules live under `definitions/rule_sets/`. They:

- match allowlisted structural facts from `reactive_taxonomy`;
- never route using source reaction names, display labels, or legacy rule IDs;
- reference explicit recipe templates owned by `condition_registry`;
- apply chemistry scope and compatibility checks before returning a recipe;
- report complete match evidence and failed predicates;
- distinguish reviewed production rules from review-only drafts.

Legacy rulebook identifiers may appear only in provenance. They do not determine
whether a rule matches.

### Phase I C-N scope

The first rule set proposes palladium conditions for verified, single-event,
intermolecular `sp2_c_n_substitution` reactions.

The current hierarchy is:

1. combined hindered-Ar-Cl override;
2. Ar-Cl, hindered-Ar-Br, or alpha-branched-amine override;
3. nucleophile-class default;
4. bounded free-amine Ar-Br fallback;
5. typed abstention when no supported rule matches.

Defaults currently distinguish:

- primary alkyl amines;
- secondary alkyl amines, including cyclic secondary amines;
- primary aryl and heteroaryl amines;
- secondary aryl-containing amines;
- aromatic N-H partners;
- primary carboxamides.

Structural overrides and defaults are selected deterministically. Within the
first matching tier, only the highest-priority rule is selected. A matching
draft override blocks a lower-priority production default, preventing unsafe
fallback around a recognized special case.

The Phase I rules abstain from unsupported handles and chemistry, including
Ar-F, uncurated pseudohalides, intramolecular reactions, multi-event reactions,
and unsupported nitrogen classes. An electron-poor heteroaryl halide may also
support SNAr; proposing palladium conditions does not assert that the product
proves a Buchwald-Hartwig mechanism.

### Production status

The primary-amide/aryl-chloride rule is currently the only active rule. It uses
a located primary-literature protocol and returns a fully specified
tBuBrettPhos Pd G3/K3PO4/t-BuOH recipe.

```python
from condition_recommender import recommend_rule_conditions

result = recommend_rule_conditions(
    "Clc1ccccc1.CC(N)=O>>CC(=O)Nc1ccccc1"
)
```

Production mode excludes every draft rule and draft recipe variant.

### Review the draft C-N rules

The draft C-N defaults and structural overrides are standalone clean-system
rules driven by structural taxonomy facts. Their screening recipes specify
identities, loadings, equivalents, partner stoichiometry, temperature, time,
concentration, and atmosphere, but remain drafts until checked against located
primary procedures.

For chemistry review, export the concise rulebook CSV:

```powershell
python -m condition_recommender.rule_review_cli export `
  results/rule_review/pd_sp2_cn_rulebook.csv
```

The export has one row per rule and explicit recipe variant, with 20 columns
chosen for chemistry review. It includes a readable structural-match summary,
condition names and quantities, operating parameters, rationale, cautions,
source locator, and editable review fields. Alternative condition sets are
separate rows; slash-separated combinations are not generated.

Internal schema fields, individual structural predicates, registry IDs,
definition versions, and provenance JSON are deliberately not flattened into
this sheet. Those belong to the validated canonical rule and recipe-template
definitions, not the chemist review view.

The final columns are blank review fields:

- `review_decision`;
- `reviewer`;
- `review_date`;
- `review_notes`;
- `proposed_change`.

The CSV is a generated review view, not a runtime definition. Editing it does
not change recommendations. Record chemistry feedback in those columns, then
curate and validate the canonical rule JSON and registry recipe-template JSON.
Regenerate the CSV after canonical changes. The `rule_id` and
`recipe_variant_id` link each review row back to those canonical definitions.

The exporter includes both active and draft definitions because its purpose is
review. It does not need `--include-draft`. That flag belongs to reaction
recommendation, where it explicitly opts into returning unreviewed recipes.

Enable them only for review:

```powershell
# Primary-alkyl-amine Ar-Br default
python -m condition_recommender.rule_recommend_cli `
  "Brc1ccccc1.CN>>CNc1ccccc1" `
  --include-draft

# Ar-Cl structural override
python -m condition_recommender.rule_recommend_cli `
  "Clc1ccccc1.CN>>CNc1ccccc1" `
  --include-draft

# Ortho-hindered Ar-Br structural override
python -m condition_recommender.rule_recommend_cli `
  "Cc1cccc(C)c1Br.CN>>Cc1cccc(C)c1NC" `
  --include-draft

# Alpha-branched-amine structural override
python -m condition_recommender.rule_recommend_cli `
  "Brc1ccccc1.CC(C)(C)N>>CC(C)(C)Nc1ccccc1" `
  --include-draft
```

The result reports:

- the selected rule and whether it is a default, structural override, or
  fallback;
- canonical recipe components with quantities;
- temperature, time, concentration, atmosphere, and partner stoichiometry;
- compatibility score, rationale, cautions, and review flags.

This concise review view is the default. Add `--json` to emit the complete
audit payload, including structural match evidence, failures for nonmatching
rules, excluded variants, provenance, and definition versions:

```powershell
python -m condition_recommender.rule_recommend_cli `
  "Clc1ccccc1.CN>>CNc1ccccc1" `
  --include-draft `
  --json
```

Without `--include-draft`, structural draft matches remain visible in the
complete JSON `match_traces`, but their recipes are not returned.

Load and validate the rule definitions directly:

```python
from condition_recommender.rules import load_condition_rule_set

rule_set = load_condition_rule_set()
```

Loading performs strict schema, taxonomy-vocabulary, environment, role,
substance-identity, template-reference, status, and provenance validation.

### Rule path limitations and remaining work

The following work is still required before the C-N drafts can be activated:

- verify every draft against a located primary procedure;
- complete an independent C-N applicability and recipe benchmark;
- use the generated rulebook in blind chemist review and track disagreements;
- calibrate which coordination and functional-group risks should be hard
  exclusions rather than cautions;
- add separately reviewed regimes for currently unsupported handles,
  intramolecular C-N formation, and additional deactivated nitrogen classes.

## Generic structure-backed recommendation

This is the intended type-agnostic precedent path. It operates on canonical
converted records rather than reaction-name partitions.

### Build reference-safe pilot samples

Develop and calibrate against deterministic samples before converting the full
corpus:

```powershell
python -m condition_recommender.sample_cli `
  data-processor/reaction_dataset `
  results/generic_sampling/v1
```

The sampler writes smoke, development, validation, and untouched-test CSVs plus
a versioned manifest. Rows connected by the same normalized publication
reference or exact normalized reaction text remain in one primary partition.
The smoke sample is a subset of development. Original source file and row
provenance survive conversion of the sampled CSVs.

Convert the smoke sample while developing:

```powershell
python -m condition_recommender.generic_conversion_cli `
  results/generic_sampling/v1/smoke.csv `
  results/generic_conversion/smoke_v1
```

Once the smoke path works, use the same restartable converter intended for the
full corpus to prepare development and validation artifacts:

```powershell
python -m condition_recommender.sharded_conversion_cli `
  results/generic_sampling/v1/development.csv `
  results/generic_conversion/v2/development `
  --mode development --shard-size 100 --workers 4

python -m condition_recommender.sharded_conversion_cli `
  results/generic_sampling/v1/validation.csv `
  results/generic_conversion/v2/validation `
  --mode validation --shard-size 100 --workers 4

python -m condition_recommender.generic_index_cli `
  results/generic_conversion/v2/development/records.jsonl.gz `
  results/generic_conversion/v2/development/generic_index.sqlite

python -m condition_recommender.generic_index_cli `
  results/generic_conversion/v2/validation/records.jsonl.gz `
  results/generic_conversion/v2/validation/generic_index.sqlite
```

Rerunning a sharded conversion with unchanged sources and definitions reuses
checksum-valid completed shards. A source, schema, taxonomy, registry, or
converter-definition change invalidates reuse instead of silently mixing
versions.

### Audit and convert a structure-rich dataset

Audit source data without modifying it:

```powershell
python -m condition_recommender.audit_cli `
  data-processor/reaction_dataset `
  results/reaction_dataset_audit `
  --chemistry-sample-per-file 100
```

Convert a CSV file or directory:

```powershell
python -m condition_recommender.generic_conversion_cli `
  data-processor/reaction_dataset `
  results/generic_conversion
```

For unresolved or ambiguous reactions, the standard converter can optionally
run the installed RXNMapper provider:

```powershell
python -m pip install -r requirements-mapping.txt

python -m condition_recommender.generic_conversion_cli `
  data-processor/reaction_dataset/Fischer_indole_synthesis.csv `
  results/fischer_indole_rxnmapper_conversion `
  --use-rxnmapper
```

The generated mapping, provider/model provenance, confidence, coverage,
internal-hypothesis match, and disposition are stored in canonical JSONL and
review CSV. Mapper-derived records remain `review_only` and are excluded from
the default index.

The converter writes canonical nested `records.jsonl`, tiered CSV review views,
and JSON/Markdown coverage reports. Exact reconstructed signatures and valid
mapped signatures may be admitted even when `named_family` is absent. Source
family labels are provenance only. Chemistry, condition, outcome, and index
eligibility are assessed independently. The legacy verified/review/rejected
tier remains a compatibility summary, not the indexing contract.

For rapid reaction-family review, export the compact chemistry-facing columns:

```powershell
python -m condition_recommender.concise_review_cli `
  results/generic_conversion/smoke/records.jsonl.gz `
  results/generic_conversion/smoke/concise_reaction_review.csv
```

The output contains canonical reaction SMILES, the detailed structural label,
the original source reaction type, the detected optional reaction family,
transformation class, signature and evidence identifiers, completeness
diagnostics, admission/index statuses, warning codes, unchanged spectator
groups, and a compact steric/electronic summary for each reactive partner. In
the compact fields,
`d` is graph-bond distance from the nearest reactive site, `S` means steric
context, `E` means electronic context, and `q` is the qualitative local
electronic index—not a measured physical quantity. Blank detected families are
retained and are not treated as errors.

The desktop app can now build both recommendation-ready data and the simple
review CSV:

```powershell
python -m app.reaction_converter_gui
```

Choose any mixture of source files and folders. The output defaults to
`datasets/literature`, but it remains editable. Select **Full** to convert every
record or **Compact** for a faster development library. Compact retains every
record in files with at most 200 records; for larger files it retains the first
200 plus a content-seeded random 15% sample of the remainder, rounded up. The
seed makes reruns and resumes deterministic. Full and Compact artifacts are
isolated under `full/` and `compact/`, respectively. If the selected output
already contains a root-level pre-mode library, that location remains the Full
library so its expensive completed batches do not need to be rebuilt; Compact
still uses `compact/`.

**Convert & Save Batch** stores the selection as an independent, restartable
folder under the selected mode's `batches/`. You can repeat this for later
groups of files. **Combine Saved Batches / Build Index** then streams every
saved batch in the selected mode, removes identical duplicate observations, and
rebuilds that mode's active recommender artifacts. If two saved
records have the same observation ID but different content, combining stops
instead of silently choosing one. Combining also refuses any cancelled or
otherwise incomplete batch, so partial source coverage cannot replace a
previously complete recommender index. In the desktop converter, the combine
button offers to resume those stored source selections first, reusing valid
checkpointed shards, and proceeds to indexing only after coverage is complete.

Each mode-specific batch library creates:

- `batches/<batch-name>/shard_manifest.json` and compressed shards: canonical
  data for each independently reusable batch. A blank batch name produces a
  stable name from the selected inputs, so cancellation or restart reuses
  completed shards.
- `combined_records.jsonl.gz`: the deduplicated canonical corpus represented by
  the active index.
- `reaction_review.csv`: the combined compact file for quick human review, including
  structural-evidence and admission diagnostics, spectators, and
  reactive-partner steric/electronic context. When RXNMapper is enabled it also
  shows mapping disposition, provider, and confidence. Machine-only hash
  identifiers and internal grouping keys are omitted; they remain available in
  the canonical shards. The CSV is not used as recommendation input.
- `generic_index.sqlite`: the trusted-precedent runtime index used by default.
  Rows and lookup values are materialized lazily for fast application startup.
- `generic_review_index.sqlite`: the paired expert-use runtime index containing
  the trusted rows plus qualified review-core precedents. The combined workflow
  writes this explicitly even when it currently contains the same reactions as
  the trusted index, keeping its expert-use scope unambiguous across batches.
- `combined_batch_manifest.json` and `combined_recommendation_report.json`:
  included batches, checksums, duplicate counts, output paths, and index counts.

The checkbox beside the conversion settings refreshes the combined index
automatically after a batch is saved. It is unchecked by default; enable it for
a one-step save-and-refresh, or use the combine button after saving all batches.

The default settings are intended for large datasets:

- `1,000` rows per shard gives useful restart points without creating too many
  small files.
- `Use RXNMapper` is unchecked by default, leaving parallel workers available.
  Enabling it uses one conversion worker so only one model copy is loaded.
- Mapper-only or conflicting records never become precedents. Full-coverage
  external mappings that agree with an internal hypothesis may enter only the
  expert review-core index. The mapper setting and model hash participate in
  shard reuse identity.
- Automatic combining is unchecked by default. Save all intended batches, then
  combine once to avoid repeated index builds.

Cancellation is safe: the active shard finishes, the manifest is checkpointed,
and choosing the same source and output folders later resumes by reusing valid
completed shards.

### Build an index and recommend

```powershell
python -m condition_recommender.generic_index_cli `
  results/generic_conversion/records.jsonl `
  results/generic_conversion/generic_index.sqlite

python -m condition_recommender.generic_index_cli `
  results/generic_conversion/records.jsonl `
  results/generic_conversion/generic_review_index.sqlite `
  --include-review-core

python -m condition_recommender.generic_recommend_cli `
  "<reaction_smiles>" `
  --records results/generic_conversion/generic_index.sqlite
```

Add `--use-rxnmapper` for an unmapped query that ordinary featurization cannot
resolve. A mapper-supported query can retrieve independently verified
precedents, but the result reports the external-mapping status and adds
expert-review cautions:

```powershell
python -m condition_recommender.generic_recommend_cli `
  "O=C1CCCCC1.Cl.NNc1ccc(F)cc1>>Fc1ccc2[nH]c3c(c2c1)CCCC3" `
  --records results/generic_conversion/generic_index.sqlite `
  --use-rxnmapper
```

For output made by the desktop app, use the fast index directly:

```powershell
python -m condition_recommender.generic_recommend_cli `
  "<reaction_smiles>" `
  --records datasets/literature/full/generic_index.sqlite

# Expert mode resolves the paired generic_review_index.sqlite automatically.
python -m condition_recommender.generic_recommend_cli `
  "<reaction_smiles>" `
  --records datasets/literature/full/generic_index.sqlite `
  --unrestricted
```

If the fast index was not built, use the canonical manifest instead:

```powershell
python -m condition_recommender.generic_recommend_cli `
  "<reaction_smiles>" `
  --records datasets/literature/full/shard_manifest.json
```

### Desktop recommender

Launch the simpler Qt6 interface:

```powershell
python -m app.reaction_recommender_gui
```

Use the **Full / Compact** selector to switch between
`datasets/literature/full/` and `datasets/literature/compact/`. Existing
root-level artifacts remain a temporary Full fallback until a new Full library
is built. The app prefers `generic_index.sqlite`, then `shard_manifest.json` as
the rebuild fallback. Persisted JSON runtime indexes are no longer accepted. Paste a complete
`reactants>>product` reaction SMILES and choose how many recipes to return.
`Use RXNMapper` is checked by default; it is invoked only when supplied mapping
and ordinary internal analysis do not already resolve the query. Clear the
checkbox to use internal evidence alone.

The window shows:

- the detected optional reaction family and generic transformation;
- the retrieval fallback level and compatible precedent counts;
- ranked catalysts, ligands, bases, solvents, additives, temperature, and time;
- similarity, compatibility, expected yield, and independent support;
- explanations, cautions, and precedent IDs for the selected recipe;
- external-mapping disposition, provider, confidence, and mandatory review
  cautions when RXNMapper participates;
- the displayed precedent hit's reaction SMILES, structural label, spectator
  groups, and local steric/electronic partner analysis.

Recommendation runs in a background thread. The validated index is cached in
memory after its first load and is automatically reloaded if the selected file
changes. Results can be exported as the complete typed JSON response.

Retrieval follows a chemistry hierarchy: exact signature, relaxed handle
signature, high-confidence family, generic transformation, then compatible bond
edits narrowed by local-environment neighbors, followed by broad compatible
bond edits, exact reaction-core shape, and anonymous edit-graph neighbors. Hard
chemistry and recipe-compatibility gates run before similarity and aggregation.
Minimum support is counted by independent
publication, or by canonical reaction when no reference is available; repeated
scope examples from one paper cannot satisfy the threshold by themselves.

If the query has no signature, this hierarchy is not entered. Typed edit
hypotheses and eligible reaction cores have their own query-only review routes;
otherwise the unverified structure fallback uses its own similarity gate and
trace. None may be presented as a verified bond-edit-compatible result.

Results aggregate by canonical `RCORE1` recipe core and expose the observed
`RCR1` operating-condition variants. Support distinguishes raw observations,
canonical reactions, reference-local condition series, independent references,
and datasets. Repeated rows from one publication contribute one evidence unit
to similarity and yield aggregation. Expected yield is omitted when no usable
outcomes support the recipe; it is never invented from missing values. In that
case the unavailable yield term receives zero applied weight and the available
ranking terms are renormalized.

Generic configuration is split by responsibility:

- `definitions/generic_retrieval.v1.json` selects the narrowest adequately
  supported candidate pool;
- `definitions/generic_similarity.v1.json` weights interpretable molecular
  comparison features;
- `definitions/generic_ranking.v1.json` combines similarity, compatibility,
  condition certainty, outcome evidence, and reference-aware support.
- `definitions/fragment_source_capabilities.v1.json` states which curated
  reported condition identities can satisfy typed product-fragment
  requirements;
- `definitions/fallback_retrieval.v1.json` controls the isolated unresolved-query
  structure fallback and its stricter support gates;
- `definitions/reaction_core_retrieval.v3.json` controls unsigned reaction-core
  evidence status, blocking warnings, and independent support.

Ranking and retrieval parameters are selected on development splits and
promoted only after a separate validation gate. Similarity weights remain an
explicit chemistry prior rather than a statistical fit.

Every result includes a typed retrieval trace. Each attempted tier reports raw
rows, independent support units, compatibility-filtered rows, exclusions, the
configured minimum support, and why the tier was selected or skipped.
Each recommendation also includes a typed score trace with similarity feature
values and contributions, ranking components and contributions, applied
weights, definition versions, and evidence/outcome counts.

For an application or batch process, load the validated index once:

```python
from condition_recommender import GenericConditionRecommender

recommender = GenericConditionRecommender.from_path(
    "results/generic_conversion/generic_index.sqlite"
)
first = recommender.recommend("<reaction_smiles>")
second = recommender.recommend("<another_reaction_smiles>")
```

`recommend_generic_conditions()` remains the convenient one-shot entry point.

Reaction topology participates in the more specific signature tiers. A
topology-agnostic fallback is allowed only with an explicit
`REACTION_TOPOLOGY_FALLBACK_USED` warning. Multi-event queries and precedents
are compared using their complete event and net-edit sets.

### Evaluate recommendations

Think of evaluation as a closed-book test:

1. The system receives 80% of the reactions as examples.
2. The other 20% are hidden.
3. The system recommends conditions for each hidden reaction.
4. We compare its suggestions with the conditions actually reported.

Reactions from the same paper are always kept on the same side. Otherwise the
system could see a nearly identical example from the same paper and the score
would look better than it really is.

Run the basic test:

```powershell
python -m condition_recommender.evaluation_cli `
  results/generic_conversion/v2/development/generic_index.sqlite `
  results/generic_evaluation/v2/development/basic `
  --test-fraction 0.2 --seed 17 --top-k 5
```

Read `evaluation_report.md` for the short report or
`evaluation_report.json` for all details. Start with these numbers:

| Result | Plain meaning | What we want |
| --- | --- | --- |
| Coverage | How often the system returned any suggestion | Higher is useful, but not at the cost of unsafe chemistry |
| Seen-recipe top-1 | How often the first suggestion matched a known recipe that was available in the examples | Higher is better |
| Seen-recipe top-5 | How often the known recipe appeared anywhere in the first five suggestions | Higher is better |
| Hard-incompatible count | Suggestions that violate a chemistry compatibility rule | Must be zero |
| Yield error | Average difference between predicted and reported yield, in percentage points | Lower is better; use only when enough yields are available |
| Abstention | Cases where the system correctly said it had insufficient support | Review these; abstaining can be safer than guessing |

Always check `seen_recipe_query_count` beside the top-1 and top-5 values. A high
percentage based on only a few reactions is weak evidence. Also compare a new
version with the current version on exactly the same split and seed.

Recipe recovery is deliberately strict: a chemically reasonable alternative
still counts as a miss when it differs from the recorded recipe. This is why
the numerical test and independent chemist check are both needed.

Calibrate the reaction-core path separately:

```powershell
python -m condition_recommender.core_evaluation_cli `
  results/generic_conversion/v2/development/generic_index.sqlite `
  results/reaction_core_evaluation/v1/development `
  --test-fraction 0.2 --seed 17 --top-k 5
```

This writes exact/local/context and pass/review metrics,
mapping-equivalence deduplication statistics, per-query JSONL, and
`reaction_core_chemist_review.csv`. The review CSV deliberately leaves the
chemist assessment and notes blank and omits the expected recipe.

A simple acceptance checklist is:

- no paper or reaction appears in both the example and hidden sets;
- no hard-incompatible recipe is returned;
- top-1 or top-5 recovery improves without a large loss of coverage;
- errors and abstentions have understandable explanations;
- an independent chemist finds no repeated chemistry problem.

There is no single “good accuracy” threshold for every reaction class. Report
the number of tested reactions and results by transformation class instead of
only one overall percentage.

#### Optional stronger tests

After the basic test works, test reactions with different scaffolds, different
source datasets, and newer publication years:

```powershell
python -m condition_recommender.evaluation_cli `
  results/generic_conversion/v2/validation/generic_index.sqlite `
  results/generic_evaluation/v2/validation/scaffold `
  --split-mode scaffold_disjoint --seed 71

python -m condition_recommender.evaluation_cli `
  results/generic_conversion/v2/validation/generic_index.sqlite `
  results/generic_evaluation/v2/validation/source `
  --split-mode source_disjoint --seed 71

python -m condition_recommender.evaluation_cli `
  results/generic_conversion/v2/validation/generic_index.sqlite `
  results/generic_evaluation/v2/validation/time `
  --split-mode forward_time
```

Compare the available retrieval approaches on the same hidden reactions:

```powershell
python -m condition_recommender.baseline_cli `
  results/generic_conversion/v2/validation/generic_index.sqlite `
  results/generic_evaluation/v2/validation/baselines `
  --seed 71
```

Calibration uses development data to choose settings and validation data only
to decide whether those settings are safe to promote:

```powershell
python -m condition_recommender.calibration_cli `
  results/generic_conversion/v2/development/generic_index.sqlite `
  results/generic_conversion/v2/validation/generic_index.sqlite `
  results/generic_evaluation/v2/calibration
```

#### Chemist check

Generate a blind review packet before using the untouched test:

```powershell
python -m condition_recommender.chemist_review_cli <index> <output>
```

The chemist opens `review_packet.html`, judges each condition set, and records
the answer in `review_form.csv`. They must not open `answer_key.jsonl` first.
Allowed answers are `compatible`, `compatible_with_caution`, `incompatible`,
and `uncertain`.

After every row is complete, create the signed review summary:

```powershell
python -m condition_recommender.chemist_review_adjudication_cli `
  <review-output> <review-output>/review_summary.json `
  --reviewer "<name>" --independent-reviewer --sign-off
```

Do not sign off while a repeated chemistry problem remains unresolved.

### Restartable full conversion

After chemist review and untouched-test gates pass:

```powershell
python -m condition_recommender.sharded_conversion_cli `
  data-processor/reaction_dataset results/generic_conversion/v2/full `
  --shard-size 100 --workers 6

python -m condition_recommender.conversion_integrity_cli `
  results/generic_conversion/v2/full/shard_manifest.json

python -m condition_recommender.generic_index_cli `
  results/generic_conversion/v2/full/records.jsonl.gz `
  results/generic_conversion/v2/full/generic_index.sqlite

python -m condition_recommender.generic_index_integrity_cli `
  results/generic_conversion/v2/full/generic_index.sqlite
```

Canonical shards and catalogs are deterministic compressed nested JSONL. The
manifest binds each shard to source and output checksums plus every
participating schema and chemistry/registry definition.

Compose the final machine and independent-human gates with
`condition_recommender.release_validation_cli`. A review summary is accepted
only when it has the adjudication schema, independent reviewer sign-off, no
unresolved systematic defect, and hashes for the packet, form, report, and
answer key.

### Test and validate the new system

Run the complete deterministic test suite before handing off any change:

```powershell
pytest -q
```

Run package-focused suites while developing:

```powershell
# Molecular observations, edits, environments, and reaction signatures
pytest -q tests/reactive_taxonomy

# Condition identities, roles, compatibility vocabulary, and recipes
pytest -q tests/condition_registry

# Conversion, indexing, retrieval, ranking, evaluation, and release gates
pytest -q tests/condition_recommender
```

Useful focused regressions for the generic path are:

```powershell
pytest -q tests/condition_recommender/test_generic_conversion.py
pytest -q tests/condition_recommender/test_generic_retrieval.py
pytest -q tests/condition_recommender/test_generic_evaluation.py
pytest -q tests/reactive_taxonomy/test_reaction_signatures.py
```

Dataset testing has four distinct levels:

1. **Mechanical smoke:** convert one 100-row shard, build an index, run index
   integrity, and request at least one recommendation using the quick-start
   commands.
2. **Development:** use the reference-safe development sample for
   implementation, baseline comparison, and candidate selection.
3. **Validation:** run grouped, scaffold-disjoint, source-disjoint, and
   forward-time reports without tuning directly against them.
4. **Untouched and full corpus:** run only after independent chemist review and
   adjudication pass.

Validate sharded artifacts and their persisted indices independently:

```powershell
python -m condition_recommender.conversion_integrity_cli `
  results/generic_conversion/v2/development/shard_manifest.json

python -m condition_recommender.generic_index_integrity_cli `
  results/generic_conversion/v2/development/generic_index.sqlite `
  --output-path results/generic_conversion/v2/development/index_integrity.json
```

A machine-valid sampled run should satisfy all of the following:

- shard input and output row counts agree and no shard failed;
- source and output checksums match the manifest;
- index schema and all chemistry/registry definition versions are current;
- duplicate observation count is zero;
- reference and canonical-reaction overlap counts are zero;
- strict scaffold overlap is zero in the scaffold-disjoint report;
- hard-incompatible recommendation count is zero;
- recovery metrics state their seen-recipe denominator and yield metrics state
  their usable-outcome count;
- invalid or unsupported query chemistry returns a typed error or abstention.

Do not treat successful unit tests, high coverage, or a checksum-valid index as
evidence that a recipe is chemically suitable. Chemistry review and held-out
metrics are separate required gates.

### Generic path limitations

- Performance depends on the quality and coverage of the supplied converted
  records.
- Current evaluation datasets are pilots and do not establish broad production
  validity.
- Registry coverage and role normalization are still incomplete for some
  source datasets.
- Calibration is sample-level and is not broad production-accuracy evidence.
- Production release remains gated by independent chemist review and the
  untouched-test report.

## Removed weak-label recommendation (historical)

> Historical design record only. The weak-label recommender and CLI described
> in this section have been removed. Label-only observations may still be
> ingested and condition-normalized for audit, but they are structurally
> ineligible as recommendation precedents and none of the following commands is
> part of the current public workflow.

The weak-label path uses a cleaned flat CSV whose precedent rows contain source
reaction labels and reactive-site labels rather than reaction structures.

```powershell
python -m condition_recommender.recommend_cli `
  "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
```

Options:

```powershell
python -m condition_recommender.recommend_cli "<reaction_smiles>" --top-k 10

python -m condition_recommender.recommend_cli "<reaction_smiles>" `
  --records datasets/reaction_label/v2.1_cleaned.csv
```

The query reaction is structure-verified with `reactive_taxonomy`, but precedent
rows are not. The implementation:

1. assigns a query reaction interpretation annotation;
2. restricts source reaction types through a declarative crosswalk;
3. matches two reactive participants as an unordered pair;
4. ranks by functional-group signature and structural qualifiers;
5. uses yield, z-score, and support as lower-weight signals;
6. aggregates canonical `condition_registry` recipes.

Weights live in `definitions/label_retrieval.v1.json`.

Weak handle signatures may be family-level prefixes when the source label lacks
structural detail. Known incompatible handle families or `XH` center classes
are hard mismatches; an unknown/blank participant may still contribute a
partial match through the other participant. This allows `XH|Csp3` precedents
to support an exact activated-carbon query without allowing N–H or O–H rows to
leak into the candidate recipes.

For the cleaned HTE data, `CH-Activation` is routed to the intermolecular
`Ar–X + Ar–H → Ar–Ar` interpretation annotation. The same source reaction type also contains
`Ar–H + alkyne` records, so the generic `alkyne` label is normalized to
`PI|Alkyne`; those records then fail participant compatibility instead of
entering direct-arylation recommendations. Friedel–Crafts remains a separate
structural interpretation and is not backed by this CSV because the dataset has no
corresponding acylation source type.

Source reaction context disambiguates the same `alkyne` label in Sonogashira
rows as the terminal reactive handle `XH|Csp|H1|Alkynyl`. Two-partner
retrieval requires both precedent signatures, and transfer handles from
different element families (for example Zn, Mg, Sn, and B) are hard
incompatible even when the generic transformation pattern is shared.

Important limitations:

- precedent reactions are not structure-verified;
- source reactive-site order is unreliable and treated as unordered;
- unresolved condition names remain explicit with identity uncertainty;
- reaction topology cannot be recovered from the precedent rows;
- only structurally verified intermolecular queries are accepted;
- recommendations are weak-label precedents, not verified literature matches.

The path returns explicit warnings including
`WEAK_LABEL_PRECEDENTS_NOT_STRUCTURE_VERIFIED` and
`UNRESOLVED_CONDITION_IDENTITIES_RETAINED_WITH_UNCERTAINTY`.

The cleaned CSV uses a component-oriented condition contract instead of fixed
`catalyst`/`base`/`solvent` slots:

- `condition_recipe_id` is the canonical `RCR1` aggregation key;
- `condition_display` is a compact human-review label;
- `temperature_c`, `time_h`, `concentration_m`, and `atmosphere` are convenient
  single-stage review fields;
- `procedure_text` preserves source text and never contributes opaque text to
  recipe identity.

The paired `<csv-stem>.condition_recipes.jsonl.gz` catalog stores each unique
full recipe exactly once, including resolved identities, contextual roles,
provenance, stages, declared absences, and uncertainty. It is losslessly
compressed and read transparently by the recommender. This keeps the CSV small
without making the recommendation path parse display strings.

Multi-stage procedures remain ordered in the recipe catalog; their
top-level temperature/time cells stay blank rather than collapsing multiple
stages. Source columns are role hints only and cannot override contradictory
registry chemistry.

Regenerate the cleaned CSV with:

```powershell
python -m scripts.clean_reaction_label_csv `
  datasets/reaction_label/v2.1.csv `
  datasets/reaction_label/v2.1_cleaned.csv
```

## Shared compatibility behavior

`definitions/compatibility.v1.json` compares unchanged query spectator groups
with contextual recipe roles and condition families.

Hard conflicts currently include explicit oxidation, reduction, acid, moisture,
and atmosphere conflicts. Hydrolysis, acid/base, metal coordination, and
catalyst-poisoning risks generally remain eligible with auditable penalties and
cautions. Within one penalty category, only the strongest matching penalty is
applied. The definition loader validates query tags against
`reactive_taxonomy`, and recipe buckets and family IDs against
`condition_registry`, so a future taxonomy extension fails explicitly if the
matching compatibility vocabulary is stale.

Compatibility evidence does not repair an unsupported transformation or make a
review-qualified recipe production-ready.

## Input contract

Recommendation APIs require a complete product-specified reaction SMILES:

```text
reactant1.reactant2>>product
```

Invalid, unsupported, ambiguous, or chemically conflicting inputs return a
typed error or abstention rather than an unrelated recommendation.

## Development status

The package exposes one canonical structure-backed recommendation workflow.
Verified-signature retrieval, review-qualified query routes, source-supported
partial transformations, and unsigned structure fallback are evidence modes of
that workflow rather than separate recommenders. Application code must preserve
their `recommendation_mode`, retrieval trace, warnings, and provenance.

Production release remains gated by full-corpus reproducibility, leakage-safe
evaluation, independent chemist review, and removal or archival of the
historical documentation above.
