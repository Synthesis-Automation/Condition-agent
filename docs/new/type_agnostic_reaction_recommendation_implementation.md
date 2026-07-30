# Type-Agnostic Reaction Recommendation: Implementation Status and Roadmap

**Status:** Consolidated implementation reference  
**Reviewed against code, definitions, tests, and local artifacts:** 2026-07-30  
**Companion design:** 
[`reaction_condition_recommendation_design_for_chemists.md`](reaction_condition_recommendation_design_for_chemists.md)

## 1. Purpose and authority

This document replaces the fragmented status narrative spread across the older
files in `docs/`. It states what is implemented now, what is still provisional,
and which gates must pass before production release.

For conflicts, use this order:

1. executable contracts, versioned definitions, validation, and tests;
2. this consolidated status and roadmap;
3. the companion chemist design guide;
4. older detailed design and migration documents.

Older documents remain useful for rationale and design history. Their
`proposed`, `in progress`, phase, schema-version, artifact-path, and
`recommended immediate work` statements are not automatically current.

The package boundary remains mandatory:

```text
reactive_taxonomy       condition_registry
          \                 /
           condition_recommender
                    |
               app / CLI / API
```

- `reactive_taxonomy` owns molecular and reaction chemistry.
- `condition_registry` owns condition-substance identity, roles, and recipes.
- `condition_recommender` owns conversion, admission, retrieval, compatibility,
  ranking, evaluation, and explanation.
- The standalone packages do not depend on legacy `chemtools`.
- Application layers compose packages but do not own chemistry rules.

## 2. Current implementation snapshot

### 2.1 What is implemented

| Area | Current state | Principal implementation |
| --- | --- | --- |
| Molecular features | Implemented | Functional groups, reactive sites, canonical connectivity interfaces, typed reactivity profiles |
| Reaction parsing | Implemented | Two- and three-part reaction SMILES with component, map, and source preservation |
| Connectivity execution | V2 implemented | One bounded declarative rewrite path for all registered grammars; no legacy reaction-operator fallback |
| Edit evidence | Implemented | Mapping, exact single/multi-event reconstruction, conservative scaffold correspondence, bounded global and fragmented-scaffold correspondence, conflicts, H/charge/stereo observations |
| Product completeness | Implemented | Verified/incomplete/unresolved accounting plus observation-only product-origin gaps |
| Reaction signatures | Implemented | Deterministic RS3 L0–L4 signatures, events, topology, profiles, spectators, unknown-family support |
| Reactivity descriptors | Implemented and active | Typed context-aware profiles are the sole active environment path |
| Condition registry | Implemented, curation incomplete | Conservative identity resolution, contextual roles, RCORE1/RCR1 recipes, stages, provenance |
| Generic conversion | Implemented | Nested canonical records, independent quality dimensions, review exports, sharding, restart/integrity checks |
| Generic index | Implemented | Version-checked persisted index with signature, environment, family, fallback, recipe, and reference keys |
| Generic retrieval | Implemented pilot | Explicit ladder, independent-support threshold, hard compatibility, similarity, reference-aware recipe aggregation |
| Unverified-query fallback | Implemented, conservative | Separate structure-derived fallback; not represented as verified edit retrieval |
| Evaluation and calibration | Implemented tooling | Grouped, scaffold-, source-, and time-disjoint modes; baselines; calibration; blind review/adjudication |
| Release validation | Implemented tooling | Machine artifact checks plus hash-bound independent chemist sign-off |
| Expert rules | Separate limited path | One active C–N protocol; other templates are draft/review-only |
| Weak-label retrieval | Transitional separate path | Works with structure-poor precedents but carries explicit weaker-evidence warnings |

The generic structure-backed path is the intended cross-family direction.
There is no automatic fallback among generic retrieval, expert rules, and
weak-label retrieval. A caller must select the path and preserve its provenance.

### 2.2 Active contract versions

The current code declares:

| Contract | Version |
| --- | --- |
| Reaction analysis | `3.1` |
| Reaction signature / ID namespace | `3.0` / `RS3` |
| Taxonomy identity manifest | `3.0` |
| Connectivity site interface | `2.0` |
| Connectivity rewrite | `2.0` |
| Typed reactivity profile | `1.0` |
| Reaction fallback descriptor | `1.2` |
| Resolved condition recipe | `1.2` |
| Recommendation record | `3.4` |
| Generic converter definition | `generic_conversion.v2.4` |
| Generic persisted index | `2.3` |
| Generic recommendation result | `2.1` |
| Reaction correspondence definitions | `2.3` |
| Generic retrieval definition | `1.6` |
| Generic admission policy | `generic_admission.v1.7` |

Do not copy this table into executable code. The constants and definition files
remain authoritative, and stale artifacts must fail validation rather than
silently mixing chemistry identities.

### 2.3 Current local literature artifact

The checked local artifact under `datasets/literature/` is a complete
conversion of the bounded `examples/all_dataset_300` collection, not the full
source corpus. Its report records:

| Measurement | Count |
| --- | ---: |
| Source files / converted rows | 118 / 30,100 |
| Verified product-provenance signatures | 16,118 |
| Index-eligible records | 12,816 |
| Chemistry verified / review / rejected | 7,417 / 13,793 / 8,890 |
| Records without an assigned named family | 28,351 |
| Generic-index rows | 12,816 |
| Unique references in the index | 5,330 |
| Unique recipe cores in the index | 3,820 |

The sharded conversion report records zero failed shards and zero duplicate
observations. The current generic-index integrity command also validates the
artifact as schema `2.3`.

These counts demonstrate that the type-agnostic path works on a broad bounded
sample and that unnamed reactions are retained. They are coverage and integrity
measurements, not recommendation-accuracy claims. Because the source is capped
per dataset, this artifact must not be described as a full-corpus production
index.

### 2.4 Grammar-independent Fischer POC

The executable benchmark
`benchmarks/fischer_indole_edit_poc.py` evaluates the 542-row Fischer source
cohort without passing its source reaction name into featurization, signature
identity, or retrieval. The current report is
`results/fischer_indole_edit_poc.json`.

The 2026-07-30 run found 356 unique reaction SMILES:

| Measurement | Count |
| --- | ---: |
| Existing global-correspondence signatures | 54 |
| New fragmented-scaffold signatures | 4 |
| Total signatures | 58 |
| Ambiguous correspondences retained as abstentions | 205 |
| Anonymous edit prototypes among signatures | 3 |
| Reference-disjoint signed queries with exact L3 support | 54 / 58 |
| Reference-disjoint signed queries with edit-graph support | 58 / 58 |
| Reference-disjoint queries linked across precursor modes | 58 / 58 |

This is a successful retrieval POC but not a solved atom-mapping POC. The
approximate edit graph connects the separate-reactant and preformed-hydrazone
views for every currently signed query. The conservative correspondence
provider admits only four additional reactions because competing nitrogen
origins or symmetry-related mappings remain genuinely ambiguous from the
unmapped source structures. Those 205 cases must not be bulk-admitted without
curated mapping or an independently validated correspondence provider.

The normal generic converter was also run over all 542 source rows. It produced
59 signed observations and 45 index-eligible precedents; the compact report is
`results/fischer_indole_poc_conversion/conversion_report.json`. A representative
accepted preformed-hydrazone query used the persisted index and selected
`edit_graph_neighbors`: 41 compatible candidates from 18 independent
references. Its exact-signature through L3 pools contained only one independent
reference, so the new tier supplied genuinely broader support. The result
retained `named_family=None`, disclosed the intramolecular/intermolecular scope
mismatch, and warned that the exact edit signatures differed.

## 3. Implemented chemistry contracts

### 3.1 Observation

`reactive_taxonomy` currently provides:

- immutable atom-provenanced bond and schema-level hydrogen edits;
- stronger internal before/after connectivity observations with explicit
  observed, projected, reconstructed, inferred, or unresolved scope;
- formal-charge and explicit stereochemical observations;
- canonical molecular reactive links, bond capacities, and connection
  endpoints;
- exact single-event and composable multi-event reconstruction;
- reaction topology and event relationships;
- typed product-atom completeness and partial product-origin gaps;
- graph-derived local reactivity profiles and unchanged spectators.

The internal `CEG1` connectivity-edit graph remains evaluation/shadow output.
The public and persisted chemistry identity is the RS3 reaction signature.
Connectivity V2 refers to the active rewrite executor, not promotion of CEG1
into public identity.

### 3.2 Interpretation

Structural archetypes are derived from edits. Transformation classes and named
families are optional interpretations. Mapped or otherwise verified unknown
chemistry can receive:

```text
transformation_class = generic_graph_transformation
named_family = None
```

Source-declared family is stored as provenance outside
`featurize_reaction()`. It cannot determine the structural output.

Display labels are built after the structural analysis. Their wording and
rendering versions do not affect signature identity.

### 3.3 Recommendation records

The converter serializes:

- source row, file, reaction, and reference provenance;
- reaction analysis, signature, completeness, fallback, and conflicts;
- chemistry, condition, stage, outcome, and index-eligibility statuses;
- resolved recipe core and variants;
- reference-local condition-series identity;
- outcome values and their usability;
- every participating schema and definition version.

Compressed nested JSONL shards are the current canonical implementation
artifact. CSV is a chemist review/export view. The longer-term preference for
partitioned Parquet at full-corpus scale remains reasonable but is not the
current implementation.

## 4. Current retrieval and ranking behavior

The verified-signature ladder in
`condition_recommender/definitions/generic_retrieval.v1.json` is:

```text
exact signature
  -> handle signature
  -> high-confidence named family
  -> transformation signature
  -> environment neighbors within compatible bond edits
  -> broad compatible bond edits
  -> chemistry-gated anonymous edit-graph neighbors
  -> abstain
```

The current minimum adequate pool is two independent support units. Family use
requires high confidence and never bypasses edit compatibility. The anonymous
edit-graph tier may cross an exact L3 key only after shared formed/broken bond
and ring-direction gates pass, and its use is disclosed in recommendation
cautions.

Compatibility runs before similarity. Ranking definitions are separated from
retrieval definitions. Recipes are aggregated by `RCORE1`, while `RCR1`
variants, independent references, reference-local series, reaction breadth,
dataset breadth, usable yields, missingness, and uncertainty remain visible.

The fallback descriptor and retrieval policy are separate from RS3 retrieval.
They must return `recommendation_mode`, fallback evidence, stricter thresholds,
and cautions. A structure-neighbor fallback must never be described as a
verified reaction-center match.

## 5. What is not complete

### 5.1 Production validation

The software pipeline is substantially ahead of the early phase documents, but
production validation is not complete:

- current performance measurements are pilot/sample measurements;
- the repository does not contain the previously described current blind
  chemist-review packet and signed adjudication summary;
- the untouched-test gate has not been demonstrated against the current
  contracts and definitions;
- the checked literature artifact is based on a per-dataset bounded sample,
  not the full source corpus;
- no broad production-accuracy claim is justified.

### 5.2 Chemistry and data coverage

Remaining limitations include:

- incomplete or missing participating reactants;
- unresolved product provenance and atom correspondence;
- unstructured multi-stage procedures;
- incomplete condition-substance and contextual-role curation;
- sparse independent support for many transformation/environment combinations;
- complex insertion, extrusion, metathesis, migration, annulation,
  cycloaddition, rearrangement, cascade, and protection chemistry without
  fully validated general contracts;
- coordination, ionic association, isotope exchange, and some radical or
  charge-only chemistry outside the ordinary covalent-edit model.

The template registry can preserve and execute selected mapped reference
transformations, including explicit multiplicity and curated role constraints.
It is an interpretation/reconstruction aid, not permission to invent
correspondence or issue an RS3 signature for incomplete chemistry.

### 5.3 Architectural cleanup

The desired end state is one canonical structure-backed recommendation path.
That cleanup is not finished:

- expert rule recommendation remains a separate deliberately limited path;
- weak-label retrieval remains a transitional lower-evidence path;
- legacy `chemtools` and legacy application paths still exist outside the new
  package dependency graph;
- duplicate or transitional paths should be removed only after parity,
  evaluation, migration, and user-facing replacement are demonstrated.

## 6. Required next sequence

Do not add large new chemistry families or tune against the untouched test
before completing the following gates.

### Gate 1: Freeze a current validation baseline

1. Run the full deterministic suite.
2. Validate taxonomy and condition definitions.
3. Validate the current sharded conversion and persisted index.
4. Record current schemas, definition hashes, source checksums, index ID,
   coverage, admission reasons, retrieval distribution, and performance.
5. Regenerate any evaluation artifact made with older schemas or definitions.

**Pass condition:** all machine artifacts are reproducible, mutually compatible,
and tied to the exact source and code state.

### Gate 2: Generate a current blind chemist review

Build a stratified packet from the current development/validation artifact.
Include:

- verified, review, rejected, and fallback cases;
- common and rare transformation classes;
- mapped, reconstructed, global-correspondence, incomplete, and conflict
  evidence;
- intramolecular and multi-event examples;
- high-confidence family and unnamed reactions;
- exact, relaxed, environment, broad-edit, and abstaining retrieval results;
- deliberately unsuitable condition controls.

The reviewer should see structures, highlighted reaction centers, concise
edits, recipes, precedents, support, and cautions, but not the answer key.

**Pass condition:** an independent chemist signs the hash-bound adjudication
summary with no unresolved systematic defect.

### Gate 3: Resolve disagreements

Classify every disagreement as:

- implementation defect;
- chemistry-definition defect;
- source-data defect;
- benchmark/review defect;
- ambiguous chemistry;
- supported limitation or required abstention.

Confirmed defects require a focused code/definition change, named regression
test, and regenerated affected artifacts. Do not change an expected snapshot
merely to match new output.

**Pass condition:** all systematic disagreements are resolved or explicitly
converted into supported abstention/review behavior.

### Gate 4: Run the untouched evaluation

Use reference- and canonical-reaction-connected partitions, then report strict
scaffold-disjoint, source-disjoint, and forward-time results where possible.
Compare family-only, generic-only, hybrid, transformation-prior, and current
baseline behavior.

Minimum reporting:

- query and recommendation coverage;
- abstention rate and reviewed correctness;
- seen-recipe top-1 and top-k recovery with denominators;
- variant recovery where the variant exists in training;
- fallback-level distribution;
- hard-incompatible recommendation count;
- independent-reference support;
- yield error only with usable-outcome counts;
- explanation completeness;
- results by transformation and evidence class.

**Pass condition:** zero unexplained hard-incompatible recommendations; hybrid
retrieval expands justified coverage without material established-family
regression; the review and numerical evidence support release.

### Gate 5: Convert and validate the full source corpus

1. Freeze contract and definition versions.
2. Convert the full source corpus in restartable deterministic shards.
3. Preserve source checksums and all artifact provenance.
4. Build the production index.
5. Run conversion and index integrity checks.
6. Re-run coverage, performance, and release validation.

**Pass condition:** no failed or mixed-version shards, no duplicate
observations, current schemas only, acceptable runtime/memory, and all machine
plus human release gates pass.

### Gate 6: Consolidate public paths

After parity and release evidence:

1. make the generic structure-backed workflow the canonical application path;
2. retain expert rules only as explicitly identified reviewed protocols or
   integrate them as provenance-preserving overlays;
3. retire weak-label or legacy paths when their remaining use cases have a
   validated replacement;
4. remove migrated duplicate logic rather than maintaining parallel behavior;
5. update application and package READMEs to point to these two consolidated
   documents.

## 7. Chemistry expansion after the release gates

Add new chemistry by definite edit shape and validated site requirements, not
by dataset filename or source reaction name.

A chemistry increment must include:

1. molecular site/interface definitions;
2. a bounded connectivity rewrite or a demonstrated need for a new instruction;
3. product projection and atom-conservation rules;
4. positive, negative, ambiguous, conflict, and invariance tests;
5. mapped and exact-reconstruction evidence where applicable;
6. condition compatibility and explanation updates;
7. a stratified dataset-impact report and chemist review.

Suggested progression:

1. close high-frequency unsupported simple substitutions, cleavages, and
   localized redox gaps;
2. validate condensation and olefination with explicit omitted-fragment
   handling;
3. add insertion, extrusion, exchange/metathesis, and migration only after
   their atom accounting is explicit;
4. add annulation, cycloaddition, rearrangement, and cascades only after
   connected multi-edit validation is robust.

Family assignment remains optional. A new family must not create a separate
converter or recommender.

## 8. Testing and operational checks

Run the complete suite before handing off any change:

```powershell
pytest -q
```

Useful focused suites:

```powershell
pytest -q tests/reactive_taxonomy
pytest -q tests/condition_registry
pytest -q tests/condition_recommender
```

Validate definitions:

```powershell
python -m reactive_taxonomy.cli validate
python -m condition_registry.cli validate
```

Validate current persisted artifacts:

```powershell
python -m condition_recommender.conversion_integrity_cli `
  datasets/literature/shard_manifest.json

python -m condition_recommender.generic_index_integrity_cli `
  datasets/literature/generic_index.json.gz
```

For a new release candidate, generate and adjudicate the chemist packet using
the current CLI help and repository README. Paths should be release-specific;
do not reuse an old answer key or signed summary after any identity-bearing
schema or definition changes.

Required chemistry regressions include:

- Suzuki C–C and sp2 C–N/O/S;
- acyl N/O/S formation;
- additions, eliminations, and localized redox;
- mapped and unmapped equivalents;
- mapped unknown-family reactions;
- valid, invalid, partial, and conflicting atom maps;
- incomplete product provenance;
- single- and multi-event reactions;
- intramolecular/intermolecular distinctions;
- partner-order and SMILES-serialization invariance;
- deterministic signature and recipe IDs;
- hard compatibility before similarity;
- missing features not treated as matches;
- stale artifact rejection.

## 9. Documentation reconciliation

| Older document | How to use it now |
| --- | --- |
| [`../type_agnostic_reaction_recommendation_implementation.md`](../type_agnostic_reaction_recommendation_implementation.md) | Primary historical architecture and phase plan. Its Phase A “immediate work” is implemented; current status and versions are superseded by this document. |
| [`../reaction_featurization_workflow.md`](../reaction_featurization_workflow.md) | Detailed near-current technical walkthrough. Retain for code-level flow; read CEG1 as an internal shadow observation and RS3 as public identity. |
| [`../connectivity_first_reaction_grammar_v2.md`](../connectivity_first_reaction_grammar_v2.md) | Current summary of the V2 connectivity executor. Its own note correctly says environment/signature identity advanced to RS3. |
| [`../connectivity_first_reaction_grammar_design.md`](../connectivity_first_reaction_grammar_design.md) | Design rationale and migration history. Shadow/authority/legacy-operator stages are historical; V2 is the active executor state. |
| [`../context_aware_reactivity_descriptors_implementation_plan.md`](../context_aware_reactivity_descriptors_implementation_plan.md) | Completed descriptor migration record and useful detailed contract. The phase sequence is history; typed profiles are already active. |
| [`../aromatic_reactivity_descriptor_proposal.md`](../aromatic_reactivity_descriptor_proposal.md) | Historical rationale for the aromatic slice, now generalized and implemented by typed profiles. |
| [`../recommend_generic_conditions_implementation_plan.md`](../recommend_generic_conditions_implementation_plan.md) | Most detailed recommendation/evaluation plan and historical pilot metrics. Its implementation-slice versions and artifact paths may lag current code. |
| [`../phased_implementation_and_evaluation_plan.md`](../phased_implementation_and_evaluation_plan.md) | Early evaluation record. Phase 1/2 status is not the current whole-project status. |
| [`../reactive_handle_taxonomy_design.md`](../reactive_handle_taxonomy_design.md) | Foundational V1 reactive-site design. Its site-first principle remains; old API, mapping, family-routing, and milestone proposals are superseded. |

The old detailed files should remain available until any still-useful chemistry
details are either migrated or intentionally archived. New status changes
should be made here first, with links to machine-generated reports rather than
copying transient metrics into several plans.

Some older files still link to removed
`connectivity_first_reaction_grammar_phase*_implementation*.md` documents.
Those links describe migration history, not missing active contracts; the V2
summary and current code supersede them. Do not restore the deleted phase files
only to preserve an obsolete implementation sequence.

## 10. Release definition of done

The type-agnostic system is ready for production use only when:

1. structure-backed unknown-family reactions are represented and retrievable;
2. conflicting or incomplete evidence cannot enter verified retrieval;
3. hard chemistry compatibility always precedes similarity;
4. recipes preserve identity, role, stage, quantity, and provenance;
5. support distinguishes rows, reactions, series, references, and corpora;
6. unsupported queries produce a typed abstention or clearly isolated fallback;
7. recommendations cite precedents and explain matches, mismatches, fallback,
   missing data, and uncertainty;
8. current artifacts pass version, checksum, duplicate, and integrity checks;
9. leakage-safe untouched evaluation passes with zero unexplained hard
   incompatibilities;
10. an independent chemist signs the current blind-review adjudication with no
    unresolved systematic defect;
11. the full source corpus is converted and indexed reproducibly;
12. the full deterministic test suite passes;
13. obsolete duplicate recommendation/conversion paths have a documented
    removal or containment decision.
