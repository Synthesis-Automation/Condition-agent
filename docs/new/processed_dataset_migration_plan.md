# Unified processed dataset migration plan

Status: implementation plan; runtime changes and corpus regeneration have not started.
Prepared: 2026-10-05.

Deliver one complete processed corpus under `datasets/processed_datasets`, with
one build workflow for condition retrieval, shared-core retrieval, fragment
search, route context, and forward/retrosynthesis artifacts. Apps and agents
will resolve a single release manifest. Remove production Full/Compact modes,
retain all observations and their evidence, and make routine query responses
concise. Fragment search will cover starting materials and products, with the
matched side explicit in every result.

This plan follows the [primary implementation roadmap](type_agnostic_reaction_recommendation_implementation.md)
and the [intermediate migration](../../raw_datasets/intermediate_dataset_migration.md).
It proposes one coordinated migration, implemented through resumable stages.
It does not declare chemistry review or untouched-evaluation gates complete.

## Scope and decisions

- Convert all physical observation inputs from `datasets/intermediate_datasets`.
  Keep verified, review, and rejected observations with separate capability and
  admission statuses. Retention never implies eligibility for every tool.
- Remove production sampling, Full/Compact selectors, directory conventions,
  mode parameters, and automatic fallbacks to old libraries. Preserve explicit
  development samples and leakage-controlled evaluation splits outside the
  production release.
- Keep canonical nested JSON as the primary artifact, partitioned into shards
  and shared objects. SQLite artifacts are derived query stores. CSV remains an
  optional export.
- Preserve full evidence while reducing repeated serialization and default
  agent output. Avoid parallel abbreviated and full corpus products.
- Include both reaction sides in fragment discovery. Preserve the current
  product default for synthesis discovery; expose starting-material and either-side
  choices in the GUI, API, CLI, and agent tools.
- Keep stock availability and weak-label datasets explicitly separate external
  capabilities. A literature occurrence does not establish purchasability.
- Keep `reactive_taxonomy` responsible for chemistry, `condition_registry` for
  substance/recipe identity, and `condition_recommender` for corpus conversion,
  storage, indexing, and retrieval. Existing planning packages own their operator
  and route contracts. The application orchestrates these packages without
  moving chemistry rules into GUI code.

## Investigation baseline

The following are observations from local files on 2026-10-05, not performance
targets or projected savings.

| Item | Observed value |
| --- | --- |
| Existing full converted corpus | 660,190 observations |
| Existing compact converted corpus | 117,232 observations |
| New intermediate physical observations | 2,054,437 |
| Preserved source routes | 976,597 |
| Algorithmic abstractions outside physical conversion | 1,549,008 |
| Full combined canonical gzip | 8.29 GB |
| Full recommendation SQLite index | 5.43 GB |
| Full shared-core SQLite index | 4.17 GB |
| Entire literature directory, including copies and backups | Approximately 80.2 GB |

Sizes use decimal GB. A sample of the first 60 HiTEA canonical rows averaged
85.6 KB of uncompressed JSON; all 60 repeated the same core at the top level and
inside `reaction_observation`. Sixty positions spread across the recommendation
index averaged 46.1 KB per serialized indexed row. These bounded samples do not
establish whole-corpus savings.

The existing converter keeps the first 200 source rows plus 15 percent of the
remainder in Compact mode. Saved shards and merged records coexist. Catalog and
index construction reread canonical records. The agent precedent accessor
returns complete indexed payloads, while procedure/reference accessors scan gzip
catalogs. Fragment search currently indexes product components only and embeds
discovery details and optional procedures in its own database.

## Published layout and release identity

Use a small active manifest at the stable root and immutable release directories.
This adds versioning for publication and rollback, not multiple coverage modes.

```text
datasets/processed_datasets/
  manifest.json                         active release reference and digest
  releases/<release_id>/
    manifest.json                       immutable artifact inventory
    records/
      observations-*.jsonl.gz
      chemistry-*.jsonl.gz
      recipes-*.jsonl.gz
      evidence-*.jsonl.gz
      route_memberships-*.jsonl.gz
      routes-*.jsonl.gz                  validated route trees, where supported
    indexes/
      generic_index.sqlite
      generic_index.shared_core.sqlite
      catalogs.sqlite
      fragment_index.sqlite
    operators/
      retrosynthesis/
      forward/
      composite/
    reports/
      coverage.json
      validation.json
      performance.json

results/dataset_builds/<build_id>/        checkpoints and temporary outputs
```

The release manifest records source fingerprints, schema and definition versions,
RDKit and mapper identity where relevant, artifact hashes, observation counts,
capability scopes, dependencies, and validation status. Store relocatable paths
relative to the release. Reject paths escaping the release for internal artifacts;
declare external capabilities separately. Derive identity from stable content
and contracts, not machine-specific absolute paths or build timestamps.

Build in a staging directory on the destination volume. Validate and close files,
finalize the immutable release, then atomically replace the small active manifest.
A process resolves and pins one release for its lifetime or investigation; it
must not mix files across a publication. Rollback changes the active manifest
to a retained compatible release. Partial builds never become active.

## Canonical storage contracts

| Object | Content and identity requirements |
| --- | --- |
| Observation | Stable observation ID, source aliases, reaction identity, raw and canonical structures, experiment-specific outcomes/stages, chemistry and recipe references, admission/capability statuses, warnings, evidence references |
| Chemistry analysis | Parsed observations, validated atom correspondence, signatures, cores, descriptors, interpretations, conflicts and candidates; share only for equivalent complete analysis inputs and versions |
| Recipe | Canonical resolved identities and process structure; preserve observation-specific raw identifiers, role evidence, amounts, stages and conflicts through separate evidence links |
| Evidence | Procedures, reference identities, original values needed for audit, supplied mappings, diagnostic detail and source locators; retrievable by stable IDs |
| Molecule occurrence | Molecular identity plus observation ID, reaction side, component index, original structure and atom-order mapping |
| Route occurrence | Source route/tree/subtree membership, source reaction aliases and observation links, connectivity status, source evidence and any validated tree reference |

Do not merge experiments because their reaction SMILES or recipe match. Deduplicate
identical source observations under explicit rules while retaining all aliases.
Report overlap between route-derived patent steps and other USPTO sources; do not
count duplicate publications or route repetitions as independent evidence.

A shared chemistry object must not erase source-specific atom indices, atom maps,
warnings, or conflicting interpretations. Keep occurrence mappings alongside
shared identities. Do not change canonical IDs merely to compact storage. Where
a schema change requires new IDs, produce an explicit migration map.

Define typed serialization and hydration contracts with versioned schemas.
Hydration reconstructs all required semantic fields from canonical objects;
redundant legacy fields may disappear only after reader migration and parity
checks. Record the intentional removals. Definition bundles are stored once
with version references from records. Keep full raw source archives upstream,
but include the evidence needed for ordinary app/agent inspection in the release.

Canonical shards are the sole persistent corpus export. Replace mandatory merged
gzip files with manifest iteration. Derived indexes retain only fields needed
for filtering, scoring, rendering summaries, and efficient evidence lookup.
The catalog SQLite store provides indexed retrieval of observations, shared
chemistry, recipes, publications, procedures, and route memberships. Detailed
compressed JSON objects may be materialized there once for random access; avoid
copying them into every search index. A gzip line number alone is not random access.

## Fragment discovery on both reaction sides

Index valid connected components from the first and third fields of
`reactants>agents>products`. Do not silently classify the middle agents field as
starting materials. The first field is a reported reactant-side occurrence;
chemical consumption and partner roles require additional evidence.

Store each normalized searchable molecule once, with occurrences keyed by
`(observation_id, side, component_index)`. Preserve stereochemistry, isotope and
charge distinctions under a versioned normalization policy. Store mappings from
searchable atom order back to the original side/component atom references.
The same molecule on both sides has two occurrences, not one overwritten link.

The proposed public request adds `search_side` with values `product`, `reactant`,
or `either`; default `product`. UI labels are Products, Starting materials, and
Either side. Required hit fields include:

- `matched_side`, `component_index`, and occurrence/observation IDs;
- `matched_molecule_smiles` and original source structure;
- `match_extent`: `whole_molecule` or `substructure`, derived from structural
  matching under the declared stereo/isotope/charge policy;
- query-to-source atom mapping and the highlighted side;
- reference/source links, reported reaction, and available procedure references;
- relationship evidence status, warnings, and separately qualified admission;
- matched sides and occurrence counts when multiple occurrences are grouped.

Examples of display wording are "Reported as starting material — exact compound
match" and "Product contains queried fragment — substructure match". A fragment
match does not establish that the standalone queried compound was made or used.
Even an exact starting-material occurrence establishes reported use, not how the
compound was prepared. Expose the cited source so its preparation or procurement
can be investigated. Missing references remain explicitly missing.

Retain the established product-relative relationship vocabulary and evidence
rules. For starting-material hits, expose reported use first; add transformation
relationships only where supplied mapping and validated edits support them.
Never label a reactant match "constructed" merely because a matching product
exists. Unresolved correspondence must not prevent valid structure discovery.

Fragment-guided retrosynthesis and transfer remain restricted to eligible
product-side construction evidence. Revalidate selected occurrences on the
server. Observation IDs alone cannot authorize transfer when the selected match
was on the reactant side. Agent summaries, exports, cached results, UI highlights,
and transfer requests must retain side information end to end.

Use one unique-molecule RDKit library with side-indexed occurrence links, or
equivalent side-specific search partitions over shared molecules. Apply side
selection before query caps so hits on the other side do not exhaust the search
budget. Count molecules, occurrences, observations, and publications separately.
Version the product-specific result names and limits into side-neutral contracts.
Keep deterministic ordering, explicit truncation and timeout status, and no
unsupported absence claims after partial searches.

## Unified processing workflow

```text
Validate and fingerprint intermediate sources
                 |
                 v
Convert all observations and persist canonical shared objects
                 |
       +---------+----------------+-------------------+
       |                          |                   |
       v                          v                   v
Condition and shared-core   Molecule occurrence   Route membership and
indexes and catalogs       and fragment index    qualified route conversion
       |                          |                   |
       +--------------------------+-------------------+
                                  |
                                  v
                    Eligible forward and retro operators
                    and supported composite artifacts
                                  |
                                  v
                  Validate coherent release and publish
```

One application action and one CLI workflow orchestrate this dependency graph.
Internal stages remain restartable. No stage owns a private conversion path or
recomputes chemistry merely to produce another artifact format. Persist shared
features and mapping assessments once with their full input/definition keys.
Use bounded batches and disk-backed deduplication; benchmark the in-memory RDKit
library and molecule dictionary before selecting a large-corpus partition strategy.

Checkpoint completed shards and artifact partitions by dependency hashes. Resume
only compatible checkpoints. Invalidate affected descendants when chemistry,
conditions, fragment policy, source data, or builder versions change. For example,
a condition-only change need not reparse molecules, but must rebuild recipe-bearing
outputs; any reused artifact must declare and validate its actual dependencies.
Check available disk against measured staging and publication requirements.

The production build covers every declared physical observation without
`max_records` or Compact sampling. Each capability reports indexed, unavailable,
review-only, and failed counts as applicable. Missing conditions can leave an
observation searchable by structure while excluding it from condition transfer.
An unexpected stage failure blocks publication of the requested release; an
explicit, justified unsupported capability is reported rather than fabricated.

## Multistep routes and other inputs

Read every released route wrapper for membership/source context, joining steps
through original reaction aliases. Keep this distinct from deduplicated physical
step observations. Preserve conflicting or ambiguous joins for review.

Reuse the existing [route conversion contract](route_tree_contract_and_conversion.md)
for routes that pass its curation requirements. Adapt the new wrapper source to
that contract rather than creating another route schema. Reconstruct connectivity
from validated structures and the declared target; never trust step-array order.
Record observed steps separately from inferred route connections. Every input
route receives a validated-tree reference or an explicit unresolved/rejection
status. Retain source memberships for routes that cannot form a validated tree.
Composite operators consume only appropriately qualified routes and step evidence.

Keep algorithmic abstractions archived upstream and excluded from physical
reaction counts, fragment occurrence evidence, and observed operator training.
Record their exclusion in coverage reports. Do not reingest the superseded
`uspto_original.csv` alongside extracted route steps.

Reconcile old source IDs and duplicate aliases against the new corpus. Resolve the
old Organic Syntheses example explicitly: include its existing example source
under the raw-source contract or record its retirement. Do not silently discard
unmatched legacy records. Stock and weak-label sources keep their own manifests
and evidence scopes; the release references them as external capabilities when used.

## Application and agent output

Default precedent results contain structures, side and match scope where relevant,
reported conditions/outcome, essential evidence quality, source/publication IDs,
material warnings, and references to details. Preserve the distinction between
an observed condition and a proposed recipe or predicted product.

Extend existing workspace summaries and accessors with explicit field selection,
pagination, and byte budgets. Detailed chemistry, raw source values, and procedures
are fetched by ID from the shared catalog. A bounded response must report omitted
sections and continuation references. Keep blocking conflicts visible even when
verbose diagnostics are omitted. Do not use generated prose as canonical evidence.

The web app and agent must resolve the same release and verify dependency hashes.
Workspace investigations pin immutable release identities and external dependencies.
Preserve old saved evidence as historical artifacts without maintaining automatic
runtime fallback to `datasets/literature`. Document that new investigations use
the new release; do not rewrite old conclusions under changed evidence.

## Implementation sequence and exit criteria

Each stage is part of the same migration. Checkboxes are intentionally incomplete.

| Stage | Work | Exit criterion |
| --- | --- | --- |
| 1 | Freeze current code/definitions/artifact inventory, counts, baseline queries and release-gate status; select representative development pilots | Reproducible baseline with old/new source reconciliation inputs and no untouched-set tuning |
| 2 | Define release, normalized storage, hydration, occurrence, and query contracts; specify ID/schema migrations | Validated schemas and fixtures preserve evidence, identity and uncertainty |
| 3 | Refactor canonical writer/readers, shared catalogs, and index payloads; remove repeated corpus serialization | Semantic reconstruction and unchanged-input chemistry/retrieval parity pass |
| 4 | Implement both-side fragment indexing, explicit hit roles, side-aware limits, lookup joins and transfer guards | Product regressions pass; reactant/either queries and UI/export labeling pass |
| 5 | Add the resumable build orchestrator, route adapter, and operator dependencies | One pilot build produces all declared artifacts; cancellation/resume and changed-input invalidation pass |
| 6 | Migrate GUI, CLI, API, agent, web clients and planning readers; remove Full/Compact contracts | All active consumers use one manifest and no old-path automatic fallback remains |
| 7 | Run representative performance and chemistry parity checks, full test suite and required independent evaluation gates | Engineering checks pass; roadmap prerequisites for full-corpus conversion are satisfied or an explicit user-authorized exception is recorded |
| 8 | Build and validate all current physical inputs and route records, then publish the release atomically | Coverage reconciles, all required artifacts agree, app and agent smoke checks pass |
| 9 | Remove obsolete conversion paths and document retirement/retention of old files and checkpoints | One canonical runtime path remains, rollback is documented, no necessary historical evidence is deleted |

Do not infer an exception to the primary roadmap from this planning request.
Freeze the current validation baseline, generate the blind chemist review packet,
resolve disagreements without tuning on untouched data, and run the untouched
evaluation in the prescribed order before the full-corpus release work. Engineering
pilots support migration decisions but do not satisfy independent scientific gates.

## Files and ownership to change

This is a starting inventory; stage 1 must search for remaining active consumers.

| Area | Existing files or paths |
| --- | --- |
| Canonical schema and conversion | `condition_recommender/models.py`, `conversion/generic.py`, `conversion/sharded.py`, `conversion/artifacts.py`, `corpus_io.py` |
| Condition and core indexes | `condition_recommender/generic_indexing.py`, `sqlite_indexing.py`, `shared_core_index.py` |
| Fragment chemistry and policy | `reactive_taxonomy/fragment_search.py`, `definitions/fragment_search.v1.json` |
| Fragment storage and queries | `condition_recommender/fragment_index.py`, `fragment_search.py` |
| Planning inputs and builders | `core_retrosynthesis/sources.py`, existing operator builders, route conversion/curation, fragment transfer and composite consumers |
| Desktop workflow | `app/reaction_converter_gui.py`; inspect preprocessor and featurizer integration for shared configuration without duplicating domain rules |
| Recommendation entry points | `condition_recommender/generic_api.py`, `generic_recommend_cli.py`, `chem_coworker/service.py` |
| Web API and UI | `app/web_api/contracts.py`, `runtime.py`, `references.py`, CLI options; `web/reaction_recommender/src` API types, hooks and fragment/planning components |
| Scientific workspace | `chem_coworker/scientific_workspace/adapters/operations.py`, `source_catalogs.py`, `fragment_search.py`, views and artifact validation |
| Configuration and examples | `examples/ai_native/artifacts.local.example.json`, active build/query examples, environment defaults |
| Documentation | Primary roadmap, intermediate migration guide, `docs/AI-native/readme.md`, CLI/API and fragment search guidance |

New orchestration and release-access code belongs in the owning standalone package
or thin application layer as appropriate. Do not add dependencies from chemistry
or registry packages onto application/planning code. Audit exported APIs, request
models and tests for `library_mode`, not just directory strings.

## Validation and acceptance

Storage tests cover deterministic IDs and ordering, hydration of all retained
evidence, conflicting duplicate IDs, multi-stage conditions, provenance aliases,
relocation, schema mismatch, missing objects, corrupt artifacts, and release
dependency verification. Retain the existing Suzuki, C-N, C-O, C-S, ordering,
unknown-family mapping, invalid-map and ambiguity/conflict regressions.

Fragment tests cover product-only and reactant-only compounds, the same molecule
on both sides, repeated components, multiple products, salts/disconnected sides,
stereochemistry, isotopes, charges, invalid components, missing maps, ambiguous
correspondence, exact versus substructure matches, ring-preserving queries,
side-specific atom highlights, and middle-field exclusion. Confirm that an
unresolved reactant hit remains discoverable and cannot silently become product
construction evidence. Test side filtering before result caps and truthful
partial-search counts, missing references, and source retrieval by observation ID.

Route tests cover linear and branched routes, repeated molecule occurrences,
ambiguous joins, cycles, unreachable steps, and preservation of unresolved source
memberships. Keep generated abstractions out of observed evidence.

Integration tests cover interrupted/resumed builds, hash invalidation, old/new
schema rejection, immutable reader pinning, publication and rollback, all public
clients, and fragment transfer guards. Run `pytest -q` before handing off runtime
changes, plus the relevant web type/build and end-to-end checks. Verify API startup
and desktop conversion workflow where their contracts change.

Use a frozen development panel spanning HiTEA, USPTO conditions, combined RDF
literature, and route-derived steps, including large/ambiguous/rejected records.
Compare old and new processing of the same inputs under the same definitions;
measure new corpus coverage separately so source expansion does not mask a storage
regression. Preserve publication/route split isolation and keep untouched data out
of implementation tuning.

Measure disk bytes by artifact, peak build RAM and temporary disk, elapsed stage
times, cold/warm lookup latency, and default response bytes at representative
percentiles. Establish numeric budgets from this pilot before the production run;
do not claim a reduction percentage from the small investigation samples. Require
smaller representative default responses with preserved warnings, indexed detail
lookups, and no loss of eligible results on unchanged-input regression fixtures.

At full build time, reconcile the current manifest's expected 2,054,437 physical
observations, including 1,414,143 extracted route steps, with every conversion
outcome. Reconfirm counts if raw inputs change. Route records must likewise
reconcile to the current 976,597 inputs. Report parse failures and per-capability
exclusions explicitly; no silent skips. A complete build is not a claim that all
records are chemically verified or all routes reconstruct successfully.

## Cutover and completion

Publish only a validated, internally consistent release. Switch default consumers
together and create a new scientific-workspace baseline. Preserve a rollback
release and any artifacts pinned by historical investigations. Remove obsolete
runtime code and mode selectors after parity; do not leave a permanent compatibility
conversion path. Retire old intermediate/literature data only after reconciliation
and retention checks, with deletion performed as a separately explicit operation.

The migration is complete when one command/GUI action builds the complete declared
corpus, every required app/agent capability uses the same manifest, starting-material
hits are visibly and structurally distinguished from product hits, concise responses
can retrieve full evidence by ID, all required validation is recorded, and active
runtime paths no longer depend on Full/Compact or `datasets/literature`.
