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
| Edit evidence | Implemented | Typed provider candidates from mapping, exact single/multi-event reconstruction, conservative scaffold, bounded global, and fragmented-scaffold correspondence; distinct ambiguous edit hypotheses; conflicts and H/charge/stereo observations |
| External atom mapping | Optional review/query integration | RXNMapper proposals are structure-validated, reconciled against internal hypotheses, persisted with model provenance, and never admitted as verified precedents |
| Product completeness | Implemented | Verified/incomplete/unresolved accounting, observation-only product-origin gaps, and typed fragment-source requirements |
| Reaction signatures | Implemented | Deterministic RS3 L0–L4 signatures, events, topology, profiles, spectators, unknown-family support |
| Reaction minimization | V2 integrated, calibration pending | Template-free mapped-edit core projection with atom transitions, remote subgraphs, typed attachment ports, robust shape keys, review export, verified-precedent index maps, and conservative query routing |
| Reactivity descriptors | Implemented and active | Typed context-aware profiles are the sole active environment path |
| Condition registry | Implemented, curation incomplete | Conservative identity resolution, contextual roles, RCORE1/RCR1 recipes, stages, provenance |
| Generic conversion | Implemented | Nested canonical records, independent quality dimensions, review exports, sharding, restart/integrity checks |
| Generic index | Implemented | Version-checked persisted index with signature, reaction-core, exact partial-transformation, environment, family, fallback, recipe, and reference keys |
| Generic retrieval | Implemented pilot | Explicit verified-signature ladder, robust reaction-core shape tier, conservative unsigned-query core and all-hypothesis routes, independent-support thresholds, hard compatibility, similarity, and reference-aware recipe aggregation |
| Source-supported partial transformations | Implemented, review-qualified | Product-observed attachment replacement may retrieve only exact partial transformations whose precedent conditions have a curated source capability |
| Unverified-query fallback | Implemented, conservative | Separate structure-derived fallback for other unsigned queries; not represented as verified edit retrieval |
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
| Reaction analysis | `3.5` |
| Reaction signature / ID namespace | `3.0` / `RS3` |
| Reaction core projection / algorithm | `2.2` / `reaction_core_projection.v8` |
| Taxonomy identity manifest | `3.0` |
| Connectivity site interface | `2.0` |
| Connectivity rewrite | `2.0` |
| Typed reactivity profile | `1.0` |
| Reaction fallback descriptor | `1.3` |
| Resolved condition recipe | `1.2` |
| Recommendation record | `3.9` |
| Generic converter definition | `generic_conversion.v3.3` |
| Generic sharded converter definition | `generic_sharded_conversion.v2.0` |
| Concise reaction review | `2.7` |
| Recommendation artifact workflow | `1.1` |
| Generic persisted index | `2.5` |
| Generic recommendation result | `2.4` |
| Reaction correspondence definitions | `2.3` |
| Generic retrieval definition | `1.8` |
| Reaction-core retrieval policy | `reaction_core_retrieval.v2@1.0` |
| Generic admission policy | `generic_admission.v2.0` |
| Fragment-source capabilities | `fragment_source_capabilities.v1@1.0` |
| Fallback retrieval definition | `1.4` |

Do not copy this table into executable code. The constants and definition files
remain authoritative, and stale artifacts must fail validation rather than
silently mixing chemistry identities.

### 2.3 Current local literature artifact

The checked local artifact under `datasets/literature/` was most recently built
from the bounded `examples/small_300` collection, not the full source corpus.
Its pre-source-support report records:

| Measurement | Count |
| --- | ---: |
| Source files / converted rows | 19 / 5,510 |
| Generic-index rows | 2,603 |

The sharded conversion report records zero failed shards and zero duplicate
observations. The artifact predates reaction-core projection `2.0`, converter
`v2.9`, and index `2.5`; current code therefore rejects it until it is
regenerated. This is deliberate because reaction-core identity and index maps
change artifact identity and retrieval coverage.

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
| Ambiguous queries with typed edit hypotheses | 205 / 205 |
| Queries whose hypotheses share one anonymous edit prototype | 116 / 205 |
| Reference-disjoint queries with a verified all-hypothesis consensus neighbor | 195 / 205 |
| Consensus links across precursor modes | 117 / 205 |
| Anonymous edit prototypes among signatures | 3 |
| Reference-disjoint signed queries with exact L3 support | 54 / 58 |
| Reference-disjoint signed queries with edit-graph support | 58 / 58 |
| Reference-disjoint queries linked across precursor modes | 58 / 58 |

This is a successful retrieval POC but not a solved atom-mapping POC. The
approximate edit graph connects the separate-reactant and preformed-hydrazone
views for every currently signed query. All 205 ambiguous cases now preserve
their distinct edit alternatives as typed `ReactionEditHypothesis` objects.
For query use only, each hypothesis is independently compared with verified
precedents; a precedent survives only when it is compatible with every
hypothesis, and ranking uses the worst hypothesis match. This recovers a
reference-disjoint consensus neighbor for 195 cases without promoting any
ambiguous dataset row into the verified index. The remaining ten abstain.

Competing nitrogen origins or symmetry-related mappings remain genuinely
ambiguous from the unmapped source structures. The query-consensus result is
therefore explicitly review-required; those 205 cases must not be bulk-admitted
without curated mapping or an independently validated correspondence provider.

The representative cyclohexanone/4-fluorophenylhydrazine query now retains two
`REH1` alternatives. Against the 45-row verified Fischer POC index, both
hypotheses converge on 45 compatible precedents from 19 independent support
units and the API returns five review-required recipes through
`edit_hypothesis_consensus`. This is an executable end-to-end demonstration of
query coverage; it does not resolve which atom correspondence is correct.

The normal generic converter was also run over all 542 source rows. It produced
59 signed observations and 45 index-eligible precedents; the compact report is
`results/fischer_indole_poc_conversion/conversion_report.json`. A representative
accepted preformed-hydrazone query used the persisted index and selected
`edit_graph_neighbors`: 41 compatible candidates from 18 independent
references. Its exact-signature through L3 pools contained only one independent
reference, so the new tier supplied genuinely broader support. The result
retained `named_family=None`, disclosed the intramolecular/intermolecular scope
mismatch, and warned that the exact edit signatures differed.

### 2.5 Template-free reaction-minimization POC

`benchmarks/reaction_core_sample_report.py` evaluates
`examples/sample_reactions.csv` without passing source reaction types into
mapping, minimization, keys, or labels. The current V2 artifacts are
`results/reaction_core_sample_report.csv` and
`results/reaction_core_sample_report_summary.json`.

The 2026-07-30 run reports:

| Measurement | Count |
| --- | ---: |
| Rows / valid reactions | 446 / 445 |
| Baseline RS3 signatures | 341 |
| Available reaction-core projections | 445 |
| Core coverage where the baseline had no signature | 104 |
| One-event / two-event cores | 421 / 24 |
| Unique center-transition keys | 258 |
| Unique robust shape keys | 351 |
| Rows in repeated center clusters | 241 |
| Rows in repeated shape clusters | 141 |
| Rows in mixed-source-label center clusters | 146 |
| Rows in mixed-source-label shape clusters | 58 |
| No-op primary-center warnings | 6 |
| Unresolved remote-continuity warnings | 3 |

The shape key substantially reduces conflation from the center-only key by
including participant handle/site evidence, retained remote shape, and event
count. A hard-negative regression confirms that a vinyl Suzuki and a Heck
reaction may share `RCS2` while receiving different `RSH2` keys. This is why
`RCS2` is diagnostic only.

The sample contains an acetonide-protection reaction for which the baseline has
no signature. V2 produces the grammar-independent minimized label
`C(R)2(=O) → C(R)2(O-R)2` and derives the remote classes `alkyl` and
`ring_aliphatic` from the cut molecular graph. The supplied mapped
benzaldehyde/methanol acetal regression similarly produces
`C(H)(Ar)(=O) + O(H)(R) → C(H)(Ar)(O-R)2`; mapper atom-origin alternatives
keep the same shape and center keys while exact and typed identities preserve
provenance. Single-center labels include a distinct external reactant when its
atom forms a bond to the primary center.

These measurements demonstrate coverage and improved cluster discrimination,
not recommendation accuracy. Of the 445 cores, 434 are external-mapper
proposals and therefore remain review evidence. The core retrieval policy is
still marked `poc_review_only_not_calibrated`; verified-precedent coverage,
false-neighbor review, and blind condition review remain release gates.

#### Graphical-core feasibility POC

`benchmarks/click_reaction_core_graphic_poc.py` exercises graphical
minimization on a large alkyne/azide cycloaddition without passing the source
reaction name into featurization. The mapped projection contains five active
atoms and one connected event. Two unchanged retained substituents are cut at
their attachment bonds, paired across the reaction by mapped identity, and
rendered as `R1` and `R2`. The resulting scheme is an alkyne plus azide forming
a substituted five-membered triazole.

The versioned renderer lives in `visualization`, while atom transitions,
continuity, remote classes, exact fragments, and typed attachment ports remain
owned by `reactive_taxonomy`. This keeps drawing policy out of core chemistry
and out of the application layer. Only retained remote subgraphs may be
abstracted; the graphic result always includes a placeholder-to-fragment
legend. Singleton retained heteroatoms remain explicit. Multi-port ring
remainders surrounding one active atom absorb that atom into a scaffold
placeholder while retaining its active bonds; other multi-port topologies use
one boundary node per port. This avoids parallel dummy bonds and supports
transformations at atoms embedded in fused rings.

Artifacts are:

- `results/click_minimized_reaction_poc.svg`;
- `results/click_minimized_reaction_poc.png`; and
- `results/click_minimized_reaction_poc.json`.

The current external mapping confidence is `0.29134308160116135`, with
`EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW`. Therefore this POC establishes
drawing feasibility and UI usefulness, not validated mapper accuracy. The
featurizer keeps resolved-reaction shadow mapping off by default and exposes it
as an explicit option for generating this review graphic.

### 2.6 External RXNMapper Fischer POC

The optional offline dependency is pinned in `requirements-mapping.txt`.
`reactive_taxonomy.external_atom_mapping` provides an instance-scoped
`RxnMapperProvider`; it does not import or depend on the legacy `chemtools`
adapter. Every generated mapping is reparsed, checked for exact preservation of
the original reactant/product structures, and passed through the typed mapped
edit normalizer. Provider and model version, model SHA-256, confidence,
coverage, warnings, and the mapped reaction are retained.

RXNMapper normally leaves atoms absent from the reported product unnumbered.
The adapter identifies only an unmapped reactant atom directly attached to a
mapped retained atom as a projected boundary. This recovers attachment losses
such as C=O and N–N cleavage in Fischer cyclization without mapping all bonds in
an unreported fragment or treating an unattached chloride as a broken bond.

`benchmarks/fischer_rxnmapper_poc.py` evaluated all 542 rows (356 unique
reaction SMILES). The compact report is
`results/fischer_rxnmapper_poc/report.json`; auditable per-reaction mappings
and normalized edits are in
`results/fischer_rxnmapper_poc/mapped_reactions.jsonl.gz`.

| Measurement | Result |
| --- | ---: |
| Valid, structure-preserving mappings | 356 / 356 |
| Full reported-product atom coverage | 356 / 356 |
| Exact agreement with existing signed edit sets | 58 / 58 |
| Ambiguous cases selecting exactly one internal hypothesis | 188 / 205 |
| Ambiguous cases matching no internal hypothesis | 17 / 205 |
| Mapper-only cases retained for review | 93 |
| Reactant-order stable edit profiles | 24 / 24 |
| Alternate-SMILES stable edit profiles | 24 / 24 |
| Median / mean mapper confidence | 0.661 / 0.628 |
| Batch mapping time, excluding existing analysis | 24.4 s |

The representative cyclohexanone/4-fluorophenylhydrazine reaction receives the
expected two formed aromatic bonds, C=O and N–N attachment losses, two bond
order changes, and three hydrogen losses. RXNMapper selects exactly one of the
two existing edit hypotheses at confidence 0.656.

This evidence is now available through an optional converter and query-time
integration, while retaining the POC's trust boundary:

- the provider runs only when supplied mapping and resolved internal evidence
  are absent;
- one external profile matching exactly one retained internal hypothesis
  receives `external_mapping_internal_consensus`;
- a valid external profile with no internal hypothesis receives
  `external_atom_mapping`;
- zero or multiple hypothesis matches preserve the original unresolved
  analysis and record an explicit conflict;
- every converted mapper-derived signature is `review_only` and therefore
  excluded from the default precedent index;
- the recommender may use a mapper-supported query signature to retrieve only
  already-verified indexed precedents, and adds mapper provenance and mandatory
  expert-review cautions.

Enable this mode explicitly:

```powershell
python -m condition_recommender.generic_conversion_cli `
  data-processor/reaction_dataset/Fischer_indole_synthesis.csv `
  results/fischer_indole_rxnmapper_conversion `
  --use-rxnmapper

python -m condition_recommender.generic_index_cli `
  results/fischer_indole_rxnmapper_conversion/records.jsonl `
  results/fischer_indole_rxnmapper_conversion/generic_index.json

python -m condition_recommender.generic_recommend_cli `
  "O=C1CCCCC1.Cl.NNc1ccc(F)cc1>>Fc1ccc2[nH]c3c(c2c1)CCCC3" `
  --records results/fischer_indole_rxnmapper_conversion/generic_index.json `
  --use-rxnmapper
```

The converter writes the assessment status, provider/model version and hash,
confidence, mapping coverage, mapped reaction, matched hypothesis IDs,
warnings, and error into `external_atom_mapping`. Review CSVs expose the
high-value fields plus the complete nested JSON. Mapper confidence is not an
admission or tie-breaking authority. The 188 selected hypotheses and 93
mapper-only POC results still require broader cross-chemistry or independent
chemist/provider validation before any future verified-precedent promotion;
the 17 disagreements remain explicit review cases.

The integrated converter was run over all 542 Fischer rows. Its report is
`results/fischer_indole_rxnmapper_conversion/conversion_report.json`:

| Integrated conversion measurement | Rows |
| --- | ---: |
| Mapper/internal single-hypothesis consensus | 324 |
| Mapper-only signatures | 111 |
| External/internal hypothesis conflicts | 38 |
| External signature unavailable | 10 |
| Existing internally resolved rows; mapper skipped | 59 |
| Rows with a reaction signature | 494 / 542 |
| Default-index eligible rows | 45 |
| Review-only / ineligible rows | 430 / 67 |

The 45-row default index is unchanged from the internal-only Fischer
precedent set; none of the mapper-participating rows entered it. The
representative cyclohexanone/4-fluorophenylhydrazine query then ran through the
persisted index with `external_mapping_internal_consensus` at mapper confidence
0.656. It selected `handle_signature`, found 11 compatible observations from
two independent support units, and returned two recipes with external-mapping
and expert-review cautions.

The featurizer, dataset-converter, and recommender desktop applications now
show a checked-by-default `Use RXNMapper` control. The converter runs one worker
while checked, records provider/model identity in the shard-reuse contract, and
exports mapping columns in the concise chemist review. The recommender caches
mapper-enabled and internal-only instances separately. The featurizer reports
the mapper disposition and confidence beside the structural analysis.

### 2.7 Knorr/Paal–Knorr transfer audit

The same name-free benchmark was run against the 652-row
`Knorr_pyrrole_synthesis.csv` source. The report is
`results/knorr_pyrrole_edit_poc.json`; the source reaction type selects the
cohort but is not passed to featurization or retrieval.

Among 519 unique reaction SMILES:

| Measurement | Count |
| --- | ---: |
| Verified name-free signatures | 310 (59.7%) |
| Global / fragmented correspondence signatures | 306 / 4 |
| Reference-disjoint queries with exact L3 support | 261 / 310 |
| Reference-disjoint queries with anonymous edit-graph support | 286 / 310 |
| Ambiguous queries with typed edit hypotheses | 5 / 5 |
| Ambiguous queries with a reference-disjoint all-hypothesis neighbor | 2 / 5 |
| Ambiguous queries with no consensus neighbor | 3 / 5 |

The hypothesis feature therefore helps Knorr/Paal–Knorr coverage, but its
increment is smaller than in Fischer because most usable Knorr records already
obtain a unique global-correspondence signature. The 24 signed reactions
without an independent edit-graph neighbor and the three non-consensus
ambiguous reactions remain abstentions in this cohort. These are neighbor
coverage measurements, not validation of condition suitability.

### 2.8 High-ROI canonical-site expansion

The current identity definitions add five conservative observation classes
without adding a named-family route or automatically registering a reaction
grammar:

| New observation | Canonical representation | Connectivity interface |
| --- | --- | --- |
| Explicit C⁻, N⁻, and O⁻ nucleophiles | `nucleophile_anion` | Connection endpoint requiring formal-charge neutralization |
| Imine, oxime, and hydrazone C=N | `PI\|PolarizedC=N` | Localized bond capacity and two connection endpoints |
| Epoxide and aziridine carbons | `EC\|StrainedRing\|...` | One C–O or C–N release link per possible opening carbon |
| Carbon-bound silyl ethers | `LG\|O\|SiR3` | Conditional O–Si release link |
| Bonded organolithium, organocopper, and H-free organoaluminium | `TM\|...\|Li/Cu/Al` | Carbon–metal transfer link |

The definitions deliberately do not infer anions from neutral molecules.
Nitro and azide resonance anions are not exposed as free nucleophiles, silanols
and siloxanes are not treated as silyl ethers, cyclopropanes are not strained
heterocycle centers, and DIBAL alkyl groups are not marked as carbon-transfer
partners.

`benchmarks/high_roi_site_coverage.py` scanned ten rows from each of the 118
bounded mixed datasets. The report
`results/high_roi_reactive_site_coverage.json` records 1,168 valid analyzed
rows and 120 rows containing at least one new observation:

| Observation | Rows | Datasets |
| --- | ---: | ---: |
| Explicit C/N/O anion | 28 | 4 |
| Polarized C=N | 53 | 14 |
| Silyl-ether O–Si link | 30 | 14 |
| Strained heterocycle | 10 | 6 |
| Li/Cu/Al transfer link | 4 | 2 |

Rows can contain more than one observation, so category counts are not
additive. This measures site incidence, not verified transformation or
recommendation coverage. New grammars and rewrites require their own chemistry
validation before they may consume these observations.

## 3. Implemented chemistry contracts

### 3.1 Observation

`reactive_taxonomy` currently provides:

- immutable atom-provenanced bond and schema-level hydrogen edits;
- typed `ReactionEvidenceCandidate` records for each attempted evidence
  provider, so failed interpretation does not erase the observation trail;
- deterministic `ReactionEditHypothesis` alternatives when correspondence
  supports several chemically distinct minimum-edit explanations;
- stronger internal before/after connectivity observations with explicit
  observed, projected, reconstructed, inferred, or unresolved scope;
- formal-charge and explicit stereochemical observations;
- canonical molecular reactive links, bond capacities, and connection
  endpoints;
- exact single-event and composable multi-event reconstruction;
- reaction topology and event relationships;
- typed product-atom completeness, partial product-origin gaps, deterministic
  `PTS1` partial-transformation keys, and `FSR1` fragment-source requirements;
- graph-derived local reactivity profiles and unchanged spectators;
- a serialized, template-free `ReactionCoreProjection` for mapped edit
  observations. The version 2 schema keeps every edit-participating atom as an atom transition,
  chooses a smaller set of primary centers only for explanation, removes each
  unchanged remote connected component once, and records every cut as a typed
  attachment port on that remote subgraph;
- four purpose-specific reaction-core identities: exact `RCX2`, typed `RCT2`,
  mapping-robust shape `RSH2`, and diagnostic center transition `RCS2`.

The remote classes (`aryl`, `heteroaryl`, `alkyl`, and related classes) are
derived from the removed molecular graph. They are not selected from a reaction
template or inferred from a source reaction name. Exact fragment SMILES and
functional-group evidence remain available beside the concise class.

`RSH2` is the only minimized identity used for retrieval. It includes generic
primary-center transitions, participant handle/site tokens, retained remote
shape, and event count. This prevents the broad center-only conflation seen in
V1, such as assigning the same retrieval identity to Suzuki and Heck reactions
that share a coarse carbon-center transition. `RCS2` remains explanation and
analysis output and is never a retrieval key.

Reaction-core projection is a sibling observation to RS3, not a substitute for
it. It does not change reaction-signature identity, named-family assignment, or
converter admission. Its schema and algorithm versions are serialized
explicitly, and stale core artifacts are rejected. The CLI and desktop
featurizer expose the minimized label, evidence, `RSH2`, diagnostic `RCS2`,
core size, remote subgraphs, attachment counts, and warnings; reaction batch
CSV exports the corresponding audit fields.

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
- reaction-core projection, evidence status, exact/typed/shape/center keys,
  atom transitions, events, remote subgraphs, attachment ports, and warnings;
- optional external-mapping disposition, mapped proposal, coverage, confidence,
  provider/model identity, matched internal hypotheses, and warnings;
- chemistry, condition, stage, outcome, and index-eligibility statuses;
- condition-source capability matches for each typed fragment requirement;
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
  -> exact reaction-core shape
  -> chemistry-gated anonymous edit-graph neighbors
  -> abstain
```

The current minimum adequate pool is two independent support units. Family use
requires high confidence and never bypasses edit compatibility. The
reaction-core tier uses exact `RSH2` equality and matching event count; it
never uses `RCS2` center equality. It is below verified bond-edit retrieval
because RS3 is more specific, and above anonymous edit-graph neighbors because
`RSH2` retains mapped center-state and participant evidence. The anonymous
edit-graph tier may cross an exact L3 key only after shared formed/broken bond
and ring-direction gates pass. Both fallback tiers disclose their use in
recommendation cautions.

Compatibility runs before similarity. Ranking definitions are separated from
retrieval definitions. Recipes are aggregated by `RCORE1`, while `RCR1`
variants, independent references, reference-local series, reaction breadth,
dataset breadth, usable yields, missingness, and uncertainty remain visible.

The fallback descriptor and retrieval policy are separate from RS3 retrieval.
They must return `recommendation_mode`, fallback evidence, stricter thresholds,
and cautions. A structure-neighbor fallback must never be described as a
verified reaction-center match.

Source-incomplete attachment replacements use an isolated exact-transformation
branch within that fallback path:

```text
exact PTS1 partial transformation
  -> hard FSR1 condition-source support
  -> compatibility and independent-support gates
  -> local-environment ranking
  -> review-qualified recommendation or abstention
```

The condition capability registry is generic over fragment graph, composition,
attachment element, resolved substance/family identity, and reported recipe
role. Its initial F/Cl/Br/I/CN/N3 entries are coverage seeds, not
reaction-family routing rules. A condition match supports reagent plausibility
but never upgrades the incomplete observation to verified atom mapping.

An unsigned query with two or more typed edit hypotheses takes a separate
query-only branch before ordinary structure fallback. The policy in
`condition_recommender/definitions/edit_hypothesis_retrieval.v1.json` requires:

- an anonymous edit prototype for every retained hypothesis;
- compatibility and minimum edit similarity for every hypothesis;
- intersection of the per-hypothesis verified-precedent pools;
- at least two independent support units;
- worst-case edit similarity for ranking; and
- explicit ambiguity and expert-review cautions.

This branch consumes only verified indexed precedents. It does not issue an
RS3 identity for the query and does not alter converter admission.

When optional external mapping produces a review-qualified query signature,
the normal verified-signature ladder is used against the same verified index.
The result reports `external_mapping_status`, provider, confidence,
`recommendation_mode`, and mandatory expert-review warnings. External
mapper-derived converted rows do not enter that index by default.

An unsigned query may also have a usable reaction core but no RS3 signature.
After edit-hypothesis consensus and before ordinary structure fallback, the
query-only core branch applies
`condition_recommender/definitions/reaction_core_retrieval.v2.json`:

```text
eligible mapped core with passing or review-qualified controls
  -> exact RCX2 identity when independently supported
  -> typed-local RCT2 identity when exact support is sparse
  -> mapping-robust RSH2 context when local support is sparse
  -> hard condition compatibility
  -> independent-support gate
  -> fallback-environment ranking
  -> review-qualified recommendation or abstention
```

This branch returns `recommendation_mode = reaction_core_review`,
`retrieval_strategy = reaction_core_ladder`, and mandatory expert-review
cautions. It never creates an RS3 signature for the query, upgrades its
admission, or admits external-mapper records as precedents. A core-only failure
continues to the existing structure fallback with both attempts preserved in
the retrieval trace.

The index stores exact, typed, shape, and center maps. Production retrieval
uses the first three as a narrow-to-broad ladder; the center-transition map is
audit-only because it is too coarse for condition transfer. Each core also
records map-number-independent equivalence, explicit H/charge/radical/isotope/
aromaticity/hybridization/stereo changes, graph-edit consistency, mapping
coverage, and pass/review/blocked quality. Presentation strings are generated
after identity construction and never participate in retrieval keys. Indexed
rows remain subject to the verified-record admission gate. Consequently,
external mapping may expand query coverage but does not manufacture trusted
training evidence.

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
- sparse mapped, internally verified precedent coverage for many `RSH2`
  reaction-core shapes and no production calibration of the core-only route;
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
- reaction-core shape hits, center-only hard negatives, core-only query hits,
  and core-to-structure-fallback abstentions;
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
- reaction-core shape coverage, independent support, mixed-label cluster audit,
  and false-neighbor review versus the center-only key;
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
