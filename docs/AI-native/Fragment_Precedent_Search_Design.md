# Fragment precedent search for scientific agents

Status: initial local implementation, 2026-09-29; independent chemistry and
agent-effectiveness evaluation remain pending. The [canonical-system roadmap](../new/type_agnostic_reaction_recommendation_implementation.md)
remains the authority for chemistry validation and release gates.

Implementation entry points are `reactive_taxonomy.fragment_search`,
`condition_recommender.fragment_index`, `condition_recommender.fragment_search`,
and the workspace's `search_fragment_precedents` operation. See the
[workspace instructions](readme.md#prepare-and-search-fragment-precedents) for
offline building and recorded calls.

The initial implementation makes these scoped choices:

- One immutable SQLite file contains the RDKit library, manifest, observation
  projections, and procedures. Existing baseline file hashing pins every component.
- Relationship evidence currently uses validated supplied maps, checked against
  the recorded structures. Inferred/reconstructed correspondence, conflicts, and
  observations without validated mapping remain unresolved; they still match
  product queries. The index reports this coverage explicitly.
- Procedure text is stored as pageable chunks with offsets and hashes, capped
  at 60,000 characters per procedure. Longer text is explicitly marked truncated;
  the source catalog remains identified in the index manifest. Automatic separate
  full-source snapshots are deferred.
- A prefix pilot is labeled `prefix_pilot`; a completed full input is labeled
  `full_source`. Neither label claims coverage of all published chemistry.
- The library checks deadlines between bounded batches. Workspace calls also use
  a killable subprocess covering imports and native matching, with persistent
  stage/error logs. A hard timeout can return no retained hits with unknown counts.

The following sections retain the design rationale and remaining evaluation plan.

## 1. Decision

Build one operation, `search_fragment_precedents`, that answers:

> Which recorded reactions contain this core in their products, and which
> provide evidence for constructing or modifying that core?

The agent chooses the fragment, interprets the evidence, and decides whether to
refine the query or pursue a route. The tool performs deterministic graph search,
projects recorded reaction changes onto the match, and returns a small, diverse
set of inspectable records. No proposed disconnection is required.

Start with substitution-tolerant substructure search over the existing canonical
corpus. Preserve the selected core's connectivity and ring topology by default.
Support explicit SMARTS constraints for deliberate broadening. Optional target-derived
candidate generation is now implemented as described below. Final core selection
remains with the agent. Unrestricted graph-edit similarity, new web-search
orchestration and an internal route planner remain deferred.

### Optional fragment suggestions

`reactive_taxonomy.search_fragments.suggest_search_fragments` provides a small
structural helper beside precedent search. It does not depend on the k-way
partition pipeline or validated retrosynthetic candidates. The workspace exposes
the same recorded operation; the regular research workbench exposes
`POST /api/v1/fragments/suggest` and **Suggest fragments** in the existing editor.

The versioned `search_fragments.v1` definition limits the target to 200 atoms,
materializes at most 64 proposals and returns at most five distinct queries.
Candidates include connected ring-bond systems (including fused, bridged and
spiro systems), one-bond contextual variants, rings with their connecting paths,
local functional/stereo regions, and the whole target when it fits query limits.
Candidates overlap and need not cover every target atom. Equivalent query strings
are deduplicated with one representative embedding; these are not all embeddings.

Extraction retains complete ring systems, multiple-bond context, charged valence
and specified stereochemistry. Each query is compiled through the precedent-search
compiler and verified against the exact selected target atom set. Extractions
that lose specified stereo after capping are rejected. Counts expose rejected
and truncated generation; no valid candidate is an acceptable result. The query
compiler now uses atomic-number hydrogen queries for explicit D/T atoms, and
reports `fragment_query_compiler.v2` in query identity/metadata. Stored product
graphs and index policy are unchanged, so existing indexes remain usable.

Results contain canonical target identity, source/input atom correspondence,
query-atom correspondence, omitted atoms, boundary bonds, structural descriptors,
reasons and cautions. An explicit `selected_atom_ids` request extracts exactly the
agent's chosen connected region or rejects it; it never silently expands it.
Boundary bonds do not assert synthetic disconnections or precursor feasibility.

Ordering first separates structurally simple queries, then uses declared candidate
kind and feature weights (ring junctions, heteroatoms, specified stereo and
nonaromatic multiple bonds), with diversity across proposal origins. These features
do not measure synthetic difficulty or corpus rarity. Common/simple motifs carry
a broad-search caution. Only the existing search can measure actual corpus breadth.

The agent may ignore suggestions, choose a different subgraph, or refine a query.
Generation invokes no corpus search, mapping, retro or forward prediction. The UI
shows each selected region on the target and requires a separate search action.
Development validation covers graph/stereo/provenance contracts and end-to-end
search compatibility. Whether suggestions improve route quality over agent-only
selection still requires a controlled agent/chemist evaluation.

The reproducible development panel is
[`examples/ai_native/fragment_suggestions_pilot.py`](../../examples/ai_native/fragment_suggestions_pilot.py).
It records 24 targets and optionally searches just the chosen cyclic-ether core
when an index is supplied. The initial run returned 68 graph-valid candidates;
the local [validation report](../../results/ai_native/fragment_suggestions_pilot/validation.md)
records limits, timings, test results and the remaining agent-evaluation gap.

## 2. What existing systems teach us

These are design references, not claims that we have their data coverage or
licensed access. Sources were checked on 2026-09-29.

| System | Relevant capability | Design consequence |
|---|---|---|
| [Open Reaction Database query implementation](https://github.com/open-reaction-database/ord-interface/blob/main/ord_interface/api/queries.py) | Input/output predicates support exact, substructure/SMARTS, and similarity searches using indexed RDKit chemistry. | Use a typed product-query contract over a prebuilt index. |
| [Pistachio query interface](https://nextmovesoftware.com/pistachio.html) | Accepts SMILES/SMARTS with component-role and search-type constraints, including product, substructure, and synthesis queries. | Keep structural predicates and product roles explicit even when an agent generates the query. |
| [RDKit SubstructLibrary](https://www.rdkit.org/docs/source/rdkit.Chem.rdSubstructLibrary.html) | Cached molecule storage, pattern-fingerprint screening, threaded graph matching, and result limits. | Reuse an existing local search engine; a new database service is unnecessary for the first version. |
| [Reaxys similarity search](https://www.elsevier.support/reaxys/answer/how-do-i-perform-a-similarity-search) | Product-only similarity searches compare molecules; reaction similarity considers reaction centers and their surroundings. | Keep product matching separate from evidence that a reaction constructed the core. |
| [CAS SciFinder similar reactions](https://cas-product-help.zendesk.com/hc/en-us/articles/19669803324301-Get-Similar-Reactions) | Narrow, medium, and broad searches vary the environment retained around a reaction center. | Make relaxation explicit. Center-based similarity becomes useful after a candidate transformation is identified. |
| [CASREACT user guide](https://cas-stnext.zendesk.com/hc/en-us/article_attachments/37936073871885) | Distinguishes groups present in a product from groups formed in a reaction. | Report presence, construction, modification, and unresolved evidence separately. |
| [SmallWorld](https://nextmovesoftware.com/smallworld.html) | Describes molecular searches using graph-edit distance and maximum common subgraphs. | Consider explainable analogue searching later; molecular distance alone does not establish a synthetic precedent. |

The useful addition for our agent is the combination of **core search, a local
reaction-change witness, source provenance, and bounded execution**. General web
search can then follow a returned publication or patent identifier. It does not
need to reproduce the structure index.

## 3. Corpus and ownership

The current [conversion report](../../datasets/literature/full/conversion_report.json)
contains 660,190 canonical observations, including
132,160 admitted verified records, 445,784 review records, and 82,246 rejected
records. These are observation counts, not distinct molecules or publications.
Only 21,518 procedure records are currently available. Counts are a local snapshot
and must be recomputed for each build.

Build the discovery index from `datasets/literature/full/combined_records.jsonl.gz`
or its canonical shard manifest. The condition-recommendation index excludes
records for reasons such as unresolved recipes that need not prevent product
substructure discovery. Index every valid product component, preserving source
admission status and exclusion reasons. A product with unresolved reactants or
mapping can be a structural lead with an unresolved reaction relationship.

Search never upgrades a record's chemistry verification or recommendation
admission. Retrieved conditions remain reported conditions; applying them to a
different target requires the existing compatibility and evidence checks.

| Owner | Responsibility |
|---|---|
| `reactive_taxonomy` | Query compilation, topology constraints, graph matching, atom-reference projection, and local reaction-change evidence. |
| `condition_recommender` | Derived corpus index, retrieval, deterministic grouping/ranking, and observation/reference/procedure joins. |
| `condition_registry` | Existing condition identities and resolved recipes; no new fragment logic. |
| `chem_coworker.scientific_workspace` | Thin operation adapter, baseline identity, deadlines, saved artifacts, and concise summaries. |

This is a new projection of the canonical corpus, not another converter. Promote
the existing streaming record reader from `conversion/concise_review.py` into
neutral corpus I/O and use it consistently. Index publication must validate all
expected shards and checksums; silently skipping an incomplete shard is invalid.

Existing `fragment_source_support` concerns the source of atoms or groups in a
reaction, and existing reaction-core indexes describe edits. Neither provides
arbitrary product-fragment search.

## 4. Small agent-facing contract

Proposed operation arguments:

```python
search_fragment_precedents(
    query: str,
    query_format: Literal["smiles", "smarts"] = "smiles",
    topology: Literal["preserve_rings", "subgraph"] = "preserve_rings",
    limit: int = 5,
    timeout_seconds: int = 10,
)
```

These are API defaults, not performance guarantees. Accept
1–10 returned observations and a deadline capped at 30 seconds. Resource limits
belong to a versioned policy rather than agent-controlled arbitrary scan sizes.

### Query meaning

- **SMILES:** retain atom identity, charge, aromaticity, bond connectivity, and
  specified isotope/stereochemistry. Ordinary implicit hydrogens allow peripheral
  substitution; explicitly constrained hydrogens remain constraints. Unspecified
  stereochemistry is unconstrained, not evidence of a stereochemical match.
- **SMARTS:** accept a documented bounded subset covering atom alternatives,
  bond alternatives, aromaticity, charge, hydrogen/degree, ring, and stereo
  constraints. Reject unsupported constructs with an actionable error. Compile
  through the centralized `compile_smarts()` cache.
- **`preserve_rings`:** preserve query atom/bond ring membership and the complete
  ring systems represented in the query. Reject additional fusion, bridging, or
  spiro rings incident to those ring atoms; permit peripheral branches and rings
  attached through acyclic bonds. Implement this as explicit graph constraints,
  not an assumption about RDKit's ordinary substructure match or SSSR numbering.
- **`subgraph`:** permit the fragment inside a larger ring system, while retaining
  explicit SMARTS restrictions. Echo this choice in every result.

For SMARTS, `preserve_rings` requires a concrete ring topology. A partial ring
query such as `[C;R]-[O;R]` or `c:c`, or alternatives allowing both ring and chain
atoms, must explicitly use `subgraph`; otherwise return a correction hint. Never infer
non-ring constraints that contradict an explicit ring predicate merely because
the SMARTS does not draw a closed cycle.

V1 searches one connected fragment in individual product components. It does not
silently remove salts from a disconnected query, neutralize charges, enumerate
tautomers, delete atoms, open rings, or relax stereo. Invalid/disconnected input
returns a correction hint rather than a guessed core.

Echo the original expression, compiled constraints, query atom references, query
policy version, and index identity. A query ID hashes its representation and
policy; it must not claim that arbitrary equivalent SMARTS expressions have been
semantically canonicalized.

### Progressive disclosure

The normal workspace `call_summary` returns the query interpretation, search
status, count scope, relationship groups, at most five compact hit cards, and an
artifact reference. Each hit contains:

- Stable hit, observation, reaction, and available reference identifiers.
- Matched product SMILES and a short explanation of local changes.
- Relationship evidence status, mapping ambiguity, and admission cautions.
- Citation availability and whether a procedure is linked to this observation.
- Literal `inspect_artifact` paths for the matched atoms, change witnesses,
  source record, reported conditions, and available procedure text.

Save bounded selected-hit details in the search artifact so the existing
`inspect_artifact` operation can reveal them without rerunning search. Store
procedure/source text as pageable chunks of at most 300 characters with exact
source offsets and a source hash. The current inspector truncates long scalar
strings and cannot page within them; a list of chunks fits its existing literal
path and list-pagination interface. Text over the initial implementation's
60,000-character per-procedure cap has an explicit truncation marker and source
hash. A later extension can save a separate full-text artifact without changing
the search operation. No second fragment-specific inspection operation is needed
initially.

Current `get_precedents` only hydrates condition-index records. Do not route these
broader discovery hits through it. Hydrate selected observations directly from
the derived index by `observation_id`. Build procedure lookups offline: exact
observation matches first; explicitly unassigned reaction-level procedures may
be shown with that scope. Never attach another observation's yield or conditions
merely because it has the same reaction ID.

### Search outcomes

| Outcome | Agent meaning |
|---|---|
| `complete` | All eligible candidates in the declared indexed scope were examined; zero hits means none in this corpus and query scope. |
| `too_broad` | Confirmed structural matches exceeded the breadth budget. Return count bounds and narrowing hints; no arbitrary globally ranked top list. |
| `partial` | A time, candidate, embedding, or evidence budget stopped work. Any returned hits are verified within the examined subset; ranking is subset-only. |

Keep invalid input, unavailable/stale index, and worker failures as distinct
typed execution errors. They are not successful zero-hit searches. Report
separate counts for distinct products, observations, and references, each with
`exact`, `at_least`, or `unknown` precision. A fingerprint shortlist is an upper
bound on possible hits, not a confirmed breadth count.

Relationship-group counts are distinct observations per group and are
non-additive: one observation may contain constructed and retained embeddings.
Global totals deduplicate IDs; embedding/witness counts are separate. Exact
relationship counts require completed classification, not just completed product
matching. If classification is interrupted, its counts retain lower-bound or
unknown precision even when the product count is exact.

## 5. Evidence about the matched core

Classify each query-to-product **embedding** using existing
`ReactionAtomReference`, `ReactionEdit`, `BondState`, correspondence hypotheses,
and completeness evidence. A match alone supports only product presence.

| Local evidence | Relationship |
|---|---|
| A verified newly formed bond is an edge of the matched query, with both atom origins accounted for and a verified prior `no_bond` state | `constructed`: some queried core connectivity was formed. |
| An internal bond order, atom state, or stereochemical state changed without new queried connectivity | `modified`. |
| An edit crosses from the matched core to atoms outside it | `boundary_changed`: attachment or functionalization. |
| Sufficient local correspondence establishes unchanged core and boundary states | `carried_through`. |
| Missing origins, insufficient correspondence, or conflicting admissible hypotheses prevent classification | `unresolved`. |

Construction means at least one core bond was built, not that the entire core was
made from simple starting materials in this step. A bond between matched atoms
that is not a queried edge is not a construction witness. Oxidation or
aromatization can be useful `modified` precedents without claiming connectivity
construction.

Return nonexclusive local-change flags and their atom/bond witnesses. A reaction
can construct a core and change its boundary in the same step. If a product
contains constructed and retained copies, retain both embeddings and report a
mixed relationship. Deduplicate equivalent witnesses without choosing a
convenient atom mapping. Mapping hypotheses that disagree remain visible.

Preserve the exact indexed atom ordering and its mapping back to each original
product component. Reusing old atom indices after canonicalizing a SMILES would
invalidate evidence. A minimized reaction-edit core is not complete product
correspondence: absence of an edit in that projection cannot establish retention.

An unrelated global completeness warning need not erase a verified local witness,
but the record's original warnings and admission remain visible. Reaching an
embedding cap prevents a definitive whole-record retention claim. Do not run a
new atom mapper or infer missing correspondence during search.

## 6. Retrieval, diversity, and execution cost

```text
Agent selects core
        |
Validate explicit query constraints
        |
Prebuilt product index: fingerprint screen -> exact graph confirmation
        |
Project saved edits onto each matched embedding
        |
Group diverse observations and attach source evidence
        |
Compact summary -> inspect saved detail -> agent decides next action
```

Use a derived immutable index with one row per unique indexed product component
and links to all original observations. Store canonical source identities,
original/indexed atom correspondence, saved edit evidence, and keyed source and
procedure records. Deduplication of molecules must not merge experimental
observations or condition recipes.

Use RDKit `SubstructLibrary` with a pattern-fingerprint screen and read-only
SQLite metadata. The initial layout stores the serialized library inside the same
SQLite file, so existing workspace file hashing covers the whole artifact. The
internal manifest binds source identities and versions. Benchmark cold loading
before considering a separate fingerprint store or persistent worker. No automatic
runtime rebuild.

The manifest records source/shard hashes, source and indexed counts, exclusions,
RDKit version, schema/definition versions, parsing/aromaticity/stereo policies,
fingerprint settings, relationship rules, procedure catalog identity, and a
successful-build marker. Missing, incompatible, or incomplete indexes fail with
an offline build instruction.

Only use a fingerprint screen whose necessary-feature behavior is validated for
the supported query syntax. Confirm every result by graph matching. Unsupported
screening must not create false negatives; reject that query feature in V1 or
use a separately bounded exact path. Morgan similarity is not an admission gate.

Rank verified construction evidence first, then useful modifications/boundary
changes, then retained/unresolved leads with their limitations. Group by observed
local transformation and reference so five hits can expose different approaches.
Use stable IDs for ties, preserve group counts even if omitted from the preview,
and allow artifact inspection of saved alternatives. Repeated examples from one
document are not independent support. Unknown publication families must remain
unknown rather than being deduplicated by a guessed identity.

Budgets cover loading, candidate checks, embeddings, linked observations,
relationship interpretation, ranking, hydration, and serialization. `limit=5`
alone does not bound the work. Use a cancellable worker with an end-to-end
deadline; record stage timings and the stopping reason in existing progress logs.
Bound thread count to avoid exhausting the workspace host.

Common biphenyl is not blacklisted. Stop after a configured number of confirmed
distinct-product matches and return `too_broad`, with suggestions such as retaining
the heteroatom pattern, substitution position, complete ring system, or a specific
functional handle. Do not invent counts for the unsearched corpus. If a compute
limit is reached without confirmed broadness, return `partial` instead.

Measure cold subprocess startup as well as warm search: current agent scripts
can start fresh Python processes. Initial performance targets are p95 under
5 seconds cold and 2 seconds warm for selective development queries on a stated
machine, with a 10-second default deadline. These are acceptance targets to test,
not current measurements; broad-query refusal and process cleanup are separate
latency checks. Add a persistent worker only if measured cold-start cost warrants
the extra lifecycle complexity.

## 7. Example: the cyclic ether target

For the recent target:

```text
CC(c1cc(COc2c3cccc2)c3cc1)=O
```

The parsed canonical target is `CC(=O)c1ccc2c(c1)COc1ccccc1-2`. An initial core
query is `c1ccc2c(c1)COc1ccccc1-2`: two benzene rings and the central six-membered
oxygen-containing ring, retaining the CH2–O sequence and ring connections.

Search that complete core with `preserve_rings`, allowing peripheral substitution.
The acetyl substituent remains part of the agent's target context; omitting it
from the retrieval query does not establish its compatibility with a precedent's
conditions. Inspect hits that actually construct a queried C–O or C–C edge, as
well as clearly labeled core modifications. Do not assume either category is
present in the current corpus before querying it.

If there are no hits, the agent can submit an explicitly revised SMARTS or use
returned structural/source clues for literature search. Each call records its
own constraints; no hidden broadening runs inside the first call. If the exact
core still lacks coverage, return that fact and let the agent continue with
single-step retrosynthesis. The previous proposed route is a development case,
not a validated reference answer.

## 8. Implementation plan and acceptance gates

| Phase | Deliverable | Exit gate |
|---|---|---|
| 1. Freeze contracts and audit coverage | Typed query/result/evidence models, versioned matching policy, small fixtures, and a stratified audit of product validity, full correspondence, and procedure joins. | SMILES/SMARTS and ring-policy semantics agreed; missing evidence represented honestly; cold-load prototype measured before storage layout is frozen. |
| 2. Build offline index and graph retrieval | Strict canonical reader, deduplicated product library, observation/source tables, complete hashed manifest, atomic publication, and bounded substructure search. | Screened results equal exhaustive graph search on fixtures; source counts reconcile; incomplete or stale artifacts cannot be published/used. |
| 3. Add local-change evidence and useful ranking | Per-embedding witnesses, ambiguity handling, source diversity, selected-hit hydration, and exact procedure joins. | Construction, modification, boundary, retention, mixed, and unresolved cases pass; no cross-observation evidence leakage. |
| 4. Expose one workspace operation | Capability reporting, baseline pinning, compact summaries, artifact paths, cancellation, progress timings, and optional task advice. | A fresh workspace process searches and inspects evidence within budgets; unavailable/broad/partial results remain distinguishable. |
| 5. Evaluate with agents and chemists | Equal-budget comparison of current web/single-step tools versus those tools plus fragment search. | Demonstrated retrieval usefulness, honest evidence, and acceptable cold latency before default exposure; add broader analogue search only for measured gaps. |

Keep the guide advisory: an agent may use web search, core precedent search,
single-step disconnection, or its own investigation order. Calling this tool must
not trigger forward prediction, condition recommendation, recursive searching,
or automatic route expansion.

Required deterministic regression cases include:

- Cyclization forming a queried edge; coupling inside a biaryl query versus
  coupling outside a phenyl-only query; oxidation and unchanged-core controls.
- Multiple core copies, symmetry, conflicting mappings, missing core atom
  origins, and unrelated remote completeness warnings.
- Extra fusion/spiro/bridging, changed ring size or oxygen position, substituted
  cores, aromatic/Kekule forms, specified stereo, and explicit SMARTS alternatives.
- Common-fragment breadth, fingerprint false positives, embedding/observation/time
  limits, complete zero hits, invalid inputs, and stale/missing index artifacts.
- Product atom-order changes, duplicate observations/publications, unavailable
  procedures, and wrong-observation procedure joins.

For the agent evaluation, include named drugs, unnamed substituted cores, rare
ring systems, and deliberately broad common fragments. Have chemists label
whether retrieved steps construct or usefully modify the requested core. Measure
precision at five, recall on curated known precedents, source/procedure fidelity,
useful transformation diversity, runtime, calls, and token use. Assess whether
the evidence improves proposed route ideas separately from whether a route is
experimentally feasible. Freeze a held-out panel before tuning.

Run owning-package tests and the full `pytest -q` suite for implementation
changes. New definition files require schema/loader validation. Development
pilots do not replace the independent chemistry review or untouched evaluations
required by the primary roadmap.

## 9. Explicitly deferred

- Automatic core extraction/ranking: first measure whether agent-selected cores
  cause retrieval failures.
- General analogue retrieval by graph-edit distance, MCS, tautomer expansion, or
  embeddings: first establish exact-query coverage and explainable constraints.
- Licensed external corpus connectors: add only through available authorized
  interfaces, preserving their provenance and access conditions.
- A separate search agent, mandatory workflow, automatic learning of chemistry
  facts from search success, and rebuilding the existing recommendation system.

The initial success criterion is practical: an agent can submit an unfamiliar
core and quickly receive a few source-linked synthetic leads, understand what
each match proves, and decide the next investigation step.
