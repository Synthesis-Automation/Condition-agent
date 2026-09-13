# Shared reaction core v2: general edit-graph beta

Date: 2026-09-13

## Status and scope

V2 extends the shared-core implementation from a qualified single-join subset
to general observed before/after edit graphs. It is an integrated, opt-in beta;
it is not a default replacement or a claim that recommended conditions transfer
successfully between substrates.

The previous [rollout report](shared_reaction_core_rollout_20260913.md) identified
the main architectural gap: reactions outside single intermolecular joins often
had no broad core. V2 addresses that gap in `reactive_taxonomy`, using the same
derived-index, compatibility, recipe-ranking and application pipeline.

## Chemistry contract

One canonical colored graph contains observed center atom correspondence,
before/after bond states, explicit edit nodes, local molecular context and
functional-unit closure over double/triple bonds. Context atoms remain attached
to their observed side; the implementation never invents their correspondence.

Center labels retain element, charge, aromaticity, hybridization, isotope,
radical state, hydrogen count, absolute stereochemistry and ring membership.
Before/after center distances retain intramolecular geometry and ring-closure
size when local shells no longer visibly connect the centers. Edit topology
therefore remains meaningful for multiple events and multiple products.

| Level | Preserved | Relaxed |
| --- | --- | --- |
| Whole reaction | Canonical supplied reactants and products | Serialization and component order |
| L0 `observed_local` | Actual edits, center states, functional units and radius-one before/after context | Remote structure beyond the core |
| L1 `retained_local` | Protected edits and radius-one context | Independently qualified departure/source ports |
| L2 `retained_typed` | Protected edit graph, center states, functional units and center topology | Additional surrounding context, using radius-zero functional-unit closure |

For transformations without qualified source ports, L1 may preserve the same
chemistry as L0. L2 can still relax remote context. Reductions, oxidations,
cleavages and cyclizations no longer need to resemble substitution reactions
to receive broad projections.

Departure/source abstraction remains an optional graph-validated operator for
single intermolecular joins. Its previous protections remain: aromatic C–H
substitution is not silently equated with Ar–halide substitution, broken bonds
are not discarded indiscriminately, and literal HCN remains a different input
from CuCN. Multievent departure/source abstraction is not asserted merely
because the general multievent graph can now be represented.

No named reaction family is required. Conflicting atom references, bond states,
hydrogen direction or event counts reject the projection. Invalid or ambiguous
correspondence does not gain permission to broaden. Appearing centers without
reactant correspondence remain unavailable. Pure stereo/state-only changes
that the upstream featurizer does not yet admit as signatures also remain
outside this contract.

## Retrieval and recipe evidence

Direct whole/L0/L1/L2 keys, product-side keys and explicit retro precedent links
still feed one comparator against the original query, followed by condition
compatibility and canonical recipe ranking. A matching product reached through
different edits remains a rejected candidate. No recursive retro planner call
or hypothetical reactant rewriting was added.

Recipe aggregation now uses an anchored observed source-port graph where that
context is qualified. A recipe supported on several remote scaffolds can occupy
one result slot, retaining each observation's actual input requirements and
citations. Different leaving/source attachments remain separate. Unqualified
source contexts retain exact input grouping. Evidence from different abstraction
levels is also kept separate.

Top-k remains an output limit, not a requirement to fill every slot. Automatic
search stops when narrower matches provide the requested independent support.
Broad search may add more distant graph-qualified analogues. For example,
Compact's ketone-reduction check returns three narrow recommendations in automatic
mode and five in broad mode. Candidate-budget warnings remain explicit.

## Versioning and dataset migration

The projection schema and namespace are `2.0` / `shared_reaction_core.v2`.
Retrieval definitions are `shared_core_retrieval.v2@2.0`. The anchored source-port
identity is part of the versioned chemistry definitions. V1 definition files
were replaced; stale v1 artifacts fail binding checks instead of silently loading.

Canonical observations, source admission and registry recipes are unchanged.
Only the derived `generic_index.shared_core.sqlite` artifacts need rebuilding.
The existing bounded parallel/checkpoint builder publishes completed artifacts
atomically and verifies their source identities. Current conversion builds that
select shared cores automatically produce v2 projections through this same path.

## Completed corpus builds

| Library | Canonical observations | V1 broad projections | V2 L0/L1/L2 projections | V2 artifact size |
| --- | ---: | ---: | ---: | ---: |
| Full | 571,157 | 178,347 | 561,283 | 3,972,395,008 bytes |
| Compact | 94,643 | 28,171 | 92,597 | 666,800,128 bytes |

The remaining observations have explicit missing-evidence or contradiction
reasons. No admission tiers or source labels changed. The definition hash is
`SCD1:e168a158aaf8fcdd2a34186528e6aefd62421d9ba38267f86c02075714904e1e`.
Full completed in 1,216.42 seconds with eight workers; Compact in 805.97 seconds
with two workers. Builds overlapped with tests/evaluation, so these timings
are not controlled performance measurements.

Both artifacts pass SQLite integrity, foreign-key and source-binding checks,
plus row/payload verification of 128 sampled projections each. Full and Compact
query audits cover cyanation, aldehyde/ketone reduction, alcohol oxidation and
amide dehydration. Full returns five recommendations for the user's bromide
and iodide cyanation queries; Compact returns four. Frequent transformations
can reach the explicit 512-candidate budget, so reported candidate counts are
bounded search results rather than exhaustive corpus counts.

## Fresh evaluation

A new 2,000-record panel excludes all previously used development/holdout
observations, publications, canonical reactions and mapping-equivalent records.
Selection preserves random draw order while filtering; sorting occurs only after
selection, correcting the earlier cohort's positional sampling bias.

The source IDs, selected positions, split membership, code/definition hashes,
seed and acceptance gates were frozen before retrieval outcomes were inspected.
Shared-core indices contain training rows only. The grouped split uses 1,602
training and 398 test records; the scaffold-disjoint split uses 1,726 and 274.

| Split | Existing coverage | V2 coverage | Existing top-5 recipe recovery | V2 top-5 recipe recovery |
| --- | ---: | ---: | ---: | ---: |
| Publication/reaction grouped | 79.65% | 67.84% | 7.79% | 6.28% |
| Scaffold disjoint | 77.01% | 67.52% | 11.31% | 8.76% |

Both paths report zero hard incompatibilities under the existing automated
rules and complete explanations for covered queries. Publication/reaction
leakage is zero; the scaffold split additionally has zero scaffold overlap.
These are sparse-training engineering evaluations, not Full-training performance
estimates or independent proofs of chemical precision. The two splits reuse the
same panel and are not independent trials. Historical recipe recovery is not
experimental success prediction.

The frozen cutover limits were at most five percentage points of coverage loss
and two points of top-five recovery loss, plus independent review of the new
chemistry cases. Coverage misses the limit in both splits; scaffold recovery
also misses its limit. Default cutover therefore remains disabled. No chemistry
rules were tuned against these outcomes.

A diagnostic replay inspected the top baseline precedent for the 53 grouped
queries and 29 scaffold-split queries covered only by the existing path. In
52 + 29 comparisons the candidate's protected graph differed; one grouped query
lacked a usable observation. None of these examined top precedents qualified
under v2 but was missed by lookup. This is not an independent judgment that all
excluded analogues are chemically unsuitable; it identifies comparison-policy
disagreements for review. Details are in `disagreement_audit.json`.

The old, already consumed diagnostic cohort showed coverage rising from v1's
30.85% to an intermediate v2 build's 63.68%. That is development evidence only,
not a substitute for the fresh results above.

## Review and use

The expanded blind packet is
`results/shared_core_v2/chemist_review/review_packet.html`. It is generated from
the older development cohort, not the fresh panel. Decisions are left blank.
It contains 30 cases, 63 candidate entries and 30 negative controls, using
797 training records and 203 held-out development records with zero group leakage.
The raw results, frozen manifest and reproducible scripts are under
`results/shared_core_v2/`.

To use the completed Full/Compact v2 artifacts, stop the previous server in its
terminal, then launch without a development-sample `--index` argument:

```powershell
$env:CONDITION_SHARED_CORE_EXPERIMENTAL = '1'
python -m app.web_api --workbench
```

The workbench displays the actual selected dataset size. The same opt-in backend
also works with the condition-only UI. Unsetting the environment flag and
restarting selects the existing default retrieval path.

The next release work is independent adjudication of broader matches and the
remaining disagreements, followed by a new frozen evaluation if chemistry is
changed. The current panel is consumed evidence now. Broader coverage alone is
not a reason to weaken the graph or compatibility gates.

## Final automated validation

- Complete `pytest -q`: **1,454 passed in 709.73 seconds**.
- Frontend build: passed in 59.18 seconds; existing bundle-size warning.
- Development browser regression: **1 passed in 12.8 seconds**, covering original
  query preservation, product-side evidence and strict-scope abstention.
- Full browser regression: **1 passed in 14.0 seconds**, covering actual library
  counts, five cyanation recommendations and non-join aldehyde reduction.
- Both browser checks reported zero page JavaScript errors; screenshots were
  inspected. The expanded review packet also renders in headless Edge.
- Modified Python modules pass Ruff name checks; `git diff --check` passes.

Final source hashes, artifact hashes and validation results are recorded in
`results/shared_core_v2/implementation_manifest.json`. The temporary port-8026
servers were stopped; restart the user's existing server to load v2.
