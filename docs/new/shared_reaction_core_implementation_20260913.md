# Shared reaction core: experimental implementation

Date: 2026-09-13

Status: an executable development implementation of the shared graph contract,
derived dataset artifact, candidate union and application integration. Independent
chemist adjudication, untouched evaluation, Full/Compact migration and removal
of the superseded default retrieval paths remain release gates.

Design reference: [shared reaction core proposal](shared_reaction_core_design_proposal_20260913.md).
Release authority: [primary roadmap](type_agnostic_reaction_recommendation_implementation.md).

## Implemented contracts

`reactive_taxonomy.shared_reaction_core` owns immutable projection and comparison
contracts. It parses original structures, validates stored atom references and
observed bonds against those graphs, and retains the existing observation core
ID, evidence quality and warnings. Conflicting or ambiguous evidence cannot
produce a qualified projection. No hypothetical reactants are enumerated.

The new namespace is `shared_reaction_core.v1`, separate from the existing
condition-signature and retro-template L-number namespaces:

| View | Implementation |
| --- | --- |
| Whole reaction | Normalized reactants and products agree; atom-map numbering is excluded from identity |
| L0 / `observed_local` | Existing detailed observed core plus the product-side local graph and protected multiple-bond units |
| L1 / `retained_local` | Qualified transformation and before/after center states, with one-shell product environment; departing/source port identities generalized |
| L2 / `retained_typed` | Same protected join and functional units, retaining reacting-center classes while generalizing surrounding environment |

Atom/edge labels are canonicalized as a colored incidence graph. The label
vocabulary is part of the payload, avoiding graph-local color collisions.
Correspondence provenance remains separate from canonical labels. Definition,
algorithm, upstream core/schema and RDKit versions bind the derived keys.

The first generalization operator covers one intermolecular single-bond join
with one graph-qualified departing or hydrogen port on each joining endpoint.
It supports the tested C-C, C-N, C-O and C-S examples through one implementation.
Nitrile and other connected nonaromatic multiple-bond units remain protected.
Versioned definitions qualify departing roots; they do not assign a universal
leaving-group reactivity score or establish condition transfer.

Hydrogen loss at supported N/O/S and sp-carbon sources is recorded as an
explicit realization difference. Literal `C#N` remains HCN input and is never
relabelled as an unspecified cyanide source. Aromatic C-H edits remain protected.

## Retrieval and aggregation

`condition_recommender.shared_core_index` builds one atomic SQLite projection
artifact over an already admitted generic index. It stores both direct and
product-side lookup keys, unavailable-level reasons and per-observation hashes.
It verifies source-index identity and checks candidate row and payload hashes.
Existing molecular observations are reused; no Full re-featurization is needed
to build this derived artifact when the stored evidence is sufficient.

`recommend_indexed_signature()` selects the new path only when an explicit
shared-core artifact is provided. The ordinary `GenericConditionRecommender`
remains the application entry point. The experimental path does not invoke the
old Br/I hypothetical-query path or fall back silently to another comparison
contract when its query is unavailable.

Candidate channels:

1. Whole-reaction and L0/L1/L2 direct lookup.
2. Product identity and protected product-core lookup when support is sparse or
   broad exploration is requested.
3. Precedent IDs supplied by an existing retro operator through the normal
   `preferred_reaction_ids` bridge.

All channels use the same original-query graph comparison and original-query
recipe/constraint compatibility checks. A shared product with a different
transformation is excluded with `DIFFERENT_TRANSFORMATION_SAME_PRODUCT` in the
audit trace. A missing operator-to-condition link is reported explicitly.
Product-side retrieval does not run the retrosynthesis planner or infer missing
experimental conditions. It uses the observed product projections in the
condition artifact; operator seeds are an additional input to that same path.

Queries have fixed literal inputs in this initial contract. Results are
`same_setup` or `analogue_evidence`; automatic input substitution and generated
alternative setups are not implemented. Each analogue retains its actual
precedent reaction and source-port requirements. Recipes with different required
input identities are aggregated separately, even when their recipe core IDs agree.

Independent support controls broadening. `top_k` controls the display limit.
Candidate union is capped at 512 with an explicit truncation warning. One
observation found by several channels is counted once. Existing reference-aware
support and recipe scoring are reused after chemistry qualification.

The recommendation response schema is **4.2**. Additive fields include
`shared_core_trace`, `match_namespace`, `evidence_relation`, `candidate_channels`
and `source_input_requirements`. Existing default result fields remain intact.
The frontend shows the new L labels without adding a conflicting ordinal label.

## Run the prepared development workbench

A 68-record development artifact is available locally under
`results/shared_core_validation/`. Its rows were selected from the current
571,157-row trusted Full index: 64 deterministic sample positions plus the
previously investigated cyanation references. It is not a representative
accuracy benchmark or the Full library.

```powershell
$env:CONDITION_SHARED_CORE_EXPERIMENTAL = '1'
python -m app.web_api --workbench --index results/shared_core_validation/generic_index.sqlite
```

The workbench's Full selector points to that explicitly supplied development
index in this command. It does not mean all Full observations are loaded.
Remove the environment flag to return to the current default retrieval path:

```powershell
Remove-Item Env:CONDITION_SHARED_CORE_EXPERIMENTAL
```

## Build from another canonical development index

The canonical conversion workflow can build the artifact in the same run:

```python
from condition_recommender.conversion import build_recommendation_artifacts

report = build_recommendation_artifacts(
    "path/to/source.csv", "results/development_build",
    build_fast_index=True, experimental_shared_core=True,
)
```

The workflow report includes the projection artifact checksum and manifest.
Cancellation checks and progress callbacks cover projection building, and an
interrupted build preserves the previous projection artifact. Existing index
rebuilds without this option cannot silently activate an older projection set:
binding validation rejects incompatible artifacts.

An already converted canonical index can be backfilled independently:

```powershell
python -m condition_recommender.shared_core_cli build `
  path/to/generic_index.sqlite path/to/generic_index.shared_core.sqlite

python -m condition_recommender.shared_core_cli query `
  path/to/generic_index.sqlite path/to/generic_index.shared_core.sqlite `
  'C#N.IC1=NNC2=C1C=CC=C2>>N#CC1=NNC2=CC=CC=C12' --scope broad
```

The web runtime expects a sibling named `<index-stem>.shared_core.sqlite` when
the experimental flag is set. Missing or incompatible artifacts fail explicitly.
Different Full, Compact or review-admission indices need matching artifacts.
The Python API also accepts `GenericConditionRecommender.from_path(index_path,
shared_core_path=projection_path)` without environment configuration.

## Development evidence and review

The pre-change complete suite passed **1,391 tests in 354.77 seconds**.
New deterministic regressions cover leaving-group/source variation, C-C/N/O/S
examples, graph and map-order invariance, negatives and conflicts, source-bound
artifacts, atomic build failure, support, budgets, compatibility, retro seed
parity, HTTP integration and training-only review artifact construction.

The 68-record backfill produced 67 eligible observed-local projections and 24
eligible L1/L2 projections. The remaining records retain explicit missing-evidence
or unsupported-generalization reasons; no source admission tiers were changed.
Full and Compact artifacts were not modified.

For the user's literal-HCN indazole query, the development artifact retrieves:

- `US05801183:417039_0` at L1;
- `US07166621B2:767236_0` at L2;
- `US08450363B2:1326394_0` at L2.

These are related cyanation observations, not identical reactions or validated
recipes for literal HCN. The before/after audit is saved in
`results/shared_core_validation/development_audit.json`; per-query responses
retain complete candidate and compatibility explanations. Timings in that
development report use different cache states and are not a performance claim.

The existing blind-review generator now accepts `--experimental-shared-core`:

```powershell
python -m condition_recommender.chemist_review_cli `
  results/shared_core_validation/generic_index.sqlite `
  results/shared_core_validation/chemist_review `
  --experimental-shared-core --max-cases 12 --minimum-pool-size 2
```

The generated development packet contains 12 cases, 19 candidate entries and
12 negative controls. Training contains 54 rows, held-out development contains
14, and the reported shared reference/reaction group overlap is zero. Product
and direct projections are built from training rows only. The manifest binds
the source and projection artifacts by hash. This packet is ready for review;
it has not been independently adjudicated and is not the untouched release set.

## Remaining gates and removal criteria

The broad operator intentionally abstains on multi-event changes, most charge or
stereochemical changes, ring transformations, unqualified departing fragments,
unsupported source representations and unresolved correspondence. A detailed
observed-local view may remain available. General source-class query semantics,
qualified alternative-input proposals and general multievent L1/L2 alignment
need further contracts and chemistry validation.

Before default cutover: complete independent review, resolve disagreements,
freeze broader acceptance thresholds and evaluation partitions, run the untouched
evaluation, and measure Full-scale lookup/coverage/performance. Then backfill
the Full/Compact artifacts and remove the superseded hypothetical-query and
duplicate retrieval paths. Keeping the experimental selection is temporary,
with removal tied to those parity and release gates rather than indefinite
backward compatibility.

## Final automated validation

- Complete `pytest -q`: **1,433 passed in 359.80 seconds**, including 42 added
  regressions beyond the pre-change baseline. The final log is
  `results/shared_core_validation/pytest_final.log`.
- Frontend `npm run build`: passed in 21.64 seconds. The existing large Ketcher
  bundle warning remains.
- Workbench Playwright check: passed; verifies literal-query preservation,
  three cyanation precedents, L1/L2 display semantics, product-side evidence,
  strict-scope abstention and no page JavaScript errors.
- HTTP activation, canonical conversion with projections, cancelled/failed
  atomic builds, artifact binding and training-only review generation are
  covered by the complete suite.
- New modules pass Ruff's undefined/unused-name checks; `git diff --check`
  passes. Source hashes are recorded in
  `results/shared_core_validation/implementation_manifest.json`.

The temporary browser-test server was stopped. No Full or Compact dataset was
rebuilt, and no independent adjudication or untouched release result is claimed.
