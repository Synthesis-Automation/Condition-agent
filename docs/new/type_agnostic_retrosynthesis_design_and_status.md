# Type-Agnostic Retrosynthesis: Design and Implementation Status

**Status date:** 2026-08-14
**Primary implementation:** [`core_retrosynthesis`](../../core_retrosynthesis/)  
**Current product priority:** reliable, diverse single-step disconnections

## Purpose and authority

This document describes the retrosynthesis system that exists in the repository
today: its intended design, implemented behavior, evidence boundaries, measured
results, and remaining work. It is the current design-and-status reference for
retrosynthesis. Package READMEs remain useful as command references and experiment
logs, but some of their historical terminology and results predate the generic
operator system described here.

`core_retrosynthesis` is now the canonical package owner. This is not a claim
that the system is production-ready: broader held-out evaluation and integration
of the strategy-level result contract into public interfaces are still required.

The condition-recommendation architecture is documented separately in
[`type_agnostic_reaction_recommendation_implementation.md`](type_agnostic_reaction_recommendation_implementation.md).

## Executive status

The strongest implemented path is a deterministic, data-derived single-step
pipeline:

1. Compile graph operators from atom-mapped reactions.
2. Admit an operator only when it reconstructs its source precursors.
3. Retrieve operators using target-observable features and product SMARTS.
4. Apply each operator in reverse with atom correspondence retained.
5. Forward-analyze every proposed precursor set.
6. Keep only proposals whose forward signature agrees with the operator.
7. Rank the verified proposals and group alternate realizations of the same
   synthetic strategy.

The system is type-agnostic in the intended sense: named reaction families may be
reported as annotations, but they are not required routing keys and they do not
define operator identity.

| Area | State | Practical meaning |
| --- | --- | --- |
| Generic graph-operator compilation | Implemented | Learns observed edit operators without a seven-family gate. |
| Source round-trip admission | Implemented | Every stored template must reproduce its source precursor set. |
| L2/L1/L0 retrieval ladder | Implemented | Falls back from local context to more general operators explicitly. |
| Reverse application and forward validation | Implemented | Returned candidates must reproduce the expected reaction signature. |
| Deterministic ranking and diversity | Implemented | Evidence-aware ranking with explicit policy versions and score-band guards. |
| Strategy identity and grouping | Implemented in Python API | `STRAT1` groups verified alternate precursor realizations after validation. |
| Full-scale sharded build | Implemented | Restartable builds, ledgers, manifests, exact deduplication, and merge. |
| Condition-aware reranking | Optional/partial | Can reorder close structural candidates; cannot create or rescue them. |
| Selectivity competition model | Review-only POC | Emits cautions; it is not calibrated for production decisions. |
| Shared route-tree contract | Implemented | Observed and planned routes use the same evidence-neutral, occurrence-preserving schema. |
| Observed route corpus | 5,000-route POC | Patent-disjoint curated routes are converted with zero rejected records; strategy annotations remain weak labels. |
| Route-core projection | Implemented | 18,644/18,647 steps have generic cores, signatures and minimized drawings; cross-step symmetry remains explicit. |
| Observed route-action labels | 5,000-route POC | 18,280/18,647 steps have retained-edit labels, 18,078 have synthon labels, 18,307 have exact precursors, and 13,099 have STRAT1; the independent executable-operator facet is 1,201/18,647. |
| Multistep planning | Bounded POC | Available for experiments, but not the current product priority. |
| Strategy-aware validation scheduling | Not implemented | Grouping currently improves output quality, not validation cost. |
| Public CLI/API for grouped strategies | Not implemented | Current grouped result contract is Python-only. |
| Broad release benchmark and chemist review | Review artifact ready | Full labels are audited and 120 promoted examples are rendered; human decisions and leakage-controlled replay remain. |

## Design principles

### Molecular evidence is authoritative

The source of truth is the molecular graph: parsed structures, atom
correspondence, bond edits, reactive atoms, local environments, and reconstruction.
Reaction names and source labels can provide supporting evidence, but cannot
override contradictory structure-derived evidence.

Evidence is used in this order:

1. valid supplied atom mapping and observed bond edits;
2. successful graph-operator reconstruction;
3. compatible interpretation annotations;
4. unresolved or conflicting evidence retained as warnings or review records.

The system does not invent atom correspondence or silently turn an ambiguous
transformation into a named reaction.

### Type-agnostic does not mean chemistry-agnostic

Candidate generation is generic, but chemistry checks remain hard gates. A
proposal must be structurally applicable and must survive forward analysis. Named
families are optional annotations; graph transformation identity is mandatory.

### Generation, validation, and preference are separate

- **Generation** asks which graph operators can disconnect the target.
- **Validation** asks whether the proposed precursors reconstruct the intended
  transformation.
- **Preference** ranks chemically valid alternatives using evidence, context,
  support, simplicity, compatibility, and optional condition evidence.

This separation prevents a ranking heuristic or learned preference from being
mistaken for chemical validity.

### Deterministic and versioned behavior

Public identities and serialized artifacts are based on normalized chemistry and
versioned definitions, not display labels or reaction names. Ranking,
compatibility, strategic-complexity, and competition policies are stored in
validated definitions where practical.

## System boundary

The current code has the following dependency direction:

```text
reactive_taxonomy          condition_registry
        |                         |
        +------> core_retrosynthesis <------+
                         |
                condition_recommender
                 (optional adapter)
```

`reactive_taxonomy` owns reaction observation and structural chemistry.
`condition_registry` owns condition-substance identity and contextual roles.
`condition_recommender` can add condition evidence after a structural candidate
exists. `core_retrosynthesis` composes these capabilities and owns the generic
operator, search, ranking, strategy, and route-planning contracts.

The package organization is:

```text
core_retrosynthesis/
  chemistry.py, mapping.py, row_io.py   shared deterministic foundations
  generic_*                            generic operator library and search
  strategy_*                           strategy identity and grouped results
  ranking / compatibility / multistep  downstream preference and planning
  baselines/cx_rdchiral/                evaluation-only historical baseline
```

The former top-level `retrosynthesis_poc` package and its standalone CLI have
been removed. Its shared chemistry and mapping utilities now belong to
`core_retrosynthesis`; the still-useful RDChiral C-X implementation is isolated
under `baselines` so comparison code does not look like an alternative product
path. The family-constrained compiler remains only for controlled regression
comparisons. The data-driven generic compiler is the canonical direction.

Historical result-directory names and serialized v2/v3 definition IDs are left
unchanged so existing operator libraries and benchmark artifacts remain readable;
they are data compatibility identifiers, not importable package boundaries.

No new-system package should depend on legacy `chemtools`.

## Single-step data flow

```text
mapped source reactions
        |
        v
mapping + reaction analysis + graph edits
        |
        v
generic operator compilation at L2 / L1 / L0
        |
        v
source round-trip admission
        |
        v
versioned operator library + retrieval index

target molecule
        |
        v
necessary-feature retrieval -> product SMARTS match
        |
        v
reverse operator application with atom correspondence
        |
        v
forward reaction analysis and signature agreement
        |
        v
compatibility + context + strategic-complexity evidence
        |
        v
ranked candidates -> STRAT1 groups -> alternate realizations
```

### Compilation and admission

The generic compiler recomputes mapping and reaction analysis from the source
reaction. It requires:

- a valid analysis and materialized reaction core;
- verified product completeness;
- one product;
- a compilable operator at the requested abstraction level; and
- exact source-precursor recovery when the operator is applied backward to its
  own source product.

Two admission modes exist:

- `supported` retains the historical seven-archetype constraint for controlled
  comparisons;
- `data_driven` admits any structurally verified operator and is the intended
  generic direction.

Rejections are recorded with explicit reasons rather than being silently dropped.

### Abstraction ladder

Operators are compiled at three specificity levels:

- **L2:** context-rich local environment;
- **L1:** reaction core and directly relevant handle context;
- **L0:** the most generalized fallback.

Search proceeds from L2 to L1 to L0. The abstraction level used is retained in
the result so fallback is visible rather than hidden inside a score.

### Retrieval and validation

The retrieval index uses necessary product-observable features only to prune the
library. It is not a chemistry validator. Product SMARTS matching and graph
application remain authoritative.

For each reverse-applied candidate, the system runs forward reaction analysis and
rejects proposals with unresolved analysis, missing identity, or an operator
signature inconsistent with the compiled operator. A retained candidate therefore
has `forward_validation_status="verified_signature"`.

## Identity hierarchy

Identity is deliberately split so chemistry, target site, strategy, and purchasable
realization are not conflated.

| Namespace | Meaning | Main normalized inputs |
| --- | --- | --- |
| `OP1` | Generic graph transformation | Bond-edit/operator signature |
| `SITE1` | Target-specific disconnection site | Canonical target and edited product bonds |
| `SYN1` | Synthon-level precursor concept | Operator plus normalized precursor skeletons with precursor-only handles removed |
| `STRAT1` | One synthetic strategy for this target | `OP1 + SITE1 + SYN1` |
| `REAL2` | Concrete precursor realization | Operator, handle signature, and precursor SMARTS |
| `GRT3` | Stored generic template | Versioned template chemistry and metadata |
| `COMP2` | Handle-completion evidence group | Operator/synthon completion observations |

All digest inputs are canonicalized. Identities do not use source reaction names or
display labels.

`STRAT1` is intentionally target-specific. It groups candidates that make the same
edit at the same target site and produce the same synthons, while retaining
different leaving groups, protecting groups, or other concrete precursor choices
as alternate realizations.

## Candidate ranking and strategy grouping

The flat ranker is deterministic and controlled by
[`retrosynthesis_ranking.v3.json`](../../core_retrosynthesis/definitions/retrosynthesis_ranking.v3.json).
It combines validation evidence, operator support, local context, compatibility,
strategic complexity, and precursor realism. Diversity is enforced using operator,
site, and synthon identities. More generalized abstraction levels do not silently
jump ahead of more specific levels.

The hierarchical ranker, defined by
[`hierarchical_ranking.v2.json`](../../core_retrosynthesis/definitions/hierarchical_ranking.v2.json),
orders verified results at three levels:

1. target disconnection site;
2. synthon strategy;
3. concrete realization.

Completion priors use explicit support backoff from operator-plus-synthon, to
operator, to global evidence. They may reorder candidates only within guarded
abstraction and score bands.

The current `STRAT1` API adds a simpler immutable strategy view on top of the
existing validated search:

- only `verified_signature` candidates can enter a strategy;
- all alternates must have the same target and `STRAT1` identity;
- identical precursor sets are deduplicated;
- the first existing ranked candidate remains the representative;
- support is conservatively aggregated as independent maximum support rather than
  naively summing correlated template observations.

This grouping is post-validation. It improves the diversity and readability of
the returned result but does not yet avoid validating many equivalent realizations.

## Compatibility, conditions, and selectivity

Hard graph and forward-signature validation occur before softer evidence.
Precursor compatibility can demote or warn, but cannot create a candidate or
rescue an invalid transformation. The current policy is
[`precursor_compatibility_policy.v1.json`](../../core_retrosynthesis/definitions/precursor_compatibility_policy.v1.json).

Condition recommendation is optional and downstream of structural validation.
Only verified candidates are submitted. Condition evidence can distinguish direct,
fallback, and insufficient support, and may reorder only candidates in the same
abstraction level and a narrow score band. Missing conditions do not invalidate a
structurally valid disconnection.

The endpoint-competition/selectivity model is a review-only POC. It currently
covers a narrow class of single-event reactions and relies on weak labels. It can
emit warnings, but it is deliberately fail-open and cannot affect admission or
ranking. It should not be presented as a calibrated yield or major-product model.

## Artifacts and scale

`GenericTemplateLibrary` currently writes schema `3.0`; the reader supports older
schemas needed for experiment recovery. A full-scale build is sharded and
restartable. Each shard writes:

- a partial operator library;
- a compressed admission ledger;
- a source/config manifest.

Merge uses exact disk-backed deduplication and produces a versioned operator
library plus build report.

The local compact full-scale build at
[`full_scale_v3/compact/build_report.json`](../../results/operator_retrosynthesis_poc/full_scale_v3/compact/build_report.json)
contains:

| Metric | Value |
| --- | ---: |
| Source rows | 113,923 |
| Accepted observations | 76,225 (66.9%) |
| Operators | 3,691 |
| Completion groups | 4,926 |
| Realizations | 9,420 |
| Templates | 31,798 |

The largest recorded rejection classes are mapping unavailable (20,343), reaction
core not verified (13,268), and source round-trip failure (2,378). These counts are
useful coverage diagnostics for this local artifact; they are not a production
release statement.

## Interfaces

The CLI currently exposes library construction, generic and operator search,
ensemble comparisons, coverage audits, bounded route planning, and report
rendering. See [`cli.py`](../../core_retrosynthesis/cli.py) for the exact
arguments.

The strategy-grouped single-step interface is currently Python-only:

```python
from core_retrosynthesis import disconnect_strategies

strategies = disconnect_strategies(
    target_smiles,
    library,
    top_k_strategies=10,
)
for strategy in strategies:
    print(strategy.strategy_id, strategy.representative.precursor_smiles)
    print(strategy.alternate_realizations)
```

It is not yet wired into the CLI, API, or application UI. That omission is one of
the highest-value integration tasks because grouped strategies are the intended
user-facing unit.

## Measured evidence

The following repository artifacts are useful checkpoints. They were produced by
different historical experiment configurations and should not be combined into a
single accuracy estimate.

| Experiment | Targets | Main result |
| --- | ---: | --- |
| Balanced held-out comparison | 86 | Baseline-first L1 ensemble exact precursor recovery: 37.2% @1, 80.2% @5, 96.5% @10. |
| Diverse round-robin | 27 | Core-context path: 100% coverage, 74.1% exact @1, 77.8% exact @10. |
| C-C/Michael stress set | 60 | Core-context path: 56.7% exact @1, 88.3% exact @10; site recovery 91.7% @1. |
| Operator specificity audit | 24 | Data-driven L2/L1 reached 100% coverage and 87.5% site @25; L0 added no gain on this small sample. |

Artifacts:

- [`balanced_500_hydrogen_normalized/comparison.json`](../../results/core_retrosynthesis_comparison/balanced_500_hydrogen_normalized/comparison.json)
- [`balanced_round_robin_100/comparison.json`](../../results/diverse_retrosynthesis_poc/balanced_round_robin_100/comparison.json)
- [`stress_cc_michael_250/comparison.json`](../../results/diverse_retrosynthesis_poc/stress_cc_michael_250/comparison.json)
- [`all_shards_400/specificity_comparison.json`](../../results/operator_retrosynthesis_poc/all_shards_400/specificity_comparison.json)

Exact source-precursor recovery is a strict and useful regression metric, but it is
not a complete measure of synthetic usefulness: a different precursor realization
or synthon strategy may be equally valid. The new `STRAT1` contract exists partly
to enable evaluation at the correct level.

## Known limitations

- Strategy grouping is not yet benchmarked on a broad untouched corpus using
  `STRAT1@k`, site coverage, alternate-realization recall, and invalid-rate metrics.
- Candidate validation is still realization-first, so equivalent leaving-group or
  protection variants can consume redundant compute before grouping.
- Ranking scores are deterministic priorities, not calibrated probabilities of
  reaction success or yield.
- There is no learned target-to-operator relevance model; retrieval and ranking
  use explicit structural features and empirical support.
- Stereo, multi-event, and unusual mapped reactions are limited by mapping quality,
  reaction-core materialization, operator compilation, and source round-trip
  success.
- Condition and selectivity evidence is incomplete and intentionally subordinate
  to structural chemistry.
- The historical RDChiral baseline and family-constrained compiler remain in an
  explicitly evaluation-only boundary until generic-path parity tests make them
  removable.
- The shared route-tree contract represents one concrete synthesis tree, not the
  alternative choices in an AND/OR search network. Search state should remain a
  separate contract rather than overloading observed route data.
- Converted route steps do not yet carry `ReactionSignature` or canonical
  condition-recipe references. The schema reserves these links, but their owning
  packages must populate them from molecular and condition evidence.
- Multistep planning uses bounded search and simple terminal/stock evidence. It
  should remain secondary until single-step strategy quality and interfaces are
  stable.

## Recommended next work for single-step retrosynthesis

The clean, high-ROI sequence is:

1. **Add strategy-level evaluation.** Extend held-out evaluation with `STRAT1@k`,
   `SITE1@k`, alternate-realization recall, candidates validated per returned
   strategy, and runtime. Preserve exact-precursor metrics for regression.
2. **Schedule validation by provisional strategy.** Interleave candidates across
   distinct operator/site/synthon groups before spending work on a second
   realization of the same strategy. Keep final `STRAT1` assignment dependent on
   forward verification. This should improve both latency and top-k diversity
   without a new model.
3. **Expose the grouped contract.** Add `disconnect-strategies` to the CLI and use
   the same typed result in API/application integration. Do not create a second
   serialization model in the app layer.
4. **Improve realization ordering inside a verified strategy.** Use existing
   completion support, compatibility, stock/realism, and optional condition
   evidence under strict score-band guards. Keep the representative and alternates
   explainable.
5. **Run a release-quality audit.** Freeze source/library versions, evaluate an
   untouched and chemistry-stratified set, report rejection and fallback reasons,
   and conduct blinded chemist review of useful non-reference strategies.
6. **Retire evaluation-only paths.** Once the benchmark and interface are stable,
   remove the isolated RDChiral/family-constrained paths after parity tests pass;
   keep `core_retrosynthesis` as the single package owner.

This sequence keeps the system simple: it improves the unit returned to users and
the allocation of existing validation work before adding learned ranking or a
larger multistep architecture.

## Verification expectations

Changes to this system should preserve deterministic IDs, partner-order
invariance, source round-trip admission, forward-signature validation, and the
L2/L1/L0 fallback contract. Tests should include positive, negative, ambiguous,
and conflicting evidence, with representative C-C, C-N, C-O, C-S, Suzuki, mapped
unknown-family, invalid-map, and deterministic-identity cases.

The full repository suite remains the handoff gate:

```powershell
pytest -q
```
