# Coupled Route-Strategy Mining POC

**Status date:** 2026-08-18

## Decision

Do not treat arbitrary adjacent reactions as one synthetic strategy. A strict
two-event strategy must have cross-step structural evidence that the second
reaction consumes or further transforms the site created or edited by the first
reaction.

The POC is an interpretation and retrieval artifact. It does not collapse two
physical reactions into one executable reaction, and it does not replace either
constituent reaction signature or route observation.

## Structural classifier

`core_retrosynthesis.coupled_route_strategy` evaluates each producer-consumer
lineage link in a `RouteCoreProjection`:

1. collect mapped active atoms from both generic reaction signatures;
2. restrict the first step to active atoms retained in the carried intermediate;
3. translate those atoms into the second reaction's map frame through every
   stored lineage candidate;
4. calculate active-atom overlap and minimum graph distance in the intermediate;
5. preserve incomplete or relationship-dependent lineage as review evidence;
6. assign content-derived exact, typed, shape, occurrence, and overall-reaction
   identities.

Version 2 preserves the version-1 structural relationship as comparison
evidence and adds a causal `dependency_class`:

| Relationship | Evidence | Admission |
| --- | --- | --- |
| `handle_progression` | Step 2 edits an atom in a bond formed by step 1 for every complete lineage candidate | strict |
| `same_site_coupled` | Both steps edit the same retained atom for every complete lineage candidate | strict |
| `shared_local_environment` | Active sites do not overlap but remain within two bonds | review |
| `lineage_ambiguous` | Coupling depends on unresolved or bounded correspondence | review |
| `unresolved` | One or both reaction cores lack admissible chemistry evidence | review |
| `independent_sites` | The carried active sites are absent or separated by more than two bonds | rejected |

The refined dependency classes are:

| Dependency | Evidence | Admission |
| --- | --- | --- |
| `created_handle_consumed` | Step 2 transforms an atom in a bond installed by step 1 | strict |
| `activation_then_conversion` | Step 2 breaks a step-1-installed bond and forms or changes another heavy-atom bond at that handle | strict |
| `continued_site_transformation` | Both steps transform the same lineage-traced site without consuming a newly installed bond | strict |
| `temporary_group_removed` | Step 2 removes a step-1-installed bond without a mapped heavy-atom replacement | review |
| `shared_local_environment` | Non-overlapping active sites are within two bonds | review |
| `lineage_ambiguous` | The dependency type changes across correspondence candidates, or lineage is incomplete | review |
| `unresolved` | One or both step cores lack admissible chemistry evidence | review |
| `independent_sites` | No carried site continuity is present | rejected |

Strict admission proves structural causal continuity, not experimental
viability or strategic usefulness. Chemist review remains the next gate.

## Full 5,000-route result

Source:
`datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.core.v1.jsonl.gz`

| Measure | Result |
| --- | ---: |
| Source routes | 5,000 |
| Producer-consumer lineage pairs | 13,647 |
| Strict pairs, version 1 | 6,589 |
| Strict pairs, refined version 2 | 6,495 |
| Review pairs, version 1 | 2,933 |
| Review pairs, refined version 2 | 3,027 |
| Rejected independent-site pairs | 4,125 |
| `handle_progression` pairs | 2,500 |
| `same_site_coupled` pairs | 4,089 |
| `created_handle_consumed` pairs | 1,535 |
| `activation_then_conversion` pairs | 871 |
| `temporary_group_removed` pairs | 91 |
| Refined dependency-ambiguous pairs | 881 |
| Shared-local-environment pairs | 1,589 |
| Unresolved pairs | 466 |
| Unique strict typed strategies, version 2 | 4,451 |
| Version-2 strict strategies supported by at least two patents | 931 |

Of the 2,500 version-1 `handle_progression` pairs, version 2 resolves 1,535
as created-handle progression, 871 as activation followed by conversion, 91 as
temporary-group removal requiring review, and 3 as lineage-dependent.

The requested aromatic sequence was recovered from route
`US04011323_913997`, source reactions `21308_0 → 21309_0`:

```text
Ar-H → Ar-NO2 → Ar-NH2
```

It retains the baseline `handle_progression` relationship and is refined to
`created_handle_consumed`, has coupling score `1.0`, and its typed strategy
occurs in 21 distinct patents in this corpus.

## Review artifact

The new deterministic 70-pair panel contains:

- 30 strict pairs;
- 20 local, ambiguous, or unresolved review pairs; and
- 20 rejected independent-site controls.

The HTML starts with the complete version-1-to-version-2 count comparison. Each
card then shows the baseline relationship beside the refined dependency, the
logical overall transformation, both physical reactions, the carried
intermediate, step-specific edit tokens, transient/replacement bond evidence,
strategy identities, and review controls. Version-2 review decisions use a
separate browser-storage key and can be exported as JSON.

Artifacts:

```text
results/core_retrosynthesis/coupled_route_strategy/
  coupled_route_strategy_poc.v1.json
  coupled_route_strategy_poc.v1.html
  coupled_route_strategy_poc.v2.json
  coupled_route_strategy_poc.v2.html
```

Run the POC with:

```powershell
python -m core_retrosynthesis mine-coupled-route-strategies `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.core.v1.jsonl.gz `
  results/core_retrosynthesis/coupled_route_strategy/coupled_route_strategy_poc.v2.json `
  results/core_retrosynthesis/coupled_route_strategy/coupled_route_strategy_poc.v2.html `
  --strict-sample-size 30 `
  --review-sample-size 20 `
  --rejected-sample-size 20 `
  --include-reaction-pair 21308_0,21309_0
```

## Evaluation gate

The first decision metric is chemist-reviewed precision among strict pairs.
Report precision separately for `handle_progression` and `same_site_coupled`,
and use the rejected controls to measure false-negative rate. Do not promote a
typed strategy to planning evidence solely because it is frequent.

A later executable composite action must additionally require:

1. patent-disjoint recurrent support;
2. acceptable chemist-review precision;
3. both physical operators independently applicable and forward-validatable;
4. explicit intermediate and two-step condition compatibility; and
5. improved held-out strategy recovery or search efficiency without reducing
   chemically valid fallback coverage.

Protection/activation sequences at different sites remain outside this strict
POC. They may be strategically related, but route adjacency alone cannot prove
that causal relationship.
