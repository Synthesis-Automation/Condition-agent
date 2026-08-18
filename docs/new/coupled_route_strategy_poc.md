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

The relationship classes are:

| Relationship | Evidence | Admission |
| --- | --- | --- |
| `handle_progression` | Step 2 edits an atom in a bond formed by step 1 for every complete lineage candidate | strict |
| `same_site_coupled` | Both steps edit the same retained atom for every complete lineage candidate | strict |
| `shared_local_environment` | Active sites do not overlap but remain within two bonds | review |
| `lineage_ambiguous` | Coupling depends on unresolved or bounded correspondence | review |
| `unresolved` | One or both reaction cores lack admissible chemistry evidence | review |
| `independent_sites` | The carried active sites are absent or separated by more than two bonds | rejected |

Strict admission proves structural site continuity, not experimental viability
or strategic usefulness. Chemist review remains the next gate.

## Full 5,000-route result

Source:
`datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.core.v1.jsonl.gz`

| Measure | Result |
| --- | ---: |
| Source routes | 5,000 |
| Producer-consumer lineage pairs | 13,647 |
| Strict pairs | 6,589 |
| Review pairs | 2,933 |
| Rejected independent-site pairs | 4,125 |
| `handle_progression` pairs | 2,500 |
| `same_site_coupled` pairs | 4,089 |
| Shared-local-environment pairs | 1,589 |
| Lineage-ambiguous pairs | 878 |
| Unresolved pairs | 466 |
| Unique strict typed strategies | 4,520 |
| Strict typed strategies supported by at least two patents | 942 |

The requested aromatic sequence was recovered from route
`US04011323_913997`, source reactions `21308_0 → 21309_0`:

```text
Ar-H → Ar-NO2 → Ar-NH2
```

It is classified as `handle_progression`, has coupling score `1.0`, and its
typed strategy occurs in 21 distinct patents in this corpus.

## Review artifact

The deterministic 70-pair panel contains:

- 30 strict pairs;
- 20 local, ambiguous, or unresolved review pairs; and
- 20 rejected independent-site controls.

Each card shows the logical overall transformation, both physical reactions,
the carried intermediate, step-specific edit tokens, lineage overlap and graph
distance, strategy identities, and review controls. Decisions persist in local
browser storage and can be exported as JSON.

Artifacts:

```text
results/core_retrosynthesis/coupled_route_strategy/
  coupled_route_strategy_poc.v1.json
  coupled_route_strategy_poc.v1.html
```

Run the POC with:

```powershell
python -m core_retrosynthesis mine-coupled-route-strategies `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.core.v1.jsonl.gz `
  results/core_retrosynthesis/coupled_route_strategy/coupled_route_strategy_poc.v1.json `
  results/core_retrosynthesis/coupled_route_strategy/coupled_route_strategy_poc.v1.html `
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
