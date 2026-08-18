# Two-Step Precedent-Route Expansion POC

**Status date:** 2026-08-17

## Question

Given one established two-step route, what nested chemical space can be built
without losing the structural connection to that precedent?

This is intentionally different from the optional route-action policy. The
action policy asks whether route precedents can rank retrosynthetic search
choices. This POC executes the observed transformations forward and measures the
products reachable from a bounded, declared building-block panel.

## Expansion levels

The levels are cumulative. A product visible at a narrower level remains visible
at every broader level.

| Level | Inputs | Operators | Claim |
| --- | --- | --- | --- |
| R0 | Observed inputs only | Route-derived L2 | Exact observed-route replay |
| R1 | Observed plus close declared building blocks | Route-derived L2 | Expansion under the observed local reaction context |
| R2 | All declared POC building blocks | Route-derived L2 and L1 | Locally relaxed, explicitly extrapolative expansion |

R2 is not a yield or feasibility guarantee. It states that both graph
transformations apply, the products sanitize, the observed edit signature is
preserved, and the reverse operator recovers the proposed inputs.

## Implementation

`core_retrosynthesis.precedent_route_expansion`:

1. validates the declared linear two-step route;
2. compiles each observed reaction into L2 and L1 generic operators;
3. admits only source-forward-round-tripped operators;
4. enumerates the first-step input sets allowed at each level;
5. propagates every validated intermediate into the second step;
6. enumerates the second-step partners allowed at that level; and
7. retains complete operator, input, intermediate, product, validation, level,
   and source-reference provenance for every pathway.

The bounded panel is versioned in
`core_retrosynthesis/definitions/two_step_precedent_route_expansion_poc.v1.json`.
It contains N-, O-, and S-linked routes:

- alkylation followed by amine acylation;
- alkylation followed by alcohol acylation; and
- thiol formation followed by thioether alkylation.

These are intentionally small graph-chemistry probes, not literature-backed
scope claims. The `POC:` references make that evidence boundary explicit.

## Result

All three routes replayed their exact intermediate and final target at R0.

| Route | R0 products | R1 cumulative products | R2 cumulative products |
| --- | ---: | ---: | ---: |
| Amine acylation | 1 | 9 | 20 |
| Alcohol acylation | 1 | 9 | 20 |
| Thiol alkylation | 1 | 6 | 15 |
| **Total** | **3** | **24** | **55** |

The deterministic report is
`results/core_retrosynthesis/precedent_route_expansion/two_step_poc.v1.json`.
Each pathway records both steps' operator IDs, directional operator IDs,
abstraction levels, reverse-round-trip results, edit-agreement results, and any
warnings.

Run it with:

```powershell
python -m core_retrosynthesis expand-precedent-routes `
  core_retrosynthesis/definitions/two_step_precedent_route_expansion_poc.v1.json `
  results/core_retrosynthesis/precedent_route_expansion/two_step_poc.v1.json
```

## Evidence boundary and next gate

This POC establishes deterministic route-conditioned enumeration. It does not
yet establish experimental viability, condition transfer, yield, selectivity,
building-block availability, or useful drug-like coverage.

The next evidence gate should replace the synthetic panel with several clean,
literature-backed two-step routes and a curated purchasable building-block set.
For each level, report step attrition, incompatibility reasons, condition
warnings, novelty, scaffold diversity, and chemist accept/question/reject review.
R2 should remain visibly separate from precedent-local R1 results.

