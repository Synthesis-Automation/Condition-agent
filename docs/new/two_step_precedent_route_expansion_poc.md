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
4. verifies observed segment IDs and unique lineage against the source route-core
   artifact;
5. verifies every curated substitute by exact identity against the supplier-stock
   snapshot;
6. enumerates the first-step input sets allowed at each level;
7. propagates every validated intermediate into the second step;
8. enumerates the second-step partners allowed at that level; and
9. retains complete operator, input, intermediate, product, validation, level,
   and source-reference provenance for every pathway.

### Observed-route panel

The primary evidence panel is versioned in
`core_retrosynthesis/definitions/two_step_observed_route_expansion_poc.v1.json`.
It contains three adjacent, uniquely linked, chemistry-resolved segments from the
5,000-route patent corpus:

- `US09238655B2_37209`: directed aryl addition followed by alcohol oxidation;
- `US06673800B2_537307`: nitrone cycloaddition followed by alcohol oxidation;
  and
- `US07612113B2_830855`: ring reduction followed by urea formation.

The runner verifies the declared patent, route tree, route core, step reactions,
source reaction IDs, and lineage link against
`routes.poc.core.v1.jsonl.gz`. The verified source SHA-256 is
`7396640f621c27666b8ac4dd065b3d51f0c6c7ddaa709224bdf646cd72ee7f51`.

All substituted inputs are exact matches in the local Mcule In Stock snapshot
dated 2026-08-02. The report retains supplier, collection, snapshot date,
availability, evidence level, region, and supplier record ID. A generated path
is rejected unless all declared structural components participate; this excludes
self-reaction and ignored-component artifacts.

### Synthetic regression panel

The original graph-smoke panel remains versioned in
`core_retrosynthesis/definitions/two_step_precedent_route_expansion_poc.v1.json`.
It contains N-, O-, and S-linked routes:

- alkylation followed by amine acylation;
- alkylation followed by alcohol acylation; and
- thiol formation followed by thioether alkylation.

These are intentionally small graph-chemistry regression probes. Their `POC:`
references make that evidence boundary explicit; they are no longer the primary
evidence panel.

## Result

All three observed segments replayed their exact intermediate and final product
at R0.

| Route | R0 products | R1 cumulative products | R2 cumulative products |
| --- | ---: | ---: | ---: |
| Aryl addition → oxidation | 1 | 7 | 16 |
| Nitrone cycloaddition → oxidation | 1 | 3 | 5 |
| Ring reduction → urea formation | 1 | 6 | 12 |
| **Total** | **3** | **16** | **33** |

The deterministic report is
`results/core_retrosynthesis/precedent_route_expansion/two_step_observed_poc.v1.json`.
Each pathway records both steps' operator IDs, directional operator IDs,
abstraction levels, reverse-round-trip results, edit-agreement results, and any
warnings.

Run it with:

```powershell
python -m core_retrosynthesis expand-precedent-routes `
  core_retrosynthesis/definitions/two_step_observed_route_expansion_poc.v1.json `
  results/core_retrosynthesis/precedent_route_expansion/two_step_observed_poc.v1.json `
  --stock-index cas_tools/data/stock_portfolio.sqlite `
  --route-core-source datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.core.v1.jsonl.gz
```

## Evidence boundary and next gate

This iteration establishes deterministic route-conditioned enumeration from
source-verified patent routes and stock-verified inputs. It does not establish
experimental viability, condition transfer, yield, selectivity, or useful
drug-like coverage. A patent-reported reaction is evidence for the exact R0
route, not experimental evidence for its R1/R2 analogues.

The next gate is chemist review and condition-aware attrition. Render a
stratified R1/R2 sample with both steps and source precedent side by side, attach
condition-compatibility cautions, and record accept/question/reject decisions.
After that review, scale to more route motifs selected without looking at their
expansion outcomes.
