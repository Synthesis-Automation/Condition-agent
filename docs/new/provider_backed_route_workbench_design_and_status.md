# Provider-backed route workbench design and status

## Goal

The workbench is the first executable vertical slice from an agent action to a
deterministically audited route. It is intended for real examples and replay,
not as a claim that a returned route is synthesis-ready.

The LLM may choose a registered provider, a leaf ID, and a bounded proposal
budget. It cannot submit SMILES, edit a reaction, invent an operator, override a
compatibility result, or mark a route solved.

## Implemented flow

```text
target + explicit provider choice
              |
              v
identifier-only provider action
              |
              v
deterministic transition admission
              |
              v
bounded multistep beam search
              |
              v
precedents + route verification + typed issues
              |
              v
weakest issue + supported repair proposals / unresolved leaf
```

The implementation reuses the canonical components rather than introducing a
second route representation:

- `core_retrosynthesis.provider_multistep` adapts admitted provider output to
  `OneStepExpansionBatch` and projects route leaves to stable agent IDs.
- `core_retrosynthesis.route_workbench` runs the existing multistep planner,
  precedent lookup, route verifier, issue collector, and repair enumerator.
- `benchmarks.provider_route_workbench_examples` runs a CSV panel against an
  actual operator library and stock database and writes full JSON plus a compact
  case-summary CSV.
- `benchmarks.render_provider_route_workbench_review` turns that JSON into a
  self-contained chemist review with inline RDKit SVG structures and reactions.

Render a completed panel report with:

```powershell
python -m benchmarks.render_provider_route_workbench_review `
  results/provider_route_workbench_examples_20260830/report.json `
  results/provider_route_workbench_examples_20260830/review.html
```

Provider-local ranks are provenance, not cross-provider scores. A route step is
not strengthened merely because it was returned early by a provider.

## Result contract

For every retained route, the workbench records:

- the canonical route tree and terminal evidence;
- stable leaf IDs and whether each leaf is currently expandable;
- provider ID, provider-local rank, and admitted transition ID for every step;
- operator, strategy, strategic class, compatibility, and selectivity evidence;
- admitted precedents and local similarity values;
- route-verification gates;
- typed route issues and one deterministic weakest issue; and
- repair proposals when a deterministic repair tool supports that issue.

`solved` continues to mean that all leaves satisfy the configured terminal
predicate. It does not mean that a route has been experimentally demonstrated.
The molecular-weight terminal rule is a search heuristic and should not be
confused with catalog availability.

## Real-example panel, 2026-08-30

The initial panel contains umeclidinium cation, vilanterol, and a
thiol-containing stereochemical target. It was run against the validated
generic operator library and the local stock portfolio with depth 2, three
candidates per expansion, beam width 6, and four maximum expansions.

Observed result:

| Case | Retained result | Steps | Admitted transitions | Strongest issue |
|---|---:|---:|---:|---|
| Umeclidinium cation | partial | 2 | 9 | unresolved leaf |
| Vilanterol | partial | 2 | 12 | unresolved leaf |
| Thiol-containing target | partial | 2 | 12 | unresolved leaf |

All three runs completed without execution errors. Every retained best route
passed tree integrity, target identity, step graph consistency, and the current
forward-validation gate. None was labeled solved: each retained one unresolved
leaf and therefore failed terminal completion and overall route verification.
The workbench also preserved selectivity or tactical-step advisories.

This is a useful conservative result. It proves that provider admission,
multistep assembly, real stock lookup, precedent lookup, verification, and issue
selection can run together. It does not yet demonstrate automatic repair or a
chemically acceptable complete route.

## What the examples expose

The main remaining gap is no longer data plumbing. It is action coverage and
route-level chemical judgment:

1. `unresolved_leaf` is a valid next search action, not a protecting-group
   repair. The current repair enumerator correctly returns no fabricated repair.
2. Several returned steps are graph-valid and precedent-bearing but remain
   strategically questionable. Precedent similarity and operator validity do
   not by themselves establish route suitability.
3. The example runner does not invoke condition recommendation, so the verifier
   reports missing condition evidence.
4. The current forward gate checks evidence carried by admitted operators; a
   broader independent forward assessment remains a separate release gate.

## Next release gate

Add one deterministic loop controller over the current contracts:

1. choose the strongest unresolved issue;
2. expand an unresolved leaf, or select one supported repair action;
3. rerun transition admission and route verification;
4. accept the successor only when the targeted issue is removed without adding
   a stronger issue; and
5. stop on solved, no progress, repeated state, or budget exhaustion.

Test that controller first on curated single-step and two-step fixtures with
known expected actions. Only after those regressions should an LLM be permitted
to choose among the enumerated IDs.
