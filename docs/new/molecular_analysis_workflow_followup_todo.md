# Molecular analysis follow-up to-do list

Created: 2026-10-11

Status: reusable engineering changes implemented after completion of the walkthrough
and the user's implementation request. See
[strategy-led planning implementation](strategy_led_planning_implementation_20261011.md).
Independent review and multi-target workflow evaluation remain open.

## Context

Target input:

```text
Cn(cc(c1ncc(CC)cc1)c2c3CCC2)c3=O
```

Canonical representation from the recorded analysis:

```text
CCc1ccc(-c2cn(C)c(=O)c3c2CCC3)nc1
```

The automatic motif summary did not explicitly identify the pyridone/lactam
environment. The existing amide definition in
`reactive_taxonomy/definitions/molecular_motifs.v1.json` uses
`[C:1](=[O:2])[N:3]`, which requires non-aromatic carbon and nitrogen. In this
target, RDKit represents the corresponding ring atoms as aromatic `c` and `n`,
and that pattern returns no matches. The molecular graph still exposes the
ring nitrogen adjacent to the carbonyl carbon.

Local investigation evidence, when available:
[step 1 report](../../results/ai_native/stepwise_pyridone_20261011/step1.html)
and [graph inventory](../../results/ai_native/stepwise_pyridone_20261011/structure.json).
Generated investigation artifacts are local and may not exist in another checkout.

## Intended analysis workflow

Use the existing general analysis tool as the reusable foundation, focused
inspection for specific planning questions, and recorded custom RDKit code when
the existing interfaces cannot answer a question. Promote recurring, validated
custom checks into the owning package.

## Implementation checklist

- [x] Review the completed planning walkthrough and consolidate additional
  analysis gaps before choosing implementation scope and priorities.
- [x] Extend the existing molecular motif definitions to recognize the target's
  pyridone environment. Assess whether a distinct subtype or annotation is
  appropriate; do not blindly broaden the ordinary amide rule or inherit all its
  reactivity tags.
- [x] Keep structural recognition separate from reactivity interpretation.
  Preserve the observed carbonyl–nitrogen relationship even when interpretation
  is incomplete or ambiguous.
- [x] Add regression coverage for this target, representative pyridones and
  ordinary amides/lactams, and negative or overlapping heteroaromatic cases.
  Test equivalent representations without silently treating distinct tautomers
  as identical.
- [x] Audit atom-reference consistency across molecular summaries, site
  inspection, drawings and focused disconnection. Make canonical versus
  input/component atom scopes explicit and retain correspondence where needed.
- [x] Make analysis coverage limitations explicit: no registered motif match
  must not imply that the functionality is absent. Do not claim exhaustive
  functional-group coverage.
- [x] Review the general analysis output for planning relevance: identity,
  charge, stereochemistry, ring systems, motifs, stable atom references and
  limitations. Keep the default summary compact and provide details on demand.
- [x] Check whether existing focused inspection APIs answer bond, ring,
  neighbourhood and target–precedent comparison questions. Extend those APIs
  only for demonstrated gaps; avoid a duplicate molecular-analysis system.
- [x] Define when an agent should use custom RDKit analysis: state the unresolved
  question, reuse canonical structures and existing helpers, and record code,
  inputs, versions, outputs and evidence references in the scientific workspace.
  Keep custom observations distinct from validated domain rules and experimental
  feasibility claims.
- [x] Assess existing caching before adding any new cache. Reuse unchanged
  analyses with keys that include molecular identity, settings, schema and
  definition versions, and relevant chemistry-engine version. Preserve atom
  reference semantics when reusing results.
- [x] Establish a promotion process for recurring custom checks: observed gap,
  focused investigation, reusable definition or implementation, positive and
  negative tests, then integration into `reactive_taxonomy`.

## Implementation and validation boundaries

- [x] Put molecular chemistry in `reactive_taxonomy`; keep application and
  workspace layers thin. Do not reintroduce legacy dependencies or parallel
  public chemistry paths.
- [x] Use versioned, validated definitions and the centralized SMARTS cache.
  Review effects on downstream interpretation, compatibility and retrieval when
  motif definitions or tags change.
- [x] Validate definition loaders and schemas, run targeted chemistry
  regressions, and run the complete `pytest -q` suite before handing off an
  implementation change. Explain chemistry-driven snapshot changes.
- [ ] Evaluate whether the changes reduce repeated custom work and improve
  planning decisions without introducing false motif assignments. Do not infer
  route feasibility from improved molecular annotation.

The original target is covered by deterministic regression tests. Comparative
route-discovery performance and independent chemistry review remain unmeasured.
No new result cache was added; existing workspace evidence and SMARTS caches are
reused. The detailed report documents schema, indexing and migration boundaries.

Validation completed: `pytest -q` — **2980 passed, 4 skipped** (803.54 seconds).
No chemistry snapshots were regenerated. Recorded target smoke checks passed.
