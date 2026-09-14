# Advisory planning alongside condition recommendations

The Workbench now offers **Reaction context: possible products and alternative
precursors**, below the original condition results. Open the section and click
**Explore reaction context**. This explicitly runs the existing forward and
single-step retro engines on the complete reaction's respective sides.

Ordinary condition retrieval, chemistry gates, recipe scores, and independent
support counts are unchanged. The new action is on demand because loading and
applying operator libraries costs substantially more than stored side lookups.

## Chemistry and application boundaries

- `condition_recommender/reaction_context.py` owns the advisory relation contract.
  It uses taxonomy parsing, featurization, and the existing full shared-core
  comparison. It retains unresolved and conflicting graph evidence.
- `chem_coworker/reaction_context.py` composes the existing forward engine,
  single-step retro engine, and canonical condition recommender. No planner is
  called recursively from condition retrieval.
- The HTTP runtime selects artifacts; React displays the results. Neither owns
  new chemical matching rules.

Forward analysis is target-blind enumeration followed by comparison with the
requested product. It retains generated products, site/operator competition
groups, checks, precedent IDs, and search diagnostics. No operator hint or
recipe is supplied. Thus it describes structural possibilities, **not
recipe-specific selectivity, product probability, or experimental feasibility**.
A missing target is inconclusive; an out-of-scope engine result remains visible.

Retro proposals are independently graph-validated by the existing strategy
search, then compared with the original reaction:

| Relation | Condition behavior |
| --- | --- |
| Supplied reaction recovered | Refer to the original recommendations; add no independent evidence. |
| Shared transformation with changed inputs | Query conditions for the explicitly displayed precursor alternative. |
| Same product through a different graph transformation | Keep as an alternative synthesis route; query its own conditions. |
| Different product | Do not request alternative condition evidence. |
| Unresolved/conflicting core relationship | Retain the proposal with uncertainty; do not request alternative condition evidence. |

Generated structures always carry `evidence_kind: generated_hypothesis`.
Operator source IDs identify supporting observations, not observations of the
generated reaction itself. Alternative recipes keep the canonical recommender's
matching level, source requirements, cautions, and precedent IDs. A graph match
does not establish that changing substrates or sources preserves a procedure.

Duplicate molecular proposals merge source/operator provenance and trigger only
one condition query. Unverified proposals are excluded. Forward/retro agreement
does not increase recipe scores or experimental support.

## Limits and failure behavior

Validated `reaction_context.v1.json` sets the limits: 40 forward operators,
32 assignments and 64 outcomes per operator, 10 displayed forward products;
40 retro templates and 12 validations **per attempted tier**, four strategies,
two realizations per strategy; four unique displayed proposals and up to three
recipes per qualified alternative. Representatives receive display slots before
additional realizations. These are work budgets, not a wall-clock deadline.

The endpoint uses complete reactants and products. Multiple product components
receive `SINGLE_PRODUCT_REQUIRED`; it never silently chooses one product.
Standalone forward/retro tools remain the entry points for partial queries.
Core relationship analysis uses deterministic taxonomy evidence without adding
an external mapping dependency. If that evidence is insufficient, the
relationship remains unresolved even if another workflow obtained a mapping.

Unavailable forward artifacts do not prevent retro analysis, and unavailable
condition evidence does not erase valid precursor proposals. A stale or absent
prepared forward artifact is reported explicitly; this action never builds an
operator library during a request.

## API and artifact impact

`POST /api/v1/recommendations/context`, Workbench profile only:

```json
{
  "reaction_smiles": "Ic1ccccc1.N#C[Cu]>>N#Cc1ccccc1",
  "library_mode": "full"
}
```

The existing API envelope contains a context response with schema `1.0`, policy
ID, search limits, advisory flag, separate reactant/product analyses, warnings,
and proposal-specific condition results. Existing recommendation contracts are
unchanged. The focused Condition Desk deployment does not expose this research
endpoint.

**No dataset rebuild is required.** This consumes the current prepared forward,
retro, and condition artifacts without modifying their identities or admission.
The frontend was rebuilt. Restart `python -m app.web_api --workbench` to load the
new backend and browser bundle.

## Verification and practical limits

Evidence is under `results/reaction_context_20260915/` (local generated output).
Unit/integration tests cover real operator application, source/departure
analogues, partner/map-order invariance, different constructions of the same
product, wrong products, unavailable/conflicting correspondence, unverified
proposals, deduplication, artifact failures, and HTTP deployment boundaries.

Live Full and Compact browser checks use aryl-iodide cyanation. Both return one
precursor alternative with three condition recipes and three separate synthesis
routes without qualified recipes. The original condition table remains
unchanged. Both forward audits return no generated product within the bounded
search. In Full, the coarse precursor index nominates 5,087 operators and the
40 attempted operators yield no outcomes. The existing diagnostic
`applied_operator_count` counts operators producing outcomes, so its zero does
not mean the action skipped forward search.

A separate Compact ethyl-bromide/ammonia query also finds no forward product
under these limits; retro returns alternative routes, three with condition
recipes. The small deterministic fixture confirms a valid forward/retro engine
path independently of the large-library coverage limitation.

This establishes operational integration and separation of evidence, **not
better recall or more accurate conditions**. The next evaluation should measure
target recovery, competing-site coverage, useful precursor alternatives, and
qualified condition support separately on a frozen diverse query set. Forward
applicability filtering before the application budget is a concrete remaining
limitation: coarse retrieval plus support ordering can spend the bounded search
on operators that do not match the actual input structures. Increasing budgets
alone should not be treated as a demonstrated chemistry improvement.

Final verification:

- Complete Python suite: **1,535 passed in 380.38 seconds** (`full_suite.log`).
- Production frontend build passed in 22.07 seconds; the existing large-bundle
  warning remains (`frontend_build.log`).
- Full/Compact Playwright checks passed in 12.7 seconds against the final bundle
  (`browser.log`). Both screenshots were inspected; no page JavaScript errors.
- A real single-observation operator fixture recovered ethylamine exactly and
  identified the supplied reaction, without requesting duplicate condition
  evidence (`positive_fixture.json`).
- Ruff checks and `git diff --check` passed.

The temporary Workbench server on port 8027 was used only for validation.
