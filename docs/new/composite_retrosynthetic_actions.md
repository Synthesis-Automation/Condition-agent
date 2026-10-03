# Composite retrosynthetic actions

Implemented 2026-10-03; experimental chemistry scope. Independent chemistry review
and untouched evaluation remain release gates.

The Workbench **Composite two-step strategies**, Python API, CLI and scientific
agent use `core_retrosynthesis.composite_actions.search_composite_actions`.
One logical action expands into two predicted physical reactions. The intermediate,
additional partners, precedents, validation and compatibility evidence stay explicit.
Both reactions count toward physical route cost. Conditions, one-pot execution and
material supply remain `not_assessed`; graph consistency does not establish feasibility.

## Chemistry gates

Existing generic single-step operators perform the graph transformations. Before
composite scoring, both steps must pass structural validation and compatibility.
The predicted `ReactionRouteTree` is projected through the existing route-core analysis
and coupled-route classifier. Query-site coupling must agree with the catalogue's
structural relationship. Complete symmetry alternatives must give invariant evidence.
Independent sites, invalid maps, unresolved chemistry, ambiguous lineage and conflicting
relationships cannot be admitted by frequency or reaction labels.

`dependency_reviews` preserves rejected hypotheses and evidence (up to 20 entries);
`diagnostics.dependency_rejected_count` counts all rejected hypotheses. Rejection
does not erase ordinary one-step fallbacks. Accepted actions retain both generic
reaction signatures, atom references and signature schema/definition versions in
`dependency.physical_step_signatures`. Each physical step retains its source IDs,
selectivity warnings and precursor/reaction compatibility assessments.

`definitions/composite_actions.v1.json` owns validated gates, scoring weights,
support saturation, review bounds and physical cost. Executable graph behavior stays
in Python. No arbitrary executable JSON hooks or removed chemistry imports exist.

Regressions cover alcohol mesylation/substitution and aromatic nitration/reduction.
These are scoped operators, not universal ROH-to-RNu or ArH-to-ArNH2 rules.
Mesylation/substitution is classified as `created_handle_consumed`: the installed
O–S bond stays in the departing group while substitution breaks C–O. The current
`activation_then_conversion` classifier requires breaking a bond installed in step 1.
Chemical names do not override this distinction. Broader scope requires admitted
constituent operators and supported strategy definitions.

## Standalone catalogue and entry points

`CompositeStrategyCatalog` contains reusable `CompositeStrategyDefinition` entries,
training-patent support, versions and provenance. It contains no held-out cases or
evaluation target/intermediate structures. Content identity is validated on load and
is invariant to strategy/support ordering.

Export existing experimental panel definitions once:

```powershell
python -m core_retrosynthesis export-composite-catalog `
  results/core_retrosynthesis/coupled_strategy_evaluation/route_only_fixed_panel.v1.json `
  results/core_retrosynthesis/coupled_strategy_evaluation/composite_strategy_catalog.v1.json
```

The default local catalogue starts with the same 12 strategies. For broader coverage,
`build-composite-catalog SOURCE_ROUTE_CORES LIBRARY OUTPUT_JSON` exports all recurrent
training-supported pairs covered by the library, rather than sampling a panel.
Training support stays distinct from review approval. Inspect source reviews before
expanding admitted scope; neither command upgrades chemistry-review status.

```python
from core_retrosynthesis import (
    load_composite_strategy_catalog, load_generic_library, search_composite_actions,
)

catalog = load_composite_strategy_catalog("catalog.json")
library = load_generic_library("operators.json.gz")
result = search_composite_actions(target_smiles, library, catalog.strategies)
for action in result.actions:
    print(action.action_id, action.terminal_precursor_smiles, action.physical_step_cost)
    for step in action.retrosynthetic_steps:
        print(step.reaction_smiles, step.precedent_reaction_ids)
```

```powershell
python -m core_retrosynthesis disconnect-composite `
  PATH_TO_OPERATOR_LIBRARY PATH_TO_COMPOSITE_CATALOG TARGET_SMILES --top-k 5
```

See `disconnect-composite --help` for budgets, context, L0 and one-step fallback
options. Missing operators are capability gaps, not chemical impossibility. The
planner must check cycles, leaf supply and full physical cost when selecting actions.

The agent operation is `disconnect_composite(target_smiles=..., top_k=...)`. It uses
separately baseline-pinned `composite_library` and `composite_catalog` artifacts,
included in `examples/ai_native/artifacts.local.example.json`. The experimental
library stays separate from ordinary `retro_library`. Saved calls and summaries
retain both physical steps, dependency evidence and deterministic replay. The agent
owns branch selection/stopping; the tool performs bounded pair expansion. Assess
and inspect selected physical steps individually before finalizing a route.

## Workbench migration

Restart the server and start a new scientific investigation for the new code/baseline.
The endpoint remains `POST /api/v1/retrosynthesis/coupled-strategies`. It calls the
shared service and reads the standalone catalogue via
`CORE_RETROSYNTHESIS_COMPOSITE_CATALOG` or `coupled_strategy_catalog_path`.
The former `CORE_RETROSYNTHESIS_COUPLED_PANEL` / `coupled_strategy_panel_path`
configuration is removed; export a catalogue before switching an existing deployment.

Query schema is 1.1 / `composite_actions.v1`. `catalog_id` replaces `panel_id`, and
capabilities expose `coupled_strategy_catalog_name` rather than the panel name.
The existing workbench artifact discriminator remains for saved-result reading.
Flat first/second reaction fields are compatibility views of `physical_steps`, with
parity regressions. Remove those views and the v1 discriminator when Workbench
render/export readers migrate entirely to the physical-step list.

The former `search_promoted_v1_strategies` function is removed. Evaluation calls the
shared service; there is no separate promoted-query implementation. Exact observed
replay remains an audit operation. Evaluation advances to 1.1 /
`coupled_strategy_eval.v2` because dependency gates can change admission. Historical
metrics are unchanged and do not establish the new version's performance.

## Validation scope

Tests cover both positive sequences, partners/intermediates, independent sites,
symmetry ambiguity, invalid maps, structural/relationship conflicts, failed
compatibility, deterministic identities, serialization invariance, fallbacks,
catalogue isolation, agent baseline/replay and real Workbench/CLI parity. Source
round trips establish graph consistency, not prospective selectivity or experimental
success. Existing held-out records are not used to tune this feature.
