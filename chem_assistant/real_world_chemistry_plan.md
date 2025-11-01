# Roadmap Toward Real-World Chemistry Recommendations

## Current Gaps

- Precedent-only suggestions lack broad datasets, uncertainty estimates, and lab feedback.
- No integration with ELNs, reagent inventory, or literature updates.
- Agent workflow does not validate feasibility (toxicity, solubility, selectivity, availability).
- Limited support for iterative optimisation or multi-objective trade-offs (yield, cost, sustainability).

## Tooling Additions

1. **Data Connectors**

   - ELN/HTE ingestion pipelines with metadata normalisation.
   - Reagent inventory + pricing APIs; vendor availability checks.
   - Literature/patent summariser (e.g., SciFinder, Reaxys) with citation capture.
   - Live precedent fetchers (open databases, statutory corpuses).

2. **Property & Feasibility Models**

   - Solubility/miscibility predictors for solvent choice.
   - pKa, logP, redox, and stability estimators for substrates/intermediates.
   - Reaction risk flags (peroxide formation, over-alkylation, exotherms).
   - Selectivity conflict detector for protecting groups and functional handles.

3. **Advanced Recommenders**

   - Catalyst/ligand selector balancing performance, availability, and cost.
   - Solvent swapper with green-chemistry scoring.
   - Base/temperature/pressure optimiser respecting user constraints.
   - Multi-objective planner (yield × cost × sustainability × tox).

4. **Uncertainty & Validation**

   - Similarity dispersion + precedent diversity metrics.
   - Ensemble confidence intervals for recommendations.
   - Rule-based plausibility checks (functional-group incompatibilities, safety).
   - “What-if” scenario simulator leveraging kinetic heuristics/HTE data.

5. **Feedback & Continuous Learning**
   - Experiment result ingester (structured + free text) to weight precedents.
   - Active-learning experiment queue generator.
   - Lab notebook summariser for human annotations.

6. **Registry Automation**
   - Agent-accessible API for dry-run and persistent reagent additions.
   - LLM-guarded workflows that confirm intent before writing to the taxonomy.
   - GUI hooks (e.g., `reagent_addition_ui.py`) wired through the shared service layer.

## Workflow Enhancements

1. **Goal-Oriented Agent Loop**

   1. Normalise inputs → detect family → parse constraints.
   2. Generate candidate conditions via multi-tool planner.
   3. Validate (rule checks, inventory, property screens).
   4. Rank with confidence + trade-off explanations.
   5. Present structured output (tables, rationale, caveats).

2. **Interactive Session Memory**

   - Persist catalyst/base locks, cost ceilings, sustainability targets across turns.
   - Offer alternative scenarios (Pd-free, low temp, low cost) on demand.

3. **Human-in-the-Loop Review**

   - Provide concise rationale, uncertainty badges, and precedent snapshots.
   - Bundle export-ready summaries for ELN insertion.

4. **Automation Hooks**
   - Optional export to robo-lab schedulers or DoE planners.
   - Alerts when new relevant precedents or literature updates appear.

## Initial Implementation Steps

1. Prototype inventory connector, property predictor, and uncertainty scorer.
2. Extend LangChain workflow with validation and richer response formatting.
3. Pilot on a high-impact family (e.g., C–N coupling) using historical ELN data to evaluate accuracy uplift and gather feedback.

## Success Metrics

- Reduction in manual triage time for condition selection.
- Increase in first-pass experimental success rate.
- Positive chemist feedback on explanation clarity and constraint handling.
- Demonstrated ability to learn from new experiment data within one sprint.
