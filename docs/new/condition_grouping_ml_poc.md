# Weak-label condition grouping ML proof of concept

## Purpose

This proof of concept tests whether the structured conditions in
`datasets/reaction_label/v2.1_cleaned.csv` can support role-aware statistical
condition groups. It does not change production recipe identity, admission,
retrieval, compatibility, or ranking.

The POC preserves complete normalized material systems as `CMAT1` identities,
factors exact decisive systems into `CCOREPOC1` identities, and generates
experimental, snapshot-bound `CGPOC2` groups. These groups are review proposals,
not claims that their recipes are interchangeable.

The factorization follows the chemist-facing presentation illustrated by
`knowledge_base/roche_duck_reagent_conditions_summary.html`: an amidation is
commonly summarized by coupling reagent/base (for example HATU/DIPEA), while a
Suzuki condition is commonly summarized by catalyst/base with an explicit
ligand subsystem where present. Solvent and operating details remain protocol
context rather than fragmenting the core identity.

## Feature boundary

The grouping model consumes only the decisive core:

- catalyst and ligand;
- coupling or condensation reagent;
- acid, base, oxidant, and reductant;
- decisive-core role topology;
- catalyst-ligand pair tokens;
- catalyst-base and activation-reagent/base combination tokens;
- explicit `No ligand` only when a catalyst is actually present.

Solvents, additives, other contextual components, temperature, and time are not
clustering features. They remain attached to groups and exact material
assignments as cross-references. Reaction type, reactive-site labels, yield,
z-score, and procedure text are excluded and may appear only in post-grouping
diagnostics.

The implementation uses TF-IDF, truncated SVD, and mini-batch k-means. Duplicate
support is square-root weighted so common observations contribute evidence
without completely dominating rare material systems.

## Run

```powershell
python -m condition_recommender.condition_grouping_poc_cli `
  datasets/reaction_label/v2.1_cleaned.csv `
  results/condition_grouping_poc/v2_core `
  --clusters 256 `
  --latent-dimensions 48 `
  --seed 17
```

The output contains:

- `condition_cores.jsonl`: exact decisive cores with observed protocol variants,
  solvent/additive frequencies, operating summaries, and learned-group links;
- `condition_groups.jsonl`: a statistical medoid, the most prevalent exact core,
  role-consensus coverage, core support, solvent/additive cross-references,
  operating summaries, and representative cores;
- `condition_group_assignments.csv`: exact material-to-group assignments with
  core identity, solvent/additive cross-references, centroid similarity,
  assignment margin, and an explicitly uncalibrated
  `supported`/`provisional`/`review` status;
- `condition_grouping_model.joblib`: the snapshot-bound POC model;
- `condition_grouping_report.json`: machine-readable run and quality report; and
- `condition_grouping_report.md`: concise human review summary.

## Complete protocol context from a core

The completion API first looks for an observed exact `CCOREPOC1`. If found, it
returns solvent/additive choices and complete material variants reported with
that exact core. If the core is unseen, it may return context from the nearest
`CGPOC2`, explicitly labeled as a broader learned prior rather than an observed
protocol.

```powershell
python -m condition_recommender.condition_completion_poc_cli `
  results/condition_grouping_poc/v2_core `
  "K3PO4 [base]; dppfPdCl2 [catalyst]; dppf [ligand]" `
  --top-k 5
```

Completion does not bypass reaction/recipe compatibility filters. Frequency is
reported as precedent support and must not be presented as proof of optimality.

## Promotion boundary

`CCOREPOC1` and `CGPOC2` must not become production identity namespaces.
Promotion requires
chemist review, explicit split/merge decisions, validated condition-family
definitions, stable definition-derived IDs, and regression coverage. A future
production condition-system identity belongs in `condition_registry`; the
statistical model may propose assignments but must preserve ambiguity and
confidence.
