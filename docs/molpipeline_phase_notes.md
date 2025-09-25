# MolPipeline Due Diligence (Phase 0)

- **License:** MIT (compatible with ChemTools) - upstream notice must accompany redistributed components.
- **Core dependencies:** joblib>=1.3, loguru>=0.7, matplotlib>=3.10, numpy>=1.26, pandas>=2.2, rdkit>=2023.9, scikit-learn>=1.6, scipy>=1.15, shap>=0.47, typing-extensions>=4.13.
- **Optional extras:** `chemprop` (deep models, via lightning) and notebook tooling; not required for baseline integration.
- **Python support:** upstream targets Python 3.11 (`pyproject.toml`); aligns with ChemTools runtime matrix as long as RDKit 2023.09+ is used.
- **Runtime considerations:** SHAP adds gradient-based explainability (GPU optional). Pipelines leverage joblib for multiprocessing; respect ChemTools determinism via explicit `n_jobs` configuration.
- **Recommendation:** expose MolPipeline features behind optional extras flag so base installs avoid the heavier scientific stack.

# Prototype Notes (Phase 1)

- Added `chemtools.integrations.molpipeline` providing:
  - Role-aware aggregators for combining MolPipeline outputs per reagent role.
  - Environment probe (`environment_snapshot`) to surface versions/availability.
  - Standardized builders for MolPipeline physicochemical and Morgan pipelines using ChemTools-friendly defaults.
  - Input coercion helpers that accept raw SMILES strings or ChemTools reagent dicts and recover normalized SMILES via `chemtools.smiles` when needed.
  - Thin wrappers (`transform_smiles`, `fit_pipeline`) keeping notebooks/CLI consistent without importing MolPipeline unless available.
- Exposed package via `chemtools.integrations` export for downstream imports.
- Smoke-tested descriptor and fingerprint pipelines against benzene to confirm operational wiring.
- Optional extras (`pip install chemtools-project[molpipeline]`) wire MolPipeline into the install flow.

- The `precedent.knn` relax dict now accepts a `molpipeline` block (e.g., `{"molpipeline": {"include_role_features": True}}`) to attach aggregated per-role MolPipeline vectors to each returned precedent and an optional `molpipeline_query_vector` for the input reaction.
