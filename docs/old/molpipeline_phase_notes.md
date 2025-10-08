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

# Reaction-Centric Feature Engineering (Phase 2)

- Hooked MolPipeline aggregation into `chemtools.precedent._attach_molpipeline_features`, wiring per-role feature vectors and concatenated embeddings onto each precedent row when a `molpipeline` relax block is present.
- Normalise incoming config (roles, aggregator strategy, job counts, fingerprint sizes) and propagate `query_role_smiles` so the query reaction can emit its own feature vector alongside precedent candidates.
- Surface integration warnings instead of hard failures by default (`suppress_errors`), keeping precedent search deterministic even when MolPipeline raises.

# Optional Dependency Packaging (Phase 3)

- Registered a `molpipeline` optional extra in `pyproject.toml` so `pip install chemtools-project[molpipeline]` pulls the RDKit/MolPipeline stack on demand.
- Added `tests/test_molpipeline_integration.py` and `tests/test_precedent_molpipeline.py` behind `pytest.importorskip("molpipeline")` to exercise the new aggregator while staying green when the extra is absent.
- Gate every public entry point through `chemtools.integrations.molpipeline.environment_snapshot()` so downstream callers can report availability/versions without importing optional binaries eagerly.

# UI & CLI Exposure (Phase 4)

- Updated `app/ui_gradio.py` precedent tab with a MolPipeline accordion: toggle attachment, edit JSON config/query snippets, and render a live summary (env metadata, warnings, feature lengths) returned from `ui_precedent_search`.
- Extended `ui_precedent_search` to accept MolPipeline toggles/config from the UI, sanitise them, forward into the precedent relax dict, and return a structured status payload for display.
- Expanded `scripts/precedent_from_rxn.py` with `--molpipeline`, `--molpipeline-config`, and `--molpipeline-query` flags; CLI runs now inject the config into the relax dict and emit a `molpipeline` summary block (availability, config echo, vector lengths, warnings) alongside the precedent pack.
