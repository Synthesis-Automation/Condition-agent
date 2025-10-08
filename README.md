# chemtools-project (v2)

Deterministic chemistry tools (FastAPI + RDKit-friendly) to support condition recommendation. Focuses on Ullmann first, with simple hooks for more families.

## 🆕 Recent Updates (October 2025)

**Recommendation Module Refactoring - Complete** ✅

We've successfully refactored the recommendation system for better maintainability and performance:

- ✅ **59% code reduction**: Main recommendation engine reduced from 1,454 → ~600 lines
- ✅ **DRFP-first approach**: Focus on reaction similarity (DRFP fingerprinting) instead of complex substrate featurization
- ✅ **Modular structure**: New `chemtools/recommend/` package with focused modules
- ✅ **Better ML organization**: ML code grouped in `chemtools/ml/` folder
- ✅ **100% backward compatible**: Old imports still work (with deprecation warnings)

**Migration**: No action required! Your existing code continues working. For details, see:

- `RECOMMEND_REFACTORING_SUCCESS.md` - Executive summary and test results
- `MIGRATION_GUIDE.md` - Migration scenarios and best practices

## Quickstart

- Create a virtualenv and install deps:

  - macOS/Linux:

    ```bash
    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    ```

  - Windows (PowerShell):

    ```powershell

 Where it lives

- Core function: `chemtools.recommend.recommend_from_reaction` (legacy API: `POST /api/v1/recommend`).
  - **New modular structure**: Refactored into `chemtools/recommend/core.py` (59% smaller, DRFP-first)
  - **Backward compatible**: Old import `from chemtools.recommend import recommend_from_reaction` still works (shows deprecation warning)
- Structured API: `chemtools.recommend.recommend_conditions_structured` (`POST /api/v1/recommend/conditions`) returns top-N condition sets (default 5) with metadata for clients.
- ML-enhanced: `chemtools.ml.recommender.hybrid_recommend` for yield-prediction-based re-ranking
  - **Previously**: `chemtools.recommend_ml.hybrid_recommend` (now deprecated, but backward compatible)
- UI: Gradio tab "Recommend Conditions" now surfaces the structured JSON, a table of the top 5 core/base/solvent combos, detection notes, and accepts an optional reaction-type override.thon -m venv .venv
    .\.venv\Scripts\Activate.ps1
    pip install -r requirements.txt

    ```

- Run API server:

  ```bash
  uvicorn app.main:app --reload --port 8000
  # open http://127.0.0.1:8000/docs
  ```

- Run tests:

  ```bash
  pip install -r requirements-dev.txt
  pytest -q
  ```

  Or with Makefile shortcuts: `make install`, `make run`, `make test`.

## Repository Layout

`data/` and `data-processor/` hold raw and processed datasets for demos and ETL pipelines; they remain unchanged and are omitted from the listing below so the focus stays on the runtime code.

```text
.
├── app/
│   ├── main.py                 # FastAPI service wiring and HTTP endpoints
│   └── ui_gradio.py            # Gradio UI surface for interactive testing
├── chemtools/                  # Core chemistry tools package
│   ├── features/               # Role-aware descriptors and fingerprint builders
│   │   └── role/               # SMARTS patterns and role-specific features
│   ├── featurizers/
│   │   ├── molecular.py        # Primary molecular featurizer API
│   │   └── ullmann.py          # Backwards-compatible Ullmann-specific featurizer
│   ├── integrations/
│   │   ├── molpipeline.py      # Optional MolPipeline feature hooks
│   │   └── mcp/
│   │       └── tools/          # Schema-stamped wrappers exposed as MCP tools
│   ├── ml/                     # Machine learning modules (refactored Oct 2025)
│   │   ├── drfp_predictor.py   # DRFP-based yield predictor (LightGBM)
│   │   ├── evaluation.py       # ML evaluation metrics and plotting
│   │   └── recommender.py      # ML-enhanced recommendation (hybrid ML + k-NN)
│   ├── precedent/              # Precedent search and similarity (refactored)
│   │   ├── search.py           # Main k-NN search with DRFP integration
│   │   ├── similarity.py       # DRFP similarity computation
│   │   ├── loader.py           # Dataset loading and caching
│   │   ├── catalyst.py         # Catalyst system utilities
│   │   ├── core_utils.py       # Core condition helpers
│   │   └── integrations.py     # External integrations (MolPipeline, etc.)
│   ├── recommend/              # Modular recommendation package (NEW Oct 2025)
│   │   ├── __init__.py         # Public API exports
│   │   ├── core.py             # DRFP-first recommendation engine (~600 lines)
│   │   ├── utils.py            # Shared utilities (canonical_family, median)
│   │   ├── substrate_analysis.py  # Optional functional group detection
│   │   ├── plate_design.py     # Plate design generation
│   │   └── families/           # Family-specific recommendation logic
│   ├── rule_scdb_matcher/      # Rule-based matcher (SchemeConditionDB)
│   │   ├── cli.py              # CLI for rule-based matching
│   │   ├── loader.py           # Database loader
│   │   ├── matcher.py          # Core matching engine
│   │   └── types.py            # Type definitions
│   ├── util/                   # Shared helpers (RDKit guards, IO utilities)
│   ├── condition_core.py       # Condition-core parsing and normalization
│   ├── constraints.py          # Post-processing filters for recommendations
│   ├── context.py              # Global context and resource management
│   ├── contracts.py            # Pydantic request/response contracts
│   ├── explain.py              # Human-readable justification builders
│   ├── output_formatter.py     # Output formatting utilities
│   ├── properties.py           # Name/CAS/token property lookup
│   ├── reaction_similarity.py  # DRFP/Tanimoto similarity utilities
│   ├── reaction_type_detector.py  # Reaction family inference
│   ├── reagent_lookup.py       # Reagent database lookup
│   ├── recommend_DEPRECATED.py # Backward-compat wrapper (→ recommend/ package)
│   ├── recommend_ml_DEPRECATED.py  # Backward-compat wrapper (→ ml/recommender.py)
│   ├── router.py               # Reaction family detection heuristics
│   ├── selector_payloads.py    # Selector payload builders
│   └── smiles.py               # SMILES normalization and reaction parsing
├── docs/                       # Documentation
├── scripts/                    # Developer CLIs (DRFP caches, ML training, etc.)
│   └── ml/                     # ML-specific scripts (train, evaluate, prepare)
├── tests/                      # Test suite
│   ├── chemtools/              # Package-focused unit tests
│   └── test_*.py               # Integration tests (API, precedent, router, ML)
├── results/                    # ML evaluation results and plots
├── AGENTS.md                   # Agent integration design notes
├── Makefile                    # Convenience targets (install/run/test)
├── pyproject.toml              # Project metadata and entry points
├── requirements.txt            # Core dependencies
├── requirements-ml.txt         # ML dependencies (DRFP, LightGBM, etc.)
├── requirements-dev.txt        # Development dependencies
├── MIGRATION_GUIDE.md          # Migration guide for refactored modules
├── RECOMMEND_REFACTORING_SUCCESS.md  # Refactoring summary
└── README.md                   # This guide
```

**Key Changes (October 2025 Refactoring)**:
- ✅ `chemtools/recommend/` is now a modular package (was single `recommend.py`)
- ✅ `chemtools/ml/recommender.py` for ML-enhanced recommendations (was `recommend_ml.py`)
- ✅ `chemtools/precedent/` refactored into focused modules (was single `precedent.py`)
- ✅ Deprecated files kept for backward compatibility with warnings

**Supporting data directories**:
- `data/`: Reaction datasets, compound taxonomy, and SchemeConditionDB snapshots
- `data-processor/`: ETL utilities for building/updating datasets

## How to Test

### 1) Unit tests (fast, deterministic)

- Run all: `pytest -q`
- Notes: tests disable heavy RDKit paths automatically; no extra setup needed.

### 2) Exercise the API via Swagger

- Start: `uvicorn app.main:app --reload --port 8000`
- Open: <http://127.0.0.1:8000/docs>
- Try endpoints:
  - POST `/api/v1/smiles/normalize` with `{ "smiles": "c1ccccc1O" }`
  - POST `/api/v1/router/detect-family` with `{ "reactants": ["Clc1ccccc1","Nc1ccccc1"] }`
  - POST `/api/v1/featurize/molecular` with `{ "electrophile": "Clc1ccccc1", "nucleophile": "Nc1ccccc1" }` (alias: `/api/v1/featurize/ullmann`)

### 3) Quick Python one-liners

- Normalize: `python -c "from chemtools.smiles import normalize; print(normalize('c1ccccc1O')['smiles_norm'])"`
- Detect family: `python -c "from chemtools.router import detect_family; print(detect_family(['Clc1ccccc1','Nc1ccccc1']))"`
- Featurize: `python -c "from chemtools.featurizers.molecular import featurize; print(featurize('Clc1ccccc1','Nc1ccccc1'))"`

### 4) Gradio UI (no-code testing)

- Launch: `python app/ui_gradio.py`
- Open: <http://127.0.0.1:7860>
- Tabs included:
  - SMILES Normalize: normalize a single SMILES.
  - Detect Family: infer reaction family from reactants (dot or newline separated).
  - Single Molecule (basic): globals-only vector by default; optionally add roles (amine/alcohol/aryl_halide).
  - Properties Lookup: query by name/CAS/token (e.g., "Water", "7778-53-2"). Uses the compound taxonomy under `data/compound_taxonomy` as the default source (no dependency on a merged CAS registry).
  - Recommend Conditions: runs the recommender over a reaction SMILES.
  - Design Plate: proposes a plate CSV across top cores.
  - Precedent Search: retrieves similar precedents (optionally DRFP re-ranking).
  - DRFP Similarity: Tanimoto between two reaction SMILES (if DRFP installed).
  - Core Search: find dataset reactions by condition core (e.g., `Pd/XPhos` or `XPhos`). Enter a core token, optionally select a family (e.g., "Suzuki_CC", "Ullmann C–N"), keep "Fuzzy" on to match ligand names in catalyst systems, then click Search.

Optional: install DRFP to enable similarity re-ranking and the similarity tab

```
pip install drfp==0.4.0 numpy
```

## MolPipeline Integration

- Optional: `pip install chemtools-project[molpipeline]` to pull MolPipeline, RDKit, and friends on demand. Without it the app keeps running, but MolPipeline-specific toggles stay disabled.
- The Precedent Search tab exposes a MolPipeline accordion; enable it to attach per-role fingerprints/physchem descriptors to each precedent and to inspect environment metadata or warnings returned by `ui_precedent_search`.
- API callers can pass a `molpipeline` block inside the `/api/v1/precedent/knn` relax payload (roles, aggregation strategy, query role overrides) to get the same feature vectors programmatically.
- The CLI helper `scripts/precedent_from_rxn.py` mirrors the UI via `--molpipeline`, `--molpipeline-config`, and `--molpipeline-query`, echoing availability, applied config, and vector lengths in the output.
- All entry points rely on `chemtools.integrations.molpipeline.environment_snapshot()` so they degrade gracefully and can report version info when the extra stack is present.

Gradio Core Search quick demo

- Start Gradio: `python app/ui_gradio.py`
- Open the "Core Search" tab.
- Try core `Pd/XPhos` (or just `XPhos`) and click Search.
- Optionally set Family to `Suzuki_CC` or `Amide_Coupling` to narrow results.
- To browse cores, pick a Family (optional) and click "List cores" to see unique cores with counts.

Tips

- Speed up or run without RDKit: set `CHEMTOOLS_DISABLE_RDKIT=1`
  - macOS/Linux: `export CHEMTOOLS_DISABLE_RDKIT=1`
  - Windows (PowerShell): `$env:CHEMTOOLS_DISABLE_RDKIT='1'`
- If running `python app/ui_gradio.py` outside the repo root, prefer module mode: `python -m app.ui_gradio`

### Data Sources (Registry/Properties)

- By default, chemtools loads compound roles, names, and tokens from the taxonomy under `data/compound_taxonomy` (ligands, bases, solvents, coupling reagents, catalyst precursors). The merged CAS registry JSONL is not required for Properties Lookup or name resolution in the UI/API.
- To point to an alternate taxonomy location, set `CHEMTOOLS_TAXONOMY_DIR`.

## Endpoints

- `GET /health`: health probe.
- `POST /api/v1/smiles/normalize` 鈥?Normalize SMILES.
  - In: `{ "smiles": "Brc1ccc(F)cc1.O" }`
  - Out: `{ "input", "fragments", "largest_smiles", "smiles_norm" }`
- `POST /api/v1/router/detect-family` 鈥?Detect reaction family.
  - In: `{ "reactants": ["Brc1ccc(F)cc1", "Nc1ccccc1"] }`
  - Out: `{ "family", "confidence", "hits" }`
- `POST /api/v1/featurize/molecular` (alias: `/api/v1/featurize/ullmann`) 鈥?Substrate features for Ullmann C 鈥揘.
  - In: `{ "electrophile": "...", "nucleophile": "..." }`
  - Out: feature dict including `LG`, `nuc_class`, `bin`, 鈥?
- `POST /api/v1/condition-core/parse` 鈥?Parse ConditionCore from reagents/text.
  - In: `{ "reagents": [{"uid":"7681-65-4","role":"CATALYST","name":"CuI"}, {"uid":"72-52-8","role":"LIGAND","name":"Phenanthroline"}], "text": "optional free text" }`
  - Out: `{ "core": "Cu/Phenanthroline", ... }`
- `POST /api/v1/properties/lookup` 鈥?Minimal properties by UID/token.
  - In: `{ "query": "K3PO4" }` or `{ "query": "7778-53-2" }`
- `POST /api/v1/precedent/knn` 鈥?Retrieve similar precedents.
  - In: `{ "family": "Ullmann_CN", "features": {"bin":"LG:Br|NUC:aniline"}, "k": 50 }`
- `POST /api/v1/constraints/filter` 鈥?Inventory/blacklist filtering of candidate IDs.
- `POST /api/v1/explain/precedents` 鈥?Short reasons and example precedents.

- `POST /api/v1/core/search` – Find reactions by condition core (exact or fuzzy match).
  - In: `{ "core": "Pd/XPhos", "family": null, "fuzzy": true, "limit": 25 }`
  - Out: `{ "query": {..}, "count": 2, "results": [ {"reaction_id": "...", "condition_core": "Pd/XPhos", ...}, ... ] }`

Notes:

- RDKit is optional at runtime; if unavailable, SMILES normalization falls back to simple heuristics.
- Sample data lives under `data/` and is used by the precedent demo.

## Project Structure

- `app/main.py` 鈥?FastAPI app and routes.
- `chemtools/` 鈥?Deterministic tool implementations:
  - `smiles.py`, `router.py`, `properties.py`, `precedent.py`, `constraints.py`, `explain.py`
  - `condition_core.py` 鈥?normalizes ConditionCore
  - `featurizers/molecular.py` (alias kept: `featurizers/ullmann.py`) 鈥?Ullmann C 鈥揘 featurizer
  - `util/rdkit_helpers.py` 鈥?RDKit helpers (safe import)
- `tests/` 鈥?Lightweight unit tests
- `data/` 鈥?Small JSONL samples for precedents
- `chemistry_tool_bulild_plan.md` 鈥?Build plan and roadmap

## Development Notes

- API surfaces match the contracts in `chemtools/contracts.py`.
- Performance-sensitive code paths are simple and deterministic; add caching where beneficial.
- See the build plan for roadmap and design details: `chemistry_tool_bulild_plan.md`.

## Machine Learning Yield Prediction (Phase 1)

**NEW**: Predictive yield modeling using DRFP fingerprints and LightGBM, avoiding hand-crafted features.

### Quick Start

1. Install ML dependencies:

   ```bash
   pip install -r requirements-ml.txt
   ```

2. Prepare dataset (extract from `data/reaction_dataset/*.jsonl`):

   ```bash
   python scripts/ml/prepare_dataset.py
   ```

   This creates:
   - `data/ml_train.jsonl` (70%)
   - `data/ml_val.jsonl` (15%)
   - `data/ml_test.jsonl` (15%)

3. Train DRFP yield predictor:

   ```bash
   python scripts/ml/train_drfp_model.py
   ```

   Model saved to: `models/drfp_yield_v1.pkl`

4. Evaluate on test set:

   ```bash
   python scripts/ml/evaluate_models.py
   ```

   Results: `results/evaluation_results.png`, `results/test_metrics.json`

### Hybrid ML + k-NN Recommendation

Use ML predictions to re-rank condition recommendations:

```python
# New import (recommended)
from chemtools.ml.recommender import hybrid_recommend

# Or backward-compatible import (shows deprecation warning)
# from chemtools.recommend_ml import hybrid_recommend

results = hybrid_recommend(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    model_path="models/drfp_yield_v1.pkl",
)

for cond in results['recommended_conditions'][:3]:
    print(f"Core: {cond['core']}, Predicted Yield: {cond['predicted_yield_pct']:.1f}%")
```

**Strategy**:

- If `n_precedents >= 50`: Use ML to predict yields and re-rank
- If `n_precedents < 50`: Fall back to k-NN voting (more robust for sparse data)

### Architecture (Phase 1)

- **Features**: DRFP (4096-bit differential reaction fingerprints) + one-hot encoded conditions (core, base, solvent, T, time)
- **Model**: LightGBM gradient boosting (1000 estimators, depth 8)
- **Target**: MAE < 15%, inference < 100ms (CPU-friendly)
- **No hand-crafted features**: Generalizes to all reaction families

### ML Module Structure

```text
chemtools/
  ml/                      # ML modules (refactored for better organization)
    __init__.py            # ML module exports
    drfp_predictor.py      # DRFPYieldPredictor class
    recommender.py         # Hybrid ML + k-NN recommender (moved from recommend_ml.py)
    evaluation.py          # Metrics and plotting
  
  recommend/               # Modular recommendation package (NEW)
    __init__.py            # Public API exports
    core.py                # DRFP-first recommendation engine (~600 lines, 59% smaller)
    utils.py               # Shared utilities
    substrate_analysis.py  # Optional functional group detection
    plate_design.py        # Plate design generation

scripts/ml/
  prepare_dataset.py       # Data extraction and splitting
  train_drfp_model.py      # Training script
  evaluate_models.py       # Test set evaluation

tests/ml/
  test_drfp_predictor.py   # Unit tests

# Deprecated (backward compatible):
chemtools/recommend_DEPRECATED.py     # Old import wrapper (use chemtools.recommend instead)
chemtools/recommend_ml_DEPRECATED.py  # Old import wrapper (use chemtools.ml.recommender instead)
```

**Refactoring Summary** (October 2025):

- ✅ **59% code reduction**: `recommend.py` (1,454 lines) → `recommend/core.py` (~600 lines)
- ✅ **DRFP-first approach**: Removed 538 lines of complex substrate featurization
- ✅ **Better organization**: ML code grouped in `ml/` folder, recommendation split into modular package
- ✅ **100% backward compatible**: Old imports work with deprecation warnings
- See `RECOMMEND_REFACTORING_SUCCESS.md` and `MIGRATION_GUIDE.md` for details

### Performance Targets (Phase 1)

- MAE: 12-15% (vs 20-25% from k-NN baseline)
- Inference: <100ms per reaction (CPU)
- Training: ~5-10 min on 10K reactions

### Future Phases

- **Phase 2** (3-4 weeks): Neural network with learned condition embeddings → MAE 10-12%
- **Phase 3** (3-4 weeks): Graph Neural Network (GNN) end-to-end → MAE 8-10%

See `docs/ML_Updated_Plan_End_to_End_Learning.md` for full implementation plan.

## License

Internal project scaffold. Add a LICENSE if distributing.

## Role-Aware Featurization (API)

- POST `/api/v1/featurize/role-aware/molecule`
  - In: `{ "smiles": "Clc1ccccc1", "roles": ["aryl_halide"] }` (roles optional; default `["amine","alcohol","aryl_halide"]`)
  - Out: `{ vector: [float...], fields: [str...], masks: {amine|alcohol|aryl_halide}, meta: {centers, roles} }`
- POST `/api/v1/featurize/role-aware/reaction`
  - In: `{ "reaction": "Brc1ccccc1.Nc1ccccc1>>" }`
  - Out: `{ reactants: [ {smiles, vector, fields, masks, meta}, ... ] }`
- GET `/api/v1/featurize/role-aware/fields?roles=amine,alcohol,aryl_halide`
  - Out: `{ roles: [...], fields: [...], counts: {global, by_role, fingerprints, total}, registry: {...} }`

### Examples

- Molecule (aryl halide):
  - Request:

    ```bash
    curl -s -X POST http://127.0.0.1:8000/api/v1/featurize/role-aware/molecule \
      -H "Content-Type: application/json" \
      -d '{"smiles":"Clc1ccccc1","roles":["aryl_halide"]}' | jq '{fields: .fields[0:8], vector: .vector[0:8], masks: .masks}'
    ```

  - Output (sample):
    `{"fields":["global.MW","global.logP",...], "vector":[112.56,2.93,...], "masks":{"aryl_halide":1,"amine":0,"alcohol":0}}`
- Reaction (default roles for all reactants):
  - Request:

    ```bash
    curl -s -X POST http://127.0.0.1:8000/api/v1/featurize/role-aware/reaction \
      -H "Content-Type: application/json" \
      -d '{"reaction":"Brc1ccccc1.Nc1ccccc1>>"}' | jq '.reactants[0] | {smiles, len: (.vector|length), masks}'
    ```

- Field list for specific roles:
  - Request:

    ```bash
    curl -s "http://127.0.0.1:8000/api/v1/featurize/role-aware/fields?roles=amine,aryl_halide" | jq '.counts, (.fields[0:10])'
    ```

  - Use this to align downstream models to a stable field order.

## Performance & Tuning

- First run can be slow (60–120s) because it loads and parses the reaction dataset under `data/reaction_dataset/`, runs RDKit parsing/heuristics, and warms in‑memory caches. Subsequent calls are much faster.
- Keep the server/UI running to reuse caches instead of restarting between calls.

Speed tips

- Disable heavy RDKit paths (fast demo mode):
  - macOS/Linux: `export CHEMTOOLS_DISABLE_RDKIT=1`
  - Windows (PowerShell): `$env:CHEMTOOLS_DISABLE_RDKIT='1'`
- Skip attaching role‑aware vectors during dataset loading (faster Core Search/Precedents):
  - Default is off. Enable only when needed with:
    - macOS/Linux: `export CHEMTOOLS_ATTACH_ROLE_AWARE=1`
    - Windows (PowerShell): `$env:CHEMTOOLS_ATTACH_ROLE_AWARE='1'`
- Skip loading the large dataset (for quick dev/testing):
  - macOS/Linux: `export CHEMTOOLS_LOAD_DATASET=0`
  - Windows (PowerShell): `$env:CHEMTOOLS_LOAD_DATASET='0'`
  - Note: recommendation/precedent features depend on precedents; disabling the dataset may reduce results.
- Turn off DRFP re‑ranking if installed (saves fingerprinting time):
  - API: send `{"relax": {"use_drfp": false}}` in `/api/v1/recommend` or Precedent Search.
  - UI: uncheck “Use DRFP re‑ranking” (Precedent Search tab).
- Avoid “precompute all” for DRFP: keep precompute scope to `candidates` (default) or disable precompute to prevent warming all rows.
- Role‑aware vectors are attached for inspection but not used for candidate caching/scoring; computation still costs time. If you only need basic features, use the Single Molecule (basic) tab or the convenience endpoint `/api/v1/featurize/molecule` (roles omitted → globals only).

## Condition Recommendation Workflow

- Input: reaction SMILES (e.g., `Brc1ccccc1.Nc1ccccc1>>`).
- Normalize: split and canonicalize reactants/products (`chemtools.smiles.normalize_reaction`).
- Family: detect reaction family from reactants (`chemtools.router.detect_family`).
- Featurize: build substrate features with `chemtools.featurizers.molecular.featurize` and a coarse `bin` label; drop `role_aware` from features for caching.
- Precedents: retrieve neighbors with `chemtools.precedent.knn` (optionally DRFP re‑ranking controlled by `relax.use_drfp`, `drfp_weight`, `drfp_n_bits`, `drfp_radius`).
- Core: choose condition core by Laplace‑smoothed vote over precedents; compute a simple confidence from the vote share and support.
- Base/Solvent: pick most frequent candidates conditioned on the chosen core, filtered by `chemtools.constraints.apply_filter` using optional rules (e.g., `no_chlorinated`, `no_HMPA`, `aqueous_only`, `min_bp_C`).
- Temperature/Time: median of values from precedents for the chosen core with fallback to all precedents.
- Reasons: short textual explanations from `chemtools.explain.for_pack` (e.g., precedent support, alternative cores/bases/solvents).
- Output: a dict with `recommendation` (core/base/solvent/T/time/confidence), `alternatives`, `precedent_pack`, `reasons`, `filters`, and a human‑friendly `formatted` block.

Where it lives

- Core function: `chemtools/recommend.py:recommend_from_reaction` (legacy API: `POST /api/v1/recommend`).
- Structured API: `chemtools/recommend.py:recommend_conditions_structured` (`POST /api/v1/recommend/conditions`) returns top‑N condition sets (default 5) with metadata for clients.
- UI: Gradio tab “Recommend Conditions” now surfaces the structured JSON, a table of the top 5 core/base/solvent combos, detection notes, and accepts an optional reaction-type override.

Tuning knobs (relax/constraints)

- `relax`: `{use_drfp, drfp_weight, drfp_n_bits, drfp_radius, precompute_drfp, precompute_scope, reaction_smiles}`.
- `constraints`: inventory/green‑chemistry style filters (see `chemtools/constraints.py`).

Performance notes

- Dataset is cached in‑memory on first use. DRFP precompute is limited to candidates by default.
- RDKit is enabled by default; set `CHEMTOOLS_DISABLE_RDKIT=1` to disable heavy paths. Role‑aware vectors are not attached by default during recommendation.

### Convenience: Single-Molecule Featurization

- POST `/api/v1/featurize/molecule`
  - Default (globals-only): roles omitted → roles=[]
  - In: `{ "smiles": "Clc1ccccc1" }` → returns only global descriptors
  - With types: include role blocks
    - In: `{ "smiles": "Nc1ccccc1", "roles": ["amine"] }`
    - In: `{ "smiles": "Clc1ccccc1", "roles": ["aryl_halide"] }`
- Field schema: `GET /api/v1/featurize/role-aware/fields?roles=amine` or `?roles=` (globals only)

Examples

- Globals-only curl:

  ```bash
  curl -s -X POST http://127.0.0.1:8000/api/v1/featurize/molecule \
    -H "Content-Type: application/json" \
    -d '{"smiles":"Clc1ccccc1"}' | jq '{len: (.vector|length), fields: .fields[0:8], masks: .masks}'
  ```

- Amine-only curl:

  ```bash
  curl -s -X POST http://127.0.0.1:8000/api/v1/featurize/molecule \
    -H "Content-Type: application/json" \
    -d '{"smiles":"Nc1ccccc1","roles":["amine"]}' | jq '{len: (.vector|length), masks: .masks}'
  ```
