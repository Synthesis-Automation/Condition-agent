# ChemTools Project (v2)

ChemTools is a deterministic toolkit for reaction condition recommendation, centered on Ullmann C–N couplings with hooks for additional reaction families. The project ships as a FastAPI application (`app/main.py`) backed by reusable Python modules in `chemtools/`, command-line helpers, and a pytest suite.

## Feature Highlights

- Reaction parsing and normalization via `chemtools.smiles` with typed request/response contracts defined in `chemtools/contracts.py`.
- Rule-based family detection and router (`chemtools.router`) that drives the downstream recommendation flow.
- Precedent search (`chemtools.precedent`) powered by DRFP similarity, Laplace-smoothed voting, and constraint filters (`chemtools.constraints`).
- Deterministic condition recommendation core (`chemtools.recommend.core`) with optional ML re-ranking (`chemtools.ml`), fusion weighting, and plate design helpers.
- Explanation packs (`chemtools.explain`) that summarize precedents, alternatives, and rule hits for API and UI consumption.
- CLI utilities in `chemtools.rule_scdb_matcher.cli` for querying the Scheme Condition DB (SCDB).

## Quickstart

### Environment setup

```bash
# macOS / Linux
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Windows (PowerShell)
python -m venv .venv
.\\.venv\\Scripts\\Activate.ps1
pip install -r requirements.txt
```

Install optional dev tooling (linters, pytest plugins) with:

```bash
pip install -r requirements-dev.txt
```

### Run the API

```bash
uvicorn app.main:app --reload --port 8000
# Open http://127.0.0.1:8000/docs for the OpenAPI UI
```

Make targets wrap the common tasks: `make install`, `make run`, `make test`, `make registry`.

### Run the CLI against the Scheme Condition DB

```bash
python -m chemtools.rule_scdb_matcher.cli "Toluene"
python -m chemtools.rule_scdb_matcher.cli --jsonl -f queries.txt
```

Pass `--demo` to load sample data from `data/registry_sample.jsonl`.

### Run the Gradio demo UI

```bash
python scripts/ui_gradio.py
```

### Execute the test suite

```bash
pytest -q
```

The tests rely on lightweight fixtures stored under `tests/` and do not require RDKit unless explicitly enabled.

## Repository Layout

```text
.
|-- app/                     # FastAPI application wiring
|   |-- main.py              # REST endpoints and request orchestration
|   `-- ui_gradio.py         # Optional Gradio front-end
|-- chemtools/               # Core deterministic libraries
|   |-- contracts.py         # Typed request/response dataclasses
|   |-- condition_core.py    # Legacy core selection logic (compat shim)
|   |-- constraints.py       # Constraint rule definitions
|   |-- context.py           # Shared runtime configuration helpers
|   |-- explain.py           # Explanation and rationale builders
|   |-- features/            # Role-aware feature engineering
|   |-- featurizers/         # Molecule and reaction featurizers (incl. Ullmann)
|   |-- ml/                  # Optional ML models (DRFP predictor, fusion)
|   |-- precedent/           # Precedent search, similarity, and loaders
|   |-- recommend/           # Modular recommendation engine and plate design
|   |-- rule_scdb_matcher/   # Scheme Condition DB tooling
|   |-- router.py            # Reaction family detection and dispatch
|   `-- smiles.py            # Reaction SMILES utilities
|-- data/                    # Sample datasets (JSONL) used in demos/tests
|-- data-processor/          # QT utilities for reagent taxonomy demos
|-- docs/                    # Project documentation (this README lives here)
|-- scripts/                 # Developer helpers and DRFP precompute scripts
|-- tests/                   # Pytest suite covering API and core modules
`-- Makefile                 # Task shortcuts (install, run, test, registry, drfp-index)
```

## API Surface

All request/response schemas live in `chemtools/contracts.py`. The FastAPI server exposes the following groups of endpoints:

| Category | Endpoints | Notes |
| --- | --- | --- |
| Health | `GET /health` | Lightweight readiness probe. |
| Reaction utilities | `POST /normalize`, `POST /detect_family`, `POST /detect_type` | Normalize SMILES strings and infer reaction families/types. |
| Featurization | `POST /featurize/ullmann`, `POST /featurize/role-aware` | Deterministic role-aware feature generation with RDKit opt-in. |
| Precedents | `POST /precedent/knn`, `POST /precedent/filters`, `POST /precedent/explain` | DRFP-backed precedent lookup and explanation packs. |
| Recommendation | `POST /api/v1/recommend`, `POST /api/v1/recommend/conditions`, `POST /api/v1/recommend/fusion` | Core recommendation, structured condition sets, and fusion-mode outputs. |
| Plate design | `POST /api/v1/design_plate` | Builds plate-ready experiment grids from recommendation output. |
| Core search | `POST /api/v1/core/search` | Lookup precedents by catalyst core and family. |
| Registry | `POST /match`, `POST /registry` | Scheme Condition DB matcher and reagent registry queries. |

See `docs/API_DOCUMENTATION.md` for detailed request examples and troubleshooting guides.

## CLI and Scripts

- `python -m chemtools.rule_scdb_matcher.cli`: search the Scheme Condition DB, supports CSV/JSONL output and demo mode.
- `make registry`: convenience wrapper for the CLI with flags such as `Q`, `FILE`, `JSONL`, `PRETTY`, and `DEMO`.
- `scripts/precompute_drfp.py`: build DRFP NPZ bundles (`make drfp-index`, `make drfp-index-4096`) to warm caches at API startup.
- `scripts/ui_gradio.py`: launch the browser UI for manual testing of the recommendation workflow.

## Data and Sample Assets

The project ships with small JSONL samples in `data/` to keep demos deterministic:

- `registry_sample.jsonl`: reagent registry snippets for CLI demo mode.
- `reactions_sample.jsonl`: compact reaction precedent data set for smoke testing.

Larger proprietary datasets should not be committed. Configure their paths via environment variables or local overrides.

## Testing and QA

- Run `pytest -q` locally before submitting changes. The suite covers API contracts (`tests/test_*_api.py`) and domain logic (`tests/test_*.py`).
- Add new tests near the affected modules and use fixtures defined in `tests/conftest.py`.
- For API changes, regenerate or update examples in `docs/API_DOCUMENTATION.md` and the Suzuki/Fusion test guides under `docs/`.

## Environment Toggles

ChemTools can run without RDKit for quick demos. Key environment variables:

- `CHEMTOOLS_DISABLE_RDKIT=1` to skip RDKit-heavy paths.
- `CHEMTOOLS_LOAD_DATASET=0` to avoid loading large precedent bundles.
- `CHEMTOOLS_ATTACH_ROLE_AWARE=1` to include role-aware vectors during dataset loading.
- `CHEMTOOLS_DRFPPATH=artifacts/drfp_index.npz` to preload a DRFP cache generated by `make drfp-index`.

## Additional References

- `docs/API_DOCUMENTATION.md`: endpoint details and example payloads.
- `docs/QUICKSTART_API_TEST.md`: step-by-step instructions for exercising the API.
- `docs/TEST_SUZUKI_API_README.md` and related guides: focused testing playbooks for the Suzuki workflow.
- `chemtools/contracts.py`: authoritative source for all request and response schemas.

For questions or missing docs, open an issue with reproduction steps and the endpoint/module involved.
