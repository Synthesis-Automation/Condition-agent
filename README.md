# chemtools-project (v2)

Deterministic chemistry tools (FastAPI + RDKit-friendly) to support condition recommendation. Focuses on Ullmann C–N first, with simple hooks for more families.

## Quickstart
- Create a virtualenv and install deps:
  - macOS/Linux:
    ```bash
    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    ```
  - Windows (PowerShell):
    ```powershell
    python -m venv .venv
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
  pytest -q
  ```
  Or with Makefile shortcuts: `make install`, `make run`, `make test`.

## How to Test

### 1) Unit tests (fast, deterministic)
- Run all: `pytest -q`
- Notes: tests disable heavy RDKit paths automatically; no extra setup needed.

### 2) Exercise the API via Swagger
- Start: `uvicorn app.main:app --reload --port 8000`
- Open: http://127.0.0.1:8000/docs
- Try endpoints:
  - POST `/api/v1/smiles/normalize` with `{ "smiles": "c1ccccc1O" }`
  - POST `/api/v1/router/detect-family` with `{ "reactants": ["Clc1ccccc1","Nc1ccccc1"] }`
  - POST `/api/v1/featurize/ullmann` with `{ "electrophile": "Clc1ccccc1", "nucleophile": "Nc1ccccc1" }`

### 3) Quick Python one-liners
- Normalize: `python -c "from chemtools.smiles import normalize; print(normalize('c1ccccc1O')['smiles_norm'])"`
- Detect family: `python -c "from chemtools.router import detect_family; print(detect_family(['Clc1ccccc1','Nc1ccccc1']))"`
- Featurize: `python -c "from chemtools.featurizers.ullmann import featurize; print(featurize('Clc1ccccc1','Nc1ccccc1'))"`

### 4) Gradio UI (no-code testing)
- Launch: `python app/ui_gradio.py`
- Open: http://127.0.0.1:7860
- Tabs included:
  - SMILES Normalize: normalize a single SMILES.
  - Detect Family: infer reaction family from reactants (dot or newline separated).
  - Ullmann Featurizer: inspect LG/nucleophile class/bin features.
  - Properties Lookup: query by name/CAS/token (e.g., "Water", "7778-53-2").
  - Recommend Conditions: runs the recommender over a reaction SMILES.
  - Design Plate: proposes a plate CSV across top cores.
  - Precedent Search: retrieves similar precedents (optionally DRFP re-ranking).
  - DRFP Similarity: Tanimoto between two reaction SMILES (if DRFP installed).

Optional: install DRFP to enable similarity re-ranking and the similarity tab
```
pip install drfp==0.4.0 numpy
```

Tips
- Speed up or run without RDKit: set `CHEMTOOLS_DISABLE_RDKIT=1`
  - macOS/Linux: `export CHEMTOOLS_DISABLE_RDKIT=1`
  - Windows (PowerShell): `$env:CHEMTOOLS_DISABLE_RDKIT='1'`
- If running `python app/ui_gradio.py` outside the repo root, prefer module mode: `python -m app.ui_gradio`

## Endpoints
- `GET /health`: health probe.
- `POST /api/v1/smiles/normalize` — Normalize SMILES.
  - In: `{ "smiles": "Brc1ccc(F)cc1.O" }`
  - Out: `{ "input", "fragments", "largest_smiles", "smiles_norm" }`
- `POST /api/v1/router/detect-family` — Detect reaction family.
  - In: `{ "reactants": ["Brc1ccc(F)cc1", "Nc1ccccc1"] }`
  - Out: `{ "family", "confidence", "hits" }`
- `POST /api/v1/featurize/ullmann` — Substrate features for Ullmann C–N.
  - In: `{ "electrophile": "...", "nucleophile": "..." }`
  - Out: feature dict including `LG`, `nuc_class`, `bin`, …
- `POST /api/v1/condition-core/parse` — Parse ConditionCore from reagents/text.
  - In: `{ "reagents": [{"uid":"7681-65-4","role":"CATALYST","name":"CuI"}, {"uid":"72-52-8","role":"LIGAND","name":"Phenanthroline"}], "text": "optional free text" }`
  - Out: `{ "core": "Cu/Phenanthroline", ... }`
- `POST /api/v1/properties/lookup` — Minimal properties by UID/token.
  - In: `{ "query": "K3PO4" }` or `{ "query": "7778-53-2" }`
- `POST /api/v1/precedent/knn` — Retrieve similar precedents.
  - In: `{ "family": "Ullmann_CN", "features": {"bin":"LG:Br|NUC:aniline"}, "k": 50 }`
- `POST /api/v1/constraints/filter` — Inventory/blacklist filtering of candidate IDs.
- `POST /api/v1/explain/precedents` — Short reasons and example precedents.

Notes:
- RDKit is optional at runtime; if unavailable, SMILES normalization falls back to simple heuristics.
- Sample data lives under `data/` and is used by the precedent demo.

## Project Structure
- `app/main.py` — FastAPI app and routes.
- `chemtools/` — Deterministic tool implementations:
  - `smiles.py`, `router.py`, `properties.py`, `precedent.py`, `constraints.py`, `explain.py`
  - `condition_core.py` — normalizes ConditionCore
  - `featurizers/ullmann.py` — Ullmann C–N featurizer
  - `util/rdkit_helpers.py` — RDKit helpers (safe import)
- `tests/` — Lightweight unit tests
- `data/` — Small JSONL samples for precedents
- `chemistry_tool_bulild_plan.md` — Build plan and roadmap

## Development Notes
- API surfaces match the contracts in `chemtools/contracts.py`.
- Performance-sensitive code paths are simple and deterministic; add caching where beneficial.
- See the build plan for roadmap and design details: `chemistry_tool_bulild_plan.md`.

## License
Internal project scaffold. Add a LICENSE if distributing.
