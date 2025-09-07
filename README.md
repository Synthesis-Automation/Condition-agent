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
