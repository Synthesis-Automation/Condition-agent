# Repository Guidelines

## Project Structure & Module Organization
- `app/main.py`: FastAPI app and route wiring (OpenAPI at `/docs`).
- `chemtools/`: Core deterministic libraries (`smiles.py`, `router.py`, `properties.py`, `precedent.py`, `constraints.py`, `explain.py`).
  - `condition_core.py` and `featurizers/ullmann.py` for Ullmann C–N.
  - `cli/registry.py`: CLI entrypoint (`chem-registry`).
- `tests/`: Pytest suite (`test_*.py`, fixtures in `conftest.py`).
- `data/`: Small JSONL samples for demos (`registry_sample.jsonl`, `reactions_sample.jsonl`).
- `scripts/`: Lightweight dev helpers.

## Build, Test, and Development Commands
- Create env and install:
  - macOS/Linux: `python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt`
  - Windows (PowerShell): `python -m venv .venv; .\.venv\Scripts\Activate.ps1; pip install -r requirements.txt`
- Run API (dev): `make run` or `uvicorn app.main:app --reload --port 8000`
- Run tests: `make test` or `pytest -q`
- Registry CLI:
  - Make: `make registry Q="Toluene" PRETTY=1`
  - Module: `python -m chemtools.cli.registry --jsonl -f queries.txt`

## Coding Style & Naming Conventions
- Python ≥ 3.10, PEP 8, 4-space indents; keep modules deterministic (no global state).
- Type hints for all public functions; docstrings for modules and complex functions.
- Naming: `snake_case` (functions/vars), `PascalCase` (classes), `UPPER_SNAKE` (constants).
- Keep API contracts in sync with `chemtools/contracts.py`; prefer simple, explicit data models.

## Testing Guidelines
- Framework: `pytest`; tests live in `tests/test_*.py` mirroring module names.
- Add unit tests with meaningful edge cases; use fixtures from `tests/conftest.py`.
- Keep tests fast and deterministic; aim to maintain or improve coverage for changed code.

## Commit & Pull Request Guidelines
- Conventional Commits: `feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`; scopes like `api`, `cli`, `chemtools` are encouraged.
- Subject in imperative mood (≤72 chars) with optional body for context.
- PRs must include: summary, motivation, linked issues, testing notes, and API impact (screenshots of `/docs` if routes change).
- Ensure `pytest -q` passes and `make run` boots locally before request for review.

## Security & Configuration Tips
- Do not commit secrets; prefer environment variables if configuration is introduced.
- RDKit is optional at runtime; ensure graceful degradation where used (see `util/rdkit_helpers.py`).
- Sample data in `data/` is for demos only—avoid adding large/proprietary datasets.

