# Repository Guidelines

## System Principles

- **Taxonomy-driven**: Reaction typing, reactant classification, and recommendation logic must flow from `chemtools/taxonomy` and chemistry rules, not ad-hoc heuristics.
- **Chemistry-first**: Validate transformations against known reaction families and functional group constraints before scoring or ranking.
- **Keep it clean**: Remove legacy code paths and deprecated modules when refactoring; backward compatibility is not a goal unless explicitly requested.

## Project Structure & Module Organization

- `app/`: FastAPI app (`app/main.py`) plus GUI/CLI entry scripts and service wiring.
- `chemtools/`: Core deterministic libraries (routing, detection, constraints, explainers).
  - `taxonomy/`: Reaction taxonomy data, loaders, validation, and catalog logic.
  - `recommend/`: HTE-based recommendation engine and CLI.
  - `formatters/`: Output normalization and protocol/rule formatting.
  - `precedent/`, `protocol/`, `reagent/`, `featurizers/`, `util/`, `visualization/`.
  - `util/functional_groups.py`: Comprehensive functional group detection using SMARTS patterns.
- `llmtools/`: LLM integration for advanced operations (`clients.py`, `agents.py`, `prompts.py`).
  - Multi-provider support (OpenAI, Aliyun/DeepSeek).
  - Chemistry-specific agents combining LLM reasoning with chemtools.
  - See `llmtools/README.md` for details.
- `data-processor/`: ETL and data preparation utilities.
- `data/`: Curated datasets (HTE, protocols, reaction datasets, reagent DB, knowledge base).
  - `HTE_db/`, `protocol_db_v2/`, `reaction_dataset/`, `reagent_db/`, `knowledge_base/`.
- `tests/`: Pytest suite (`test_*.py`, fixtures in `conftest.py`).
- `scripts/`: Lightweight dev helpers and CLI tools.
- `docs/`: Documentation for features and usage.
- `results/`: Local analysis outputs (do not commit large artifacts).

## Build, Test, and Development Commands

- Create env and install:
  - macOS/Linux: `python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt`
  - Windows (PowerShell): `python -m venv .venv; .\.venv\Scripts\Activate.ps1; pip install -r requirements.txt`
- Run API (dev): `uvicorn app.main:app --reload --port 8000`
- Run tests: `pytest -q`
- HTE CLI:
  - `python -m chemtools.recommend.cli --help`

## Coding Style & Naming Conventions

- Python >= 3.10, PEP 8, 4-space indents; keep modules deterministic (no global state).
- Type hints for all public functions; docstrings for modules and complex functions.
- Naming: `snake_case` (functions/vars), `PascalCase` (classes), `UPPER_SNAKE` (constants).
- Keep API contracts in sync with `chemtools/contracts.py`; prefer simple, explicit data models.
- Align any new reaction rules, labels, or schemas with `chemtools/taxonomy` and its validation tools.

### SMARTS Pattern Compilation

**Always use the centralized cache** for SMARTS pattern compilation:

```python
from chemtools.util.smarts_cache import compile_smarts

# Compile a single pattern (cached automatically)
pattern = compile_smarts("[CX4][Cl,Br,I]", validate=False)
if pattern and mol.HasSubstructMatch(pattern):
    # Match found
    pass
```

**Best Practices:**

1. **Module-level pattern definitions**: Define SMARTS strings as module-level constants:

   ```python
   _MY_PATTERNS = {
       "aryl_halide": "[$(c[Cl,Br,I]),$(c-[Cl,Br,I])]",
       "alcohol": "[OX2H]",
   }
   ```

2. **Lazy compilation**: Compile patterns only when needed, not at module import:

   ```python
   # Good: Lazy compilation via centralized cache
   def detect_feature(mol):
       pattern = compile_smarts(_MY_PATTERNS["aryl_halide"], validate=False)
       return mol.HasSubstructMatch(pattern) if pattern else False
   
   # Bad: Eager compilation at module load
   _PATTERN = Chem.MolFromSmarts("[CX4][Cl,Br,I]")  # Don't do this
   ```

3. **Validation**: Use `validate=True` for critical patterns, `validate=False` for speed:

   ```python
   # Development: Validate new patterns
   pattern = compile_smarts(smarts, validate=True)
   
   # Production: Skip validation for known-good patterns
   pattern = compile_smarts(smarts, validate=False)
   ```

4. **Never call `Chem.MolFromSmarts()` directly** - always use `compile_smarts()` to benefit from global caching.

5. **Batch compilation**: For multiple patterns at startup:

   ```python
   from chemtools.util.smarts_cache import compile_smarts_batch
   
   patterns = compile_smarts_batch(_MY_PATTERNS, skip_invalid=True)
   ```

**Benefits**: 1024-entry global LRU cache provides 10-100x speedup for repeated patterns, shared across all modules.

## Testing Guidelines

- Framework: `pytest`; tests live in `tests/test_*.py` mirroring module names.
- Add unit tests with meaningful edge cases; use fixtures from `tests/conftest.py`.
- Keep tests fast and deterministic; aim to maintain or improve coverage for changed code.

## Commit & Pull Request Guidelines

- Conventional Commits: `feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`; scopes like `api`, `cli`, `chemtools` are encouraged.
- Subject in imperative mood (<=72 chars) with optional body for context.
- PRs must include: summary, motivation, linked issues, testing notes, and API impact (screenshots of `/docs` if routes change).
- Ensure `pytest -q` passes and `uvicorn app.main:app --reload --port 8000` boots locally before request for review.

## Security & Configuration Tips

- Do not commit secrets; prefer environment variables if configuration is introduced.
- RDKit is required whenever molecular operations are needed; do not implement optional fallbacks unless explicitly requested.
- Sample data in `data/` is for demos only - avoid adding large/proprietary datasets.
