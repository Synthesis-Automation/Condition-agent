# Repository Guidelines

## Project Structure & Module Organization

- `app/main.py`: FastAPI app and route wiring (OpenAPI at `/docs`).
- `chemtools/`: Core deterministic libraries (`smiles.py`, `router.py`, `properties.py`, `precedent.py`, `constraints.py`, `explain.py`).
  - `condition_core.py` and `featurizers/molecular.py` for C-N coupling substrate features.
  - `util/functional_groups.py`: Comprehensive functional group detection (80+ groups) using SMARTS patterns.
  - `protocol/`: Protocol-based recommendation using DRFP similarity.
- `chemtools/recommend/`: Primary HTE-based recommendation interface.
  - `formatters/rule_to_protocol.py`: Converts rule conditions to protocol-compatible format with ordered addition sequences.
  - `cli/registry.py`: CLI entrypoint (`chem-registry`).
- `llmtools/`: LLM integration for advanced operations (`clients.py`, `agents.py`, `prompts.py`).
  - Multi-provider support (OpenAI, Aliyun/DeepSeek).
  - Chemistry-specific agents combining LLM reasoning with chemtools.
  - See `llmtools/README.md` for details.
- `tests/`: Pytest suite (`test_*.py`, fixtures in `conftest.py`).
- `data/`: Small JSONL samples for demos (`registry_sample.jsonl`, `reactions_sample.jsonl`).
  - `protocol_db_v2/`: Protocol database for protocol-based recommendations (validated schema).
  - `rule_db_v2/`: Rule database with standardized types (9 reaction families).
- `scripts/`: Lightweight dev helpers and CLI tools.
- `docs/`: Documentation for features and usage.

## Build, Test, and Development Commands

- Create env and install:
  - macOS/Linux: `python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt`
  - Windows (PowerShell): `python -m venv .venv; .\.venv\Scripts\Activate.ps1; pip install -r requirements.txt`
- Run API (dev): `make run` or `uvicorn app.main:app --reload --port 8000`
- Run tests: `make test` or `pytest -q`
- Registry CLI:
  - Make: `make registry Q="Toluene" PRETTY=1`

## Coding Style & Naming Conventions

- Python ≥ 3.10, PEP 8, 4-space indents; keep modules deterministic (no global state).
- Type hints for all public functions; docstrings for modules and complex functions.
- Naming: `snake_case` (functions/vars), `PascalCase` (classes), `UPPER_SNAKE` (constants).
- Keep API contracts in sync with `chemtools/contracts.py`; prefer simple, explicit data models.

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
   _PATTERN = Chem.MolFromSmarts("[CX4][Cl,Br,I]")  # ❌ Don't do this
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

**Benefits**: 1024-entry global LRU cache provides 10-100× speedup for repeated patterns, shared across all modules.

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
