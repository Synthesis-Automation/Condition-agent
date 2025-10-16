# ChemTools Project (v2)

ChemTools is a deterministic toolkit for reaction condition recommendation, centered on C-N coupling reactions (Ullmann and Buchwald-Hartwig) with hooks for additional reaction families. The project ships as a FastAPI application (`app/main.py`) backed by reusable Python modules in `chemtools/`, command-line helpers, and a pytest suite.

## Feature Highlights

- Reaction parsing and normalization via `chemtools.smiles` with typed request/response contracts defined in `chemtools/contracts.py`.
- Rule-based family detection and router (`chemtools.router`) that drives the downstream recommendation flow.
- **Functional group detection (`chemtools.util.functional_groups`)** with 80+ detectable groups using SMARTS patterns (RDKit) or text-based fallbacks, accessible via unified Context API (`chem.functional_groups.*`).
- Precedent search (`chemtools.precedent`) powered by DRFP similarity, Laplace-smoothed voting, and constraint filters (`chemtools.constraints`).
- Deterministic condition recommendation core (`chemtools.recommend.core`) with optional ML re-ranking (`chemtools.ml`), fusion weighting, and plate design helpers.
- **Protocol recommendation (`chemtools.protocol`)** powered by DRFP similarity search over experimental protocols with automatic condition extraction and standard JSON output.
- Explanation packs (`chemtools.explain`) that summarize precedents, alternatives, and rule hits for API and UI consumption.
- **LLM integration (`llmtools`)** with multi-provider support (OpenAI, Aliyun/DeepSeek) for chemistry-specific agents and natural language processing.
- CLI utilities in `chemtools.rule_scdb_matcher.cli` for querying the Scheme Condition DB (SCDB) and `chemtools.protocol.cli` for protocol index management.

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
|   |-- context.py           # Unified ChemTools API (chem.* namespace)
|   |-- explain.py           # Explanation and rationale builders
|   |-- output_formatter.py  # Standard JSON output formatting for all recommendation modes
|   |-- features/            # Role-aware feature engineering
|   |-- featurizers/         # Molecule and reaction featurizers for C-N coupling
|   |   `-- molecular.py     # C-N coupling substrate featurization (Ullmann, Buchwald, etc.)
|   |-- ml/                  # Optional ML models (DRFP predictor, fusion)
|   |-- precedent/           # Precedent search, similarity, and loaders
|   |-- protocol/            # Protocol recommendation (DRFP-based protocol search)
|   |   |-- indexer.py       # Protocol indexing with DRFP fingerprints
|   |   |-- recommend.py     # DRFP similarity-based protocol recommendation
|   |   `-- cli.py           # CLI for protocol index management
|   |-- recommend/           # Modular recommendation engine and plate design
|   |-- rule_scdb_matcher/   # Scheme Condition DB tooling
|   |-- router.py            # Reaction family detection and dispatch
|   |-- smiles.py            # Reaction SMILES utilities
|   `-- util/                # Shared utilities
|       |-- functional_groups.py  # Comprehensive functional group detection (80+ groups)
|       |-- rdkit_helpers.py      # RDKit integration with graceful fallbacks
|       `-- drfp_storage.py       # DRFP fingerprint caching
|-- data/                    # Sample datasets (JSONL) used in demos/tests
|   |-- protocol_db/         # Experimental protocol JSON files (Structured_Output_schema.json format)
|   |   `-- .protocol_index.json  # Protocol index with precomputed DRFP fingerprints
|   |-- registry_sample.jsonl     # Reagent registry snippets for CLI demo mode
|   `-- reactions_sample.jsonl    # Compact reaction precedent data set for smoke testing
|-- data-processor/          # QT utilities for reagent taxonomy demos
|   `-- requirements.txt     # Data processor specific dependencies
|-- docs/                    # Project documentation (this README lives here)
|   |-- API_DOCUMENTATION.md      # Endpoint details and example payloads
|   |-- PROTOCOL_MODULE.md        # Protocol module technical documentation
|   |-- PROTOCOL_OUTPUT_FORMAT.md # Protocol standard JSON output specification
|   |-- PROTOCOL_QUICKSTART.md    # Protocol module 2-minute quick start
|   |-- PROTOCOL_CLI_GUIDE.md     # Interactive CLI user guide
|   |-- PROTOCOL_README.md        # Protocol module overview
|   `-- requirements.txt          # Documentation build dependencies
|-- llmtools/                # LLM integration for advanced operations
|   |-- clients.py           # Multi-provider LLM client (OpenAI, Aliyun/DeepSeek)
|   |-- agents.py            # Chemistry-specific LLM agents
|   `-- prompts.py           # Prompt templates for chemistry tasks
|-- scripts/                 # Developer helpers and DRFP precompute scripts
|-- tests/                   # Pytest suite covering API and core modules
|   |-- test_protocol_recommendation.py  # Protocol module tests
|   `-- test_protocol_cli.py             # Interactive protocol CLI tester
|-- requirements.txt         # Main production dependencies
|-- requirements-dev.txt     # Development and testing tools (NEW)
|-- FUNCTIONAL_GROUPS_GUIDE.md   # Functional groups detection user guide (NEW)
|-- REQUIREMENTS_GUIDE.md    # Installation and dependency guide (NEW)
|-- run_protocol_cli.ps1     # PowerShell launcher for protocol CLI (Windows)
|-- run_protocol_cli.bat     # Batch launcher for protocol CLI (Windows)
`-- Makefile                 # Task shortcuts (install, run, test, registry, drfp-index)
```

## API Surface

All request/response schemas live in `chemtools/contracts.py`. The FastAPI server exposes the following groups of endpoints:

| Category | Endpoints | Notes |
| --- | --- | --- |
| Health | `GET /health` | Lightweight readiness probe. |
| Reaction utilities | `POST /normalize`, `POST /detect_family`, `POST /detect_type` | Normalize SMILES strings and infer reaction families/types. |
| Featurization | `POST /featurize/molecular`, `POST /featurize/role-aware` | Deterministic role-aware feature generation with RDKit opt-in. |
| Functional groups | Via Context API: `chem.functional_groups.*` | Detect 80+ functional groups in molecules (SMARTS or text patterns). |
| Precedents | `POST /precedent/knn`, `POST /precedent/filters`, `POST /precedent/explain` | DRFP-backed precedent lookup and explanation packs. |
| Recommendation | `POST /api/v1/recommend`, `POST /api/v1/recommend/conditions`, `POST /api/v1/recommend/fusion` | Core recommendation, structured condition sets, and fusion-mode outputs. |
| Plate design | `POST /api/v1/design_plate` | Builds plate-ready experiment grids from recommendation output. |
| Core search | `POST /api/v1/core/search` | Lookup precedents by catalyst core and family. |
| Registry | `POST /match`, `POST /registry` | Scheme Condition DB matcher and reagent registry queries. |

See `docs/API_DOCUMENTATION.md` for detailed request examples and troubleshooting guides.

## CLI and Scripts

- `python -m chemtools.rule_scdb_matcher.cli`: search the Scheme Condition DB, supports CSV/JSONL output and demo mode.
- `python -m chemtools.protocol.cli`: manage protocol index (build, stats, list-families, show-family, show-tag).
- `make registry`: convenience wrapper for the CLI with flags such as `Q`, `FILE`, `JSONL`, `PRETTY`, and `DEMO`.
- `scripts/precompute_drfp.py`: build DRFP NPZ bundles (`make drfp-index`, `make drfp-index-4096`) to warm caches at API startup.
- `scripts/ui_gradio.py`: launch the browser UI for manual testing of the recommendation workflow.
- `test_protocol_cli.py`: interactive CLI for testing protocol recommendations with reaction SMILES input.

## Data and Sample Assets

The project ships with small JSONL samples in `data/` to keep demos deterministic:

- `registry_sample.jsonl`: reagent registry snippets for CLI demo mode.
- `reactions_sample.jsonl`: compact reaction precedent data set for smoke testing.
- `protocol_db/*.json`: experimental protocol files following the Structured_Output_schema.json format (16 protocols currently indexed).
- `protocol_db/.protocol_index.json`: precomputed protocol index with DRFP fingerprints for fast similarity search.

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

## Protocol Recommendation Module

ChemTools now includes a **protocol recommendation system** that finds the most relevant experimental protocol for a given reaction using DRFP similarity search.

### Key Features

- **DRFP-based similarity**: Uses molecular fingerprints to match query reactions to experimental protocols
- **Standard JSON output**: Returns the same format as ML-based and Rule-based recommendations (`meta`, `input`, `detection`, `recommended_conditions`)
- **Fast indexing**: Precomputes DRFP fingerprints for instant search (~100ms per query)
- **Flexible filtering**: Filter protocols by reaction family and tags
- **Condition extraction**: Automatically extracts catalyst, ligand, base, solvent, temperature, time, and atmosphere
- **CLI tools**: Command-line interface for index management and interactive testing

### Quick Start

```bash
# Build the protocol index (one-time setup)
python -m chemtools.protocol.cli build

# View index statistics
python -m chemtools.protocol.cli stats

# Test interactively (Windows)
.\run_protocol_cli.ps1

# Or use Python directly
python test_protocol_cli.py
```

### Python API

```python
from chemtools.protocol import ProtocolRecommender

# Initialize (loads index)
recommender = ProtocolRecommender()

# Get top-3 similar protocols (standard JSON format)
results = recommender.recommend(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=3
)

# Access results using standard format
print(f"Model: {results['meta']['model']}")  # "Protocol-DRFP"
print(f"Detected family: {results['detection']['family']}")
print(f"Confidence: {results['detection']['confidence']:.3f}")

for rec in results['recommended_conditions']:
    protocol = rec['protocol']
    print(f"Rank {rec['rank']}: {protocol['title']}")
    print(f"  Similarity: {rec['similarity']:.3f}")
    print(f"  DOI: {protocol['doi']}")
```

### Output Format

The protocol module returns **standard JSON format** compatible with ML-based and Rule-based recommendations:

```json
{
  "meta": {
    "model": "Protocol-DRFP",
    "status": "success",
    "processing_time_ms": 1702.6
  },
  "input": {
    "reaction_smiles": "...",
    "options": {"k": 3}
  },
  "detection": {
    "family": "Suzuki_Cu_alkyl_halide+aryl_boron",
    "confidence": 0.8018,
    "method": "protocol-similarity"
  },
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 0.8018,
      "protocol": {
        "filename": "Suzuki_Cu_C(sp3)-C(sp2).json",
        "title": "Copper-Catalyzed Suzuki-Miyaura Coupling...",
        "doi": "10.15227/orgsyn.102.0086",
        ...
      },
      "similarity": 0.8018,
      "source": "protocol_database"
    }
  ],
  "extras": {
    "num_candidates": 16,
    "num_total_protocols": 16
  }
}
```

### Documentation

See `docs/PROTOCOL_MODULE.md` for complete documentation, or `docs/PROTOCOL_QUICKSTART.md` for a quick introduction.

## Additional References

- `docs/API_DOCUMENTATION.md`: endpoint details and example payloads.
- `docs/QUICKSTART_API_TEST.md`: step-by-step instructions for exercising the API.
- `docs/TEST_SUZUKI_API_README.md` and related guides: focused testing playbooks for the Suzuki workflow.
- `docs/PROTOCOL_MODULE.md`: complete technical documentation for the protocol recommendation module.
- `docs/PROTOCOL_OUTPUT_FORMAT.md`: detailed specification of the standard JSON output format used by protocol recommendations.
- `docs/PROTOCOL_QUICKSTART.md`: 2-minute quick start guide for the protocol module.
- `docs/PROTOCOL_CLI_GUIDE.md`: user guide for the interactive protocol CLI tester.
- `FUNCTIONAL_GROUPS_GUIDE.md`: comprehensive guide for functional group detection API (80+ groups).
- `REQUIREMENTS_GUIDE.md`: installation guide and dependency troubleshooting.
- `chemtools/contracts.py`: authoritative source for all request and response schemas.
- `chemtools/output_formatter.py`: standard JSON output formatting utilities shared across ML, Rule, and Protocol recommendation modes.

For questions or missing docs, open an issue with reproduction steps and the endpoint/module involved.
