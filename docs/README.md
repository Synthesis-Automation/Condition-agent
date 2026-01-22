# ChemTools Toolkit

ChemTools is a deterministic toolkit for reaction-condition design and analysis, focused on (but not limited to) C–N coupling chemistry. It powers API services, CLI tools, and automation workflows used to normalize reactions, detect functional motifs, search precedent data, and surface recommended conditions.

---

## Key Capabilities

- **Reaction normalization** – canonicalize SMILES/reaction SMILES with graceful fallbacks when RDKit is absent.
- **Molecular feature extraction** – fingerprint/descriptor generation for electrophile–nucleophile pairs.
- **Unified recommendations** – DRFP + feature-tag similarity over reaction datasets, protocols, and HTE screens.
- **Protocol & taxonomy utilities** – curated reagent registries, taxonomy schemas, and protocol discovery/SMARTS tooling for experimental planning.
- **LLM extensions** – chemistry-aware agent orchestration with multi-provider support (OpenAI, Aliyun/DeepSeek) via `llmtools`.

---

## Getting Started

```bash
# create & activate virtual environment
python -m venv .venv
. .venv/Scripts/activate    # Windows (PowerShell): .\.venv\Scripts\Activate.ps1

# install core dependencies
pip install -r requirements.txt

# optional dev tooling
pip install -r requirements-dev.txt

# run tests
pytest -q
```

## Repository Layout

```text
.
├─ app/                      # FastAPI application wiring & Gradio demo UI
├─ chemtools/                # Core deterministic libraries
│  ├─ analysis/              # SMILES & reaction normalization helpers
│  ├─ featurizers/           # Molecular featurization
│  ├─ precedent/             # DRFP-based precedent search
│  ├─ protocol/              # Protocol recommendation & SMARTS tools
│  ├─ recommend/             # Condition recommendation engines
│  └─ util/                  # Functional group detection, RDKit helpers, etc.
├─ chem_assistant/           # Assistant-facing wrappers and CLI helpers
├─ data/                     # Sample datasets (registry, reactions, HTE)
├─ docs/                     # Project documentation (this README)
├─ llmtools/                 # LLM client integrations and prompt assets
├─ scripts/                  # Developer utilities & showcases
├─ tests/                    # Pytest suite with fixtures
└─ requirements*.txt         # Dependency manifests
```

---

## Useful Commands

```bash
# run API locally
uvicorn app.main:app --reload --port 8000

```

ChemTools is fully deterministic by default, with optional higher-level integrations (ML re-ranking, LLM agents) activated as your workflow requires.
