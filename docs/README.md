# Condition Agent Documentation

This folder collects product docs, API notes, and developer references. This README provides a short index and a current map of the codebase layout so new contributors can orient quickly.

## Quick links

- API entry: `app/main.py` (FastAPI, OpenAPI at `/docs`)
- Core deterministic library: `chemtools/`
- LLM integration: `llmtools/`
- Agent wrapper and UX helpers: `chem_assistant/`
- HTE recommender docs: `docs/HTE_RECOMMENDER.md`
- API usage notes: `docs/HTE_RECOMMENDER_API.md`

## Repository layout (current)

```
app/
  main.py                         FastAPI wiring
  HTE_recommender_gui.py          GUI for condition recommendations
  HTE_analytics_cli.py            Analytics CLI
  A_convert_to_hte_format*.py     HTE format converters
  reagent_*_gui.py                Reagent analysis/ordering UI
chemtools/
  smiles.py                       SMILES helpers
  router.py                       Reaction routing
  context.py                      Reaction context assembly
  reaction_similarity.py          Similarity utilities
  reaction_type_detection.py      Reaction family detection
  precedent/                      Precedent search and similarity
  constraints.py                  Rule constraints
  explain.py                      Explanation helpers
  featurizers/                    Molecular and reaction featurizers
  recommend/                      Primary HTE recommendation interface
  formatters/                     Output formatting and normalization
  protocol/                       Protocol DB indexing/recommender
  reagent/                        Reagent registry + taxonomy tools
  taxonomy/                       Reaction taxonomy + specs
  util/                           SMARTS cache, RDKit helpers, etc.
chem_assistant/
  chemtools_agent.py              LLM-based agent entry
  chemtools_wrapper.py            Tooling surface for agents
  rag.py                          Lightweight RAG utilities
  gui/                            Local UI helpers
llmtools/
  clients.py                      Multi-provider LLM clients
  agents.py                       LLM agent utilities
  prompts.py                      Prompt templates
  recommendation_llm.py           LLM-assisted recommendation
  reagent_classifier.py           Reagent classification LLM
  README.md                       LLM tools overview
examples/
  sample_reactions.csv            Example reaction inputs
  rule_db_v2/                     Example rules
  sample500/                      Small sample dataset
data-processor/
  process_reactions.py            Dataset normalization pipeline
  v2_processor_core.py            Canonicalization core
  Scifinder_rdf_processer.py      RDF ingestion for source data
scripts/
  index_knowledge_base.py         Indexing utilities for RAG
  validate_taxonomy_chemistry.py  Taxonomy validation
  analyze_taxonomy.py             Taxonomy analysis
  rag_eval.py                     RAG evaluation
  eval_spectator_groups.py        Spectator group scoring
  check_kb_conditions.py          KB sanity checks
data/
  protocol_db_v2/                 Protocol database (schema-validated)
  rule_db_v2/                     Rule database (9 reaction families)
  registry_sample.jsonl           Example registry data
  reactions_sample.jsonl          Example reactions
results/
  hte_cache/                      Cached HTE outputs
  visualization/                  Plots and summaries
  test_out.csv                    Sample output

tests/
  test_*.py                       Pytest suite
  conftest.py                     Fixtures

app/services/                     Service helpers for the API
chemtools/visualization/          Plotting helpers
chemtools/util/functional_groups.py  80+ functional group SMARTS
```

## Codebase structure and responsibilities

### 1) Deterministic chemistry toolkit (`chemtools/`)

- **Core routing and context**: `router.py`, `context.py`, `reaction_type_detection.py`
- **Precedent search**: `precedent/` loaders, search, and similarity utilities
- **Featurization**: `featurizers/` and `featurizers/analysis/` for C-N coupling and reaction pairs
- **Recommendation**: `recommend/` provides the HTE-based recommender and adapters
- **Protocols**: `protocol/` handles DRFP indexing, matching, and protocol recommendation
- **Reagents**: `reagent/` provides registry, taxonomy, validation, and lookup utilities
- **Utilities**: `util/` for SMARTS caching, RDKit guards, DRFP storage, and helpers

### 2) API and user-facing entry points (`app/`)

- **FastAPI server**: `app/main.py` exposes the REST API and OpenAPI docs
- **GUIs and CLIs**: `HTE_recommender_gui.py`, `HTE_analytics_cli.py`, and format conversion tools

### 3) LLM integrations (`llmtools/` + `chem_assistant/`)

- **Provider clients**: `llmtools/clients.py` (OpenAI + Aliyun/DeepSeek)
- **Agent logic**: `llmtools/agents.py` and `chem_assistant/chemtools_agent.py`
- **Tool surface**: `chem_assistant/chemtools_wrapper.py` exposes deterministic tools to LLMs
- **RAG**: `chem_assistant/rag.py` + `scripts/index_knowledge_base.py` for indexing

### 4) Data, results, and examples

- **Data**: `data/` stores small JSONL and structured protocol/rule DBs
- **Examples**: `examples/` contains sample inputs and rule DB snippets
- **Results**: `results/` holds cache and visualizations (safe to delete caches)

## Build and test quickstart

```
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt

pytest -q

uvicorn app.main:app --reload --port 8000
```

## SMARTS compilation rule (important)

Always use the centralized cache:

```
from chemtools.util.smarts_cache import compile_smarts

pattern = compile_smarts("[CX4][Cl,Br,I]", validate=False)
```

Never call `Chem.MolFromSmarts()` directly.

## Documentation index

- `docs/HTE_RECOMMENDER.md` - System behavior and scoring
- `docs/HTE_RECOMMENDER_API.md` - API payloads and usage
- `docs/HTE_ANALYTICS.md` - Analytics and evaluation
- `docs/REACTION_DETECTION_METHODS.md` - Detection heuristics
- `docs/UNIFIED_FEATURIZERS.md` - Feature alignment
- `docs/TAXONOMY_*` - Taxonomy strategy and alignment
- `docs/WEB_UI_IMPLEMENTATION_PLAN.md` - UI plan

## Notes

- RDKit is optional at runtime; see `chemtools/util/rdkit_helpers.py` for graceful degradation.
- Keep API contracts in sync with `chemtools/contracts.py` when changing request/response formats.
- The authoritative contribution rules live in `AGENTS.md` at repo root.
