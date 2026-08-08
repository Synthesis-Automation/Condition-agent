# Reaction Condition Recommender Web UI

Local React/Ketcher client for the repository's canonical structure-first
condition recommender. The browser only owns interaction and presentation; all
reaction analysis, chemistry filtering, retrieval, scoring, discovery, recipe
aggregation, and SVG rendering stay in the standalone Python packages behind a
versioned FastAPI boundary.

Prerequisites are Python 3.10+ with the repository chemistry dependencies and
Node.js 24.14.1+.

## Development

From the repository root, start the API:

```powershell
python -m pip install -r requirements-web.txt
python -m app.web_api
```

In another terminal, start the client:

```powershell
cd web/reaction_recommender
npm install
npm run dev
```

Open `http://127.0.0.1:5173/`. Vite proxies `/api` to the local service on port
8000. The default recommendation index is
`datasets/literature/generic_index.sqlite`; override it with
`CONDITION_RECOMMENDER_INDEX` or `python -m app.web_api --index <path>`.

## Single-port local build

```powershell
cd web/reaction_recommender
npm run build
cd ../..
python -m app.web_api
```

Open `http://127.0.0.1:8000/`. FastAPI serves the compiled client and the
versioned `/api/v1` routes from the same local process. Interactive API
documentation is available at `http://127.0.0.1:8000/api/docs`.

## Supported workflow

- Draw, clear, load, paste, and export reaction SMILES with Ketcher.
- Validate product-fragment source requirements before recommendation.
- Retrieve and rank chemically compatible canonical condition recipes.
- Apply declarative ranking profiles or transparent custom weights.
- Browse structural precedents through the four discovery views.
- Analyze a molecule or reaction with the same deterministic featurization used
  by the Qt tool, including motifs, reactive sites, reaction-core evidence,
  partner roles, mapping provenance, and the canonical nested analysis.
- Inspect score traces, structural matches and mismatches, cautions, conditions,
  yields, fallback levels, and precedent provenance.
- Export the complete versioned recommendation, discovery, or feature result as
  JSON.
- Preview and download the selected recommendation's versioned synthesis
  protocol JSON, including registry substance IDs, CAS numbers, quantities,
  operating conditions, observed operations, and execution-readiness gaps.

The UI intentionally exposes no arbitrary file paths or upload endpoints. Local
dataset identity and access remain server configuration concerns.
