# HTE Web Interface Implementation Plan (FastAPI + JS)

## Goals

- Provide a cloud-friendly web interface for the HTE recommender.
- Keep chemistry logic in the backend; frontend only renders results.
- Match current GUI inputs/outputs (reaction SMILES, filters, top_k, min_exp, source group, weighting).

## Scope

- Backend: FastAPI endpoints for HTE recommendation + stats.
- Frontend: lightweight static HTML/JS/CSS UI.
- Deployment: Linux server with uvicorn/gunicorn, optional Nginx reverse proxy.

## Assumptions

- The HTE database is available on the server filesystem (e.g., `data/HTE_db`).
- RDKit is optional; degrade gracefully if visualization unavailable.
- Users will not be allowed to browse arbitrary paths on the server.

## Deliverables

- `app/main.py` routes for HTE recommendations and stats.
- Request/response schemas (Pydantic) aligned to GUI fields.
- Static frontend under `app/static/` (or `docs/` if preferred).
- Basic README updates for running locally and on Linux.

## API Design

### POST `/hte/recommend`

Request JSON:

- `reaction_smiles` (string, required)
- `db_path` (string, optional; default: `data/HTE_db`)
- `top_k` (int, default 10)
- `min_experiments` (int, default 1)
- `reaction_type_filter` (string, optional)
- `catalyst_filter` (string, optional)
- `source_group` (string, optional: literature|experiments|rules)
- `use_aryl_weighting` (bool, default false)

Response JSON:

- `summary` (reactant types, predicted reaction type, coverage, matches)
- `recommendations` (list of condition recommendations)
- `recommendations_by_source` (optional; keyed by source group)
- `stats` (HTE database summary)
- `warnings` (optional list of strings)

### GET `/hte/stats`

Returns HTE database summary stats.

## Frontend UI

- Form inputs mirroring the GUI:
  - Reaction SMILES
  - Top K
  - Min experiments
  - Reaction filter
  - Catalyst filter
  - Source group
  - Aryl weighting
- Result view:
  - Reaction image panel (rendered from backend or client-side if available)
  - Summary block
  - Results table with uniform columns
  - Tabs/filters by source group (optional)

## Implementation Steps

1) Backend wiring
   - Add new FastAPI routes in `app/main.py` (or a dedicated router module).
   - Create Pydantic schemas for request/response payloads.
   - Reuse `chemtools.recommend.HTERecommender` for recommendation logic.
   - Normalize results to JSON (consistent with GUI fields).

2) Reaction image support
   - Option A: Return image as base64 PNG in the API response.
   - Option B: Create a dedicated `/hte/render` endpoint that returns a PNG.
   - Graceful fallback if visualization libraries are missing.

3) Frontend
   - Build a static HTML page with a simple form and a results table.
   - Use vanilla JS `fetch` to call `/hte/recommend`.
   - Display reaction image + summary + results table.

4) Security/Validation
   - Restrict `db_path` to allowed directories (e.g., `data/HTE_db`).
   - Validate `reaction_smiles` is non-empty.
   - Set max limits on `top_k` / `min_experiments`.

5) Deployment
   - Add instructions for Linux deployment (uvicorn + gunicorn).
   - Optional Nginx reverse proxy and static file hosting.

## Testing

- Unit test for API response structure (pytest + FastAPI TestClient).
- Smoke test: request with known example, ensure non-empty recommendations.

## Open Questions

- Do you want to allow user-supplied DB paths or lock to a single server path?
   No
- Should the reaction image be server-rendered or skipped if RDKit is missing?
- Do you need auth (basic token) for the endpoint?
