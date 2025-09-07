# chemtools-project (v2)

Minimal scaffold for deterministic **chemistry tools** (FastAPI + RDKit-friendly) to power condition recommendation.

## Quickstart
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
# open http://127.0.0.1:8000/docs
```
