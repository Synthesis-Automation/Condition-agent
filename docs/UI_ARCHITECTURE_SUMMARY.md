# UI Architecture Summary

**Date:** October 6, 2025

## Quick Answer

**NO**, `ui_simple.py` does **NOT** use FastAPI endpoints. It is a **standalone Gradio application** that directly imports and calls `chemtools` functions.

## Architecture Overview

Your repository has **two separate application types**:

### 1. **Standalone Gradio UI** (No FastAPI)

**Files:**
- `app/ui_simple.py` (1209 lines) - Simplified recommendation-focused UI
- `app/ui_gradio.py` (2705 lines) - Full feature testing UI

**Architecture:**
```
User Browser
    ↓
Gradio Web Interface (port 7861 or 7860)
    ↓
Direct Python imports
    ↓
chemtools functions (recommend.py, scdb_matcher, router, etc.)
```

**How it works:**
```python
# app/ui_simple.py
import gradio as gr
from chemtools import recommend, router
from chemtools.scdb_matcher import load_db, match

# Direct function calls - NO HTTP requests
def get_ml_recommendations(reaction_smiles, reaction_type, top_k=3):
    result = recommend.recommend_conditions_structured(
        reaction=reaction_smiles,
        reaction_type=reaction_type,
        limit=top_k
    )
    return result

# Gradio interface
demo = gr.Interface(...)
demo.launch(server_port=7861)  # Gradio's own web server
```

**Key characteristics:**
- ✅ Standalone Python process
- ✅ Direct function calls (in-memory)
- ✅ Gradio's built-in web server
- ✅ Runs on port 7861 (`ui_simple.py`) or 7860 (`ui_gradio.py`)
- ❌ No FastAPI involved
- ❌ No HTTP requests between UI and backend

### 2. **FastAPI REST API** (Separate Service)

**File:**
- `app/main.py` (473 lines) - HTTP API endpoints

**Architecture:**
```
HTTP Client (Postman, curl, robot, etc.)
    ↓
FastAPI Server (port 8000)
    ↓
HTTP endpoints (@app.post, @app.get)
    ↓
chemtools functions
```

**Available endpoints (20+ endpoints):**

```python
# app/main.py
from fastapi import FastAPI

app = FastAPI(title="Chemistry Tools API")

# Health check
@app.get("/health")
async def health(): ...

# Recommendation endpoints
@app.post("/api/v1/recommend/from-reaction")
async def api_recommend_from_reaction(req: RecommendFromReactionRequest): ...

@app.post("/api/v1/recommend/conditions")
async def api_recommend_conditions(req: RecommendConditionsRequest): ...

# Rule-based matching
@app.post("/match")
async def match_scheme(req: SchemeMatchRequest): ...

# Utility endpoints
@app.post("/api/v1/smiles/normalize")
async def api_normalize(req: NormalizeRequest): ...

@app.post("/api/v1/router/detect-family")
async def api_detect_family(req: DetectFamilyRequest): ...

@app.post("/api/v1/reaction/detect-type")
async def api_detect_type(req: DetectTypeRequest): ...

# Feature extraction
@app.post("/api/v1/featurize/ullmann")
async def api_featurize_ullmann(req: FeaturizeUllmannRequest): ...

@app.post("/api/v1/featurize/molecular")
async def api_featurize_molecular(req: FeaturizeMolecularRequest): ...

# And 15+ more endpoints...
```

**How to run:**
```bash
# Start FastAPI server
uvicorn app.main:app --reload --port 8000

# Access API docs
http://127.0.0.1:8000/docs
```

**Key characteristics:**
- ✅ HTTP REST API
- ✅ OpenAPI/Swagger documentation at `/docs`
- ✅ Pydantic request/response validation
- ✅ Can be called from any HTTP client
- ✅ Separate from Gradio UIs
- ✅ Runs on port 8000

## Comparison

| Feature | `ui_simple.py` / `ui_gradio.py` | `main.py` (FastAPI) |
|---------|--------------------------------|---------------------|
| **Type** | Gradio Web UI | FastAPI REST API |
| **Port** | 7860/7861 | 8000 |
| **Access** | Browser UI | HTTP endpoints |
| **Backend** | Direct imports | HTTP requests |
| **Documentation** | UI help text | OpenAPI/Swagger |
| **Use Case** | Interactive human use | Programmatic access |
| **Robot Integration** | ❌ Not suitable | ✅ Perfect for robots |
| **Protocol** | WebSocket (Gradio) | HTTP REST |
| **Startup** | `python app/ui_simple.py` | `uvicorn app.main:app` |

## Current Running Processes

Based on terminal history, you tried to run:
```bash
python app/ui_simple.py  # Exit Code: 1 (failed)
```

**This runs:**
- Standalone Gradio UI
- No FastAPI involvement
- Direct chemtools imports

**If you want FastAPI endpoints, you would run:**
```bash
uvicorn app.main:app --reload --port 8000
```

## For Robotic System Integration

### Option 1: Use Existing FastAPI Server ✅ RECOMMENDED

The **FastAPI server (`app/main.py`) is already perfect for robot integration**:

```python
# Robot code (any language with HTTP support)
import requests

# Normalize reaction
response = requests.post(
    "http://localhost:8000/api/v1/smiles/normalize",
    json={"smiles_rxn": "Brc1ccccc1.Nc1ccccc1>>product"}
)
normalized = response.json()

# Get recommendations
response = requests.post(
    "http://localhost:8000/api/v1/recommend/conditions",
    json={
        "reaction": "Brc1ccccc1.Nc1ccccc1>>product",
        "reaction_type": "C_N_Coupling_Pd",
        "limit": 3
    }
)
recommendations = response.json()

# Execute on robot with recommendations
robot.execute(recommendations)
```

**Advantages:**
- ✅ Already implemented
- ✅ 20+ endpoints available
- ✅ OpenAPI documentation
- ✅ Request/response validation
- ✅ Language-agnostic (HTTP)
- ✅ Easy to test (use `/docs` interface)

### Option 2: Direct Python Import (if robot uses Python)

```python
# Robot code (Python only)
from chemtools import recommend

# Direct function call (no HTTP overhead)
recommendations = recommend.recommend_conditions_structured(
    reaction="Brc1ccccc1.Nc1ccccc1>>product",
    reaction_type="C_N_Coupling_Pd",
    limit=3
)

robot.execute(recommendations)
```

**Advantages:**
- ✅ Fastest (no network overhead)
- ✅ Simplest (no server needed)
- ❌ Only works with Python robots
- ❌ Requires chemtools installed on robot

### Option 3: Create MCP Server (for AI-assisted robots)

See `docs/MCP_SERVER_IMPLEMENTATION_GUIDE.md` for details.

## Recommendation for Your Robot

**Use the existing FastAPI server (`app/main.py`)**:

1. **Start the server:**
   ```bash
   uvicorn app.main:app --host 0.0.0.0 --port 8000
   ```

2. **Browse available endpoints:**
   ```
   http://localhost:8000/docs
   ```

3. **Key endpoints for robots:**
   - `POST /api/v1/smiles/normalize` - Validate/canonicalize SMILES
   - `POST /api/v1/reaction/detect-type` - Auto-detect reaction type
   - `POST /api/v1/recommend/conditions` - Get ML recommendations
   - `POST /match` - Get rule-based recommendations
   - `POST /api/v1/featurize/molecular` - Extract molecular features

4. **Robot calls endpoints via HTTP:**
   - Any language (Python, C++, LabVIEW, etc.)
   - Standard HTTP POST requests
   - JSON payloads

## Summary

### Current State:
- ✅ **`ui_simple.py`** - Gradio UI, NO FastAPI, direct imports
- ✅ **`ui_gradio.py`** - Gradio UI, NO FastAPI, direct imports  
- ✅ **`main.py`** - FastAPI REST API, 20+ HTTP endpoints

### For Robotic Integration:
- ✅ **FastAPI server already exists** - Ready to use!
- ✅ **20+ endpoints available** - Comprehensive API
- ✅ **OpenAPI documentation** - Auto-generated at `/docs`
- ✅ **No additional work needed** - Just start the server

### To Get Started:
```bash
# Terminal 1: Start FastAPI server
uvicorn app.main:app --reload --port 8000

# Terminal 2: Test from robot or browser
curl -X POST http://localhost:8000/api/v1/recommend/conditions \
  -H "Content-Type: application/json" \
  -d '{"reaction": "Brc1ccccc1.Nc1ccccc1>>product", "reaction_type": "C_N_Coupling_Pd", "limit": 3}'
```

---

**Bottom line:** Your FastAPI server is already perfect for robot integration. Just start it and point your robot at the HTTP endpoints!
