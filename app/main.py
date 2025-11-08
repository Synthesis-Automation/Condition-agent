from fastapi import FastAPI, Request, Response
import os
import logging
import time

from chemtools.contracts import RecommendConditionsRequest

# Enable RDKit by default unless explicitly disabled
os.environ.setdefault("CHEMTOOLS_DISABLE_RDKIT", "0")

# Import service layer
from app.services import recommendation_service

# Import error handlers
from app.services.error_handlers import register_error_handlers

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("chemtools.api")

app = FastAPI(title="Chemistry Tools API", version="0.1.2")

# Register error handlers for consistent error responses
register_error_handlers(app)


@app.middleware("http")
async def logging_timing_middleware(request: Request, call_next):
    start = time.perf_counter()
    path = request.url.path
    method = request.method
    status = None
    try:
        response: Response = await call_next(request)
        status = response.status_code
        return response
    finally:
        dur = time.perf_counter() - start
        # Log without body to avoid PII
        status_str = str(status) if status is not None else "ERROR"
        logger.info("%s %s -> %s in %.3f s", method, path, status_str, dur)

@app.get("/health")
def health():
    return {"ok": True}


@app.post("/api/v1/recommend/conditions")
def api_recommend_conditions(req: RecommendConditionsRequest):
    """Get ML-based condition recommendations.
    
    Returns clean, compact output matching the ui_simple.py format.
    This is robot-friendly and excludes large feature vectors.
    """
    return recommendation_service.recommend_conditions(req)
