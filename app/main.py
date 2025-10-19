from fastapi import FastAPI, HTTPException, Request, Response
import os
from functools import lru_cache
from pathlib import Path
from chemtools.contracts import (
    NormalizeRequest, DetectFamilyRequest, FeaturizeUllmannRequest,
    ConditionCoreParseRequest, PrecedentKNNRequest,
    ConstraintsFilterRequest, ExplainPrecedentsRequest, ConditionCoreValidateRequest,
    RecommendFromReactionRequest, RecommendConditionsRequest,
    PlateDesignRequest,
    CoreSearchRequest,
    RoleAwareMolRequest, RoleAwareReactionRequest,
    DetectTypeRequest,
    SchemeMatchRequest,
)
# Enable RDKit by default unless explicitly disabled
os.environ.setdefault("CHEMTOOLS_DISABLE_RDKIT", "0")

# ChemTools v2.0 - Unified API
from chemtools import chem, output_formatter

# Import specific modules that don't have chem.* equivalents yet
from chemtools import featurizers, condition_core

# Import service layer
from app.services import (
    matching_service,
    featurization_service,
    recommendation_service,
    rule_matching_service,
    precedent_service,
)

try:
    from chemtools.reaction_type_detector import detect_reaction_type as rxn_detect_type, is_available as rxn_insight_available  # type: ignore
    _HAS_RXN_INSIGHT = True
except Exception:
    _HAS_RXN_INSIGHT = False
try:
    from chemtools.features.role import (
        featurize_mol as role_featurize_mol,
        featurize_reaction as role_featurize_reaction,
    )  # type: ignore
    from chemtools.features.role.registry import REGISTRY as ROLE_REGISTRY  # type: ignore
    _HAS_ROLE_FEATS = True
except Exception:
    _HAS_ROLE_FEATS = False
import logging, time

# Import error handlers
from app.error_handlers import register_error_handlers
from chemtools.exceptions import (
    ValidationError,
    DatabaseNotFoundError,
    ProcessingError,
    SchemeConditionDBError,
)

# SCDB integration is now part of ChemTools v2.0 via chem.rules.*
try:
    from chemtools.rule_scdb_matcher import loader as scdb_loader
    _HAS_SCDB = True
except Exception:
    _HAS_SCDB = False


_SCDB_DEFAULT_DB = (
    os.environ.get("SCDB_MATCH_DB_PATH", "cn_coupling_pd_db.json").strip()
    or "cn_coupling_pd_db.json"
)


# Optional Prometheus metrics
try:
    from prometheus_client import Counter, Histogram, generate_latest, CONTENT_TYPE_LATEST  # type: ignore
    _PROM = True
except Exception:
    _PROM = False

if _PROM:
    REQUEST_COUNT = Counter(
        'chemtools_requests_total', 'Total HTTP requests', ['method', 'path', 'status']
    )
    REQUEST_LATENCY = Histogram(
        'chemtools_request_latency_seconds', 'Latency per request', ['method', 'path']
    )
else:
    REQUEST_COUNT = None
    REQUEST_LATENCY = None

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
        if REQUEST_COUNT is not None and status is not None:
            REQUEST_COUNT.labels(method=method, path=path, status=str(status)).inc()
        if REQUEST_LATENCY is not None:
            REQUEST_LATENCY.labels(method=method, path=path).observe(dur)

@app.get("/health")
def health(): return {"ok": True}


@app.post("/match")
def api_scdb_match(req: SchemeMatchRequest):
    """Match reaction to rule-based scheme database using ChemTools v2.0."""
    return rule_matching_service.match_reaction(req)

@app.post("/api/v1/smiles/normalize")
def api_smiles_normalize(req: NormalizeRequest): 
    return matching_service.normalize_smiles(req)

@app.post("/api/v1/router/detect-family")
def api_router_detect(req: DetectFamilyRequest): 
    return matching_service.detect_family(req)

@app.post("/api/v1/reaction/detect-type")
def api_detect_type(req: DetectTypeRequest):
    """Detect reaction type using rxn-insight when available, with router fallback."""
    return matching_service.detect_reaction_type(req)

@app.post("/api/v1/featurize/molecular")
def api_featurize_molecular(req: FeaturizeUllmannRequest):
    return featurization_service.featurize_molecular(req)

@app.post("/api/v1/featurize/role-aware/molecule")
def api_featurize_role_molecule(req: RoleAwareMolRequest):
    return featurization_service.featurize_role_aware_molecule(req)

@app.post("/api/v1/featurize/role-aware/reaction")
def api_featurize_role_reaction(req: RoleAwareReactionRequest):
    return featurization_service.featurize_role_aware_reaction(req)


@app.get("/api/v1/featurize/role-aware/fields")
def api_role_aware_fields(roles: str | None = None):
    """Describe role-aware field order and registry.

    Query param `roles` can be a comma-separated list; defaults to amine,alcohol,aryl_halide.
    """
    return featurization_service.get_role_aware_fields(roles)


@app.post("/api/v1/featurize/molecule")
def api_featurize_molecule(req: RoleAwareMolRequest):
    """Convenience alias to role-aware featurizer.

    If roles is omitted/null, returns global/basic features only (roles=[]).
    Pass roles=[...] to include role-specific blocks.
    """
    # Use featurization service with default roles=[] for backwards compatibility
    req.roles = req.roles if req.roles is not None else []
    return featurization_service.featurize_role_aware_molecule(req)

@app.post("/api/v1/condition-core/parse")
def api_condition_core(req: ConditionCoreParseRequest): 
    return precedent_service.parse_condition_core(req)

@app.get("/api/v1/precedent/cores")
def api_core_list(family: str | None = None, limit: int = 200, counts: bool = True):
    return precedent_service.list_cores(family, limit, counts)

@app.post("/api/v1/precedent/knn")
def api_precedent_knn(req: PrecedentKNNRequest): 
    return precedent_service.knn_search(req)

@app.post("/api/v1/constraints/filter")
def api_constraints_filter(req: ConstraintsFilterRequest): 
    return precedent_service.filter_constraints(req)

@app.post("/api/v1/explain/precedents")
def api_explain_precedents(req: ExplainPrecedentsRequest): 
    return precedent_service.explain_precedents(req)


@app.post("/api/v1/condition-core/validate-dataset")
def api_condition_core_validate(req: ConditionCoreValidateRequest):
    return precedent_service.validate_dataset(req)


@app.get("/metrics")
def metrics():
    if _PROM:
        data = generate_latest(    )  # type: ignore
        return Response(content=data, media_type=CONTENT_TYPE_LATEST    )  # type: ignore
    # Fallback minimal metrics when prometheus_client is unavailable
    return {
        "ok": True,
        "note": "prometheus_client not installed; exposing minimal metrics only",
    }


@app.on_event("startup")
async def warm_startup_caches() -> None:
    # Preload dataset-derived aliases to reduce first-request latency
    try:
        from chemtools import condition_core as _cc
        # Touch module-level caches to force initialization
        _ = _cc._LIG_BY_CAS  # type: ignore[attr-defined]
        _ = _cc._LIG_BY_NAME  # type: ignore[attr-defined]
    except Exception:
        pass
    try:
        # Optional: load DRFP precomputed vectors is handled above; nothing else here
        pass
    except Exception:
        pass


@app.post("/api/v1/recommend/conditions")
def api_recommend_conditions(req: RecommendConditionsRequest):
    """Get ML-based condition recommendations.
    
    Returns clean, compact output matching the ui_simple.py format.
    This is robot-friendly and excludes large feature vectors.
    """
    return recommendation_service.recommend_conditions(req)


@app.post("/api/v1/recommend")
def api_recommend(req: RecommendFromReactionRequest):
    return recommendation_service.recommend_from_reaction(req)


@app.post("/api/v1/design_plate")
def api_design_plate(req: PlateDesignRequest):
    return recommendation_service.design_plate(req)


@app.post("/api/v1/core/search")
def api_core_search(req: CoreSearchRequest):
    return precedent_service.search_by_core(req)
