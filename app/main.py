from fastapi import FastAPI, HTTPException, Request, Response
from chemtools.contracts import (
    NormalizeRequest, DetectFamilyRequest, FeaturizeUllmannRequest,
    ConditionCoreParseRequest, PropertiesLookupRequest, PrecedentKNNRequest,
    ConstraintsFilterRequest, ExplainPrecedentsRequest
)
from chemtools import smiles, router, featurizers, condition_core, properties, precedent, constraints, explain
import logging, time

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


@app.middleware("http")
async def logging_timing_middleware(request: Request, call_next):
    start = time.perf_counter()
    path = request.url.path
    method = request.method
    try:
        response: Response = await call_next(request)
        status = response.status_code
        return response
    finally:
        dur = time.perf_counter() - start
        # Log without body to avoid PII
        logger.info("%s %s -> %s in %.3f s", method, path, status, dur)
        if REQUEST_COUNT is not None:
            REQUEST_COUNT.labels(method=method, path=path, status=str(status)).inc()
        if REQUEST_LATENCY is not None:
            REQUEST_LATENCY.labels(method=method, path=path).observe(dur)

@app.get("/health")
def health(): return {"ok": True}

@app.post("/api/v1/smiles/normalize")
def api_smiles_normalize(req: NormalizeRequest): return smiles.normalize(req.smiles)

@app.post("/api/v1/router/detect-family")
def api_router_detect(req: DetectFamilyRequest): return router.detect_family(req.reactants)

@app.post("/api/v1/featurize/ullmann")
def api_featurize_ullmann(req: FeaturizeUllmannRequest): return featurizers.ullmann.featurize(req.electrophile, req.nucleophile)

@app.post("/api/v1/condition-core/parse")
def api_condition_core(req: ConditionCoreParseRequest): return condition_core.parse_core(req.reagents, req.text or "")

@app.post("/api/v1/properties/lookup")
def api_properties(req: PropertiesLookupRequest): return properties.lookup(req.query)

@app.post("/api/v1/precedent/knn")
def api_precedent_knn(req: PrecedentKNNRequest): return precedent.knn(req.family, req.features, req.k, req.relax or {})

@app.post("/api/v1/constraints/filter")
def api_constraints_filter(req: ConstraintsFilterRequest): return constraints.apply_filter(req.candidates, req.rules or {})

@app.post("/api/v1/explain/precedents")
def api_explain_precedents(req: ExplainPrecedentsRequest): return explain.for_pack(req.pack, req.features)


@app.get("/metrics")
def metrics():
    if _PROM:
        data = generate_latest()  # type: ignore
        return Response(content=data, media_type=CONTENT_TYPE_LATEST)  # type: ignore
    # Fallback minimal metrics when prometheus_client is unavailable
    return {
        "ok": True,
        "note": "prometheus_client not installed; exposing minimal metrics only",
    }
