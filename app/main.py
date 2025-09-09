from fastapi import FastAPI, HTTPException, Request, Response
from chemtools.contracts import (
    NormalizeRequest, DetectFamilyRequest, FeaturizeUllmannRequest,
    ConditionCoreParseRequest, PropertiesLookupRequest, PrecedentKNNRequest,
    ConstraintsFilterRequest, ExplainPrecedentsRequest, ConditionCoreValidateRequest,
    RecommendFromReactionRequest, PlateDesignRequest,
    RoleAwareMolRequest, RoleAwareReactionRequest,
)
from chemtools import smiles, router, featurizers, condition_core, properties, precedent, constraints, explain, recommend
try:
    from chem_feats import featurize_mol as role_featurize_mol, featurize_reaction as role_featurize_reaction  # type: ignore
    from chem_feats.registry import REGISTRY as ROLE_REGISTRY  # type: ignore
    _HAS_ROLE_FEATS = True
except Exception:
    _HAS_ROLE_FEATS = False
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
def api_featurize_ullmann(req: FeaturizeUllmannRequest, response: Response):
    # Backwards-compatible alias; prefer /api/v1/featurize/molecular
    logger.warning("DEPRECATED endpoint /api/v1/featurize/ullmann; use /api/v1/featurize/molecular")
    try:
        response.headers["X-Deprecated"] = "true"
        response.headers["Link"] = "</api/v1/featurize/molecular>; rel=\"successor-version\""
    except Exception:
        pass
    return featurizers.molecular.featurize(req.electrophile, req.nucleophile)

@app.post("/api/v1/featurize/molecular")
def api_featurize_molecular(req: FeaturizeUllmannRequest):
    return featurizers.molecular.featurize(req.electrophile, req.nucleophile)

@app.post("/api/v1/featurize/role-aware/molecule")
def api_featurize_role_molecule(req: RoleAwareMolRequest):
    if not _HAS_ROLE_FEATS:
        raise HTTPException(status_code=503, detail="role-aware featurization unavailable")
    out = role_featurize_mol(req.smiles, roles=req.roles or None)
    vec = out.get("vector")
    try:
        out["vector"] = vec.tolist()  # type: ignore
    except Exception:
        pass
    return out

@app.post("/api/v1/featurize/role-aware/reaction")
def api_featurize_role_reaction(req: RoleAwareReactionRequest):
    if not _HAS_ROLE_FEATS:
        raise HTTPException(status_code=503, detail="role-aware featurization unavailable")
    out = role_featurize_reaction(req.reaction)
    # Ensure vectors are JSON-serializable lists
    try:
        for item in out.get("reactants") or []:  # type: ignore[union-attr]
            vec = item.get("vector")
            try:
                item["vector"] = vec.tolist()  # type: ignore
            except Exception:
                pass
    except Exception:
        pass
    return out


@app.get("/api/v1/featurize/role-aware/fields")
def api_role_aware_fields(roles: str | None = None):
    """Describe role-aware field order and registry.

    Query param `roles` can be a comma-separated list; defaults to amine,alcohol,aryl_halide.
    """
    if not _HAS_ROLE_FEATS:
        raise HTTPException(status_code=503, detail="role-aware featurization unavailable")
    # Parse and normalize roles
    default_roles = ["amine", "alcohol", "aryl_halide"]
    if roles is None or not str(roles).strip():
        use_roles = default_roles
    else:
        use_roles = [r.strip() for r in str(roles).split(",") if r.strip()]
        # Keep only known roles, preserve order
        known = {"amine", "alcohol", "aryl_halide"}
        use_roles = [r for r in use_roles if r in known]
        if not use_roles:
            use_roles = default_roles

    # Assemble fields: global -> role fields (in order) -> fingerprints per role
    fields: list[str] = []
    fields.extend([f.get("name", "") for f in ROLE_REGISTRY.get("global", [])])
    for r in use_roles:
        fields.extend([f.get("name", "") for f in ROLE_REGISTRY.get(r, [])])
    for r in use_roles:
        bits = int(ROLE_REGISTRY.get("fingerprints", {}).get(r, {}).get("bits", 512))
        fields.extend([f"fp.{r}.{i}" for i in range(bits)])

    counts = {
        "global": len(ROLE_REGISTRY.get("global", [])),
        "by_role": {r: len(ROLE_REGISTRY.get(r, [])) for r in use_roles},
        "fingerprints": {r: int(ROLE_REGISTRY.get("fingerprints", {}).get(r, {}).get("bits", 512)) for r in use_roles},
    }
    return {
        "roles": use_roles,
        "fields": fields,
        "counts": {**counts, "total": len(fields)},
        "registry": ROLE_REGISTRY,
    }


@app.post("/api/v1/featurize/molecule")
def api_featurize_molecule(req: RoleAwareMolRequest):
    """Convenience alias to role-aware featurizer.

    If roles is omitted/null, returns global/basic features only (roles=[]).
    Pass roles=[...] to include role-specific blocks.
    """
    if not _HAS_ROLE_FEATS:
        raise HTTPException(status_code=503, detail="role-aware featurization unavailable")
    roles = req.roles if req.roles is not None else []  # default to globals-only
    out = role_featurize_mol(req.smiles, roles=roles)
    vec = out.get("vector")
    try:
        out["vector"] = vec.tolist()  # type: ignore
    except Exception:
        pass
    return out

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


@app.post("/api/v1/condition-core/validate-dataset")
def api_condition_core_validate(req: ConditionCoreValidateRequest):
    import json, os
    path = req.path
    limit = int(req.limit or 0)
    if not os.path.exists(path):
        raise HTTPException(status_code=400, detail=f"Dataset not found: {path}")

    total = 0
    ok = 0
    mismatches = []

    def _norm_core(s: str) -> str:
        s = (s or "").strip()
        return s[:-5] if s.endswith("/none") else s

    def _metal_part(s: str) -> str:
        s = (s or "").strip()
        return s.split("/", 1)[0] if s else ""

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except Exception:
                continue
            cat = rec.get("catalyst") or {}
            reagents = []
            for item in (cat.get("full_system") or cat.get("core") or []):
                reagents.append({"name": item.get("name"), "uid": item.get("cas"), "role": "CATALYST"})
            for item in (rec.get("reagents") or []):
                reagents.append({"name": item.get("name"), "uid": item.get("cas"), "role": item.get("role") or "ADDITIVE"})
            for item in (rec.get("solvents") or []):
                reagents.append({"name": item.get("name"), "uid": item.get("cas"), "role": "SOLVENT"})

            out = condition_core.parse(reagents, "")
            truth = (rec.get("condition_core") or "").strip()
            pred = (out.get("core") or "").strip()

            t = _norm_core(truth)
            p = _norm_core(pred)
            ok_flag = (t == p) or (req.metal_only_ok and _metal_part(t) and _metal_part(t) == _metal_part(p))

            total += 1
            if ok_flag:
                ok += 1
            elif len(mismatches) < int(req.show_mismatches or 0):
                mismatches.append({
                    "reaction_id": rec.get("reaction_id"),
                    "truth": truth,
                    "pred": pred,
                })
            if limit and total >= limit:
                break

    acc = (ok / total) * 100.0 if total else 0.0
    return {"records": total, "matches": ok, "accuracy": round(acc, 2), "mismatches": mismatches}


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


@app.on_event("startup")
async def warm_startup_caches() -> None:
    # Preload registry and dataset-derived aliases to reduce first-request latency
    try:
        from chemtools import registry as _reg
        _reg._load_registry()  # type: ignore[attr-defined]
    except Exception:
        pass
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


@app.post("/api/v1/recommend")
def api_recommend(req: RecommendFromReactionRequest):
    return recommend.recommend_from_reaction(req.reaction, k=req.k, relax=req.relax or {}, constraint_rules=req.constraints or {})


@app.post("/api/v1/design_plate")
def api_design_plate(req: PlateDesignRequest):
    return recommend.design_plate_from_reaction(req.reaction, plate_size=req.plate_size, relax=req.relax or {}, constraint_rules=req.constraints or {})
    # Optionally preload DRFP fingerprints from NPZ bundle if provided via env
    try:
        import os
        from chemtools.reaction_similarity import load_precomputed_npz  # type: ignore
        path = os.environ.get("CHEMTOOLS_DRFPPATH", "").strip()
        if path and os.path.exists(path):
            res = load_precomputed_npz(path)
            # Silent on failure to avoid impacting startup
            _ = res
    except Exception:
        pass
