from fastapi import FastAPI, HTTPException, Request, Response
import os
from functools import lru_cache
from pathlib import Path
from chemtools.contracts import (
    NormalizeRequest, DetectFamilyRequest, FeaturizeUllmannRequest,
    ConditionCoreParseRequest, PropertiesLookupRequest, PrecedentKNNRequest,
    ConstraintsFilterRequest, ExplainPrecedentsRequest, ConditionCoreValidateRequest,
    RecommendFromReactionRequest, RecommendConditionsRequest, PlateDesignRequest,
    CoreSearchRequest,
    RoleAwareMolRequest, RoleAwareReactionRequest,
    DetectTypeRequest,
    SchemeMatchRequest,
)
# Enable RDKit by default unless explicitly disabled
os.environ.setdefault("CHEMTOOLS_DISABLE_RDKIT", "0")

# ChemTools v2.0 - Unified API
from chemtools import chem

# Deprecated: Keep old module imports for gradual migration
# TODO: Remove these imports once all code uses chem.*
from chemtools import smiles, router, featurizers, condition_core, properties, precedent, constraints, explain, recommend
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

# SCDB integration is now part of ChemTools v2.0 via chem.rules.*
# Keeping backward compatibility for error handling
try:
    from chemtools.scdb_matcher.loader import SchemeConditionDBError
    _HAS_SCDB = True
except Exception:
    _HAS_SCDB = False
    
    class SchemeConditionDBError(RuntimeError):
        pass


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


@app.post("/match")
def api_scdb_match(req: SchemeMatchRequest):
    """Match reaction to rule-based scheme database using ChemTools v2.0."""
    if not _HAS_SCDB:
        raise HTTPException(status_code=503, detail="SchemeConditionDB matcher unavailable")
    
    reaction = (req.reaction or "").strip()
    if not reaction:
        raise HTTPException(status_code=400, detail="Reaction must be a non-empty string")
    
    db_path = (req.db or _SCDB_DEFAULT_DB).strip()
    if not db_path:
        raise HTTPException(status_code=400, detail="No database path configured for matching")
    
    try:
        # Use ChemTools v2.0 chem.rules.* API
        db = chem.rules.load_database(db_path, cache=True)
        result = chem.rules.match(db, reaction)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=f"Database file not found: {db_path}") from exc
    except SchemeConditionDBError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Unexpected error while matching SchemeConditionDB")
        raise HTTPException(status_code=500, detail="Internal error while performing SchemeConditionDB match") from exc
    
    payload = result.to_json_dict()
    if not req.include_trace:
        payload.pop("trace", None)
    return payload

@app.post("/api/v1/smiles/normalize")
def api_smiles_normalize(req: NormalizeRequest): return chem.smiles.normalize(req.smiles)

@app.post("/api/v1/router/detect-family")
def api_router_detect(req: DetectFamilyRequest): return router.detect_family(req.reactants)  # Note: router.detect_family expects list, not reaction string

@app.post("/api/v1/reaction/detect-type")
def api_detect_type(req: DetectTypeRequest):
    """Detect reaction type using rxn-insight when available, with router fallback.

    Returns a combined payload including rxn-insight results, router fallback, and the selected family.
    """
    rxn = req.reaction
    norm = chem.smiles.normalize_reaction(rxn)
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]
    fallback = router.detect_family(reactants)  # Using list of reactants as expected
    auto = None
    if _HAS_RXN_INSIGHT:
        try:
            auto = rxn_detect_type(norm.get("normalized") or rxn)
        except Exception:
            auto = None
    selected = None
    if isinstance(auto, dict) and (auto.get("mapped_family") or auto.get("success")):
        selected = auto.get("mapped_family") or fallback.get("family")
    else:
        selected = fallback.get("family")
    return {
        "input": {"reaction_smiles": norm.get("normalized") or rxn},
        "rxn_insight_available": bool(_HAS_RXN_INSIGHT),
        "rxn_insight": auto,
        "router_fallback": fallback,
        "selected_family": selected,
    }

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
        out["vector"] = vec.tolist(    )  # type: ignore
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
                item["vector"] = vec.tolist(    )  # type: ignore
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
        out["vector"] = vec.tolist(    )  # type: ignore
    except Exception:
        pass
    return out

@app.post("/api/v1/condition-core/parse")
def api_condition_core(req: ConditionCoreParseRequest): return condition_core.parse_core(req.reagents, req.text or "")

@app.post("/api/v1/properties/lookup")
def api_properties(req: PropertiesLookupRequest): return chem.properties.lookup(req.query)

@app.get("/api/v1/cores")
def api_core_list(family: str | None = None, limit: int = 200, counts: bool = True):
    data = chem.precedent.list_cores(family=family, top_n=int(limit or 200), include_counts=bool(counts))
    return {"cores": data}

@app.post("/api/v1/precedent/knn")
def api_precedent_knn(req: PrecedentKNNRequest): return chem.precedent.knn(req.family, req.features, req.k, req.relax or {})

@app.post("/api/v1/constraints/filter")
def api_constraints_filter(req: ConstraintsFilterRequest): return chem.constraints.filter(req.candidates, req.rules or {})

@app.post("/api/v1/explain/precedents")
def api_explain_precedents(req: ExplainPrecedentsRequest): return chem.explain.for_pack(req.pack, req.features)


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
    import time
    from chemtools import output_formatter
    
    start_time = time.time()
    
    # Get raw ML recommendations
    raw_data = chem.recommend.recommend_conditions_structured(
        reaction=req.reaction,
        reaction_type=req.reaction_type,
        k=req.k,
        limit=req.limit,
        relax=req.relax or {},
        constraints=req.constraints or {},
    )
    
    # Extract detection info
    detection = raw_data.get("detection", {})
    detected_type = detection.get("reaction_type", "Unknown")
    confidence = detection.get("auto", {}).get("confidence", 0.0) if detection.get("source") == "auto" else 1.0
    
    # Format recommendations using the UI formatter
    recommendations_data = []
    for rec in raw_data.get("recommendations", [])[:req.limit]:
        summary = rec.get("summary", {})
        chemicals = rec.get("chemicals", [])
        conditions_info = rec.get("conditions", {})
        
        # Build reagents list
        reagents = []
        for chem in chemicals:
            reagents.append({
                "id": chem.get("uid", chem.get("cas")),
                "role": chem.get("role", "reagent"),
                "name": chem.get("name"),
                "abbreviation": None,
                "cas": chem.get("cas"),
                "smiles": chem.get("smiles"),
                "equivalents": None,  # Not in summary
            })
        
        # Format conditions
        conditions = {}
        if conditions_info.get("temperature"):
            conditions["temperature"] = {
                "value": conditions_info["temperature"],
                "unit": "°C"
            }
        if conditions_info.get("time"):
            conditions["time"] = {
                "value": conditions_info["time"],
                "unit": "hours"
            }
        if conditions_info.get("atmosphere"):
            conditions["atmosphere"] = conditions_info["atmosphere"]
        
        formatted_rec = {
            "rank": rec.get("rank", len(recommendations_data) + 1),
            "confidence": summary.get("confidence", 0.0) / 100.0 if summary.get("confidence") else 0.0,  # Convert % to decimal
            "reagents": reagents,
            "conditions": conditions,
            "precedent_count": summary.get("support", {}).get("count") if isinstance(summary.get("support"), dict) else summary.get("support", 0),
        }
        
        recommendations_data.append(formatted_rec)
    
    # Build compact output using UI formatter
    processing_time_ms = (time.time() - start_time) * 1000
    
    output = output_formatter.format_ml_output(
        reaction_smiles=req.reaction,
        requested_type=req.reaction_type,
        detected_type=detected_type,
        detection_confidence=confidence,
        recommendations_data=recommendations_data,
        processing_time_ms=processing_time_ms,
    )
    
    return output




@app.post("/api/v1/recommend")
def api_recommend(req: RecommendFromReactionRequest):
    return chem.recommend.recommend_from_reaction(req.reaction, k=req.k, relax=req.relax or {}, constraint_rules=req.constraints or {})


@app.post("/api/v1/design_plate")
def api_design_plate(req: PlateDesignRequest):
    return chem.recommend.design_plate_from_reaction(req.reaction, plate_size=req.plate_size, relax=req.relax or {}, constraint_rules=req.constraints or {})
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


@app.post("/api/v1/core/search")
def api_core_search(req: CoreSearchRequest):
    rows = precedent.find_reactions_by_core(
        req.core,
        family=req.family,
        fuzzy=bool(req.fuzzy),
        limit=int(req.limit or 50),
    )
    return {
        "query": {"core": req.core, "family": req.family, "fuzzy": bool(req.fuzzy), "limit": int(req.limit or 50)},
        "count": len(rows),
        "results": rows,
    }
