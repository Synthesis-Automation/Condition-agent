from fastapi import FastAPI, HTTPException
from chemtools.contracts import (
    NormalizeRequest, DetectFamilyRequest, FeaturizeUllmannRequest,
    ConditionCoreParseRequest, PropertiesLookupRequest, PrecedentKNNRequest,
    ConstraintsFilterRequest, ExplainPrecedentsRequest
)
from chemtools import smiles, router, featurizers, condition_core, properties, precedent, constraints, explain

app = FastAPI(title="Chemistry Tools API", version="0.1.1")

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
