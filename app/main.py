"""
FastAPI application for Chemistry Tools.

Provides HTE recommendation endpoint that returns the automation-ready JSON
format used by the HTE recommender GUI export.
"""

from __future__ import annotations

import time
from typing import Optional, Tuple, Dict, Any, List

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response
from pydantic import BaseModel, Field

from chemtools.core.errors import ValidationError
from chemtools.formatters import format_hte_output
from chemtools.recommend import HTERecommender
from app.services.error_handlers import register_error_handlers


app = FastAPI(title="ChemTools API", version="1.0.0")
register_error_handlers(app)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/favicon.ico")
def favicon() -> Response:
    return Response(status_code=204)


class HTERecommendationInput(BaseModel):
    reaction_smiles: str = Field(..., description="Reaction SMILES (A.B>>P or A.B or A)")
    reaction_type_filter: Optional[str] = Field(
        None, description="Optional reaction family filter (e.g., C_N_Coupling)"
    )
    catalyst_filter: Optional[str] = Field(
        None, description="Optional catalyst filter (e.g., Pd, Cu)"
    )
    reaction_key_only: bool = Field(
        False, description="Only match using reaction_key/signatures; disable reactant-type fallback"
    )
    top_k: int = Field(10, ge=1, le=200, description="Number of recommendations to return")
    min_experiments: int = Field(2, ge=1, le=200, description="Minimum experiments per condition")


class HTERecommendationRequest(BaseModel):
    input: HTERecommendationInput


def _parse_reaction_smiles(reaction_smiles: str) -> Tuple[str, Optional[str], Optional[str]]:
    text = (reaction_smiles or "").strip()
    if not text:
        return "", None, None
    if ">>" in text:
        reactants_part, product = text.split(">>", 1)
        reactants = [r for r in reactants_part.split(".") if r]
        reactant_a = reactants[0] if reactants else ""
        reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
        product_smiles = product.strip() or None
        return reactant_a, reactant_b, product_smiles
    reactants = [r for r in text.split(".") if r]
    reactant_a = reactants[0] if reactants else ""
    reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
    return reactant_a, reactant_b, None


def _normalize_source_group_label(value: Any) -> str:
    text = str(value or "").strip().lower()
    if not text or text == "nan":
        return "unknown"
    if text in ("literature", "datasets", "dataset", "lit"):
        return "literature"
    if text in ("motif", "motifs", "experiments", "experiment", "experiements"):
        return "motif"
    if text == "rules":
        return "rules"
    return text


@app.post("/hte/recommendations")
def recommend_hte(payload: HTERecommendationRequest) -> Dict[str, Any]:
    request_input = payload.input
    reaction_smiles = request_input.reaction_smiles.strip()
    if not reaction_smiles:
        raise ValidationError("reaction_smiles cannot be empty")

    reactant_a, reactant_b, product = _parse_reaction_smiles(reaction_smiles)
    if not reactant_a:
        raise ValidationError("Provide a reaction SMILES with at least one reactant.")

    start = time.perf_counter()
    recommender = HTERecommender()
    result = recommender.recommend(
        reactant_a_smiles=reactant_a,
        reactant_b_smiles=reactant_b,
        product_smiles=product,
        top_k=request_input.top_k,
        min_experiments=request_input.min_experiments,
        reaction_type_filter=request_input.reaction_type_filter or None,
        catalyst_filter=request_input.catalyst_filter or None,
        reaction_key_only=request_input.reaction_key_only,
    )

    source_map = getattr(result, "recommendations_by_source", {}) or {}
    normalized_map: Dict[str, List[Any]] = {}
    for key, items in source_map.items():
        normalized_key = _normalize_source_group_label(key)
        normalized_map.setdefault(normalized_key, []).extend(items)
    literature_recs = list(normalized_map.get("literature") or [])
    if not literature_recs:
        literature_recs = list(getattr(result, "recommendations", []) or [])

    output = format_hte_output(
        result,
        recommendations=literature_recs,
        reaction_smiles=reaction_smiles,
        reaction_type_filter=request_input.reaction_type_filter,
        catalyst_filter=request_input.catalyst_filter,
        explanation=None,
    )

    processing_ms = round((time.perf_counter() - start) * 1000, 2)
    output.setdefault("meta", {})["processing_time_ms"] = processing_ms
    return output
