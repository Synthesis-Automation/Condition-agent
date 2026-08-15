"""FastAPI composition layer for the local reaction-recommender web app."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

from .contracts import (
    API_SCHEMA_VERSION,
    DiscoveryRequest,
    FeatureAnalysisRequest,
    ForwardSynthesisRequest,
    MultistepRetrosynthesisRequest,
    PrepareReactionRequest,
    RecommendationRequest,
    RetrosynthesisConditionsRequest,
    RetrosynthesisRequest,
    RenderMoleculeRequest,
    RenderReactionRequest,
    envelope,
)
from .runtime import (
    LocalRecommendationRuntime,
    WebRuntime,
    error_payload,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FRONTEND_DIST = (
    PROJECT_ROOT / "web" / "reaction_recommender" / "dist"
)


def create_app(
    *,
    runtime: WebRuntime | None = None,
    frontend_dist: str | Path | None = None,
) -> FastAPI:
    """Create an injectable local API without importing domain logic into UI code."""

    app = FastAPI(
        title="Reaction Condition Recommender API",
        version=API_SCHEMA_VERSION,
        docs_url="/api/docs",
        openapi_url="/api/openapi.json",
    )
    app.add_middleware(
        CORSMiddleware,
        allow_origins=[
            "http://127.0.0.1:5173",
            "http://localhost:5173",
        ],
        allow_credentials=False,
        allow_methods=["GET", "POST"],
        allow_headers=["Content-Type"],
    )
    app.state.runtime = runtime or LocalRecommendationRuntime()

    def active_runtime(request: Request) -> WebRuntime:
        return request.app.state.runtime

    @app.get("/api/v1/health")
    def health() -> dict[str, Any]:
        return envelope({"status": "ok", "local_only": True})

    @app.get("/api/v1/capabilities")
    def capabilities(request: Request) -> dict[str, Any]:
        return envelope(active_runtime(request).capabilities())

    @app.get("/api/v1/ranking-profiles")
    def ranking_profiles(request: Request) -> dict[str, Any]:
        return envelope(
            {"profiles": list(active_runtime(request).ranking_profiles())}
        )

    @app.post("/api/v1/reactions/prepare")
    def prepare_reaction(
        payload: PrepareReactionRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).prepare_reaction(
                payload.reaction_smiles
            )
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=error_payload(exc)) from exc
        return envelope(data)

    @app.post("/api/v1/recommendations")
    def recommend(
        payload: RecommendationRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).recommend(payload)
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            status = 422 if isinstance(exc, ValueError) else 503
            raise HTTPException(status_code=status, detail=error_payload(exc)) from exc
        return envelope(data)

    @app.post("/api/v1/discovery")
    def discover(
        payload: DiscoveryRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).discover(payload)
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            status = 422 if isinstance(exc, ValueError) else 503
            raise HTTPException(status_code=status, detail=error_payload(exc)) from exc
        return envelope(data)

    @app.post("/api/v1/features/analyze")
    def analyze_features(
        payload: FeatureAnalysisRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).analyze_features(payload)
        except (ValueError, RuntimeError) as exc:
            status = 422 if isinstance(exc, ValueError) else 503
            raise HTTPException(
                status_code=status,
                detail=error_payload(exc),
            ) from exc
        return envelope(data)

    @app.post("/api/v1/retrosynthesis")
    def retrosynthesize(
        payload: RetrosynthesisRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).retrosynthesize(payload)
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            status = 422 if isinstance(exc, ValueError) else 503
            raise HTTPException(
                status_code=status,
                detail=error_payload(exc),
            ) from exc
        return envelope(data)

    @app.post("/api/v1/forward-synthesis")
    def forward_synthesize(
        payload: ForwardSynthesisRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).forward_synthesize(payload)
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            status = 422 if isinstance(exc, ValueError) else 503
            raise HTTPException(
                status_code=status,
                detail=error_payload(exc),
            ) from exc
        return envelope(data)

    @app.post("/api/v1/retrosynthesis/conditions")
    def retrosynthesis_conditions(
        payload: RetrosynthesisConditionsRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).retrosynthesis_conditions(payload)
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            status = 422 if isinstance(exc, ValueError) else 503
            raise HTTPException(
                status_code=status,
                detail=error_payload(exc),
            ) from exc
        return envelope(data)

    @app.post("/api/v1/retrosynthesis/routes")
    def multistep_retrosynthesize(
        payload: MultistepRetrosynthesisRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).multistep_retrosynthesize(payload)
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            status = 422 if isinstance(exc, ValueError) else 503
            raise HTTPException(
                status_code=status,
                detail=error_payload(exc),
            ) from exc
        return envelope(data)

    @app.post("/api/v1/render/reaction")
    def render_reaction(
        payload: RenderReactionRequest,
        request: Request,
    ) -> Response:
        try:
            drawing = active_runtime(request).render_reaction(
                payload.reaction_smiles,
                width=payload.width,
                height=payload.height,
            )
        except (ValueError, RuntimeError) as exc:
            raise HTTPException(status_code=422, detail=error_payload(exc)) from exc
        return Response(
            content=drawing,
            media_type="image/svg+xml",
            headers={"Cache-Control": "no-store"},
        )

    @app.post("/api/v1/render/molecule")
    def render_molecule(
        payload: RenderMoleculeRequest,
        request: Request,
    ) -> Response:
        try:
            drawing = active_runtime(request).render_molecule(
                payload.molecule_smiles,
                width=payload.width,
                height=payload.height,
            )
        except (ValueError, RuntimeError) as exc:
            raise HTTPException(
                status_code=422,
                detail=error_payload(exc),
            ) from exc
        return Response(
            content=drawing,
            media_type="image/svg+xml",
            headers={"Cache-Control": "no-store"},
        )

    @app.exception_handler(Exception)
    async def unexpected_error(
        _request: Request,
        exc: Exception,
    ) -> JSONResponse:
        return JSONResponse(
            status_code=500,
            content={"detail": error_payload(exc)},
        )

    dist = Path(frontend_dist) if frontend_dist is not None else DEFAULT_FRONTEND_DIST
    if dist.is_dir():
        app.mount("/", StaticFiles(directory=dist, html=True), name="frontend")

    return app


app = create_app()


__all__ = ["app", "create_app"]
