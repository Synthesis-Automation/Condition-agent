"""FastAPI composition layer for the local reaction-recommender web app."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles

from .contracts import (
    API_SCHEMA_VERSION,
    AssistanceConfirmationRequest,
    AssistanceSessionRequest,
    CoupledStrategyRetrosynthesisRequest,
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
from .conditions import router as conditions_router


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FRONTEND_DIST = PROJECT_ROOT / "web" / "reaction_recommender" / "dist"


def create_app(
    *,
    runtime: WebRuntime | None = None,
    assistance_service: Any | None = None,
    frontend_dist: str | Path | None = None,
    recommendation_only: bool = True,
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
    app.state.assistance_service = assistance_service
    profile = "recommendation_only" if recommendation_only else "research_workbench"
    app.state.deployment_profile = profile
    app.include_router(conditions_router)

    def active_runtime(request: Request) -> WebRuntime:
        return request.app.state.runtime

    @app.get("/api/v1/health")
    def health() -> dict[str, Any]:
        return envelope(
            {"status": "ok", "local_only": True, "deployment_profile": profile}
        )

    @app.get("/api/v1/capabilities")
    def capabilities(request: Request) -> dict[str, Any]:
        data = active_runtime(request).capabilities()
        if recommendation_only:
            data = {
                key: value
                for key, value in data.items()
                if key
                in {
                    "service",
                    "index_name",
                    "index_available",
                    "loaded_runtime_variants",
                    "rxnmapper_available",
                    "recommendation",
                    "weak_label_recommendation",
                    "weak_label_dataset_name",
                    "reaction_rendering",
                    "local_only",
                }
            }
        data = {**data, "deployment_profile": profile}
        return envelope(data)

    @app.get("/api/v1/ranking-profiles")
    def ranking_profiles(request: Request) -> dict[str, Any]:
        return envelope({"profiles": list(active_runtime(request).ranking_profiles())})

    @app.post("/api/v1/reactions/prepare")
    def prepare_reaction(
        payload: PrepareReactionRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).prepare_reaction(payload.reaction_smiles)
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

    @app.post("/api/v1/experimental/assistance")
    def start_assistance(
        payload: AssistanceSessionRequest,
        request: Request,
    ) -> dict[str, Any]:
        service = request.app.state.assistance_service
        if service is None:
            raise HTTPException(
                status_code=503,
                detail={
                    "error": "ASSISTANCE_NOT_CONFIGURED",
                    "message": "Experimental advisory assistance is disabled.",
                },
            )
        try:
            data = service.start(
                objective=payload.objective,
                mode=payload.mode,
                structure_input=payload.structure_input,
                provider_settings=payload.provider.model_dump(),
            )
        except (ValueError, RuntimeError) as exc:
            raise HTTPException(status_code=422, detail=error_payload(exc)) from exc
        return envelope(data)

    @app.post("/api/v1/experimental/assistance/confirm-condition")
    def confirm_assistance_condition(
        payload: AssistanceConfirmationRequest,
        request: Request,
    ) -> dict[str, Any]:
        service = request.app.state.assistance_service
        if service is None:
            raise HTTPException(
                status_code=503,
                detail={
                    "error": "ASSISTANCE_NOT_CONFIGURED",
                    "message": "Experimental advisory assistance is disabled.",
                },
            )
        try:
            data = service.confirm_condition_constraint(
                state=payload.state,
                raw_value=payload.raw_value,
            )
        except (ValueError, RuntimeError) as exc:
            raise HTTPException(status_code=422, detail=error_payload(exc)) from exc
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

    @app.get("/api/v1/forward-synthesis/condition-profiles")
    def forward_condition_profiles(request: Request) -> dict[str, Any]:
        return envelope(active_runtime(request).forward_condition_profiles())

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

    @app.post("/api/v1/retrosynthesis/coupled-strategies")
    def coupled_strategy_retrosynthesize(
        payload: CoupledStrategyRetrosynthesisRequest,
        request: Request,
    ) -> dict[str, Any]:
        try:
            data = active_runtime(request).coupled_strategy_retrosynthesize(payload)
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

    if recommendation_only:
        # Explicit allowlist: research capabilities cannot be invoked through
        # the focused deployment, even when their optional libraries exist.
        allowed = {
            "/api/docs",
            "/docs/oauth2-redirect",
            "/api/openapi.json",
            "/api/v1/health",
            "/api/v1/capabilities",
            "/api/v1/conditions/recommend",
            "/api/v1/reactions/prepare",
            "/api/v1/render/reaction",
            "/api/v1/render/molecule",
        }
        app.router.routes[:] = [route for route in app.routes if route.path in allowed]
        app.title = "Condition Recommendations"

    dist = Path(frontend_dist) if frontend_dist is not None else DEFAULT_FRONTEND_DIST
    entry = dist / ("index.html" if recommendation_only else "workbench.html")

    @app.get("/", include_in_schema=False)
    def frontend() -> Response:
        if not entry.is_file():
            return JSONResponse(
                status_code=503,
                content={
                    "detail": "Frontend build is missing. Run python -m app.web_api --build"
                },
                headers={"Cache-Control": "no-store"},
            )
        return FileResponse(entry, headers={"Cache-Control": "no-store"})

    @app.get("/index.html", include_in_schema=False)
    @app.get("/workbench.html", include_in_schema=False)
    def canonical_frontend() -> Response:
        return RedirectResponse(
            "/", status_code=307, headers={"Cache-Control": "no-store"}
        )

    if (dist / "assets").is_dir():
        app.mount(
            "/assets", StaticFiles(directory=dist / "assets"), name="frontend_assets"
        )

    return app


app = create_app()


__all__ = ["app", "create_app"]
