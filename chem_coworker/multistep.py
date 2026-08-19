"""Application orchestration for bounded multistep retrosynthesis."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional, Protocol, Sequence

from cas_tools import open_stock_lookup
from core_retrosynthesis import (
    GenericTemplateLibrary,
    MultistepRetrosynthesisRoute,
    plan_multistep_routes,
    recommend_retrosynthesis_conditions,
)

from .contracts import (
    ConditionReviewSettings,
    MultistepRetrosynthesisRequest,
    MultistepRetrosynthesisResponse,
    MultistepReview,
)
from .multistep_rendering import render_multistep_retrosynthesis


_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STOCK_PORTFOLIO_PATH = (
    _REPOSITORY_ROOT / "cas_tools" / "data" / "stock_portfolio.sqlite"
)
DEFAULT_LITERATURE_INDEX_PATH = (
    _REPOSITORY_ROOT / "results" / "literature_molecule_index.sqlite"
)


class _ConditionRecommender(Protocol):
    def recommend(self, reaction_smiles: str, **kwargs):
        """Return deterministic condition evidence."""


class _MultistepReviewer(Protocol):
    def review(
        self,
        routes: Sequence[MultistepRetrosynthesisRoute],
        route_kind: str,
        settings: ConditionReviewSettings,
        guidance: str = "",
    ) -> MultistepReview:
        """Return bounded annotations for existing route candidates."""


@dataclass(frozen=True)
class MultistepRetrosynthesisCoworker:
    """Compose deterministic route search, conditions, and advisory review."""

    library: GenericTemplateLibrary
    stock_path: Path
    library_path: Optional[Path] = None
    condition_recommender: Optional[_ConditionRecommender] = None
    reviewer: Optional[_MultistepReviewer] = None

    @classmethod
    def from_retrosynthesis_coworker(
        cls,
        coworker,
        *,
        stock_path: str | Path | None = None,
        reviewer: Optional[_MultistepReviewer] = None,
    ) -> "MultistepRetrosynthesisCoworker":
        """Reuse the already loaded one-step operator library."""

        selected_stock = cls.default_stock_path() if stock_path is None else Path(stock_path)
        selected_stock = selected_stock.expanduser().resolve()
        if not selected_stock.is_file():
            raise FileNotFoundError(f"stock lookup is unavailable: {selected_stock}")
        return cls(
            library=coworker.library,
            library_path=coworker.library_path,
            stock_path=selected_stock,
            condition_recommender=coworker.condition_recommender,
            reviewer=reviewer,
        )

    @staticmethod
    def default_stock_path() -> Path:
        """Return the strongest available explicit terminal-evidence index."""

        for candidate in (
            DEFAULT_STOCK_PORTFOLIO_PATH,
            DEFAULT_LITERATURE_INDEX_PATH,
        ):
            if candidate.is_file():
                return candidate
        raise FileNotFoundError(
            "No stock lookup is available; expected "
            f"{DEFAULT_STOCK_PORTFOLIO_PATH} or {DEFAULT_LITERATURE_INDEX_PATH}"
        )

    def plan(
        self,
        request: MultistepRetrosynthesisRequest,
    ) -> MultistepRetrosynthesisResponse:
        """Search bounded routes, then optionally review the fixed candidates."""

        search_top_k = request.top_k
        if request.review.mode != "off":
            search_top_k = min(
                10,
                max(request.top_k * 2, request.review.max_candidates),
            )

        condition_evaluator = None
        warnings: list[str] = []
        if request.include_conditions:
            if self.condition_recommender is None:
                warnings.append(
                    "Condition evidence was requested but no condition recommender "
                    "is configured."
                )
            else:

                def condition_evaluator(reaction_smiles: str):
                    return recommend_retrosynthesis_conditions(
                        reaction_smiles,
                        self.condition_recommender,
                        condition_top_k=request.condition_top_k,
                    )

        try:
            with open_stock_lookup(self.stock_path) as stock_index:
                result = plan_multistep_routes(
                    request.target_smiles.strip(),
                    self.library,
                    stock_index,
                    max_depth=request.max_depth,
                    molecular_weight_threshold=request.molecular_weight_threshold,
                    top_k_routes=search_top_k,
                    per_step_top_k=request.per_step_top_k,
                    beam_width=request.beam_width,
                    max_expansions=request.max_expansions,
                    max_templates_to_apply=request.max_templates_to_apply,
                    max_candidates_to_validate=request.max_candidates_to_validate,
                    use_context=request.use_context,
                    include_l0=request.include_l0,
                    condition_evidence_evaluator=condition_evaluator,
                )
        except (OSError, RuntimeError, ValueError) as exc:
            response = MultistepRetrosynthesisResponse(
                request=request,
                valid=False,
                error=str(exc),
                warnings=tuple(warnings),
                library_path=str(self.library_path) if self.library_path else None,
                stock_path=str(self.stock_path),
            )
            return replace(response, answer=render_multistep_retrosynthesis(response))

        routes = result.routes if result.routes else result.partial_routes
        route_kind = "solved" if result.routes else "partial"
        review = None
        if self.reviewer is not None:
            review = self.reviewer.review(
                routes,
                route_kind,
                request.review,
                request.strategic_guidance,
            )
        elif request.review.mode != "off":
            warnings.append("LLM route review was requested but is not configured.")
        if not result.routes:
            warnings.append(
                "No fully solved route was found within the configured search bounds; "
                "partial routes remain explicitly unresolved."
            )
        response = MultistepRetrosynthesisResponse(
            request=request,
            valid=True,
            result=result,
            review=review,
            warnings=tuple(warnings),
            library_path=str(self.library_path) if self.library_path else None,
            stock_path=str(self.stock_path),
        )
        return replace(response, answer=render_multistep_retrosynthesis(response))


__all__ = [
    "DEFAULT_LITERATURE_INDEX_PATH",
    "DEFAULT_STOCK_PORTFOLIO_PATH",
    "MultistepRetrosynthesisCoworker",
]
