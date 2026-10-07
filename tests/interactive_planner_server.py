"""Local browser-test server: real planning API with deterministic search fixtures.

Run only for UI tests; the candidate fixtures are not chemistry validations.
"""

from __future__ import annotations

import argparse
import time
from typing import Any

from app.web_api.main import create_app
from app.web_api.runtime import LocalRecommendationRuntime
from tests.test_web_recommender_api import FakeRuntime
from tests.core_retrosynthesis_tests.test_interactive_planning import search_result


class BrowserPlanningRuntime(FakeRuntime):
    """Keep chemistry searches predictable while exercising real route edits."""

    render_molecule = LocalRecommendationRuntime.render_molecule
    render_reaction = LocalRecommendationRuntime.render_reaction

    def capabilities(self) -> dict[str, Any]:
        return {
            **super().capabilities(),
            "planner_fixture": True,
            "stock_portfolio_available": True,
            "retrosynthesis_library_modes": {
                "full": {"library_available": True},
                "compact": {"library_available": True},
            },
        }

    def retrosynthesize(self, request: Any) -> dict[str, Any]:
        """Return two root alternatives, repeat occurrences, cycles, or empty hits."""
        if request.target_smiles == "CCN":
            time.sleep(2)
            raise RuntimeError("Fixture search unavailable; retry is safe")
        choices = {
            "CCOC": ("CCO.CI", "CCBr.CO"),
            "CCO": ("CCBr.O",),
            "CCBr": ("C.CBr", "CCOC"),
            "CC": ("C.C",),
        }
        # The first CCBr proposal must be acyclic: methyl bromide is CBr, not CCBr.
        result = search_result(
            request.target_smiles, *choices.get(request.target_smiles, ())
        )
        for strategy in result["strategies"]:
            candidate = strategy["representative"]
            candidate.update(
                abstraction_level="L2",
                transformation_kind="fixture proposal",
                forward_assessment=None,
                supporting_precedents=[],
            )
        if request.target_smiles == "CCOC":
            alternate = search_result("CCOC", "CCCl.CO")["strategies"][0][
                "representative"
            ]
            alternate["abstraction_level"] = "L1"
            result["strategies"][0]["alternate_realizations"] = [alternate]
        result.update(
            library_mode=request.library_mode,
            valid=bool(result["strategies"]),
            search_diagnostics={"budget_limited": False, "levels_attempted": ["L2"]},
        )
        return result

    def retrosynthesis_conditions(self, request: Any) -> dict[str, Any]:
        """Return explicit absence of condition support, without fabrication."""
        return dict(
            status="insufficient_evidence",
            recommendations=[],
            warnings=["No compatible recipe in this test fixture"],
            query_reaction_smiles=request.reaction_smiles,
        )

    def planning_stock(self, smiles: str) -> dict[str, Any]:
        """Provide explicit saved-stock evidence for the browser interaction."""
        return dict(
            smiles=smiles,
            status="verified_stock_match",
            checked_at="2026-10-07T00:00:00+00:00",
            source_records=[
                {"supplier": "Fixture supplier", "terminal_eligible": "true"}
            ],
        )


if __name__ == "__main__":
    import uvicorn

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=5184)
    args = parser.parse_args()
    uvicorn.run(
        create_app(runtime=BrowserPlanningRuntime(), recommendation_only=False),
        host="127.0.0.1",
        port=args.port,
    )
