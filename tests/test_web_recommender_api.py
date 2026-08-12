"""Contract tests for the local web application boundary."""

from __future__ import annotations

import gzip
import json
from types import SimpleNamespace
from typing import Any, Dict

from fastapi.testclient import TestClient

from app.web_api.main import create_app
from app.web_api.contracts import (
    MultistepRetrosynthesisRequest,
    RetrosynthesisConditionsRequest,
    RetrosynthesisRequest,
)
import app.web_api.runtime as runtime_module
from app.web_api.runtime import LocalRecommendationRuntime


class FakeRuntime:
    def capabilities(self) -> Dict[str, Any]:
        return {
            "index_available": True,
            "rxnmapper_available": False,
            "featurization": True,
            "local_only": True,
        }

    def ranking_profiles(self) -> tuple[Dict[str, Any], ...]:
        return (
            {
                "profile_id": "default",
                "label": "Balanced",
                "description": "Balanced evidence",
                "weights": {"similarity": 1.0},
            },
        )

    def prepare_reaction(self, reaction_smiles: str) -> Dict[str, Any]:
        if reaction_smiles == "invalid":
            raise ValueError("INVALID_REACTION")
        return {
            "valid": True,
            "input_reaction_smiles": reaction_smiles,
            "completion_proposal": {
                "status": "not_required",
                "requirements": [],
            },
            "warnings": [],
        }

    def recommend(self, request: Any) -> Dict[str, Any]:
        return {
            "query_reaction_smiles": request.reaction_smiles,
            "library_mode": request.library_mode,
            "valid": True,
            "recommendations": [
                {
                    "rank": 1,
                    "recipe_id": "recipe:1",
                    "score": 0.9,
                    "synthesis_protocol": {
                        "protocol_id": "CP1:test",
                        "materials": [
                            {
                                "category": "condition",
                                "canonical_name": "Potassium carbonate",
                                "cas": "584-08-7",
                            }
                        ],
                        "execution_readiness": "review_required",
                    },
                }
            ],
            "schema_version": "test",
        }

    def discover(self, request: Any) -> Dict[str, Any]:
        return {
            "query_reaction_smiles": request.reaction_smiles,
            "library_mode": request.library_mode,
            "valid": True,
            "hits": [{"rank": 1, "reaction_id": "rxn:1"}],
            "schema_version": "test",
        }

    def analyze_features(self, request: Any) -> Dict[str, Any]:
        input_kind = "reaction" if ">" in request.input_smiles else "molecule"
        return {
            "input_kind": input_kind,
            "input_smiles": request.input_smiles,
            "valid": True,
            "overview": {"atom_count": 6},
            "motifs": [],
            "reactive_sites": [],
            "reaction_core": None,
            "partners": [],
            "mapping": None,
            "core_graphic_svg": None,
            "core_projection": None,
            "warnings": [],
            "error": None,
            "analysis": {"schema_version": "test"},
        }

    def retrosynthesize(self, request: Any) -> Dict[str, Any]:
        return {
            "target_smiles": request.target_smiles,
            "library_mode": request.library_mode,
            "valid": True,
            "error": None,
            "schema_version": "1.0",
            "candidate_count": 1,
            "library_operator_count": 12,
            "library_template_count": 34,
            "warnings": [],
            "candidates": [
                {
                    "rank": 1,
                    "precursor_smiles": "CCBr.N",
                    "proposed_reaction_smiles": "CCBr.N>>CCN",
                    "score": 0.91,
                    "abstraction_level": "L2",
                    "forward_validation_status": "verified_signature",
                }
            ],
        }

    def retrosynthesis_conditions(self, request: Any) -> Dict[str, Any]:
        return {
            "status": "recommended_direct",
            "query_reaction_smiles": request.reaction_smiles,
            "recommender_valid": True,
            "recommendation_mode": "verified_signature",
            "retrieval_level": "L2",
            "uses_type_agnostic_fallback": False,
            "candidate_count": 1,
            "independent_candidate_count": 1,
            "compatible_candidate_count": 1,
            "independent_compatible_candidate_count": 1,
            "excluded_candidate_count": 0,
            "best_recipe_score": 0.9,
            "best_recipe_compatibility_score": 1.0,
            "best_recipe_reference_support": 1,
            "recommendations": [],
            "warnings": [],
            "error": None,
        }

    def multistep_retrosynthesize(self, request: Any) -> Dict[str, Any]:
        return {
            "target_smiles": request.target_smiles,
            "library_mode": request.library_mode,
            "valid": True,
            "error": None,
            "schema_version": "1.0",
            "max_depth": request.max_depth,
            "molecular_weight_threshold": (request.molecular_weight_threshold),
            "route_count": 1,
            "partial_route_count": 0,
            "routes": [
                {
                    "route_id": "ROUTE1:test",
                    "solved": True,
                    "route_cost": 1.2,
                    "reaction_count": 1,
                    "maximum_depth": 1,
                    "steps": [],
                    "leaves": [],
                    "warnings": [],
                }
            ],
            "partial_routes": [],
            "diagnostics": {"expanded_states": 1},
            "warnings": [],
        }

    def render_reaction(
        self, reaction_smiles: str, *, width: int, height: int
    ) -> bytes:
        assert reaction_smiles
        assert width == 760
        assert height == 220
        return b"<svg xmlns='http://www.w3.org/2000/svg'></svg>"

    def render_molecule(
        self, molecule_smiles: str, *, width: int, height: int
    ) -> bytes:
        assert molecule_smiles
        assert width == 760
        assert height == 220
        return b"<svg xmlns='http://www.w3.org/2000/svg'></svg>"


def client() -> TestClient:
    return TestClient(create_app(runtime=FakeRuntime(), frontend_dist="missing"))


def test_local_runtime_reports_isolated_full_and_compact_indexes(tmp_path) -> None:
    for mode in ("full", "compact"):
        mode_dir = tmp_path / mode
        mode_dir.mkdir()
        (mode_dir / "generic_index.sqlite").touch()

    capabilities = LocalRecommendationRuntime(library_root=tmp_path).capabilities()

    assert capabilities["default_library_mode"] == "full"
    assert capabilities["library_modes"]["full"]["index_available"] is True
    assert capabilities["library_modes"]["compact"]["index_available"] is True


def test_local_runtime_reports_retrosynthesis_library_modes(tmp_path) -> None:
    compact = tmp_path / "compact"
    compact.mkdir()
    (compact / "operator_library_v3.json.gz").touch()

    capabilities = LocalRecommendationRuntime(
        retrosynthesis_library_root=tmp_path
    ).capabilities()

    assert capabilities["retrosynthesis"] is True
    assert capabilities["default_retrosynthesis_library_mode"] == "compact"
    assert (
        capabilities["retrosynthesis_library_modes"]["compact"]["library_available"]
        is True
    )
    assert (
        capabilities["retrosynthesis_library_modes"]["full"]["library_available"]
        is False
    )


def test_local_runtime_reports_multistep_index_availability(tmp_path) -> None:
    compact = tmp_path / "operators" / "compact"
    compact.mkdir(parents=True)
    (compact / "operator_library_v3.json.gz").touch()
    literature_index = tmp_path / "literature.sqlite"
    literature_index.touch()

    capabilities = LocalRecommendationRuntime(
        retrosynthesis_library_root=tmp_path / "operators",
        literature_index_path=literature_index,
    ).capabilities()

    assert capabilities["multistep_retrosynthesis"] is True
    assert capabilities["literature_molecule_index_available"] is True
    assert capabilities["literature_molecule_index_name"] == "literature.sqlite"


def test_local_runtime_runs_multistep_planner_with_web_limits(
    monkeypatch,
    tmp_path,
) -> None:
    literature_index = tmp_path / "literature.sqlite"
    literature_index.touch()
    library = SimpleNamespace(operators=("operator",), templates=("template",))
    runtime = LocalRecommendationRuntime(
        literature_index_path=literature_index,
    )
    monkeypatch.setattr(
        runtime,
        "_get_retrosynthesis_library",
        lambda mode: library,
    )

    opened_paths = []

    class FakeMoleculeIndex:
        def __init__(self, path) -> None:
            opened_paths.append(path)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback) -> None:
            return None

    captured = {}

    def fake_plan(target_smiles, selected_library, stock_index, **options):
        captured.update(
            {
                "target_smiles": target_smiles,
                "library": selected_library,
                "stock_index": stock_index,
                **options,
            }
        )
        return SimpleNamespace(
            routes=("route",),
            partial_routes=(),
            to_dict=lambda: {
                "target_smiles": target_smiles,
                "routes": [],
                "partial_routes": [],
                "diagnostics": {},
                "warnings": [],
            },
        )

    monkeypatch.setattr(runtime_module, "CanonicalMoleculeIndex", FakeMoleculeIndex)
    monkeypatch.setattr(runtime_module, "plan_multistep_routes", fake_plan)

    payload = runtime.multistep_retrosynthesize(
        MultistepRetrosynthesisRequest(
            target_smiles=" Cc1ccccc1 ",
            library_mode="compact",
            top_k_routes=3,
            max_depth=2,
            molecular_weight_threshold=140.0,
            include_l0=False,
            use_context=True,
            diversify=False,
        )
    )

    assert opened_paths == [literature_index]
    assert captured["target_smiles"] == "Cc1ccccc1"
    assert captured["library"] is library
    assert captured["max_depth"] == 2
    assert captured["molecular_weight_threshold"] == 140.0
    assert captured["top_k_routes"] == 3
    assert captured["per_step_top_k"] == 5
    assert captured["beam_width"] == 12
    assert captured["max_expansions"] == 4
    assert captured["max_templates_to_apply"] == 40
    assert captured["max_candidates_to_validate"] == 10
    assert captured["include_l0"] is False
    assert captured["diversify"] is False
    assert payload["valid"] is True
    assert payload["route_count"] == 1


def test_local_runtime_loads_paired_catalogs_from_string_index_path(
    tmp_path,
) -> None:
    index_path = tmp_path / "generic_index.sqlite"
    index_path.touch()
    reference = {"reference_id": "REF1:paper", "raw_reference": "Example"}
    detail = {
        "observation_id": "observation:1",
        "reaction_id": "reaction:1",
        "procedure_text": "Stir for 1 h.",
    }
    with gzip.open(
        tmp_path / "reference_catalog.jsonl.gz", "wt", encoding="utf-8"
    ) as handle:
        handle.write(json.dumps(reference) + "\n")
    with gzip.open(
        tmp_path / "experimental_detail_catalog.jsonl.gz",
        "wt",
        encoding="utf-8",
    ) as handle:
        handle.write(json.dumps(detail) + "\n")

    runtime = LocalRecommendationRuntime(index_path=index_path)

    assert runtime._get_reference_catalog(str(index_path)) == {"REF1:paper": reference}
    assert runtime._get_experimental_detail_catalog(str(index_path)) == {
        "observation:observation:1": detail,
        "reaction:reaction:1": detail,
    }


def test_local_retrosynthesis_returns_hits_before_condition_lookup(
    monkeypatch, tmp_path
) -> None:
    candidate = SimpleNamespace(
        target_smiles="CCN",
        template_id="template:1",
        proposed_reaction_smiles="CCBr.N>>CCN",
        to_dict=lambda: {
            "target_smiles": "CCN",
            "template_id": "template:1",
            "precursor_smiles": "CCBr.N",
            "proposed_reaction_smiles": "CCBr.N>>CCN",
            "precedent_reaction_ids": ["internal-row-id"],
            "selectivity_warnings": [
                {
                    "code": "POSSIBLE_FUNCTIONAL_GROUP_COMPETITION",
                    "message": "Review a competing endpoint.",
                    "ranking_impact": "none",
                    "competing_outcomes": [],
                }
            ],
        },
    )
    precedent = SimpleNamespace(
        precursor_smiles="CCBr.N",
        product_smiles="CCN",
        reference_id="REF1:paper",
    )
    library = SimpleNamespace(
        operators=(object(),),
        templates=(
            SimpleNamespace(
                template_id="template:1",
                precedents=(precedent,),
            ),
        ),
    )
    reference = {
        "reference_id": "REF1:paper",
        "raw_reference": (
            "Example synthesis | A. Chemist | Journal of Examples (2025), 1, 1"
        ),
    }
    runtime = LocalRecommendationRuntime(
        library_root=tmp_path,
        retrosynthesis_library_root=tmp_path,
    )
    monkeypatch.setattr(runtime, "_get_retrosynthesis_library", lambda mode: library)
    monkeypatch.setattr(
        runtime,
        "_get_recommender",
        lambda **options: (_ for _ in ()).throw(
            AssertionError("condition lookup must not block retrosynthesis")
        ),
    )
    monkeypatch.setattr(
        runtime,
        "_get_reference_catalog",
        lambda path: {"REF1:paper": reference},
    )
    search_options = {}

    def fake_disconnect(*args, **kwargs):
        search_options.update(kwargs)
        return (candidate,)

    monkeypatch.setattr(
        runtime_module,
        "disconnect_operator_ladder",
        fake_disconnect,
    )

    payload = runtime.retrosynthesize(
        RetrosynthesisRequest(target_smiles="CCN", library_mode="compact")
    )

    result = payload["candidates"][0]
    assert payload["schema_version"] == "1.3"
    assert search_options["max_templates_to_apply"] == 100
    assert search_options["max_candidates_to_validate"] == 30
    assert "precedent_reaction_ids" not in result
    assert result["supporting_precedents"] == [
        {
            "reaction_smiles": "CCBr.N>>CCN",
            "reference_record": reference,
        }
    ]
    assert result["condition_evidence"]["status"] == "pending"
    assert result["condition_evidence"]["recommendations"] == []
    assert result["selectivity_warnings"][0]["ranking_impact"] == "none"


def test_local_retrosynthesis_condition_lookup_attaches_publication_records(
    monkeypatch, tmp_path
) -> None:
    recommendation = {
        "rank": 1,
        "recipe_id": "recipe:1",
        "score": 0.9,
        "resolved_recipe": {"bases": []},
        "precedent_reaction_ids": ["reaction:1"],
        "precedent_reference_ids": ["REF1:paper"],
    }
    evidence = SimpleNamespace(
        to_dict=lambda: {
            "status": "recommended_direct",
            "recommendations": [recommendation],
            "warnings": [],
        }
    )
    reference = {"reference_id": "REF1:paper", "raw_reference": "Example"}
    experimental = {
        "observation_id": "observation:1",
        "reaction_id": "reaction:1",
        "procedure_text": "Stir for one hour.",
    }
    recommender = SimpleNamespace(source_path=str(tmp_path / "generic_index.sqlite"))
    runtime = LocalRecommendationRuntime(library_root=tmp_path)
    monkeypatch.setattr(runtime, "_get_recommender", lambda **options: recommender)
    monkeypatch.setattr(
        runtime_module,
        "recommend_retrosynthesis_conditions",
        lambda *args, **kwargs: evidence,
    )
    monkeypatch.setattr(
        runtime,
        "_get_reference_catalog",
        lambda path: {"REF1:paper": reference},
    )
    monkeypatch.setattr(
        runtime,
        "_get_experimental_detail_catalog",
        lambda path: {"reaction:reaction:1": experimental},
    )

    payload = runtime.retrosynthesis_conditions(
        RetrosynthesisConditionsRequest(
            reaction_smiles="CCBr.N>>CCN",
            library_mode="compact",
        )
    )

    condition = payload["recommendations"][0]
    assert condition["precedent_references"] == [reference]
    assert condition["precedent_experimental_details"] == [experimental]


def test_health_and_capabilities_are_versioned() -> None:
    web = client()
    health = web.get("/api/v1/health")
    capabilities = web.get("/api/v1/capabilities")

    assert health.status_code == 200
    assert health.json()["api_schema_version"] == "1.0"
    assert capabilities.json()["data"]["index_available"] is True


def test_multistep_retrosynthesis_contract_forwards_bounded_options() -> None:
    response = client().post(
        "/api/v1/retrosynthesis/routes",
        json={
            "target_smiles": "Cc1ccccc1",
            "library_mode": "compact",
            "top_k_routes": 3,
            "max_depth": 2,
            "molecular_weight_threshold": 140.0,
            "include_l0": False,
            "use_context": True,
            "diversify": True,
        },
    )

    assert response.status_code == 200
    result = response.json()["data"]
    assert result["route_count"] == 1
    assert result["max_depth"] == 2
    assert result["molecular_weight_threshold"] == 140.0


def test_multistep_retrosynthesis_rejects_depth_above_three() -> None:
    response = client().post(
        "/api/v1/retrosynthesis/routes",
        json={
            "target_smiles": "Cc1ccccc1",
            "max_depth": 4,
        },
    )

    assert response.status_code == 422


def test_prepare_reaction_returns_typed_error() -> None:
    response = client().post(
        "/api/v1/reactions/prepare",
        json={"reaction_smiles": "invalid"},
    )

    assert response.status_code == 422
    assert response.json()["detail"] == {
        "code": "INVALID_REACTION",
        "message": "INVALID_REACTION",
    }


def test_recommendation_contract_forwards_validated_options() -> None:
    response = client().post(
        "/api/v1/recommendations",
        json={
            "reaction_smiles": "CCBr.N>>CCN",
            "library_mode": "compact",
            "top_k": 4,
            "use_rxnmapper": False,
            "ranking_preferences": {
                "profile_id": "default",
                "weights": {},
            },
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["api_schema_version"] == "1.0"
    assert payload["data"]["query_reaction_smiles"] == "CCBr.N>>CCN"
    assert payload["data"]["library_mode"] == "compact"
    assert payload["data"]["recommendations"][0]["recipe_id"] == "recipe:1"
    assert (
        payload["data"]["recommendations"][0]["synthesis_protocol"]["materials"][0][
            "cas"
        ]
        == "584-08-7"
    )


def test_discovery_and_svg_rendering_contracts() -> None:
    web = client()
    discovery = web.post(
        "/api/v1/discovery",
        json={
            "reaction_smiles": "CCBr.N>>CCN",
            "library_mode": "compact",
            "view": "closest_chemistry",
            "use_rxnmapper": False,
        },
    )
    drawing = web.post(
        "/api/v1/render/reaction",
        json={"reaction_smiles": "CCBr.N>>CCN"},
    )

    assert discovery.status_code == 200
    assert discovery.json()["data"]["library_mode"] == "compact"
    assert discovery.json()["data"]["hits"][0]["reaction_id"] == "rxn:1"
    assert drawing.status_code == 200
    assert drawing.headers["content-type"].startswith("image/svg+xml")
    assert drawing.content.startswith(b"<svg")


def test_feature_analysis_and_molecule_rendering_contracts() -> None:
    web = client()
    feature_result = web.post(
        "/api/v1/features/analyze",
        json={"input_smiles": "Brc1ccccc1", "use_rxnmapper": False},
    )
    drawing = web.post(
        "/api/v1/render/molecule",
        json={"molecule_smiles": "Brc1ccccc1"},
    )

    assert feature_result.status_code == 200
    assert feature_result.json()["data"]["input_kind"] == "molecule"
    assert feature_result.json()["data"]["overview"]["atom_count"] == 6
    assert drawing.status_code == 200
    assert drawing.headers["content-type"].startswith("image/svg+xml")
    assert drawing.content.startswith(b"<svg")


def test_retrosynthesis_contract_forwards_operator_options() -> None:
    response = client().post(
        "/api/v1/retrosynthesis",
        json={
            "target_smiles": "CCN",
            "library_mode": "compact",
            "top_k": 7,
            "include_l0": False,
            "use_context": True,
            "diversify": True,
        },
    )

    assert response.status_code == 200
    payload = response.json()["data"]
    assert payload["target_smiles"] == "CCN"
    assert payload["library_mode"] == "compact"
    assert payload["candidates"][0]["precursor_smiles"] == "CCBr.N"


def test_retrosynthesis_conditions_contract_is_independent() -> None:
    response = client().post(
        "/api/v1/retrosynthesis/conditions",
        json={
            "reaction_smiles": "CCBr.N>>CCN",
            "library_mode": "compact",
            "top_k": 3,
        },
    )

    assert response.status_code == 200
    payload = response.json()["data"]
    assert payload["status"] == "recommended_direct"
    assert payload["query_reaction_smiles"] == "CCBr.N>>CCN"


def test_retrosynthesis_contract_requires_one_target() -> None:
    response = client().post(
        "/api/v1/retrosynthesis",
        json={"target_smiles": "", "library_mode": "compact"},
    )

    assert response.status_code == 422


def test_request_contract_rejects_unknown_fields() -> None:
    response = client().post(
        "/api/v1/recommendations",
        json={
            "reaction_smiles": "CCBr.N>>CCN",
            "unknown_option": True,
        },
    )

    assert response.status_code == 422
