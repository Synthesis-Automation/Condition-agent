"""Process-local runtime that composes the standalone chemistry packages."""

from __future__ import annotations

import os
import time
from pathlib import Path
from threading import RLock
from typing import Any, Dict, Protocol

from cas_tools import CanonicalMoleculeIndex
from condition_recommender import (
    ChemistRankingPreferences,
    GenericConditionRecommender,
    ReactionDiscoveryExplorer,
    available_ranking_profiles,
    build_completion_selection,
    propose_reaction_completion,
    resolve_ranking_preferences,
)
from condition_recommender.reaction_completion import (
    validate_completion_selections,
)
from core_retrosynthesis_poc import (
    GenericTemplateLibrary,
    disconnect_operator_ladder,
    load_generic_library,
    plan_multistep_routes,
    recommend_retrosynthesis_conditions,
)
from reactive_taxonomy import RxnMapperProvider
from visualization import (
    render_molecule_image_bytes,
    render_reaction_image_bytes,
)

from .contracts import (
    DiscoveryRequest,
    FeatureAnalysisRequest,
    MultistepRetrosynthesisRequest,
    RecommendationRequest,
    RetrosynthesisConditionsRequest,
    RetrosynthesisRequest,
)
from .experimental_details import (
    EXPERIMENTAL_DETAIL_CATALOG_FILENAME,
    attach_discovery_experimental_details,
    attach_recommendation_experimental_details,
    load_experimental_detail_catalog,
)
from .features import analyze_features, detect_input_kind
from .references import (
    REFERENCE_CATALOG_FILENAME,
    attach_discovery_references,
    attach_recommendation_references,
    load_reference_catalog,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LIBRARY_ROOT = PROJECT_ROOT / "datasets" / "literature"
DEFAULT_INDEX_PATH = DEFAULT_LIBRARY_ROOT / "generic_index.sqlite"
DEFAULT_RETROSYNTHESIS_LIBRARY_ROOT = (
    PROJECT_ROOT
    / "results"
    / "operator_retrosynthesis_poc"
    / "full_scale_v3"
)
DEFAULT_LITERATURE_MOLECULE_INDEX = (
    PROJECT_ROOT / "results" / "literature_molecule_index.sqlite"
)
WEB_RETROSYNTHESIS_BASE_TEMPLATE_BUDGET = 100
WEB_RETROSYNTHESIS_BASE_VALIDATION_BUDGET = 30
WEB_MULTISTEP_TEMPLATE_BUDGET = 40
WEB_MULTISTEP_VALIDATION_BUDGET = 10


def _unavailable_retrosynthesis_conditions(
    reaction_smiles: str,
    error: Exception | None = None,
) -> Dict[str, Any]:
    """Preserve a disconnection when condition evidence is unavailable."""

    return {
        "status": "insufficient_evidence",
        "query_reaction_smiles": reaction_smiles,
        "recommender_valid": False,
        "recommendation_mode": "unavailable",
        "retrieval_level": None,
        "uses_type_agnostic_fallback": False,
        "candidate_count": 0,
        "independent_candidate_count": 0,
        "compatible_candidate_count": 0,
        "independent_compatible_candidate_count": 0,
        "excluded_candidate_count": 0,
        "best_recipe_score": None,
        "best_recipe_compatibility_score": None,
        "best_recipe_reference_support": 0,
        "recommendations": [],
        "warnings": ["CONDITION_RECOMMENDATION_UNAVAILABLE"],
        "error": type(error).__name__ if error is not None else None,
    }


def _pending_retrosynthesis_conditions(
    reaction_smiles: str,
) -> Dict[str, Any]:
    """Represent condition evidence that the browser will fetch progressively."""

    value = _unavailable_retrosynthesis_conditions(reaction_smiles)
    value.update(
        {
            "status": "pending",
            "recommendation_mode": "pending",
            "warnings": [],
        }
    )
    return value


class WebRuntime(Protocol):
    """Narrow application runtime consumed by the HTTP routes."""

    def capabilities(self) -> Dict[str, Any]: ...

    def ranking_profiles(self) -> tuple[Dict[str, Any], ...]: ...

    def prepare_reaction(self, reaction_smiles: str) -> Dict[str, Any]: ...

    def recommend(self, request: RecommendationRequest) -> Dict[str, Any]: ...

    def discover(self, request: DiscoveryRequest) -> Dict[str, Any]: ...

    def analyze_features(
        self, request: FeatureAnalysisRequest
    ) -> Dict[str, Any]: ...

    def retrosynthesize(
        self, request: RetrosynthesisRequest
    ) -> Dict[str, Any]: ...

    def multistep_retrosynthesize(
        self, request: MultistepRetrosynthesisRequest
    ) -> Dict[str, Any]: ...

    def retrosynthesis_conditions(
        self, request: RetrosynthesisConditionsRequest
    ) -> Dict[str, Any]: ...

    def render_reaction(
        self, reaction_smiles: str, *, width: int, height: int
    ) -> bytes: ...

    def render_molecule(
        self, molecule_smiles: str, *, width: int, height: int
    ) -> bytes: ...


class LocalRecommendationRuntime:
    """Load each immutable index/mapping configuration at most once."""

    def __init__(
        self,
        index_path: str | Path | None = None,
        *,
        library_root: str | Path | None = None,
        retrosynthesis_library_root: str | Path | None = None,
        literature_index_path: str | Path | None = None,
    ) -> None:
        configured = index_path or os.environ.get("CONDITION_RECOMMENDER_INDEX")
        self._configured_index_path = Path(configured) if configured else None
        configured_library = library_root or os.environ.get(
            "CONDITION_RECOMMENDER_LIBRARY"
        )
        if configured_library is not None:
            self.library_root = Path(configured_library)
        elif self._configured_index_path is not None:
            configured_parent = self._configured_index_path.parent
            self.library_root = (
                configured_parent.parent
                if configured_parent.name.casefold() == "full"
                else configured_parent
            )
        else:
            self.library_root = DEFAULT_LIBRARY_ROOT
        self.index_path = self._index_path("full")
        configured_retrosynthesis = (
            retrosynthesis_library_root
            or os.environ.get("CORE_RETROSYNTHESIS_LIBRARY_ROOT")
        )
        self.retrosynthesis_library_root = Path(
            configured_retrosynthesis or DEFAULT_RETROSYNTHESIS_LIBRARY_ROOT
        )
        configured_literature_index = (
            literature_index_path
            or os.environ.get("CORE_RETROSYNTHESIS_LITERATURE_INDEX")
        )
        self.literature_index_path = Path(
            configured_literature_index or DEFAULT_LITERATURE_MOLECULE_INDEX
        )
        self._recommenders: Dict[
            tuple[str, int, int, bool, bool], GenericConditionRecommender
        ] = {}
        self._lock = RLock()
        self._feature_mapping_provider: RxnMapperProvider | None = None
        self._reference_catalogs: Dict[
            tuple[str, int, int], Dict[str, Dict[str, Any]]
        ] = {}
        self._experimental_detail_catalogs: Dict[
            tuple[str, int, int], Dict[str, Dict[str, Any]]
        ] = {}
        self._retrosynthesis_libraries: Dict[
            tuple[str, int, int], GenericTemplateLibrary
        ] = {}

    def _index_path(self, library_mode: str) -> Path:
        """Resolve one isolated Full or Compact runtime index."""
        mode = library_mode.strip().casefold()
        if mode not in {"full", "compact"}:
            raise ValueError(f"unsupported library mode: {library_mode}")
        if self._configured_index_path is not None and mode == "full":
            return self._configured_index_path
        candidate = self.library_root / mode / "generic_index.sqlite"
        if mode == "full" and not candidate.is_file():
            legacy = self.library_root / "generic_index.sqlite"
            if legacy.is_file() or self._configured_index_path is None:
                return legacy
        return candidate

    def _retrosynthesis_library_path(self, library_mode: str) -> Path:
        """Resolve the immutable operator library for one web mode."""

        mode = library_mode.strip().casefold()
        if mode not in {"full", "compact"}:
            raise ValueError(f"unsupported library mode: {library_mode}")
        return (
            self.retrosynthesis_library_root
            / mode
            / "operator_library_v3.json.gz"
        )

    def _get_retrosynthesis_library(
        self,
        library_mode: str,
    ) -> GenericTemplateLibrary:
        """Load and cache one versioned executable operator library."""

        path = self._retrosynthesis_library_path(library_mode)
        if not path.is_file():
            raise FileNotFoundError(
                "retrosynthesis operator library is unavailable for "
                f"{library_mode} mode"
            )
        stat = path.stat()
        key = (str(path.resolve()), stat.st_size, stat.st_mtime_ns)
        with self._lock:
            cached = self._retrosynthesis_libraries.get(key)
            if cached is not None:
                return cached
            for old_key in tuple(self._retrosynthesis_libraries):
                if old_key[0] == key[0] and old_key[1:] != key[1:]:
                    self._retrosynthesis_libraries.pop(old_key, None)
            library = load_generic_library(path)
            self._retrosynthesis_libraries[key] = library
            return library

    def _get_reference_catalog(
        self, index_path: str | Path
    ) -> Dict[str, Dict[str, Any]]:
        """Load and cache the reference artifact paired with the active index."""

        normalized_index_path = Path(index_path)
        catalog_path = normalized_index_path.parent / REFERENCE_CATALOG_FILENAME
        if not catalog_path.is_file():
            return {}
        stat = catalog_path.stat()
        key = (str(catalog_path.resolve()), stat.st_size, stat.st_mtime_ns)
        with self._lock:
            if key not in self._reference_catalogs:
                self._reference_catalogs[key] = load_reference_catalog(
                    normalized_index_path
                )
            return self._reference_catalogs[key]

    def _get_experimental_detail_catalog(
        self, index_path: str | Path
    ) -> Dict[str, Dict[str, Any]]:
        """Load and cache observed procedures paired with the active index."""

        normalized_index_path = Path(index_path)
        catalog_path = (
            normalized_index_path.parent / EXPERIMENTAL_DETAIL_CATALOG_FILENAME
        )
        if not catalog_path.is_file():
            return {}
        stat = catalog_path.stat()
        key = (str(catalog_path.resolve()), stat.st_size, stat.st_mtime_ns)
        with self._lock:
            if key not in self._experimental_detail_catalogs:
                self._experimental_detail_catalogs[key] = (
                    load_experimental_detail_catalog(normalized_index_path)
                )
            return self._experimental_detail_catalogs[key]

    def _cache_key(
        self,
        index_path: Path,
        *,
        use_rxnmapper: bool,
        include_review: bool,
    ) -> tuple[str, int, int, bool, bool]:
        resolved = index_path.resolve()
        stat = resolved.stat()
        return (
            str(resolved),
            stat.st_size,
            stat.st_mtime_ns,
            use_rxnmapper,
            include_review,
        )

    def _get_recommender(
        self,
        *,
        library_mode: str,
        use_rxnmapper: bool,
        include_review: bool,
    ) -> GenericConditionRecommender:
        index_path = self._index_path(library_mode)
        if not index_path.is_file():
            raise FileNotFoundError(
                "recommendation index is unavailable for "
                f"{library_mode} mode: {index_path}"
            )
        if use_rxnmapper and not RxnMapperProvider.is_available():
            raise RuntimeError("RXNMAPPER_UNAVAILABLE")
        key = self._cache_key(
            index_path,
            use_rxnmapper=use_rxnmapper,
            include_review=include_review,
        )
        with self._lock:
            cached = self._recommenders.get(key)
            if cached is not None:
                return cached
            for old_key in tuple(self._recommenders):
                if old_key[0] == key[0] and old_key[1:3] != key[1:3]:
                    self._recommenders.pop(old_key, None)
            recommender = GenericConditionRecommender.from_path(
                index_path,
                mapping_provider=(
                    RxnMapperProvider() if use_rxnmapper else None
                ),
                include_review=include_review,
            )
            self._recommenders[key] = recommender
            return recommender

    def capabilities(self) -> Dict[str, Any]:
        """Report local feature availability without exposing absolute paths."""

        mode_paths = {
            mode: self._index_path(mode) for mode in ("full", "compact")
        }
        retrosynthesis_paths = {
            mode: self._retrosynthesis_library_path(mode)
            for mode in ("full", "compact")
        }
        default_retrosynthesis_mode = (
            "full" if retrosynthesis_paths["full"].is_file() else "compact"
        )
        return {
            "service": "reaction-condition-recommender",
            "index_name": mode_paths["full"].name,
            "index_available": mode_paths["full"].is_file(),
            "default_library_mode": "full",
            "library_modes": {
                mode: {
                    "label": mode.title(),
                    "index_name": path.name,
                    "index_available": path.is_file(),
                }
                for mode, path in mode_paths.items()
            },
            "loaded_runtime_variants": len(self._recommenders),
            "rxnmapper_available": RxnMapperProvider.is_available(),
            "recommendation": True,
            "discovery": True,
            "featurization": True,
            "reaction_rendering": True,
            "retrosynthesis": any(
                path.is_file() for path in retrosynthesis_paths.values()
            ),
            "multistep_retrosynthesis": (
                self.literature_index_path.is_file()
                and any(path.is_file() for path in retrosynthesis_paths.values())
            ),
            "literature_molecule_index_available": (
                self.literature_index_path.is_file()
            ),
            "literature_molecule_index_name": self.literature_index_path.name,
            "default_retrosynthesis_library_mode": (
                default_retrosynthesis_mode
            ),
            "retrosynthesis_library_modes": {
                mode: {
                    "label": mode.title(),
                    "library_name": path.name,
                    "library_available": path.is_file(),
                }
                for mode, path in retrosynthesis_paths.items()
            },
            "local_only": True,
        }

    def ranking_profiles(self) -> tuple[Dict[str, Any], ...]:
        """Expose declarative chemist ranking profiles."""

        profiles = []
        for value in available_ranking_profiles():
            payload = dict(value)
            resolved = resolve_ranking_preferences(
                ChemistRankingPreferences(profile_id=payload["profile_id"])
            )
            payload["weights"] = dict(resolved.weights)
            profiles.append(payload)
        return tuple(profiles)

    def prepare_reaction(self, reaction_smiles: str) -> Dict[str, Any]:
        """Validate a reaction and return any source-completion proposal."""

        proposal = propose_reaction_completion(reaction_smiles.strip())
        return {
            "valid": True,
            "input_reaction_smiles": proposal.query_reaction_smiles,
            "completion_proposal": proposal.to_dict(),
            "warnings": list(proposal.warnings),
        }

    @staticmethod
    def _completion_selections(
        reaction_smiles: str,
        choices: tuple[Any, ...],
    ) -> tuple[Any, ...]:
        proposal = propose_reaction_completion(reaction_smiles)
        if not proposal.requirements:
            if choices:
                raise ValueError("REACTION_COMPLETION_NOT_REQUIRED")
            return ()
        expected = {value.requirement_id for value in proposal.requirements}
        supplied = [value.requirement_id for value in choices]
        if set(supplied) != expected or len(supplied) != len(expected):
            raise ValueError("REACTION_COMPLETION_CHOICES_INCOMPLETE")
        selections = tuple(
            build_completion_selection(
                proposal,
                choice.requirement_id,
                option_id=choice.option_id,
                custom_identifier=choice.custom_identifier,
            )
            for choice in choices
        )
        validate_completion_selections(proposal, selections)
        return selections

    def recommend(self, request: RecommendationRequest) -> Dict[str, Any]:
        """Execute the canonical generic recommendation use case."""

        reaction_smiles = request.reaction_smiles.strip()
        selections = self._completion_selections(
            reaction_smiles,
            request.completion_choices,
        )
        ranking = request.ranking_preferences
        preferences = ChemistRankingPreferences(
            profile_id=ranking.profile_id,
            weights=dict(ranking.weights),
            customized=bool(ranking.weights),
        )
        recommender = self._get_recommender(
            library_mode=request.library_mode,
            use_rxnmapper=request.use_rxnmapper,
            include_review=request.unrestricted_fallback,
        )
        result = recommender.recommend(
            reaction_smiles,
            top_k=request.top_k,
            minimum_pool_size=request.minimum_pool_size,
            unrestricted_fallback=request.unrestricted_fallback,
            ranking_preferences=preferences,
            completion_selections=selections,
        )
        payload = attach_recommendation_references(
            result.to_dict(),
            self._get_reference_catalog(recommender.source_path),
        )
        return attach_recommendation_experimental_details(
            payload,
            self._get_experimental_detail_catalog(recommender.source_path),
        )

    def discover(self, request: DiscoveryRequest) -> Dict[str, Any]:
        """Execute exploratory precedent discovery over the shared index."""

        recommender = self._get_recommender(
            library_mode=request.library_mode,
            use_rxnmapper=request.use_rxnmapper,
            include_review=request.include_review,
        )
        explorer = ReactionDiscoveryExplorer(
            recommender.index,
            recommender.source_path,
            recommender.mapping_provider,
        )
        result = explorer.discover(
            request.reaction_smiles.strip(),
            top_k=request.top_k,
            view=request.view,
            include_low_yield=request.include_low_yield,
            include_unreported_outcomes=request.include_unreported_outcomes,
        )
        payload = attach_discovery_references(
            result.to_dict(),
            self._get_reference_catalog(recommender.source_path),
        )
        return attach_discovery_experimental_details(
            payload,
            self._get_experimental_detail_catalog(recommender.source_path),
        )

    def analyze_features(
        self,
        request: FeatureAnalysisRequest,
    ) -> Dict[str, Any]:
        """Auto-detect and featurize one molecule or reaction."""

        value = request.input_smiles.strip()
        mapping_provider = None
        if detect_input_kind(value) == "reaction" and request.use_rxnmapper:
            if not RxnMapperProvider.is_available():
                raise RuntimeError("RXNMAPPER_UNAVAILABLE")
            with self._lock:
                if self._feature_mapping_provider is None:
                    self._feature_mapping_provider = RxnMapperProvider()
                mapping_provider = self._feature_mapping_provider
        return analyze_features(
            value,
            mapping_provider=mapping_provider,
            force_resolved_mapping=request.force_resolved_mapping,
        )

    def retrosynthesize(
        self,
        request: RetrosynthesisRequest,
    ) -> Dict[str, Any]:
        """Apply operators and return hits before slower condition lookup."""

        library = self._get_retrosynthesis_library(request.library_mode)
        template_budget = min(
            300,
            max(
                WEB_RETROSYNTHESIS_BASE_TEMPLATE_BUDGET,
                request.top_k * 10,
            ),
        )
        validation_budget = min(
            100,
            max(
                WEB_RETROSYNTHESIS_BASE_VALIDATION_BUDGET,
                request.top_k * 2,
            ),
        )
        candidates = disconnect_operator_ladder(
            request.target_smiles.strip(),
            library,
            top_k=request.top_k,
            use_context=request.use_context,
            include_l0=request.include_l0,
            diversify=request.diversify,
            max_templates_to_apply=template_budget,
            max_candidates_to_validate=validation_budget,
        )
        index_path = self._index_path(request.library_mode)
        reference_catalog = self._get_reference_catalog(index_path)
        templates = {
            template.template_id: template for template in library.templates
        }
        serialized = []
        for rank, candidate in enumerate(candidates, start=1):
            value = {
                "rank": rank,
                **candidate.to_dict(),
            }
            condition_query = getattr(
                candidate,
                "condition_query_reaction_smiles",
                "",
            ) or candidate.proposed_reaction_smiles
            value["condition_evidence"] = _pending_retrosynthesis_conditions(
                condition_query
            )
            value.pop("precedent_reaction_ids", None)
            template = templates.get(candidate.template_id)
            supporting_precedents = []
            seen_support = set()
            for precedent in template.precedents if template is not None else ():
                reaction_smiles = (
                    f"{precedent.precursor_smiles}>>{precedent.product_smiles}"
                )
                support_key = (reaction_smiles, precedent.reference_id)
                if support_key in seen_support:
                    continue
                seen_support.add(support_key)
                reference = reference_catalog.get(precedent.reference_id)
                supporting_precedents.append(
                    {
                        "reaction_smiles": reaction_smiles,
                        "reference_record": (
                            dict(reference) if reference is not None else None
                        ),
                    }
                )
            value["supporting_precedents"] = supporting_precedents
            serialized.append(value)
        return {
            "target_smiles": (
                candidates[0].target_smiles
                if candidates
                else request.target_smiles.strip()
            ),
            "library_mode": request.library_mode,
            "valid": bool(candidates),
            "error": None if candidates else "NO_RETROSYNTHESIS_CANDIDATES",
            "schema_version": "1.3",
            "candidate_count": len(candidates),
            "library_operator_count": len(library.operators),
            "library_template_count": len(library.templates),
            "warnings": [
                "Single-step proposals are generated from source-round-tripped "
                "graph operators and still require chemist review.",
                "Interactive search uses a bounded operator-validation budget; "
                "condition evidence loads progressively for each hit.",
            ],
            "candidates": serialized,
        }

    def retrosynthesis_conditions(
        self,
        request: RetrosynthesisConditionsRequest,
    ) -> Dict[str, Any]:
        """Load condition evidence independently of disconnection generation."""

        recommender = self._get_recommender(
            library_mode=request.library_mode,
            use_rxnmapper=False,
            include_review=False,
        )
        evidence = recommend_retrosynthesis_conditions(
            request.reaction_smiles.strip(),
            recommender,
            condition_top_k=request.top_k,
        ).to_dict()
        recommendation_payload = {
            "recommendations": list(evidence.get("recommendations") or ())
        }
        attach_recommendation_references(
            recommendation_payload,
            self._get_reference_catalog(recommender.source_path),
        )
        attach_recommendation_experimental_details(
            recommendation_payload,
            self._get_experimental_detail_catalog(recommender.source_path),
        )
        evidence["recommendations"] = recommendation_payload["recommendations"]
        return evidence

    def multistep_retrosynthesize(
        self,
        request: MultistepRetrosynthesisRequest,
    ) -> Dict[str, Any]:
        """Search short routes using explicit literature/MW terminal rules."""

        library = self._get_retrosynthesis_library(request.library_mode)
        started_at = time.perf_counter()
        with CanonicalMoleculeIndex(self.literature_index_path) as stock_index:
            result = plan_multistep_routes(
                request.target_smiles.strip(),
                library,
                stock_index,
                max_depth=request.max_depth,
                molecular_weight_threshold=(
                    request.molecular_weight_threshold
                ),
                top_k_routes=request.top_k_routes,
                per_step_top_k=5,
                beam_width=max(12, request.top_k_routes * 3),
                max_expansions=max(4, request.top_k_routes),
                max_templates_to_apply=WEB_MULTISTEP_TEMPLATE_BUDGET,
                max_candidates_to_validate=(
                    WEB_MULTISTEP_VALIDATION_BUDGET
                ),
                use_context=request.use_context,
                include_l0=request.include_l0,
                diversify=request.diversify,
            )
        payload = result.to_dict()
        elapsed_seconds = round(time.perf_counter() - started_at, 3)
        payload.update(
            {
                "library_mode": request.library_mode,
                "valid": bool(result.routes),
                "error": (
                    None if result.routes else "NO_MULTISTEP_ROUTES"
                ),
                "route_count": len(result.routes),
                "partial_route_count": len(result.partial_routes),
                "library_operator_count": len(library.operators),
                "library_template_count": len(library.templates),
                "search_elapsed_seconds": elapsed_seconds,
                "search_budget": {
                    "per_step_top_k": 5,
                    "beam_width": max(12, request.top_k_routes * 3),
                    "max_expansions": max(4, request.top_k_routes),
                    "max_templates_to_apply": WEB_MULTISTEP_TEMPLATE_BUDGET,
                    "max_candidates_to_validate": (
                        WEB_MULTISTEP_VALIDATION_BUDGET
                    ),
                },
            }
        )
        return payload

    def render_reaction(
        self,
        reaction_smiles: str,
        *,
        width: int,
        height: int,
    ) -> bytes:
        """Render trusted RDKit SVG without making it chemistry evidence."""

        return render_reaction_image_bytes(
            reaction_smiles.strip(),
            size=(width, height),
            image_format="svg",
            render_preset="web_consistent",
        )

    def render_molecule(
        self,
        molecule_smiles: str,
        *,
        width: int,
        height: int,
    ) -> bytes:
        """Render a molecule with the same web drawing preset."""

        return render_molecule_image_bytes(
            molecule_smiles.strip(),
            size=(width, height),
            image_format="svg",
            render_preset="web_consistent",
        )


def error_payload(exc: Exception) -> Dict[str, str]:
    """Return a stable error shape while retaining deterministic error codes."""

    message = str(exc).strip() or type(exc).__name__
    code = (
        message
        if message.isupper() and " " not in message
        else type(exc).__name__.upper()
    )
    return {"code": code, "message": message}


__all__ = [
    "DEFAULT_INDEX_PATH",
    "LocalRecommendationRuntime",
    "WebRuntime",
    "error_payload",
]
