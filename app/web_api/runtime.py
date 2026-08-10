"""Process-local runtime that composes the standalone chemistry packages."""

from __future__ import annotations

import os
from pathlib import Path
from threading import RLock
from typing import Any, Dict, Protocol

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
from reactive_taxonomy import RxnMapperProvider
from visualization import (
    render_molecule_image_bytes,
    render_reaction_image_bytes,
)

from .contracts import (
    DiscoveryRequest,
    FeatureAnalysisRequest,
    RecommendationRequest,
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
