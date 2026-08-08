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
from .features import analyze_features, detect_input_kind


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INDEX_PATH = (
    PROJECT_ROOT / "datasets" / "literature" / "generic_index.sqlite"
)


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

    def __init__(self, index_path: str | Path | None = None) -> None:
        configured = index_path or os.environ.get("CONDITION_RECOMMENDER_INDEX")
        self.index_path = Path(configured) if configured else DEFAULT_INDEX_PATH
        self._recommenders: Dict[
            tuple[str, int, int, bool, bool], GenericConditionRecommender
        ] = {}
        self._lock = RLock()
        self._feature_mapping_provider: RxnMapperProvider | None = None

    def _cache_key(
        self, *, use_rxnmapper: bool, include_review: bool
    ) -> tuple[str, int, int, bool, bool]:
        resolved = self.index_path.resolve()
        stat = resolved.stat()
        return (
            str(resolved),
            stat.st_size,
            stat.st_mtime_ns,
            use_rxnmapper,
            include_review,
        )

    def _get_recommender(
        self, *, use_rxnmapper: bool, include_review: bool
    ) -> GenericConditionRecommender:
        if not self.index_path.is_file():
            raise FileNotFoundError(
                f"recommendation index is unavailable: {self.index_path.name}"
            )
        if use_rxnmapper and not RxnMapperProvider.is_available():
            raise RuntimeError("RXNMAPPER_UNAVAILABLE")
        key = self._cache_key(
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
                self.index_path,
                mapping_provider=(
                    RxnMapperProvider() if use_rxnmapper else None
                ),
                include_review=include_review,
            )
            self._recommenders[key] = recommender
            return recommender

    def capabilities(self) -> Dict[str, Any]:
        """Report local feature availability without exposing absolute paths."""

        return {
            "service": "reaction-condition-recommender",
            "index_name": self.index_path.name,
            "index_available": self.index_path.is_file(),
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
        return result.to_dict()

    def discover(self, request: DiscoveryRequest) -> Dict[str, Any]:
        """Execute exploratory precedent discovery over the shared index."""

        recommender = self._get_recommender(
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
        return result.to_dict()

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
