"""Process-local runtime that composes the standalone chemistry packages."""

from __future__ import annotations

import os
import time
from pathlib import Path
from threading import RLock
from typing import Any, Dict, Protocol

from cas_tools import (
    CanonicalMoleculeIndex,
    MoleculeIdentity,
    PrecursorEvidence,
    StockPortfolio,
    assess_precursor_components,
    molecule_identity,
)
from condition_registry.loader import load_substances
from condition_recommender import (
    ChemistRankingPreferences,
    GenericConditionRecommender,
    available_ranking_profiles,
    build_completion_selection,
    propose_reaction_completion,
    resolve_ranking_preferences,
)
from condition_recommender.reaction_completion import (
    validate_completion_selections,
)
from core_retrosynthesis import (
    GenericTemplateLibrary,
    disconnect_operator_ladder,
    load_generic_library,
    plan_multistep_routes,
    recommend_retrosynthesis_conditions,
)
from forward_synthesis import (
    ForwardOperatorLibrary,
    assess_proposed_step,
    build_forward_library,
    condition_profile_catalog,
    load_forward_library,
    predict_products,
)
from reactive_taxonomy import (
    STRATEGIC_COMPLEXITY_DEFINITION_ID,
    RxnMapperProvider,
    complex_target_requires_strategic_candidate,
)
from visualization import (
    render_molecule_image_bytes,
    render_reaction_image_bytes,
)

from .contracts import (
    FeatureAnalysisRequest,
    ForwardSynthesisRequest,
    MultistepRetrosynthesisRequest,
    RecommendationRequest,
    RetrosynthesisConditionsRequest,
    RetrosynthesisRequest,
)
from .condition_precedents import attach_condition_precedents
from .experimental_details import (
    EXPERIMENTAL_DETAIL_CATALOG_FILENAME,
    attach_recommendation_experimental_details,
    load_experimental_detail_catalog,
)
from .features import analyze_features, detect_input_kind
from .references import (
    REFERENCE_CATALOG_FILENAME,
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
DEFAULT_STOCK_PORTFOLIO = (
    PROJECT_ROOT / "cas_tools" / "data" / "stock_portfolio.sqlite"
)
WEB_RETROSYNTHESIS_BASE_TEMPLATE_BUDGET = 100
WEB_RETROSYNTHESIS_BASE_VALIDATION_BUDGET = 30
WEB_MULTISTEP_TEMPLATE_BUDGET = 40
WEB_MULTISTEP_VALIDATION_BUDGET = 10
WEB_FORWARD_AUDIT_OPERATOR_BUDGET = 40
WEB_FORWARD_AUDIT_ASSIGNMENT_BUDGET = 32
WEB_FORWARD_AUDIT_OUTCOME_BUDGET = 64


def _web_forward_audit_limits() -> Dict[str, int]:
    """Return bounded enumeration limits for interactive forward audits."""

    return {
        "max_operators_to_apply": WEB_FORWARD_AUDIT_OPERATOR_BUDGET,
        "max_assignments_per_operator": WEB_FORWARD_AUDIT_ASSIGNMENT_BUDGET,
        "max_outcomes_per_operator": WEB_FORWARD_AUDIT_OUTCOME_BUDGET,
    }


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


def _compact_forward_assessment(assessment: Any) -> Dict[str, Any]:
    """Project a full forward audit into a bounded retrosynthesis payload."""

    payload = assessment.to_dict()
    blind = payload.pop("blind_prediction")
    intended_rank = payload.get("intended_product_rank")
    candidates = list(blind.get("candidates") or ())
    payload.update(
        {
            "evaluated": True,
            "blind_prediction_summary": {
                "valid": bool(blind.get("valid")),
                "status": blind.get("status"),
                "conditions_supplied": bool(blind.get("conditions_supplied")),
                "condition_profile_supplied": bool(
                    blind.get("condition_profile_supplied")
                ),
                "candidate_count": len(candidates),
                "valid_pathway_count": int(
                    (blind.get("diagnostics") or {}).get(
                        "valid_pathway_count", 0
                    )
                ),
                "top_products": [
                    {
                        "rank": candidate.get("rank"),
                        "product_smiles": candidate.get("product_smiles"),
                        "score": candidate.get("score"),
                        "is_intended": candidate.get("rank") == intended_rank,
                    }
                    for candidate in candidates[:5]
                ],
                "warnings": list(blind.get("warnings") or ()),
            },
        }
    )
    return payload


class WebRuntime(Protocol):
    """Narrow application runtime consumed by the HTTP routes."""

    def capabilities(self) -> Dict[str, Any]: ...

    def ranking_profiles(self) -> tuple[Dict[str, Any], ...]: ...

    def prepare_reaction(self, reaction_smiles: str) -> Dict[str, Any]: ...

    def recommend(self, request: RecommendationRequest) -> Dict[str, Any]: ...

    def analyze_features(
        self, request: FeatureAnalysisRequest
    ) -> Dict[str, Any]: ...

    def forward_synthesize(
        self, request: ForwardSynthesisRequest
    ) -> Dict[str, Any]: ...

    def forward_condition_profiles(self) -> Dict[str, Any]: ...

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
        stock_portfolio_path: str | Path | None = None,
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
        configured_stock_portfolio = (
            stock_portfolio_path or os.environ.get("RETROSYNTHESIS_STOCK_PORTFOLIO")
        )
        self.stock_portfolio_path = Path(
            configured_stock_portfolio or DEFAULT_STOCK_PORTFOLIO
        )
        # An explicitly supplied legacy index is an intentional test/runtime
        # override unless a stock portfolio was also explicitly configured.
        self._prefer_stock_portfolio = bool(configured_stock_portfolio) or (
            literature_index_path is None
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
        self._forward_libraries: Dict[
            tuple[str, int, int], ForwardOperatorLibrary
        ] = {}
        self._compound_registry_identities: (
            tuple[frozenset[str], frozenset[str]] | None
        ) = None

    def _registered_compound_identities(
        self,
    ) -> tuple[frozenset[str], frozenset[str]]:
        """Return exact canonical SMILES and full InChIKeys from the registry."""

        with self._lock:
            if self._compound_registry_identities is None:
                identities = tuple(
                    molecule_identity(substance.smiles)
                    for substance in load_substances()
                    if substance.smiles
                )
                self._compound_registry_identities = (
                    frozenset(
                        identity.canonical_smiles
                        for identity in identities
                        if identity is not None
                    ),
                    frozenset(
                        identity.inchi_key
                        for identity in identities
                        if identity is not None and identity.inchi_key is not None
                    ),
                )
            return self._compound_registry_identities

    @staticmethod
    def _is_terminal_stock_match(match: Any) -> bool:
        """Return whether a stock match has current terminal-grade evidence."""

        if match is None:
            return False
        return any(
            str(record.get("terminal_eligible") or "").casefold() == "true"
            for record in match.source_records
        )

    def _precursor_realism_scorer(self):
        """Open configured evidence stores and return a scorer plus close hook."""

        stock = (
            StockPortfolio(self.stock_portfolio_path)
            if self._prefer_stock_portfolio and self.stock_portfolio_path.is_file()
            else None
        )
        literature = (
            CanonicalMoleculeIndex(self.literature_index_path)
            if self.literature_index_path.is_file()
            else None
        )
        registry_smiles, registry_inchi_keys = (
            self._registered_compound_identities()
        )

        def evidence(identity: MoleculeIdentity) -> PrecursorEvidence:
            stock_match = stock.lookup(identity) if stock is not None else None
            literature_match = (
                literature.lookup(identity) if literature is not None else None
            )
            return PrecursorEvidence(
                buyable=self._is_terminal_stock_match(stock_match),
                in_compound_registry=(
                    identity.canonical_smiles in registry_smiles
                    or (
                        identity.inchi_key is not None
                        and identity.inchi_key in registry_inchi_keys
                    )
                ),
                in_literature=(literature_match is not None),
            )

        def score(precursors: str):
            return assess_precursor_components(precursors, evidence)

        def close() -> None:
            if stock is not None:
                stock.close()
            if literature is not None:
                literature.close()

        return score, close

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

    def _forward_library_path(self, library_mode: str) -> Path:
        """Resolve an optional prebuilt forward-operator artifact."""

        mode = library_mode.strip().casefold()
        if mode not in {"full", "compact"}:
            raise ValueError(f"unsupported library mode: {library_mode}")
        return (
            self.retrosynthesis_library_root
            / mode
            / "forward_operator_library_v1.json.gz"
        )

    def _get_forward_library(self, library_mode: str) -> ForwardOperatorLibrary:
        """Load a prebuilt forward library or derive and cache it once."""

        prepared_path = self._forward_library_path(library_mode)
        source_path = (
            prepared_path
            if prepared_path.is_file()
            else self._retrosynthesis_library_path(library_mode)
        )
        if not source_path.is_file():
            raise FileNotFoundError(
                "forward operator source library is unavailable for "
                f"{library_mode} mode"
            )
        stat = source_path.stat()
        key = (str(source_path.resolve()), stat.st_size, stat.st_mtime_ns)
        with self._lock:
            cached = self._forward_libraries.get(key)
            if cached is not None:
                return cached
            for old_key in tuple(self._forward_libraries):
                if old_key[0] == key[0] and old_key[1:] != key[1:]:
                    self._forward_libraries.pop(old_key, None)
            if prepared_path.is_file():
                library = load_forward_library(prepared_path)
            else:
                library = build_forward_library(
                    self._get_retrosynthesis_library(library_mode)
                )
            self._forward_libraries[key] = library
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
        forward_paths = {
            mode: self._forward_library_path(mode)
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
            "featurization": True,
            "reaction_rendering": True,
            "forward_synthesis": any(
                path.is_file() for path in retrosynthesis_paths.values()
            ),
            "retrosynthesis": any(
                path.is_file() for path in retrosynthesis_paths.values()
            ),
            "multistep_retrosynthesis": (
                (
                    (
                        self._prefer_stock_portfolio
                        and self.stock_portfolio_path.is_file()
                    )
                    or self.literature_index_path.is_file()
                )
                and any(path.is_file() for path in retrosynthesis_paths.values())
            ),
            "literature_molecule_index_available": self.literature_index_path.is_file(),
            "literature_molecule_index_name": self.literature_index_path.name,
            "stock_portfolio_available": (
                self._prefer_stock_portfolio
                and self.stock_portfolio_path.is_file()
            ),
            "stock_portfolio_name": self.stock_portfolio_path.name,
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
            "forward_library_modes": {
                mode: {
                    "label": mode.title(),
                    "library_name": (
                        forward_paths[mode].name
                        if forward_paths[mode].is_file()
                        else path.name
                    ),
                    "library_available": path.is_file(),
                    "prepared": forward_paths[mode].is_file(),
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
        reference_catalog = self._get_reference_catalog(recommender.source_path)
        experimental_catalog = self._get_experimental_detail_catalog(
            recommender.source_path
        )
        payload = attach_recommendation_references(
            result.to_dict(),
            reference_catalog,
        )
        payload = attach_recommendation_experimental_details(
            payload,
            experimental_catalog,
        )
        return attach_condition_precedents(
            payload,
            reference_catalog,
            experimental_catalog,
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
        realism_scorer = None

        def close_realism_sources() -> None:
            return None

        if request.use_precursor_realism:
            realism_scorer, close_realism_sources = (
                self._precursor_realism_scorer()
            )
        try:
            candidates = disconnect_operator_ladder(
                request.target_smiles.strip(),
                library,
                top_k=request.top_k,
                use_context=request.use_context,
                include_l0=request.include_l0,
                diversify=request.diversify,
                max_templates_to_apply=template_budget,
                max_candidates_to_validate=validation_budget,
                precursor_realism_scorer=realism_scorer,
            )
        finally:
            close_realism_sources()
        index_path = self._index_path(request.library_mode)
        reference_catalog = self._get_reference_catalog(index_path)
        templates = {
            template.template_id: template for template in library.templates
        }
        forward_library = None
        forward_setup_warning = None
        if request.use_forward_validation:
            try:
                forward_library = self._get_forward_library(
                    request.library_mode
                )
            except (FileNotFoundError, RuntimeError, ValueError):
                forward_setup_warning = "FORWARD_VALIDATION_UNAVAILABLE"
        forward_validity_counts: Dict[str, int] = {}
        forward_partial_failure = False
        forward_levels = ("L4", "L3", "L2", "L1", "RDCHIRAL")
        if request.include_l0:
            forward_levels += ("L0",)
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
            if forward_library is not None:
                try:
                    forward_assessment = assess_proposed_step(
                        value["precursor_smiles"],
                        value["target_smiles"],
                        forward_library,
                        operator_hint=value.get("operator_id") or None,
                        levels=forward_levels,
                        top_k=max(10, min(20, request.top_k)),
                        **_web_forward_audit_limits(),
                    )
                except (RuntimeError, ValueError):
                    value["forward_assessment"] = None
                    forward_partial_failure = True
                else:
                    value["forward_assessment"] = _compact_forward_assessment(
                        forward_assessment
                    )
                    validity = forward_assessment.validity
                    forward_validity_counts[validity] = (
                        forward_validity_counts.get(validity, 0) + 1
                    )
            else:
                value["forward_assessment"] = None
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
        warnings = [
            "Single-step proposals are generated from source-round-tripped "
            "graph operators and still require chemist review.",
            "Interactive search uses a bounded operator-validation budget; "
            "condition evidence loads progressively for each hit.",
        ]
        if forward_setup_warning is not None:
            warnings.append(forward_setup_warning)
        if forward_partial_failure:
            warnings.append("FORWARD_VALIDATION_PARTIALLY_UNAVAILABLE")
        if request.use_precursor_realism:
            warnings.append(
                "Precursor realism is an experimental heuristic that uses "
                "versioned evidence bonuses and penalty bands; it is not a "
                "probability."
            )
            if not (
                self._prefer_stock_portfolio
                and self.stock_portfolio_path.is_file()
            ):
                warnings.append(
                    "No supplier stock portfolio was available; buyable "
                    "evidence was treated as absent."
                )
            if not self.literature_index_path.is_file():
                warnings.append(
                    "No literature molecule index was available; literature "
                    "evidence was treated as absent."
                )
        strategic_assessments = tuple(
            getattr(candidate, "strategic_complexity", None)
            for candidate in candidates
            if getattr(candidate, "strategic_complexity", None) is not None
        )
        if (
            strategic_assessments
            and complex_target_requires_strategic_candidate(
                strategic_assessments[0]
            )
            and not any(
                assessment.is_strategic
                for assessment in strategic_assessments
            )
        ):
            warnings.append("NO_SCAFFOLD_SIMPLIFYING_CANDIDATE_GENERATED")
        return {
            "target_smiles": (
                candidates[0].target_smiles
                if candidates
                else request.target_smiles.strip()
            ),
            "library_mode": request.library_mode,
            "valid": bool(candidates),
            "error": None if candidates else "NO_RETROSYNTHESIS_CANDIDATES",
            "schema_version": "1.8",
            "forward_validation_enabled": request.use_forward_validation,
            "forward_validity_counts": dict(
                sorted(forward_validity_counts.items())
            ),
            "precursor_realism_enabled": request.use_precursor_realism,
            "strategic_complexity_definition_id": (
                STRATEGIC_COMPLEXITY_DEFINITION_ID
            ),
            "strategic_candidate_count": sum(
                bool(getattr(candidate, "strategic_candidate", False))
                for candidate in candidates
            ),
            "precursor_realism_sources": {
                "buyable": bool(
                    self._prefer_stock_portfolio
                    and self.stock_portfolio_path.is_file()
                ),
                "compound_registry": True,
                "literature": self.literature_index_path.is_file(),
            },
            "candidate_count": len(candidates),
            "library_operator_count": len(library.operators),
            "library_template_count": len(library.templates),
            "warnings": warnings,
            "candidates": serialized,
        }

    def forward_synthesize(
        self,
        request: ForwardSynthesisRequest,
    ) -> Dict[str, Any]:
        """Predict products and optionally audit a proposed route product."""

        library = self._get_forward_library(request.library_mode)
        levels = ("L4", "L3", "L2", "L1", "RDCHIRAL")
        if request.include_l0:
            levels += ("L0",)
        starting_materials = request.starting_materials.strip()
        condition_profile = (
            request.condition_profile.model_dump()
            if request.condition_profile is not None
            else None
        )
        if request.intended_product and request.intended_product.strip():
            assessment = assess_proposed_step(
                starting_materials,
                request.intended_product.strip(),
                library,
                recipe=request.recipe,
                operator_hint=(
                    request.operator_hint.strip()
                    if request.operator_hint and request.operator_hint.strip()
                    else None
                ),
                levels=levels,
                include_self_reactions=request.include_self_reactions,
                condition_profile=condition_profile,
                top_k=request.top_k,
            )
            assessment_payload = assessment.to_dict()
            prediction = assessment_payload.pop("blind_prediction")
            return {
                "analysis_mode": "step_assessment",
                "library_mode": request.library_mode,
                "valid": bool(prediction["valid"]),
                "schema_version": assessment.schema_version,
                "forward_library_operator_count": len(library.operators),
                "prediction": prediction,
                "assessment": assessment_payload,
            }
        prediction_result = predict_products(
            starting_materials,
            library,
            recipe=request.recipe,
            levels=levels,
            include_self_reactions=request.include_self_reactions,
            condition_profile=condition_profile,
            top_k=request.top_k,
        )
        return {
            "analysis_mode": "blind_prediction",
            "library_mode": request.library_mode,
            "valid": prediction_result.valid,
            "schema_version": prediction_result.schema_version,
            "forward_library_operator_count": len(library.operators),
            "prediction": prediction_result.to_dict(),
            "assessment": None,
        }

    def forward_condition_profiles(self) -> Dict[str, Any]:
        """Expose the versioned chemist-facing condition profile catalog."""

        return condition_profile_catalog()

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
        reference_catalog = self._get_reference_catalog(recommender.source_path)
        experimental_catalog = self._get_experimental_detail_catalog(
            recommender.source_path
        )
        recommendation_payload = {
            "recommendations": list(evidence.get("recommendations") or ())
        }
        attach_recommendation_references(
            recommendation_payload,
            reference_catalog,
        )
        attach_recommendation_experimental_details(
            recommendation_payload,
            experimental_catalog,
        )
        attach_condition_precedents(
            recommendation_payload,
            reference_catalog,
            experimental_catalog,
        )
        evidence["recommendations"] = recommendation_payload["recommendations"]
        evidence["forward_assessment"] = None
        if request.use_forward_validation and evidence["recommendations"]:
            recipe = evidence["recommendations"][0].get("resolved_recipe")
            try:
                forward_library = self._get_forward_library(request.library_mode)
                levels = ("L4", "L3", "L2", "L1", "RDCHIRAL")
                if request.include_l0:
                    levels += ("L0",)
                assessment = assess_proposed_step(
                    request.starting_materials or "",
                    request.intended_product or "",
                    forward_library,
                    operator_hint=request.operator_hint or None,
                    recipe=recipe,
                    levels=levels,
                    top_k=10,
                    **_web_forward_audit_limits(),
                )
            except (FileNotFoundError, RuntimeError, ValueError):
                evidence["warnings"].append(
                    "CONDITIONED_FORWARD_VALIDATION_UNAVAILABLE"
                )
            else:
                evidence["forward_assessment"] = _compact_forward_assessment(
                    assessment
                )
        return evidence

    def multistep_retrosynthesize(
        self,
        request: MultistepRetrosynthesisRequest,
    ) -> Dict[str, Any]:
        """Search short routes using explicit stock and MW terminal rules."""

        library = self._get_retrosynthesis_library(request.library_mode)
        started_at = time.perf_counter()
        stock_path = (
            self.stock_portfolio_path
            if (
                self._prefer_stock_portfolio
                and self.stock_portfolio_path.is_file()
            )
            else self.literature_index_path
        )
        stock_type = (
            "supplier_stock_portfolio"
            if (
                self._prefer_stock_portfolio
                and self.stock_portfolio_path.is_file()
            )
            else "literature_molecule_index"
        )
        stock_context = (
            StockPortfolio(stock_path)
            if stock_type == "supplier_stock_portfolio"
            else CanonicalMoleculeIndex(stock_path)
        )
        realism_scorer = None

        def close_realism_sources() -> None:
            return None

        if request.use_precursor_realism:
            realism_scorer, close_realism_sources = (
                self._precursor_realism_scorer()
            )
        condition_recommender = None
        condition_setup_warning = None
        if request.use_condition_availability:
            try:
                condition_recommender = self._get_recommender(
                    library_mode=request.library_mode,
                    use_rxnmapper=False,
                    include_review=False,
                )
            except FileNotFoundError:
                condition_setup_warning = "CONDITION_INDEX_UNAVAILABLE"

        def condition_evidence(reaction_smiles: str):
            if condition_recommender is None:
                raise RuntimeError("condition recommender is unavailable")
            return recommend_retrosynthesis_conditions(
                reaction_smiles,
                condition_recommender,
                condition_top_k=3,
            )

        try:
            with stock_context as stock_index:
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
                    precursor_realism_scorer=realism_scorer,
                    condition_evidence_evaluator=(
                        condition_evidence
                        if condition_recommender is not None
                        else None
                    ),
                )
        finally:
            close_realism_sources()
        payload = result.to_dict()
        forward_setup_warning = None
        forward_partial_failure = False
        forward_validity_counts: Dict[str, int] = {}
        if request.use_forward_validation:
            try:
                forward_library = self._get_forward_library(
                    request.library_mode
                )
            except (FileNotFoundError, RuntimeError, ValueError):
                forward_library = None
                forward_setup_warning = "FORWARD_VALIDATION_UNAVAILABLE"
            if forward_library is not None:
                audit_cache: Dict[
                    tuple[str, str, str, str], Dict[str, Any]
                ] = {}
                forward_levels = ("L4", "L3", "L2", "L1", "RDCHIRAL")
                if request.include_l0:
                    forward_levels += ("L0",)
                for route in (*payload["routes"], *payload["partial_routes"]):
                    route_counts: Dict[str, int] = {}
                    for step in route["steps"]:
                        candidate = step["candidate"]
                        recommendations = list(
                            (step.get("condition_evidence") or {}).get(
                                "recommendations"
                            )
                            or ()
                        )
                        recipe = (
                            recommendations[0].get("resolved_recipe")
                            if recommendations
                            else None
                        )
                        recipe_key = (
                            str(recommendations[0].get("recipe_id") or recipe)
                            if recommendations
                            else ""
                        )
                        cache_key = (
                            str(candidate["precursor_smiles"]),
                            str(step["product_smiles"]),
                            str(candidate.get("operator_id") or ""),
                            recipe_key,
                        )
                        audit = audit_cache.get(cache_key)
                        if audit is None:
                            try:
                                assessment = assess_proposed_step(
                                    cache_key[0],
                                    cache_key[1],
                                    forward_library,
                                    operator_hint=cache_key[2] or None,
                                    recipe=recipe,
                                    levels=forward_levels,
                                    top_k=10,
                                    **_web_forward_audit_limits(),
                                )
                            except (RuntimeError, ValueError):
                                step["forward_assessment"] = None
                                forward_partial_failure = True
                                continue
                            audit = _compact_forward_assessment(assessment)
                            audit_cache[cache_key] = audit
                        step["forward_assessment"] = audit
                        validity = str(audit["validity"])
                        route_counts[validity] = route_counts.get(validity, 0) + 1
                        forward_validity_counts[validity] = (
                            forward_validity_counts.get(validity, 0) + 1
                        )
                    route["forward_validity_counts"] = dict(
                        sorted(route_counts.items())
                    )
            else:
                for route in (*payload["routes"], *payload["partial_routes"]):
                    for step in route["steps"]:
                        step["forward_assessment"] = None
        if condition_recommender is not None:
            reference_catalog = self._get_reference_catalog(
                condition_recommender.source_path
            )
            experimental_catalog = self._get_experimental_detail_catalog(
                condition_recommender.source_path
            )
            for route in (*payload["routes"], *payload["partial_routes"]):
                for step in route["steps"]:
                    evidence = step.get("condition_evidence")
                    if not evidence:
                        continue
                    recommendation_payload = {
                        "recommendations": list(
                            evidence.get("recommendations") or ()
                        )
                    }
                    attach_recommendation_references(
                        recommendation_payload,
                        reference_catalog,
                    )
                    attach_recommendation_experimental_details(
                        recommendation_payload,
                        experimental_catalog,
                    )
                    attach_condition_precedents(
                        recommendation_payload,
                        reference_catalog,
                        experimental_catalog,
                    )
                    evidence["recommendations"] = recommendation_payload[
                        "recommendations"
                    ]
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
                "precursor_realism_requested": request.use_precursor_realism,
                "condition_availability_requested": (
                    request.use_condition_availability
                ),
                "forward_validation_requested": request.use_forward_validation,
                "forward_validity_counts": dict(
                    sorted(forward_validity_counts.items())
                ),
                "strategic_complexity_definition_id": (
                    STRATEGIC_COMPLEXITY_DEFINITION_ID
                ),
                "precursor_realism_sources": {
                    "buyable": bool(
                        self._prefer_stock_portfolio
                        and self.stock_portfolio_path.is_file()
                    ),
                    "compound_registry": True,
                    "literature": self.literature_index_path.is_file(),
                },
                "terminal_stock_source": {
                    "type": stock_type,
                    "name": stock_path.name,
                },
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
        if condition_setup_warning is not None:
            payload["warnings"].append(condition_setup_warning)
        if forward_setup_warning is not None:
            payload["warnings"].append(forward_setup_warning)
        if forward_partial_failure:
            payload["warnings"].append(
                "FORWARD_VALIDATION_PARTIALLY_UNAVAILABLE"
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
    "DEFAULT_STOCK_PORTFOLIO",
    "LocalRecommendationRuntime",
    "WebRuntime",
    "error_payload",
]
