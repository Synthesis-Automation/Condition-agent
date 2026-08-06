"""Simple Qt6 interface for type-agnostic condition recommendation."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA_FOLDER = PROJECT_ROOT / "datasets" / "literature"
DEFAULT_SQLITE_INDEX_PATH = DEFAULT_DATA_FOLDER / "generic_index.sqlite"
DEFAULT_MANIFEST_PATH = DEFAULT_DATA_FOLDER / "shard_manifest.json"

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from condition_recommender import (  # noqa: E402
    ChemistRankingPreferences,
    GenericConditionRecommender,
    GenericRecommendationResult,
    ReactionDiscoveryExplorer,
    ReactionDiscoveryHit,
    ReactionDiscoveryResult,
    ReactionCompletionProposal,
    ReactionCompletionSelection,
    available_ranking_profiles,
    build_completion_selection,
    propose_reaction_completion,
    resolve_ranking_preferences,
)
from reactive_taxonomy import (  # noqa: E402
    RxnMapperProvider,
    render_reactivity_profile,
)
from visualization import render_reaction_image_bytes  # noqa: E402
from visualization.qt_widgets import StructureImageLabel  # noqa: E402

_RECOMMENDER_CACHE: Dict[
    Tuple[str, int, int, bool, bool],
    GenericConditionRecommender,
] = {}

QUERY_REACTION_IMAGE_SIZE = (680, 168)
PRECEDENT_REACTION_IMAGE_SIZE = (760, 240)

_RECIPE_ROLE_LABELS = (
    ("catalysts", "Catalyst"),
    ("ligands", "Ligand"),
    ("bases", "Base"),
    ("condensation_agents", "Condensation agent"),
    ("oxidants", "Oxidant"),
    ("reductants", "Reductant"),
    ("acids", "Acid"),
    ("additives", "Additive"),
    ("solvents", "Solvent"),
    ("other_components", "Other"),
)

_RANKING_COMPONENT_LABELS = {
    "similarity": "Structural similarity",
    "partner_category": "Reactant category",
    "functional_group_tolerance": "Functional-group tolerance",
    "yield": "Expected yield",
    "independent_support": "Independent support",
    "reaction_breadth": "Reaction breadth",
    "dataset_diversity": "Dataset diversity",
    "compatibility": "Condition compatibility",
    "condition_certainty": "Procedure completeness",
}


def default_recommendation_data_path() -> Path:
    """Return the fastest available default recommendation artifact."""
    if DEFAULT_SQLITE_INDEX_PATH.is_file():
        return DEFAULT_SQLITE_INDEX_PATH
    return DEFAULT_MANIFEST_PATH


def _cache_key(
    path: Path,
    use_rxnmapper: bool,
    include_review: bool,
) -> Tuple[str, int, int, bool, bool]:
    resolved = path.resolve()
    stat = resolved.stat()
    return (
        str(resolved),
        stat.st_size,
        stat.st_mtime_ns,
        use_rxnmapper,
        include_review,
    )


def _get_cached_recommender(
    path: str | Path,
    *,
    use_rxnmapper: bool = True,
    include_review: bool = False,
) -> GenericConditionRecommender:
    """Load a validated index once and invalidate it when the file changes."""
    source = Path(path)
    key = _cache_key(source, use_rxnmapper, include_review)
    recommender = _RECOMMENDER_CACHE.get(key)
    if recommender is not None:
        return recommender
    resolved = key[0]
    for old_key in tuple(_RECOMMENDER_CACHE):
        if old_key[0] == resolved and old_key[1:3] != key[1:3]:
            _RECOMMENDER_CACHE.pop(old_key, None)
    recommender = GenericConditionRecommender.from_path(
        source,
        mapping_provider=RxnMapperProvider() if use_rxnmapper else None,
        include_review=include_review,
    )
    _RECOMMENDER_CACHE[key] = recommender
    return recommender


def _get_cached_explorer(
    path: str | Path,
    *,
    use_rxnmapper: bool = True,
    include_review: bool = False,
) -> ReactionDiscoveryExplorer:
    """Share the GUI's cached validated index with the discovery API."""
    recommender = _get_cached_recommender(
        path,
        use_rxnmapper=use_rxnmapper,
        include_review=include_review,
    )
    return ReactionDiscoveryExplorer(
        recommender.index,
        recommender.source_path,
        recommender.mapping_provider,
    )


def _component_name(component: Mapping[str, Any]) -> str:
    return str(
        component.get("canonical_name")
        or component.get("raw_identifier")
        or component.get("substance_id")
        or "Unknown"
    )


def _component_names(recipe: Mapping[str, Any], field: str) -> str:
    values = recipe.get(field) or ()
    names = [
        _component_name(component)
        for component in values
        if isinstance(component, Mapping)
    ]
    return ", ".join(names)


def format_recipe_summary(recipe: Mapping[str, Any]) -> str:
    """Render the most useful recipe fields in one compact table cell."""
    parts = []
    for field, label in _RECIPE_ROLE_LABELS:
        names = _component_names(recipe, field)
        if names:
            parts.append(f"{label}: {names}")
    temperature = recipe.get("temperature_c")
    time_h = recipe.get("time_h")
    if temperature is not None:
        parts.append(f"{temperature:g} °C")
    if time_h is not None:
        parts.append(f"{time_h:g} h")
    atmosphere = str(recipe.get("atmosphere") or "").strip()
    if atmosphere:
        parts.append(f"Atmosphere: {atmosphere}")
    return "; ".join(parts) or "No resolved components"


def _display_name(value: Any) -> str:
    text = str(value or "").strip()
    return text.replace("_", " ").title() if text else "Unassigned"


def _spectator_summary(groups: Tuple[Dict[str, Any], ...]) -> str:
    """Render unchanged spectator groups with reactant and distance context."""
    spectator_labels = []
    for group in groups:
        if not isinstance(group, Mapping):
            continue
        label = str(
            group.get("chemist_label") or group.get("group_id") or "Unclassified group"
        )
        group_id = str(group.get("group_id") or "").strip()
        if group_id and group_id != label:
            label = f"{label} [{group_id}]"
        component_index = group.get("component_index")
        distance = group.get("graph_distance")
        context = []
        if isinstance(component_index, int):
            context.append(f"reactant {component_index + 1}")
        if distance is not None:
            context.append(f"d={distance}")
        spectator_labels.append(f"{label} ({', '.join(context)})" if context else label)
    return "; ".join(spectator_labels) if spectator_labels else "None detected"


def _partner_analysis_summaries(
    partners: Tuple[Dict[str, Any], ...],
) -> Tuple[str, ...]:
    """Render the canonical context-aware profile for reaction partners."""
    summaries = []
    for partner in partners:
        if not isinstance(partner, Mapping):
            continue
        role = _display_name(partner.get("role"))
        chemist_label = str(partner.get("chemist_label") or "").strip()
        identity = chemist_label or "Unclassified partner"
        summaries.append(
            f"{role} — {identity}: "
            f"{render_reactivity_profile(partner.get('reactivity_profile'))}"
        )
    return tuple(summaries)


def format_query_summary(result: GenericRecommendationResult) -> str:
    """Render query identity, spectators, and local reaction-partner context."""
    reaction_label = result.reaction_label or {}
    label_evidence = _display_name(str(reaction_label.get("status") or "unavailable"))
    lines = [
        (
            f"Reaction: {reaction_label.get('text') or 'Unresolved'}  •  "
            f"Label evidence: {label_evidence}"
        ),
        (
            f"Family: {_display_name(result.named_family)}  •  "
            f"Transformation: {_display_name(result.transformation_class)}  •  "
            f"Mode: {_display_name(result.recommendation_mode)}  •  "
            f"Retrieval: {_display_name(result.retrieval_level)}"
        ),
    ]
    preferences = result.ranking_preferences or {}
    proposal = result.completion_proposal or {}
    requirements = tuple(proposal.get("requirements") or ())
    if requirements:
        fragments = ", ".join(
            str(value.get("rooted_fragment_smiles") or value.get("fragment_key") or "")
            for value in requirements
            if isinstance(value, Mapping)
        )
        lines.append(f"Missing source requirement: {fragments}")
    if result.completion_selections:
        selections = ", ".join(
            f"{value.get('display_name') or 'Unresolved'} "
            f"[{_display_name(value.get('provenance'))}]"
            for value in result.completion_selections
            if isinstance(value, Mapping)
        )
        lines.append(f"Source completion: {selections}")
    lines.append(
        "Ranking profile: "
        f"{_display_name(preferences.get('profile_id'))}"
        + (" (customized)" if preferences.get("customized") else "")
    )
    if result.query_edit_hypothesis_ids:
        lines.append(
            "Ambiguous edit alternatives (all must agree): "
            + ", ".join(result.query_edit_hypothesis_ids)
        )
    if result.external_mapping_status:
        confidence = (
            f"{result.external_mapping_confidence:.3f}"
            if result.external_mapping_confidence is not None
            else "unavailable"
        )
        lines.append(
            "External mapping: "
            f"{_display_name(result.external_mapping_status)} • "
            f"provider {result.external_mapping_provider or 'unavailable'} • "
            f"confidence {confidence}"
        )
    lines.append(f"Spectator groups: {_spectator_summary(result.spectator_groups)}")
    partner_summaries = _partner_analysis_summaries(result.reaction_partners)
    if partner_summaries:
        lines.append("Reactivity profile:")
        lines.extend(f"  {summary}" for summary in partner_summaries)
    warnings = ", ".join(result.warnings) if result.warnings else "None"
    lines.append(
        f"Candidates: {result.candidate_count}  •  "
        f"Compatible: {result.compatible_candidate_count}  •  "
        f"Independent: {result.independent_compatible_candidate_count}  •  "
        f"Excluded: {result.excluded_candidate_count}  •  "
        f"Warnings: {warnings}"
    )
    return "\n".join(lines)


def format_discovery_summary(result: ReactionDiscoveryResult) -> str:
    """Render the query evidence and exploratory-result scope."""
    reaction_label = result.reaction_label or {}
    lines = [
        (
            f"Reaction: {reaction_label.get('text') or 'Unresolved'}  •  "
            f"Label evidence: {_display_name(reaction_label.get('status'))}"
        ),
        (
            f"Transformation: {_display_name(result.transformation_class)}  •  "
            f"Discovery view: {_display_name(result.discovery_view)}"
        ),
        "Discovery ranks structural analogues, not condition recipes. "
        "All displayed conditions are observations from individual precedents.",
    ]
    if result.query_edit_hypothesis_ids:
        lines.append(
            "Ambiguous edit hypotheses (searched separately): "
            + ", ".join(result.query_edit_hypothesis_ids)
        )
    lines.append(f"Spectator groups: {_spectator_summary(result.spectator_groups)}")
    partner_summaries = _partner_analysis_summaries(result.reaction_partners)
    if partner_summaries:
        lines.append("Reactivity profile:")
        lines.extend(f"  {summary}" for summary in partner_summaries)
    relation_text = (
        ", ".join(
            f"{_display_name(key)} {value}"
            for key, value in sorted(result.relation_counts.items())
        )
        or "None"
    )
    lines.append(
        f"Candidates: {result.candidate_count}  •  Relations: {relation_text}  •  "
        f"Warnings: {', '.join(result.warnings) if result.warnings else 'None'}"
    )
    return "\n".join(lines)


def _format_discovery_factor(value: Optional[float]) -> str:
    return "—" if value is None else f"{value:.3f}"


def _friendly_error(error: Any) -> str:
    code = str(error or "RECOMMENDATION_FAILED")
    messages = {
        "EMPTY_GENERIC_INDEX": "The selected recommendation index is empty.",
        "QUERY_HAS_NO_USABLE_REACTION_SIGNATURE": (
            "The molecular graphs did not provide enough evidence for a "
            "recommendable reaction transformation."
        ),
        "NO_CHEMICALLY_COMPATIBLE_PRECEDENT": (
            "No chemically compatible precedent was found in this dataset."
        ),
        "NO_COMPATIBLE_CONDITION_PRECEDENT": (
            "Matching reactions were found, but their condition recipes were "
            "not compatible with the query."
        ),
        "INCOMPATIBLE_REACTION_SIGNATURE_SCHEMA": (
            "The query and saved index use incompatible signature versions. "
            "Rebuild the converted data."
        ),
        "INCOMPATIBLE_REACTION_TAXONOMY_DEFINITIONS": (
            "The saved index uses older chemistry definitions. Rebuild it."
        ),
    }
    return messages.get(code, code.replace("_", " ").title())


ReactionImageLabel = StructureImageLabel


class RankingPreferencesDialog(QtWidgets.QDialog):
    """Edit transparent ranking priorities without changing chemistry gates."""

    def __init__(
        self,
        weights: Mapping[str, float],
        *,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Customize ranking priorities")
        self.setModal(True)
        self.spins: Dict[str, QtWidgets.QDoubleSpinBox] = {}
        layout = QtWidgets.QVBoxLayout(self)
        note = QtWidgets.QLabel(
            "Priorities are normalized automatically. Structural validity, "
            "precedent admission, and hard compatibility gates remain locked."
        )
        note.setWordWrap(True)
        layout.addWidget(note)
        form = QtWidgets.QFormLayout()
        for component, label in _RANKING_COMPONENT_LABELS.items():
            spin = QtWidgets.QDoubleSpinBox()
            spin.setObjectName(f"rankingWeight_{component}")
            spin.setRange(0.0, 100.0)
            spin.setDecimals(1)
            spin.setSingleStep(1.0)
            spin.setValue(100.0 * float(weights.get(component, 0.0)))
            spin.setSuffix(" %")
            form.addRow(f"{label}:", spin)
            self.spins[component] = spin
        layout.addLayout(form)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept_if_valid)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    @property
    def weights(self) -> Dict[str, float]:
        return {name: spin.value() for name, spin in self.spins.items()}

    def _accept_if_valid(self) -> None:
        if sum(self.weights.values()) <= 0.0:
            QtWidgets.QMessageBox.warning(
                self,
                "Ranking priorities required",
                "At least one ranking priority must be greater than zero.",
            )
            return
        self.accept()


class CompletionSourceDialog(QtWidgets.QDialog):
    """Collect explicit user decisions for proposed missing condition sources."""

    def __init__(
        self,
        proposal: ReactionCompletionProposal,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.proposal = proposal
        self.setWindowTitle("Confirm missing reaction source")
        self.setModal(True)
        self.source_combos: Dict[str, QtWidgets.QComboBox] = {}

        layout = QtWidgets.QVBoxLayout(self)
        note = QtWidgets.QLabel(
            "The product contains atoms not present in the submitted reactants. "
            "Confirm a condition source for recommendation filtering. The "
            "reaction SMILES will remain unchanged."
        )
        note.setWordWrap(True)
        layout.addWidget(note)
        form = QtWidgets.QFormLayout()
        for requirement in proposal.requirements:
            combo = QtWidgets.QComboBox()
            combo.setObjectName(
                f"completionSource_{requirement.requirement_id}"
            )
            combo.setEditable(True)
            combo.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
            for option in requirement.options:
                combo.addItem(option.display_name, option.option_id)
            combo.setToolTip(
                "Choose a curated source class or substance, type another "
                "identifier, or leave the source unresolved."
            )
            form.addRow(
                f"Installed fragment {requirement.rooted_fragment_smiles}:",
                combo,
            )
            self.source_combos[requirement.requirement_id] = combo
        layout.addLayout(form)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setText(
            "Use selection"
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    @property
    def selections(self) -> Tuple[ReactionCompletionSelection, ...]:
        selections = []
        for requirement in self.proposal.requirements:
            combo = self.source_combos[requirement.requirement_id]
            current_text = combo.currentText().strip()
            current_index = combo.currentIndex()
            selected_text = (
                combo.itemText(current_index).strip()
                if current_index >= 0
                else ""
            )
            if current_index >= 0 and current_text == selected_text:
                selections.append(
                    build_completion_selection(
                        self.proposal,
                        requirement.requirement_id,
                        option_id=str(combo.currentData()),
                    )
                )
            else:
                selections.append(
                    build_completion_selection(
                        self.proposal,
                        requirement.requirement_id,
                        custom_identifier=current_text,
                    )
                )
        return tuple(selections)


class GenericRecommendationWorker(QtCore.QObject):
    """Load the index and recommend outside the Qt event loop."""

    progress = QtCore.pyqtSignal(str)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(
        self,
        data_path: str,
        reaction_smiles: str,
        *,
        top_k: int,
        minimum_pool_size: Optional[int],
        unrestricted_fallback: bool = False,
        use_rxnmapper: bool = True,
        ranking_preferences: ChemistRankingPreferences | None = None,
        completion_selections: Tuple[ReactionCompletionSelection, ...] = (),
    ) -> None:
        super().__init__()
        self.data_path = data_path
        self.reaction_smiles = reaction_smiles
        self.top_k = top_k
        self.minimum_pool_size = minimum_pool_size
        self.unrestricted_fallback = unrestricted_fallback
        self.use_rxnmapper = use_rxnmapper
        self.ranking_preferences = ranking_preferences
        self.completion_selections = completion_selections

    @QtCore.pyqtSlot()
    def run(self) -> None:
        """Load or reuse the index and execute one recommendation."""
        try:
            self.progress.emit("Loading recommendation index…")
            recommender = _get_cached_recommender(
                self.data_path,
                use_rxnmapper=self.use_rxnmapper,
                include_review=self.unrestricted_fallback,
            )
            self.progress.emit(
                "Analyzing reaction with RXNMapper and ranking conditions…"
                if self.use_rxnmapper
                else "Analyzing reaction and ranking conditions…"
            )
            result = recommender.recommend(
                self.reaction_smiles,
                top_k=self.top_k,
                minimum_pool_size=self.minimum_pool_size,
                unrestricted_fallback=self.unrestricted_fallback,
                ranking_preferences=self.ranking_preferences,
                completion_selections=self.completion_selections,
            )
        except Exception as exc:
            self.finished.emit(
                False,
                None,
                f"{type(exc).__name__}: {exc}",
            )
        else:
            self.finished.emit(True, result, "")


class ReactionDiscoveryWorker(QtCore.QObject):
    """Find exploratory analogues outside the Qt event loop."""

    progress = QtCore.pyqtSignal(str)
    finished = QtCore.pyqtSignal(bool, object, str)

    def __init__(
        self,
        data_path: str,
        reaction_smiles: str,
        *,
        top_k: int,
        view: str,
        include_low_yield: bool,
        include_unreported_outcomes: bool,
        use_rxnmapper: bool = True,
        include_review: bool = False,
    ) -> None:
        super().__init__()
        self.data_path = data_path
        self.reaction_smiles = reaction_smiles
        self.top_k = top_k
        self.view = view
        self.include_low_yield = include_low_yield
        self.include_unreported_outcomes = include_unreported_outcomes
        self.use_rxnmapper = use_rxnmapper
        self.include_review = include_review

    @QtCore.pyqtSlot()
    def run(self) -> None:
        try:
            self.progress.emit("Loading shared reaction index…")
            explorer = _get_cached_explorer(
                self.data_path,
                use_rxnmapper=self.use_rxnmapper,
                include_review=self.include_review,
            )
            self.progress.emit("Comparing bond edits and local structures…")
            result = explorer.discover(
                self.reaction_smiles,
                top_k=self.top_k,
                view=self.view,
                include_low_yield=self.include_low_yield,
                include_unreported_outcomes=self.include_unreported_outcomes,
            )
        except Exception as exc:
            self.finished.emit(False, None, f"{type(exc).__name__}: {exc}")
        else:
            self.finished.emit(True, result, "")


class GenericRecommenderWindow(QtWidgets.QWidget):
    """Simple desktop interface for the clean generic recommender."""

    def __init__(self) -> None:
        super().__init__()
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Reaction Condition Recommender")
        self.resize(1180, 760)
        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[QtCore.QObject] = None
        self.last_result: Optional[
            GenericRecommendationResult | ReactionDiscoveryResult
        ] = None

        self.data_path_edit = QtWidgets.QLineEdit(
            str(default_recommendation_data_path())
        )
        self.data_path_edit.setObjectName("recommendationDataPath")
        self.data_path_edit.setPlaceholderText(
            "generic_index.sqlite or shard_manifest.json; "
            "review mode uses the paired review-core index"
        )
        self.data_summary = QtWidgets.QLabel()
        self.data_summary.setObjectName("dataSummary")
        self.data_summary.setStyleSheet("color: #52606d;")

        self.reaction_edit = QtWidgets.QLineEdit()
        self.reaction_edit.setObjectName("reactionSmiles")
        self.reaction_edit.setPlaceholderText(
            "Enter reaction SMILES with product, for example: "
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
        )

        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.setObjectName("analysisMode")
        self.mode_combo.addItem("Condition recommendation", "recommendation")
        self.mode_combo.addItem("Reaction discovery", "discovery")
        self.mode_combo.setToolTip(
            "Recommendation ranks compatible recipes. Discovery searches "
            "structural analogues and shows their observed conditions."
        )
        self.discovery_view_combo = QtWidgets.QComboBox()
        self.discovery_view_combo.setObjectName("discoveryView")
        for label, value in (
            ("Closest chemistry", "closest_chemistry"),
            ("Diverse strategies", "diverse_strategies"),
            ("Successful precedents", "successful_precedents"),
            ("Failure-informed", "failure_informed"),
        ):
            self.discovery_view_combo.addItem(label, value)
        self.include_low_yield_check = QtWidgets.QCheckBox(
            "Include low-yield precedents"
        )
        self.include_low_yield_check.setObjectName("includeLowYield")
        self.include_low_yield_check.setChecked(True)
        self.include_unreported_check = QtWidgets.QCheckBox(
            "Include unreported outcomes"
        )
        self.include_unreported_check.setObjectName("includeUnreportedOutcomes")
        self.include_unreported_check.setChecked(True)

        self.top_k_spin = QtWidgets.QSpinBox()
        self.top_k_spin.setObjectName("topK")
        self.top_k_spin.setRange(1, 50)
        self.top_k_spin.setValue(5)
        self.minimum_pool_spin = QtWidgets.QSpinBox()
        self.minimum_pool_spin.setObjectName("minimumPoolSize")
        self.minimum_pool_spin.setRange(0, 100)
        self.minimum_pool_spin.setSpecialValueText("Default")
        self.minimum_pool_spin.setValue(0)
        self.minimum_pool_spin.setToolTip(
            "Default uses the versioned retrieval definition. A larger value "
            "may force a broader chemistry fallback."
        )
        self.unrestricted_fallback_check = QtWidgets.QCheckBox(
            "Review-core and unrestricted fallback mode (expert use)"
        )
        self.unrestricted_fallback_check.setObjectName("unrestrictedFallback")
        self.unrestricted_fallback_check.setChecked(False)
        self.unrestricted_fallback_check.setToolTip(
            "Also load qualified review-core precedents. For unresolved "
            "fallbacks, bypass eligibility, similarity, independent-support, "
            "and condition-compatibility gates. Expert review is required."
        )
        self.ranking_profile_combo = QtWidgets.QComboBox()
        self.ranking_profile_combo.setObjectName("rankingProfile")
        self.ranking_profiles = {
            value["profile_id"]: value for value in available_ranking_profiles()
        }
        for value in self.ranking_profiles.values():
            self.ranking_profile_combo.addItem(value["label"], value["profile_id"])
        self.ranking_profile_combo.setCurrentIndex(
            max(0, self.ranking_profile_combo.findData("default"))
        )
        self.ranking_profile_combo.setToolTip(
            self.ranking_profiles["default"]["description"]
        )
        self.customize_ranking_button = QtWidgets.QPushButton("Customize priorities…")
        self.customize_ranking_button.setObjectName("customizeRanking")
        self.ranking_profile_status = QtWidgets.QLabel("Preset weights")
        self.ranking_profile_status.setObjectName("rankingProfileStatus")
        self.ranking_profile_status.setStyleSheet("color: #52606d;")
        self.custom_ranking_weights: Optional[Dict[str, float]] = None
        self.use_rxnmapper_check = QtWidgets.QCheckBox(
            "Use RXNMapper for unresolved or ambiguous queries"
        )
        self.use_rxnmapper_check.setObjectName("useRxnMapper")
        self.use_rxnmapper_check.setChecked(True)
        self.use_rxnmapper_check.setToolTip(
            "Checked by default. Generated correspondence is clearly reported "
            "and recommendations remain expert-review required."
        )

        self.run_button = QtWidgets.QPushButton("Recommend Conditions")
        self.run_button.setObjectName("recommendButton")
        self.run_button.setMinimumHeight(34)
        self.run_button.setStyleSheet(
            "QPushButton { background: #2b6cb0; color: white; "
            "font-weight: 600; border-radius: 4px; padding: 5px 12px; }"
            "QPushButton:hover { background: #2c5282; }"
            "QPushButton:disabled { background: #9bb7d6; }"
        )
        self.example_button = QtWidgets.QPushButton("Load Example")
        self.clear_button = QtWidgets.QPushButton("Clear")
        self.export_button = QtWidgets.QPushButton("Export JSON")
        self.export_button.setObjectName("exportButton")
        self.export_button.setEnabled(False)

        self.summary_box = QtWidgets.QPlainTextEdit()
        self.summary_box.setObjectName("recommendationSummary")
        self.summary_box.setReadOnly(True)
        self.summary_box.setFixedHeight(QUERY_REACTION_IMAGE_SIZE[1])
        self.summary_box.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        self.reaction_image_label = ReactionImageLabel(
            object_name="queryReactionGraph",
            minimum_height=QUERY_REACTION_IMAGE_SIZE[1],
        )
        self.reaction_image_label.setFixedHeight(QUERY_REACTION_IMAGE_SIZE[1])

        self.results_table = QtWidgets.QTableWidget()
        self.results_table.setObjectName("recommendationTable")
        self.results_table.setColumnCount(10)
        self.results_table.setHorizontalHeaderLabels(
            [
                "Rank",
                "Default rank",
                "Score",
                "Similarity",
                "Compatibility",
                "Expected yield",
                "Rxn support",
                "Ref support",
                "Conditions",
                "Cautions",
            ]
        )
        self.results_table.setAlternatingRowColors(True)
        self.results_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.results_table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.results_table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.SingleSelection
        )
        self.results_table.setSortingEnabled(True)
        header = self.results_table.horizontalHeader()
        header.setSectionResizeMode(
            8,
            QtWidgets.QHeaderView.ResizeMode.Stretch,
        )

        self.details_box = QtWidgets.QPlainTextEdit()
        self.details_box.setObjectName("recommendationDetails")
        self.details_box.setReadOnly(True)
        self.selected_reaction_image_label = ReactionImageLabel(
            placeholder="Select a recipe to view its first precedent reaction.",
            object_name="selectedPrecedentReactionGraph",
            minimum_height=180,
        )

        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setObjectName("recommendationProgress")
        self.progress_bar.setVisible(False)
        self.status_label = QtWidgets.QLabel("Ready")
        self.status_label.setObjectName("recommendationStatus")

        self._build_layout()
        self._bind_signals()
        self._update_data_summary()
        self._mode_changed(0)
        QtGui.QShortcut(
            QtGui.QKeySequence("Ctrl+Return"),
            self,
            activated=self.start_recommendation,
        )

    def _build_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(8)

        self.data_row_layout = QtWidgets.QHBoxLayout()
        self.data_label = QtWidgets.QLabel("Recommendation data")
        self.data_label.setObjectName("recommendationDataLabel")
        self.data_row_layout.addWidget(self.data_label)
        self.data_row_layout.addWidget(self.data_path_edit, stretch=1)
        browse = QtWidgets.QPushButton("Browse…")
        browse.clicked.connect(self.choose_data_path)
        self.data_row_layout.addWidget(browse)
        self.data_row_layout.addWidget(self.data_summary)
        layout.addLayout(self.data_row_layout)

        form = QtWidgets.QFormLayout()
        self.reaction_row_label = QtWidgets.QWidget()
        reaction_row_label_layout = QtWidgets.QHBoxLayout(self.reaction_row_label)
        reaction_row_label_layout.setContentsMargins(0, 0, 0, 0)
        reaction_row_label_layout.setSpacing(8)
        reaction_row_label_layout.addWidget(self.mode_combo)
        reaction_row_label_layout.addWidget(QtWidgets.QLabel("Reaction SMILES:"))
        form.addRow(self.reaction_row_label, self.reaction_edit)
        options = QtWidgets.QHBoxLayout()
        options.addWidget(QtWidgets.QLabel("Top results"))
        options.addWidget(self.top_k_spin)
        options.addSpacing(18)
        options.addWidget(QtWidgets.QLabel("Minimum precedent pool"))
        options.addWidget(self.minimum_pool_spin)
        options.addSpacing(18)
        options.addWidget(self.unrestricted_fallback_check)
        options.addSpacing(18)
        options.addWidget(self.use_rxnmapper_check)
        options.addStretch()
        form.addRow("Options:", options)
        ranking_options = QtWidgets.QHBoxLayout()
        ranking_options.addWidget(self.ranking_profile_combo)
        ranking_options.addWidget(self.customize_ranking_button)
        ranking_options.addWidget(self.ranking_profile_status)
        ranking_options.addStretch()
        form.addRow("Ranking priorities:", ranking_options)
        self.discovery_options_widget = QtWidgets.QWidget()
        discovery_options = QtWidgets.QHBoxLayout(self.discovery_options_widget)
        discovery_options.setContentsMargins(0, 0, 0, 0)
        discovery_options.addWidget(self.discovery_view_combo)
        discovery_options.addWidget(self.include_low_yield_check)
        discovery_options.addWidget(self.include_unreported_check)
        discovery_options.addStretch()
        self.discovery_options_label = QtWidgets.QLabel("Discovery options:")
        form.addRow(
            self.discovery_options_label,
            self.discovery_options_widget,
        )
        layout.addLayout(form)

        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.example_button)
        buttons.addWidget(self.clear_button)
        buttons.addWidget(self.export_button)
        buttons.addStretch()
        layout.addLayout(buttons)

        splitter = QtWidgets.QSplitter(
            QtCore.Qt.Orientation.Vertical,
        )
        upper = QtWidgets.QWidget()
        upper_layout = QtWidgets.QVBoxLayout(upper)
        upper_layout.setContentsMargins(0, 0, 0, 0)
        self.query_summary_panel = QtWidgets.QWidget()
        self.query_summary_panel.setObjectName("querySummaryPanel")
        self.query_summary_layout = QtWidgets.QHBoxLayout(self.query_summary_panel)
        self.query_summary_layout.setContentsMargins(0, 0, 0, 0)
        self.query_summary_layout.setSpacing(12)

        text_column = QtWidgets.QWidget()
        text_layout = QtWidgets.QVBoxLayout(text_column)
        text_layout.setContentsMargins(0, 0, 0, 0)
        text_layout.setSpacing(4)
        text_layout.addWidget(QtWidgets.QLabel("Query summary"))
        text_layout.addWidget(self.summary_box)

        graph_column = QtWidgets.QWidget()
        graph_layout = QtWidgets.QVBoxLayout(graph_column)
        graph_layout.setContentsMargins(0, 0, 0, 0)
        graph_layout.setSpacing(4)
        graph_layout.addWidget(QtWidgets.QLabel("Reaction graph"))
        graph_layout.addWidget(self.reaction_image_label)

        self.query_summary_layout.addWidget(text_column, stretch=1)
        self.query_summary_layout.addWidget(graph_column, stretch=1)
        upper_layout.addWidget(self.query_summary_panel)
        self.results_heading = QtWidgets.QLabel("Recommended recipes")
        upper_layout.addWidget(self.results_heading)
        upper_layout.addWidget(self.results_table)
        splitter.addWidget(upper)

        lower = QtWidgets.QWidget()
        lower_layout = QtWidgets.QVBoxLayout(lower)
        lower_layout.setContentsMargins(0, 0, 0, 0)
        self.selected_details_layout = QtWidgets.QHBoxLayout()
        self.selected_details_layout.setContentsMargins(0, 0, 0, 0)
        self.selected_details_layout.setSpacing(12)

        details_column = QtWidgets.QWidget()
        details_layout = QtWidgets.QVBoxLayout(details_column)
        details_layout.setContentsMargins(0, 0, 0, 0)
        details_layout.setSpacing(4)
        self.details_heading = QtWidgets.QLabel("Selected recipe details")
        details_layout.addWidget(self.details_heading)
        details_layout.addWidget(self.details_box)

        precedent_column = QtWidgets.QWidget()
        precedent_layout = QtWidgets.QVBoxLayout(precedent_column)
        precedent_layout.setContentsMargins(0, 0, 0, 0)
        precedent_layout.setSpacing(4)
        self.precedent_heading = QtWidgets.QLabel("First precedent reaction")
        precedent_layout.addWidget(self.precedent_heading)
        precedent_layout.addWidget(self.selected_reaction_image_label)

        self.selected_details_layout.addWidget(details_column, stretch=1)
        self.selected_details_layout.addWidget(precedent_column, stretch=1)
        lower_layout.addLayout(self.selected_details_layout)
        splitter.addWidget(lower)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        layout.addWidget(splitter, stretch=1)

        status_row = QtWidgets.QHBoxLayout()
        self.progress_bar.setMaximumWidth(220)
        status_row.addWidget(self.progress_bar)
        status_row.addWidget(self.status_label)
        status_row.addStretch()
        layout.addLayout(status_row)

    def _bind_signals(self) -> None:
        self.data_path_edit.editingFinished.connect(self._update_data_summary)
        self.run_button.clicked.connect(self.start_recommendation)
        self.example_button.clicked.connect(self.load_example)
        self.clear_button.clicked.connect(self.clear_results)
        self.export_button.clicked.connect(self.export_json)
        self.ranking_profile_combo.currentIndexChanged.connect(
            self._ranking_profile_changed
        )
        self.mode_combo.currentIndexChanged.connect(self._mode_changed)
        self.customize_ranking_button.clicked.connect(self.customize_ranking_priorities)
        self.results_table.itemSelectionChanged.connect(
            self._synchronize_results_current_cell
        )
        self.results_table.itemSelectionChanged.connect(self._show_selected_details)

    @QtCore.pyqtSlot(int)
    def _mode_changed(self, _index: int) -> None:
        """Switch presentation and controls without changing backend contracts."""
        discovery = self.mode_combo.currentData() == "discovery"
        self.discovery_options_label.setVisible(discovery)
        self.discovery_options_widget.setVisible(discovery)
        self.ranking_profile_combo.setEnabled(not discovery)
        self.customize_ranking_button.setEnabled(not discovery)
        self.minimum_pool_spin.setEnabled(not discovery)
        self.run_button.setText(
            "Discover Related Reactions" if discovery else "Recommend Conditions"
        )
        self.results_heading.setText(
            "Related reaction precedents" if discovery else "Recommended recipes"
        )
        self.details_heading.setText(
            "Selected precedent evidence" if discovery else "Selected recipe details"
        )
        self.precedent_heading.setText(
            "Precedent reaction" if discovery else "First precedent reaction"
        )
        headers = (
            [
                "Rank",
                "Relation",
                "Score",
                "Edit",
                "Center",
                "Environment",
                "Reactant category",
                "Yield",
                "Observed conditions",
                "Cautions",
            ]
            if discovery
            else [
                "Rank",
                "Default rank",
                "Score",
                "Similarity",
                "Compatibility",
                "Expected yield",
                "Rxn support",
                "Ref support",
                "Conditions",
                "Cautions",
            ]
        )
        self.results_table.setHorizontalHeaderLabels(headers)
        self.results_table.horizontalHeader().setSectionResizeMode(
            8,
            QtWidgets.QHeaderView.ResizeMode.Stretch,
        )
        self.clear_results()

    @QtCore.pyqtSlot()
    def _synchronize_results_current_cell(self) -> None:
        """Keep the current-cell highlight inside the selected recipe row.

        Selecting a row through its vertical header changes the selection but
        some Qt styles retain and paint the previous current cell. Moving the
        current index without changing selection prevents that stale blue cell.
        """
        selection_model = self.results_table.selectionModel()
        selected_rows = selection_model.selectedRows()
        if len(selected_rows) != 1:
            return
        selected_row = selected_rows[0].row()
        current = self.results_table.currentIndex()
        if current.isValid() and current.row() == selected_row:
            return
        self.results_table.setCurrentCell(
            selected_row,
            0,
            QtCore.QItemSelectionModel.SelectionFlag.NoUpdate,
        )

    @QtCore.pyqtSlot()
    def choose_data_path(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Choose recommendation index or manifest",
            self.data_path_edit.text() or str(DEFAULT_DATA_FOLDER),
            ("Recommendation data (*.sqlite shard_manifest.json);;All files (*)"),
        )
        if path:
            self.data_path_edit.setText(path)
            self._update_data_summary()

    @QtCore.pyqtSlot(int)
    def _ranking_profile_changed(self, _index: int) -> None:
        profile_id = str(self.ranking_profile_combo.currentData() or "default")
        profile = self.ranking_profiles.get(profile_id, {})
        self.ranking_profile_combo.setToolTip(str(profile.get("description") or ""))
        self.custom_ranking_weights = None
        self.ranking_profile_status.setText("Preset weights")

    @QtCore.pyqtSlot()
    def customize_ranking_priorities(self) -> None:
        profile_id = str(self.ranking_profile_combo.currentData() or "default")
        base = resolve_ranking_preferences(
            ChemistRankingPreferences(profile_id=profile_id)
        )
        dialog = RankingPreferencesDialog(
            self.custom_ranking_weights or base.weights,
            parent=self,
        )
        if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        self.custom_ranking_weights = dialog.weights
        self.ranking_profile_status.setText("Custom normalized weights")

    def _ranking_preferences(self) -> ChemistRankingPreferences:
        profile_id = str(self.ranking_profile_combo.currentData() or "default")
        return ChemistRankingPreferences(
            profile_id=(
                f"{profile_id}:custom" if self.custom_ranking_weights else profile_id
            ),
            weights=dict(self.custom_ranking_weights or {}),
            customized=bool(self.custom_ranking_weights),
        )

    @QtCore.pyqtSlot()
    def _update_data_summary(self) -> None:
        path = Path(self.data_path_edit.text().strip())
        if not path.is_file():
            self.data_summary.setText("Data file not found.")
            return
        report_path = path.parent / "recommendation_artifacts_report.json"
        if report_path.is_file():
            try:
                report = json.loads(report_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                report = {}
            total = report.get("record_count")
            trusted = report.get(
                "trusted_precedent_count",
                report.get("eligible_index_record_count"),
            )
            review_core = report.get("review_core_precedent_count")
            query_core = report.get("query_core_eligible_count")
            if total is not None:
                details = []
                if trusted is not None:
                    details.append(f"{trusted} trusted precedents")
                if review_core is not None:
                    details.append(f"{review_core} review-core precedents")
                if query_core is not None:
                    details.append(f"{query_core} query-core eligible")
                self.data_summary.setText(
                    f"{total} converted reactions"
                    + (f"; {'; '.join(details)}" if details else "")
                    + "."
                )
                return
        self.data_summary.setText(
            f"Using {path.name} ({path.stat().st_size / 1024:.1f} KB)."
        )

    @QtCore.pyqtSlot()
    def load_example(self) -> None:
        self.reaction_edit.setText("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")

    @QtCore.pyqtSlot()
    def clear_results(self) -> None:
        self.summary_box.clear()
        self.reaction_image_label.clear_image()
        self.selected_reaction_image_label.clear_image()
        self.details_box.clear()
        self.results_table.setSortingEnabled(False)
        self.results_table.setRowCount(0)
        self.results_table.setSortingEnabled(True)
        self.last_result = None
        self.export_button.setEnabled(False)
        self.status_label.setText("Ready")

    @QtCore.pyqtSlot()
    def start_recommendation(self) -> None:
        if self.thread is not None:
            return
        data_path = Path(self.data_path_edit.text().strip())
        reaction_smiles = self.reaction_edit.text().strip()
        if not data_path.is_file():
            QtWidgets.QMessageBox.warning(
                self,
                "Recommendation data required",
                "Choose an existing generic index or shard manifest.",
            )
            return
        if not reaction_smiles:
            QtWidgets.QMessageBox.warning(
                self,
                "Reaction required",
                "Enter a reaction SMILES including its product.",
            )
            return
        if ">>" not in reaction_smiles:
            QtWidgets.QMessageBox.warning(
                self,
                "Product required",
                "Use reaction SMILES in reactants>>product form.",
            )
            return
        if (
            self.use_rxnmapper_check.isChecked()
            and not RxnMapperProvider.is_available()
        ):
            QtWidgets.QMessageBox.warning(
                self,
                "RXNMapper unavailable",
                "RXNMapper is checked but not installed. Run "
                "'python -m pip install -r requirements-mapping.txt' or "
                "clear the RXNMapper checkbox.",
            )
            return

        completion_selections: Tuple[ReactionCompletionSelection, ...] = ()
        if self.mode_combo.currentData() == "recommendation":
            try:
                completion_proposal = propose_reaction_completion(reaction_smiles)
            except ValueError as exc:
                QtWidgets.QMessageBox.warning(
                    self,
                    "Reaction could not be analyzed",
                    str(exc),
                )
                return
            if completion_proposal.requirements:
                completion_dialog = CompletionSourceDialog(
                    completion_proposal,
                    parent=self,
                )
                if (
                    completion_dialog.exec()
                    != QtWidgets.QDialog.DialogCode.Accepted
                ):
                    return
                completion_selections = completion_dialog.selections

        self.clear_results()
        self._render_reaction_graph(reaction_smiles)
        self.run_button.setEnabled(False)
        self.status_label.setText("Starting…")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)

        minimum = self.minimum_pool_spin.value()
        thread = QtCore.QThread(self)
        if self.mode_combo.currentData() == "discovery":
            worker = ReactionDiscoveryWorker(
                str(data_path),
                reaction_smiles,
                top_k=self.top_k_spin.value(),
                view=str(self.discovery_view_combo.currentData()),
                include_low_yield=self.include_low_yield_check.isChecked(),
                include_unreported_outcomes=(self.include_unreported_check.isChecked()),
                use_rxnmapper=self.use_rxnmapper_check.isChecked(),
                include_review=self.unrestricted_fallback_check.isChecked(),
            )
        else:
            worker = GenericRecommendationWorker(
                str(data_path),
                reaction_smiles,
                top_k=self.top_k_spin.value(),
                minimum_pool_size=minimum or None,
                unrestricted_fallback=(self.unrestricted_fallback_check.isChecked()),
                use_rxnmapper=self.use_rxnmapper_check.isChecked(),
                ranking_preferences=self._ranking_preferences(),
                completion_selections=completion_selections,
            )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self.status_label.setText)
        worker.finished.connect(self._on_finished)
        worker.finished.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_thread)
        self.thread = thread
        self.worker = worker
        thread.start()

    @QtCore.pyqtSlot(bool, object, str)
    def _on_finished(
        self,
        success: bool,
        result: Optional[object],
        error: str,
    ) -> None:
        self.run_button.setEnabled(True)
        self.progress_bar.setVisible(False)
        if not success or result is None:
            self.status_label.setText("Error")
            self.summary_box.setPlainText(error)
            return
        self.last_result = result
        self.export_button.setEnabled(True)
        if isinstance(result, ReactionDiscoveryResult):
            self._render_discovery_result(result)
        elif isinstance(result, GenericRecommendationResult):
            self._render_result(result)
        else:
            self.status_label.setText("Error")
            self.summary_box.setPlainText("Unsupported result contract")

    def _render_discovery_result(
        self,
        result: ReactionDiscoveryResult,
    ) -> None:
        self._render_reaction_graph(result.query_reaction_smiles)
        self.summary_box.setPlainText(
            format_discovery_summary(result)
            + (
                f"\nDiscovery: {_friendly_error(result.error)}"
                if not result.valid
                else ""
            )
        )
        hits = tuple(result.hits) if result.valid else ()
        self.results_table.setSortingEnabled(False)
        self.results_table.setRowCount(len(hits))
        for row, hit in enumerate(hits):
            components = hit.score_trace.components
            values = (
                str(hit.rank),
                _display_name(hit.relation_class),
                f"{hit.discovery_score:.3f}",
                _format_discovery_factor(components.get("edit_similarity")),
                _format_discovery_factor(components.get("reaction_center")),
                _format_discovery_factor(components.get("local_environment")),
                _format_discovery_factor(components.get("partner_category")),
                (f"{hit.yield_pct:.1f}%" if hit.yield_pct is not None else "—"),
                format_recipe_summary(hit.resolved_recipe),
                "; ".join(hit.cautions),
            )
            for column, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                if column == 0:
                    item.setData(QtCore.Qt.ItemDataRole.UserRole, hit)
                self.results_table.setItem(row, column, item)
        self.results_table.setSortingEnabled(True)
        self.results_table.resizeColumnsToContents()
        self.results_table.horizontalHeader().setSectionResizeMode(
            8,
            QtWidgets.QHeaderView.ResizeMode.Stretch,
        )
        self.status_label.setText(
            f"Done — {len(hits)} analogue(s)" if hits else "No structural analogue"
        )
        if hits:
            self.results_table.selectRow(0)

    def _render_result(self, result: GenericRecommendationResult) -> None:
        self._render_reaction_graph(result.query_reaction_smiles)
        if not result.valid:
            self.status_label.setText("No recommendation")
            summary = format_query_summary(result)
            self.summary_box.setPlainText(
                f"{summary}\nRecommendation: {_friendly_error(result.error)}"
            )
            self.results_table.setRowCount(0)
            return
        self.summary_box.setPlainText(format_query_summary(result))
        recommendations = tuple(result.recommendations)
        self.results_table.setSortingEnabled(False)
        self.results_table.setRowCount(len(recommendations))
        for row, recommendation in enumerate(recommendations):
            values = (
                str(recommendation.rank),
                (
                    str(recommendation.default_rank)
                    if recommendation.default_rank is not None
                    else "—"
                ),
                f"{recommendation.score:.3f}",
                f"{recommendation.similarity_score:.3f}",
                f"{recommendation.compatibility_score:.3f}",
                (
                    f"{recommendation.expected_yield_pct:.1f}%"
                    if recommendation.expected_yield_pct is not None
                    else "—"
                ),
                str(recommendation.support),
                str(recommendation.reference_support),
                format_recipe_summary(recommendation.resolved_recipe),
                "; ".join(recommendation.cautions),
            )
            for column, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                if column == 0:
                    item.setData(
                        QtCore.Qt.ItemDataRole.UserRole,
                        recommendation,
                    )
                self.results_table.setItem(row, column, item)
        self.results_table.setSortingEnabled(True)
        self.results_table.resizeColumnsToContents()
        self.results_table.horizontalHeader().setSectionResizeMode(
            8,
            QtWidgets.QHeaderView.ResizeMode.Stretch,
        )
        self.status_label.setText(f"Done — {len(recommendations)} recipe(s)")
        if recommendations:
            self.results_table.selectRow(0)

    def _render_reaction_graph(self, reaction_smiles: str) -> None:
        """Render the query graph without making visualization outcome-critical."""
        try:
            drawing = render_reaction_image_bytes(
                reaction_smiles,
                size=QUERY_REACTION_IMAGE_SIZE,
                image_format="svg",
            )
        except (RuntimeError, ValueError) as exc:
            self.reaction_image_label.clear_image("Reaction graph unavailable.")
            self.reaction_image_label.setToolTip(str(exc))
            return
        if not self.reaction_image_label.set_image_bytes(drawing):
            self.reaction_image_label.setToolTip(
                "The renderer returned an unsupported image."
            )
            return
        self.reaction_image_label.setToolTip(reaction_smiles)

    def _render_selected_reaction_graph(self, reaction_smiles: str) -> None:
        """Render the first precedent associated with the selected recipe."""
        if not reaction_smiles:
            self.selected_reaction_image_label.clear_image(
                "No precedent reaction structure is available."
            )
            return
        try:
            drawing = render_reaction_image_bytes(
                reaction_smiles,
                size=PRECEDENT_REACTION_IMAGE_SIZE,
                image_format="svg",
            )
        except (RuntimeError, ValueError) as exc:
            self.selected_reaction_image_label.clear_image(
                "Precedent reaction graph unavailable."
            )
            self.selected_reaction_image_label.setToolTip(str(exc))
            return
        if not self.selected_reaction_image_label.set_image_bytes(drawing):
            self.selected_reaction_image_label.setToolTip(
                "The renderer returned an unsupported image."
            )
            return
        self.selected_reaction_image_label.setToolTip(reaction_smiles)

    @QtCore.pyqtSlot()
    def _show_selected_details(self) -> None:
        row = self.results_table.currentRow()
        if row < 0:
            self.details_box.clear()
            self.selected_reaction_image_label.clear_image()
            return
        item = self.results_table.item(row, 0)
        if item is None:
            return
        recommendation = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if recommendation is None:
            return
        if isinstance(recommendation, ReactionDiscoveryHit):
            self._show_discovery_details(recommendation)
            return
        precedent_reaction_smiles = tuple(recommendation.precedent_reaction_smiles)
        self._render_selected_reaction_graph(
            precedent_reaction_smiles[0] if precedent_reaction_smiles else ""
        )
        recipe = recommendation.resolved_recipe
        precedent_contexts = tuple(recommendation.precedent_reaction_contexts)
        selected_context = (
            precedent_contexts[0]
            if precedent_contexts and isinstance(precedent_contexts[0], Mapping)
            else {}
        )
        spectator_groups = tuple(
            dict(group)
            for group in (selected_context.get("spectator_groups") or ())
            if isinstance(group, Mapping)
        )
        reaction_partners = tuple(
            dict(partner)
            for partner in (selected_context.get("reaction_partners") or ())
            if isinstance(partner, Mapping)
        )
        fragment_source_support = tuple(
            dict(value)
            for value in (selected_context.get("fragment_source_support") or ())
            if isinstance(value, Mapping)
        )
        source_support_summary = (
            "; ".join(
                (
                    f"{_display_name(value.get('status'))}: "
                    + ", ".join(
                        str(identifier)
                        for identifier in (value.get("component_raw_identifiers") or ())
                    )
                ).rstrip(": ")
                for value in fragment_source_support
            )
            or "Not required"
        )
        partner_summaries = _partner_analysis_summaries(reaction_partners)
        lines = [
            f"Rank {recommendation.rank}",
            (
                f"Default rank {recommendation.default_rank}; "
                f"rank change {recommendation.rank_change:+d}"
                if recommendation.default_rank is not None
                else "Default rank unavailable"
            ),
            "",
            "Displayed precedent context",
            (
                "Reaction label: "
                f"{(selected_context.get('reaction_label') or {}).get('text') or 'Unresolved'}"
            ),
            (
                f"Selected hit reaction SMILES: {precedent_reaction_smiles[0]}"
                if precedent_reaction_smiles
                else "Selected hit reaction SMILES: Unavailable"
            ),
            f"Fragment source support: {source_support_summary}",
            f"Spectator groups: {_spectator_summary(spectator_groups)}",
            "Reactivity profile:",
        ]
        if partner_summaries:
            lines.extend(f"• {summary}" for summary in partner_summaries)
        else:
            lines.append("None available")
        lines.extend(("", "Conditions"))
        for field, label in _RECIPE_ROLE_LABELS:
            names = _component_names(recipe, field)
            if names:
                lines.append(f"{label}: {names}")
        for field, label, unit in (
            ("temperature_c", "Temperature", "°C"),
            ("time_h", "Time", "h"),
            ("concentration_m", "Concentration", "M"),
            ("pressure_bar", "Pressure", "bar"),
        ):
            value = recipe.get(field)
            if value is not None:
                lines.append(f"{label}: {value:g} {unit}")
        atmosphere = str(recipe.get("atmosphere") or "").strip()
        if atmosphere:
            lines.append(f"Atmosphere: {atmosphere}")
        lines.extend(
            (
                "",
                "Evidence",
                (
                    f"Score: {recommendation.score:.3f}; similarity: "
                    f"{recommendation.similarity_score:.3f}; compatibility: "
                    f"{recommendation.compatibility_score:.3f}"
                ),
                (
                    f"Default score: {recommendation.default_score:.3f}"
                    if recommendation.default_score is not None
                    else "Default score: unavailable"
                ),
                (
                    f"Support: {recommendation.support} reaction(s), "
                    f"{recommendation.reference_support} reference(s), "
                    f"{recommendation.dataset_support} dataset(s)"
                ),
            )
        )
        if recommendation.expected_yield_pct is not None:
            lines.append(
                f"Evidence-weighted expected yield: "
                f"{recommendation.expected_yield_pct:.1f}%"
            )
        trace = recommendation.score_trace
        lines.extend(("", f"Ranking profile: {trace.ranking_profile}", "Score factors"))
        for component, weight in trace.applied_ranking_weights.items():
            value = trace.ranking_components.get(component)
            contribution = trace.ranking_contributions.get(component, 0.0)
            value_text = "unavailable" if value is None else f"{value:.3f}"
            lines.append(
                f"• {_RANKING_COMPONENT_LABELS.get(component, _display_name(component))}: "
                f"value {value_text}; weight {weight:.3f}; "
                f"contribution {contribution:.3f}"
            )
        partner_evidence = recommendation.factor_evidence.get("partner_category", {})
        query_categories = tuple(partner_evidence.get("query_categories") or ())
        if query_categories:
            lines.extend(("", "Reactant-category evidence"))
            lines.append("• Query: " + ", ".join(query_categories))
            for assessment in partner_evidence.get("precedent_assessments") or ():
                precedent_categories = tuple(
                    assessment.get("precedent_categories") or ()
                )
                lines.append(
                    "• Precedent "
                    f"{assessment.get('evidence_id')}: "
                    + (", ".join(precedent_categories) or "unavailable")
                    + (
                        f" (score {assessment.get('score'):.3f})"
                        if assessment.get("score") is not None
                        else ""
                    )
                )
        tolerance = recommendation.factor_evidence.get("functional_group_tolerance", {})
        tolerance_groups = tuple(tolerance.get("groups") or ())
        if tolerance_groups:
            lines.extend(("", "Functional-group tolerance evidence"))
            for group in tolerance_groups:
                lines.append(
                    f"• {group.get('chemist_label') or group.get('group_id')}: "
                    f"{_display_name(group.get('status'))}; "
                    f"{group.get('matched_independent_count', 0)} independent "
                    "precedent(s)"
                )
        if recommendation.explanation:
            lines.extend(("", "Why this recipe"))
            lines.extend(f"• {value}" for value in recommendation.explanation)
        if recommendation.compatibility_evidence:
            lines.extend(("", "Compatibility"))
            lines.extend(
                f"• {value}" for value in recommendation.compatibility_evidence
            )
        if recommendation.cautions:
            lines.extend(("", "Cautions"))
            lines.extend(f"• {value}" for value in recommendation.cautions)
        if recommendation.precedent_reference_ids:
            lines.extend(
                (
                    "",
                    "Precedent references",
                    *recommendation.precedent_reference_ids,
                )
            )
        if recommendation.precedent_reaction_ids:
            lines.extend(
                (
                    "",
                    "Precedent reactions",
                    *recommendation.precedent_reaction_ids,
                )
            )
        self.details_box.setPlainText("\n".join(lines))

    def _show_discovery_details(self, hit: ReactionDiscoveryHit) -> None:
        """Explain one analogue without presenting its conditions as advice."""
        self._render_selected_reaction_graph(hit.reaction_smiles)
        recipe = hit.resolved_recipe
        lines = [
            f"Rank {hit.rank}: {_display_name(hit.relation_class)}",
            f"Discovery score: {hit.discovery_score:.3f}",
            f"Precedent reaction: {hit.reaction_smiles}",
            f"Evidence tier: {_display_name(hit.evidence_tier)}",
            f"Chemistry status: {_display_name(hit.chemistry_status)}",
            f"Source: {hit.source_dataset or 'Unavailable'}",
            f"Reference: {hit.reference_id or 'Unavailable'}",
            (
                f"Observed yield: {hit.yield_pct:.1f}%"
                if hit.yield_pct is not None
                else "Observed yield: Unreported"
            ),
        ]
        if hit.hypothesis_id:
            lines.append(f"Query edit hypothesis: {hit.hypothesis_id}")
        lines.extend(("", "Why it is related"))
        lines.extend(f"• {value}" for value in hit.score_trace.matches)
        if hit.score_trace.mismatches:
            lines.extend(("", "Structural differences"))
            lines.extend(f"• {value}" for value in hit.score_trace.mismatches)
        lines.extend(("", "Discovery score factors"))
        for name, configured_weight in hit.score_trace.configured_weights.items():
            value = hit.score_trace.components.get(name)
            effective = hit.score_trace.effective_weights.get(name)
            contribution = hit.score_trace.contributions.get(name)
            lines.append(
                f"• {_display_name(name)}: "
                f"value {_format_discovery_factor(value)}; "
                f"configured weight {configured_weight:.3f}; "
                + (
                    f"effective weight {effective:.3f}; contribution {contribution:.3f}"
                    if effective is not None and contribution is not None
                    else "unavailable for this comparison"
                )
            )
        lines.extend(
            (
                "",
                "Observed conditions (precedent evidence, not a recommendation)",
            )
        )
        for field, label in _RECIPE_ROLE_LABELS:
            names = _component_names(recipe, field)
            if names:
                lines.append(f"{label}: {names}")
        for field, label, unit in (
            ("temperature_c", "Temperature", "°C"),
            ("time_h", "Time", "h"),
            ("concentration_m", "Concentration", "M"),
            ("pressure_bar", "Pressure", "bar"),
        ):
            value = recipe.get(field)
            if value is not None:
                lines.append(f"{label}: {value:g} {unit}")
        if hit.insights:
            lines.extend(("", "Insights"))
            lines.extend(f"• {value}" for value in hit.insights)
        if hit.cautions:
            lines.extend(("", "Cautions"))
            lines.extend(f"• {value}" for value in hit.cautions)
        self.details_box.setPlainText("\n".join(lines))

    @QtCore.pyqtSlot()
    def export_json(self) -> None:
        if self.last_result is None:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export analysis result",
            str(
                PROJECT_ROOT
                / "results"
                / (
                    "reaction_discovery.json"
                    if isinstance(self.last_result, ReactionDiscoveryResult)
                    else "generic_recommendation.json"
                )
            ),
            "JSON files (*.json)",
        )
        if not path:
            return
        destination = Path(path)
        if destination.suffix.casefold() != ".json":
            destination = destination.with_suffix(".json")
        try:
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(
                json.dumps(
                    self.last_result.to_dict(),
                    ensure_ascii=False,
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
        except OSError as exc:
            QtWidgets.QMessageBox.critical(
                self,
                "Export failed",
                str(exc),
            )
            return
        self.status_label.setText(f"Exported {destination}")

    @QtCore.pyqtSlot()
    def _clear_thread(self) -> None:
        self.thread = None
        self.worker = None

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        if self.thread is not None and self.thread.isRunning():
            QtWidgets.QMessageBox.information(
                self,
                "Recommendation running",
                "Wait for the current recommendation to finish before closing.",
            )
            event.ignore()
            return
        event.accept()


def _show_main_window(window: GenericRecommenderWindow) -> None:
    """Show the application window in its requested initial state."""
    window.showMaximized()


def main() -> None:
    """Launch the generic condition recommender."""
    application = QtWidgets.QApplication(sys.argv)
    window = GenericRecommenderWindow()
    _show_main_window(window)
    raise SystemExit(application.exec())


if __name__ == "__main__":
    main()


__all__ = [
    "GenericRecommendationWorker",
    "GenericRecommenderWindow",
    "default_recommendation_data_path",
    "format_recipe_summary",
    "main",
]
