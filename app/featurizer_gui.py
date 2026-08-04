"""Small Qt6 wrapper for molecule and reaction featurization."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Literal, Sequence

from PyQt6 import QtCore, QtGui, QtWidgets

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy import (  # noqa: E402
    AtomMappingProvider,
    ExternalMappingAssessment,
    RxnMapperProvider,
    analyze_reaction_with_external_mapping,
    build_reaction_display_projection,
    analyze_molecule,
    featurize_reaction,
    reaction_render_context_from_analysis,
)
from reactive_taxonomy.cli import format_concise_analysis  # noqa: E402
from visualization import (  # noqa: E402
    available_render_presets,
    build_reaction_display_graphic,
    render_molecule_image_bytes,
    render_reaction_image_bytes,
)
from visualization.qt_widgets import StructureImageLabel  # noqa: E402


InputKind = Literal["molecule", "reaction"]

REACTION_EXAMPLE = (
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
)
MOLECULE_EXAMPLE = "Brc1ccc(N)cc1C#N"
REACTION_IMAGE_SIZE = (680, 168)
MOLECULE_IMAGE_SIZE = (480, 300)
DEFAULT_RENDER_PRESET = "current"


def detect_input_kind(text: str) -> InputKind:
    """Classify input from reaction-SMILES delimiters.

    The ``>`` character is reserved for reaction separators and is not part of
    ordinary molecular SMILES. Detecting any separator also routes malformed
    reaction input to the reaction parser, which can return the useful error.
    """
    return "reaction" if ">" in text else "molecule"


def featurize_text(
    text: str,
    *,
    mapping_provider: AtomMappingProvider | None = None,
    force_resolved_mapping: bool = False,
) -> tuple[InputKind, object, ExternalMappingAssessment | None]:
    """Featurize stripped text through the appropriate public taxonomy API."""
    value = text.strip()
    if not value:
        raise ValueError("Enter a molecule or reaction SMILES.")
    kind = detect_input_kind(value)
    if kind == "reaction":
        base_analysis = featurize_reaction(value)
        assessment = (
            analyze_reaction_with_external_mapping(
                value,
                mapping_provider,
                base_analysis=base_analysis,
                force_resolved_shadow=force_resolved_mapping,
            )
            if mapping_provider is not None
            else None
        )
        return (
            kind,
            assessment.analysis if assessment is not None else base_analysis,
            assessment,
        )
    return kind, analyze_molecule(value), None


def _readable(value: object) -> str:
    return str(value or "").replace("_", " ").strip().lower()


def _evidence_text(value: object) -> str:
    labels = {
        "validated_atom_mapping": "validated atom mapping",
        "validated_mapping": "validated atom mapping",
        "global_atom_correspondence": "whole-reaction atom correspondence",
        "fragmented_scaffold_correspondence": "fragment-aware atom correspondence",
        "scaffold_correspondence": "scaffold atom correspondence",
        "ambiguous_correspondence_core_consensus": (
            "consensus across atom-correspondence alternatives"
        ),
    }
    raw = str(value or "")
    return labels.get(raw, _readable(raw))


def _annotation_text(value: object) -> str:
    labels = {
        "sp2_c_n_substitution_like": "sp² C–N substitution",
        "sp2_c_o_substitution_like": "sp² C–O substitution",
        "sp2_c_s_substitution_like": "sp² C–S substitution",
        "organoboron_c_c_coupling_like": "organoboron C–C coupling",
        "alkene_hydrogenation_like": "alkene hydrogenation",
        "generic_graph_transformation": "graph-defined transformation",
        "generic_multi_event_graph_transformation": (
            "multi-event graph transformation"
        ),
        "suzuki_miyaura": "Suzuki–Miyaura",
        "buchwald_hartwig_c_n": "Buchwald–Hartwig C–N coupling",
        "snar_amination": "SNAr amination",
        "ullmann_c_n": "Ullmann C–N coupling",
    }
    raw = str(value or "")
    return labels.get(raw, _readable(raw))


def _quality_note(value: object) -> str:
    labels = {
        "not_all_edits_graph_checked": (
            "some bond changes could not be independently checked"
        ),
        "partial_active_atom_mapping": (
            "atom provenance is incomplete at the reaction center"
        ),
        "remote_continuity_unresolved": "R-group continuity is not fully resolved",
        "active_core_size_exceeds_review_limit": "reaction core is unusually large",
        "event_count_exceeds_review_limit": "reaction contains many connected events",
        "primary_center_has_no_state_change": (
            "the proposed primary center has no recorded state change"
        ),
    }
    raw = str(value or "")
    return labels.get(raw, _readable(raw))


def _warning_text(value: object) -> str:
    labels = {
        "INFERRED_GLOBAL_ATOM_CORRESPONDENCE": "atom correspondence was inferred",
        "CORE_SHARED_BY_ALL_EDIT_HYPOTHESES": (
            "all mapping alternatives support the same reaction core"
        ),
        "REACTION_COMPLETENESS_UNRESOLVED": (
            "individual product-atom origins are not fully resolved"
        ),
    }
    raw = str(value or "")
    return labels.get(raw, _readable(raw))


def _mapping_status_text(value: object) -> str:
    labels = {
        "not_requested_invalid_reaction": "not run; reaction is invalid",
        "not_requested_supplied_mapping": "not needed; supplied mapping used",
        "not_requested_resolved_internal_evidence": (
            "not needed; internal correspondence resolved"
        ),
        "external_mapping_failed": "failed",
        "external_mapping_internal_consensus": (
            "supports the internal correspondence"
        ),
        "external_mapping_hypothesis_conflict": (
            "conflicts with an internal mapping alternative"
        ),
        "external_mapping_signature_conflict": (
            "conflicts with the internal reaction signature"
        ),
        "external_mapping_ambiguous_hypothesis_match": (
            "matches more than one internal mapping alternative"
        ),
        "external_mapping_only": "external mapping only; expert review required",
        "external_mapping_signature_unavailable": (
            "mapping produced, but no verified signature is available"
        ),
    }
    raw = str(value or "")
    return labels.get(raw, _readable(raw))


def _bond_text(value: str) -> str:
    return value.replace("-", "–")


def _edit_text(token: str) -> str:
    """Translate one normalized edit token into concise chemical language."""
    try:
        edit_type, atoms, transition = str(token).split(":", 2)
        before, after = transition.split(">", 1)
    except ValueError:
        return _readable(token)
    bond = _bond_text(atoms)
    if edit_type == "formed":
        return f"form {bond}"
    if edit_type == "broken":
        return f"break {bond}"
    if edit_type == "hydrogen_change":
        return f"remove {bond}" if after == "NONE" else f"add {bond}"
    if edit_type == "order_changed":
        return (
            f"change {bond} bond from {_readable(before)} "
            f"to {_readable(after)}"
        )
    return f"{_readable(edit_type)} at {bond}"


def _attachment_text(profile: Any) -> str:
    hybridization = {
        "SP3": "sp³",
        "SP2": "sp²",
        "SP": "sp",
    }.get(str(profile.attachment_hybridization).upper(), "")
    atom = str(profile.attachment_element)
    center = f"{hybridization} {atom}" if hybridization else atom
    order = _readable(profile.attachment_bond_order)
    return f"attached through {center} by a {order} bond"


def _site_profile_text(profile: Any) -> str:
    """Render the most decision-relevant site context in chemical language."""
    context = profile.context
    center = profile.reactive_center
    kind = str(profile.context_kind)
    values = []
    if kind == "heteroatom":
        substitution = _readable(context.substitution_class)
        values.append(f"{substitution} {context.element}")
        resonance = _readable(context.resonance_class)
        if resonance:
            values.append(resonance)
        if center.lone_pair_availability:
            values.append(
                f"{_readable(center.lone_pair_availability)} "
                "lone-pair availability"
            )
    elif kind == "aromatic":
        values.append(f"{_readable(context.ring_family)} aromatic {center.element}")
        values.append(f"{_readable(context.ortho_burden_class)} ortho crowding")
    elif kind == "alkyl":
        values.append(
            f"{_readable(context.carbon_substitution)} sp³ carbon"
        )
        for name in ("benzylic", "allylic", "propargylic", "cyclic"):
            if getattr(context, name, False):
                values.append(name)
    elif kind in {"alkenyl", "alkynyl"}:
        values.append(_readable(kind))
        values.append(_readable(context.conjugation_class))
    elif kind in {"acyl", "sulfonyl", "phosphoryl"}:
        values.append(_readable(context.center_class))
    else:
        values.append(f"{_readable(kind)} {center.element}")

    access = _readable(profile.steric.accessibility_class)
    values.append(f"{access} access")
    electronic = profile.electronic
    if electronic.activation_axis != "lone_pair_availability":
        electronic_class = _readable(electronic.activation_class).replace(
            "electron poor", "electron-poor"
        ).replace("electron rich", "electron-rich")
        values.append(electronic_class)
    if any(
        modifier.modifier_type == "coordination"
        for modifier in profile.modifiers
    ):
        values.append("coordination risk")
    return "; ".join(value for value in values if value)


def _port_profile_text(profile: Any) -> str:
    """Render one attachment-port profile without interpreting its R label."""
    class_name = {
        "ring_aliphatic": "saturated ring",
        "generic_R": "organic substituent",
    }.get(str(profile.base_class), _readable(profile.base_class))
    values = []
    if profile.carbon_substitution not in {
        "not_carbon",
        "not_applicable",
        "unresolved",
    }:
        values.append(f"{profile.carbon_substitution} carbon attachment")
    values.append(class_name)
    values.extend(
        name
        for name, present in (
            ("benzylic", profile.benzylic),
            ("allylic", profile.allylic),
            ("propargylic", profile.propargylic),
        )
        if present
    )
    if profile.ring_sizes:
        values.append(
            "/".join(str(size) for size in profile.ring_sizes)
            + "-membered ring"
        )
    if profile.alpha_branch_count == 0 and profile.beta_branch_count == 0:
        values.append("unbranched near the attachment")
    else:
        if profile.alpha_branch_count:
            values.append(f"α branching={profile.alpha_branch_count}")
        if profile.beta_branch_count:
            values.append(f"β branching={profile.beta_branch_count}")
    if profile.radius_1_heteroatoms:
        values.append(
            "/".join(profile.radius_1_heteroatoms) + " one bond away"
        )
    if profile.radius_2_heteroatoms:
        values.append(
            "/".join(profile.radius_2_heteroatoms) + " two bonds away"
        )
    if profile.aromatic_substituent_relations:
        relations = sorted(
            {
                f"{_readable(item.positional_relation)} "
                f"{item.substituent_fragment_smiles or item.substituent_element}"
                for item in profile.aromatic_substituent_relations
            }
        )
        values.append("ring substituents: " + ", ".join(relations))
    values.append(_attachment_text(profile))
    return "; ".join(values)


def _r_group_functional_context_lines(
    analysis: Any,
    ports: Sequence[Any],
) -> list[str]:
    """Render motif overlays already associated with remote core subgraphs."""
    interpretation = getattr(analysis, "interpretation", None)
    contexts = tuple(
        getattr(interpretation, "r_group_functional_contexts", ()) or ()
    )
    lines = []
    for context in contexts:
        matching_ports = tuple(
            value
            for value in ports
            if value.source_component_index == context.component_index
            and tuple(value.atom_indices) == tuple(context.remote_atom_indices)
        )
        if not matching_ports:
            continue
        matching_ports = tuple(
            sorted(
                matching_ports,
                key=lambda value: int(value.placeholder_index or 0),
            )
        )
        labels = tuple(str(value.display_label) for value in matching_ports)
        label_text = "/".join(labels)
        label_by_attachment = {
            int(value.attachment_atom_index): str(value.display_label)
            for value in matching_ports
        }
        lines.append(
            f"Unchanged functional groups on the {label_text} scaffold:"
        )
        for group in context.functional_groups:
            port_values = tuple(
                (
                    label_by_attachment[distance.attachment_atom_index],
                    distance.bond_distance,
                )
                for distance in group.port_distances
                if distance.attachment_atom_index in label_by_attachment
            )
            distance_notes = []
            if port_values:
                unique_distances = {distance for _, distance in port_values}
                if len(port_values) > 1 and len(unique_distances) == 1:
                    distance_notes.append(
                        f"{next(iter(unique_distances))} bonds from each R attachment"
                    )
                elif len(port_values) == 1:
                    label, distance = port_values[0]
                    distance_notes.append(
                        f"{distance} bonds from the {label} attachment"
                    )
                else:
                    distance_notes.append(
                        ", ".join(
                            f"{label}: {distance} bonds"
                            for label, distance in port_values
                        )
                    )
            if group.distance_to_reactive_site is not None:
                distance_notes.append(
                    f"{group.distance_to_reactive_site} bonds from the reaction center"
                )
            suffix = (
                f" ({'; '.join(distance_notes)})" if distance_notes else ""
            )
            lines.append(f"  {_annotation_text(group.motif_id)}{suffix}")
    return lines


def _r_group_analysis_lines(analysis: Any) -> list[str]:
    """Return display-labeled R-port profiles from the canonical projection."""
    try:
        projection = build_reaction_display_projection(
            reaction_render_context_from_analysis(analysis)
        )
    except (RuntimeError, ValueError) as exc:
        return [f"R-group attachment profiles: unavailable ({exc})"]

    candidates = tuple(
        value
        for value in projection.substituents
        if value.boundary_action == "r_placeholder"
        and value.placeholder_index is not None
        and value.display_label
    )
    by_index = {}
    for value in candidates:
        index = int(value.placeholder_index)
        current = by_index.get(index)
        if current is None or (
            value.side == "reactant" and current.side != "reactant"
        ):
            by_index[index] = value
    ports = tuple(by_index[index] for index in sorted(by_index))
    if not ports:
        return ["R-group attachment profiles: none"]

    lines = ["R-group attachment profiles:"]
    continuity_labels = {
            "retained": "retained",
            "departing": "departing group",
            "appearing": "product-only group",
            "changed": "changed during reaction",
            "unresolved": "continuity unresolved",
    }
    profile_groups = {}
    for value in ports:
        key = (
            value.substituent_profile.profile_id,
            str(value.continuity),
        )
        profile_groups.setdefault(key, []).append(value)
    for (_, continuity), values in sorted(
        profile_groups.items(),
        key=lambda item: min(
            int(value.placeholder_index or 0) for value in item[1]
        ),
    ):
        ordered = sorted(
            values,
            key=lambda item: int(item.placeholder_index or 0),
        )
        labels = "/".join(str(value.display_label) for value in ordered)
        profile = ordered[0].substituent_profile
        continuity_text = continuity_labels.get(
            continuity,
            _readable(continuity),
        )
        lines.append(
            f"  {labels} — {continuity_text}; "
            f"{_port_profile_text(profile)}"
        )

    label_indices = {
        str(value.display_label): int(value.placeholder_index)
        for value in ports
    }
    shared: dict[tuple[int, tuple[int, ...], str], list[str]] = {}
    for value in ports:
        key = (
            int(value.source_component_index),
            tuple(value.atom_indices),
            str(value.fragment_smiles),
        )
        shared.setdefault(key, []).append(str(value.display_label))
    for _, labels in sorted(
        shared.items(),
        key=lambda item: min(label_indices[label] for label in item[1]),
    ):
        ordered_labels = sorted(
            set(labels),
            key=label_indices.__getitem__,
        )
        if len(ordered_labels) > 1:
            label_text = " and ".join(ordered_labels)
            lines.append(
                f"  {label_text} are connected through the same omitted scaffold."
            )
    lines.extend(_r_group_functional_context_lines(analysis, ports))
    connectors = {}
    for value in projection.connectors:
        key = tuple(value.placeholder_indices)
        current = connectors.get(key)
        if current is None or (
            value.side == "reactant" and current.side != "reactant"
        ):
            connectors[key] = value
    if connectors:
        lines.append("Hidden core connectors:")
        for value in sorted(
            connectors.values(),
            key=lambda item: item.placeholder_indices,
        ):
            path = (
                f"{value.shortest_path_bond_count} omitted bonds"
                if value.shortest_path_bond_count is not None
                else f"{value.hidden_atom_count} omitted atoms"
            )
            lines.append(
                f"  {' and '.join(value.port_display_labels)} remain connected "
                f"through {path}."
            )
    return lines


def _active_site_profile_lines(analysis: Any) -> list[str]:
    """Render annotations for reactant sites participating in core edits."""
    core = getattr(analysis, "reaction_core", None)
    if core is None:
        return ["Active-site steric/electronic context: unavailable"]
    roles = {}
    for transition in core.atom_transitions:
        state = transition.before_state
        if state is not None:
            roles[(state.component_index, state.atom_index)] = transition.role
    values = []
    for component in getattr(analysis, "reactants", ()):
        molecule_analysis = getattr(component, "molecule_analysis", None)
        if molecule_analysis is None:
            continue
        for environment in molecule_analysis.reactive_site_environments:
            key = (component.component_index, environment.center_atom_index)
            if key not in roles or environment.reactivity_profile is None:
                continue
            values.append(
                (
                    component.component_index,
                    roles[key],
                    _site_profile_text(environment.reactivity_profile),
                )
            )
    if not values:
        return ["Active-site steric/electronic context: unavailable"]
    lines = ["Active-site steric/electronic context:"]
    for component_index, role, profile in sorted(values):
        role_label = {
            "primary_center": "Primary center",
            "participant": "Participating site",
        }.get(str(role), _readable(role).capitalize())
        lines.append(
            f"  {role_label} in reactant {component_index + 1}: {profile}"
        )
    return lines


def format_core_graph_analysis(analysis: Any) -> str:
    """Render the canonical reaction core, R ports, and site annotations."""
    reaction_label = getattr(analysis, "reaction_label", None)
    lines = [
        f"Reaction: {getattr(reaction_label, 'text', None) or '-'}",
    ]
    core = getattr(analysis, "reaction_core", None)
    if core is None:
        lines.append("Core: unavailable (verified or consensus edits required)")
        return "\n".join(lines)
    evidence_status = {
        "verified": "Verified",
        "inferred": "Inferred",
        "hypothesis": "Consensus hypothesis",
        "external": "External mapping; review required",
    }.get(str(core.evidence_status), _readable(core.evidence_status).capitalize())
    lines.append(
        f"Evidence: {evidence_status} from {_evidence_text(core.evidence)} "
        f"({core.confidence:.0%} confidence)"
    )
    quality = core.quality
    cautions = tuple(quality.review_reasons) + tuple(quality.blocking_reasons)
    if cautions:
        lines.append(
            "Review note: "
            + "; ".join(_quality_note(value) for value in cautions)
        )
    lines.append("Bond changes:" if len(core.events) == 1 else "Reaction events:")
    for index, event in enumerate(core.events, start=1):
        changes = "; ".join(_edit_text(token) for token in event.edit_tokens)
        prefix = f"{index}. " if len(core.events) > 1 else ""
        lines.append(f"  {prefix}{changes}")
    lines.extend(_r_group_analysis_lines(analysis))
    lines.extend(_active_site_profile_lines(analysis))
    return "\n".join(lines)


def format_reaction_additional_analysis(analysis: Any) -> str:
    """Render a compact current-contract summary without legacy duplicates."""
    observation = getattr(analysis, "observation", None)
    signature = getattr(analysis, "reaction_signature", None)
    interpretation = getattr(analysis, "interpretation", None)
    completeness = getattr(analysis, "reaction_completeness", None)
    lines = [
        f"Status: {'valid' if getattr(analysis, 'valid', False) else 'invalid'}",
        "Observation evidence: "
        f"{_evidence_text(getattr(observation, 'evidence_quality', None)) or '-'}",
        f"Generic signature: {'available' if signature is not None else 'unavailable'}",
        "Transformation: "
        f"{_annotation_text(getattr(signature, 'transformation_class', None)) or 'unknown'}",
        "Completeness: "
        f"{_readable(getattr(completeness, 'status', None)) or 'unknown'}",
    ]
    primary_pattern = getattr(interpretation, "primary_pattern_id", None)
    if primary_pattern:
        lines.append(f"Primary pattern: {_annotation_text(primary_pattern)}")
    compatible = tuple(
        getattr(interpretation, "compatible_named_families", ()) or ()
    )
    if compatible:
        lines.append(
            "Compatible annotations: "
            + ", ".join(_annotation_text(value) for value in compatible)
        )
    named_family = getattr(interpretation, "named_family", None)
    if named_family:
        lines.append(f"Resolved annotation: {_annotation_text(named_family)}")
    warnings = tuple(
        value
        for value in (getattr(analysis, "warnings", ()) or ())
        if value != "INFERRED_GLOBAL_ATOM_CORRESPONDENCE"
    )
    if warnings:
        lines.append(
            "Warnings: "
            + ", ".join(_warning_text(value) for value in warnings)
        )
    error = getattr(analysis, "error", None)
    if error:
        lines.append(f"Error: {error}")
    return "\n".join(lines)


class ReactiveTaxonomyWindow(QtWidgets.QMainWindow):
    """Auto-detecting molecule and reaction featurization window."""

    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("reactiveTaxonomyWindow")
        self.setFont(QtGui.QFont("Segoe UI", 9))
        self.setWindowTitle("Reactive Taxonomy Featurizer")
        self.resize(900, 620)
        self._mapping_provider: RxnMapperProvider | None = None
        self._last_analysis: object | None = None
        self._last_kind: InputKind | None = None
        self._last_input_text = ""
        self._build_ui()
        self._connect_signals()
        self._apply_style()
        self._update_detected_kind("")

    def _build_ui(self) -> None:
        central = QtWidgets.QWidget()
        central.setObjectName("centralPanel")
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)
        layout.setContentsMargins(28, 24, 28, 24)
        layout.setSpacing(14)

        title = QtWidgets.QLabel("Molecule & reaction featurizer")
        title.setObjectName("title")
        layout.addWidget(title)

        description = QtWidgets.QLabel(
            "Paste a molecular SMILES or reaction SMILES. Input containing a "
            "reaction separator (>) is analyzed as a reaction automatically."
        )
        description.setObjectName("description")
        description.setWordWrap(True)
        layout.addWidget(description)

        input_row = QtWidgets.QHBoxLayout()
        input_row.setSpacing(10)
        self.input_edit = QtWidgets.QLineEdit()
        self.input_edit.setObjectName("smilesInput")
        self.input_edit.setClearButtonEnabled(True)
        self.input_edit.setPlaceholderText(
            "SMILES or reactants>>product"
        )
        self.input_edit.setAccessibleName("Molecule or reaction SMILES")
        input_row.addWidget(self.input_edit, 1)

        self.analyze_button = QtWidgets.QPushButton("Analyze")
        self.analyze_button.setObjectName("analyzeInput")
        self.analyze_button.setDefault(True)
        input_row.addWidget(self.analyze_button)
        layout.addLayout(input_row)

        controls = QtWidgets.QHBoxLayout()
        controls.setSpacing(8)
        self.kind_label = QtWidgets.QLabel()
        self.kind_label.setObjectName("detectedKind")
        controls.addWidget(self.kind_label)
        self.use_rxnmapper_check = QtWidgets.QCheckBox(
            "Use RXNMapper for unresolved or ambiguous reactions"
        )
        self.use_rxnmapper_check.setObjectName("useRxnMapper")
        self.use_rxnmapper_check.setChecked(True)
        self.use_rxnmapper_check.setToolTip(
            "Checked by default. Supplied maps and resolved internal evidence "
            "still take precedence; generated mapping remains review evidence."
        )
        controls.addWidget(self.use_rxnmapper_check)
        self.force_core_mapping_check = QtWidgets.QCheckBox(
            "Map resolved reactions for minimized graphic"
        )
        self.force_core_mapping_check.setObjectName("forceCoreMapping")
        self.force_core_mapping_check.setChecked(False)
        self.force_core_mapping_check.setToolTip(
            "Optional and slower. Runs RXNMapper even when internal evidence "
            "already resolves the reaction, so a mapped graphical core can be "
            "drawn. The mapped result remains review evidence."
        )
        controls.addWidget(self.force_core_mapping_check)
        controls.addStretch(1)

        self.reaction_example_button = QtWidgets.QPushButton("Reaction example")
        self.reaction_example_button.setObjectName("reactionExample")
        controls.addWidget(self.reaction_example_button)

        self.molecule_example_button = QtWidgets.QPushButton("Molecule example")
        self.molecule_example_button.setObjectName("moleculeExample")
        controls.addWidget(self.molecule_example_button)

        self.copy_button = QtWidgets.QPushButton("Copy result")
        self.copy_button.setObjectName("copyResult")
        self.copy_button.setEnabled(False)
        controls.addWidget(self.copy_button)
        layout.addLayout(controls)

        result_panel = QtWidgets.QWidget()
        result_layout = QtWidgets.QHBoxLayout(result_panel)
        result_layout.setContentsMargins(0, 0, 0, 0)
        result_layout.setSpacing(12)

        analysis_column = QtWidgets.QWidget()
        analysis_layout = QtWidgets.QVBoxLayout(analysis_column)
        analysis_layout.setContentsMargins(0, 0, 0, 0)
        analysis_layout.setSpacing(4)
        self.core_analysis_heading = QtWidgets.QLabel("Core graph analysis")
        self.core_analysis_heading.setObjectName("coreGraphAnalysisHeading")
        analysis_layout.addWidget(self.core_analysis_heading)

        self.core_analysis_output = QtWidgets.QPlainTextEdit()
        self.core_analysis_output.setObjectName("coreGraphAnalysisOutput")
        self.core_analysis_output.setReadOnly(True)
        self.core_analysis_output.setPlaceholderText(
            "Reaction-core events, R-group attachment profiles, and active-site "
            "steric/electronic context appear here."
        )
        self.core_analysis_output.setLineWrapMode(
            QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth
        )
        self.core_analysis_output.setMinimumHeight(240)
        self.core_analysis_output.setMaximumHeight(340)
        analysis_layout.addWidget(self.core_analysis_output)

        self.additional_analysis_heading = QtWidgets.QLabel(
            "Additional analysis"
        )
        self.additional_analysis_heading.setObjectName(
            "additionalAnalysisHeading"
        )
        analysis_layout.addWidget(self.additional_analysis_heading)

        self.output = QtWidgets.QPlainTextEdit()
        self.output.setObjectName("analysisOutput")
        self.output.setReadOnly(True)
        self.output.setPlaceholderText("Featurization results appear here.")
        self.output.setLineWrapMode(
            QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth
        )
        fixed_font = QtGui.QFontDatabase.systemFont(
            QtGui.QFontDatabase.SystemFont.FixedFont
        )
        fixed_font.setPointSize(10)
        self.output.setFont(fixed_font)
        analysis_layout.addWidget(self.output)

        graph_column = QtWidgets.QWidget()
        graph_layout = QtWidgets.QVBoxLayout(graph_column)
        graph_layout.setContentsMargins(0, 0, 0, 0)
        graph_layout.setSpacing(4)
        graph_header = QtWidgets.QHBoxLayout()
        self.graph_heading = QtWidgets.QLabel("Structure graph")
        self.graph_heading.setObjectName("structureGraphHeading")
        graph_header.addWidget(self.graph_heading)
        graph_header.addStretch(1)
        graph_header.addWidget(QtWidgets.QLabel("Drawing style"))
        self.render_style_combo = QtWidgets.QComboBox()
        self.render_style_combo.setObjectName("renderStylePreset")
        for preset_id, label in available_render_presets():
            self.render_style_combo.addItem(label, preset_id)
        default_index = self.render_style_combo.findData(DEFAULT_RENDER_PRESET)
        if default_index >= 0:
            self.render_style_combo.setCurrentIndex(default_index)
        self.render_style_combo.setToolTip(
            "Current is the default project drawing style. Other presets can "
            "be selected without rerunning the chemistry analysis."
        )
        graph_header.addWidget(self.render_style_combo)
        graph_layout.addLayout(graph_header)
        self.full_structure_panel = QtWidgets.QGroupBox("Full structure")
        self.full_structure_panel.setObjectName("fullStructurePanel")
        full_structure_layout = QtWidgets.QVBoxLayout(
            self.full_structure_panel
        )
        full_structure_layout.setContentsMargins(6, 10, 6, 6)
        self.structure_image_label = StructureImageLabel(
            placeholder="Reaction or compound graph will appear here.",
            object_name="featurizedStructureGraph",
            minimum_height=220,
        )
        full_structure_layout.addWidget(self.structure_image_label)
        graph_layout.addWidget(self.full_structure_panel, 1)

        self.minimized_panel = QtWidgets.QGroupBox("Minimized reaction")
        self.minimized_panel.setObjectName("minimizedReactionPanel")
        minimized_layout = QtWidgets.QVBoxLayout(self.minimized_panel)
        minimized_layout.setContentsMargins(6, 10, 6, 6)
        minimized_layout.setSpacing(4)
        self.core_image_label = StructureImageLabel(
            placeholder="Mapped minimized reaction will appear here.",
            object_name="minimizedReactionGraphic",
            minimum_height=190,
        )
        minimized_layout.addWidget(self.core_image_label, 1)
        self.core_graphic_note = QtWidgets.QLabel(
            "A mapped reaction core is required."
        )
        self.core_graphic_note.setObjectName("coreGraphicNote")
        self.core_graphic_note.setWordWrap(True)
        minimized_layout.addWidget(self.core_graphic_note)
        graph_layout.addWidget(self.minimized_panel, 1)

        result_layout.addWidget(analysis_column, stretch=9)
        result_layout.addWidget(graph_column, stretch=11)
        layout.addWidget(result_panel, 1)

        self.status_label = QtWidgets.QLabel("Ready")
        self.status_label.setObjectName("statusLabel")
        layout.addWidget(self.status_label)

        QtGui.QShortcut(
            QtGui.QKeySequence("Ctrl+Return"),
            self,
            activated=self.analyze,
        )
        QtGui.QShortcut(
            QtGui.QKeySequence("Ctrl+Enter"),
            self,
            activated=self.analyze,
        )

    def _connect_signals(self) -> None:
        self.input_edit.textChanged.connect(self._update_detected_kind)
        self.use_rxnmapper_check.toggled.connect(
            self.force_core_mapping_check.setEnabled
        )
        self.input_edit.returnPressed.connect(self.analyze)
        self.analyze_button.clicked.connect(self.analyze)
        self.reaction_example_button.clicked.connect(
            lambda: self._load_example(REACTION_EXAMPLE)
        )
        self.molecule_example_button.clicked.connect(
            lambda: self._load_example(MOLECULE_EXAMPLE)
        )
        self.copy_button.clicked.connect(self.copy_result)
        self.render_style_combo.currentIndexChanged.connect(
            self._rerender_last_structure
        )

    def _apply_style(self) -> None:
        """Add converter-style accents while preserving the native dark palette."""
        self.setStyleSheet(
            """
            QLabel#title {
                font-size: 22px;
                font-weight: 600;
            }
            QLabel#description {
                color: #d0d7de;
                font-size: 13px;
            }
            QPushButton#analyzeInput {
                background-color: #0078d7;
                color: white;
                font-weight: 700;
                border: none;
                border-radius: 6px;
                padding: 8px 13px;
                padding-left: 20px;
                padding-right: 20px;
            }
            QPushButton#analyzeInput:hover {
                background-color: #1689e5;
            }
            QPushButton#analyzeInput:disabled {
                background-color: #355b78;
                color: #aebdcc;
            }
            QLabel#detectedKind {
                color: #6cb6ff;
                font-weight: 600;
            }
            QLabel#statusLabel {
                color: #c5ced8;
            }
            """
        )

    @QtCore.pyqtSlot(str)
    def _update_detected_kind(self, text: str) -> None:
        value = text.strip()
        if not value:
            self.kind_label.setText("Detected: waiting for input")
            return
        kind = detect_input_kind(value)
        self.kind_label.setText(f"Detected: {kind}")

    def _load_example(self, text: str) -> None:
        self.input_edit.setText(text)
        self.input_edit.setFocus()
        self.input_edit.selectAll()

    @QtCore.pyqtSlot()
    def analyze(self) -> None:
        """Analyze the current input and display the concise CLI-equivalent view."""
        self.analyze_button.setEnabled(False)
        self.status_label.setText("Analyzing…")
        QtWidgets.QApplication.processEvents(
            QtCore.QEventLoop.ProcessEventsFlag.ExcludeUserInputEvents
        )
        try:
            input_text = self.input_edit.text()
            requested_kind = detect_input_kind(input_text.strip())
            mapping_provider = None
            if (
                requested_kind == "reaction"
                and self.use_rxnmapper_check.isChecked()
            ):
                if not RxnMapperProvider.is_available():
                    raise RuntimeError(
                        "RXNMapper is not installed. Run "
                        "'python -m pip install -r requirements-mapping.txt' "
                        "or clear the RXNMapper checkbox."
                    )
                if self._mapping_provider is None:
                    self._mapping_provider = RxnMapperProvider()
                mapping_provider = self._mapping_provider
            kind, analysis, assessment = featurize_text(
                input_text,
                mapping_provider=mapping_provider,
                force_resolved_mapping=(
                    requested_kind == "reaction"
                    and self.force_core_mapping_check.isChecked()
                ),
            )
            heading = f"{kind.upper()} FEATURIZATION"
            mapping_summary = ""
            if kind == "reaction":
                if assessment is None:
                    mapping_summary = "\nRXNMapper: disabled"
                else:
                    result = assessment.mapping_result
                    confidence = (
                        f"{result.mapper_confidence:.3f}"
                        if result is not None
                        and result.mapper_confidence is not None
                        else "not run"
                    )
                    mapping_summary = (
                        f"\nRXNMapper: {_mapping_status_text(assessment.status)}"
                    )
                    if result is not None:
                        mapping_summary += f"\nMapping confidence: {confidence}"
            if kind == "reaction":
                self.core_analysis_output.setPlainText(
                    format_core_graph_analysis(analysis)
                )
                additional = format_reaction_additional_analysis(analysis)
            else:
                self.core_analysis_output.setPlainText(
                    "Core graph analysis applies to reaction SMILES."
                )
                additional = format_concise_analysis(analysis)
            self.output.setPlainText(
                f"{heading}{mapping_summary}\n\n{additional}"
            )
            self._last_analysis = analysis
            self._last_kind = kind
            self._last_input_text = self.input_edit.text().strip()
            self._render_structure(
                kind,
                self.input_edit.text().strip(),
                analysis=analysis,
            )
            valid = bool(getattr(analysis, "valid", False))
            state = "valid" if valid else "invalid"
            self.status_label.setText(
                f"Complete · {kind} input · {state}"
                + (
                    " · RXNMapper on"
                    if kind == "reaction"
                    and self.use_rxnmapper_check.isChecked()
                    else ""
                )
            )
            self.copy_button.setEnabled(True)
        except Exception as exc:
            self._last_analysis = None
            self._last_kind = None
            self._last_input_text = ""
            self.output.setPlainText(f"Unable to analyze input.\n\n{exc}")
            self.core_analysis_output.setPlainText(
                "Core graph analysis unavailable."
            )
            self.graph_heading.setText("Structure graph")
            self.structure_image_label.clear_image(
                "Structure graph unavailable."
            )
            self.core_image_label.clear_image(
                "Minimized reaction graphic unavailable."
            )
            self.core_graphic_note.setText(str(exc))
            self.status_label.setText("Analysis failed")
            self.copy_button.setEnabled(True)
        finally:
            self.analyze_button.setEnabled(True)

    def _render_structure(
        self,
        kind: InputKind,
        text: str,
        *,
        analysis: object,
    ) -> None:
        """Render the full graph and, when available, its minimized core."""
        self.graph_heading.setText(
            "Reaction graph" if kind == "reaction" else "Compound graph"
        )
        render_preset = str(self.render_style_combo.currentData() or "current")
        try:
            if kind == "reaction":
                drawing = render_reaction_image_bytes(
                    text,
                    size=REACTION_IMAGE_SIZE,
                    image_format="svg",
                    render_preset=render_preset,
                )
            else:
                drawing = render_molecule_image_bytes(
                    text,
                    size=MOLECULE_IMAGE_SIZE,
                    image_format="svg",
                    render_preset=render_preset,
                )
        except (RuntimeError, ValueError) as exc:
            self.structure_image_label.clear_image(
                f"{self.graph_heading.text()} unavailable."
            )
            self.structure_image_label.setToolTip(str(exc))
            return
        if not self.structure_image_label.set_image_bytes(
            drawing,
            trim_white_space=True,
            max_upscale=(6.0 if render_preset == "acs_1996_compact" else None),
        ):
            self.structure_image_label.setToolTip(
                "The renderer returned an unsupported image."
            )
            return
        self.structure_image_label.setToolTip(text)
        if kind != "reaction":
            self.core_image_label.clear_image(
                "Reaction minimization applies only to reactions."
            )
            self.core_graphic_note.setText(
                "Enter a reaction SMILES to generate a minimized graphic."
            )
            return
        core = getattr(analysis, "reaction_core", None)
        if core is None:
            self.core_image_label.clear_image(
                "Mapped reaction core unavailable."
            )
            self.core_graphic_note.setText(
                "Supply atom mapping, use RXNMapper for an unresolved "
                "reaction, or enable resolved-reaction mapping."
            )
            return
        try:
            projection = build_reaction_display_projection(
                reaction_render_context_from_analysis(analysis)
            )
            graphic = build_reaction_display_graphic(
                projection,
                size=REACTION_IMAGE_SIZE,
                render_preset=render_preset,
            )
        except (RuntimeError, ValueError) as exc:
            self.core_image_label.clear_image(
                "Unable to render minimized reaction."
            )
            self.core_graphic_note.setText(str(exc))
            return
        if not self.core_image_label.set_image_bytes(
            graphic.image_bytes,
            trim_white_space=True,
            max_upscale=(6.0 if render_preset == "acs_1996_compact" else None),
        ):
            self.core_graphic_note.setText(
                "The minimized renderer returned an unsupported image."
            )
            return
        removed_substituent_count = sum(
            component.removed_substituent_count
            for component in projection.reactants + projection.products
        )
        r_group_count = sum(
            component.r_group_count
            for component in projection.reactants + projection.products
        )
        label_by_index = {
            int(value.placeholder_index): str(value.display_label)
            for value in projection.substituents
            if value.placeholder_index is not None and value.display_label
        }
        placeholder_labels = tuple(
            label_by_index[index] for index in sorted(label_by_index)
        )
        profile_source = tuple(
            value
            for value in projection.substituents
            if value.side == "reactant"
            and value.placeholder_index is not None
            and value.display_label
        ) or tuple(
            value
            for value in projection.substituents
            if value.side == "product"
            and value.placeholder_index is not None
            and value.display_label
        )
        profile_legend = []
        seen_profile_labels = set()
        for value in sorted(
            profile_source,
            key=lambda item: int(item.placeholder_index or 0),
        ):
            label = str(value.display_label)
            if label in seen_profile_labels:
                continue
            seen_profile_labels.add(label)
            profile = value.substituent_profile
            descriptors = []
            if profile.carbon_substitution not in {
                "not_carbon",
                "not_applicable",
                "unresolved",
            }:
                descriptors.append(profile.carbon_substitution)
            descriptors.append(profile.base_class.replace("_", " "))
            descriptors.extend(
                name
                for name, present in (
                    ("benzylic", profile.benzylic),
                    ("allylic", profile.allylic),
                    ("propargylic", profile.propargylic),
                )
                if present
            )
            profile_legend.append(f"{label}={' '.join(descriptors)}")
        hidden_aromatic = tuple(
            sorted(
                {
                    f"{relation.positional_relation} "
                    f"{value.fragment_smiles}"
                    for value in projection.substituents
                    if value.boundary_action == "aromatic_hydrogen_cap"
                    for relation in value.aromatic_relations
                }
            )
        )
        connector_source = tuple(
            value
            for value in projection.connectors
            if value.side == "reactant"
        ) or tuple(
            value
            for value in projection.connectors
            if value.side == "product"
        )
        hidden_connectors = tuple(
            f"{value.display_label} "
            f"{'⋯'.join(value.port_display_labels)} "
            + (
                f"({len(value.shortest_path_atom_indices)} hidden-path atoms)"
                if value.shortest_path_atom_indices
                else f"({len(value.placeholder_indices)} shared ports)"
            )
            for value in connector_source
        )
        note = (
            f"Display-only minimum: {projection.minimum_reaction_smiles}; "
            f"{projection.evidence_status} evidence; "
            f"confidence {projection.confidence:.3f}; "
            f"R groups: {r_group_count}; "
            f"removed unchanged aromatic substituents: "
            f"{removed_substituent_count}"
        )
        if placeholder_labels:
            note += f"; labels: {', '.join(placeholder_labels)}"
        if profile_legend:
            note += f"; R profiles: {', '.join(profile_legend)}"
        if hidden_aromatic:
            note += f"; hidden aromatic: {', '.join(hidden_aromatic)}"
        if hidden_connectors:
            note += f"; hidden connectors: {', '.join(hidden_connectors)}"
        if projection.annotation:
            note += f"; {projection.annotation}"
        if projection.evidence_status == "external":
            note += "; external mapping requires expert review"
        self.core_graphic_note.setText(note)
        self.core_image_label.setToolTip(note)

    @QtCore.pyqtSlot()
    def _rerender_last_structure(self) -> None:
        """Apply a newly selected drawing style without rerunning chemistry."""
        current_text = self.input_edit.text().strip()
        if (
            self._last_analysis is None
            or self._last_kind is None
            or current_text != self._last_input_text
        ):
            return
        self._render_structure(
            self._last_kind,
            current_text,
            analysis=self._last_analysis,
        )

    @QtCore.pyqtSlot()
    def copy_result(self) -> None:
        """Copy the visible analysis to the system clipboard."""
        sections = (
            self.core_analysis_output.toPlainText().strip(),
            self.output.toPlainText().strip(),
        )
        QtWidgets.QApplication.clipboard().setText(
            "\n\n".join(section for section in sections if section)
        )
        self.status_label.setText("Result copied to clipboard")


def _show_main_window(window: ReactiveTaxonomyWindow) -> None:
    """Show the featurizer in its requested initial state."""
    window.showMaximized()


def main(argv: Sequence[str] | None = None) -> int:
    """Launch the Qt6 featurization application."""
    application = QtWidgets.QApplication(
        list(argv) if argv is not None else sys.argv
    )
    application.setApplicationName("Reactive Taxonomy Featurizer")
    window = ReactiveTaxonomyWindow()
    _show_main_window(window)
    return application.exec()


if __name__ == "__main__":
    raise SystemExit(main())
