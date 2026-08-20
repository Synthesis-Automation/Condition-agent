"""Qt-free presentation service for molecule and reaction featurization."""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, Literal

from reactive_taxonomy import (
    AtomMappingProvider,
    analyze_molecule,
    analyze_reaction_with_external_mapping,
    build_reaction_display_projection,
    featurize_reaction,
    reaction_render_context_from_analysis,
)
from visualization import build_reaction_display_graphic


InputKind = Literal["molecule", "reaction"]


def detect_input_kind(text: str) -> InputKind:
    """Route any reaction separator to the reaction parser."""

    return "reaction" if ">" in text else "molecule"


def _molecule_overview(analysis: Any) -> Dict[str, Any]:
    components = tuple(analysis.components)
    atoms = tuple(
        atom for component in components for atom in component.atoms
    )
    bonds = tuple(
        bond for component in components for bond in component.bonds
    )
    elements = Counter(str(atom.element) for atom in atoms)
    return {
        "canonical_smiles": analysis.canonical_smiles,
        "component_count": len(components),
        "atom_count": len(atoms),
        "bond_count": len(bonds),
        "element_counts": dict(sorted(elements.items())),
        "motif_count": len(analysis.motifs),
        "reactive_site_count": len(analysis.reactive_site_hypotheses),
    }


def _molecule_motifs(analysis: Any) -> list[Dict[str, Any]]:
    return [
        {
            "motif_id": value.motif_id,
            "label": value.chemist_label,
            "component_index": value.component_index,
            "atom_indices": list(value.atom_indices),
            "tags": list(value.tags),
            "confidence": value.confidence,
        }
        for value in analysis.motifs
    ]


def _molecule_sites(analysis: Any) -> list[Dict[str, Any]]:
    environment_by_hypothesis = {
        value.hypothesis_id: value
        for value in analysis.reactive_site_environments
    }
    sites = []
    for value in analysis.reactive_site_hypotheses:
        environment = environment_by_hypothesis.get(value.hypothesis_id)
        profile = (
            environment.reactivity_profile
            if environment is not None
            else None
        )
        sites.append(
            {
                "hypothesis_id": value.hypothesis_id,
                "label": value.chemist_label,
                "site_type": value.site_type,
                "availability": value.availability,
                "component_index": value.component_index,
                "atom_indices": list(value.atom_indices),
                "context_kind": (
                    profile.context_kind if profile is not None else None
                ),
            }
        )
    return sites


def _reaction_motifs(analysis: Any) -> list[Dict[str, Any]]:
    values = []
    seen = set()
    for side_name, components in (
        ("reactant", analysis.reactants),
        ("agent", analysis.agents),
        ("product", analysis.products),
    ):
        for component in components:
            for motif in component.molecule_analysis.motifs:
                key = (
                    side_name,
                    component.component_index,
                    motif.motif_id,
                    tuple(motif.atom_indices),
                )
                if key in seen:
                    continue
                seen.add(key)
                values.append(
                    {
                        "side": side_name,
                        "component_index": component.component_index,
                        "motif_id": motif.motif_id,
                        "label": motif.chemist_label,
                        "atom_indices": list(motif.atom_indices),
                        "tags": list(motif.tags),
                        "confidence": motif.confidence,
                    }
                )
    return values


def _reaction_sites(analysis: Any) -> list[Dict[str, Any]]:
    """Project reactant-side site hypotheses with their local contexts."""

    sites = []
    for component in analysis.reactants:
        molecule = component.molecule_analysis
        environment_by_hypothesis = {
            value.hypothesis_id: value
            for value in molecule.reactive_site_environments
        }
        for value in molecule.reactive_site_hypotheses:
            environment = environment_by_hypothesis.get(value.hypothesis_id)
            profile = (
                environment.reactivity_profile
                if environment is not None
                else None
            )
            sites.append(
                {
                    "hypothesis_id": value.hypothesis_id,
                    "label": value.chemist_label,
                    "site_type": value.site_type,
                    "availability": value.availability,
                    "side": "reactant",
                    "component_index": component.component_index,
                    "atom_indices": list(value.atom_indices),
                    "context_kind": (
                        profile.context_kind if profile is not None else None
                    ),
                }
            )
    return sites


def _reaction_core_summary(
    core: Any,
    completeness: Any = None,
) -> Dict[str, Any] | None:
    if core is None:
        return None
    incomplete = (
        completeness is not None and completeness.status == "incomplete"
    )
    blocking_reasons = set(core.quality.blocking_reasons)
    warnings = set(core.warnings)
    if incomplete:
        blocking_reasons.add("incomplete_product_atom_provenance")
        warnings.add("REACTION_CORE_DISPLAY_BLOCKED_INCOMPLETE_REACTION")
    return {
        "core_id": core.core_id,
        "evidence": core.evidence,
        "evidence_status": core.evidence_status,
        "confidence": core.confidence,
        "active_atom_count": core.active_atom_count,
        "event_count": core.event_count,
        "events": [
            {
                "event_id": event.event_id,
                "edit_tokens": list(event.edit_tokens),
                "reactant_component_indices": list(
                    event.reactant_component_indices
                ),
                "product_component_indices": list(
                    event.product_component_indices
                ),
            }
            for event in core.events
        ],
        "quality": {
            "status": "blocked" if incomplete else core.quality.status,
            "review_reasons": list(core.quality.review_reasons),
            "blocking_reasons": sorted(blocking_reasons),
            "checked_edit_fraction": core.quality.checked_edit_fraction,
            "active_atom_mapping_coverage": (
                core.quality.active_atom_mapping_coverage
            ),
        },
        "warnings": sorted(warnings),
    }


def _reaction_overview(analysis: Any) -> Dict[str, Any]:
    signature = analysis.reaction_signature
    completeness = analysis.reaction_completeness
    label = analysis.reaction_label
    return {
        "reaction_label": label.text if label is not None else None,
        "label_confidence": label.confidence if label is not None else None,
        "reactant_count": len(analysis.reactants),
        "agent_count": len(analysis.agents),
        "product_count": len(analysis.products),
        "edit_archetype": analysis.edit_archetype,
        "transformation_class": analysis.transformation_class,
        "named_family": analysis.named_family,
        "compatible_named_families": list(
            analysis.compatible_named_families
        ),
        "evidence_quality": analysis.evidence_quality,
        "signature_id": (
            signature.signature_id if signature is not None else None
        ),
        "completeness_status": (
            completeness.status if completeness is not None else None
        ),
    }


def _reaction_partners(analysis: Any) -> list[Dict[str, Any]]:
    signature = analysis.reaction_signature
    if signature is None:
        return []
    return [
        {
            "role": value.role,
            "component_index": value.component_index,
            "label": value.chemist_label,
            "role_confidence": value.role_confidence,
            "anchor_contexts": list(value.anchor_contexts),
        }
        for value in signature.partners
    ]


def analyze_features(
    text: str,
    *,
    mapping_provider: AtomMappingProvider | None = None,
    force_resolved_mapping: bool = False,
) -> Dict[str, Any]:
    """Return a compact UI projection plus the canonical nested analysis."""

    value = str(text or "").strip()
    if not value:
        raise ValueError("Enter a molecule or reaction SMILES.")
    kind = detect_input_kind(value)
    if kind == "molecule":
        analysis = analyze_molecule(value)
        return {
            "input_kind": kind,
            "input_smiles": value,
            "valid": analysis.valid,
            "schema_version": analysis.schema_version,
            "overview": _molecule_overview(analysis),
            "motifs": _molecule_motifs(analysis),
            "reactive_sites": _molecule_sites(analysis),
            "reaction_core": None,
            "partners": [],
            "mapping": None,
            "core_graphic_svg": None,
            "core_projection": None,
            "warnings": list(analysis.structure.warnings),
            "error": analysis.structure.error,
            "analysis": analysis.to_dict(),
        }

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
    analysis = assessment.analysis if assessment is not None else base_analysis
    core_graphic_svg = None
    core_projection = None
    core_graphic_error = None
    if (
        analysis.reaction_core is not None
        and analysis.reaction_core.quality.status != "blocked"
        and (
            analysis.reaction_completeness is None
            or analysis.reaction_completeness.status != "incomplete"
        )
    ):
        try:
            projection = build_reaction_display_projection(
                reaction_render_context_from_analysis(analysis)
            )
            graphic = build_reaction_display_graphic(
                projection,
                size=(760, 190),
                image_format="svg",
                render_preset="web_consistent",
            )
            core_graphic_svg = graphic.image_bytes.decode("utf-8")
            core_projection = {
                "minimum_reaction_smiles": projection.minimum_reaction_smiles,
                "reaction_scope": projection.reaction_scope,
                "evidence_status": projection.evidence_status,
                "confidence": projection.confidence,
                "annotation": projection.annotation,
                "warnings": list(projection.warnings),
            }
        except (RuntimeError, ValueError) as exc:
            core_graphic_error = str(exc)
    elif analysis.reaction_core is not None:
        core_graphic_error = (
            "Reaction core graphic withheld because the structural evidence "
            "is blocked; review the listed core reasons and reaction warnings."
        )

    return {
        "input_kind": kind,
        "input_smiles": value,
        "valid": analysis.valid,
        "schema_version": analysis.schema_version,
        "overview": _reaction_overview(analysis),
        "motifs": _reaction_motifs(analysis),
        "reactive_sites": _reaction_sites(analysis),
        "reaction_core": _reaction_core_summary(
            analysis.reaction_core,
            analysis.reaction_completeness,
        ),
        "partners": _reaction_partners(analysis),
        "mapping": (
            assessment.to_provenance_dict()
            if assessment is not None
            else None
        ),
        "core_graphic_svg": core_graphic_svg,
        "core_graphic_error": core_graphic_error,
        "core_projection": core_projection,
        "warnings": list(analysis.warnings),
        "error": analysis.error,
        "analysis": analysis.to_dict(),
    }


__all__ = ["InputKind", "analyze_features", "detect_input_kind"]
