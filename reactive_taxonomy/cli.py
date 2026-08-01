"""Integrated command-line tester for standalone reactive-taxonomy features."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from . import (
    featurize_molecule,
    featurize_reaction,
    render_reactivity_profile,
    validate_taxonomy,
)


_MOLECULE_SITE_TYPES = (
    "leaving_group",
    "pronucleophile_XH",
    "transfer_group",
    "electrophilic_center",
    "aromatic_CH",
    "unsaturated_bond",
    "dipolar_group",
    "heteroatom_bond",
)

_SELF_TEST_MOLECULES = (
    "Brc1ccc(N)cc1C#N",
    "CC(=O)Cl",
    "CS(=O)(=O)Cl",
    "OB(O)c1ccccc1",
    "CNC",
)

_SELF_TEST_REACTIONS = (
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
    "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1",
    "Brc1ccccc1.Sc1ccccc1>>c1ccc(Sc2ccccc2)cc1",
)

_REACTION_CORE_CSV_FIELDS = (
    "reaction_core_available",
    "reaction_core_id",
    "reaction_core_evidence_status",
    "reaction_core_exact_key",
    "reaction_core_typed_key",
    "reaction_core_shape_key",
    "reaction_core_center_transition_key",
    "reaction_core_mapping_equivalence_key",
    "reaction_core_quality_status",
    "reaction_core_quality_reasons",
    "reaction_core_bond_changes",
    "reaction_core_state_changes",
    "reaction_core_retained_context",
    "reaction_core_departing_context",
    "reaction_core_appearing_context",
    "reaction_core_motif_key",
    "reaction_core_general_equation",
    "reaction_core_limiter",
    "reaction_core_event_count",
    "reaction_core_primary_center_count",
    "reaction_core_remote_classes",
    "reaction_core_remote_subgraphs",
    "reaction_core_warnings",
)

_REACTION_RING_CSV_FIELDS = (
    "reaction_display_label_detailed",
    "reaction_display_source",
    "reaction_display_status",
    "reaction_display_confidence",
    "formed_ring_sizes",
    "ring_count_delta",
    "ring_change_count",
    "ring_changes_json",
)


def _json_dump(value: Any, *, compact: bool = False) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        indent=None if compact else 2,
    )


def _joined(values: Iterable[Any]) -> str:
    return ", ".join(str(value) for value in values) or "none"


def _reaction_partners(result: Any) -> tuple[Any, ...]:
    """Return interpreted partners, with generic partners as a fallback."""
    interpretation = getattr(result, "interpretation", None)
    interpreted_partners = (
        getattr(interpretation, "partners", ())
        if interpretation is not None
        else ()
    )
    environment = (
        getattr(interpretation, "family_environment", None)
        if interpretation is not None
        else None
    )
    if environment is None:
        environment = getattr(result, "family_environment", None)
    signature = getattr(result, "reaction_signature", None)
    if interpreted_partners:
        partners = interpreted_partners
    elif environment and environment.partners:
        partners = environment.partners
    elif signature and signature.partners:
        partners = signature.partners
    else:
        partners = ()
    role_order = {"electrophile": 0, "nucleophile": 1, "transfer_partner": 2}
    return tuple(sorted(
        partners,
        key=lambda partner: (
            role_order.get(getattr(partner, "role", None), 99),
            getattr(partner, "role", None) or "unassigned",
            int(getattr(partner, "component_index", -1)),
        ),
    ))


def _partner_value(partner: Any, name: str, default: Any = None) -> Any:
    if isinstance(partner, Mapping):
        return partner.get(name, default)
    return getattr(partner, name, default)


def _partner_analysis(result: Any) -> str:
    """Render canonical context-aware profiles in one compact field."""
    partners = _reaction_partners(result)

    def role(partner: Any) -> str:
        return str(_partner_value(partner, "role") or "unassigned")

    return "; ".join(
        f"{role(partner)}="
        f"{render_reactivity_profile(_partner_value(partner, 'reactivity_profile'))}"
        for partner in partners
    )


def _counted_elements(values: Mapping[str, Any]) -> str:
    """Render a deterministic elemental count for chemist-facing summaries."""
    return ", ".join(
        f"{element} × {int(count)}"
        for element, count in sorted(values.items())
        if int(count)
    ) or "none"


def _bond_inventory_summary(tokens: Iterable[str]) -> str:
    """Render non-corresponded bond-inventory deltas without claiming edits."""
    counts = Counter(str(token) for token in tokens)
    action_order = {"gained": 0, "lost": 1, "changed": 2}

    def sort_key(item: tuple[str, int]) -> tuple[Any, ...]:
        token, _ = item
        parts = token.split(":", 2)
        action = parts[0] if parts else ""
        return (action_order.get(action, 99), token)

    rendered = []
    for token, count in sorted(counts.items(), key=sort_key):
        parts = token.split(":", 2)
        if len(parts) != 3:
            rendered.append(f"{count} × {token}" if count > 1 else token)
            continue
        action, elements, order = parts
        bond = elements.replace("-", "–")
        multiplicity = f"{count} × " if count > 1 else ""
        rendered.append(f"{action} {multiplicity}{bond} {order.lower()}")
    return "; ".join(rendered)


def _group_inventory_summary(tokens: Iterable[str]) -> str:
    """Render fallback functional-group tokens as readable names."""
    names = sorted(
        {
            str(token).removeprefix("group:").replace("_", " ")
            for token in tokens
            if str(token).startswith("group:")
        }
    )
    return _joined(names)


_INTERPRETATION_DISPLAY_NAMES = {
    "carbonyl_amine_reductive_coupling": "reductive carbonyl–amine coupling",
    "carbonyl_reduction": "carbonyl reduction",
    "sp2_c_activated_c_substitution": "activated sp² C–C substitution",
    "sp2_c_aromatic_ch_substitution": "aromatic C–H substitution",
    "sp2_c_n_substitution": "sp² C–N substitution",
}


def _interpretation_display_name(annotation_id: str) -> str:
    """Return a compact chemist-facing interpretation description."""
    return _INTERPRETATION_DISPLAY_NAMES.get(
        annotation_id,
        annotation_id.replace("_", " "),
    )


def _ambiguity_count(warnings: Iterable[str]) -> int | None:
    """Read the public ambiguity count from a structured warning code."""
    prefix = "AMBIGUOUS_SCAFFOLD_CORRESPONDENCE:"
    for warning in warnings:
        value = str(warning)
        if not value.startswith(prefix):
            continue
        try:
            return int(value.removeprefix(prefix))
        except ValueError:
            return None
    return None


def _atom_reference_label(atom: Any) -> str:
    """Render one provenance-bearing atom compactly for review."""
    if atom is None:
        return "H"
    side = "r" if str(getattr(atom, "side", "")) == "reactant" else "p"
    return (
        f"{getattr(atom, 'element', '?')}"
        f"({side}{int(getattr(atom, 'component_index', 0))}:"
        f"a{int(getattr(atom, 'atom_index', 0))})"
    )


def _edit_hypothesis_lines(result: Any) -> list[str]:
    """Render retained alternatives without promoting them to observations."""
    lines = []
    for index, hypothesis in enumerate(
        getattr(result, "edit_hypotheses", ()) or (),
        start=1,
    ):
        lines.append(
            f"Edit hypothesis {index}: {hypothesis.hypothesis_id} "
            f"({hypothesis.provider}; "
            f"{hypothesis.correspondence_count} correspondences; unverified)"
        )
        hydrogen_change_count = 0
        for edit in hypothesis.edits:
            if edit.edit_type == "hydrogen_change":
                hydrogen_change_count += 1
                continue
            left = _atom_reference_label(edit.atom_1)
            right = _atom_reference_label(edit.atom_2)
            lines.append(
                f"  {edit.edit_type}: {left}–{right}, "
                f"{edit.old_order or 'none'} → {edit.new_order or 'none'}"
            )
        if hydrogen_change_count:
            lines.append(
                f"  hydrogen-count changes: {hydrogen_change_count}"
            )
    return lines


def _reaction_diagnostic_lines(result: Any) -> list[str]:
    """Explain unresolved reaction evidence using already-public observations."""
    lines: list[str] = []
    completeness = getattr(result, "reaction_completeness", None)
    warnings = tuple(getattr(result, "warnings", ()) or ())

    if completeness is not None:
        lines.append(
            "Atom accounting: "
            f"{completeness.reactant_heavy_atom_count} reactant → "
            f"{completeness.product_heavy_atom_count} product heavy atoms"
        )
        lines.append(
            "Unaccounted product atoms: "
            f"{_counted_elements(completeness.product_element_excess)}"
        )
        lines.append(
            "Atoms not in the main product: "
            f"{_counted_elements(completeness.reactant_element_excess)}"
        )
        coverage = completeness.product_heavy_atom_coverage
        coverage_text = (
            f"{100.0 * float(coverage):.1f}%"
            if coverage is not None
            else "not verifiable without atom provenance"
        )
        lines.append(
            f"Completeness: {completeness.status} "
            f"({completeness.evidence}; coverage {coverage_text})"
        )

    ambiguity_count = _ambiguity_count(warnings)
    if ambiguity_count is not None:
        lines.append(
            "Correspondence ambiguity: "
            f"{ambiguity_count} distinct edit hypotheses; "
            "verified edits withheld"
        )
        lines.extend(_edit_hypothesis_lines(result))

    fallback = getattr(result, "fallback_descriptor", None)
    if fallback is not None:
        lines.append(
            "Detected functional groups: reactants "
            f"[{_group_inventory_summary(fallback.reactant_group_tokens)}] → "
            "product "
            f"[{_group_inventory_summary(fallback.product_group_tokens)}]"
        )
    if fallback is not None and fallback.bond_inventory_delta_tokens:
        bond_summary = _bond_inventory_summary(
            fallback.bond_inventory_delta_tokens
        )
        if bond_summary:
            lines.append(
                "Net bond inventory (unmapped, not verified edits): "
                f"{bond_summary}"
            )

    candidates = tuple(getattr(result, "candidates", ()) or ())
    selected = getattr(result, "selected_candidate", None)
    selected_events = tuple(getattr(result, "selected_events", ()) or ())
    if candidates and selected is None and not selected_events:
        verification_counts = Counter(
            str(candidate.verification) for candidate in candidates
        )
        outcomes = ", ".join(
            f"{count} {verification.replace('_', ' ')}"
            for verification, count in sorted(verification_counts.items())
        )
        lines.append(
            f"Interpretation checks: {len(candidates)} candidates "
            f"({outcomes})"
        )
        interpretation_names = sorted(
            {
                _interpretation_display_name(str(candidate.annotation_id))
                for candidate in candidates
            }
        )
        heading = (
            "Rejected interpretations"
            if all(
                candidate.verification
                in {"construction_failed", "product_mismatch"}
                for candidate in candidates
            )
            else "Interpretations checked"
        )
        lines.append(f"{heading}: {_joined(interpretation_names)}")

    signature = getattr(result, "reaction_signature", None)
    if signature is not None:
        lines.append(f"Retrieval: verified signature {signature.signature_id}")
    elif fallback is not None:
        mode = str(fallback.evidence_mode).replace("_", " ")
        if fallback.retrieval_eligible:
            lines.append(
                f"Retrieval: eligible via {mode} "
                f"(confidence {fallback.confidence:.2f})"
            )
        else:
            reasons = ", ".join(
                str(reason).replace("_", " ")
                for reason in fallback.ineligibility_reasons
            ) or "insufficient structural evidence"
            lines.append(
                f"Retrieval: not eligible — {reasons} "
                f"({mode}, confidence {fallback.confidence:.2f})"
            )
    else:
        lines.append("Retrieval: unavailable")

    if "PRODUCT_CONTRADICTED_INTERPRETATION_CANDIDATES" in warnings:
        lines.append(
            "Interpretation note: registered annotations did not match the "
            "supplied product; this does not invalidate the product structure"
        )
    if (
        "REACTION_COMPLETENESS_UNRESOLVED" in warnings
        and completeness is not None
        and not completeness.product_element_excess
    ):
        lines.append(
            "Completeness note: the product has no unexplained elements, but "
            "individual atom origins are not uniquely verified"
        )
    return lines


def _molecule_summary(result: Any) -> str:
    lines = [
        f"valid: {result.valid}",
        f"input: {result.input_smiles}",
        f"canonical: {result.canonical_smiles or '-'}",
        f"components: {len(result.components)}",
        f"reactive sites: {len(result.sites)}",
    ]
    for site in result.sites:
        environment = (site.context_features or {}).get("environment") or {}
        profile = environment.get("reactivity_profile")
        lines.append(
            f"  {site.site_id}: {site.chemist_label} [{site.site_type}; "
            f"availability={site.availability}]"
        )
        lines.append(f"    profile: {render_reactivity_profile(profile)}")
    lines.append(f"functional groups: {len(result.functional_groups)}")
    for group in result.functional_groups:
        lines.append(
            f"  component {group.component_index}: {group.chemist_label} "
            f"[{group.group_id}; atoms={_joined(group.atom_indices)}]"
        )
    if result.warnings:
        lines.append(f"warnings: {_joined(result.warnings)}")
    if result.error:
        lines.append(f"error: {result.error}")
    return "\n".join(lines)


def _reaction_summary(result: Any) -> str:
    reaction_label = result.reaction_label
    lines = [
        f"valid: {result.valid}",
        f"input: {result.input_reaction_smiles}",
        f"evidence: {result.evidence_quality}",
        f"transformation: {result.transformation_class or '-'}",
        f"named family: {result.named_family or '-'}",
        f"compatible families: {_joined(result.compatible_named_families)}",
        f"reaction label: {reaction_label.concise if reaction_label else '-'}",
        f"label status: {reaction_label.status if reaction_label else 'unavailable'}",
        f"candidates: {len(result.candidates)}",
        f"mapped bond changes: {len(result.mapped_bond_changes)}",
    ]
    if result.selected_candidate:
        lines.append(
            "selected interpretation: "
            f"{result.selected_candidate.annotation_id}"
        )
        for change in result.selected_candidate.predicted_bond_changes:
            lines.append(
                f"  edit: {change.change_type} {change.atom_1_role}-{change.atom_2_role} "
                f"{change.old_order or '-'}->{change.new_order or '-'}"
            )
    if result.product_connection:
        lines.append(
            f"product connection: {result.product_connection.concise_label} "
            f"[{result.product_connection.connection_type}; {result.product_connection.evidence}]"
        )
    lines.append(f"spectator groups: {len(result.spectator_groups)}")
    for group in result.spectator_groups:
        distance = "-" if group.graph_distance is None else group.graph_distance
        lines.append(
            f"  component {group.component_index}: {group.chemist_label} "
            f"[{group.group_id}; distance={distance}]"
        )
    if result.family_environment:
        lines.append(f"family environment: {result.family_environment.family_id}")
        for partner in result.family_environment.partners:
            lines.append(
                f"  {partner.role}: {partner.chemist_label} "
                f"[handle={partner.handle_token or '-'}; context={partner.anchor_context or '-'}; "
                f"flags={_joined(partner.flags)}]"
            )
    partner_analysis = _partner_analysis(result)
    if partner_analysis:
        lines.append(f"partner analysis: {partner_analysis}")
    lines.extend(_reaction_core_lines(result))
    if result.warnings:
        lines.append(f"warnings: {_joined(result.warnings)}")
    if result.error:
        lines.append(f"error: {result.error}")
    return "\n".join(lines)


def _molecule_concise_summary(result: Any) -> str:
    lines = [
        f"Molecule: {result.canonical_smiles or result.input_smiles}",
        f"Status: {'valid' if result.valid else 'invalid'}",
    ]
    if result.sites:
        lines.append("Reactive sites:")
        for site in result.sites:
            lines.append(
                f"  {site.chemist_label} — {site.site_type}, {site.availability}"
            )
            environment = (site.context_features or {}).get("environment") or {}
            lines.append(
                "    "
                + render_reactivity_profile(
                    environment.get("reactivity_profile")
                )
            )
    else:
        lines.append("Reactive sites: none")
    if result.functional_groups:
        groups = sorted({group.chemist_label for group in result.functional_groups})
        lines.append(f"Functional groups: {_joined(groups)}")
    else:
        lines.append("Functional groups: none")
    if result.warnings:
        lines.append(f"Warnings: {_joined(result.warnings)}")
    if result.error:
        lines.append(f"Error: {result.error}")
    return "\n".join(lines)


def _reaction_core_lines(result: Any) -> list[str]:
    """Render the minimized mapped-edit observation without overstating it."""
    core = getattr(result, "reaction_core", None)
    if core is None:
        return [
            "Reaction minimization: unavailable "
            "(mapped edit correspondence required)"
        ]
    primary_center_count = sum(
        transition.role == "primary_center"
        for transition in core.atom_transitions
    )
    lines = [
        "Reaction minimization:",
        f"  Minimized reaction: {core.generic_label}",
        (
            f"  Core evidence: {core.evidence_status} "
            f"({core.evidence}; confidence {core.confidence:.3f})"
        ),
        f"  Core shape (retrieval): {core.shape_core_key}",
        (
            "  Center transition (diagnostic only): "
            f"{core.center_transition_key}"
        ),
        (
            f"  Core size: {core.event_count} event(s); "
            f"{primary_center_count} primary center(s); "
            f"{core.active_atom_count} active atom(s)"
        ),
    ]
    if core.abstraction is not None:
        lines[1:2] = [
            f"  General motif: {core.abstraction.general_label}",
            f"  Specific limiter: {core.abstraction.limiter_label}",
            f"  Atom-level reaction: {core.generic_label}",
            f"  Motif key (future retrieval tier): {core.abstraction.motif_key}",
        ]
    if core.remote_subgraphs:
        remote = []
        for subgraph in core.remote_subgraphs:
            fragment = (
                f" [{subgraph.fragment_smiles}]"
                if subgraph.fragment_smiles
                else ""
            )
            port_label = (
                "port"
                if len(subgraph.attachment_ports) == 1
                else "ports"
            )
            remote.append(
                f"{subgraph.side} {subgraph.continuity} "
                f"{subgraph.remote_class}{fragment} "
                f"({len(subgraph.attachment_ports)} {port_label})"
            )
        lines.append(f"  Remote subgraphs: {'; '.join(remote)}")
    else:
        lines.append("  Remote subgraphs: none")
    if core.warnings:
        lines.append(f"  Core warnings: {_joined(core.warnings)}")
    return lines


def _reaction_concise_summary(result: Any) -> str:
    reaction_label = result.reaction_label
    lines = [
        f"Reaction: {reaction_label.concise if reaction_label else '-'}",
        f"Status: {'valid' if result.valid else 'invalid'}",
        f"Evidence: {result.evidence_quality}",
        f"Transformation: {result.transformation_class or 'unknown'}",
        f"Named family: {result.named_family or 'unknown'}",
    ]
    if result.product_connection:
        lines.append(
            f"Product connection: {result.product_connection.concise_label} "
            f"({result.product_connection.connection_type})"
        )
    else:
        lines.append("Product connection: not verified")
    partners = _reaction_partners(result)
    if partners:
        partner_labels = [
            f"{getattr(partner, 'role', None) or 'unassigned'}={partner.chemist_label}"
            for partner in partners
        ]
        lines.append(f"Reactive partners: {_joined(partner_labels)}")
        lines.append(f"Partner analysis: {_partner_analysis(result)}")
    if result.spectator_groups:
        spectators = sorted({group.chemist_label for group in result.spectator_groups})
        lines.append(f"Spectator groups: {_joined(spectators)}")
    lines.extend(("", *_reaction_core_lines(result)))
    diagnostic_lines = _reaction_diagnostic_lines(result)
    if diagnostic_lines:
        lines.extend(("", "Structural evidence:"))
        lines.extend(f"  {line}" for line in diagnostic_lines)
    if result.warnings:
        lines.append(f"Warnings: {_joined(result.warnings)}")
    if result.error:
        lines.append(f"Error: {result.error}")
    return "\n".join(lines)


def format_concise_analysis(result: Any) -> str:
    """Render a molecule or reaction analysis in the concise CLI format."""
    if hasattr(result, "input_reaction_smiles"):
        return _reaction_concise_summary(result)
    return _molecule_concise_summary(result)


def _print_result(result: Any, output_format: str, *, concise: bool = False) -> None:
    if output_format == "json":
        print(_json_dump(result.to_dict()))
    elif concise:
        print(format_concise_analysis(result))
    elif hasattr(result, "input_reaction_smiles"):
        print(_reaction_summary(result))
    else:
        print(_molecule_summary(result))


def _command_validate(args: argparse.Namespace) -> int:
    errors = validate_taxonomy()
    payload = {"valid": not errors, "error_count": len(errors), "errors": errors}
    if args.format == "json":
        print(_json_dump(payload))
    elif errors:
        print(f"taxonomy invalid: {len(errors)} error(s)")
        for error in errors:
            print(f"  - {error}")
    else:
        print("taxonomy valid: all definitions and SMARTS passed validation")
    return 0 if not errors else 1


def _command_molecule(args: argparse.Namespace) -> int:
    result = featurize_molecule(
        args.smiles,
        site_types=args.site_type or None,
        include_context_features=not args.no_context,
        label_style=args.label_style,
    )
    _print_result(result, args.format, concise=args.concise)
    return 0 if result.valid else 1


def _command_reaction(args: argparse.Namespace) -> int:
    result = featurize_reaction(
        args.reaction_smiles,
        label_style=args.label_style,
        max_candidates=args.max_candidates,
    )
    _print_result(result, args.format, concise=args.concise)
    return 0 if result.valid else 1


def _detect_column(fieldnames: Sequence[str], mode: str) -> str | None:
    candidates = (
        ("reaction_smiles", "rxn_smiles", "rxn_smiles_clean", "reaction")
        if mode == "reaction"
        else ("smiles", "canonical_smiles", "molecule_smiles")
    )
    normalized = {name.strip().lower(): name for name in fieldnames}
    return next((normalized[name] for name in candidates if name in normalized), None)


def _batch_summary(records: Sequence[dict[str, Any]], mode: str) -> dict[str, Any]:
    valid_count = sum(bool(record["analysis"].get("valid")) for record in records)
    summary: dict[str, Any] = {
        "mode": mode,
        "total": len(records),
        "valid": valid_count,
        "invalid": len(records) - valid_count,
    }
    if mode == "molecule":
        summary["site_type_counts"] = dict(sorted(Counter(
            site["site_type"]
            for record in records
            for site in record["analysis"].get("sites", [])
        ).items()))
        summary["functional_group_counts"] = dict(sorted(Counter(
            group["group_id"]
            for record in records
            for group in record["analysis"].get("functional_groups", [])
        ).items()))
    else:
        summary["evidence_counts"] = dict(sorted(Counter(
            record["analysis"].get("evidence_quality", "unavailable")
            for record in records
        ).items()))
        summary["family_counts"] = dict(sorted(Counter(
            record["analysis"].get("named_family") or "unknown"
            for record in records
        ).items()))
        summary["label_status_counts"] = dict(sorted(Counter(
            (record["analysis"].get("reaction_label") or {}).get(
                "status", "unavailable"
            )
            for record in records
        ).items()))
    return summary


def _molecule_csv_columns(
    source_fields: Sequence[str],
    smiles_column: str,
) -> list[str]:
    columns = ["source_row"]
    for field in source_fields:
        if field == "source_row":
            continue
        columns.append(field)
        if field == smiles_column:
            columns.extend(("reactive_site_labels", "functional_group_labels"))
    columns.extend(["valid", "canonical_smiles", "component_count", "reactive_site_count"])
    for site_type in _MOLECULE_SITE_TYPES:
        columns.extend((f"{site_type}_count", f"{site_type}_labels"))
    columns.extend((
        "canonical_signatures", "functional_group_count", "functional_group_ids",
        "warnings", "error",
    ))
    return columns


def _molecule_csv_row(record: dict[str, Any]) -> dict[str, Any]:
    analysis = record["analysis"]
    sites = list(analysis.get("sites") or [])
    groups = list(analysis.get("functional_groups") or [])
    row: dict[str, Any] = {
        "source_row": record["source_row"],
        **record["source"],
        "valid": bool(analysis.get("valid")),
        "canonical_smiles": analysis.get("canonical_smiles") or "",
        "component_count": len(analysis.get("components") or []),
        "reactive_site_count": len(sites),
        "reactive_site_labels": "; ".join(str(site.get("chemist_label") or "") for site in sites),
        "canonical_signatures": "; ".join(str(site.get("canonical_signature") or "") for site in sites),
        "functional_group_count": len(groups),
        "functional_group_ids": "; ".join(str(group.get("group_id") or "") for group in groups),
        "functional_group_labels": "; ".join(str(group.get("chemist_label") or "") for group in groups),
        "warnings": "; ".join(str(value) for value in analysis.get("warnings") or []),
        "error": analysis.get("error") or "",
    }
    for site_type in _MOLECULE_SITE_TYPES:
        selected = [site for site in sites if site.get("site_type") == site_type]
        row[f"{site_type}_count"] = len(selected)
        row[f"{site_type}_labels"] = "; ".join(str(site.get("chemist_label") or "") for site in selected)
    return row


def _reaction_csv_columns(
    source_fields: Sequence[str],
    reaction_smiles_column: str,
) -> list[str]:
    columns = ["source_row"]
    for field in source_fields:
        if field == "source_row":
            continue
        columns.append(field)
        if field == reaction_smiles_column:
            columns.extend(
                (
                    "reaction_core_label",
                    "reaction_display_label",
                    "spectator_groups",
                )
            )
    columns.extend((
        "valid", "evidence_quality", "transformation_class", "named_family",
        "partner_analysis", *_REACTION_RING_CSV_FIELDS,
        *_REACTION_CORE_CSV_FIELDS,
        "candidate_count", "warnings", "error",
    ))
    return columns


def _reaction_csv_row(record: dict[str, Any]) -> dict[str, Any]:
    analysis = record["analysis"]
    reaction_core = analysis.get("reaction_core") or {}
    core_abstraction = reaction_core.get("abstraction") or {}
    spectator_groups = analysis.get("spectator_groups") or []
    interpretation = analysis.get("interpretation") or {}
    family_environment = (
        interpretation.get("family_environment")
        or analysis.get("family_environment")
        or {}
    )
    signature_partners = (
        interpretation.get("partners")
        or family_environment.get("partners")
        or (analysis.get("reaction_signature") or {}).get("partners")
        or []
    )

    def dict_role(partner: dict[str, Any]) -> str:
        return str(partner.get("role") or "unassigned")

    role_order = {"electrophile": 0, "nucleophile": 1, "transfer_partner": 2}
    signature_partners = sorted(
        signature_partners,
        key=lambda partner: (
            role_order.get(partner.get("role"), 99),
            dict_role(partner),
            int(partner.get("component_index", -1)),
        ),
    )
    atom_transitions = reaction_core.get("atom_transitions") or []
    core_quality = reaction_core.get("quality") or {}
    core_presentation = reaction_core.get("presentation") or {}
    topology = analysis.get("reaction_topology") or {}
    ring_changes = topology.get("ring_changes") or []
    reaction_label = analysis.get("reaction_label") or {}
    remote_subgraphs = reaction_core.get("remote_subgraphs") or []
    remote_classes = sorted(
        {
            str(subgraph.get("remote_class") or "")
            for subgraph in remote_subgraphs
            if subgraph.get("remote_class")
        }
    )
    remote_summary = sorted(
        (
            f"{subgraph.get('side') or 'unknown'}:"
            f"{subgraph.get('continuity') or 'unknown'}:"
            f"{subgraph.get('remote_class') or 'generic_R'}:"
            f"{subgraph.get('fragment_smiles') or '-'}:"
            f"{len(subgraph.get('attachment_ports') or [])}"
        )
        for subgraph in remote_subgraphs
    )

    return {
        "source_row": record["source_row"],
        **record["source"],
        "valid": bool(analysis.get("valid")),
        "evidence_quality": analysis.get("evidence_quality") or "",
        "transformation_class": analysis.get("transformation_class") or "",
        "named_family": analysis.get("named_family") or "",
        "reaction_core_label": reaction_core.get("generic_label") or "",
        "reaction_display_label": reaction_label.get("concise") or "",
        "spectator_groups": "; ".join(
            str(group.get("group_id") or "") for group in spectator_groups
        ),
        "partner_analysis": "; ".join(
            f"{dict_role(partner)}="
            f"{render_reactivity_profile(partner.get('reactivity_profile'))}"
            for partner in signature_partners
        ),
        "reaction_display_label_detailed": reaction_label.get("detailed") or "",
        "reaction_display_source": reaction_label.get("source") or "",
        "reaction_display_status": reaction_label.get("status") or "",
        "reaction_display_confidence": (
            reaction_label.get("confidence")
            if reaction_label.get("confidence") is not None
            else ""
        ),
        "formed_ring_sizes": "; ".join(
            str(value) for value in topology.get("formed_ring_sizes") or []
        ),
        "ring_count_delta": (
            topology.get("ring_count_delta")
            if topology.get("ring_count_delta") is not None
            else ""
        ),
        "ring_change_count": len(ring_changes),
        "ring_changes_json": _json_dump(ring_changes, compact=True)
        if ring_changes
        else "",
        "reaction_core_available": bool(reaction_core),
        "reaction_core_id": reaction_core.get("core_id") or "",
        "reaction_core_evidence_status": (
            reaction_core.get("evidence_status") or ""
        ),
        "reaction_core_exact_key": reaction_core.get("exact_core_key") or "",
        "reaction_core_typed_key": reaction_core.get("typed_core_key") or "",
        "reaction_core_shape_key": reaction_core.get("shape_core_key") or "",
        "reaction_core_center_transition_key": (
            reaction_core.get("center_transition_key") or ""
        ),
        "reaction_core_mapping_equivalence_key": (
            reaction_core.get("mapping_equivalence_key") or ""
        ),
        "reaction_core_quality_status": core_quality.get("status") or "",
        "reaction_core_quality_reasons": "; ".join(
            tuple(core_quality.get("review_reasons") or ())
            + tuple(core_quality.get("blocking_reasons") or ())
        ),
        "reaction_core_bond_changes": "; ".join(
            core_presentation.get("bond_changes") or ()
        ),
        "reaction_core_state_changes": "; ".join(
            core_presentation.get("atom_state_changes") or ()
        ),
        "reaction_core_retained_context": "; ".join(
            core_presentation.get("retained_context") or ()
        ),
        "reaction_core_departing_context": "; ".join(
            core_presentation.get("departing_context") or ()
        ),
        "reaction_core_appearing_context": "; ".join(
            core_presentation.get("appearing_context") or ()
        ),
        "reaction_core_motif_key": core_abstraction.get("motif_key") or "",
        "reaction_core_general_equation": (
            core_abstraction.get("general_label") or ""
        ),
        "reaction_core_limiter": (
            core_abstraction.get("limiter_label") or ""
        ),
        "reaction_core_event_count": reaction_core.get("event_count") or "",
        "reaction_core_primary_center_count": (
            sum(
                transition.get("role") == "primary_center"
                for transition in atom_transitions
            )
            if reaction_core
            else ""
        ),
        "reaction_core_remote_classes": "; ".join(remote_classes),
        "reaction_core_remote_subgraphs": "; ".join(remote_summary),
        "reaction_core_warnings": "; ".join(
            str(value) for value in reaction_core.get("warnings") or []
        ),
        "candidate_count": len(analysis.get("candidates") or []),
        "warnings": "; ".join(str(value) for value in analysis.get("warnings") or []),
        "error": analysis.get("error") or "",
    }


def _write_batch_csv(
    records: Sequence[dict[str, Any]],
    output_path: Path,
    mode: str,
    source_fields: Sequence[str],
    input_column: str,
    *,
    concise: bool = False,
) -> None:
    if concise:
        columns = [
            "reaction_smiles",
            "reaction_core_label",
            "reaction_display_label",
            "partner_analysis",
            "spectator_groups",
            *_REACTION_RING_CSV_FIELDS,
            *_REACTION_CORE_CSV_FIELDS,
        ]
        row_builder = _reaction_csv_row
    else:
        columns = (
            _molecule_csv_columns(source_fields, input_column)
            if mode == "molecule"
            else _reaction_csv_columns(source_fields, input_column)
        )
        row_builder = _molecule_csv_row if mode == "molecule" else _reaction_csv_row
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = row_builder(record)
            if concise:
                row["reaction_smiles"] = record["source"].get(input_column) or ""
                row = {column: row[column] for column in columns}
            writer.writerow(row)


def _command_batch(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = tuple(reader.fieldnames or ())
        column = args.column or _detect_column(fieldnames, args.mode)
        if not column or column not in fieldnames:
            print(
                f"error: could not find {args.mode} column; available columns: "
                f"{_joined(fieldnames)}",
                file=sys.stderr,
            )
            return 2
        rows = list(reader)

    records: list[dict[str, Any]] = []
    for source_row, row in enumerate(rows, start=2):
        value = str(row.get(column) or "").strip()
        analysis = (
            featurize_reaction(value, label_style=args.label_style, max_candidates=args.max_candidates)
            if args.mode == "reaction"
            else featurize_molecule(value, label_style=args.label_style)
        )
        records.append({
            "source_row": source_row,
            "source": row,
            "analysis": analysis.to_dict(),
        })

    summary = _batch_summary(records, args.mode)
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_format = args.output_format or ("csv" if output_path.suffix.lower() == ".csv" else "jsonl")
        if args.concise and (args.mode != "reaction" or output_format != "csv"):
            print("error: --concise requires reaction mode and CSV output", file=sys.stderr)
            return 2
        if output_format == "csv":
            _write_batch_csv(
                records,
                output_path,
                args.mode,
                fieldnames,
                column,
                concise=args.concise,
            )
        else:
            with output_path.open("w", encoding="utf-8", newline="") as handle:
                for record in records:
                    handle.write(_json_dump(record, compact=True) + "\n")
        summary["output"] = str(output_path)
        summary["output_format"] = output_format
    if args.format == "json":
        print(_json_dump(summary))
    else:
        print(f"mode: {summary['mode']}")
        print(f"records: {summary['total']}")
        print(f"valid: {summary['valid']}")
        print(f"invalid: {summary['invalid']}")
        for key, value in summary.items():
            if key not in {"mode", "total", "valid", "invalid", "output", "output_format"}:
                print(f"{key}: {_json_dump(value, compact=True)}")
        if "output" in summary:
            print(f"{str(summary['output_format']).upper()} output: {summary['output']}")
    return 0 if summary["invalid"] == 0 else 1


def _command_self_test(args: argparse.Namespace) -> int:
    taxonomy_errors = validate_taxonomy()
    molecule_results = [featurize_molecule(value, label_style=args.label_style) for value in _SELF_TEST_MOLECULES]
    reaction_results = [featurize_reaction(value, label_style=args.label_style) for value in _SELF_TEST_REACTIONS]
    payload = {
        "taxonomy_valid": not taxonomy_errors,
        "taxonomy_errors": taxonomy_errors,
        "molecules": [result.to_dict() for result in molecule_results],
        "reactions": [result.to_dict() for result in reaction_results],
    }
    passed = not taxonomy_errors and all(result.valid for result in molecule_results + reaction_results)
    if args.format == "json":
        print(_json_dump(payload))
    else:
        print(f"taxonomy: {'PASS' if not taxonomy_errors else 'FAIL'}")
        print("\nmolecule feature checks")
        for result in molecule_results:
            print(f"  {'PASS' if result.valid else 'FAIL'} {result.input_smiles}: "
                  f"{len(result.sites)} sites, {len(result.functional_groups)} groups")
        print("\nreaction feature checks")
        for result in reaction_results:
            print(f"  {'PASS' if result.valid else 'FAIL'} "
                  f"{result.reaction_label.concise if result.reaction_label else result.input_reaction_smiles}: "
                  f"{result.evidence_quality}, family={result.named_family or 'unknown'}")
        print(f"\noverall: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m reactive_taxonomy.cli",
        description="Validate and exercise standalone compound/reaction featurization.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser("validate", help="validate all taxonomy definitions")
    validate_parser.add_argument("--format", choices=("text", "json"), default="text")
    validate_parser.set_defaults(func=_command_validate)

    molecule_parser = subparsers.add_parser("molecule", help="featurize one molecule SMILES")
    molecule_parser.add_argument("smiles")
    molecule_parser.add_argument("--site-type", action="append", choices=_MOLECULE_SITE_TYPES)
    molecule_parser.add_argument("--no-context", action="store_true")
    molecule_parser.add_argument("--concise", action="store_true", help="show only key chemist-readable features")
    molecule_parser.add_argument("--label-style", choices=("unicode", "ascii", "hte_legacy"), default="unicode")
    molecule_parser.add_argument("--format", choices=("text", "json"), default="text")
    molecule_parser.set_defaults(func=_command_molecule)

    reaction_parser = subparsers.add_parser("reaction", help="featurize one reaction SMILES")
    reaction_parser.add_argument("reaction_smiles")
    reaction_parser.add_argument("--max-candidates", type=int, default=500)
    reaction_parser.add_argument("--concise", action="store_true", help="show only key chemist-readable features")
    reaction_parser.add_argument("--label-style", choices=("unicode", "ascii", "hte_legacy"), default="unicode")
    reaction_parser.add_argument("--format", choices=("text", "json"), default="text")
    reaction_parser.set_defaults(func=_command_reaction)

    batch_parser = subparsers.add_parser("batch", help="test every molecule or reaction in a CSV")
    batch_parser.add_argument("input")
    batch_parser.add_argument("--mode", choices=("molecule", "reaction"), required=True)
    batch_parser.add_argument("--column", help="input column; auto-detected when omitted")
    batch_parser.add_argument("--output", help="optional JSONL or CSV result path")
    batch_parser.add_argument("--output-format", choices=("jsonl", "csv"), help="infer from --output suffix when omitted")
    batch_parser.add_argument("--max-candidates", type=int, default=500)
    batch_parser.add_argument(
        "--concise",
        action="store_true",
        help=(
            "for reaction CSV output, write only reaction_smiles, "
            "reaction_display_label, and spectator_groups"
        ),
    )
    batch_parser.add_argument("--label-style", choices=("unicode", "ascii", "hte_legacy"), default="unicode")
    batch_parser.add_argument("--format", choices=("text", "json"), default="text")
    batch_parser.set_defaults(func=_command_batch)

    self_test_parser = subparsers.add_parser("self-test", help="exercise representative built-in features")
    self_test_parser.add_argument("--label-style", choices=("unicode", "ascii", "hte_legacy"), default="unicode")
    self_test_parser.add_argument("--format", choices=("text", "json"), default="text")
    self_test_parser.set_defaults(func=_command_self_test)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the integrated reactive-taxonomy tester."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
