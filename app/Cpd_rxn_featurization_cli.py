"""
CLI for testing compound and reaction featurizations.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
import argparse
from functools import lru_cache
from typing import Any, Dict, Iterable, List

# Ensure repo root is on sys.path for direct execution.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_molecule, featurize_reaction


def _print_json(payload: Dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True))


def _format_list(items: Iterable[str], *, indent: int = 2) -> str:
    prefix = " " * indent
    return "\n".join(f"{prefix}- {item}" for item in items)


def _format_value(value: Any) -> str:
    if value is None:
        return "none"
    if isinstance(value, float):
        return f"{value:.3f}".rstrip("0").rstrip(".")
    if isinstance(value, list):
        if not value:
            return "none"
        if all(isinstance(v, (str, int, float, bool)) for v in value):
            # Clean up compound IDs in lists
            return ", ".join(_clean_compound_id(str(v)) for v in value)
        return json.dumps(value, sort_keys=True)
    if isinstance(value, dict):
        if not value:
            return "none"
        return json.dumps(value, sort_keys=True)
    return _clean_compound_id(str(value))


def _format_kv_inline(data: Dict[str, Any]) -> str:
    if not data:
        return "none"
    parts = []
    for key in sorted(data):
        parts.append(f"{key}={_format_value(data[key])}")
    return ", ".join(parts)


def _format_mapping_lines(data: Dict[str, Any], *, sort_keys: bool = True) -> List[str]:
    if not data:
        return []
    keys = sorted(data) if sort_keys else data.keys()
    lines = []
    for key in keys:
        value = data.get(key)
        if value is None or value == "" or value == [] or value == {}:
            continue
        lines.append(f"{key}: {_format_value(value)}")
    return lines


_STERIC_SUMMARY_KEYS = (
    "score_0_10",
    "group_bulk_0_10",
    "beta_hydrogens",
    "description",
    "classification",
)
_ELECTRONIC_SUMMARY_KEYS = ("score_0_10", "including_group_score_0_10", "description")


def _format_summary(data: Dict[str, Any], keys: Iterable[str]) -> str:
    summary: Dict[str, Any] = {}
    for key in keys:
        if key in data:
            value = data.get(key)
            if value is None or value == "" or value == [] or value == {}:
                continue
            summary[key] = value
    return _format_kv_inline(summary) if summary else "none"


def _clean_compound_id(cid: str) -> str:
    """Clean up double dashes in compound IDs for display."""
    if isinstance(cid, str):
        return cid.replace("--", "-")
    return str(cid)


def _format_compound_id(entry: Dict[str, Any]) -> str:
    cid = entry.get("compound_id", "unknown")
    cid = _clean_compound_id(cid)
    alt_ids = entry.get("alt_compound_ids")
    if alt_ids:
        if isinstance(alt_ids, (set, list)):
            # Use the provided order (pre-sorted by complexity/specificity)
            all_ids = [cid] + [_clean_compound_id(str(alt)) for alt in alt_ids]
            cid = ", ".join(all_ids)
    if entry.get("undocumented"):
        cid += " [UNDOCUMENTED]"
    return cid


def _format_primary_compound_id(entry: Dict[str, Any]) -> str:
    cid = entry.get("compound_id", "unknown")
    if entry.get("undocumented"):
        cid += " [UNDOCUMENTED]"
    return cid


@lru_cache(maxsize=1)
def _load_compound_registry() -> Dict[str, Any]:
    from chemtools.featurizers.motifs import registry

    return registry.build_compound_registry(registry._default_registry_paths())


def _compound_specificity_key(compound_id: str) -> tuple[int, int, int, str] | None:
    if not compound_id:
        return None
    try:
        registry = _load_compound_registry()
    except Exception:
        return None
    pattern = registry.get("compound_map", {}).get(compound_id)
    if not pattern:
        return None
    return (pattern.complexity, pattern.priority, len(compound_id), compound_id)


def _filter_most_specific_analyses(
    analyses: Iterable[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    analysis_list = [entry for entry in analyses if isinstance(entry, dict)]
    if len(analysis_list) <= 1:
        return analysis_list

    scored: List[tuple[tuple[int, int, int, str], Dict[str, Any]]] = []
    for entry in analysis_list:
        key = _compound_specificity_key(entry.get("compound_id", ""))
        if key is None:
            continue
        scored.append((key, entry))

    if not scored:
        return analysis_list

    max_key = max(key for key, _ in scored)
    winners = [entry for key, entry in scored if key == max_key]
    return winners or analysis_list


def _get_molecule_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    # Handle new core format (v2) - no nested wrapper
    if payload.get("schema_version") == "v2" and payload.get("kind") == "molecule":
        return payload
    # Handle legacy format (v1) - nested under 'molecule' key
    if "molecule" in payload and isinstance(payload.get("molecule"), dict):
        return payload["molecule"]
    return payload


def _get_reaction_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    # Handle new core format (v2) - no nested wrapper
    if payload.get("schema_version") == "v2" and payload.get("kind") == "reaction":
        return payload
    # Handle legacy format (v1) - nested under 'reaction' key
    if "reaction" in payload and isinstance(payload.get("reaction"), dict):
        return payload["reaction"]
    return payload


def _print_meta_section(title: str, meta: Dict[str, Any], *, indent: int = 0) -> None:
    if not meta:
        return
    prefix = " " * indent
    print(f"{prefix}{title}:")
    lines = _format_mapping_lines(meta)
    if lines:
        print(_format_list(lines, indent=indent + 2))


def _print_rdkit_props(props: Dict[str, Any], *, indent: int = 0) -> None:
    if not props:
        return
    prefix = " " * indent
    print(f"{prefix}RDKit Properties:")
    lines = _format_mapping_lines(props)
    if lines:
        print(_format_list(lines, indent=indent + 2))


def _print_motifs(motifs: Iterable[Dict[str, Any]], *, indent: int = 0) -> None:
    motifs_list = list(motifs or [])
    if not motifs_list:
        return
    prefix = " " * indent
    print(f"{prefix}Motifs ({len(motifs_list)}):")
    lines = []
    for motif in motifs_list:
        # v2 format: {"id": "Ar-Br", "rank": 450, "a_atom_idx": 5}
        # v1 format: {"compound_id": "Ar-Br", "rank_score": 450, "a_atom_idx": 5}
        if "id" in motif:
            # v2 format
            display_name = _clean_compound_id(motif.get("id", "unknown"))
            rank = motif.get("rank", 0)
        else:
            # v1 format
            display_name = _format_compound_id(motif)
            rank = motif.get("rank_score", 0)
        
        a_idx = motif.get("a_atom_idx")
        b_idx = motif.get("b_atom_idx")
        bond = motif.get("bond")
        
        details = []
        if a_idx is not None:
            details.append(f"a={a_idx}")
        if b_idx is not None:
            details.append(f"b={b_idx}")
        if bond is not None:
            details.append(f"bond={bond}")
        if isinstance(rank, (int, float)) and rank > 0:
            details.append(f"rank={float(rank):.2f}")
        
        if details:
            lines.append(f"{display_name} ({', '.join(details)})")
        else:
            lines.append(display_name)
    print(_format_list(lines, indent=indent + 2))


def _print_nearby_groups(groups: Iterable[Dict[str, Any]], *, indent: int = 0) -> None:
    prefix = " " * indent
    group_list = list(groups or [])
    if not group_list:
        print(f"{prefix}nearby_groups: none")
        return
    print(f"{prefix}nearby_groups:")
    for entry in group_list:
        if not isinstance(entry, dict):
            print(f"{prefix}  - {entry}")
            continue
        name = entry.get("name") or "unknown"
        details = [name]
        group_id = entry.get("group_id")
        if group_id:
            details.append(f"group_id={group_id}")
        tags = entry.get("tags") or []
        if tags:
            details.append(f"tags={_format_value(tags)}")
        print(f"{prefix}  - {'; '.join(details)}")


def _print_motif_analyses(analyses: Iterable[Dict[str, Any]], *, indent: int = 0) -> None:
    analysis_list = _filter_most_specific_analyses(analyses or [])
    if not analysis_list:
        return
    prefix = " " * indent
    print(f"{prefix}Per-Motif Analysis (most specific, {len(analysis_list)}):")
    for entry in analysis_list:
        cid = _format_primary_compound_id(entry)
        print(f"{prefix}  - {cid}")
        center = entry.get("center") or {}
        if center:
            print(f"{prefix}    center: {_format_kv_inline(center)}")
        steric = entry.get("steric")
        if steric:
            if isinstance(steric, dict):
                summary = _format_summary(steric, _STERIC_SUMMARY_KEYS)
                print(f"{prefix}    steric: {summary}")
            else:
                print(f"{prefix}    steric: {_format_value(steric)}")
        electronic = entry.get("electronic")
        if electronic is not None:
            if isinstance(electronic, list):
                print(f"{prefix}    electronic:")
                for idx, e_entry in enumerate(electronic, start=1):
                    if isinstance(e_entry, dict):
                        formatted = _format_summary(e_entry, _ELECTRONIC_SUMMARY_KEYS)
                    else:
                        formatted = _format_value(e_entry)
                    print(f"{prefix}      - [{idx}] {formatted}")
            elif isinstance(electronic, dict):
                summary = _format_summary(electronic, _ELECTRONIC_SUMMARY_KEYS)
                print(f"{prefix}    electronic: {summary}")
            else:
                print(f"{prefix}    electronic: {_format_value(electronic)}")
        if "nearby_groups" in entry:
            _print_nearby_groups(entry.get("nearby_groups") or [], indent=indent + 4)


def _print_snar_feasibility(items: Iterable[Dict[str, Any]], *, indent: int = 0) -> None:
    snar_list = list(items or [])
    if not snar_list:
        return
    prefix = " " * indent
    print(f"{prefix}SNAr Feasibility:")
    for item in snar_list:
        status = "YES" if item.get("feasible") else "NO"
        motif = item.get("motif") or "unknown"
        confidence = item.get("confidence")
        score = item.get("score")
        details = [status]
        if confidence is not None:
            details.append(f"confidence={confidence}")
        if score is not None:
            details.append(f"score={score}")
        print(f"{prefix}  - {motif}: {', '.join(details)}")
        reason = item.get("reason")
        if reason:
            print(f"{prefix}    Reason: {reason}")


def _print_molecule_detail(
    payload: Dict[str, Any],
    *,
    indent: int = 0,
    show_rdkit: bool = False,
) -> None:
    molecule = _get_molecule_payload(payload)
    prefix = " " * indent
    print(f"{prefix}SMILES: {molecule.get('smiles')}")

    if show_rdkit:
        _print_rdkit_props(molecule.get("rdkit_props") or {}, indent=indent)
    _print_motifs(molecule.get("motifs") or [], indent=indent)
    ranked = molecule.get("ranked_motifs") or []
    if ranked:
        print(f"{prefix}Ranked Motifs:")
        print(_format_list([_clean_compound_id(str(item)) for item in ranked], indent=indent + 2))
    _print_motif_analyses(molecule.get("analyses") or [], indent=indent)
    _print_snar_feasibility(molecule.get("snar_feasibility") or [], indent=indent)


def _print_molecule_summary(payload: Dict[str, Any], *, show_rdkit: bool = False) -> None:
    print("\nSummary (Molecule)")
    print("-" * 72)
    _print_molecule_detail(payload, indent=0, show_rdkit=show_rdkit)


def _print_normalized_entries(title: str, entries: Iterable[Dict[str, Any]], *, indent: int = 0) -> None:
    entry_list = list(entries or [])
    prefix = " " * indent
    print(f"{prefix}{title} ({len(entry_list)}):")
    if not entry_list:
        print(f"{prefix}  - none")
        return
    for idx, entry in enumerate(entry_list, start=1):
        label = (
            entry.get("smiles_norm")
            or entry.get("largest_smiles")
            or entry.get("input")
            or "unknown"
        )
        print(f"{prefix}  - [{idx}] {label}")
        for key in ("input", "smiles_norm", "largest_smiles", "fragments", "error"):
            value = entry.get(key)
            if value is None or value == "" or value == []:
                continue
            formatted = _format_value(value) if key == "fragments" else value
            print(f"{prefix}    {key}: {formatted}")


def _print_reaction_type(
    reaction_type: Dict[str, Any],
    *,
    indent: int = 0,
    reaction_key: str | None = None,
) -> None:
    if not reaction_type:
        return
    prefix = " " * indent
    print(f"{prefix}Reaction Type:")
    lines = []
    if reaction_key:
        lines.append(f"reaction_key: {reaction_key}")
    for key in ("reaction_type", "name", "category", "confidence"):
        value = reaction_type.get(key)
        if value is not None:
            if key == "confidence":
                lines.append(f"{key}: {value:.4f}")
            else:
                lines.append(f"{key}: {value}")
    if lines:
        print(_format_list(lines, indent=indent + 2))
    
    alternatives = reaction_type.get("alternatives")
    if alternatives:
        print(f"{prefix}  Alternatives:")
        for alt in alternatives:
            alt_name = alt.get("name") or alt.get("reaction_type")
            alt_conf = alt.get("confidence", 0.0)
            print(f"{prefix}    - {alt_name} [Conf: {alt_conf:.4f}]")

    slot_evidence = reaction_type.get("slot_evidence") or {}
    if slot_evidence:
        print(f"{prefix}  Slot Evidence:")
        for slot in sorted(slot_evidence):
            motifs = slot_evidence.get(slot) or []
            print(f"{prefix}    - {slot}: {_format_value(motifs)}")


def _print_detection(detection: Dict[str, Any], *, indent: int = 0) -> None:
    if not detection:
        return
    prefix = " " * indent
    error = detection.get("error")
    if error:
        print(f"{prefix}Detection Error: {error}")
    matches = detection.get("matches") or []
    # Show top 3 matches by confidence
    top_matches = matches[:3]
    print(f"{prefix}Detection Matches (Showing top {len(top_matches)} of {len(matches)}):")
    if not matches:
        print(f"{prefix}  - none")
        return
    for match in top_matches:
        name = match.get("name") or match.get("reaction_type") or "unknown"
        category = match.get("category")
        confidence = match.get("confidence")
        
        header = name
        if category:
            header += f" ({category})"
        if confidence is not None:
            header += f" [Confidence: {confidence:.2f}]"
            
        print(f"{prefix}  - {header}")
        slot_evidence = match.get("slot_evidence") or {}
        for slot in sorted(slot_evidence):
            motifs = slot_evidence.get(slot) or []
            print(f"{prefix}    {slot}: {_format_value(motifs)}")


def _print_roles_summary(roles: Dict[str, Any], *, indent: int = 0) -> None:
    if not roles:
        return
    prefix = " " * indent
    print(f"{prefix}Role Classification:")
    lines = []
    for key in (
        "reaction_type",
        "confidence",
        "detection_method",
        "num_reactants",
        "has_multi_functional_substrates",
    ):
        value = roles.get(key)
        if value is not None:
            lines.append(f"{key}: {value}")
    if lines:
        print(_format_list(lines, indent=indent + 2))
    reactants = roles.get("reactants") or []
    if reactants:
        print(f"{prefix}  Reactants:")
        for reactant in reactants:
            pos = reactant.get("position")
            name = reactant.get("name") or reactant.get("category") or "unknown"
            label = f"pos {pos}: {name}" if pos is not None else name
            print(f"{prefix}    - {label}")
            detail_parts = []
            for key in ("category", "member_type", "role", "is_expected", "has_alternatives"):
                value = reactant.get(key)
                if value is not None:
                    detail_parts.append(f"{key}={value}")
            if detail_parts:
                print(f"{prefix}      {', '.join(detail_parts)}")
            alternatives = reactant.get("alternative_functional_groups") or []
            if alternatives:
                print(f"{prefix}      alternatives:")
                for alt in alternatives:
                    alt_label = alt.get("name") or alt.get("category") or "unknown"
                    alt_parts = []
                    for key in ("category", "member_type"):
                        value = alt.get(key)
                        if value is not None:
                            alt_parts.append(f"{key}={value}")
                    if alt_parts:
                        alt_label = f"{alt_label} ({', '.join(alt_parts)})"
                    print(f"{prefix}        - {alt_label}")


def _print_agent_roles(agent_roles: Dict[str, Any], *, indent: int = 0) -> None:
    if not agent_roles:
        return
    prefix = " " * indent
    print(f"{prefix}Agent Roles:")
    agent_count = agent_roles.get("agent_count")
    if agent_count is not None:
        print(f"{prefix}  - agent_count: {agent_count}")
    for label in ("role_counts", "family_counts", "role_flags", "flags"):
        data = agent_roles.get(label) or {}
        if data:
            print(f"{prefix}  - {label}:")
            lines = _format_mapping_lines(data)
            if lines:
                print(_format_list(lines, indent=indent + 4))
    agents = agent_roles.get("agents") or []
    if agents:
        print(f"{prefix}  Agents:")
        for agent in agents:
            label = agent.get("input") or agent.get("smiles") or "unknown"
            print(f"{prefix}    - {label}")
            details = []
            for key in ("role_name", "family_name", "match_kind"):
                value = agent.get(key)
                if value:
                    details.append(f"{key}={value}")
            if details:
                print(f"{prefix}      {', '.join(details)}")


def _print_reaction_feasibility(feasibility: Dict[str, Any], *, indent: int = 0) -> None:
    if not feasibility:
        return
    prefix = " " * indent
    print(f"{prefix}Feasibility Analysis:")
    status = "YES" if feasibility.get("feasible") else "NO"
    confidence = feasibility.get("confidence")
    score = feasibility.get("score")
    details = [f"possible={status}"]
    if confidence is not None:
        details.append(f"confidence={confidence}")
    if score is not None:
        details.append(f"score={score}")
    print(f"{prefix}  - {', '.join(details)}")
    reason = feasibility.get("reason")
    if reason:
        print(f"{prefix}  - reason: {reason}")
    extra = feasibility.get("details") or {}
    if extra:
        print(f"{prefix}  - details:")
        lines = _format_mapping_lines(extra)
        if lines:
            print(_format_list(lines, indent=indent + 4))


def _print_reaction_summary(
    payload: Dict[str, Any],
    *,
    show_roles: bool = False,
    show_rdkit: bool = False,
) -> None:
    reaction = _get_reaction_payload(payload)
    print("\nSummary (Reaction)")
    print("-" * 72)
    print(f"Reaction SMILES: {reaction.get('reaction_smiles')}")
    reaction_key = reaction.get("reaction_key")
    if reaction_key:
        print(f"Reaction Key: {_clean_compound_id(reaction_key)}")

    # normalized = reaction.get("normalized") or {}
    # if normalized:
    #     print("Normalization:")
    #     normalization_lines = _format_mapping_lines(
    #         {
    #             "input": normalized.get("input"),
    #             "normalized": normalized.get("normalized"),
    #         }
    #     )
    #     if normalization_lines:
    #         print(_format_list(normalization_lines, indent=2))
    #     errors = normalized.get("errors") or []
    #     if errors:
    #         print(_format_list([f"errors: {_format_value(errors)}"], indent=2))
    #     _print_normalized_entries("Reactants", normalized.get("reactants") or [], indent=2)
    #     _print_normalized_entries("Agents", normalized.get("agents") or [], indent=2)
    #     _print_normalized_entries("Products", normalized.get("products") or [], indent=2)

    _print_reaction_type(
        reaction.get("reaction_type") or {},
        reaction_key=reaction.get("reaction_key"),
    )
    _print_detection(reaction.get("detection") or {})

    aggregates = reaction.get("aggregates") or {}
    if aggregates:
        print("Aggregates:")
        lines = _format_mapping_lines(aggregates)
        if lines:
            print(_format_list(lines))

    if show_roles:
        _print_roles_summary(reaction.get("roles") or {})
        _print_agent_roles(reaction.get("agent_roles") or {})

    intramolecular = reaction.get("intramolecular")
    if intramolecular:
        print("Intramolecular:")
        if isinstance(intramolecular, dict):
            lines = _format_mapping_lines(intramolecular)
            if lines:
                print(_format_list(lines, indent=2))
        else:
            print(f"  - {intramolecular}")

    reactants = reaction.get("reactants") or []
    print(f"Reactants ({len(reactants)}):")
    if reactants:
        for idx, reactant in enumerate(reactants, start=1):
            print(f"  Reactant {idx}:")
            _print_molecule_detail(reactant, indent=4, show_rdkit=show_rdkit)
    else:
        print("  - none")

    products = reaction.get("products") or []
    print(f"Products ({len(products)}):")
    if products:
        for idx, product in enumerate(products, start=1):
            print(f"  Product {idx}:")
            _print_molecule_detail(product, indent=4, show_rdkit=show_rdkit)
    else:
        print("  - none")

    _print_reaction_feasibility(reaction.get("feasibility") or {})


def _print_readable(
    payload: Dict[str, Any],
    *,
    show_roles: bool = False,
    show_rdkit: bool = False,
) -> None:
    kind = payload.get("kind")
    if kind == "molecule":
        _print_molecule_summary(payload, show_rdkit=show_rdkit)
    elif kind == "reaction":
        _print_reaction_summary(payload, show_roles=show_roles, show_rdkit=show_rdkit)
    else:
        print("\nSummary")
        print("-" * 72)
        print(f"Kind: {kind or 'unknown'}")


def _prompt(text: str) -> str:
    try:
        return input(text)
    except EOFError:
        return "q"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Interactive CLI for ChemTools featurization.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Formats:\n"
            "  summary     - readable summary only (default)\n"
            "  full-json   - full JSON payload only\n"
            "  both        - summary + full JSON\n"
        ),
    )
    parser.add_argument(
        "--format",
        default="summary",
        choices=["summary", "full-json", "both"],
        help="Output format (default: summary).",
    )
    parser.add_argument(
        "--include-ar-h",
        action="store_true",
        help="Include Ar-H motifs even if other motifs are present.",
    )
    parser.add_argument(
        "--no-discovery",
        action="store_false",
        dest="discovery",
        help="Disable discovery mode for undocumented motifs.",
    )
    parser.set_defaults(discovery=True)
    parser.add_argument(
        "--target-groups",
        help="Focus on specific motifs (e.g., 'Br', 'H', 'CN'). Comma-separated.",
    )
    parser.add_argument(
        "--show-roles",
        action="store_true",
        help="Show agent roles and role classification in reaction summaries.",
    )
    parser.add_argument(
        "--show-rdkit",
        action="store_true",
        help="Show RDKit properties in molecule summaries.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    
    target_groups = None
    if args.target_groups:
        target_groups = [tg.strip() for tg in args.target_groups.split(",")]

    options = {
        "include_ar_h": args.include_ar_h,
        "target_groups": target_groups,
        "discovery_mode": args.discovery,
        "confirm_coupling_products": True,
    }
    print("ChemTools Featurization CLI")
    print("Enter 'q' to quit.")

    while True:
        mode = "auto"  # Default mode
        if mode in {"q", "quit", "exit"}:
            return 0
        if mode not in {"compound", "reaction", "auto", "c", "r"}:
            print("Invalid mode. Choose compound, reaction, or auto.")
            continue

        smiles = _prompt("SMILES: ").strip()
        if not smiles:
            print("Empty input. Try again.")
            continue
        if smiles.lower() in {"q", "quit", "exit"}:
            return 0

        targets = _prompt("Target Groups (optional, e.g. Br,CN): ").strip()
        current_options = dict(options)
        if targets:
            current_options["target_groups"] = [tg.strip() for tg in targets.split(",")]

        try:
            if mode in {"reaction", "r"} or (mode == "auto" and ">" in smiles):
                reaction_options = dict(current_options)
                reaction_options.setdefault("motif_site_filter", "substituent")
                payload = featurize_reaction(smiles, options=reaction_options)
            else:
                payload = featurize_molecule(smiles, options=current_options)
        except Exception as exc:
            print(f"Error: {exc}")
            continue

        if args.format in {"summary", "both"}:
            _print_readable(
                payload,
                show_roles=args.show_roles,
                show_rdkit=args.show_rdkit,
            )
        if args.format in {"full-json", "both"}:
            print("\nFull Payload")
            print("-" * 72)
            _print_json(payload)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
