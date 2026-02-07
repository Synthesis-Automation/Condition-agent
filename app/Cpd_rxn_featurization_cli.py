"""
CLI for testing compound and reaction featurizations.
"""

from __future__ import annotations

import json
import os
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


def _extract_crk_section(crk_key: str, label: str) -> str:
    if not crk_key:
        return ""
    prefix = f"{label}: "
    parts = [p.strip() for p in crk_key.split(" | ") if p.strip()]
    for part in parts:
        if part.startswith(prefix):
            return part[len(prefix):].strip()
    return ""


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


def _supports_color() -> bool:
    if not sys.stdout.isatty():
        return False
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("TERM") == "dumb":
        return False
    return True


def _highlight(text: str) -> str:
    if not text:
        return text
    if not _supports_color():
        return text
    return f"\x1b[96m{text}\x1b[0m"


def _highlight_blue(text: str) -> str:
    if not text:
        return text
    if not _supports_color():
        return text
    return f"\x1b[94m{text}\x1b[0m"


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
    """Extract molecule payload (v2 format only)."""
    return payload


def _get_reaction_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Extract reaction payload (v2 format only)."""
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
        display_name = _clean_compound_id(motif.get("id", "unknown"))
        rank = motif.get("rank", 0)
        
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


def _print_extended_section(extended: Dict[str, Any], *, indent: int = 0) -> None:
    """Print extended section with detailed analysis for conditions recommendation."""
    if not extended:
        return
    prefix = " " * indent
    print(f"{prefix}Extended Analysis (for Conditions Recommendation):")
    
    # Molecule extended info
    per_motif = extended.get("per_motif_analysis")
    if per_motif:
        print(f"{prefix}  Per-Motif Analysis ({len(per_motif)}):")
        for analysis in per_motif:
            motif_id = analysis.get("motif_id", "unknown")
            print(f"{prefix}    - {motif_id}")
            
            steric = analysis.get("steric")
            if steric:
                print(f"{prefix}      steric: score={steric.get('score')}, classification={steric.get('classification')}")
            
            electronic = analysis.get("electronic")
            if electronic:
                print(f"{prefix}      electronic: score={electronic.get('score')}, description={electronic.get('description')}")
            
            nearby = analysis.get("nearby_groups")
            if nearby:
                print(f"{prefix}      nearby_groups: {', '.join(nearby)}")
    
    snar = extended.get("snar_feasibility")
    if snar:
        _print_snar_feasibility(snar, indent=indent + 2)
    
    # Reaction extended info (focus on conditions recommendation data)
    # Note: Detection info is not shown here as it duplicates the main Reaction Type field
    
    aggregates = extended.get("aggregates")
    if aggregates:
        print(f"{prefix}  Aggregates:")
        for key in sorted(aggregates):
            value = aggregates[key]
            print(f"{prefix}    - {key}: {_format_value(value)}")
    
    role_classification = extended.get("role_classification")
    if role_classification:
        reactant_roles = role_classification.get("reactants")
        if reactant_roles:
            # Note: Role classification reaction type is NOT displayed to avoid confusion
            # with the main reaction_type field. Only show the role assignments.
            _print_roles_summary(reactant_roles, indent=indent + 2)


def _print_molecule_detail(
    payload: Dict[str, Any],
    *,
    indent: int = 0,
    show_rdkit: bool = False,
    show_extended: bool = True,
) -> None:
    molecule = _get_molecule_payload(payload)
    prefix = " " * indent
    print(f"{prefix}SMILES: {molecule.get('smiles')}")

    if show_rdkit:
        _print_rdkit_props(molecule.get("rdkit_props") or {}, indent=indent)
    
    # Core information
    _print_motifs(molecule.get("motifs") or [], indent=indent)
    
    # Extended information (useful for conditions recommendation)
    if show_extended:
        extended = molecule.get("extended")
        if extended:
            _print_extended_section(extended, indent=indent)


def _print_molecule_summary(payload: Dict[str, Any], *, show_rdkit: bool = False, show_extended: bool = True) -> None:
    print("\nSummary (Molecule)")
    print("-" * 72)
    _print_molecule_detail(payload, indent=0, show_rdkit=show_rdkit, show_extended=show_extended)


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


def _print_reaction_type_summary(
    reaction_type_str: str | None,
    confidence: float | None,
    *,
    indent: int = 0,
    reaction_key: str | None = None,
) -> None:
    """Print reaction type summary from v2 core format."""
    if not reaction_type_str:
        return
    prefix = " " * indent
    print(f"{prefix}Reaction Type:")
    lines = []
    lines.append(f"reaction_type: {reaction_type_str}")
    if confidence is not None:
        lines.append(f"confidence: {confidence:.4f}")
    if lines:
        print(_format_list(lines, indent=indent + 2))


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
    print(f"{prefix}Reactant Roles:")
    
    # Only show useful metadata (skip redundant reaction_type/confidence)
    num_reactants = roles.get("num_reactants")
    has_multi_functional = roles.get("has_multi_functional_substrates")
    if num_reactants is not None:
        print(f"{prefix}  - num_reactants: {num_reactants}")
    if has_multi_functional is not None:
        print(f"{prefix}  - has_multi_functional_substrates: {has_multi_functional}")
    
    reactants = roles.get("reactants") or []
    if reactants:
        # Count actual reactants (positions) vs total items shown
        actual_count = len(reactants)
        total_alt_count = sum(len(r.get("alternative_functional_groups") or []) for r in reactants)
        
        if total_alt_count > 0:
            print(f"{prefix}  Assignments ({actual_count} reactants, {total_alt_count} alternatives shown):")
        else:
            print(f"{prefix}  Assignments:")
            
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
                print(f"{prefix}      alternatives (other groups in same molecule):")
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
    show_extended: bool = True,
) -> None:
    reaction = _get_reaction_payload(payload)
    def _collect_entry_motifs(entry: Dict[str, Any]) -> List[str]:
        seen = set()
        items: List[str] = []
        for motif in entry.get("motifs") or []:
            if not isinstance(motif, dict):
                continue
            cid = motif.get("id") or motif.get("compound_id")
            if not cid:
                continue
            cid = _clean_compound_id(str(cid))
            if cid in seen:
                continue
            seen.add(cid)
            items.append(cid)
        return items

    print(f"Reaction SMILES: {reaction.get('reaction_smiles')}")
    
    reaction_key = reaction.get("reaction_key")
    product_broad_tags = reaction.get("product_broad_tags") or []
    if reaction_key:
        print(f"CRK: {_highlight_blue(_clean_compound_id(reaction_key))}")
    if product_broad_tags:
        print(f"Broad Product Tags: {_format_value(product_broad_tags)}")
    product_motifs_reactive = reaction.get("product_motifs_reactive") or []
    if product_motifs_reactive:
        print(f"Product Motifs (Reactive): {_format_value(product_motifs_reactive)}")
        
        # Explain the reaction key components from aggregates or extended section
        aggregates = reaction.get("aggregates")
        if not aggregates:
            extended = reaction.get("extended", {})
            aggregates = extended.get("aggregates", {})
        
        if aggregates:
            reacted = aggregates.get("reacted_motifs", [])
            formed = aggregates.get("formed_motifs", [])
            formed_center = aggregates.get("formed_motifs_center", [])
            formed_context = aggregates.get("formed_motifs_context", [])
            spectators = aggregates.get("spectator_motifs", [])
            spectator_groups = (
                aggregates.get("spectator_groups_ranked")
                or aggregates.get("spectator_groups_combined")
                or []
            )
            
            print("  Reaction Key Generation:")
            if reacted:
                if isinstance(reacted, list):
                    reacted_str = "|".join(str(m) for m in reacted)
                else:
                    reacted_str = str(reacted)
                print(f"    Reacted motifs: {reacted_str}")
            
            if formed:
                if isinstance(formed, list):
                    formed_str = ", ".join(str(m) for m in formed)
                else:
                    formed_str = str(formed)
                print(f"    Formed motifs (all): {formed_str}")
            if formed_center:
                if isinstance(formed_center, list):
                    formed_center_str = ", ".join(str(m) for m in formed_center)
                else:
                    formed_center_str = str(formed_center)
                print(f"    Formed motifs (reaction-center): {formed_center_str}")
            if formed_context:
                if isinstance(formed_context, list):
                    formed_context_str = ", ".join(str(m) for m in formed_context)
                else:
                    formed_context_str = str(formed_context)
                print(f"    Formed motifs (context): {formed_context_str}")

            bond_formed = _extract_crk_section(reaction_key or "", "bond_formed")
            bond_broken = _extract_crk_section(reaction_key or "", "bond_broken")
            if bond_formed:
                print(f"    Bond formed: {bond_formed}")
            if bond_broken:
                print(f"    Bond broken: {bond_broken}")
            
            if spectator_groups:
                if isinstance(spectator_groups, list):
                    groups_str = ", ".join(str(m) for m in spectator_groups)
                else:
                    groups_str = str(spectator_groups)
                print(f"    Spectator groups: {groups_str}")
            if spectators:
                if isinstance(spectators, list):
                    spectators_str = ", ".join(str(m) for m in spectators)
                else:
                    spectators_str = str(spectators)
                print(f"    Spectator motifs: {spectators_str}")
            
            print("    Format: |Reactants -> Product | bond_formed: ... | bond_broken: ... | spectators: ...")

    reactants = reaction.get("reactants") or []
    products = reaction.get("products") or []
    print(f"Reactant motifs ({len(reactants)}):")
    if reactants:
        for idx, entry in enumerate(reactants, start=1):
            smiles = entry.get("smiles") or ""
            motifs = _collect_entry_motifs(entry)
            motif_str = ", ".join(motifs) if motifs else "none"
            print(f"  - [{idx}] {smiles}: {motif_str}")
    else:
        print("  - none")

    print(f"Product motifs ({len(products)}):")
    if products:
        for idx, entry in enumerate(products, start=1):
            smiles = entry.get("smiles") or ""
            motifs = _collect_entry_motifs(entry)
            motif_str = ", ".join(motifs) if motifs else "none"
            print(f"  - [{idx}] {smiles}: {motif_str}")
    else:
        print("  - none")

    _print_reaction_type_summary(
        reaction.get("reaction_type"),
        reaction.get("confidence"),
        reaction_key=reaction.get("reaction_key"),
    )


def _print_readable(
    payload: Dict[str, Any],
    *,
    show_roles: bool = False,
    show_rdkit: bool = False,
    show_extended: bool = True,
) -> None:
    kind = payload.get("kind")
    if kind == "molecule":
        _print_molecule_summary(payload, show_rdkit=show_rdkit, show_extended=show_extended)
    elif kind == "reaction":
        _print_reaction_summary(payload, show_roles=show_roles, show_rdkit=show_rdkit, show_extended=show_extended)
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
        "--discovery",
        action="store_true",
        dest="discovery",
        help="Enable discovery mode for undocumented motifs.",
    )
    parser.add_argument(
        "--no-discovery",
        action="store_false",
        dest="discovery",
        help="Disable discovery mode for undocumented motifs (default).",
    )
    parser.set_defaults(discovery=False)
    parser.add_argument(
        "--target-groups",
        help="Focus on specific motifs (e.g., 'Br', 'H', 'CN'). Comma-separated.",
    )
    parser.add_argument(
        "--no-reactant-coverage-guard",
        action="store_false",
        dest="reactant_coverage_guard",
        help="Disable coverage guard that labels reactants without taxonomy motifs.",
    )
    parser.set_defaults(reactant_coverage_guard=True)
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
        "reactant_coverage_guard": args.reactant_coverage_guard,
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

        current_options = dict(options)
        is_reaction = ">" in smiles
        if not is_reaction:
            targets = _prompt("Target Groups (optional, e.g. Br,CN): ").strip()
            if targets:
                current_options["target_groups"] = [tg.strip() for tg in targets.split(",")]

        try:
            if mode in {"reaction", "r"} or (mode == "auto" and is_reaction):
                reaction_options = dict(current_options)
                reaction_options.setdefault("motif_site_filter", "substituent")
                reaction_options["detailed"] = True  # Enable extended output
                payload = featurize_reaction(smiles, options=reaction_options)
            else:
                molecule_options = dict(current_options)
                molecule_options["detailed"] = True  # Enable extended output
                payload = featurize_molecule(smiles, options=molecule_options)
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
