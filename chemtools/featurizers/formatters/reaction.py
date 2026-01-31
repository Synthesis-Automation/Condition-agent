"""
Reaction-level feature formatting and bundling.

Handles reaction type detection, reactant/product processing, and reaction key generation.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from chemtools.util import rdkit_helpers
from chemtools.smiles import normalize_reaction

from ..analysis.reaction_context import classify_reactants_with_context, get_reactant_summary
from ..analysis.feasibility import analyze_snar_feasibility

from .molecule import build_molecule_bundle, to_bool
from .aggregation import aggregate_reaction_features, infer_intramolecular
from .utils import extract_motif_ids
from .simplified import build_core_reaction, build_extended_reaction


@lru_cache(maxsize=1)
def _load_motif_sets() -> Dict[str, Set[str]]:
    path = Path(__file__).resolve().parent.parent.parent / "taxonomy" / "data" / "compound_logic.json"
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)
    except Exception:
        return {}
    motif_sets: Dict[str, Set[str]] = {}
    for set_name, set_data in (payload.get("motif_sets", {}) or {}).items():
        members = set_data.get("members", []) or []
        motif_sets[set_name] = set(str(m) for m in members if m)
    return motif_sets


def _expand_pattern(pattern: Any, motif_sets: Dict[str, Set[str]]) -> Set[str]:
    if pattern is None:
        return set()
    if isinstance(pattern, str):
        if pattern.startswith("@"):
            return motif_sets.get(pattern[1:], set())
        return {pattern}
    if isinstance(pattern, list):
        expanded: Set[str] = set()
        for item in pattern:
            expanded.update(_expand_pattern(item, motif_sets))
        return expanded
    return set()


@lru_cache(maxsize=1)
def _load_reactant_slot_sets() -> tuple[Set[str], Set[str]]:
    """Load electrophile/nucleophile motif sets from the reaction taxonomy."""
    path = Path(__file__).resolve().parent.parent.parent / "taxonomy" / "data" / "reaction_types.v4.0.json"
    if not path.exists():
        return set(), set()
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)
    except Exception:
        return set(), set()
    motif_sets = _load_motif_sets()
    electrophiles: Set[str] = set()
    nucleophiles: Set[str] = set()
    for reaction_def in payload.get("reaction_types", []) or []:
        reactants = reaction_def.get("reactants", {}) or {}
        for slot_name, slot_pattern in reactants.items():
            if slot_name == "nucleophile":
                nucleophiles.update(_expand_pattern(slot_pattern, motif_sets))
            elif slot_name in {"electrophile", "substrate"}:
                electrophiles.update(_expand_pattern(slot_pattern, motif_sets))
    return electrophiles, nucleophiles


@lru_cache(maxsize=1)
def _load_product_to_nucleophile_map() -> Dict[str, Set[str]]:
    """Map product motifs to allowed nucleophile motifs from the taxonomy."""
    path = Path(__file__).resolve().parent.parent.parent / "taxonomy" / "data" / "reaction_types.v4.0.json"
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)
    except Exception:
        return {}
    motif_sets = _load_motif_sets()
    mapping: Dict[str, Set[str]] = {}
    for reaction_def in payload.get("reaction_types", []) or []:
        reactants = reaction_def.get("reactants", {}) or {}
        products = reaction_def.get("products", {}) or {}
        nucleophile_pattern = reactants.get("nucleophile")
        if not nucleophile_pattern:
            continue
        nucleophiles = _expand_pattern(nucleophile_pattern, motif_sets)
        if not nucleophiles:
            continue
        for slot_pattern in products.values():
            product_set = _expand_pattern(slot_pattern, motif_sets)
            for motif_id in product_set:
                mapping.setdefault(motif_id, set()).update(nucleophiles)
    return mapping


def format_reaction_type_summary(detection: Any) -> Dict[str, Any]:
    """Extract reaction type information with alternatives."""
    matches = detection.matches if detection else []
    if not matches:
        return {"reaction_type": "Unknown", "confidence": 0.0, "slot_evidence": {}}
    
    best = matches[0]
    # Construct slot_evidence from new API structure
    slot_evidence = {}
    if best.electrophile:
        slot_evidence["electrophile"] = best.electrophile
    if best.nucleophile:
        slot_evidence["nucleophile"] = best.nucleophile
    if best.product:
        slot_evidence["product"] = best.product
    
    result = {
        "reaction_type": best.reaction_type,
        "name": best.reaction_type,  # Use reaction_type as name for compatibility
        "category": "coupling",  # Default category - could be improved
        "confidence": best.confidence,
        "slot_evidence": slot_evidence,
    }
    
    # Add alternatives if available (top 3 total)
    if len(matches) > 1:
        alts = []
        for alt in matches[1:3]:
            alts.append({
                "reaction_type": alt.reaction_type,
                "name": alt.reaction_type,
                "confidence": alt.confidence
            })
        result["alternatives"] = alts
        
    return result


def format_reaction_key(
    reacted: Iterable[str],
    formed: Iterable[str],
    spectators: Iterable[str]
) -> str:
    """
    Format motif sets into standardized Reaction_Key string.
    
    Format: Reacted -> Formed || Spectators
    Example: Ar-Br|R_acidic-H -> Ar-Alkyl || Ar-COR|RCH2-COR|R_acidic-H
    
    Args:
        reacted: Motifs consumed in the reaction
        formed: Motifs created in the product
        spectators: Motifs present but unchanged
    
    Returns:
        Formatted Reaction_Key string
    """
    reacted_list = sorted(reacted) if reacted else []
    formed_list = sorted(formed) if formed else []
    spectators_list = sorted(spectators) if spectators else []
    
    reacted_str = "|".join(reacted_list) if reacted_list else "[]"
    formed_str = "|".join(formed_list) if formed_list else "[]"
    spectators_str = "|".join(spectators_list) if spectators_list else "[]"
    
    return f"{reacted_str} -> {formed_str} || {spectators_str}"


def _parse_reaction_key_parts(key: str) -> Tuple[List[str], List[str], List[str]]:
    """Parse a Reaction_Key into reacted/formed/spectator motif lists."""
    if not key or " -> " not in key:
        return [], [], []
    main_part, *rest = key.split(" || ")
    reacted_part, formed_part = (main_part.split(" -> ") + [""])[:2]
    spectators_part = rest[0] if rest else ""

    def parse_part(part: str) -> List[str]:
        part = part.strip()
        if not part or part in {"[]", "None"}:
            return []
        if part.startswith("[") and part.endswith("]"):
            part = part[1:-1]
        return [token.strip() for token in part.split("|") if token.strip() and token.strip() != "[]"]

    return parse_part(reacted_part), parse_part(formed_part), parse_part(spectators_part)


def _atom_signature(atom: Any) -> Tuple[int, bool, int, int, Tuple[int, ...]]:
    """Build a simple signature for mapped-atom change detection."""
    return (
        atom.GetAtomicNum(),
        bool(atom.GetIsAromatic()),
        int(atom.GetTotalDegree()),
        int(atom.GetTotalNumHs()),
        tuple(sorted(n.GetAtomicNum() for n in atom.GetNeighbors())),
    )


def _collect_mapped_signatures(mols: List[Any]) -> Tuple[Dict[int, Tuple[int, bool, int, int, Tuple[int, ...]]], bool]:
    """Collect atom map signatures; return (map -> signature, conflict_found)."""
    signatures: Dict[int, Tuple[int, bool, int, int, Tuple[int, ...]]] = {}
    conflict = False
    for mol in mols:
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if not map_num:
                continue
            sig = _atom_signature(atom)
            if map_num in signatures and signatures[map_num] != sig:
                conflict = True
            signatures[map_num] = sig
    return signatures, conflict


def _detect_changed_map_numbers(
    reactant_mols: List[Any],
    product_mols: List[Any],
    *,
    min_common: int = 3,
) -> Set[int]:
    """Return map numbers with changed local environments, or empty if low-confidence."""
    reactant_maps, reactant_conflict = _collect_mapped_signatures(reactant_mols)
    product_maps, product_conflict = _collect_mapped_signatures(product_mols)
    if reactant_conflict or product_conflict:
        return set()
    common = set(reactant_maps) & set(product_maps)
    if len(common) < min_common:
        return set()
    changed = {m for m in common if reactant_maps[m] != product_maps[m]}
    return changed


def _select_primary_formed_by_mapping(
    product_bundles: List[Dict[str, Any]],
    product_mols: List[Any],
    formed_set: Set[str],
    changed_maps: Set[int],
) -> Optional[str]:
    """Select primary formed motif when its attachment atom is mapped as changed."""
    if not changed_maps:
        return None
    candidates: List[str] = []
    for idx, bundle in enumerate(product_bundles):
        mol = product_mols[idx] if idx < len(product_mols) else None
        if mol is None:
            continue
        for motif in bundle.get("motifs", []):
            if not isinstance(motif, dict):
                continue
            cid = motif.get("compound_id") or motif.get("id")
            if not cid or cid not in formed_set:
                continue
            a_idx = motif.get("a_atom_idx") or motif.get("a_idx")
            if a_idx is None:
                continue
            try:
                atom = mol.GetAtomWithIdx(int(a_idx))
            except Exception:
                continue
            map_num = atom.GetAtomMapNum()
            if map_num and map_num in changed_maps:
                candidates.append(str(cid))
    if not candidates:
        return None
    return select_primary_formed_motif(candidates, formed_set)


def _select_primary_reaction_key(
    reaction_key: Optional[str],
    reaction_keys_alt: List[str],
) -> Tuple[Optional[str], List[str]]:
    """Promote a matching alt key when the primary key has no reaction-type match."""
    if not reaction_key or not reaction_keys_alt:
        return reaction_key, reaction_keys_alt

    try:
        from chemtools.reaction_key_matcher_v2 import detect_from_reaction_key_v2
    except Exception:
        return reaction_key, reaction_keys_alt

    try:
        reacted_primary, formed_primary, spectators_primary = _parse_reaction_key_parts(reaction_key)
        primary_match, _ = detect_from_reaction_key_v2(
            reacted_primary, formed_primary, spectators_primary
        )
        if primary_match is not None:
            return reaction_key, reaction_keys_alt

        best_key = None
        best_score: Optional[Tuple[float, int]] = None
        for alt_key in reaction_keys_alt:
            reacted_alt, formed_alt, spectators_alt = _parse_reaction_key_parts(alt_key)
            alt_match, _ = detect_from_reaction_key_v2(
                reacted_alt, formed_alt, spectators_alt
            )
            if not alt_match:
                continue
            score = (
                alt_match.confidence,
                len(alt_match.matched_reacted) + len(alt_match.matched_formed),
            )
            if best_score is None or score > best_score:
                best_score = score
                best_key = alt_key

        if not best_key:
            return reaction_key, reaction_keys_alt

        new_alt = [k for k in reaction_keys_alt if k != best_key]
        new_alt.insert(0, reaction_key)
        return best_key, new_alt
    except Exception:
        return reaction_key, reaction_keys_alt


def select_primary_reacted_motifs(
    reactant_entries: Iterable[Any],
    reacted_set: set[str],
    formed_set: Optional[set[str]] = None,
) -> List[str]:
    """Select primary reacted motif from each reactant.

    Prioritize nucleophile motifs (taxonomy-driven) over electrophiles when both
    are present on the same reactant (e.g., Ar-SH vs Ar-Cl).
    """
    primary: List[str] = []
    electrophiles, nucleophiles = _load_reactant_slot_sets()
    preferred_nucleophiles: Set[str] = set()
    if formed_set:
        product_map = _load_product_to_nucleophile_map()
        for motif_id in formed_set:
            preferred_nucleophiles.update(product_map.get(motif_id, set()))
    for entry in reactant_entries or []:
        motifs = []
        if isinstance(entry, dict):
            motifs = extract_motif_ids(entry.get("motifs", []))
        else:
            motifs = extract_motif_ids(entry)
        reacted_here = [m for m in motifs if m in reacted_set]
        if reacted_here:
            idx_map = {m: idx for idx, m in enumerate(motifs)}
            if preferred_nucleophiles:
                preferred = [m for m in reacted_here if m in preferred_nucleophiles]
                if preferred:
                    preferred.sort(key=lambda m: idx_map.get(m, 0))
                    primary.append(preferred[0])
                    continue
            if electrophiles or nucleophiles:
                def _score(motif_id: str) -> tuple[int, int]:
                    if preferred_nucleophiles and motif_id in nucleophiles and motif_id not in preferred_nucleophiles:
                        return (0, idx_map.get(motif_id, 0))
                    if motif_id in nucleophiles:
                        return (2, idx_map.get(motif_id, 0))
                    if motif_id in electrophiles:
                        return (1, idx_map.get(motif_id, 0))
                    return (0, idx_map.get(motif_id, 0))
                best = sorted(reacted_here, key=lambda m: (-_score(m)[0], _score(m)[1]))[0]
                primary.append(best)
            else:
                primary.append(reacted_here[0])
    return primary


def select_primary_formed_motif(
    product_motifs: Iterable[Any],
    formed_set: set[str],
) -> Optional[str]:
    """
    Select the primary formed motif from products.
    
    Prioritizes compound motifs (scaffold-substituent like Ar-Br, HeteroAr-OH)
    over generic functional group motifs to identify the key transformation.
    
    Args:
        product_motifs: Product motif IDs
        formed_set: Set of motifs that were formed
        
    Returns:
        Primary formed motif ID, or None if no formed motifs
    """
    candidates = [m for m in extract_motif_ids(product_motifs) if m in formed_set]
    if not candidates:
        return None
    
    # Prioritize scaffold-substituent patterns (Ar-X, HeteroAr-X, Alkenyl-X, etc.)
    # over generic functional groups (R3C-X, RCH2-X, Alkyl-X, etc.)
    scaffold_substituent_motifs = []
    for m in candidates:
        if "-" not in m:
            continue
        scaffold = m.split("-", 1)[0]
        # Include only named scaffolds, exclude generic R-groups
        if not any(scaffold.startswith(prefix) for prefix in ["R", "Alkyl"]):
            scaffold_substituent_motifs.append(m)
    
    if scaffold_substituent_motifs:
        return scaffold_substituent_motifs[0]
    
    # Fall back to first formed motif
    return candidates[0]
    return None


def classify_agent_roles(agents: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
    """Classify reagents/solvents by role using reagent taxonomy."""
    from functools import lru_cache
    from chemtools.reagent.reagent_v2 import ReagentTaxonomyV2
    
    @lru_cache(maxsize=1)
    def load_taxonomy() -> Optional[ReagentTaxonomyV2]:
        try:
            return ReagentTaxonomyV2.from_path()
        except Exception:
            return None
    
    def get_agent_smiles(agent: Dict[str, Any]) -> str:
        """Extract SMILES from agent dict."""
        for key in ("smiles", "smiles_norm", "largest_smiles"):
            value = agent.get(key)
            if value:
                return str(value)
        return ""
    
    taxonomy = load_taxonomy()
    if not taxonomy:
        return {"agent_count": 0, "role_counts": {}, "family_counts": {}, "role_flags": {}, "flags": {}}
    
    role_flags = [
        "metal_catalyst", "organo_catalyst", "enzyme", "ligand", "base",
        "acid", "solvent", "additive", "oxidant", "reductant",
        "condensation_agent", "other_reagent"
    ]
    
    entries: List[Dict[str, Any]] = []
    role_counts: Dict[str, int] = {}
    family_counts: Dict[str, int] = {}
    flags: Dict[str, bool] = {role: False for role in role_flags}
    
    for agent in agents or []:
        smiles = get_agent_smiles(agent)
        if not smiles:
            continue
        
        reagent = taxonomy.lookup_reagent(smiles)
        if not reagent:
            continue
        
        role = reagent.get("role") or "other_reagent"
        family = reagent.get("family") or "Unknown"
        
        role_counts[role] = role_counts.get(role, 0) + 1
        family_counts[family] = family_counts.get(family, 0) + 1
        
        if role in flags:
            flags[role] = True
        
        entries.append({
            "smiles": smiles,
            "role": role,
            "family": family,
            "name": reagent.get("name"),
        })
    
    return {
        "entries": entries,
        "agent_count": len(entries),
        "role_counts": role_counts,
        "family_counts": family_counts,
        "role_flags": {k: v for k, v in flags.items() if v},
        "flags": flags,
    }


def featurize_reaction(
    reaction_smiles: str,
    *,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Return a reaction feature bundle (core or extended format).
    
    Args:
        reaction_smiles: Reaction SMILES with >> separator
        registry_paths: Custom taxonomy paths
        options: Featurization options
            - detailed (bool): If True, return extended format with all analysis.
                             If False (default), return simplified core format.
        
    Returns:
        Core bundle (7 fields) or extended bundle (8 fields with extended section)
    """
    options = options or {}
    use_detailed = to_bool(options.get("detailed"), default=False)
    
    include_roles = to_bool(options.get("include_roles"), default=True)
    include_agent_roles = to_bool(options.get("include_agent_roles"), default=True)
    
    # General coupling confirmation (supports 9+ coupling reaction types)
    # Backward compatibility: map old Suzuki-specific parameter to general one
    confirm_coupling = options.get("confirm_coupling_products")
    if confirm_coupling is None:
        # Check for legacy Suzuki-specific parameter
        confirm_coupling = options.get("confirm_suzuki_products")
    confirm_coupling_products = to_bool(confirm_coupling, default=True)

    # NOTE: Reaction type detection removed to avoid circular dependency
    # Detection now happens separately via chemtools.detection.detect_reaction_type()
    # which calls featurize_reaction() to extract motifs
    reaction_type = {"reaction_type": "Unknown", "confidence": 0.0, "slot_evidence": {}}
    detection_payload = {}

    # Normalize reaction SMILES
    normalized = normalize_reaction(reaction_smiles)
    reactant_smiles = [
        item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
        for item in (normalized.get("reactants") or [])
    ]
    reactant_smiles = [s for s in reactant_smiles if s]

    # Classify agents/reagents
    agent_roles = None
    if include_agent_roles:
        agent_roles = classify_agent_roles(normalized.get("agents") or [])

    # Featurize reactants
    reactant_bundles = [
        build_molecule_bundle(smiles, registry_paths=registry_paths, options=options)
        for smiles in reactant_smiles
    ]

    # Featurize products
    product_smiles = [
        item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
        for item in (normalized.get("products") or [])
    ]
    product_smiles = [s for s in product_smiles if s]
    product_bundles: List[Dict[str, Any]] = []
    product_motif_ids: List[str] = []
    product_motifs_full: List[Dict[str, Any]] = []  # Full motif dicts with fingerprints
    product_mols: List[Any] = []
    for smiles in product_smiles:
        try:
            bundle = build_molecule_bundle(smiles, registry_paths=registry_paths, options=options)
        except Exception:
            continue
        product_bundles.append(bundle)
        product_mols.append(rdkit_helpers.parse_smiles(smiles))
        product_motif_ids.extend(extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", [])))
        # Collect full motif dicts for fingerprint-aware comparison
        product_motifs_full.extend(bundle.get("motifs", []))

    # Aggregate features without pattern-based filtering
    # Pattern filtering is skipped to avoid removing reacted motifs based on incorrect initial detection
    # The validation step will correct the detection using unfiltered reacted motifs
    aggregates = aggregate_reaction_features(
        reactant_bundles,
        product_motif_ids=product_motif_ids,
        product_motifs=product_motifs_full,
        reaction_type=None,  # Disable pattern-based filtering
    )
    
    # Validate detection using reacted motifs patterns
    from .detection_validation import validate_detection_with_reacted_motifs
    
    validated = validate_detection_with_reacted_motifs(
        initial_detection=reaction_type.get("reaction_type", "Unknown") if isinstance(reaction_type, dict) else str(reaction_type),
        initial_confidence=reaction_type.get("confidence", 0.0) if isinstance(reaction_type, dict) else 0.0,
        reacted_motifs=aggregates.get("reacted_motifs", []),
        formed_motifs=aggregates.get("formed_motifs", []),
        spectator_motifs=aggregates.get("spectator_motifs", []),
    )
    
    # Update reaction type if corrected
    if validated.get("corrected_from"):
        # Preserve original detection info but update the main type
        if isinstance(reaction_type, dict):
            reaction_type["reaction_type"] = validated["reaction_type"]
            reaction_type["name"] = validated["reaction_type"]  # Also update 'name' field!
            reaction_type["confidence"] = validated["confidence"]
        else:
            reaction_type = validated["reaction_type"]
        
        # Add validation metadata to detection payload
        detection_payload["validation"] = {
            "original_detection": validated["corrected_from"],
            "validated_detection": validated["reaction_type"],
            "validation_method": validated["validation_method"],
            "validation_reason": validated["reason"],
            "validation_confidence": validated["confidence"],
        }

    # Classify reactant roles
    roles_summary = None
    if include_roles:
        try:
            roles = classify_reactants_with_context(reaction_smiles)
            roles_summary = get_reactant_summary(roles)
        except Exception:
            roles_summary = None

    intramolecular = infer_intramolecular(reactant_smiles, product_smiles, roles_summary)

    # Generate Reaction_Key using filtered aggregates
    reaction_key = None
    reaction_keys_alt: List[str] = []
    if product_bundles and reactant_bundles:
        # Use the pre-computed aggregates which have pattern-based filtering applied
        reacted = set(aggregates.get("reacted_motifs", []))
        formed = set(aggregates.get("formed_motifs", []))
        spectators = set(aggregates.get("spectator_motifs", []))
        
        # Select primary formed motif first to guide reacted motif selection.
        product_primary: List[str] = []
        for bundle in product_bundles:
            product_primary.extend(
                extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", []))
            )
        primary_formed = None
        # Mapping-aware primary formed selection (high-confidence only)
        if rdkit_helpers.rdkit_available():
            reactant_mols = [rdkit_helpers.parse_smiles(s) for s in reactant_smiles]
            changed_maps = _detect_changed_map_numbers(reactant_mols, product_mols)
            primary_formed = _select_primary_formed_by_mapping(
                product_bundles, product_mols, formed, changed_maps
            )
        if primary_formed is None:
            primary_formed = select_primary_formed_motif(product_primary, formed)

        formed_for_key = {primary_formed} if primary_formed else formed
        # Select primary motifs for key generation
        primary_reacted = select_primary_reacted_motifs(reactant_bundles, reacted, formed_for_key)

        reaction_key = format_reaction_key(
            primary_reacted if primary_reacted else reacted,
            [primary_formed] if primary_formed else formed,
            spectators,
        )
        # Alternate keys for multi-event reactions (one per formed motif)
        if formed:
            reacted_basis = primary_reacted if primary_reacted else sorted(reacted)
            for formed_motif in sorted(formed):
                alt_key = format_reaction_key(reacted_basis, [formed_motif], spectators)
                if alt_key != reaction_key and alt_key not in reaction_keys_alt:
                    reaction_keys_alt.append(alt_key)
        reaction_key, reaction_keys_alt = _select_primary_reaction_key(
            reaction_key, reaction_keys_alt
        )

    reaction = {
        "reaction_smiles": reaction_smiles,
        "normalized": normalized,
        "detection": detection_payload,
        "reaction_type": reaction_type,
        "reactants": reactant_bundles,
        "products": product_bundles,
        "aggregates": aggregates,
        "reaction_key": reaction_key,
        "reaction_keys_alt": reaction_keys_alt,
        "roles": roles_summary,
        "agent_roles": agent_roles,
        "intramolecular": intramolecular,
    }

    # Add feasibility analysis for specific reaction types
    rt_id = reaction_type.get("reaction_type")
    if rt_id == "snar_cn" or rt_id == "c_n_cross_coupling":
        reaction["feasibility"] = analyze_snar_feasibility(reaction)

    # Return simplified format (core or extended)
    if use_detailed:
        result = build_extended_reaction(reaction)
    else:
        result = build_core_reaction(reaction)
    
    # Add metadata
    result["kind"] = "reaction"
    result["schema_version"] = "v2"
    
    meta = {
        "rdkit_available": rdkit_helpers.rdkit_available(),
    }
    if detection_payload.get("error"):
        meta["errors"] = [detection_payload["error"]]
    if meta.get("errors") or not meta.get("rdkit_available", True):
        result["meta"] = meta
    
    return result
