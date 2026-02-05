"""
Reaction-level feature formatting and bundling.

Handles reaction type detection, reactant/product processing, and reaction key generation.
"""

from __future__ import annotations

from functools import lru_cache
import json
import re
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


def _map_atom_labels(mapped_smiles: str) -> Dict[int, str]:
    """Build a map_num -> label mapping from a mapped reaction SMILES."""
    if not mapped_smiles or ">>" not in mapped_smiles:
        return {}
    labels: Dict[int, str] = {}
    reactants_smiles, products_smiles = mapped_smiles.split(">>", 1)
    for side in (reactants_smiles, products_smiles):
        for part in side.split("."):
            mol = rdkit_helpers.parse_smiles(part)
            if mol is None:
                continue
            for atom in mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if not map_num or map_num in labels:
                    continue
                symbol = atom.GetSymbol()
                if atom.GetIsAromatic():
                    labels[map_num] = f"{symbol}(ar)"
                else:
                    labels[map_num] = symbol
    return labels


def _format_bond_tokens(
    bonds: Iterable[Any],
    labels: Dict[int, str],
) -> List[str]:
    tokens: List[str] = []
    seen: Set[Tuple[str, str]] = set()
    for bond in bonds or []:
        if not isinstance(bond, (tuple, list)) or len(bond) < 2:
            continue
        a, b = bond[0], bond[1]
        if isinstance(a, int):
            a_label = labels.get(a, f"#{a}")
        else:
            a_label = str(a).split()[0]
        if isinstance(b, int):
            b_label = labels.get(b, f"#{b}")
        else:
            b_label = str(b).split()[0]
        pair = tuple(sorted((a_label, b_label)))
        if pair in seen:
            continue
        seen.add(pair)
        tokens.append(f"{pair[0]}-{pair[1]}")
    return tokens


def _get_bond_change_analysis(reaction_smiles: str) -> Optional[Dict[str, Any]]:
    if not reaction_smiles or not rdkit_helpers.rdkit_available():
        return None
    try:
        from chemtools._atom_mapping import analyze_bond_changes_hybrid
    except Exception:
        return None
    try:
        analysis = analyze_bond_changes_hybrid(reaction_smiles, use_mcs=False)
    except Exception:
        return None
    if not analysis or not analysis.get("success"):
        return None
    result = analysis.get("recommended_result") or analysis.get("rxnmapper_result") or analysis.get("manual_result")
    if not result or not result.get("success"):
        return None
    return result


@lru_cache(maxsize=1)
def _load_compound_logic_sets() -> Dict[str, Set[str]]:
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


def format_bond_change_key(
    reaction_smiles: str,
    *,
    analysis: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """Generate a bond-change key from atom mapping (POC)."""
    result = analysis or _get_bond_change_analysis(reaction_smiles)
    if not result:
        return None
    mapped_smiles = result.get("mapped_smiles") or ""
    labels = _map_atom_labels(mapped_smiles)
    broken = _format_bond_tokens(result.get("broken_bonds"), labels)
    formed = _format_bond_tokens(result.get("formed_bonds"), labels)
    if not broken and not formed:
        return None
    parts: List[str] = []
    if broken:
        parts.append("break: " + "; ".join(broken))
    if formed:
        parts.append("form: " + "; ".join(formed))
    return " | ".join(parts)


def _extract_bond_pairs_from_key(bond_key: Optional[str]) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    formed_pairs: List[Tuple[str, str]] = []
    broken_pairs: List[Tuple[str, str]] = []
    for token in _extract_bond_section(bond_key, section="form"):
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        formed_pairs.append((left, right))
    for token in _extract_bond_section(bond_key, section="break"):
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        broken_pairs.append((left, right))
    return formed_pairs, broken_pairs


def _pair_reactants_from_bonds(
    reacted: Iterable[str],
    bond_key: Optional[str],
) -> List[str]:
    """Heuristic pairing of nucleophile/electrophile using bond-change tokens."""
    if not bond_key:
        return []
    formed_pairs, broken_pairs = _extract_bond_pairs_from_key(bond_key)
    if not formed_pairs:
        return []

    def pick_prefix(token: str, prefixes: Tuple[str, ...]) -> Optional[str]:
        for pref in prefixes:
            if token.endswith(pref):
                return pref
        return None

    halides = ("-I", "-Br", "-Cl", "-F")
    nucleophiles = ("-SH", "-NH2", "-NHR", "-OH")

    reacted_list = [str(r) for r in reacted if r]
    pairs: List[str] = []

    for formed in formed_pairs:
        # Try to detect which heteroatom formed the bond (S/N/O).
        left, right = formed
        hetero = None
        if _is_sulfur_label(left) or _is_sulfur_label(right):
            hetero = "S"
        elif left.startswith("N") or right.startswith("N"):
            hetero = "N"
        elif left.startswith("O") or right.startswith("O"):
            hetero = "O"
        if not hetero:
            continue

        nuc = None
        for motif in reacted_list:
            if hetero == "S" and motif.endswith("-SH"):
                nuc = motif
                break
            if hetero == "N" and (motif.endswith("-NH2") or motif.endswith("-NHR")):
                nuc = motif
                break
            if hetero == "O" and motif.endswith("-OH"):
                nuc = motif
                break

        if not nuc:
            continue

        # Match electrophile by halide in broken bonds.
        elec = None
        for broken in broken_pairs:
            b_left, b_right = broken
            if any(x in b_left for x in ("I", "Br", "Cl", "F")) or any(x in b_right for x in ("I", "Br", "Cl", "F")):
                for motif in reacted_list:
                    if any(motif.endswith(h) for h in halides):
                        elec = motif
                        break
            if elec:
                break

        if nuc and elec:
            pairs.append(f"{nuc}~{elec}")

    # De-dup
    out: List[str] = []
    seen: Set[str] = set()
    for p in pairs:
        if p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out


def _extract_bond_section(bond_key: Optional[str], *, section: str) -> List[str]:
    if not bond_key:
        return []
    marker = f"{section}: "
    parts = [p.strip() for p in bond_key.split(" | ") if p.strip()]
    for part in parts:
        if part.startswith(marker):
            payload = part[len(marker):].strip()
            if not payload:
                return []
            return [tok.strip() for tok in payload.split(";") if tok.strip()]
    return []


def _is_carbon_label(label: str) -> bool:
    return label == "C" or label.startswith("C(")


def _is_sulfur_label(label: str) -> bool:
    return label == "S" or label.startswith("S(")


def _parse_element(label: str) -> str:
    label = label.strip()
    if label.startswith("#"):
        return ""
    if "(" in label:
        label = label.split("(", 1)[0]
    return label


def _is_aromatic_label(label: str) -> bool:
    return label.endswith("(ar)")


def _bond_elements_from_key(bond_key: Optional[str]) -> Set[str]:
    elements: Set[str] = set()
    for section in ("break", "form"):
        for token in _extract_bond_section(bond_key, section=section):
            if "-" not in token:
                continue
            left, right = [t.strip() for t in token.split("-", 1)]
            left_el = _parse_element(left)
            right_el = _parse_element(right)
            if left_el:
                elements.add(left_el)
            if right_el:
                elements.add(right_el)
    return elements


def _match_group_id(motif_id: str, group_element_map: Dict[str, Set[str]]) -> Optional[str]:
    motif_id = str(motif_id)
    for group_id in sorted(group_element_map.keys(), key=len, reverse=True):
        if motif_id.endswith(group_id):
            return group_id
    return None


def _filter_reactants_for_crk(
    reacted: Iterable[str],
    bond_key: Optional[str],
    spectators: Optional[Iterable[str]] = None,
) -> List[str]:
    scaffold_ids: Set[str] = set()
    try:
        from .aggregation import load_scaffold_motif_ids
        scaffold_ids = load_scaffold_motif_ids()
    except Exception:
        scaffold_ids = set()

    if not bond_key:
        return sorted(str(r) for r in reacted if r and str(r) not in scaffold_ids)
    elements = _bond_elements_from_key(bond_key)
    if not elements:
        return sorted(str(r) for r in reacted if r and str(r) not in scaffold_ids)
    elements_non_c = {el for el in elements if el != "C"}

    # If N/O/S bonds are formed, make sure N/O/S nucleophiles are not dropped
    # just because they have -H suffixes (e.g., AromN-H).
    formed_tokens = _extract_bond_section(bond_key, section="form")
    formed_elements: Set[str] = set()
    for token in formed_tokens:
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        if left_el:
            formed_elements.add(left_el)
        if right_el:
            formed_elements.add(right_el)

    bond_tokens = _extract_bond_section(bond_key, section="form") + _extract_bond_section(
        bond_key, section="break"
    )

    def _parse_bond_token(token: str) -> Optional[Tuple[str, bool, str, bool]]:
        if "-" not in token:
            return None
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        if not left_el or not right_el:
            return None
        return (left_el, _is_aromatic_label(left), right_el, _is_aromatic_label(right))

    parsed_bonds = [entry for entry in (_parse_bond_token(t) for t in bond_tokens) if entry]

    def _has_c_c_bond() -> bool:
        for left_el, _left_ar, right_el, _right_ar in parsed_bonds:
            if left_el == "C" and right_el == "C":
                return True
        return False

    def _has_non_aromatic_c_bond(other_elements: Set[str]) -> bool:
        for left_el, left_ar, right_el, right_ar in parsed_bonds:
            if left_el == "C" and not left_ar and right_el in other_elements:
                return True
            if right_el == "C" and not right_ar and left_el in other_elements:
                return True
        return False

    group_element_map = {
        "-I": {"I"},
        "-Br": {"Br"},
        "-Cl": {"Cl"},
        "-F": {"F"},
        "-SH": {"S"},
        "-OH": {"O"},
        "-NH2": {"N"},
        "-NHR": {"N"},
        "-NR2": {"N"},
        "-NHCOR": {"N"},
        "-OR": {"O"},
        "-SR": {"S"},
        "-N3": {"N"},
        "-CN": {"C", "N"},
        "-COCl": {"C", "Cl"},
        "-COBr": {"C", "Br"},
        "-COI": {"C", "I"},
        "-COF": {"C", "F"},
        "-CO2H": {"C", "O"},
        "-CO2R": {"C", "O"},
        "-B(OH)2": {"B"},
        "-Bpin": {"B"},
        "-B(OR)2": {"B"},
        "-BF3K": {"B"},
        "-Sn*": {"Sn"},
        "-Zn*": {"Zn"},
        "-Mg*": {"Mg"},
        "-Si*": {"Si"},
        "-H": {"H"},
    }

    spectator_ids = {str(s) for s in (spectators or []) if s}
    filtered: List[str] = []
    for motif in reacted:
        if not motif:
            continue
        if str(motif) in scaffold_ids:
            continue
        motif_str = str(motif)
        scaffold = motif_str.split("-", 1)[0] if "-" in motif_str else motif_str
        # Keep N/O/S nucleophiles when corresponding bonds form.
        keep_by_formed = False
        if "N" in formed_elements and (
            motif_str.endswith("-NH2")
            or motif_str.endswith("-NHR")
            or motif_str.endswith("-NR2")
            or motif_str.endswith("-NHCOR")
            or motif_str == "AromN-H"
        ):
            keep_by_formed = True
        if "O" in formed_elements and (motif_str.endswith("-OH") or motif_str.endswith("-OR")):
            keep_by_formed = True
        if "S" in formed_elements and (motif_str.endswith("-SH") or motif_str.endswith("-SR")):
            keep_by_formed = True
        if scaffold.startswith("Alkenyl") or scaffold.startswith("Alkynyl"):
            if _has_c_c_bond():
                keep_by_formed = True

        group_id = _match_group_id(motif_str, group_element_map)
        keep_by_elements = False
        if group_id:
            group_elements = group_element_map.get(group_id)
            if group_elements:
                if elements_non_c:
                    if elements_non_c.intersection(group_elements):
                        keep_by_elements = True
                elif elements.intersection(group_elements):
                    keep_by_elements = True

        # For carbonyl-derived groups, element overlap alone is too permissive.
        # Require explicit nucleophile logic to keep these when they are spectators.
        if group_id in {
            "-CO2R",
            "-CO2H",
            "-COR",
            "-CONH2",
            "-CONHR",
            "-CONR2",
        }:
            # Carbonyl-derived motifs should only be kept if a non-aromatic carbonyl bond changes.
            if not _has_non_aromatic_c_bond({"O", "N", "Cl", "Br", "I", "F", "S"}):
                keep_by_elements = False

        # Drop pure spectators unless they align with bond elements or formed nucleophiles.
        if motif_str in spectator_ids and not (keep_by_formed or keep_by_elements):
            continue

        if keep_by_formed:
            filtered.append(motif_str)
            continue
        if not group_id:
            filtered.append(motif_str)
            continue
        if keep_by_elements:
            filtered.append(motif_str)
    filtered_sorted = sorted(filtered)
    if filtered_sorted:
        return filtered_sorted
    fallback = [str(r) for r in reacted if r and str(r) not in scaffold_ids]
    if fallback:
        return sorted(fallback)
    return sorted(str(r) for r in reacted if r)


def _strip_atom_mapping(smiles: str) -> str:
    return re.sub(r":\d+", "", smiles)


def _build_map_info(
    mapped_smiles: str,
    reactant_smiles: List[str],
) -> Dict[int, Dict[str, Any]]:
    """Map atom map numbers to element/aromatic/reactant index."""
    if not mapped_smiles or ">>" not in mapped_smiles:
        return {}
    reactants_side = mapped_smiles.split(">>", 1)[0]
    components = [c for c in reactants_side.split(".") if c]

    # Canonical reactant smiles for matching
    canon_reactants: List[str] = []
    for smi in reactant_smiles:
        canon = rdkit_helpers.canonical_smiles(smi) or smi
        canon_reactants.append(canon)
    remaining: Dict[str, List[int]] = {}
    for idx, smi in enumerate(canon_reactants):
        remaining.setdefault(smi, []).append(idx)

    info: Dict[int, Dict[str, Any]] = {}
    for comp in components:
        comp_unmapped = _strip_atom_mapping(comp)
        comp_canon = rdkit_helpers.canonical_smiles(comp_unmapped) or comp_unmapped
        react_idx = None
        if comp_canon in remaining and remaining[comp_canon]:
            react_idx = remaining[comp_canon].pop(0)
        mol = rdkit_helpers.parse_smiles(comp)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if not map_num:
                continue
            info[map_num] = {
                "element": atom.GetSymbol(),
                "aromatic": bool(atom.GetIsAromatic()),
                "reactant_idx": react_idx,
            }
    return info


def _find_halide_for_carbon(
    carbon_map: int,
    broken_bonds: Iterable[Any],
    map_info: Dict[int, Dict[str, Any]],
) -> Optional[str]:
    priority = {"I": 3, "Br": 2, "Cl": 1, "F": 0}
    found: List[str] = []
    for bond in broken_bonds or []:
        if not isinstance(bond, (tuple, list)) or len(bond) < 2:
            continue
        a, b = bond[0], bond[1]
        if a == carbon_map:
            other = b
        elif b == carbon_map:
            other = a
        else:
            continue

        if isinstance(other, int):
            element = map_info.get(other, {}).get("element")
        else:
            element = _parse_element(str(other))
        if element in priority:
            found.append(element)
    if not found:
        return None
    found.sort(key=lambda el: -priority.get(el, 0))
    return found[0]


def _reacted_motifs_by_reactant(
    reactant_bundles: Iterable[Dict[str, Any]],
    reacted_set: Set[str],
) -> List[List[str]]:
    buckets: List[List[str]] = []
    for bundle in reactant_bundles or []:
        motifs = extract_motif_ids(bundle.get("motifs", []))
        buckets.append([m for m in motifs if m in reacted_set])
    return buckets


def _pair_reactants_from_mapping(
    reacted_set: Set[str],
    *,
    analysis: Optional[Dict[str, Any]],
    reactant_smiles: List[str],
    reactant_bundles: List[Dict[str, Any]],
) -> List[str]:
    if not analysis:
        return []
    mapped_smiles = analysis.get("mapped_smiles") or ""
    broken_bonds = analysis.get("broken_bonds") or []
    formed_bonds = analysis.get("formed_bonds") or []
    if not mapped_smiles or not formed_bonds:
        return []

    map_info = _build_map_info(mapped_smiles, reactant_smiles)
    if not map_info:
        return []

    reacted_by_idx = _reacted_motifs_by_reactant(reactant_bundles, reacted_set)
    pairs: List[str] = []

    for bond in formed_bonds:
        if not isinstance(bond, (tuple, list)) or len(bond) < 2:
            continue
        a, b = bond[0], bond[1]
        if not isinstance(a, int) or not isinstance(b, int):
            continue
        a_info = map_info.get(a)
        b_info = map_info.get(b)
        if not a_info or not b_info:
            continue
        # Identify hetero vs carbon.
        if a_info["element"] == "C" and b_info["element"] in {"S", "N", "O"}:
            carbon_map, hetero_map = a, b
            hetero_el = b_info["element"]
        elif b_info["element"] == "C" and a_info["element"] in {"S", "N", "O"}:
            carbon_map, hetero_map = b, a
            hetero_el = a_info["element"]
        else:
            continue

        nuc_idx = map_info.get(hetero_map, {}).get("reactant_idx")
        elec_idx = map_info.get(carbon_map, {}).get("reactant_idx")
        if nuc_idx is None or elec_idx is None:
            continue
        if nuc_idx >= len(reacted_by_idx) or elec_idx >= len(reacted_by_idx):
            continue

        nuc_candidates = reacted_by_idx[nuc_idx]
        elec_candidates = reacted_by_idx[elec_idx]
        if not nuc_candidates or not elec_candidates:
            continue

        if hetero_el == "S":
            nuc = next((m for m in nuc_candidates if m.endswith("-SH")), nuc_candidates[0])
        elif hetero_el == "N":
            nuc = next((m for m in nuc_candidates if m.endswith("-NH2") or m.endswith("-NHR")), nuc_candidates[0])
        elif hetero_el == "O":
            nuc = next((m for m in nuc_candidates if m.endswith("-OH")), nuc_candidates[0])
        else:
            nuc = nuc_candidates[0]

        halide = _find_halide_for_carbon(carbon_map, broken_bonds, map_info)
        if halide:
            elec = next((m for m in elec_candidates if m.endswith(f"-{halide}")), elec_candidates[0])
        else:
            elec = elec_candidates[0]

        if nuc and elec:
            pairs.append(f"{nuc}~{elec}")

    out: List[str] = []
    seen: Set[str] = set()
    for p in pairs:
        if p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out


def _fallback_pairs_by_halide_priority(reacted: Iterable[str]) -> List[str]:
    """Fallback pairing when mapping is unavailable."""
    reacted_list = [str(r) for r in reacted if r]
    if not reacted_list:
        return []
    nucleophile = next((m for m in reacted_list if m.endswith("-SH")), None)
    if not nucleophile:
        nucleophile = next((m for m in reacted_list if m.endswith("-NH2") or m.endswith("-NHR")), None)
    if not nucleophile:
        nucleophile = next((m for m in reacted_list if m.endswith("-OH")), None)
    if not nucleophile:
        return []
    halide_priority = ["-I", "-Br", "-Cl", "-F"]
    electrophile = None
    for h in halide_priority:
        electrophile = next((m for m in reacted_list if m.endswith(h)), None)
        if electrophile:
            break
    if not electrophile:
        return []
    return [f"{nucleophile}~{electrophile}"]


def _infer_product_broad_tags(bond_key: Optional[str]) -> List[str]:
    tags: List[str] = []
    formed = _extract_bond_section(bond_key, section="form")
    if not formed:
        return tags
    for token in formed:
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        if (_is_carbon_label(left) and _is_sulfur_label(right)) or (
            _is_sulfur_label(left) and _is_carbon_label(right)
        ):
            tags.append("Product_C-S")
            if left.endswith("(ar)") or right.endswith("(ar)"):
                tags.append("Product_Aryl_S")
            break
    return sorted(set(tags))


def _select_primary_broad_tag(tags: Iterable[str]) -> str:
    if not tags:
        return "[]"
    priority = {
        "Product_Aryl_S": 2,
        "Product_C-S": 1,
    }
    candidates = [str(t) for t in tags if t]
    if not candidates:
        return "[]"
    candidates.sort(key=lambda t: (-priority.get(t, 0), t))
    return candidates[0]


def _detect_aryl_s_from_product_smiles(product_smiles: Iterable[str]) -> bool:
    if not rdkit_helpers.rdkit_available():
        return False
    # Aromatic carbon bonded to sulfur (aryl thioether/aryl sulfur link)
    # Matches c-S- and c-S(=O)- variants loosely.
    smarts = "[c][S]"
    try:
        from chemtools.util.smarts_cache import compile_smarts
    except Exception:
        return False
    patt = compile_smarts(smarts, validate=False)
    if patt is None:
        return False
    for smi in product_smiles or []:
        mol = rdkit_helpers.parse_smiles(smi)
        if mol is None:
            continue
        try:
            if mol.HasSubstructMatch(patt):
                return True
        except Exception:
            continue
    return False


def _infer_product_broad_tags_with_validation(
    *,
    bond_key: Optional[str],
    product_smiles: Iterable[str],
) -> List[str]:
    """Infer broad product tags with mapping-first + SMARTS fallback."""
    tags = _infer_product_broad_tags(bond_key)
    if "Product_Aryl_S" in tags:
        return tags
    # Fallback: product structure contains aryl-S bond.
    if _detect_aryl_s_from_product_smiles(product_smiles):
        tags.append("Product_C-S")
        tags.append("Product_Aryl_S")
    return sorted(set(tags))


def format_crk_key(
    *,
    bond_key: Optional[str],
    reacted: Iterable[str],
    spectators: Iterable[str],
    product_broad_tags: Iterable[str],
    product_motifs_reactive: Optional[Iterable[str]] = None,
    include_product: bool = False,
) -> str:
    """Build a composite Condition Recommendation Key (CRK-v1)."""
    reactants_text = _format_motif_list(reacted)
    broad_tags = sorted(str(t) for t in product_broad_tags if t)
    reactive_products = sorted(str(t) for t in (product_motifs_reactive or []) if t)

    products_primary = "[]"
    if include_product:
        if reactive_products:
            products_primary = "|".join(reactive_products)
        else:
            products_primary = _select_primary_broad_tag(broad_tags)

    summary = f"|{reactants_text} -> {products_primary}"

    sections: List[str] = [summary]

    if bond_key:
        formed = _extract_bond_section(bond_key, section="form")
        broken = _extract_bond_section(bond_key, section="break")
        if formed:
            sections.append("bond_formed: " + "; ".join(formed))
        if broken:
            sections.append("bond_broken: " + "; ".join(broken))

    if spectators:
        sections.append("spectators: " + _format_motif_list(spectators))

    return " | ".join(sections)


def _format_motif_list(items: Iterable[str]) -> str:
    values = sorted(str(item) for item in items if item)
    return "|".join(values) if values else "[]"


def _select_reactive_product_motifs(
    product_motifs: Iterable[Dict[str, Any]],
    *,
    bond_key: Optional[str],
    formed_motifs: Iterable[str],
    reacted_motifs: Iterable[str],
) -> List[str]:
    """Select product motifs that align with bond formation (reaction-center motifs)."""
    motif_ids = [
        str(m.get("compound_id") or m.get("id"))
        for m in product_motifs
        if isinstance(m, dict) and (m.get("compound_id") or m.get("id"))
    ]
    motif_set = set(motif_ids)
    formed_set = {str(m) for m in formed_motifs if m}
    reacted_set = {str(m) for m in reacted_motifs if m}

    inferred = _infer_product_motifs_from_logic(reacted_set, bond_key)

    if not bond_key:
        if inferred:
            return inferred
        if formed_set:
            return sorted(motif_set & formed_set) or sorted(formed_set)
        return []

    formed_bonds = _extract_bond_section(bond_key, section="form")
    if not formed_bonds:
        if inferred:
            return inferred
        if formed_set:
            return sorted(motif_set & formed_set) or sorted(formed_set)
        return []

    logic_sets = _load_compound_logic_sets()
    target_ids: Set[str] = set()

    def add_set(name: str) -> None:
        target_ids.update(logic_sets.get(name, set()))

    for token in formed_bonds:
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        left_ar = _is_aromatic_label(left)
        right_ar = _is_aromatic_label(right)

        elements = {left_el, right_el}
        if "C" in elements and "N" in elements:
            add_set("aryl_amines")
            add_set("aryl_amides")
        if "C" in elements and "O" in elements:
            add_set("aryl_ethers")
        if "C" in elements and "S" in elements:
            add_set("aryl_thioethers")
        if left_el == "C" and right_el == "C" and left_ar and right_ar:
            target_ids.add("Ar-Ar")

    reactive = sorted(motif_set & target_ids)
    if inferred:
        if reactive:
            overlap = sorted(set(reactive) & set(inferred))
            if overlap:
                return overlap
        return inferred
    if reactive:
        return reactive
    if formed_set:
        return sorted(motif_set & formed_set) or sorted(formed_set)
    return []


def _infer_product_motifs_from_logic(
    reacted_motifs: Iterable[str],
    bond_key: Optional[str],
) -> List[str]:
    """Infer likely product motifs from reacted motifs + bond formation."""
    if not bond_key:
        return []
    reacted_list = [str(m) for m in reacted_motifs if m]
    if not reacted_list:
        return []
    formed_bonds = _extract_bond_section(bond_key, section="form")
    if not formed_bonds:
        return []

    logic_sets = _load_compound_logic_sets()
    amines = logic_sets.get("aryl_amines", set())
    amides = logic_sets.get("aryl_amides", set())
    ethers = logic_sets.get("aryl_ethers", set())
    thioethers = logic_sets.get("aryl_thioethers", set())

    def has_suffix(suffixes: Tuple[str, ...]) -> bool:
        return any(m.endswith(sfx) for m in reacted_list for sfx in suffixes)

    def has_any(values: Set[str]) -> bool:
        if not values:
            return False
        for m in reacted_list:
            if m in values:
                return True
        return False

    inferred: Set[str] = set()
    for token in formed_bonds:
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        left_ar = _is_aromatic_label(left)
        right_ar = _is_aromatic_label(right)
        elements = {left_el, right_el}
        is_aryl_carbon = (left_el == "C" and left_ar) or (right_el == "C" and right_ar)

        if "C" in elements and "S" in elements and is_aryl_carbon:
            if has_any(logic_sets.get("thiols_sh", set())) or has_suffix(("-SH",)):
                if "Ar-SR" in thioethers:
                    inferred.add("Ar-SR")

        if "C" in elements and "O" in elements and is_aryl_carbon:
            if has_any(logic_sets.get("alcohols_oh", set())) or has_suffix(("-OH",)):
                if "Ar-OR" in ethers:
                    inferred.add("Ar-OR")

        if "C" in elements and "N" in elements and is_aryl_carbon:
            n_is_aromatic = (left_el == "N" and left_ar) or (right_el == "N" and right_ar)
            if n_is_aromatic:
                if "Ar-AromN" in amines:
                    inferred.add("Ar-AromN")
                continue

            if has_suffix(("-NHCOR",)):
                if "Ar-NRCOR" in amides:
                    inferred.add("Ar-NRCOR")
                continue
            if has_suffix(("-NH2",)):
                if "Ar-NH2" in amines:
                    inferred.add("Ar-NH2")
                continue
            if has_suffix(("-NHR",)):
                if "Ar-NHR" in amines:
                    inferred.add("Ar-NHR")
                continue
            if has_suffix(("-NR2",)):
                if "Ar-NR2" in amines:
                    inferred.add("Ar-NR2")
                continue

            if "AromN-H" in reacted_list and "Ar-AromN" in amines:
                inferred.add("Ar-AromN")

    return sorted(inferred)


def _scaffold_spectators_from_bundles(
    reactant_bundles: Iterable[Dict[str, Any]],
    product_bundles: Iterable[Dict[str, Any]],
) -> Set[str]:
    """Identify scaffold motifs present on both sides to treat as spectators in CRK."""
    try:
        from .aggregation import load_scaffold_motif_ids
    except Exception:
        return set()
    scaffold_ids = load_scaffold_motif_ids()
    if not scaffold_ids:
        return set()
    reactant_ids: Set[str] = set()
    product_ids: Set[str] = set()
    for bundle in reactant_bundles or []:
        reactant_ids.update(
            extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", []))
        )
    for bundle in product_bundles or []:
        product_ids.update(
            extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", []))
        )
    return set(m for m in reactant_ids & product_ids if m in scaffold_ids)


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
    skip_bond_analysis = to_bool(options.get("skip_bond_analysis"), default=False)
    include_product_in_crk = to_bool(options.get("include_product_in_crk"), default=True)
    
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
    for smiles in product_smiles:
        try:
            bundle = build_molecule_bundle(smiles, registry_paths=registry_paths, options=options)
        except Exception:
            continue
        product_bundles.append(bundle)
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

    rt_id = None
    if isinstance(reaction_type, dict):
        rt_id = reaction_type.get("reaction_type")
    elif reaction_type is not None:
        rt_id = str(reaction_type)
    if rt_id == "Unknown":
        rt_id = None

    # Generate CRK-v1 reaction key (single source of truth)
    reaction_key = None
    product_broad_tags: List[str] = []
    product_motifs_reactive: List[str] = []
    if product_bundles and reactant_bundles:
        reacted_full = set(aggregates.get("reacted_motifs", []))
        spectators = set(aggregates.get("spectator_motifs", []))
        scaffold_spectators = _scaffold_spectators_from_bundles(reactant_bundles, product_bundles)
        spectators_for_crk = spectators | scaffold_spectators

        bond_analysis = None
        bond_key = None
        if not skip_bond_analysis:
            bond_analysis = _get_bond_change_analysis(reaction_smiles)
            bond_key = format_bond_change_key(reaction_smiles, analysis=bond_analysis)
        product_broad_tags = _infer_product_broad_tags_with_validation(
            bond_key=bond_key,
            product_smiles=product_smiles,
        )
        product_motifs_reactive = _select_reactive_product_motifs(
            product_motifs_full,
            bond_key=bond_key,
            formed_motifs=aggregates.get("formed_motifs", []),
            reacted_motifs=reacted_full,
        )
        reacted_for_crk = _filter_reactants_for_crk(
            reacted_full,
            bond_key,
            spectators=spectators_for_crk,
        )
        if not reacted_for_crk and reacted_full:
            reacted_for_crk = sorted(str(r) for r in reacted_full if r)
        reaction_key = format_crk_key(
            bond_key=bond_key,
            reacted=reacted_for_crk,
            spectators=spectators_for_crk,
            product_broad_tags=product_broad_tags,
            product_motifs_reactive=product_motifs_reactive,
            include_product=include_product_in_crk,
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
        "product_broad_tags": product_broad_tags,
        "product_motifs_reactive": product_motifs_reactive,
        "roles": roles_summary,
        "agent_roles": agent_roles,
        "intramolecular": intramolecular,
    }

    # Add feasibility analysis for specific reaction types
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
