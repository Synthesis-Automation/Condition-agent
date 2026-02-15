"""
Deterministic reaction event extraction and reaction-key quality assessment.

This module provides a general-purpose, taxonomy-aligned evidence layer that
summarizes what transformation primitives occurred (e.g., C-N formation,
leaving-group displacement, ring closure) and flags anomalies suggesting
multi-step records or missing activation metadata.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from chemtools.util import rdkit_helpers


_HALOGENS = {"F", "Cl", "Br", "I"}
_PATTERN_SMARTS = {
    "ester": "[CX3](=O)[OX2][#6]",
    "carboxy_acid": "[CX3](=O)[OX2H]",
    "carboxylate": "[CX3](=O)[O-]",
    "phenol": "[c][OX2H]",
    "benzyl_halide": "[c][CH2][Cl,Br,I]",
    "benzyl_aryl_ether": "[c][CH2][O][c]",
}
_EVENT_SIGNATURE_PRIORITY = [
    "benzyl_o_alkylation_like",
    "ester_hydrolysis_like",
    "amidation_like",
    "ring_closure_or_annulation",
    "leaving_group_displacement",
    "c_n_bond_formation",
    "c_o_bond_formation",
    "c_s_bond_formation",
    "c_c_bond_formation",
]
_EVENT_SIGNATURE_CODE = {
    "benzyl_o_alkylation_like": "BzOAlk",
    "ester_hydrolysis_like": "EsterHyd",
    "amidation_like": "Amid",
    "ring_closure_or_annulation": "Ann",
    "leaving_group_displacement": "LGDisp",
    "c_n_bond_formation": "C-N",
    "c_o_bond_formation": "C-O",
    "c_s_bond_formation": "C-S",
    "c_c_bond_formation": "C-C",
}
_REDOX_HETERO = {"N", "O", "S", "P", "F", "Cl", "Br", "I"}
_SCAFFOLD_SP3_PREFIXES = (
    "Alkyl-",
    "CH3-",
    "RCH2-",
    "R2CH-",
    "R3C-",
    "Bn-",
    "Allyl-",
    "Propargyl-",
)
_SCAFFOLD_SP2_PREFIXES = ("Ar-", "HeteroAr-", "AromN-", "Alkenyl-")
_SCAFFOLD_SP_PREFIXES = ("Alkynyl-",)
_EWG_TOKENS = ("NO2", "CN", "CO2", "COR", "CHO", "SO2", "CF3")
_SNAR_ACTIVATING_TOKENS = (
    "HeteroAr",
    "AromN",
    "Pyridine",
    "Pyridyl",
    "Pyrimidine",
    "Pyrimidyl",
    "Pyrazine",
    "Triazine",
    "Quinoline",
    "Isoquinoline",
    "-NO2",
    "-CN",
    "-COR",
    "-CO2",
    "-SO2",
    "-CHO",
    "-CF3",
    "-N2+",
)
_SNAR_DEACTIVATING_TOKENS = (
    "Ar-OR",
    "Ar-NR2",
    "Ar-NHR",
    "Ar-NH2",
    "Ar-SR",
    "Ar-Alkyl",
)
_AMBIDENT_TOKENS = ("-SCN", "-NCS", "-NO2", "-CN", "Hydrazine")
_NUCLEOPHILE_CLASS_HINTS = {
    "amine": ("-NH2", "-NHR", "-NR2", "Hydrazine"),
    "alcohol_or_phenol": ("-OH", "-OR"),
    "thiol_or_thioether": ("-SH", "-SR"),
    "cyanide_like": ("-CN", "-NC", "-SCN", "-NCS"),
    "carbon_nucleophile_like": ("R_acidic", "Alkenyl-", "Alkynyl-"),
}
_OA_PARTNER_HINTS = ("-B(", "-BF", "Bpin", "-Zn", "-Sn", "-Mg", "-Si", "-Li", "Grignard")


def _split_reaction_sides(reaction_smiles: str) -> Tuple[List[str], List[str]]:
    if not reaction_smiles:
        return [], []
    text = str(reaction_smiles).strip()
    if ">>" in text:
        left, right = text.split(">>", 1)
        return [x for x in left.split(".") if x], [x for x in right.split(".") if x]
    parts = text.split(">")
    if len(parts) == 3:
        left = parts[0]
        right = parts[2]
        return [x for x in left.split(".") if x], [x for x in right.split(".") if x]
    return [x for x in text.split(".") if x], []


def _extract_bond_section(bond_key: Optional[str], section: str) -> List[str]:
    if not bond_key:
        return []
    marker = f"{section}: "
    parts = [p.strip() for p in str(bond_key).split(" | ") if p.strip()]
    for part in parts:
        if part.startswith(marker):
            payload = part[len(marker):].strip()
            if not payload:
                return []
            return [tok.strip() for tok in payload.split(";") if tok.strip()]
    return []


def _parse_element(label: str) -> str:
    text = str(label).strip()
    if not text:
        return ""
    if text.startswith("#"):
        return ""
    if "(" in text:
        text = text.split("(", 1)[0]
    if "-" in text and len(text) > 2:
        # tokens like "C(ar)-N" are split earlier; this is a safety guard.
        text = text.split("-", 1)[0]
    return text


def _parse_bond_token(token: str) -> Optional[Tuple[str, str]]:
    if "-" not in str(token):
        return None
    left, right = [t.strip() for t in str(token).split("-", 1)]
    left_el = _parse_element(left)
    right_el = _parse_element(right)
    if not left_el or not right_el:
        return None
    pair = tuple(sorted((left_el, right_el)))
    return pair  # type: ignore[return-value]


def _bond_pairs(tokens: Iterable[str]) -> Set[Tuple[str, str]]:
    pairs: Set[Tuple[str, str]] = set()
    for token in tokens or []:
        parsed = _parse_bond_token(token)
        if parsed is not None:
            pairs.add(parsed)
    return pairs


def _ring_count(smiles_list: Iterable[str]) -> int:
    if not rdkit_helpers.rdkit_available():
        return 0
    count = 0
    for smiles in smiles_list or []:
        mol = rdkit_helpers.parse_smiles(smiles)
        if mol is None:
            continue
        try:
            count += int(mol.GetRingInfo().NumRings())
        except Exception:
            continue
    return count


def _has_carboxylic_acid_like(motif_ids: Set[str]) -> bool:
    for mid in motif_ids:
        m = str(mid)
        if "COOH" in m:
            return True
    return False


def _has_explicit_acyl_activation_like(motif_ids: Set[str]) -> bool:
    markers = ("COCl", "COBr", "COF", "Anhydride", "MixedAnhydride", "Acyl-Imidazole")
    for mid in motif_ids:
        m = str(mid)
        if any(marker in m for marker in markers):
            return True
    return False


def _event(kind: str, confidence: float, details: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"kind": kind, "confidence": max(0.0, min(1.0, float(confidence)))}
    if details:
        payload["details"] = details
    return payload


def _quality_bucket(score: float) -> str:
    if score >= 0.75:
        return "high"
    if score >= 0.45:
        return "medium"
    return "low"


def _compute_redox_assessment(
    formed_pairs: Set[Tuple[str, str]],
    broken_pairs: Set[Tuple[str, str]],
) -> Dict[str, Any]:
    """
    Heuristic redox assessment from bond-pair deltas.

    Score convention:
    - positive: net oxidation signal on substrate carbons
    - negative: net reduction signal on substrate carbons
    """

    signals: List[str] = []
    score = 0.0

    formed_het = sum(
        1 for a, b in formed_pairs if "C" in (a, b) and ({a, b} - {"C"}) and next(iter({a, b} - {"C"})) in _REDOX_HETERO
    )
    broken_het = sum(
        1 for a, b in broken_pairs if "C" in (a, b) and ({a, b} - {"C"}) and next(iter({a, b} - {"C"})) in _REDOX_HETERO
    )
    formed_ch = int(("C", "H") in formed_pairs)
    broken_ch = int(("C", "H") in broken_pairs)

    # Forming C-hetero bonds and breaking C-H bonds are oxidation-like.
    score += float(formed_het)
    score += float(broken_ch)
    # Breaking C-hetero bonds and forming C-H bonds are reduction-like.
    score -= float(broken_het)
    score -= float(formed_ch)

    if formed_het:
        signals.append(f"formed_c_hetero={formed_het}")
    if broken_het:
        signals.append(f"broken_c_hetero={broken_het}")
    if formed_ch:
        signals.append("formed_c_h")
    if broken_ch:
        signals.append("broken_c_h")

    has_evidence = bool(formed_pairs or broken_pairs)
    if not has_evidence:
        return {
            "classification": "uncertain",
            "confidence": 0.2,
            "score": 0.0,
            "reasons": ["missing_bond_change_evidence"],
            "signals": [],
        }

    if abs(score) < 0.5:
        classification = "redox_neutral"
    elif score > 0:
        classification = "net_oxidation"
    else:
        classification = "net_reduction"

    # Confidence is higher when there are explicit C-H changes or multiple
    # oxidation-state-relevant bond changes.
    informative = formed_het + broken_het + formed_ch + broken_ch
    confidence = 0.55 + min(0.35, 0.1 * max(0, informative - 1))
    confidence = round(max(0.2, min(0.95, confidence)), 2)

    reasons: List[str] = []
    if classification == "redox_neutral":
        reasons.append("oxidation_and_reduction_signals_balanced_or_absent")
    elif classification == "net_oxidation":
        reasons.append("oxidation_signals_exceed_reduction_signals")
    else:
        reasons.append("reduction_signals_exceed_oxidation_signals")

    return {
        "classification": classification,
        "confidence": confidence,
        "score": round(score, 2),
        "reasons": reasons,
        "signals": signals,
    }


def _guess_electrophile_hybridization(reacted_motifs: Set[str]) -> str:
    for motif in reacted_motifs:
        text = str(motif)
        if text.startswith(_SCAFFOLD_SP_PREFIXES):
            return "sp"
    for motif in reacted_motifs:
        text = str(motif)
        if text.startswith(_SCAFFOLD_SP2_PREFIXES):
            if text.startswith(("Ar-", "HeteroAr-", "AromN-")):
                return "sp2_aromatic"
            return "sp2"
    for motif in reacted_motifs:
        if str(motif).startswith(_SCAFFOLD_SP3_PREFIXES):
            return "sp3"
    return "unknown"


def _electrophile_environment_tags(reacted_motifs: Set[str]) -> List[str]:
    tags: Set[str] = set()
    for motif in reacted_motifs:
        text = str(motif)
        if text.startswith("Bn-"):
            tags.add("benzylic")
        if text.startswith("Allyl-"):
            tags.add("allylic")
        if text.startswith("Propargyl-"):
            tags.add("propargylic")
        if text.startswith(("Ar-", "HeteroAr-", "AromN-")):
            tags.add("aromatic_sp2")
        if any(token in text for token in _EWG_TOKENS):
            tags.add("ewg_activated")
    return sorted(tags)


def _nucleophile_candidate_classes(
    reacted_motifs: Set[str],
    formed_pairs: Set[Tuple[str, str]],
) -> List[str]:
    classes: Set[str] = set()
    for motif in reacted_motifs:
        text = str(motif)
        for cls, tokens in _NUCLEOPHILE_CLASS_HINTS.items():
            if any(token in text for token in tokens):
                classes.add(cls)
    if ("C", "N") in formed_pairs:
        classes.add("amine")
    if ("C", "O") in formed_pairs:
        classes.add("alcohol_or_phenol")
    if ("C", "S") in formed_pairs:
        classes.add("thiol_or_thioether")
    if ("C", "C") in formed_pairs:
        classes.add("carbon_nucleophile_like")
    return sorted(classes)


def _snar_activation_supported(
    *,
    reacted_motifs: Set[str],
    formed_motifs: Set[str],
    electrophile_tags: List[str],
) -> Tuple[bool, str]:
    motif_space: Set[str] = set(reacted_motifs) | set(formed_motifs)
    if not motif_space:
        return False, "missing_motif_context"

    has_heteroaryl_context = any(
        str(motif).startswith(("HeteroAr-", "AromN-"))
        or any(tag in str(motif) for tag in ("Pyridine", "Pyridyl", "Pyrimidine", "Pyrazine", "Triazine", "Quinoline"))
        for motif in motif_space
    )
    if has_heteroaryl_context:
        return True, "heteroaryl_activation"

    if "ewg_activated" in set(electrophile_tags):
        return True, "ewg_activation"

    has_activation_token = any(
        any(token in str(motif) for token in _SNAR_ACTIVATING_TOKENS)
        for motif in motif_space
    )
    if has_activation_token:
        return True, "ewg_or_heteroaryl_activation"

    has_deactivating = any(
        any(token in str(motif) for token in _SNAR_DEACTIVATING_TOKENS)
        for motif in motif_space
    )
    if has_deactivating:
        return False, "electron_rich_aryl_context"

    return False, "no_activation_marker"


def _mechanism_shortlist(
    *,
    reacted_motifs: Set[str],
    formed_motifs: Set[str],
    event_kinds: Set[str],
    formed_pairs: Set[Tuple[str, str]],
    electrophile_hybridization: str,
    electrophile_tags: List[str],
    leaving_groups: Set[str],
) -> List[Dict[str, Any]]:
    shortlist: List[Dict[str, Any]] = []

    def add(name: str, confidence: float, reasons: List[str]) -> None:
        shortlist.append(
            {
                "name": name,
                "confidence": round(max(0.0, min(1.0, confidence)), 2),
                "reasons": reasons,
            }
        )

    has_displacement = "leaving_group_displacement" in event_kinds
    has_c_hetero = bool({("C", "N"), ("C", "O"), ("C", "S")} & formed_pairs)
    has_c_c = ("C", "C") in formed_pairs
    snar_supported, snar_reason = _snar_activation_supported(
        reacted_motifs=reacted_motifs,
        formed_motifs=formed_motifs,
        electrophile_tags=electrophile_tags,
    )

    if has_displacement and electrophile_hybridization == "sp3" and has_c_hetero:
        add("SN2", 0.72, ["sp3_electrophile", "c_hetero_bond_formation", "leaving_group_displacement"])
    if has_displacement and electrophile_hybridization == "sp2_aromatic" and has_c_hetero:
        conf = 0.78 if snar_supported else 0.52
        reasons = ["aryl_or_heteroaryl_electrophile", "leaving_group_displacement", "c_hetero_bond_formation"]
        reasons.append(snar_reason)
        add("SNAr", conf, reasons)
    if "amidation_like" in event_kinds:
        add("acyl_substitution", 0.8, ["amidation_like_event", "c_n_bond_formation"])

    has_sp2_halide = False
    for motif in reacted_motifs:
        text = str(motif)
        if text.startswith(("Ar-", "HeteroAr-", "Alkenyl-")) and any(
            token in text for token in ("-Cl", "-Br", "-I", "-OTf", "-OTs", "-OMs")
        ):
            has_sp2_halide = True
            break
    has_oa_partner = any(any(token in str(motif) for token in _OA_PARTNER_HINTS) for motif in reacted_motifs)
    if has_sp2_halide and (has_c_c or has_c_hetero) and has_oa_partner:
        add(
            "oa_based_coupling",
            0.74,
            ["sp2_electrophile_with_leaving_group", "cross_coupling_partner_like", "new_c_c_or_c_hetero_bond"],
        )

    # E2 is treated as a competing mechanism flag, not primary path assignment.
    if has_displacement and electrophile_hybridization == "sp3" and bool(leaving_groups):
        add("E2_competitor", 0.45, ["sp3_leaving_group_substrate", "possible_beta_h_elimination_competition"])

    shortlist.sort(key=lambda row: float(row.get("confidence", 0.0)), reverse=True)
    return shortlist[:4]


def _selectivity_risks(
    *,
    reacted_motifs: Set[str],
    event_kinds: Set[str],
    electrophile_hybridization: str,
    nucleophile_classes: List[str],
) -> List[str]:
    risks: Set[str] = set()
    if "leaving_group_displacement" in event_kinds and electrophile_hybridization == "sp3":
        # Practical default: substitution on sp3 C-LG often competes with elimination.
        has_non_methyl_sp3 = any(
            str(motif).startswith(prefix) for motif in reacted_motifs for prefix in ("Alkyl-", "RCH2-", "R2CH-", "R3C-", "Bn-", "Allyl-", "Propargyl-")
        )
        if has_non_methyl_sp3:
            risks.add("substitution_vs_elimination")

    ambident = any(any(token in str(motif) for token in _AMBIDENT_TOKENS) for motif in reacted_motifs)
    if ambident:
        risks.add("ambident_site_selectivity")

    # Simple overreaction signal: multiple nucleophile classes in the same record.
    if len(set(nucleophile_classes)) >= 2:
        risks.add("overreaction_or_chemoselectivity")

    return sorted(risks)


def _count_molecules_matching_smarts(
    smiles_list: Iterable[str],
    smarts_list: Iterable[str],
) -> int:
    if not rdkit_helpers.rdkit_available():
        return 0
    try:
        from chemtools.util.smarts_cache import compile_smarts
    except Exception:
        return 0
    patterns = []
    for smarts in smarts_list:
        patt = compile_smarts(smarts, validate=False)
        if patt is not None:
            patterns.append(patt)
    if not patterns:
        return 0

    hits = 0
    for smiles in smiles_list or []:
        mol = rdkit_helpers.parse_smiles(smiles)
        if mol is None:
            continue
        try:
            if any(mol.HasSubstructMatch(p) for p in patterns):
                hits += 1
        except Exception:
            continue
    return hits


def _infer_event_families(events: Iterable[Dict[str, Any]]) -> Set[str]:
    kinds = {
        str(ev.get("kind")).strip()
        for ev in (events or [])
        if isinstance(ev, dict) and str(ev.get("kind") or "").strip()
    }
    families: Set[str] = set()
    has_benzyl_o_alkylation = "benzyl_o_alkylation_like" in kinds
    if has_benzyl_o_alkylation:
        families.add("o_alkylation")

    if "ester_hydrolysis_like" in kinds:
        families.add("hydrolysis")
    if "amidation_like" in kinds:
        families.add("acyl_transfer")
    if "ring_closure_or_annulation" in kinds:
        families.add("annulation")

    bond_formation_kinds = {
        "c_n_bond_formation",
        "c_o_bond_formation",
        "c_s_bond_formation",
        "c_c_bond_formation",
    }
    has_bond_formation = bool(kinds & bond_formation_kinds)
    has_displacement = "leaving_group_displacement" in kinds
    if not has_benzyl_o_alkylation:
        if has_displacement and has_bond_formation:
            families.add("substitution")
        elif has_displacement:
            families.add("displacement")
        elif has_bond_formation:
            families.add("bond_formation")
    return families


def format_multi_event_signature(events: Iterable[Dict[str, Any]]) -> str:
    kinds = {
        str(ev.get("kind")).strip()
        for ev in (events or [])
        if isinstance(ev, dict) and str(ev.get("kind") or "").strip()
    }
    ordered_codes: List[str] = []
    for kind in _EVENT_SIGNATURE_PRIORITY:
        if kind not in kinds:
            continue
        code = _EVENT_SIGNATURE_CODE.get(kind)
        if code and code not in ordered_codes:
            ordered_codes.append(code)
    if len(ordered_codes) < 2:
        return ""
    return "+".join(ordered_codes[:4])


def summarize_reaction_events(
    *,
    reaction_smiles: str,
    bond_key: Optional[str],
    fallback_bond_key: Optional[str],
    reacted_motifs: Iterable[str],
    formed_motifs: Iterable[str],
    mapping_warning: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Build deterministic reaction-event summary and key quality diagnostics.
    """
    active_bond_key = bond_key or fallback_bond_key
    formed_tokens = _extract_bond_section(active_bond_key, section="form")
    broken_tokens = _extract_bond_section(active_bond_key, section="break")
    formed_pairs = _bond_pairs(formed_tokens)
    broken_pairs = _bond_pairs(broken_tokens)

    reactant_smiles, product_smiles = _split_reaction_sides(reaction_smiles)
    reactant_ring_count = _ring_count(reactant_smiles)
    product_ring_count = _ring_count(product_smiles)
    ring_delta = product_ring_count - reactant_ring_count

    reacted_set = {str(m) for m in reacted_motifs or [] if str(m).strip()}
    formed_set = {str(m) for m in formed_motifs or [] if str(m).strip()}

    events: List[Dict[str, Any]] = []
    anomalies: List[str] = []

    if ("C", "N") in formed_pairs:
        events.append(_event("c_n_bond_formation", 0.9))
    if ("C", "O") in formed_pairs:
        events.append(_event("c_o_bond_formation", 0.9))
    if ("C", "S") in formed_pairs:
        events.append(_event("c_s_bond_formation", 0.9))
    if ("C", "C") in formed_pairs:
        events.append(_event("c_c_bond_formation", 0.85))

    for hal in _HALOGENS:
        if tuple(sorted(("C", hal))) in broken_pairs:
            nucleophile = None
            for cand in ("N", "O", "S", "C"):
                if tuple(sorted(("C", cand))) in formed_pairs:
                    nucleophile = cand
                    break
            events.append(
                _event(
                    "leaving_group_displacement",
                    0.92,
                    {"leaving_group": hal, "nucleophile_element": nucleophile},
                )
            )

    if ring_delta > 0:
        events.append(
            _event(
                "ring_closure_or_annulation",
                0.88,
                {
                    "reactant_ring_count": reactant_ring_count,
                    "product_ring_count": product_ring_count,
                    "ring_delta": ring_delta,
                },
            )
        )

    if len(reactant_smiles) == 1 and formed_pairs:
        events.append(_event("intramolecular_likely", 0.8))
    elif len(reactant_smiles) > 1 and formed_pairs:
        events.append(_event("intermolecular_or_multi_component", 0.7))

    if ("C", "N") in formed_pairs and _has_carboxylic_acid_like(reacted_set):
        events.append(_event("amidation_like", 0.78))
        if not _has_explicit_acyl_activation_like(reacted_set):
            anomalies.append("amidation_without_explicit_activation_marker")

    reactant_ester_count = _count_molecules_matching_smarts(
        reactant_smiles,
        [_PATTERN_SMARTS["ester"]],
    )
    product_ester_count = _count_molecules_matching_smarts(
        product_smiles,
        [_PATTERN_SMARTS["ester"]],
    )
    reactant_acid_like_count = _count_molecules_matching_smarts(
        reactant_smiles,
        [_PATTERN_SMARTS["carboxy_acid"], _PATTERN_SMARTS["carboxylate"]],
    )
    product_acid_like_count = _count_molecules_matching_smarts(
        product_smiles,
        [_PATTERN_SMARTS["carboxy_acid"], _PATTERN_SMARTS["carboxylate"]],
    )
    if (
        reactant_ester_count > product_ester_count
        and product_acid_like_count > reactant_acid_like_count
    ):
        events.append(
            _event(
                "ester_hydrolysis_like",
                0.82,
                {
                    "reactant_ester_count": reactant_ester_count,
                    "product_ester_count": product_ester_count,
                    "reactant_acid_like_count": reactant_acid_like_count,
                    "product_acid_like_count": product_acid_like_count,
                },
            )
        )

    reactant_phenol_count = _count_molecules_matching_smarts(
        reactant_smiles,
        [_PATTERN_SMARTS["phenol"]],
    )
    reactant_benzyl_halide_count = _count_molecules_matching_smarts(
        reactant_smiles,
        [_PATTERN_SMARTS["benzyl_halide"]],
    )
    product_benzyl_aryl_ether_count = _count_molecules_matching_smarts(
        product_smiles,
        [_PATTERN_SMARTS["benzyl_aryl_ether"]],
    )
    if (
        reactant_phenol_count > 0
        and reactant_benzyl_halide_count > 0
        and product_benzyl_aryl_ether_count > 0
    ):
        events.append(
            _event(
                "benzyl_o_alkylation_like",
                0.84,
                {
                    "reactant_phenol_count": reactant_phenol_count,
                    "reactant_benzyl_halide_count": reactant_benzyl_halide_count,
                    "product_benzyl_aryl_ether_count": product_benzyl_aryl_ether_count,
                },
            )
        )

    event_families = _infer_event_families(events)
    event_kinds = {
        str(ev.get("kind")).strip()
        for ev in events
        if isinstance(ev, dict) and str(ev.get("kind") or "").strip()
    }
    molecularity = "unknown"
    if "intramolecular_likely" in event_kinds:
        molecularity = "intramolecular"
    elif "intermolecular_or_multi_component" in event_kinds:
        molecularity = "intermolecular_or_multi_component"

    leaving_groups: Set[str] = set()
    nucleophile_elements: Set[str] = set()
    for ev in events:
        if not isinstance(ev, dict):
            continue
        if str(ev.get("kind") or "").strip() != "leaving_group_displacement":
            continue
        details = ev.get("details") or {}
        if not isinstance(details, dict):
            continue
        lg = str(details.get("leaving_group") or "").strip()
        if lg:
            leaving_groups.add(lg)
        nuc = str(details.get("nucleophile_element") or "").strip()
        if nuc:
            nucleophile_elements.add(nuc)

    formed_bond_classes = sorted("-".join(pair) for pair in formed_pairs)
    broken_bond_classes = sorted("-".join(pair) for pair in broken_pairs)
    electrophile_hybridization = _guess_electrophile_hybridization(reacted_set)
    electrophile_tags = _electrophile_environment_tags(reacted_set)
    nucleophile_classes = _nucleophile_candidate_classes(reacted_set, formed_pairs)
    ambident_possible = any(
        any(token in str(motif) for token in _AMBIDENT_TOKENS)
        for motif in reacted_set
    )
    mechanism_shortlist = _mechanism_shortlist(
        reacted_motifs=reacted_set,
        formed_motifs=formed_set,
        event_kinds=event_kinds,
        formed_pairs=formed_pairs,
        electrophile_hybridization=electrophile_hybridization,
        electrophile_tags=electrophile_tags,
        leaving_groups=leaving_groups,
    )
    selectivity_risks = _selectivity_risks(
        reacted_motifs=reacted_set,
        event_kinds=event_kinds,
        electrophile_hybridization=electrophile_hybridization,
        nucleophile_classes=nucleophile_classes,
    )

    has_complex_family = bool(event_families & {"hydrolysis", "acyl_transfer", "annulation"})
    if len(event_families) >= 3 or (len(event_families) >= 2 and has_complex_family):
        anomalies.append("multi_transform_record_possible")

    if isinstance(mapping_warning, dict) and mapping_warning:
        anomalies.append("mapping_unreliable_fallback_used")

    quality_score = 1.0
    quality_reasons: List[str] = []
    if not formed_tokens and not formed_set:
        quality_score -= 0.45
        quality_reasons.append("missing_formed_bond_and_product_motif_evidence")
    if not active_bond_key:
        quality_score -= 0.25
        quality_reasons.append("missing_bond_key")
    if anomalies:
        quality_score -= min(0.45, 0.15 * len(anomalies))
        quality_reasons.extend(anomalies)
    quality_score = max(0.0, min(1.0, quality_score))
    redox_assessment = _compute_redox_assessment(formed_pairs, broken_pairs)

    return {
        "events": events,
        "anomalies": anomalies,
        "event_families": sorted(event_families),
        "electrophile_profile": {
            "center_element": "C",
            "hybridization_guess": electrophile_hybridization,
            "leaving_group_elements": sorted(leaving_groups),
            "environment_tags": electrophile_tags,
            "ewg_activation_likely": "ewg_activated" in set(electrophile_tags),
            "snar_activation_supported": _snar_activation_supported(
                reacted_motifs=reacted_set,
                formed_motifs=formed_set,
                electrophile_tags=electrophile_tags,
            )[0],
        },
        "nucleophile_profile": {
            "bond_forming_elements": sorted(nucleophile_elements),
            "candidate_classes": nucleophile_classes,
            "ambident_possible": ambident_possible,
        },
        "mechanism_shortlist": mechanism_shortlist,
        "selectivity_risks": selectivity_risks,
        "transformation_profile": {
            "molecularity": molecularity,
            "event_count": len(events),
            "formed_bond_classes": formed_bond_classes,
            "broken_bond_classes": broken_bond_classes,
            "leaving_groups": sorted(leaving_groups),
            "nucleophile_elements": sorted(nucleophile_elements),
            "ring_delta": ring_delta,
        },
        "bond_pairs": {
            "formed": sorted(list(formed_pairs)),
            "broken": sorted(list(broken_pairs)),
        },
        "ring_change": {
            "reactants": reactant_ring_count,
            "products": product_ring_count,
            "delta": ring_delta,
        },
        "reaction_key_quality": {
            "score_0_1": round(quality_score, 3),
            "level": _quality_bucket(quality_score),
            "reasons": quality_reasons,
        },
        "redox_assessment": redox_assessment,
    }


__all__ = [
    "format_multi_event_signature",
    "summarize_reaction_events",
]

