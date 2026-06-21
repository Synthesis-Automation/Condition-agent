"""
General taxonomy-first reaction inference with explicit evidence separation.

The analyzer is intentionally deterministic:
- Certain evidence: parsed reaction parts, motif delta, bond/event tokens, formula delta.
- Hypotheses: ranked taxonomy candidates and mechanism-family label.

This keeps the workflow general across reaction families while still supporting
taxonomy-gap signaling when no reaction type can be matched confidently.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
import re
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from chemtools.core.smiles import normalize_reaction
from chemtools.core.rdkit import parse_smiles, rdkit_available
from chemtools.reaction.feasibility import validate_detection_with_crk_key
from chemtools.reaction.typing import detect_reaction_type
from chemtools.taxonomy.reaction_catalog import get_reaction_type, motif_tokens_compatible


@dataclass(frozen=True)
class ReactionDecision:
    reaction_type: str
    confidence: float
    source: str
    mechanism_family: str
    rationale: str


@dataclass(frozen=True)
class ReactionValidation:
    passed: bool
    checks: List[Dict[str, Any]] = field(default_factory=list)
    tautomer_or_representation_issue: bool = False


@dataclass(frozen=True)
class GeneralReactionAnalysis:
    reaction_smiles: str
    reactant_smiles_list: List[str]
    product_smiles_list: List[str]
    principal_pair: Dict[str, Any]
    formula_delta: Dict[str, int]
    features: Dict[str, Any]
    evidence: Dict[str, Any]
    taxonomy_candidates: List[Dict[str, Any]]
    decision: ReactionDecision
    validation: ReactionValidation
    summary: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "reaction_smiles": self.reaction_smiles,
            "reactant_smiles_list": list(self.reactant_smiles_list),
            "product_smiles_list": list(self.product_smiles_list),
            "principal_pair": dict(self.principal_pair),
            "formula_delta": dict(self.formula_delta),
            "features": dict(self.features),
            "evidence": dict(self.evidence),
            "taxonomy_candidates": [dict(row) for row in self.taxonomy_candidates],
            "decision": {
                "reaction_type": self.decision.reaction_type,
                "confidence": self.decision.confidence,
                "source": self.decision.source,
                "mechanism_family": self.decision.mechanism_family,
                "rationale": self.decision.rationale,
            },
            "validation": {
                "passed": self.validation.passed,
                "checks": [dict(row) for row in self.validation.checks],
                "tautomer_or_representation_issue": self.validation.tautomer_or_representation_issue,
            },
            "summary": self.summary,
        }


def _preferred_side_smiles(side_payload: Iterable[Mapping[str, Any]]) -> List[str]:
    smiles_list: List[str] = []
    for entry in side_payload:
        smiles = (
            str(entry.get("smiles_norm") or "").strip()
            or str(entry.get("largest_smiles") or "").strip()
            or str(entry.get("input") or "").strip()
        )
        if smiles:
            smiles_list.append(smiles)
    return smiles_list


def _count_elements(smiles: str) -> Counter[str]:
    mol = parse_smiles(smiles)
    if mol is None:
        return Counter()
    counts: Counter[str] = Counter()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol == "H":
            continue
        counts[symbol] += 1
    return counts


def _formula_delta(reactant_smiles: str, product_smiles: str) -> Dict[str, int]:
    reactant_counts = _count_elements(reactant_smiles)
    product_counts = _count_elements(product_smiles)
    elements = sorted(set(reactant_counts) | set(product_counts))
    delta: Dict[str, int] = {}
    for element in elements:
        change = int(product_counts.get(element, 0) - reactant_counts.get(element, 0))
        if change != 0:
            delta[element] = change
    return delta


def _formula_delta_for_sets(
    reactants: Sequence[str],
    products: Sequence[str],
) -> Dict[str, int]:
    reactant_counts: Counter[str] = Counter()
    product_counts: Counter[str] = Counter()
    for smiles in reactants:
        reactant_counts.update(_count_elements(smiles))
    for smiles in products:
        product_counts.update(_count_elements(smiles))
    elements = sorted(set(reactant_counts) | set(product_counts))
    delta: Dict[str, int] = {}
    for element in elements:
        change = int(product_counts.get(element, 0) - reactant_counts.get(element, 0))
        if change != 0:
            delta[element] = change
    return delta


def _count_aromatic_ring_n(smiles: str) -> int:
    mol = parse_smiles(smiles)
    if mol is None:
        return 0
    return sum(
        1
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 7 and atom.GetIsAromatic()
    )


def _has_n_n_bond(smiles: str) -> bool:
    mol = parse_smiles(smiles)
    if mol is None:
        text = smiles.lower()
        return "nn" in text or "n=n" in text or "n-n" in text
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom().GetAtomicNum()
        b = bond.GetEndAtom().GetAtomicNum()
        if a == 7 and b == 7:
            return True
    return False


def _mcs_ratio(reactant_smiles: str, product_smiles: str) -> float:
    reactant = parse_smiles(reactant_smiles)
    product = parse_smiles(product_smiles)
    if reactant is None or product is None:
        return 0.0
    try:
        from rdkit.Chem import rdFMCS  # type: ignore
    except Exception:
        return 0.0

    reactant_heavy = max(1, reactant.GetNumHeavyAtoms())
    product_heavy = max(1, product.GetNumHeavyAtoms())
    denominator = max(reactant_heavy, product_heavy)
    if denominator <= 0:
        return 0.0

    try:
        result = rdFMCS.FindMCS(
            [reactant, product],
            timeout=5,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            matchValences=False,
        )
    except Exception:
        return 0.0

    if not getattr(result, "canceled", False):
        atoms = int(getattr(result, "numAtoms", 0) or 0)
        if atoms > 0:
            return round(min(1.0, atoms / float(denominator)), 4)
    return 0.0


def _atom_overlap_ratio(reactant_smiles: str, product_smiles: str) -> float:
    a = _count_elements(reactant_smiles)
    b = _count_elements(product_smiles)
    if not a or not b:
        return 0.0
    common = sum(min(a.get(k, 0), b.get(k, 0)) for k in set(a) | set(b))
    total = max(1, sum(b.values()))
    return round(common / float(total), 4)


def _select_principal_pair(
    reactants: Sequence[str],
    products: Sequence[str],
) -> Dict[str, Any]:
    if not reactants or not products:
        return {
            "reactant_index": 0,
            "product_index": 0,
            "reactant_smiles": reactants[0] if reactants else "",
            "product_smiles": products[0] if products else "",
            "mcs_ratio": 0.0,
            "overlap_ratio": 0.0,
            "score": 0.0,
        }

    rows: List[Dict[str, Any]] = []
    for reactant_index, reactant_smiles in enumerate(reactants):
        for product_index, product_smiles in enumerate(products):
            mcs = _mcs_ratio(reactant_smiles, product_smiles) if rdkit_available() else 0.0
            overlap = _atom_overlap_ratio(reactant_smiles, product_smiles)
            score = round(0.8 * mcs + 0.2 * overlap, 4)
            rows.append(
                {
                    "reactant_index": int(reactant_index),
                    "product_index": int(product_index),
                    "reactant_smiles": reactant_smiles,
                    "product_smiles": product_smiles,
                    "mcs_ratio": mcs,
                    "overlap_ratio": overlap,
                    "score": score,
                }
            )

    rows.sort(
        key=lambda row: (
            float(row.get("score", 0.0)),
            float(row.get("mcs_ratio", 0.0)),
            float(row.get("overlap_ratio", 0.0)),
            -int(row.get("reactant_index", 0)),
            -int(row.get("product_index", 0)),
        ),
        reverse=True,
    )
    return rows[0] if rows else {}


def _canonicalize_bond_token(token: str) -> Optional[str]:
    text = str(token or "").strip()
    if "-" not in text:
        return None
    left, right = [part.strip() for part in text.split("-", 1)]
    left_match = re.search(r"[A-Z][a-z]?", left)
    right_match = re.search(r"[A-Z][a-z]?", right)
    if not left_match or not right_match:
        return None
    left_el = left_match.group(0)
    right_el = right_match.group(0)
    return "-".join(sorted((left_el, right_el)))


def _parse_reaction_key(reaction_key: str) -> Dict[str, Any]:
    sections = [section.strip() for section in str(reaction_key or "").split(" | ") if section.strip()]
    reacted_motifs: List[str] = []
    formed_motifs: List[str] = []
    formed_bonds: List[str] = []
    broken_bonds: List[str] = []
    event_tokens: List[str] = []

    if sections:
        summary = sections[0]
        if summary.startswith("CRK-v1"):
            summary = summary[len("CRK-v1") :].strip()
        if summary.startswith("|"):
            summary = summary[1:].strip()
        if "->" in summary:
            lhs, rhs = summary.split("->", 1)
            reacted_motifs = [tok.strip() for tok in lhs.split("|") if tok.strip() and tok.strip() != "[]"]
            formed_motifs = [tok.strip() for tok in rhs.split("|") if tok.strip() and tok.strip() != "[]"]

    for section in sections[1:]:
        lower = section.lower()
        if lower.startswith("bond_formed:"):
            payload = section.split(":", 1)[1]
            tokens = [tok.strip() for tok in payload.split(";") if tok.strip()]
            formed_bonds.extend(tokens)
        elif lower.startswith("bond_broken:"):
            payload = section.split(":", 1)[1]
            tokens = [tok.strip() for tok in payload.split(";") if tok.strip()]
            broken_bonds.extend(tokens)
        elif lower.startswith("events:"):
            payload = section.split(":", 1)[1]
            for chunk in [tok.strip() for tok in payload.split(";") if tok.strip()]:
                for token in [part.strip() for part in chunk.split("+") if part.strip()]:
                    event_tokens.append(token)

    canonical_formed = sorted(
        {
            token
            for token in (_canonicalize_bond_token(value) for value in formed_bonds)
            if token
        }
    )
    canonical_broken = sorted(
        {
            token
            for token in (_canonicalize_bond_token(value) for value in broken_bonds)
            if token
        }
    )
    return {
        "reacted_motifs": reacted_motifs,
        "formed_motifs": formed_motifs,
        "formed_bonds": canonical_formed,
        "broken_bonds": canonical_broken,
        "event_tokens": event_tokens,
    }


def _score_detection_candidate(match: Any) -> float:
    slot_hits = int(bool(getattr(match, "electrophile", []))) + int(
        bool(getattr(match, "nucleophile", []))
    ) + int(bool(getattr(match, "product", [])))
    base = 0.8 * float(getattr(match, "confidence", 0.0))
    slot_bonus = 0.06 * float(slot_hits)
    source_bonus = 0.04 if bool(getattr(match, "slot_sources", {})) else 0.0
    return round(min(1.0, base + slot_bonus + source_bonus), 4)


def _score_crk_candidate(candidate: Mapping[str, Any], rank: int) -> float:
    matched_slots = len(candidate.get("matched_reactant_slots") or [])
    matched_products = int(candidate.get("matched_product_slots") or 0)
    reactant_support = int(candidate.get("reactant_support") or 0)
    product_support = int(candidate.get("product_support") or 0)
    score = (
        0.45
        + 0.08 * matched_slots
        + 0.06 * matched_products
        + 0.02 * reactant_support
        + 0.02 * product_support
        - 0.02 * max(0, rank)
    )
    return round(min(0.92, max(0.0, score)), 4)


def _merge_candidates(rows: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    merged: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        rid = str(row.get("reaction_type") or "").strip()
        if not rid:
            continue
        previous = merged.get(rid)
        if previous is None:
            merged[rid] = dict(row)
            continue
        if float(row.get("deterministic_score", 0.0)) > float(
            previous.get("deterministic_score", 0.0)
        ):
            merged[rid] = dict(row)
            continue
        evidence = list(previous.get("evidence_sources") or [])
        for item in row.get("evidence_sources") or []:
            if item not in evidence:
                evidence.append(item)
        previous["evidence_sources"] = evidence
        previous["confidence"] = max(
            float(previous.get("confidence", 0.0)),
            float(row.get("confidence", 0.0)),
        )

    ordered = list(merged.values())
    ordered.sort(
        key=lambda row: (
            float(row.get("deterministic_score", 0.0)),
            float(row.get("confidence", 0.0)),
            -len(str(row.get("reaction_type", ""))),
        ),
        reverse=True,
    )
    return ordered


def _as_set(values: Any) -> set[str]:
    if not values:
        return set()
    if isinstance(values, str):
        return {values}
    if isinstance(values, (list, tuple, set)):
        return {str(value) for value in values if str(value).strip()}
    return set()


def _constraints_hold(
    constraints: Mapping[str, Any],
    reacted_set: set[str],
    formed_set: set[str],
    formed_bonds: set[str],
    broken_bonds: set[str],
) -> bool:
    include_reacted = _as_set(constraints.get("include_reacted"))
    exclude_reacted = _as_set(constraints.get("exclude_reacted"))
    include_formed = _as_set(constraints.get("include_formed"))
    exclude_formed = _as_set(constraints.get("exclude_formed"))
    include_bond_formed = _as_set(constraints.get("include_bond_formed"))
    exclude_bond_formed = _as_set(constraints.get("exclude_bond_formed"))
    include_bond_broken = _as_set(constraints.get("include_bond_broken"))
    exclude_bond_broken = _as_set(constraints.get("exclude_bond_broken"))
    if include_reacted and not include_reacted.issubset(reacted_set):
        return False
    if exclude_reacted and (exclude_reacted & reacted_set):
        return False
    if include_formed and not include_formed.issubset(formed_set):
        return False
    if exclude_formed and (exclude_formed & formed_set):
        return False
    if include_bond_formed and formed_bonds and not include_bond_formed.issubset(formed_bonds):
        return False
    if exclude_bond_formed and (exclude_bond_formed & formed_bonds):
        return False
    if include_bond_broken and broken_bonds and not include_bond_broken.issubset(broken_bonds):
        return False
    if exclude_bond_broken and (exclude_bond_broken & broken_bonds):
        return False
    return True


def _taxonomy_consistency_check(
    reaction_type: str,
    reacted_motifs: Sequence[str],
    formed_motifs: Sequence[str],
    formed_bonds: Sequence[str],
    broken_bonds: Sequence[str],
) -> Tuple[bool, str]:
    if not reaction_type or reaction_type == "unknown":
        return False, "unknown_reaction_type"
    definition = get_reaction_type(reaction_type)
    if definition is None:
        return False, "reaction_type_missing_from_taxonomy"

    reacted_set = {str(token) for token in reacted_motifs if str(token).strip()}
    formed_set = {str(token) for token in formed_motifs if str(token).strip()}
    formed_bond_set = {str(token) for token in formed_bonds if str(token).strip()}
    broken_bond_set = {str(token) for token in broken_bonds if str(token).strip()}

    def _slot_matches(allowed_tokens: Sequence[str], observed: set[str]) -> bool:
        for observed_token in observed:
            for allowed_token in (allowed_tokens or []):
                if motif_tokens_compatible(str(observed_token), str(allowed_token)):
                    return True
        return False

    for slot_name, slot_req in definition.reactants.items():
        allowed = [str(token) for token in (slot_req.allowed or []) if str(token).strip()]
        if allowed and not _slot_matches(allowed, reacted_set):
            return False, f"missing_reactant_slot:{slot_name}"

    if definition.products:
        product_match = False
        for slot_name, slot_req in definition.products.items():
            allowed = [str(token) for token in (slot_req.allowed or []) if str(token).strip()]
            if allowed and _slot_matches(allowed, formed_set):
                product_match = True
                break
            if not allowed:
                product_match = True
                break
        if not product_match:
            return False, "missing_product_slot"

    if not _constraints_hold(
        definition.constraints or {},
        reacted_set,
        formed_set,
        formed_bond_set,
        broken_bond_set,
    ):
        return False, "constraint_mismatch"
    return True, "ok"


def _infer_mechanism_family(
    reaction_type: str,
    *,
    formed_bonds: Sequence[str],
    broken_bonds: Sequence[str],
    event_tokens: Sequence[str],
) -> str:
    rid = str(reaction_type or "").strip().lower()
    if not rid or rid == "unknown":
        events = {token.lower() for token in event_tokens}
        formed = set(formed_bonds)
        broken = set(broken_bonds)
        if "ann" in events or "cycl" in events:
            return "annulation_cyclization"
        if "amid" in events:
            return "acyl_transfer"
        if "oxid" in events:
            return "oxidation"
        if "red" in events:
            return "reduction"
        halogen_loss = bool(broken & {"C-Cl", "C-Br", "C-F", "C-I"})
        if halogen_loss and "C-N" in formed:
            return "nucleophilic_substitution"
        if "lgdisp" in events and "C-C" in formed:
            return "substitution"
        if "C-N" in formed:
            return "bond_formation_c_n"
        return "unknown"

    if rid.startswith("oxidation"):
        return "oxidation"
    if rid.startswith("reduction"):
        return "reduction"
    if "snar" in rid:
        return "SNAr"
    if "amide" in rid or "acyl" in rid:
        return "acyl_transfer"
    if "annulation" in rid or "cyclization" in rid:
        return "annulation_cyclization"
    if "coupling" in rid or rid in {
        "suzuki_miyaura",
        "sonogashira",
        "negishi",
        "stille",
        "hiyama",
        "kumada",
        "heck",
    }:
        return "cross_coupling"
    if "halogenation" in rid:
        return "halogenation"
    if "deprotection" in rid:
        return "deprotection"
    return "other"


def _carbon_bond_partners(bond_tokens: Sequence[str]) -> List[str]:
    partners: set[str] = set()
    for token in bond_tokens or []:
        text = str(token or "").strip()
        if "-" not in text:
            continue
        left, right = [part.strip() for part in text.split("-", 1)]
        if left == "C" and right and right != "H":
            partners.add(right)
        elif right == "C" and left and left != "H":
            partners.add(left)
    priority = {
        "N": 0,
        "O": 1,
        "S": 2,
        "C": 3,
        "B": 4,
        "Si": 5,
        "Sn": 6,
        "P": 7,
    }
    return sorted(partners, key=lambda el: (priority.get(el, 99), el))


def _format_general_bond_change_label(partners: Sequence[str], *, action: str) -> Optional[str]:
    normalized = [str(el).strip() for el in partners if str(el).strip()]
    if not normalized:
        return None
    pair_labels = [f"C_{el}" if el != "C" else "C_C" for el in normalized]
    joined = " / ".join(pair_labels)
    if action == "formation":
        return f"General {joined} bond formation reaction"
    if action == "cleavage":
        return f"General {joined} cleavage reaction"
    return None


def _general_bond_change_fallback_reaction_type(
    *,
    formed_bonds: Sequence[str],
    broken_bonds: Sequence[str],
) -> Optional[str]:
    formed_partners = _carbon_bond_partners(formed_bonds)
    if formed_partners:
        return _format_general_bond_change_label(formed_partners, action="formation")
    broken_partners = _carbon_bond_partners(broken_bonds)
    if broken_partners:
        return _format_general_bond_change_label(broken_partners, action="cleavage")
    return None


def _detect_tautomer_issue(
    products: Sequence[str],
    principal_product: str,
    formula_delta: Mapping[str, int],
) -> bool:
    if not products:
        return False
    if "[nH]" in principal_product:
        return True
    if formula_delta.get("H", 0) == 0 and formula_delta.get("N", 0) != 0:
        lower = principal_product.lower()
        if "n" in lower and ("=n" in lower or "[n" in lower):
            return True
    return False


def analyze_reaction_general(
    reaction_smiles: str,
    *,
    use_llm: bool = False,
    min_confidence: float = 0.55,
) -> GeneralReactionAnalysis:
    """
    Run general reaction inference from reaction SMILES.

    The function is deterministic and taxonomy-driven. `use_llm` is accepted for
    compatibility but not used in this implementation.
    """
    _ = use_llm
    normalized = normalize_reaction(reaction_smiles)
    reactants = _preferred_side_smiles(normalized.get("reactants") or [])
    products = _preferred_side_smiles(normalized.get("products") or [])
    principal_pair = _select_principal_pair(reactants, products)
    core_reactant = str(principal_pair.get("reactant_smiles") or "")
    core_product = str(principal_pair.get("product_smiles") or "")
    formula_delta = _formula_delta(core_reactant, core_product) if core_reactant and core_product else {}
    net_formula_delta = _formula_delta_for_sets(reactants, products)

    side_reactants = [
        smiles
        for idx, smiles in enumerate(reactants)
        if idx != int(principal_pair.get("reactant_index", -1))
    ]
    features: Dict[str, Any] = {
        "reactant_count": len(reactants),
        "product_count": len(products),
        "principal_pair_mcs_ratio": float(principal_pair.get("mcs_ratio", 0.0)),
        "core_reactant_aromatic_ring_n_count": _count_aromatic_ring_n(core_reactant),
        "core_product_aromatic_ring_n_count": _count_aromatic_ring_n(core_product),
        "side_reactants_have_nn": any(_has_n_n_bond(smiles) for smiles in side_reactants),
        "formula_delta_halogen_loss": sum(
            int(formula_delta.get(el, 0)) for el in ("F", "Cl", "Br", "I")
        ),
        "formula_delta_n": int(formula_delta.get("N", 0)),
        "net_formula_delta_nonzero_count": len(net_formula_delta),
    }

    detection = detect_reaction_type(reaction_smiles)
    key_parts = _parse_reaction_key(detection.reaction_key)
    reacted_motifs = sorted({str(token) for token in detection.reacted_motifs if token})
    formed_motifs = sorted({str(token) for token in detection.formed_motifs if token})

    candidate_rows: List[Dict[str, Any]] = []
    for match in detection.matches:
        candidate_rows.append(
            {
                "reaction_type": str(match.reaction_type),
                "confidence": float(match.confidence),
                "deterministic_score": _score_detection_candidate(match),
                "source": "taxonomy_detection",
                "evidence_sources": ["reaction_key", "motif_delta", "slot_matching"],
            }
        )

    crk_probe = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=detection.reaction_key,
        include_evidence=True,
    )
    crk_candidates = (
        ((crk_probe.get("evidence") or {}).get("top_candidates") or [])
        if isinstance(crk_probe, dict)
        else []
    )
    for rank, row in enumerate(crk_candidates):
        reaction_id = str(row.get("reaction_id") or "").strip()
        if not reaction_id:
            continue
        score = _score_crk_candidate(row, rank)
        candidate_rows.append(
            {
                "reaction_type": reaction_id,
                "confidence": score,
                "deterministic_score": score,
                "source": "taxonomy_crk_specificity",
                "evidence_sources": ["reaction_key", "motif_delta", "constraints"],
            }
        )

    taxonomy_candidates = _merge_candidates(candidate_rows)
    top_candidate = taxonomy_candidates[0] if taxonomy_candidates else None
    decision_source = "deterministic"
    if top_candidate and float(top_candidate.get("deterministic_score", 0.0)) >= float(min_confidence):
        decision_reaction_type = str(top_candidate.get("reaction_type") or "unknown")
        decision_confidence = float(top_candidate.get("deterministic_score", 0.0))
        decision_rationale = "top taxonomy candidate passed deterministic confidence threshold"
    else:
        fallback_reaction_type = _general_bond_change_fallback_reaction_type(
            formed_bonds=key_parts.get("formed_bonds") or [],
            broken_bonds=key_parts.get("broken_bonds") or [],
        )
        if fallback_reaction_type:
            decision_reaction_type = fallback_reaction_type
            decision_confidence = 0.35
            decision_source = "general_bond_change_fallback"
            decision_rationale = "no taxonomy candidate passed threshold; assigned bond-change fallback label"
        else:
            decision_reaction_type = "unknown"
            decision_confidence = 0.0
            decision_rationale = "no taxonomy candidate cleared deterministic confidence threshold"

    mechanism_family = _infer_mechanism_family(
        decision_reaction_type if decision_source == "deterministic" else "unknown",
        formed_bonds=key_parts.get("formed_bonds") or [],
        broken_bonds=key_parts.get("broken_bonds") or [],
        event_tokens=key_parts.get("event_tokens") or [],
    )
    if mechanism_family == "unknown" and decision_reaction_type == "unknown":
        halogen_loss = any(int(formula_delta.get(el, 0)) < 0 for el in ("F", "Cl", "Br", "I"))
        nitrogen_gain = int(formula_delta.get("N", 0)) > 0
        aromatic_n_rich = int(features.get("core_reactant_aromatic_ring_n_count") or 0) >= 2
        if halogen_loss and nitrogen_gain and aromatic_n_rich:
            mechanism_family = "nucleophilic_substitution"
    tautomer_issue = _detect_tautomer_issue(products, core_product, formula_delta)
    if decision_source == "general_bond_change_fallback":
        consistency_ok, consistency_reason = True, "not_applicable_general_bond_change_fallback"
    else:
        consistency_ok, consistency_reason = _taxonomy_consistency_check(
            decision_reaction_type,
            reacted_motifs,
            formed_motifs,
            key_parts.get("formed_bonds") or [],
            key_parts.get("broken_bonds") or [],
        )
    checks = [
        {
            "check": "parsed_reaction_components",
            "passed": bool(reactants and products),
            "detail": f"reactants={len(reactants)}, products={len(products)}",
        },
        {
            "check": "motif_delta_available",
            "passed": bool(reacted_motifs or formed_motifs),
            "detail": f"reacted={len(reacted_motifs)}, formed={len(formed_motifs)}",
        },
        {
            "check": "taxonomy_consistency",
            "passed": bool(consistency_ok),
            "detail": consistency_reason,
        },
    ]
    validation_passed = bool(decision_reaction_type != "unknown" and consistency_ok)
    validation = ReactionValidation(
        passed=validation_passed,
        checks=checks,
        tautomer_or_representation_issue=tautomer_issue,
    )

    decision = ReactionDecision(
        reaction_type=decision_reaction_type,
        confidence=round(decision_confidence, 4),
        source=decision_source,
        mechanism_family=mechanism_family,
        rationale=decision_rationale,
    )
    summary_bits = [
        f"Decision: {decision.reaction_type} ({decision.confidence:.2f}) via {decision.source}.",
        f"Mechanism family: {decision.mechanism_family}.",
    ]
    if tautomer_issue:
        summary_bits.append("Tautomer/representation ambiguity was detected in the product-side core.")
    if decision.source == "general_bond_change_fallback":
        summary_bits.append("Applied general bond-change fallback classification because no taxonomy candidate passed threshold.")
    if decision.reaction_type == "unknown":
        summary_bits.append("No taxonomy-consistent candidate passed threshold; treat as potential taxonomy gap.")

    evidence = {
        "certain": {
            "principal_pair": dict(principal_pair),
            "formula_delta": dict(formula_delta),
            "net_formula_delta": dict(net_formula_delta),
            "reaction_key": detection.reaction_key,
            "reacted_motifs": reacted_motifs,
            "formed_motifs": formed_motifs,
            "formed_bonds": list(key_parts.get("formed_bonds") or []),
            "broken_bonds": list(key_parts.get("broken_bonds") or []),
            "event_tokens": list(key_parts.get("event_tokens") or []),
        },
        "hypotheses": [dict(row) for row in taxonomy_candidates],
        "raw_detection": {
            "match_count": len(detection.matches),
            "error": detection.error,
        },
    }

    return GeneralReactionAnalysis(
        reaction_smiles=reaction_smiles,
        reactant_smiles_list=list(reactants),
        product_smiles_list=list(products),
        principal_pair=dict(principal_pair),
        formula_delta=dict(formula_delta),
        features=features,
        evidence=evidence,
        taxonomy_candidates=taxonomy_candidates,
        decision=decision,
        validation=validation,
        summary=" ".join(summary_bits),
    )


def classify_reaction(
    reaction_smiles: str,
    *,
    use_llm: bool = False,
    min_confidence: float = 0.55,
) -> ReactionDecision:
    """
    Canonical reaction classification entry point.

    This wraps ``analyze_reaction_general`` and returns only the final decision so
    application code can use one unified classification path without duplicating
    fallback logic.
    """
    analysis = analyze_reaction_general(
        reaction_smiles,
        use_llm=use_llm,
        min_confidence=min_confidence,
    )
    return analysis.decision


__all__ = [
    "ReactionDecision",
    "ReactionValidation",
    "GeneralReactionAnalysis",
    "analyze_reaction_general",
    "classify_reaction",
]
