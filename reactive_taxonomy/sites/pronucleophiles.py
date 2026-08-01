"""Pronucleophilic X-H site detection."""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Tuple

from ..context import classify_neighbor_contexts
from ..models import ContextClassification, ReactiveSiteCandidate
from ..patterns import MatchIndex


def _derived_family(center: str, h_count: int, contexts: List[str]) -> str:
    if center == "N_aromatic": return "aromatic_nh"
    if center == "Csp": return "terminal_alkyne"
    if center == "Csp3": return "activated_sp3_carbon"
    if center == "O": return "phenol" if contexts == ["Ar"] else "alcohol"
    if center == "S": return "thiophenol" if contexts == ["Ar"] else "thiol"
    if "SO2R" in contexts: return "sulfonamide"
    if any(value.startswith("C(O)") for value in contexts): return "amide_like"
    if "N" in contexts: return "hydrazine"
    if center == "N":
        if h_count == 3 and not contexts: return "ammonia"
        return "primary_amine" if h_count == 2 else "secondary_amine"
    return "pronucleophile"


def _availability(center: str, contexts: List[str]) -> str:
    if center == "Csp3":
        return "activated"
    if center == "N_aromatic" or "SO2R" in contexts or any(value.startswith("C(O)") for value in contexts):
        return "deactivated"
    if center == "O" and contexts and contexts[0].startswith("C(O)"):
        return "deactivated"
    return "free"


def _activation_contexts(
    matches: List[Tuple[Dict[str, Any], Dict[str, Tuple[int, ...]]]],
) -> tuple[
    List[str],
    List[ContextClassification],
    List[Dict[str, Any]],
    Tuple[int, ...],
]:
    """Normalize matched acidic-carbon rules into deterministic context evidence."""
    retained: Dict[Tuple[str, int], Tuple[Dict[str, Any], Dict[str, Tuple[int, ...]]]] = {}
    for definition, roles in matches:
        token = str(definition.get("activation_token") or "")
        anchors = roles.get("activation_anchor") or ()
        if not token or len(anchors) != 1:
            continue
        retained[(token, int(anchors[0]))] = (definition, roles)

    activation_records: List[Dict[str, Any]] = []
    context_records: List[ContextClassification] = []
    for (token, anchor), (definition, roles) in sorted(retained.items()):
        relationship = str(
            definition.get("activation_relationship") or "alpha_to"
        )
        fragment = tuple(
            sorted(set(roles.get("activation_fragment") or (anchor,)))
        )
        pattern_id = str(definition["id"])
        activation_records.append(
            {
                "relationship": relationship,
                "activator": token,
                "activator_atom_index": anchor,
                "fragment_atom_indices": list(fragment),
                "rule_id": pattern_id,
            }
        )
        context_records.append(
            ContextClassification(
                token=f"{relationship}:{token}",
                attachment_atom_index=anchor,
                fragment_atom_indices=fragment,
                classification_method="mapped_smarts",
                facet="activation",
                semantic_id=f"{relationship}_{token}",
                display_token=f"{relationship}:{token}",
                subtype=token,
                matched_pattern=pattern_id,
                features={
                    "relationship": relationship,
                    "activator": token,
                },
            )
        )

    counts = Counter(
        f"{record['relationship']}:{record['activator']}"
        for record in activation_records
    )
    contexts = [
        token if count == 1 else f"{token}*{count}"
        for token, count in sorted(counts.items())
    ]
    anchors = tuple(
        sorted(
            {
                int(record["activator_atom_index"])
                for record in activation_records
            }
        )
    )
    return contexts, context_records, activation_records, anchors


def detect(mol: Any, match_index: MatchIndex) -> List[ReactiveSiteCandidate]:
    sites: List[ReactiveSiteCandidate] = []
    candidate_centers = match_index.role_atoms("pronucleophile_XH", "center")
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        h_count = int(atom.GetTotalNumHs(includeNeighbors=True))
        center = symbol
        supported = symbol in {"N", "O", "S"}
        pattern_matches = match_index.role_matches_for_atom(
            "pronucleophile_XH",
            "center",
            atom.GetIdx(),
        )
        activation_matches = [
            item
            for item in pattern_matches
            if item[0].get("activation_token")
        ]
        if symbol == "C" and str(atom.GetHybridization()) == "SP" and h_count == 1:
            supported, center = True, "Csp"
        elif (
            symbol == "C"
            and str(atom.GetHybridization()) == "SP3"
            and activation_matches
        ):
            supported, center = True, "Csp3"
        if atom.GetIdx() not in candidate_centers or not supported or h_count < 1 or atom.GetFormalCharge() != 0:
            continue
        activation_records: List[Dict[str, Any]] = []
        activation_anchors: Tuple[int, ...] = ()
        if center == "Csp3":
            (
                contexts,
                context_records,
                activation_records,
                activation_anchors,
            ) = _activation_contexts(activation_matches)
            if not contexts:
                continue
        else:
            context_records = classify_neighbor_contexts(
                mol,
                atom.GetIdx(),
                match_index=match_index,
            )
            contexts = [record.token for record in context_records]
        if symbol == "N" and atom.GetIsAromatic():
            center = "N_aromatic"
            contexts = ["HeteroAr"]
            if context_records:
                merged_atoms = tuple(sorted({idx for record in context_records for idx in record.fragment_atom_indices}))
                context_records = [type(context_records[0])(
                    token="HeteroAr", attachment_atom_index=atom.GetIdx(), fragment_atom_indices=merged_atoms,
                    classification_method="aromatic_ring_system",
                    facet="scaffold",
                    semantic_id="heteroaromatic",
                    display_token="HeteroAr",
                    subtype="aromatic_nh_ring",
                    features={"ring_neighbor_count": len(context_records)},
                )]
        signature = f"XH|{center}|H{h_count}|{','.join(contexts)}"
        family = _derived_family(center, h_count, contexts)
        render_features = {}
        if center == "Csp":
            opposite = next(
                (
                    neighbor
                    for neighbor in atom.GetNeighbors()
                    if neighbor.GetAtomicNum() == 6
                    and str(neighbor.GetHybridization()) == "SP"
                ),
                None,
            )
            if opposite is not None:
                render_features["opposite_endpoint_substituent_count"] = sum(
                    neighbor.GetAtomicNum() > 1
                    and neighbor.GetIdx() != atom.GetIdx()
                    for neighbor in opposite.GetNeighbors()
                )
        elif center == "Csp3":
            render_features["activation_count"] = len(activation_records)
            render_features["activation_tokens"] = [
                str(record["activator"]) for record in activation_records
            ]
        atom_roles = {"center": (atom.GetIdx(),)}
        heavy_neighbors = tuple(
            sorted(
                neighbor.GetIdx()
                for neighbor in atom.GetNeighbors()
                if neighbor.GetAtomicNum() > 1
            )
        )
        if len(heavy_neighbors) == 1:
            atom_roles["attachment"] = heavy_neighbors
        attachment = (
            mol.GetAtomWithIdx(int(heavy_neighbors[0]))
            if len(heavy_neighbors) == 1
            else None
        )
        if activation_anchors:
            atom_roles["activation_anchor"] = activation_anchors
        sites.append(ReactiveSiteCandidate(
            site_type="pronucleophile_XH", topology="atom",
            atom_roles=atom_roles, atom_indices=(atom.GetIdx(),), bond_indices=(),
            canonical_signature=signature, render_kind="xh",
            render_data={
                "center": center,
                "h_count": h_count,
                "contexts": contexts,
                "features": render_features,
            },
            matched_patterns=tuple(
                sorted(
                    {
                        str(definition["id"])
                        for definition, _ in pattern_matches
                    }
                )
            ),
            details={
                "center_element": symbol,
                "center_token": center,
                "formal_charge": atom.GetFormalCharge(),
                "aromatic": atom.GetIsAromatic(),
                "h_count": h_count,
                "contexts": contexts,
                "derived_family": family,
                "attachment_h_count": (
                    int(attachment.GetTotalNumHs(includeNeighbors=True))
                    if attachment is not None
                    else None
                ),
                "attachment_hybridization": (
                    str(attachment.GetHybridization()).upper()
                    if attachment is not None
                    else None
                ),
                "attachment_context": (
                    str(contexts[0]) if attachment is not None and contexts else None
                ),
                "activation_relationship": (
                    "alpha_to" if activation_records else None
                ),
                "activation_records": activation_records,
                **render_features,
            },
            context_records=tuple(context_records), availability=_availability(center, contexts),
        ))
    return sites
