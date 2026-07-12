"""Pronucleophilic X-H site detection."""

from __future__ import annotations

from typing import Any, List

from ..context import classify_neighbor_contexts
from ..models import SiteCandidate
from ..patterns import MatchIndex


def _derived_family(center: str, h_count: int, contexts: List[str]) -> str:
    if center == "N_aromatic": return "aromatic_nh"
    if center == "Csp": return "terminal_alkyne"
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
    if center == "N_aromatic" or "SO2R" in contexts or any(value.startswith("C(O)") for value in contexts):
        return "deactivated"
    if center == "O" and contexts and contexts[0].startswith("C(O)"):
        return "deactivated"
    return "free"


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    sites: List[SiteCandidate] = []
    candidate_centers = match_index.role_atoms("pronucleophile_XH", "center")
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        h_count = int(atom.GetTotalNumHs(includeNeighbors=True))
        center = symbol
        supported = symbol in {"N", "O", "S"}
        if symbol == "C" and str(atom.GetHybridization()) == "SP" and h_count == 1:
            supported, center = True, "Csp"
        if atom.GetIdx() not in candidate_centers or not supported or h_count < 1 or atom.GetFormalCharge() != 0:
            continue
        context_records = classify_neighbor_contexts(mol, atom.GetIdx(), match_index=match_index)
        contexts = [record.token for record in context_records]
        if symbol == "N" and atom.GetIsAromatic():
            center = "N_aromatic"
            contexts = ["HeteroAr"]
            if context_records:
                merged_atoms = tuple(sorted({idx for record in context_records for idx in record.fragment_atom_indices}))
                context_records = [type(context_records[0])(
                    token="HeteroAr", attachment_atom_index=atom.GetIdx(), fragment_atom_indices=merged_atoms,
                    classification_method="aromatic_ring_system", subtype="aromatic_nh_ring",
                    features={"ring_neighbor_count": len(context_records)},
                )]
        signature = f"XH|{center}|H{h_count}|{','.join(contexts)}"
        family = _derived_family(center, h_count, contexts)
        sites.append(SiteCandidate(
            site_type="pronucleophile_XH", topology="atom",
            atom_roles={"center": (atom.GetIdx(),)}, atom_indices=(atom.GetIdx(),), bond_indices=(),
            canonical_signature=signature, render_kind="xh",
            render_data={"center": center, "h_count": h_count, "contexts": contexts},
            matched_patterns=tuple(item["id"] for item in match_index.patterns_for_atom("pronucleophile_XH", "center", atom.GetIdx())),
            details={"center_element": symbol, "center_token": center, "formal_charge": atom.GetFormalCharge(), "aromatic": atom.GetIsAromatic(), "h_count": h_count, "contexts": contexts, "derived_family": family},
            context_records=tuple(context_records), availability=_availability(center, contexts),
        ))
    return sites
