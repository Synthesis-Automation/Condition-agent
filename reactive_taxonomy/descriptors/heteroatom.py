"""Descriptors for N, O, S, and P reactive centers."""

from __future__ import annotations

from typing import Any, Iterable, Tuple

from .common import branch_contribution, heteroatom_substitution
from .models import (
    AttachedGroupProfile,
    ElectronicContribution,
    HeteroatomContextDescriptor,
    StericContribution,
)


def _site_contexts(site: Any) -> Tuple[str, ...]:
    details = site.details or {}
    values = details.get("contexts") or ()
    if isinstance(values, str):
        values = (values,)
    return tuple(sorted(str(value) for value in values if value))


def _attached_context(mol: Any, center: int, neighbor: Any) -> str:
    if neighbor.GetIsAromatic():
        return "Ar"
    if neighbor.GetAtomicNum() == 6:
        if any(
            bond.GetOtherAtom(neighbor).GetAtomicNum() == 8
            and str(bond.GetBondType()).upper() == "DOUBLE"
            for bond in neighbor.GetBonds()
        ):
            return "C(O)R"
        if str(neighbor.GetHybridization()) in {"SP3", "S"}:
            return "Alkyl"
        if str(neighbor.GetHybridization()) == "SP2":
            return "Alkenyl"
        if str(neighbor.GetHybridization()) == "SP":
            return "Alkynyl"
    if neighbor.GetAtomicNum() == 16 and sum(
        adjacent.GetAtomicNum() == 8
        and str(mol.GetBondBetweenAtoms(
            neighbor.GetIdx(), adjacent.GetIdx()
        ).GetBondType()).upper()
        == "DOUBLE"
        for adjacent in neighbor.GetNeighbors()
    ) >= 2:
        return "SO2R"
    return neighbor.GetSymbol()


def heteroatom_lone_pair(
    mol: Any, center: int, contexts: Iterable[str]
) -> Tuple[str, str, str | None]:
    """Return lone-pair class, availability, and acidity class."""
    atom = mol.GetAtomWithIdx(center)
    values = set(contexts)
    if atom.GetFormalCharge() > 0:
        return "cationic", "low", "not_applicable"
    if atom.GetFormalCharge() < 0:
        return "anionic", "high", "conjugate_base"
    if "SO2R" in values:
        return "sulfonamide_like", "low", "acidic"
    if values & {"C(O)R", "C(O)NR", "C(O)OR"}:
        return "amide_like", "low", "moderately_acidic"
    if values & {"Ar", "HetAr"}:
        return "aryl_delocalized", "medium", "weakly_acidic"
    if atom.GetSymbol() == "N":
        return "amine_like", "high", "weakly_acidic"
    if atom.GetSymbol() in {"O", "S"} and atom.GetTotalNumHs() > 0:
        return (
            "localized",
            "medium" if atom.GetSymbol() == "O" else "high",
            "moderately_acidic",
        )
    return "localized", "medium", None


def build_heteroatom_context(
    mol: Any,
    center: int,
    site: Any,
    *,
    excluded_atoms: Iterable[int] = (),
) -> Tuple[
    HeteroatomContextDescriptor,
    Tuple[StericContribution, ...],
    Tuple[ElectronicContribution, ...],
    str,
    str,
    str | None,
]:
    """Build attachment, resonance, accessibility, and lone-pair evidence."""
    atom = mol.GetAtomWithIdx(center)
    excluded = {int(value) for value in excluded_atoms}
    graph_contexts = tuple(
        _attached_context(mol, center, neighbor)
        for neighbor in atom.GetNeighbors()
        if neighbor.GetAtomicNum() > 1 and neighbor.GetIdx() not in excluded
    )
    contexts = tuple(sorted(set(_site_contexts(site)) | set(graph_contexts)))
    attached_groups = []
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() <= 1 or neighbor.GetIdx() in excluded:
            continue
        context_token = _attached_context(mol, center, neighbor)
        attachment_class = None
        alpha_branched = False
        beta_branch_count = 0
        if (
            neighbor.GetAtomicNum() == 6
            and not neighbor.GetIsAromatic()
            and str(neighbor.GetHybridization()) in {"SP3", "S"}
        ):
            carbon_neighbors = sum(
                adjacent.GetAtomicNum() == 6
                and adjacent.GetIdx() != center
                for adjacent in neighbor.GetNeighbors()
            )
            attachment_class = {
                0: "methyl",
                1: "primary",
                2: "secondary",
                3: "tertiary",
            }.get(carbon_neighbors, "quaternary_or_complex")
            alpha_branched = carbon_neighbors >= 2
            beta_branch_count = sum(
                max(
                    0,
                    sum(
                        branch.GetAtomicNum() == 6
                        for branch in adjacent.GetNeighbors()
                    )
                    - 1,
                )
                for adjacent in neighbor.GetNeighbors()
                if adjacent.GetAtomicNum() == 6
                and adjacent.GetIdx() != center
            )
        attached_groups.append(
            AttachedGroupProfile(
                atom_index=neighbor.GetIdx(),
                context=context_token,
                element=neighbor.GetSymbol(),
                attachment_carbon_class=attachment_class,
                alpha_branched=alpha_branched,
                beta_branch_count=beta_branch_count,
            )
        )
    lone_pair_class, availability, acidity = heteroatom_lone_pair(
        mol, center, contexts
    )
    resonance = (
        "sulfonyl_delocalized"
        if "SO2R" in contexts
        else "carbonyl_delocalized"
        if set(contexts) & {"C(O)R", "C(O)NR", "C(O)OR"}
        else "aryl_delocalized"
        if set(contexts) & {"Ar", "HetAr"}
        else "localized"
    )
    contributions = tuple(
        sorted(
            (
                branch_contribution(
                    mol,
                    neighbor.GetIdx(),
                    blocked=(center, *excluded),
                    relation="attached_group",
                    radius=3,
                )
                for neighbor in atom.GetNeighbors()
                if neighbor.GetAtomicNum() > 1
                and neighbor.GetIdx() not in excluded
            ),
            key=lambda item: item.origin_atom_index,
        )
    )
    score = {
        "high": -0.45,
        "medium": 0.0,
        "low": 0.55,
    }[availability]
    electronic = (
        ElectronicContribution(
            source_id=f"lone_pair:{lone_pair_class}",
            effect="withdrawing" if score > 0 else "donating" if score < 0 else "mixed",
            pathway="resonance" if "delocalized" in resonance else "inductive",
            positional_relation="center",
            contribution=score,
            atom_indices=(center,),
        ),
    )
    context = HeteroatomContextDescriptor(
        context_kind="heteroatom",
        element=atom.GetSymbol(),
        substitution_class=heteroatom_substitution(atom),
        attached_contexts=contexts,
        resonance_class=resonance,
        lone_pair_class=lone_pair_class,
        proton_count=int(atom.GetTotalNumHs()),
        alpha_branched_group_count=sum(
            item.alpha_branched for item in attached_groups
        ),
        attached_groups=tuple(
            sorted(attached_groups, key=lambda item: item.atom_index)
        ),
    )
    return (
        context,
        contributions,
        electronic,
        lone_pair_class,
        availability,
        acidity,
    )


__all__ = ["build_heteroatom_context", "heteroatom_lone_pair"]
