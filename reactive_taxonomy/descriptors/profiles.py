"""Context dispatch and complete site-reactivity profile construction."""

from __future__ import annotations

from typing import Any, Iterable, Tuple

from .activated_centers import (
    activated_context_kind,
    build_activated_center_context,
)
from .alkyl import build_alkyl_context
from .aromatic import aromatic_atom_role, build_aromatic_context
from .common import (
    class_from_upper_bins,
    electronic_class,
    reactive_center_profile,
    shortest_distance,
)
from .heteroatom import build_heteroatom_context
from .models import (
    DescriptorEvidence,
    ElectronicContribution,
    ElectronicProfile,
    OtherContextDescriptor,
    SiteReactivityProfile,
    StericProfile,
)
from .modifiers import build_reactivity_modifiers
from .registry import descriptor_definition_versions, descriptor_rules
from .unsaturated import build_alkenyl_context, build_alkynyl_context


def _excluded_atoms(site: Any, center: int) -> Tuple[int, ...]:
    roles = site.details.get("atom_roles") or {}
    values = []
    for role in ("handle", "leaving_or_activatable", "departing_a"):
        raw = roles.get(role) or ()
        if isinstance(raw, int):
            raw = (raw,)
        values.extend(int(value) for value in raw if int(value) != center)
    return tuple(sorted(set(values)))


def _functional_group_contributions(
    mol: Any,
    center: int,
    groups: Iterable[Any],
) -> Tuple[ElectronicContribution, ...]:
    rules = descriptor_rules()["electronic"]
    radius = int(rules["radius"])
    weights = rules["functional_group_tag_weights"]
    values = []
    for group in groups:
        distance = shortest_distance(mol, center, group.atom_indices)
        if distance is None or distance < 1 or distance > radius:
            continue
        for tag in sorted(group.tags):
            weight = float(weights.get(tag, 0.0))
            if weight == 0.0:
                continue
            contribution = max(-1.0, min(1.0, weight / distance))
            values.append(
                ElectronicContribution(
                    source_id=f"molecular_motif:{group.motif_id}:{tag}",
                    effect="withdrawing" if contribution > 0 else "donating",
                    pathway="inductive",
                    positional_relation=f"distance_{distance}",
                    contribution=round(contribution, 3),
                    atom_indices=tuple(sorted(int(value) for value in group.atom_indices)),
                )
            )
    return tuple(
        sorted(
            values,
            key=lambda item: (
                item.source_id,
                item.positional_relation,
                item.atom_indices,
            ),
        )
    )


def _formal_charge_contribution(atom: Any) -> Tuple[ElectronicContribution, ...]:
    charge = int(atom.GetFormalCharge())
    if charge == 0:
        return ()
    weight = float(
        descriptor_rules()["electronic"]["formal_charge_weight"]
    )
    value = max(-1.0, min(1.0, charge * weight))
    return (
        ElectronicContribution(
            source_id="formal_charge",
            effect="withdrawing" if value > 0 else "donating",
            pathway="charge",
            positional_relation="center",
            contribution=value,
            atom_indices=(atom.GetIdx(),),
        ),
    )


def _steric_score(context: Any, contributions: Tuple[Any, ...]) -> float:
    kind = context.context_kind
    if kind == "aromatic":
        return float(context.ortho_burden_score)
    if kind == "alkyl":
        base = {
            "methyl": 0.05,
            "primary": 0.22,
            "secondary": 0.52,
            "tertiary": 0.82,
        }[context.carbon_substitution]
        return min(1.0, base + 0.04 * context.beta_branch_count)
    if kind == "heteroatom":
        return min(
            1.0,
            sum(item.score for item in contributions) / 2.5,
        )
    if kind in {"alkenyl", "alkynyl"}:
        return min(
            1.0,
            sum(item.score for item in contributions) / 2.0,
        )
    if kind in {"acyl", "sulfonyl", "phosphoryl"}:
        return min(
            1.0,
            sum(item.score for item in contributions) / 2.25,
        )
    return min(1.0, sum(item.score for item in contributions) / 2.0)


def _context_metrics(context: Any) -> Tuple[Tuple[str, str], ...]:
    kind = context.context_kind
    values: dict[str, str] = {"context_kind": str(kind)}
    if kind == "aromatic":
        values.update(
            {
                "ring_family": context.ring_family,
                "ortho_occupancy": (
                    f"{context.ortho_occupancy_count}/{context.ortho_capacity}"
                ),
                "ortho_burden": context.ortho_burden_class,
            }
        )
    elif kind == "alkyl":
        values.update(
            {
                "carbon_substitution": context.carbon_substitution,
                "alpha_branched": str(context.alpha_branched).lower(),
                "beta_branch_count": str(context.beta_branch_count),
                "beta_hydrogen_count": str(context.beta_hydrogen_count),
            }
        )
    elif kind in {"alkenyl", "alkynyl"}:
        values["endpoint_substitution"] = "/".join(
            str(value) for value in context.endpoint_substitution
        )
    elif kind in {"acyl", "sulfonyl", "phosphoryl"}:
        values["center_class"] = context.center_class
    elif kind == "heteroatom":
        values.update(
            {
                "substitution_class": context.substitution_class,
                "lone_pair_class": context.lone_pair_class,
                "alpha_branched_group_count": str(
                    context.alpha_branched_group_count
                ),
            }
        )
    return tuple(sorted(values.items()))


def build_site_reactivity_profile(
    mol: Any,
    site: Any,
    groups: Iterable[Any],
    nearby_groups: Iterable[dict[str, Any]],
    *,
    center_atom_index: int,
) -> SiteReactivityProfile:
    """Build one deterministic context-aware profile from molecular evidence."""
    center = int(center_atom_index)
    atom = mol.GetAtomWithIdx(center)
    excluded = _excluded_atoms(site, center)
    lone_pair_class = None
    lone_pair_availability = None
    acidity_class = None
    flags = []

    if atom.GetIsAromatic():
        context, steric_contributions, intrinsic_electronic = (
            build_aromatic_context(mol, center)
        )
        if atom.GetAtomicNum() in {7, 8, 16}:
            lone_pair_class = aromatic_atom_role(atom)
            lone_pair_availability = (
                "low" if lone_pair_class == "pyrrole_like" else "medium"
            )
            acidity_class = (
                "moderately_acidic"
                if atom.GetTotalNumHs() > 0
                else "not_applicable"
            )
            lone_pair_score = 0.55 if lone_pair_availability == "low" else 0.0
            intrinsic_electronic = (
                *intrinsic_electronic,
                ElectronicContribution(
                    source_id=f"lone_pair:{lone_pair_class}",
                    effect="withdrawing" if lone_pair_score > 0 else "mixed",
                    pathway="resonance",
                    positional_relation="center",
                    contribution=lone_pair_score,
                    atom_indices=(center,),
                ),
            )
    else:
        activated_kind = activated_context_kind(atom)
        if activated_kind is not None:
            context, steric_contributions, intrinsic_electronic = (
                build_activated_center_context(
                    mol, center, excluded_atoms=excluded
                )
            )
        elif atom.GetAtomicNum() == 6 and str(atom.GetHybridization()) in {
            "SP3",
            "S",
        }:
            context, steric_contributions, intrinsic_electronic = (
                build_alkyl_context(
                    mol, center, excluded_atoms=excluded
                )
            )
        elif atom.GetAtomicNum() == 6 and str(atom.GetHybridization()) == "SP2":
            try:
                context, steric_contributions, intrinsic_electronic = (
                    build_alkenyl_context(
                        mol, center, excluded_atoms=excluded
                    )
                )
            except ValueError:
                context = OtherContextDescriptor(
                    context_kind="other",
                    center_element=atom.GetSymbol(),
                    reason="unsupported_sp2_center",
                )
                steric_contributions = ()
                intrinsic_electronic = ()
        elif atom.GetAtomicNum() == 6 and str(atom.GetHybridization()) == "SP":
            try:
                context, steric_contributions, intrinsic_electronic = (
                    build_alkynyl_context(
                        mol, center, excluded_atoms=excluded
                    )
                )
            except ValueError:
                context = OtherContextDescriptor(
                    context_kind="other",
                    center_element=atom.GetSymbol(),
                    reason="unsupported_sp_center",
                )
                steric_contributions = ()
                intrinsic_electronic = ()
        elif atom.GetAtomicNum() in {7, 8, 15, 16}:
            (
                context,
                steric_contributions,
                intrinsic_electronic,
                lone_pair_class,
                lone_pair_availability,
                acidity_class,
            ) = build_heteroatom_context(
                mol, center, site, excluded_atoms=excluded
            )
        else:
            context = OtherContextDescriptor(
                context_kind="other",
                center_element=atom.GetSymbol(),
                reason="unsupported_center_context",
            )
            steric_contributions = ()
            intrinsic_electronic = ()

    context_kind = context.context_kind
    rules = descriptor_rules()
    score = round(_steric_score(context, steric_contributions), 3)
    access = class_from_upper_bins(score, rules["steric"]["bins"])
    burden = class_from_upper_bins(
        score, rules["steric"]["ortho_burden_bins"]
    )
    electronic_contributions = tuple(
        sorted(
            (
                *intrinsic_electronic,
                *_formal_charge_contribution(atom),
                *_functional_group_contributions(mol, center, groups),
            ),
            key=lambda item: (
                item.source_id,
                item.positional_relation,
                item.atom_indices,
            ),
        )
    )
    electronic_score = round(
        max(
            -1.0,
            min(
                1.0,
                sum(item.contribution for item in electronic_contributions),
            ),
        ),
        3,
    )
    axis = str(rules["electronic"]["activation_axes"][context_kind])
    if lone_pair_availability and atom.GetAtomicNum() in {7, 8, 15, 16}:
        axis = "lone_pair_availability"
        electronic_name = lone_pair_availability
    else:
        electronic_name = electronic_class(
            electronic_score, rules["electronic"]["bins"]
        )
    center_profile = reactive_center_profile(
        mol,
        center,
        lone_pair_class=lone_pair_class,
        lone_pair_availability=lone_pair_availability,
        acidity_class=acidity_class,
    )
    beta_hydrogen_count = (
        context.beta_hydrogen_count if context_kind == "alkyl" else None
    )
    modifiers = build_reactivity_modifiers(
        mol,
        site,
        center,
        beta_hydrogen_count=beta_hydrogen_count,
        ring_sizes=center_profile.ring_sizes,
        nearby_groups=nearby_groups,
    )
    if context_kind == "alkyl":
        if context.benzylic:
            flags.append("benzylic")
        if context.allylic:
            flags.append("allylic")
        if context.propargylic:
            flags.append("propargylic")
    if any(
        modifier.modifier_type == "coordination" for modifier in modifiers
    ):
        flags.append("coordination_risk")
    return SiteReactivityProfile(
        hypothesis_id=str(site.hypothesis_id),
        center_atom_index=center,
        context_kind=context_kind,
        context=context,
        reactive_center=center_profile,
        steric=StericProfile(
            accessibility_class=access,  # type: ignore[arg-type]
            accessibility_score=score,
            approach_burden_class=burden,  # type: ignore[arg-type]
            branch_contributions=steric_contributions,
            context_metrics=_context_metrics(context),
            evidence=DescriptorEvidence(
                source="molecular_graph",
                method=f"{context_kind}_steric_graph_v1",
                confidence=1.0,
                contributing_atom_indices=tuple(
                    sorted(
                        {
                            center,
                            *(
                                item.origin_atom_index
                                for item in steric_contributions
                            ),
                        }
                    )
                ),
            ),
        ),
        electronic=ElectronicProfile(
            activation_axis=axis,
            activation_class=electronic_name,
            activation_score=electronic_score,
            contributions=electronic_contributions,
            evidence=DescriptorEvidence(
                source="molecular_graph_and_motifs",
                method=f"{context_kind}_electronic_contributions_v1",
                confidence=(
                    1.0 if electronic_contributions else 0.65
                ),
                contributing_atom_indices=tuple(
                    sorted(
                        {
                            center,
                            *(
                                index
                                for item in electronic_contributions
                                for index in item.atom_indices
                            ),
                        }
                    )
                ),
                warnings=(
                    ()
                    if electronic_contributions
                    else ("no_recognized_electronic_contributors",)
                ),
            ),
        ),
        modifiers=modifiers,
        flags=tuple(sorted(set(flags))),
        status="unresolved" if context_kind == "other" else "derived",
        definition_versions=descriptor_definition_versions(),
    )


__all__ = ["build_site_reactivity_profile"]
