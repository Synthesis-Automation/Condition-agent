"""Extract and compare non-executable reaction-core applicability context."""

from __future__ import annotations

from typing import Any, Iterable

from reactive_taxonomy import analyze_molecule

from .models import CenterReactivityContext, TemplateContext


def _member(value: Any, name: str, default: Any = None) -> Any:
    if isinstance(value, dict):
        return value.get(name, default)
    return getattr(value, name, default)


def _center_maps(observation: Any) -> tuple[tuple[int, str], ...]:
    values = []
    for edit in _member(observation, "edits", ()) or ():
        if _member(edit, "edit_type") != "formed":
            continue
        for atom_name in ("atom_1", "atom_2"):
            atom = _member(edit, atom_name)
            if atom is None:
                continue
            map_number = _member(atom, "atom_map_number")
            if map_number is not None:
                values.append((int(map_number), str(_member(atom, "element", ""))))
    return tuple(sorted(set(values)))


def _center_locations(observation: Any) -> tuple[tuple[int, int, str], ...]:
    values = []
    for edit in _member(observation, "edits", ()) or ():
        if _member(edit, "edit_type") != "formed":
            continue
        for atom_name in ("atom_1", "atom_2"):
            atom = _member(edit, atom_name)
            if atom is None:
                continue
            values.append(
                (
                    int(_member(atom, "component_index", -1)),
                    int(_member(atom, "atom_index", -1)),
                    str(_member(atom, "element", "")),
                )
            )
    return tuple(sorted(set(values)))


def _center_profiles(
    mapped_reaction_smiles: str,
    observation: Any,
) -> tuple[CenterReactivityContext, ...]:
    reactant_side = mapped_reaction_smiles.split(">>", maxsplit=1)[0]
    wanted = dict(_center_maps(observation))
    locations = {
        (component_index, atom_index): element
        for component_index, atom_index, element in _center_locations(observation)
    }
    profiles = []
    if not wanted and not locations:
        return ()
    from rdkit import Chem

    for component_index, component_smiles in enumerate(reactant_side.split(".")):
        molecule = Chem.MolFromSmiles(component_smiles)
        if molecule is None:
            continue
        indexes = {
            int(atom.GetAtomMapNum()): int(atom.GetIdx())
            for atom in molecule.GetAtoms()
            if int(atom.GetAtomMapNum()) in wanted
        }
        elements_by_index = {
            atom_index: wanted[map_number]
            for map_number, atom_index in indexes.items()
        }
        elements_by_index.update(
            {
                atom_index: element
                for (located_component, atom_index), element in locations.items()
                if located_component == component_index
            }
        )
        if not elements_by_index:
            continue
        analysis = analyze_molecule(component_smiles)
        environments = {
            int(environment.center_atom_index): environment
            for environment in analysis.reactive_site_environments
        }
        for atom_index, element in sorted(elements_by_index.items()):
            environment = environments.get(atom_index)
            if environment is None:
                continue
            profile = environment.reactivity_profile
            profiles.append(
                CenterReactivityContext(
                    role="carbon" if element == "C" else "heteroatom",
                    element=element,
                    context_kind=str(profile.context_kind),
                    accessibility_class=str(profile.steric.accessibility_class),
                    accessibility_score=float(profile.steric.accessibility_score),
                    activation_axis=str(profile.electronic.activation_axis),
                    activation_class=str(profile.electronic.activation_class),
                    activation_score=float(profile.electronic.activation_score),
                )
            )
    return tuple(
        sorted(
            profiles,
            key=lambda value: (value.role, value.element, value.context_kind),
        )
    )


def context_from_analysis(
    analysis: Any,
    mapped_reaction_smiles: str,
) -> TemplateContext:
    """Build context from core attachment profiles and site descriptors."""

    core = _member(analysis, "reaction_core")
    spectator_groups = _member(analysis, "spectator_groups", ()) or ()
    tokens = set()
    for remote in _member(core, "remote_subgraphs", ()) or ():
        for port in _member(remote, "attachment_ports", ()) or ():
            profile = _member(port, "substituent_profile")
            tokens.update(_member(profile, "feature_tokens", ()) or ())
    return TemplateContext(
        spectator_group_ids=tuple(
            sorted(
                {
                    str(_member(group, "group_id", ""))
                    for group in spectator_groups
                    if _member(group, "group_id")
                }
            )
        ),
        substituent_feature_tokens=tuple(sorted(str(token) for token in tokens)),
        center_profiles=_center_profiles(
            mapped_reaction_smiles,
            _member(analysis, "observation"),
        ),
    )


def _jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    left_set = set(left)
    right_set = set(right)
    if not left_set and not right_set:
        return 1.0
    return len(left_set.intersection(right_set)) / len(left_set.union(right_set))


def context_similarity(left: TemplateContext, right: TemplateContext) -> float:
    """Return an interpretable similarity across structural and site context."""

    spectator_score = _jaccard(
        left.spectator_group_ids,
        right.spectator_group_ids,
    )
    substituent_score = _jaccard(
        left.substituent_feature_tokens,
        right.substituent_feature_tokens,
    )
    profile_scores = []
    right_by_role = {profile.role: profile for profile in right.center_profiles}
    for profile in left.center_profiles:
        other = right_by_role.get(profile.role)
        if other is None:
            continue
        categorical = sum(
            (
                profile.context_kind == other.context_kind,
                profile.activation_axis == other.activation_axis,
                profile.activation_class == other.activation_class,
            )
        ) / 3.0
        numeric = (
            max(0.0, 1.0 - abs(profile.accessibility_score - other.accessibility_score))
            + max(0.0, 1.0 - abs(profile.activation_score - other.activation_score) / 2.0)
        ) / 2.0
        profile_scores.append((categorical + numeric) / 2.0)
    profile_score = sum(profile_scores) / len(profile_scores) if profile_scores else 0.5
    return 0.30 * spectator_score + 0.35 * substituent_score + 0.35 * profile_score


__all__ = ["context_from_analysis", "context_similarity"]
