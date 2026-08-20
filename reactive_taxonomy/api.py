"""Public molecular structure and interpretation APIs."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Iterable, List, Optional, Set, Tuple

from .chemistry.rdkit_utils import mol_to_canonical_smiles, parse_smiles, rdkit_available
from .environments import build_site_environment
from .molecular_motifs import detect_molecular_motifs
from .labels import (
    available_styles,
    render_anion,
    render_edge,
    render_named_handle,
    render_unsaturated_bond,
    render_xh,
)
from .models import (
    MolecularAtomObservation,
    MolecularBondObservation,
    MolecularComponentInterpretation,
    MolecularComponentStructure,
    MolecularInterpretation,
    MolecularStructureObservation,
    MoleculeAnalysis,
    ReactiveSiteCandidate,
    ReactiveSiteHypothesis,
    SiteType,
)
from .patterns import MatchIndex
from .resolution import resolve_candidates
from .sites import DETECTORS


def _component_molecules(molecule: Any) -> List[Any]:
    from rdkit import Chem  # type: ignore

    return list(Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True))


def _component_structure(
    molecule: Any,
    *,
    component_index: int,
    atom_offset: int,
) -> MolecularComponentStructure:
    from rdkit import Chem  # type: ignore

    canonical = mol_to_canonical_smiles(molecule) or ""
    ordered = Chem.MolToSmiles(
        molecule,
        canonical=False,
        isomericSmiles=True,
    )
    atoms = tuple(
        MolecularAtomObservation(
            atom_index=int(atom.GetIdx()),
            element=str(atom.GetSymbol()),
            isotope=int(atom.GetIsotope()),
            formal_charge=int(atom.GetFormalCharge()),
            aromatic=bool(atom.GetIsAromatic()),
            hybridization=str(atom.GetHybridization()).upper(),
            total_hydrogens=int(atom.GetTotalNumHs(includeNeighbors=True)),
            degree=int(atom.GetDegree()),
            in_ring=bool(atom.IsInRing()),
            atom_map_number=(int(atom.GetAtomMapNum()) or None),
        )
        for atom in molecule.GetAtoms()
    )
    bonds = tuple(
        MolecularBondObservation(
            bond_index=int(bond.GetIdx()),
            atom_1_index=int(bond.GetBeginAtomIdx()),
            atom_2_index=int(bond.GetEndAtomIdx()),
            order=str(bond.GetBondType()).upper(),
            aromatic=bool(bond.GetIsAromatic()),
            conjugated=bool(bond.GetIsConjugated()),
            in_ring=bool(bond.IsInRing()),
            stereo=str(bond.GetStereo()),
        )
        for bond in molecule.GetBonds()
    )
    return MolecularComponentStructure(
        component_index=component_index,
        input_smiles=str(ordered),
        canonical_smiles=canonical,
        atom_offset=atom_offset,
        atoms=atoms,
        bonds=bonds,
    )


def observe_molecular_structure(smiles: str) -> MolecularStructureObservation:
    """Parse molecular graph facts without applying annotation definitions."""
    text = str(smiles or "").strip()
    if not text:
        return MolecularStructureObservation(text, None, False, error="EMPTY_SMILES")
    if not rdkit_available():
        return MolecularStructureObservation(text, None, False, error="RDKIT_UNAVAILABLE")
    molecule = parse_smiles(text)
    if molecule is None:
        return MolecularStructureObservation(text, None, False, error="INVALID_SMILES")
    components = []
    atom_offset = 0
    for component_index, component_molecule in enumerate(_component_molecules(molecule)):
        structure = _component_structure(
            component_molecule,
            component_index=component_index,
            atom_offset=atom_offset,
        )
        components.append(structure)
        atom_offset += len(structure.atoms)
    return MolecularStructureObservation(
        input_smiles=text,
        canonical_smiles=mol_to_canonical_smiles(molecule),
        valid=True,
        components=tuple(components),
    )


def _hypothesis_id(
    component_index: int,
    site_number: int,
    candidate: ReactiveSiteCandidate,
) -> str:
    if candidate.topology == "bond" and candidate.bond_indices:
        locus = f"bond{candidate.bond_indices[0]}"
    else:
        preferred_role = "anchor" if candidate.topology == "edge" else "center"
        indices = (
            candidate.atom_roles.get(preferred_role)
            or candidate.atom_roles.get("endpoint_a")
            or candidate.atom_indices
        )
        locus = f"atom{indices[0]}"
    return f"mol{component_index}:{locus}:hypothesis{site_number}"


def _render_candidate(candidate: ReactiveSiteCandidate, style: str) -> str:
    if candidate.render_kind == "edge":
        return render_edge(style=style, **candidate.render_data)
    if candidate.render_kind == "xh":
        return render_xh(style=style, **candidate.render_data)
    if candidate.render_kind == "anion":
        return render_anion(style=style, **candidate.render_data)
    if candidate.render_kind == "named_handle":
        return render_named_handle(style=style, **candidate.render_data)
    if candidate.render_kind == "unsaturated_bond":
        return render_unsaturated_bond(style=style, **candidate.render_data)
    return candidate.canonical_signature


def interpret_molecular_reactivity(
    structure: MolecularStructureObservation,
    *,
    site_types: Optional[Iterable[SiteType]] = None,
    include_context_features: bool = True,
    label_style: str = "unicode",
) -> MolecularInterpretation:
    """Return optional motifs and reactive-site hypotheses for a structure."""
    if not structure.valid:
        return MolecularInterpretation()
    if label_style not in available_styles():
        raise ValueError(f"UNKNOWN_LABEL_STYLE:{label_style}")
    selected: Set[str] = set(site_types or DETECTORS)
    unknown = selected - set(DETECTORS)
    if unknown:
        raise ValueError(f"UNKNOWN_SITE_TYPES:{','.join(sorted(unknown))}")

    component_interpretations = []
    all_motifs = []
    all_hypotheses = []
    all_environments = []
    all_connectivity = []
    for component in structure.components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        motifs = detect_molecular_motifs(
            molecule, component.component_index, label_style=label_style
        )
        match_index = MatchIndex(molecule)
        raw: List[ReactiveSiteCandidate] = []
        for site_type, detector in DETECTORS.items():
            if site_type in selected:
                raw.extend(detector(molecule, match_index))
        candidates = resolve_candidates(raw)
        hypotheses = []
        for number, candidate in enumerate(candidates):
            details = dict(candidate.details)
            details["atom_roles"] = {
                role: list(indices) for role, indices in candidate.atom_roles.items()
            }
            for role, indices in candidate.atom_roles.items():
                field_name = (
                    f"{role}_atom_index" if len(indices) == 1 else f"{role}_atom_indices"
                )
                details.setdefault(
                    field_name,
                    indices[0] if len(indices) == 1 else list(indices),
                )
            details["context_records"] = [
                record.to_dict() for record in candidate.context_records
            ]
            hypotheses.append(
                ReactiveSiteHypothesis(
                    hypothesis_id=_hypothesis_id(component.component_index, number, candidate),
                    site_type=candidate.site_type,
                    topology=candidate.topology,
                    component_index=component.component_index,
                    atom_indices=tuple(candidate.atom_indices),
                    bond_indices=tuple(index for index in candidate.bond_indices if index >= 0),
                    canonical_signature=candidate.canonical_signature,
                    chemist_label=_render_candidate(candidate, label_style),
                    availability=candidate.availability,
                    details=details,
                    context_features=(
                        {
                            "contexts": [
                                record.to_dict() for record in candidate.context_records
                            ]
                        }
                        if include_context_features
                        else {}
                    ),
                    warnings=tuple(candidate.warnings),
                )
            )
        environments = tuple(
            build_site_environment(molecule, hypothesis, motifs)
            for hypothesis in hypotheses
        )
        by_hypothesis = {
            environment.hypothesis_id: environment for environment in environments
        }
        hypotheses = [
            replace(
                hypothesis,
                context_features={
                    **hypothesis.context_features,
                    "environment": by_hypothesis[hypothesis.hypothesis_id].to_dict(),
                },
            )
            for hypothesis in hypotheses
        ]
        from .reaction_site_interfaces import normalize_detected_site

        connectivity = tuple(
            normalize_detected_site(hypothesis, molecule)
            for hypothesis in hypotheses
        )
        component_interpretations.append(
            MolecularComponentInterpretation(
                component_index=component.component_index,
                motifs=tuple(motifs),
                reactive_site_hypotheses=tuple(hypotheses),
                reactive_site_environments=environments,
                connectivity_hypotheses=connectivity,
            )
        )
        all_motifs.extend(motifs)
        all_hypotheses.extend(hypotheses)
        all_environments.extend(environments)
        all_connectivity.extend(connectivity)
    return MolecularInterpretation(
        components=tuple(component_interpretations),
        motifs=tuple(all_motifs),
        reactive_site_hypotheses=tuple(all_hypotheses),
        reactive_site_environments=tuple(all_environments),
        connectivity_hypotheses=tuple(all_connectivity),
    )


def analyze_molecule(
    smiles: str,
    *,
    site_types: Optional[Iterable[SiteType]] = None,
    include_context_features: bool = True,
    label_style: str = "unicode",
) -> MoleculeAnalysis:
    """Compose structure observation and optional molecular interpretation."""
    structure = observe_molecular_structure(smiles)
    if label_style not in available_styles() and structure.valid:
        structure = replace(structure, valid=False, error=f"UNKNOWN_LABEL_STYLE:{label_style}")
    try:
        interpretation = interpret_molecular_reactivity(
            structure,
            site_types=site_types,
            include_context_features=include_context_features,
            label_style=label_style,
        )
    except ValueError as error:
        structure = replace(structure, valid=False, error=str(error))
        interpretation = MolecularInterpretation()
    return MoleculeAnalysis(structure=structure, interpretation=interpretation)


def detect_reactive_site_hypotheses(
    molecule: Any,
    site_types: Optional[Iterable[SiteType]] = None,
    *,
    label_style: str = "unicode",
) -> Tuple[ReactiveSiteHypothesis, ...]:
    """Return optional reactive-site hypotheses for an RDKit molecule."""
    smiles = mol_to_canonical_smiles(molecule)
    if not smiles:
        return ()
    return analyze_molecule(
        smiles,
        site_types=site_types,
        label_style=label_style,
    ).interpretation.reactive_site_hypotheses


__all__ = [
    "analyze_molecule",
    "detect_reactive_site_hypotheses",
    "interpret_molecular_reactivity",
    "observe_molecular_structure",
]
