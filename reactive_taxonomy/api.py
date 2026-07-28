"""Public compound-featurization API for the isolated reactive taxonomy."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Iterable, List, Optional, Set

from .chemistry.rdkit_utils import mol_to_canonical_smiles, parse_smiles, rdkit_available

from .labels import (
    available_styles,
    render_anion,
    render_edge,
    render_named_handle,
    render_unsaturated_bond,
    render_xh,
)
from .environments import build_site_environment
from .functional_groups import detect_functional_groups
from .models import ComponentAnalysis, CompoundAnalysis, ReactiveSite, SiteCandidate, SiteType
from .patterns import MatchIndex
from .resolution import resolve_candidates
from .sites import DETECTORS


def _component_molecules(mol: Any) -> List[Any]:
    from rdkit import Chem  # type: ignore
    return list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True))


def _site_id(component_index: int, site_number: int, candidate: SiteCandidate) -> str:
    """Build an identifier from the topology's chemically meaningful locus."""
    if candidate.topology == "bond" and candidate.bond_indices:
        return f"mol{component_index}:bond{candidate.bond_indices[0]}:site{site_number}"
    preferred_role = "anchor" if candidate.topology == "edge" else "center"
    locus = candidate.atom_roles.get(preferred_role)
    if not locus:
        locus = candidate.atom_roles.get("endpoint_a") or candidate.atom_indices
    return f"mol{component_index}:atom{locus[0]}:site{site_number}"


def featurize_molecule(
    smiles: str,
    *,
    site_types: Optional[Iterable[SiteType]] = None,
    include_context_features: bool = True,
    label_style: str = "unicode",
) -> CompoundAnalysis:
    """Detect reactive annotations and canonical v2 connectivity sites.

    Disconnected components and formal charges are preserved. Atom indices in
    each site are component-local; ``component_index`` disambiguates them.
    """
    text = str(smiles or "").strip()
    if not text:
        return CompoundAnalysis(text, None, False, error="EMPTY_SMILES")
    if not rdkit_available():
        return CompoundAnalysis(text, None, False, error="RDKIT_UNAVAILABLE")
    mol = parse_smiles(text)
    if mol is None:
        return CompoundAnalysis(text, None, False, error="INVALID_SMILES")
    if label_style not in available_styles():
        return CompoundAnalysis(text, mol_to_canonical_smiles(mol), False, error=f"UNKNOWN_LABEL_STYLE:{label_style}")

    selected: Set[str] = set(site_types or DETECTORS)
    unknown = selected - set(DETECTORS)
    if unknown:
        return CompoundAnalysis(text, mol_to_canonical_smiles(mol), False, error=f"UNKNOWN_SITE_TYPES:{','.join(sorted(unknown))}")

    components: List[ComponentAnalysis] = []
    all_sites: List[ReactiveSite] = []
    all_functional_groups = []
    all_site_environments = []
    all_connectivity_sites = []
    atom_offset = 0
    for component_index, component_mol in enumerate(_component_molecules(mol)):
        functional_groups = detect_functional_groups(component_mol, component_index)
        match_index = MatchIndex(component_mol)
        raw_sites: List[SiteCandidate] = []
        for site_type, detector in DETECTORS.items():
            if site_type in selected:
                raw_sites.extend(detector(component_mol, match_index))
        raw_sites = resolve_candidates(raw_sites)
        component_sites: List[ReactiveSite] = []
        for site_number, candidate in enumerate(raw_sites):
            details = dict(candidate.details)
            details["atom_roles"] = {role: list(indices) for role, indices in candidate.atom_roles.items()}
            for role, indices in candidate.atom_roles.items():
                field = f"{role}_atom_index" if len(indices) == 1 else f"{role}_atom_indices"
                details.setdefault(field, indices[0] if len(indices) == 1 else list(indices))
            details["context_records"] = [record.to_dict() for record in candidate.context_records]
            if candidate.render_kind == "edge":
                label = render_edge(style=label_style, **candidate.render_data)
            elif candidate.render_kind == "xh":
                label = render_xh(style=label_style, **candidate.render_data)
            elif candidate.render_kind == "anion":
                label = render_anion(style=label_style, **candidate.render_data)
            elif candidate.render_kind == "named_handle":
                label = render_named_handle(style=label_style, **candidate.render_data)
            elif candidate.render_kind == "unsaturated_bond":
                label = render_unsaturated_bond(
                    style=label_style, **candidate.render_data
                )
            else:
                label = candidate.canonical_signature
            site = ReactiveSite(
                site_id=_site_id(component_index, site_number, candidate),
                site_type=candidate.site_type, topology=candidate.topology, component_index=component_index,
                atom_indices=list(candidate.atom_indices), bond_indices=[i for i in candidate.bond_indices if i >= 0],
                canonical_signature=candidate.canonical_signature, chemist_label=label,
                availability=candidate.availability, details=details,
                context_features={"contexts": [record.to_dict() for record in candidate.context_records]} if include_context_features else {},
                warnings=list(candidate.warnings),
            )
            component_sites.append(site)
        site_environments = [
            build_site_environment(component_mol, site, functional_groups)
            for site in component_sites
        ]
        environments_by_site = {environment.site_id: environment for environment in site_environments}
        component_sites = [
            replace(
                site,
                context_features={
                    **site.context_features,
                    "environment": environments_by_site[site.site_id].to_dict(),
                },
            )
            for site in component_sites
        ]
        from .reaction_site_interfaces import normalize_detected_site

        connectivity_sites = [
            normalize_detected_site(site, component_mol)
            for site in component_sites
        ]
        canonical_component = mol_to_canonical_smiles(component_mol) or ""
        input_component = canonical_component
        components.append(
            ComponentAnalysis(
                component_index=component_index,
                input_smiles=input_component,
                canonical_smiles=canonical_component,
                atom_offset=atom_offset,
                sites=component_sites,
                functional_groups=functional_groups,
                site_environments=site_environments,
                connectivity_sites=connectivity_sites,
            )
        )
        all_sites.extend(component_sites)
        all_functional_groups.extend(functional_groups)
        all_site_environments.extend(site_environments)
        all_connectivity_sites.extend(connectivity_sites)
        atom_offset += component_mol.GetNumAtoms()

    return CompoundAnalysis(
        text, mol_to_canonical_smiles(mol), True, components, all_sites,
        functional_groups=all_functional_groups,
        site_environments=all_site_environments,
        connectivity_sites=all_connectivity_sites,
    )


def detect_sites(mol: Any, site_types: Optional[Iterable[SiteType]] = None, *, label_style: str = "unicode") -> List[ReactiveSite]:
    """Detect sites on an RDKit molecule through the public result contract."""
    smiles = mol_to_canonical_smiles(mol)
    if not smiles:
        return []
    return featurize_molecule(smiles, site_types=site_types, label_style=label_style).sites


__all__ = ["detect_sites", "featurize_molecule"]
