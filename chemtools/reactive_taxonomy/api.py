"""Public compound-featurization API for the isolated reactive taxonomy."""

from __future__ import annotations

from typing import Any, Iterable, List, Optional, Set

from chemtools.core.rdkit import mol_to_canonical_smiles, parse_smiles, rdkit_available

from .labels import available_styles, render_edge, render_xh
from .models import ComponentAnalysis, CompoundAnalysis, ReactiveSite, SiteCandidate, SiteType
from .patterns import MatchIndex
from .resolution import resolve_candidates
from .sites import DETECTORS


def _component_molecules(mol: Any) -> List[Any]:
    from rdkit import Chem  # type: ignore
    return list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True))


def featurize_molecule(
    smiles: str,
    *,
    site_types: Optional[Iterable[SiteType]] = None,
    include_context_features: bool = True,
    label_style: str = "unicode",
) -> CompoundAnalysis:
    """Detect v1 reactive handles without using the legacy taxonomy.

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
    atom_offset = 0
    for component_index, component_mol in enumerate(_component_molecules(mol)):
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
            else:
                label = candidate.canonical_signature
            site = ReactiveSite(
                site_id=f"mol{component_index}:atom{candidate.atom_indices[0]}:site{site_number}",
                site_type=candidate.site_type, topology=candidate.topology, component_index=component_index,
                atom_indices=list(candidate.atom_indices), bond_indices=[i for i in candidate.bond_indices if i >= 0],
                canonical_signature=candidate.canonical_signature, chemist_label=label,
                availability=candidate.availability, details=details,
                context_features={"contexts": [record.to_dict() for record in candidate.context_records]} if include_context_features else {},
                warnings=list(candidate.warnings),
            )
            component_sites.append(site)
        canonical_component = mol_to_canonical_smiles(component_mol) or ""
        input_component = canonical_component
        components.append(ComponentAnalysis(component_index, input_component, canonical_component, atom_offset, component_sites))
        all_sites.extend(component_sites)
        atom_offset += component_mol.GetNumAtoms()

    return CompoundAnalysis(text, mol_to_canonical_smiles(mol), True, components, all_sites)


def detect_sites(mol: Any, site_types: Optional[Iterable[SiteType]] = None, *, label_style: str = "unicode") -> List[ReactiveSite]:
    """Detect sites on an RDKit molecule through the public result contract."""
    smiles = mol_to_canonical_smiles(mol)
    if not smiles:
        return []
    return featurize_molecule(smiles, site_types=site_types, label_style=label_style).sites


__all__ = ["detect_sites", "featurize_molecule"]
