"""Public compound-featurization API for the isolated reactive taxonomy."""

from __future__ import annotations

from typing import Any, Iterable, List, Optional, Set

from chemtools.core.rdkit import mol_to_canonical_smiles, parse_smiles, rdkit_available

from .models import ComponentAnalysis, CompoundAnalysis, ReactiveSite, SiteType
from .sites import DETECTORS


def _component_molecules(mol: Any) -> List[Any]:
    from rdkit import Chem  # type: ignore
    return list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True))


def featurize_molecule(
    smiles: str,
    *,
    site_types: Optional[Iterable[SiteType]] = None,
    include_context_features: bool = True,
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

    selected: Set[str] = set(site_types or DETECTORS)
    unknown = selected - set(DETECTORS)
    if unknown:
        return CompoundAnalysis(text, mol_to_canonical_smiles(mol), False, error=f"UNKNOWN_SITE_TYPES:{','.join(sorted(unknown))}")

    components: List[ComponentAnalysis] = []
    all_sites: List[ReactiveSite] = []
    atom_offset = 0
    for component_index, component_mol in enumerate(_component_molecules(mol)):
        raw_sites = []
        for site_type, detector in DETECTORS.items():
            if site_type in selected:
                raw_sites.extend((site_type, item) for item in detector(component_mol))
        raw_sites.sort(key=lambda pair: (pair[1]["atom_indices"], pair[0], pair[1]["signature"]))
        component_sites: List[ReactiveSite] = []
        for site_number, (site_type, item) in enumerate(raw_sites):
            site = ReactiveSite(
                site_id=f"mol{component_index}:atom{item['atom_indices'][0]}:site{site_number}",
                site_type=site_type, topology=item["topology"], component_index=component_index,
                atom_indices=item["atom_indices"], bond_indices=[i for i in item["bond_indices"] if i >= 0],
                canonical_signature=item["signature"], chemist_label=item["label"],
                details=item.get("details", {}),
                context_features=item.get("context_features", {}) if include_context_features else {},
            )
            component_sites.append(site)
        canonical_component = mol_to_canonical_smiles(component_mol) or ""
        input_component = canonical_component
        components.append(ComponentAnalysis(component_index, input_component, canonical_component, atom_offset, component_sites))
        all_sites.extend(component_sites)
        atom_offset += component_mol.GetNumAtoms()

    return CompoundAnalysis(text, mol_to_canonical_smiles(mol), True, components, all_sites)


def detect_sites(mol: Any, site_types: Optional[Iterable[SiteType]] = None) -> List[ReactiveSite]:
    """Detect sites on an RDKit molecule through the public result contract."""
    smiles = mol_to_canonical_smiles(mol)
    if not smiles:
        return []
    return featurize_molecule(smiles, site_types=site_types).sites


__all__ = ["detect_sites", "featurize_molecule"]
