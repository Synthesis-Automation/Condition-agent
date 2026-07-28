"""Normalized, mechanism-neutral reactive-site interfaces.

Existing detectors retain their chemistry-specific site types and labels.
These adapters provide the smaller connectivity contracts consumed by generic
reaction rewrites without changing molecule or reaction serialization.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Literal, Mapping, Sequence, Tuple, cast

from .chemistry.rdkit_utils import parse_smiles
from .models import CompoundAnalysis, ReactiveSite
from .reaction_models import (
    ReactionComponent,
    ReactionSiteReference,
)


SITE_INTERFACE_SCHEMA_VERSION = "1.0"

EndpointKind = Literal["atom", "virtual_hydrogen"]
LinkSourceKind = Literal["explicit_bond", "implicit_hydrogen"]
EndpointSymmetry = Literal[
    "equivalent_endpoints",
    "distinguishable_endpoints",
    "virtual_hydrogen",
]

_BOND_ORDERS = {"SINGLE", "DOUBLE", "TRIPLE"}


@dataclass(frozen=True)
class ConnectivityEndpoint:
    """One real atom or schema-level virtual hydrogen endpoint."""

    endpoint_id: str
    endpoint_kind: EndpointKind
    component_index: int
    source_atom_role: str
    element: str
    atom_index: int | None
    carrier_atom_index: int | None
    formal_charge: int | None
    aromatic: bool | None
    hybridization: str | None
    current_valence: int | None
    hydrogen_count: int
    in_ring: bool | None
    local_environment: Tuple[str, ...] = ()
    context_tokens: Tuple[str, ...] = ()
    schema_version: str = SITE_INTERFACE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.endpoint_kind == "atom":
            if self.atom_index is None or self.carrier_atom_index is not None:
                raise ValueError("Atom endpoint requires only atom_index")
        elif self.endpoint_kind == "virtual_hydrogen":
            if (
                self.element != "H"
                or self.atom_index is not None
                or self.carrier_atom_index is None
            ):
                raise ValueError(
                    "Virtual hydrogen requires H and a carrier atom"
                )
        else:
            raise ValueError(f"Unsupported endpoint kind: {self.endpoint_kind}")
        if self.hydrogen_count < 0:
            raise ValueError("Endpoint hydrogen count cannot be negative")

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable endpoint view."""
        return asdict(self)


@dataclass(frozen=True)
class ReactiveLinkSite:
    """One explicit or schema-level link that a rewrite may consume."""

    site_id: str
    source_site_id: str
    source_site_type: str
    source_signature: str
    source_chemist_label: str
    endpoint_a: ConnectivityEndpoint
    endpoint_b: ConnectivityEndpoint
    before_order: str
    available_units: int
    source_kind: LinkSourceKind
    availability: str
    symmetry_class: EndpointSymmetry
    annotation_tokens: Tuple[str, ...] = ()
    schema_version: str = SITE_INTERFACE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.before_order not in _BOND_ORDERS:
            raise ValueError(f"Unsupported reactive-link order: {self.before_order}")
        if self.available_units < 1:
            raise ValueError("Reactive link requires at least one available unit")
        if self.endpoint_a.component_index != self.endpoint_b.component_index:
            raise ValueError("Reactive-link endpoints must share a component")
        if self.source_kind == "explicit_bond" and (
            self.endpoint_a.endpoint_kind != "atom"
            or self.endpoint_b.endpoint_kind != "atom"
        ):
            raise ValueError("Explicit reactive link requires two atom endpoints")
        if self.source_kind == "implicit_hydrogen" and (
            self.endpoint_a.endpoint_kind != "atom"
            or self.endpoint_b.endpoint_kind != "virtual_hydrogen"
        ):
            raise ValueError(
                "Implicit-H reactive link requires atom and virtual-H endpoints"
            )

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable link view."""
        return asdict(self)


@dataclass(frozen=True)
class BondCapacitySite:
    """One localized bond with bounded multiplicity capacity."""

    site_id: str
    source_site_id: str
    source_site_type: str
    source_signature: str
    source_chemist_label: str
    endpoint_a: ConnectivityEndpoint
    endpoint_b: ConnectivityEndpoint
    current_order: str
    maximum_decrement: int
    maximum_increment: int
    bond_class: str
    in_ring: bool
    aromatic: bool
    availability: str
    annotation_tokens: Tuple[str, ...] = ()
    schema_version: str = SITE_INTERFACE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.current_order not in _BOND_ORDERS:
            raise ValueError(f"Unsupported bond-capacity order: {self.current_order}")
        if self.maximum_decrement < 0 or self.maximum_increment < 0:
            raise ValueError("Bond capacity cannot be negative")
        if self.aromatic:
            raise ValueError("Aromatic bonds are outside the localized capacity domain")
        if self.endpoint_a.endpoint_kind != "atom" or (
            self.endpoint_b.endpoint_kind != "atom"
        ):
            raise ValueError("Bond capacity requires two atom endpoints")

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable bond-capacity view."""
        return asdict(self)


@dataclass(frozen=True)
class ConnectionEndpointSite:
    """One bounded atom endpoint that may participate in a new connection."""

    site_id: str
    source_site_id: str
    source_site_type: str
    source_signature: str
    source_chemist_label: str
    endpoint: ConnectivityEndpoint
    availability: str
    required_link_release_id: str | None = None
    required_hydrogen_delta: int = 0
    required_formal_charge_delta: int = 0
    required_bond_capacity_id: str | None = None
    required_bond_capacity_decrement: int = 0
    annotation_tokens: Tuple[str, ...] = ()
    schema_version: str = SITE_INTERFACE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.endpoint.endpoint_kind != "atom":
            raise ValueError("Connection endpoint must reference a real atom")
        if self.required_hydrogen_delta > 0:
            raise ValueError(
                "Connection endpoints may require H loss, not invented H gain"
            )
        if self.required_bond_capacity_decrement < 0:
            raise ValueError("Required capacity decrement cannot be negative")
        if bool(self.required_bond_capacity_id) != bool(
            self.required_bond_capacity_decrement
        ):
            raise ValueError(
                "Capacity ID and decrement must be declared together"
            )

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable connection-endpoint view."""
        return asdict(self)


@dataclass(frozen=True)
class NormalizedSiteInterfaces:
    """All normalized connectivity interfaces derived from one detector site."""

    source_site_id: str
    source_site_type: str
    source_signature: str
    source_chemist_label: str
    reactive_links: Tuple[ReactiveLinkSite, ...] = ()
    bond_capacities: Tuple[BondCapacitySite, ...] = ()
    connection_endpoints: Tuple[ConnectionEndpointSite, ...] = ()
    schema_version: str = SITE_INTERFACE_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable normalized site view."""
        return asdict(self)


def _context_tokens(site: ReactionSiteReference) -> Tuple[str, ...]:
    values = []
    for key in ("anchor_context", "ring_context", "retained_context"):
        if site.details.get(key):
            values.append(str(site.details[key]))
    values.extend(str(value) for value in site.details.get("contexts") or ())
    values.extend(
        str(record.get("token"))
        for record in site.details.get("context_records") or ()
        if record.get("token")
    )
    return tuple(sorted(set(values)))


def _interface_id(site: ReactionSiteReference, suffix: str) -> str:
    return (
        f"{site.side}:component{site.component_index}:"
        f"{site.site_id}:{suffix}"
    )


def _atom_endpoint(
    molecule: Any,
    site: ReactionSiteReference,
    *,
    atom_index: int,
    source_atom_role: str,
) -> ConnectivityEndpoint:
    atom = molecule.GetAtomWithIdx(atom_index)
    try:
        current_valence = int(atom.GetTotalValence())
    except (RuntimeError, ValueError):
        current_valence = int(atom.GetExplicitValence())
    local_environment = tuple(
        sorted(
            (
                f"{_order_name(bond) or 'UNSUPPORTED'}:"
                f"{bond.GetOtherAtom(atom).GetSymbol()}:"
                f"{bond.GetOtherAtom(atom).GetFormalCharge()}:"
                f"{int(bond.GetOtherAtom(atom).GetIsAromatic())}"
            )
            for bond in atom.GetBonds()
        )
    )
    return ConnectivityEndpoint(
        endpoint_id=_interface_id(site, f"endpoint:{source_atom_role}"),
        endpoint_kind="atom",
        component_index=site.component_index,
        source_atom_role=source_atom_role,
        element=atom.GetSymbol(),
        atom_index=atom_index,
        carrier_atom_index=None,
        formal_charge=int(atom.GetFormalCharge()),
        aromatic=bool(atom.GetIsAromatic()),
        hybridization=str(atom.GetHybridization()).upper(),
        current_valence=current_valence,
        hydrogen_count=int(atom.GetTotalNumHs(includeNeighbors=True)),
        in_ring=bool(atom.IsInRing()),
        local_environment=local_environment,
        context_tokens=_context_tokens(site),
    )


def _virtual_hydrogen_endpoint(
    site: ReactionSiteReference,
    *,
    carrier_atom_index: int,
    source_atom_role: str,
) -> ConnectivityEndpoint:
    return ConnectivityEndpoint(
        endpoint_id=_interface_id(site, "endpoint:virtual_hydrogen"),
        endpoint_kind="virtual_hydrogen",
        component_index=site.component_index,
        source_atom_role=source_atom_role,
        element="H",
        atom_index=None,
        carrier_atom_index=carrier_atom_index,
        formal_charge=0,
        aromatic=False,
        hybridization=None,
        current_valence=1,
        hydrogen_count=0,
        in_ring=False,
        local_environment=(),
        context_tokens=_context_tokens(site),
    )


def _single_role(
    site: ReactionSiteReference,
    role: str,
) -> int | None:
    indices = site.atom_roles.get(role) or ()
    return int(indices[0]) if len(indices) == 1 else None


def _bonded_endpoint(
    molecule: Any,
    site: ReactionSiteReference,
    *,
    anchor_index: int,
    candidate_roles: Sequence[str],
) -> tuple[int, str] | None:
    for role in candidate_roles:
        for raw_index in site.atom_roles.get(role) or ():
            atom_index = int(raw_index)
            if molecule.GetBondBetweenAtoms(anchor_index, atom_index) is not None:
                return atom_index, role
    return None


def _order_name(bond: Any) -> str | None:
    value = float(bond.GetBondTypeAsDouble())
    return {1.0: "SINGLE", 2.0: "DOUBLE", 3.0: "TRIPLE"}.get(value)


def _annotation_tokens(site: ReactionSiteReference) -> Tuple[str, ...]:
    values = [
        site.site_type,
        site.canonical_signature,
        site.chemist_label,
    ]
    for key in (
        "handle_token",
        "derived_family",
        "donor_class",
        "center_family",
    ):
        if site.details.get(key):
            values.append(str(site.details[key]))
    return tuple(dict.fromkeys(values))


def _endpoint_invariant(endpoint: ConnectivityEndpoint) -> tuple[Any, ...]:
    return (
        endpoint.endpoint_kind,
        endpoint.element,
        endpoint.formal_charge,
        endpoint.aromatic,
        endpoint.hybridization,
        endpoint.current_valence,
        endpoint.hydrogen_count,
        endpoint.in_ring,
        endpoint.local_environment,
        endpoint.context_tokens,
    )


def _link_from_site(
    site: ReactionSiteReference,
    molecule: Any,
) -> ReactiveLinkSite | None:
    endpoint_a_index: int | None = None
    endpoint_b_index: int | None = None
    endpoint_a_role = ""
    endpoint_b_role = ""
    source_kind: LinkSourceKind = "explicit_bond"
    available_units = 1

    if site.site_type == "leaving_group":
        endpoint_a_role = "anchor"
        endpoint_a_index = _single_role(site, endpoint_a_role)
        if endpoint_a_index is not None:
            bonded = _bonded_endpoint(
                molecule,
                site,
                anchor_index=endpoint_a_index,
                candidate_roles=("connector", "handle", "center"),
            )
            if bonded is not None:
                endpoint_b_index, _ = bonded
                # Keep the detector's compatibility role even when a
                # multi-atom handle uses a bonded connector internally.
                endpoint_b_role = "handle"
    elif site.site_type == "transfer_group":
        endpoint_a_role, endpoint_b_role = "anchor", "center"
        endpoint_a_index = _single_role(site, endpoint_a_role)
        endpoint_b_index = _single_role(site, endpoint_b_role)
    elif site.site_type in {"pronucleophile_XH", "aromatic_CH"}:
        endpoint_a_role = "center"
        endpoint_a_index = _single_role(site, endpoint_a_role)
        if endpoint_a_index is not None:
            source_kind = "implicit_hydrogen"
            available_units = int(site.details.get("h_count") or 0)
    elif site.site_type == "addition_donor":
        endpoint_a_role = "addend_a"
        endpoint_a_index = _single_role(site, endpoint_a_role)
        raw_source_kind = str(site.details.get("source_kind") or "")
        if raw_source_kind not in {"explicit_bond", "implicit_hydrogen"}:
            return None
        source_kind = cast(LinkSourceKind, raw_source_kind)
        if source_kind == "explicit_bond":
            endpoint_b_role = "addend_b"
            endpoint_b_index = _single_role(site, endpoint_b_role)
        elif source_kind == "implicit_hydrogen":
            endpoint_b_role = "hydrogen_carrier"
            available_units = int(site.details.get("hydrogen_count") or 0)
        else:
            return None
    elif site.site_type == "electrophilic_center":
        endpoint_a_role, endpoint_b_role = (
            "center",
            "leaving_or_activatable",
        )
        endpoint_a_index = _single_role(site, endpoint_a_role)
        endpoint_b_index = _single_role(site, endpoint_b_role)
    elif site.site_type == "eliminable_pair":
        endpoint_a_role, endpoint_b_role = "endpoint_a", "departing_a"
        endpoint_a_index = _single_role(site, endpoint_a_role)
        endpoint_b_index = _single_role(site, endpoint_b_role)
    else:
        return None

    if endpoint_a_index is None:
        return None
    endpoint_a = _atom_endpoint(
        molecule,
        site,
        atom_index=endpoint_a_index,
        source_atom_role=endpoint_a_role,
    )
    if source_kind == "implicit_hydrogen":
        if available_units < 1:
            return None
        endpoint_b = _virtual_hydrogen_endpoint(
            site,
            carrier_atom_index=endpoint_a_index,
            source_atom_role=endpoint_b_role or endpoint_a_role,
        )
        symmetry: EndpointSymmetry = "virtual_hydrogen"
    else:
        if endpoint_b_index is None:
            return None
        bond = molecule.GetBondBetweenAtoms(endpoint_a_index, endpoint_b_index)
        before_order = _order_name(bond) if bond is not None else None
        if before_order is None:
            return None
        endpoint_b = _atom_endpoint(
            molecule,
            site,
            atom_index=endpoint_b_index,
            source_atom_role=endpoint_b_role,
        )
        symmetry = (
            "equivalent_endpoints"
            if _endpoint_invariant(endpoint_a) == _endpoint_invariant(endpoint_b)
            else "distinguishable_endpoints"
        )
        return ReactiveLinkSite(
            site_id=_interface_id(site, "reactive_link"),
            source_site_id=site.site_id,
            source_site_type=site.site_type,
            source_signature=site.canonical_signature,
            source_chemist_label=site.chemist_label,
            endpoint_a=endpoint_a,
            endpoint_b=endpoint_b,
            before_order=before_order,
            available_units=1,
            source_kind=source_kind,
            availability=site.availability,
            symmetry_class=symmetry,
            annotation_tokens=_annotation_tokens(site),
        )
    return ReactiveLinkSite(
        site_id=_interface_id(site, "reactive_link"),
        source_site_id=site.site_id,
        source_site_type=site.site_type,
        source_signature=site.canonical_signature,
        source_chemist_label=site.chemist_label,
        endpoint_a=endpoint_a,
        endpoint_b=endpoint_b,
        before_order="SINGLE",
        available_units=available_units,
        source_kind=source_kind,
        availability=site.availability,
        symmetry_class=symmetry,
        annotation_tokens=_annotation_tokens(site),
    )


def _capacity_from_site(
    site: ReactionSiteReference,
    molecule: Any,
) -> BondCapacitySite | None:
    if site.site_type == "unsaturated_bond":
        endpoint_a_role, endpoint_b_role = "endpoint_a", "endpoint_b"
        bond_class = {
            "Alkene": "carbon_carbon_pi",
            "Alkyne": "carbon_carbon_pi",
            "Nitrile": "carbon_nitrogen_pi",
        }.get(
            str(site.details.get("handle_token") or ""),
            "localized_multiple_bond",
        )
    elif (
        site.site_type == "electrophilic_center"
        and "heteroatom" in site.atom_roles
    ):
        endpoint_a_role, endpoint_b_role = "center", "heteroatom"
        bond_class = "polarized_multiple_bond"
    elif (
        site.site_type == "pronucleophile_XH"
        and site.details.get("derived_family") == "alcohol"
    ):
        endpoint_a_role, endpoint_b_role = "center", "attachment"
        bond_class = "alcohol_carbon_heteroatom"
    elif site.site_type == "eliminable_pair":
        endpoint_a_role, endpoint_b_role = "endpoint_a", "endpoint_b"
        bond_class = "eliminable_backbone"
    else:
        return None
    endpoint_a_index = _single_role(site, endpoint_a_role)
    endpoint_b_index = _single_role(site, endpoint_b_role)
    if endpoint_a_index is None or endpoint_b_index is None:
        return None
    bond = molecule.GetBondBetweenAtoms(endpoint_a_index, endpoint_b_index)
    current_order = _order_name(bond) if bond is not None else None
    if current_order is None or bond.GetIsAromatic():
        return None
    order_units = {"SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3}[current_order]
    endpoint_a = _atom_endpoint(
        molecule,
        site,
        atom_index=endpoint_a_index,
        source_atom_role=endpoint_a_role,
    )
    endpoint_b = _atom_endpoint(
        molecule,
        site,
        atom_index=endpoint_b_index,
        source_atom_role=endpoint_b_role,
    )
    return BondCapacitySite(
        site_id=_interface_id(site, "bond_capacity"),
        source_site_id=site.site_id,
        source_site_type=site.site_type,
        source_signature=site.canonical_signature,
        source_chemist_label=site.chemist_label,
        endpoint_a=endpoint_a,
        endpoint_b=endpoint_b,
        current_order=current_order,
        maximum_decrement=max(0, order_units - 1),
        maximum_increment=max(0, 3 - order_units),
        bond_class=bond_class,
        in_ring=bool(bond.IsInRing()),
        aromatic=False,
        availability=site.availability,
        annotation_tokens=_annotation_tokens(site),
    )


def _connection_endpoints(
    site: ReactionSiteReference,
    molecule: Any,
    link: ReactiveLinkSite | None,
    capacity: BondCapacitySite | None,
) -> Tuple[ConnectionEndpointSite, ...]:
    endpoints = []
    annotations = _annotation_tokens(site)
    if site.site_type == "nucleophile_anion":
        center_index = _single_role(site, "center")
        if center_index is None:
            return ()
        endpoint = _atom_endpoint(
            molecule,
            site,
            atom_index=center_index,
            source_atom_role="center",
        )
        return (
            ConnectionEndpointSite(
                site_id=_interface_id(site, "connection:center"),
                source_site_id=site.site_id,
                source_site_type=site.site_type,
                source_signature=site.canonical_signature,
                source_chemist_label=site.chemist_label,
                endpoint=endpoint,
                availability=site.availability,
                required_formal_charge_delta=1,
                annotation_tokens=annotations,
            ),
        )
    if site.site_type == "eliminable_pair":
        carrier_index = _single_role(site, "hydrogen_carrier_b")
        if carrier_index is None:
            return ()
        endpoint = _atom_endpoint(
            molecule,
            site,
            atom_index=carrier_index,
            source_atom_role="hydrogen_carrier_b",
        )
        return (
            ConnectionEndpointSite(
                site_id=_interface_id(
                    site, "connection:hydrogen_carrier_b"
                ),
                source_site_id=site.site_id,
                source_site_type=site.site_type,
                source_signature=site.canonical_signature,
                source_chemist_label=site.chemist_label,
                endpoint=endpoint,
                availability=site.availability,
                required_hydrogen_delta=-1,
                annotation_tokens=annotations,
            ),
        )
    if link is not None:
        for endpoint in (link.endpoint_a, link.endpoint_b):
            if endpoint.endpoint_kind != "atom":
                continue
            endpoints.append(
                ConnectionEndpointSite(
                    site_id=(
                        _interface_id(
                            site,
                            f"connection:{endpoint.source_atom_role}",
                        )
                    ),
                    source_site_id=site.site_id,
                    source_site_type=site.site_type,
                    source_signature=site.canonical_signature,
                    source_chemist_label=site.chemist_label,
                    endpoint=endpoint,
                    availability=site.availability,
                    required_link_release_id=(
                        link.site_id
                        if link.source_kind == "explicit_bond"
                        else None
                    ),
                    required_hydrogen_delta=(
                        -1
                        if link.source_kind == "implicit_hydrogen"
                        else 0
                    ),
                    annotation_tokens=annotations,
                )
            )
    if capacity is not None and link is None:
        for endpoint in (capacity.endpoint_a, capacity.endpoint_b):
            endpoints.append(
                ConnectionEndpointSite(
                    site_id=(
                        _interface_id(
                            site,
                            f"connection:{endpoint.source_atom_role}",
                        )
                    ),
                    source_site_id=site.site_id,
                    source_site_type=site.site_type,
                    source_signature=site.canonical_signature,
                    source_chemist_label=site.chemist_label,
                    endpoint=endpoint,
                    availability=site.availability,
                    required_bond_capacity_id=capacity.site_id,
                    required_bond_capacity_decrement=1,
                    annotation_tokens=annotations,
                )
            )
    return tuple(endpoints)


def normalize_reaction_site(
    site: ReactionSiteReference,
    component: ReactionComponent,
) -> NormalizedSiteInterfaces:
    """Adapt one existing detector site to normalized connectivity interfaces."""
    if site.component_index != component.component_index:
        raise ValueError("Site and component indices do not match")
    molecule = parse_smiles(component.input_smiles)
    if molecule is None:
        raise ValueError("Cannot normalize a site from an invalid component")
    return _normalize_with_molecule(site, molecule)


def _normalize_with_molecule(
    site: ReactionSiteReference,
    molecule: Any,
) -> NormalizedSiteInterfaces:
    link = _link_from_site(site, molecule)
    capacity = _capacity_from_site(site, molecule)
    return NormalizedSiteInterfaces(
        source_site_id=site.site_id,
        source_site_type=site.site_type,
        source_signature=site.canonical_signature,
        source_chemist_label=site.chemist_label,
        reactive_links=(link,) if link is not None else (),
        bond_capacities=(capacity,) if capacity is not None else (),
        connection_endpoints=_connection_endpoints(
            site, molecule, link, capacity
        ),
    )


def normalize_detected_site(
    site: ReactiveSite,
    component_molecule: Any,
) -> NormalizedSiteInterfaces:
    """Adapt a detector site using its atom-index-identical RDKit component."""
    if not hasattr(component_molecule, "GetAtomWithIdx"):
        raise ValueError("Detected-site normalization requires an RDKit molecule")
    raw_roles = site.details.get("atom_roles") or {}
    reference = ReactionSiteReference(
        side="reactant",
        component_index=site.component_index,
        site_id=site.site_id,
        site_type=site.site_type,
        canonical_signature=site.canonical_signature,
        chemist_label=site.chemist_label,
        availability=site.availability,
        atom_roles={
            str(role): tuple(int(index) for index in indices)
            for role, indices in raw_roles.items()
        },
        details=dict(site.details),
    )
    return _normalize_with_molecule(reference, component_molecule)


def normalize_compound_sites(
    analysis: CompoundAnalysis,
) -> Tuple[NormalizedSiteInterfaces, ...]:
    """Normalize all sites while preserving detector atom-index provenance."""
    from rdkit import Chem

    molecule = parse_smiles(analysis.input_smiles)
    if molecule is None:
        raise ValueError("Cannot normalize sites from an invalid compound")
    component_molecules = tuple(
        Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)
    )
    if len(component_molecules) != len(analysis.components):
        raise ValueError("Compound components do not match detector provenance")
    if any(
        site.component_index < 0
        or site.component_index >= len(component_molecules)
        for site in analysis.sites
    ):
        raise ValueError("Detected site has invalid component provenance")
    return tuple(
        normalize_detected_site(
            site,
            component_molecules[site.component_index],
        )
        for site in analysis.sites
    )


def normalize_reaction_assignment(
    assignment: Mapping[str, ReactionSiteReference],
    components: Sequence[ReactionComponent],
) -> Dict[str, NormalizedSiteInterfaces]:
    """Normalize every grammar role assignment without changing its source site."""
    by_index = {
        component.component_index: component for component in components
    }
    missing_components = {
        site.component_index for site in assignment.values()
    } - set(by_index)
    if missing_components:
        raise ValueError("Reaction assignment has missing component provenance")
    return {
        role: normalize_reaction_site(site, by_index[site.component_index])
        for role, site in assignment.items()
    }


__all__ = [
    "BondCapacitySite",
    "ConnectionEndpointSite",
    "ConnectivityEndpoint",
    "EndpointKind",
    "EndpointSymmetry",
    "LinkSourceKind",
    "NormalizedSiteInterfaces",
    "ReactiveLinkSite",
    "SITE_INTERFACE_SCHEMA_VERSION",
    "normalize_compound_sites",
    "normalize_detected_site",
    "normalize_reaction_assignment",
    "normalize_reaction_site",
]
