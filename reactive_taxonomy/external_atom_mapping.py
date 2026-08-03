"""Optional atom-mapping providers with structure-preserving validation.

External mappers propose atom correspondence.  They do not become chemistry
truth merely because they emit atom-map numbers: every result is reparsed,
checked against the original reaction, and normalized through the same typed
edit extractor used for supplied mappings.
"""

from __future__ import annotations

import hashlib
import importlib.util
from dataclasses import dataclass, replace
from importlib import metadata, resources
from typing import Any, Iterable, Literal, Optional, Protocol, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_connectivity import build_connectivity_edit_graph
from .reaction_edits import EditNormalizationResult, normalize_mapped_edits
from .reaction_models import (
    ReactionAnalysis,
    ReactionAtomReference,
    ReactionComponent,
    ReactionCoreProjection,
    ReactionEdit,
    ReactionEvidenceCandidate,
)
from .reaction_parser import ParsedReaction, parse_reaction_smiles
from .reaction_rendering import render_reaction
from .reaction_render_context import build_reaction_render_context


EXTERNAL_ATOM_MAPPING_SCHEMA_VERSION = "1.0"
EXTERNAL_MAPPING_EVIDENCE = "external_atom_mapping"


def _attach_core_and_rerender(
    base: ReactionAnalysis,
    reaction_core: ReactionCoreProjection,
    *,
    warnings: Iterable[str],
) -> ReactionAnalysis:
    """Attach a compatible core to both contracts and refresh the label.

    External mapping may add retained-context information after the internal
    reaction interpretation was completed. Identity, signature, and family
    interpretation remain unchanged; only the evidence-aware terminal label is
    rendered again from the enriched observation.
    """
    combined_warnings = tuple(sorted(set(base.warnings).union(warnings)))
    observation = base.observation
    if observation is None:
        return replace(
            base,
            reaction_core=reaction_core,
            warnings=combined_warnings,
        )
    enriched_observation = replace(
        observation,
        core=reaction_core,
        warnings=tuple(
            sorted(set(observation.warnings).union(warnings))
        ),
    )
    return replace(
        base,
        reaction_core=reaction_core,
        observation=enriched_observation,
        reaction_label=render_reaction(build_reaction_render_context(
            observation=enriched_observation,
            reactants=base.reactants,
            agents=base.agents,
            products=base.products,
            signature=base.reaction_signature,
            partial_transformation=base.partial_product_transformation,
            style=(base.reaction_label.style if base.reaction_label else "unicode"),
            interpretation=base.interpretation,
        )),
        warnings=combined_warnings,
    )


@dataclass(frozen=True)
class AtomMappingProviderMetadata:
    """Versioned identity of one external atom-mapping implementation."""

    provider_id: str
    provider_version: str
    model_id: str
    model_sha256: Optional[str]
    schema_version: str = EXTERNAL_ATOM_MAPPING_SCHEMA_VERSION


@dataclass(frozen=True)
class ExternalAtomMappingResult:
    """One mapper proposal plus local validation and normalized edits."""

    input_reaction_smiles: str
    mapped_reaction_smiles: Optional[str]
    metadata: AtomMappingProviderMetadata
    mapper_confidence: Optional[float]
    structure_preserved: bool
    reactant_mapping_coverage: float
    product_mapping_coverage: float
    shared_mapped_heavy_atom_count: int
    normalization: Optional[EditNormalizationResult]
    valid: bool
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = EXTERNAL_ATOM_MAPPING_SCHEMA_VERSION


ExternalMappingAssessmentStatus = Literal[
    "not_requested_invalid_reaction",
    "not_requested_supplied_mapping",
    "not_requested_resolved_internal_evidence",
    "external_mapping_failed",
    "external_mapping_internal_consensus",
    "external_mapping_hypothesis_conflict",
    "external_mapping_signature_conflict",
    "external_mapping_ambiguous_hypothesis_match",
    "external_mapping_only",
    "external_mapping_signature_unavailable",
]


@dataclass(frozen=True)
class ExternalMappingAssessment:
    """Effective analysis plus auditable external-correspondence disposition."""

    input_reaction_smiles: str
    status: ExternalMappingAssessmentStatus
    analysis: ReactionAnalysis
    provider_metadata: AtomMappingProviderMetadata
    mapping_result: Optional[ExternalAtomMappingResult] = None
    matched_hypothesis_ids: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()
    schema_version: str = EXTERNAL_ATOM_MAPPING_SCHEMA_VERSION

    def to_provenance_dict(self) -> dict[str, Any]:
        """Serialize provider evidence without duplicating full reaction analysis."""
        mapping = self.mapping_result
        return {
            "schema_version": self.schema_version,
            "status": self.status,
            "provider": {
                "provider_id": self.provider_metadata.provider_id,
                "provider_version": self.provider_metadata.provider_version,
                "model_id": self.provider_metadata.model_id,
                "model_sha256": self.provider_metadata.model_sha256,
                "schema_version": self.provider_metadata.schema_version,
            },
            "mapper_confidence": (
                mapping.mapper_confidence if mapping is not None else None
            ),
            "structure_preserved": (
                mapping.structure_preserved if mapping is not None else None
            ),
            "reactant_mapping_coverage": (
                mapping.reactant_mapping_coverage if mapping is not None else None
            ),
            "product_mapping_coverage": (
                mapping.product_mapping_coverage if mapping is not None else None
            ),
            "shared_mapped_heavy_atom_count": (
                mapping.shared_mapped_heavy_atom_count
                if mapping is not None
                else None
            ),
            "mapped_reaction_smiles": (
                mapping.mapped_reaction_smiles if mapping is not None else None
            ),
            "matched_hypothesis_ids": self.matched_hypothesis_ids,
            "warnings": self.warnings,
            "error": mapping.error if mapping is not None else None,
        }


class AtomMappingProvider(Protocol):
    """Narrow injectable provider contract; implementations may be optional."""

    @property
    def metadata(self) -> AtomMappingProviderMetadata:
        """Return immutable provider and model provenance."""

    def map_reactions(
        self, reaction_smiles: Iterable[str]
    ) -> Tuple[ExternalAtomMappingResult, ...]:
        """Map and locally validate reactions in input order."""


def _canonical_without_maps(smiles: str) -> Optional[str]:
    from rdkit import Chem

    molecule = parse_smiles(smiles)
    if molecule is None:
        return None
    copy = Chem.Mol(molecule)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        # Atom mapping changes atom order, which can reverse the textual @/@@
        # tag while preserving the same stereochemical structure.  Reassign
        # stereo after removing maps so identity comparison uses normalized
        # molecular stereochemistry rather than serialization-local tags.
        Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
        return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def _side_identity(
    components: Tuple[ReactionComponent, ...],
) -> Optional[Tuple[str, ...]]:
    identities = tuple(
        _canonical_without_maps(component.input_smiles)
        for component in components
    )
    if any(identity is None for identity in identities):
        return None
    return tuple(sorted(str(identity) for identity in identities))


def _reaction_identity(
    parsed: ParsedReaction,
) -> Optional[Tuple[Tuple[str, ...], Tuple[str, ...], Tuple[str, ...]]]:
    reactants = _side_identity(parsed.reactants)
    agents = _side_identity(parsed.agents)
    products = _side_identity(parsed.products)
    if reactants is None or agents is None or products is None:
        return None
    return reactants, agents, products


def _indexed_graph(molecule: Any) -> tuple[object, ...]:
    """Return atom-order-sensitive graph data while ignoring map numbers."""
    atoms = tuple(
        (
            atom.GetSymbol(),
            int(atom.GetFormalCharge()),
            bool(atom.GetIsAromatic()),
        )
        for atom in molecule.GetAtoms()
    )
    bonds = tuple(
        sorted(
            (
                min(int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())),
                max(int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())),
                str(bond.GetBondType()).upper(),
            )
            for bond in molecule.GetBonds()
        )
    )
    return atoms, bonds


def _component_map_numbers_in_original_order(
    original: ReactionComponent,
    mapped: ReactionComponent,
) -> Optional[dict[int, int]]:
    """Project validated map numbers onto original atom coordinates."""
    from rdkit import Chem

    original_molecule = parse_smiles(original.input_smiles)
    mapped_molecule = parse_smiles(mapped.input_smiles)
    if original_molecule is None or mapped_molecule is None:
        return None
    query = Chem.Mol(mapped_molecule)
    for atom in query.GetAtoms():
        atom.SetAtomMapNum(0)
    matches = tuple(
        match
        for match in original_molecule.GetSubstructMatches(
            query,
            uniquify=False,
            useChirality=True,
            maxMatches=1000,
        )
        if len(match) == original_molecule.GetNumAtoms()
    )
    if not matches:
        return None
    match = min(matches)
    return {
        int(original_index): int(
            mapped_molecule.GetAtomWithIdx(mapped_index).GetAtomMapNum()
        )
        for mapped_index, original_index in enumerate(match)
        if int(mapped_molecule.GetAtomWithIdx(mapped_index).GetAtomMapNum()) > 0
    }


def _mapping_in_original_coordinates(
    base: ReactionAnalysis,
    mapped: ParsedReaction,
) -> Optional[dict[tuple[str, int, int], int]]:
    """Return validated map numbers keyed by base-analysis coordinates."""

    def project_side(
        side: str,
        originals: Sequence[ReactionComponent],
        mapped_components: Sequence[ReactionComponent],
    ) -> Optional[dict[tuple[str, int, int], int]]:
        available: dict[str, list[ReactionComponent]] = {}
        for component in mapped_components:
            identity = _canonical_without_maps(component.input_smiles)
            if identity is None:
                return None
            available.setdefault(identity, []).append(component)
        for values in available.values():
            values.sort(key=lambda component: component.component_index)
        projected: dict[tuple[str, int, int], int] = {}
        for original in originals:
            identity = _canonical_without_maps(original.input_smiles)
            candidates = available.get(str(identity), [])
            if identity is None or not candidates:
                return None
            mapped_component = candidates.pop(0)
            map_numbers = _component_map_numbers_in_original_order(
                original,
                mapped_component,
            )
            if map_numbers is None:
                return None
            projected.update(
                {
                    (side, original.component_index, atom_index): map_number
                    for atom_index, map_number in map_numbers.items()
                }
            )
        if any(available.values()):
            return None
        return projected

    reactants = project_side("reactant", base.reactants, mapped.reactants)
    products = project_side("product", base.products, mapped.products)
    if reactants is None or products is None:
        return None
    return {**reactants, **products}


def _mapping_coverage(
    components: Tuple[ReactionComponent, ...],
) -> tuple[int, int, set[int]]:
    total = 0
    mapped = 0
    map_numbers: set[int] = set()
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() <= 1:
                continue
            total += 1
            map_number = int(atom.GetAtomMapNum())
            if map_number > 0:
                mapped += 1
                map_numbers.add(map_number)
    return mapped, total, map_numbers


def _project_reactant_boundary_atoms(
    mapped_reaction_smiles: str,
) -> tuple[str, int]:
    """Identify only retained-to-unreported reactant attachment boundaries.

    RXNMapper normally maps atoms retained in the reported product and leaves
    departing atoms unnumbered.  Giving an otherwise unmapped boundary atom a
    reactant-only map number lets the existing normalizer represent the
    observed attachment loss without pretending to map the entire unreported
    fragment or an unmapped spectator component.
    """
    from rdkit import Chem

    if ">>" in mapped_reaction_smiles:
        left, right = mapped_reaction_smiles.split(">>", 1)
        middle = ""
        separator = ">>"
    else:
        parts = mapped_reaction_smiles.split(">")
        if len(parts) != 3:
            return mapped_reaction_smiles, 0
        left, middle, right = parts
        separator = ">"
    molecules = [
        molecule
        for token in tuple(left.split(".")) + tuple(right.split("."))
        for molecule in (parse_smiles(token),)
        if molecule is not None
    ]
    next_map = (
        max(
            (
                int(atom.GetAtomMapNum())
                for molecule in molecules
                for atom in molecule.GetAtoms()
            ),
            default=0,
        )
        + 1
    )
    projected = 0
    rendered_reactants = []
    for token in left.split("."):
        molecule = parse_smiles(token)
        if molecule is None:
            rendered_reactants.append(token)
            continue
        initially_mapped = {
            int(atom.GetIdx())
            for atom in molecule.GetAtoms()
            if int(atom.GetAtomMapNum()) > 0
        }
        boundary_indices = tuple(
            sorted(
                int(atom.GetIdx())
                for atom in molecule.GetAtoms()
                if atom.GetAtomicNum() > 1
                and int(atom.GetAtomMapNum()) == 0
                and any(
                    int(neighbor.GetIdx()) in initially_mapped
                    for neighbor in atom.GetNeighbors()
                )
            )
        )
        for atom_index in boundary_indices:
            molecule.GetAtomWithIdx(atom_index).SetAtomMapNum(next_map)
            next_map += 1
            projected += 1
        rendered_reactants.append(
            Chem.MolToSmiles(molecule, canonical=False, isomericSmiles=True)
        )
    rendered_left = ".".join(rendered_reactants)
    if separator == ">>":
        return f"{rendered_left}>>{right}", projected
    return f"{rendered_left}>{middle}>{right}", projected


def _externalize_normalization(
    normalization: EditNormalizationResult,
    confidence: float,
) -> EditNormalizationResult:
    edits = tuple(
        replace(
            edit,
            evidence=EXTERNAL_MAPPING_EVIDENCE,
            confidence=confidence,
        )
        for edit in normalization.edits
    )
    stereo_changes = tuple(
        replace(
            change,
            evidence=EXTERNAL_MAPPING_EVIDENCE,
            confidence=confidence,
        )
        for change in normalization.stereo_changes
    )
    source_graph = normalization.connectivity_edit_graph
    connectivity_graph = None
    if source_graph is not None:
        connectivity_graph = build_connectivity_edit_graph(
            bond_transitions=tuple(
                replace(
                    transition,
                    evidence=EXTERNAL_MAPPING_EVIDENCE,
                    confidence=confidence,
                )
                for transition in source_graph.bond_transitions
            ),
            hydrogen_deltas=tuple(
                replace(
                    delta,
                    evidence=EXTERNAL_MAPPING_EVIDENCE,
                    confidence=confidence,
                )
                for delta in source_graph.hydrogen_deltas
            ),
            atom_state_transitions=tuple(
                replace(
                    transition,
                    evidence=EXTERNAL_MAPPING_EVIDENCE,
                    confidence=confidence,
                )
                for transition in source_graph.atom_state_transitions
            ),
            evidence=EXTERNAL_MAPPING_EVIDENCE,
            confidence=confidence,
            warnings=source_graph.warnings,
        )
    return replace(
        normalization,
        edits=edits,
        evidence=EXTERNAL_MAPPING_EVIDENCE,
        confidence=confidence,
        warnings=tuple(
            sorted(
                set(normalization.warnings)
                | {"EXTERNAL_MAPPING_REQUIRES_INDEPENDENT_VALIDATION"}
            )
        ),
        stereo_changes=stereo_changes,
        connectivity_edit_graph=connectivity_graph,
    )


def validate_external_atom_mapping(
    input_reaction_smiles: str,
    mapped_reaction_smiles: str,
    *,
    provider_metadata: AtomMappingProviderMetadata,
    mapper_confidence: Optional[float],
) -> ExternalAtomMappingResult:
    """Validate one external proposal without promoting it to supplied evidence."""
    original = parse_reaction_smiles(input_reaction_smiles)
    mapped = parse_reaction_smiles(mapped_reaction_smiles)
    warnings: set[str] = {"EXTERNAL_MAPPING_REQUIRES_INDEPENDENT_VALIDATION"}
    if not original.valid:
        return ExternalAtomMappingResult(
            input_reaction_smiles=input_reaction_smiles,
            mapped_reaction_smiles=mapped_reaction_smiles,
            metadata=provider_metadata,
            mapper_confidence=mapper_confidence,
            structure_preserved=False,
            reactant_mapping_coverage=0.0,
            product_mapping_coverage=0.0,
            shared_mapped_heavy_atom_count=0,
            normalization=None,
            valid=False,
            warnings=tuple(sorted(warnings)),
            error="INVALID_INPUT_REACTION",
        )
    if not mapped.valid:
        warnings.update(mapped.warnings)
        return ExternalAtomMappingResult(
            input_reaction_smiles=input_reaction_smiles,
            mapped_reaction_smiles=mapped_reaction_smiles,
            metadata=provider_metadata,
            mapper_confidence=mapper_confidence,
            structure_preserved=False,
            reactant_mapping_coverage=0.0,
            product_mapping_coverage=0.0,
            shared_mapped_heavy_atom_count=0,
            normalization=None,
            valid=False,
            warnings=tuple(sorted(warnings)),
            error="INVALID_MAPPED_REACTION",
        )
    structure_preserved = _reaction_identity(original) == _reaction_identity(mapped)
    reactant_mapped, reactant_total, reactant_maps = _mapping_coverage(
        mapped.reactants
    )
    product_mapped, product_total, product_maps = _mapping_coverage(mapped.products)
    reactant_coverage = reactant_mapped / reactant_total if reactant_total else 0.0
    product_coverage = product_mapped / product_total if product_total else 0.0
    shared_count = len(reactant_maps.intersection(product_maps))
    if not structure_preserved:
        warnings.add("EXTERNAL_MAPPER_CHANGED_REACTION_STRUCTURE")
        return ExternalAtomMappingResult(
            input_reaction_smiles=input_reaction_smiles,
            mapped_reaction_smiles=mapped_reaction_smiles,
            metadata=provider_metadata,
            mapper_confidence=mapper_confidence,
            structure_preserved=False,
            reactant_mapping_coverage=reactant_coverage,
            product_mapping_coverage=product_coverage,
            shared_mapped_heavy_atom_count=shared_count,
            normalization=None,
            valid=False,
            warnings=tuple(sorted(warnings)),
            error="MAPPER_CHANGED_REACTION_STRUCTURE",
        )
    normalized_mapped_smiles, projected_boundary_count = (
        _project_reactant_boundary_atoms(mapped_reaction_smiles)
    )
    if projected_boundary_count:
        warnings.add(
            "EXTERNAL_MAPPING_PROJECTED_REACTANT_BOUNDARY_ATOMS:"
            f"{projected_boundary_count}"
        )
        mapped = parse_reaction_smiles(normalized_mapped_smiles)
    normalized = normalize_mapped_edits(mapped.reactants, mapped.products)
    confidence = (
        float(mapper_confidence)
        if mapper_confidence is not None and 0.0 <= mapper_confidence <= 1.0
        else 0.0
    )
    normalized = _externalize_normalization(normalized, confidence)
    warnings.update(normalized.warnings)
    valid = normalized.valid and bool(normalized.edits)
    return ExternalAtomMappingResult(
        input_reaction_smiles=input_reaction_smiles,
        mapped_reaction_smiles=normalized_mapped_smiles,
        metadata=provider_metadata,
        mapper_confidence=mapper_confidence,
        structure_preserved=True,
        reactant_mapping_coverage=reactant_coverage,
        product_mapping_coverage=product_coverage,
        shared_mapped_heavy_atom_count=shared_count,
        normalization=normalized,
        valid=valid,
        warnings=tuple(sorted(warnings)),
        error=None if valid else normalized.evidence,
    )


def _sha256_file(path: Any) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _component_by_index(
    components: Sequence[ReactionComponent],
    index: int,
) -> Optional[ReactionComponent]:
    return next(
        (
            component
            for component in components
            if component.component_index == index
        ),
        None,
    )


def _edit_endpoint_token(
    atom: Optional[ReactionAtomReference],
    reactants: Sequence[ReactionComponent],
    products: Sequence[ReactionComponent],
) -> tuple[Any, ...]:
    if atom is None:
        return ("H",)
    components = products if atom.side == "product" else reactants
    component = _component_by_index(components, atom.component_index)
    return (
        atom.side,
        (
            _canonical_without_maps(component.input_smiles)
            if component is not None
            else "?"
        ),
        atom.element,
        atom.formal_charge,
        atom.aromatic,
        atom.hybridization,
        atom.local_environment_id,
    )


def normalized_edit_profile(
    edits: Sequence[ReactionEdit],
    reactants: Sequence[ReactionComponent],
    products: Sequence[ReactionComponent],
) -> tuple[tuple[Any, ...], ...]:
    """Return a serialization- and component-order-independent edit profile."""
    tokens = []
    for edit in edits:
        endpoints = tuple(
            sorted(
                (
                    _edit_endpoint_token(edit.atom_1, reactants, products),
                    _edit_endpoint_token(edit.atom_2, reactants, products),
                ),
                key=lambda value: repr(value),
            )
        )
        tokens.append(
            (
                edit.edit_type,
                endpoints,
                edit.old_order or "NONE",
                edit.new_order or "NONE",
            )
        )
    return tuple(sorted(tokens, key=repr))


def _external_candidate(
    result: ExternalAtomMappingResult,
) -> ReactionEvidenceCandidate:
    normalization = result.normalization
    if result.valid and normalization is not None:
        return ReactionEvidenceCandidate(
            provider=result.metadata.provider_id,
            status="verified",
            evidence=EXTERNAL_MAPPING_EVIDENCE,
            confidence=normalization.confidence,
            edits=normalization.edits,
            stereo_changes=normalization.stereo_changes,
            warnings=result.warnings,
        )
    return ReactionEvidenceCandidate(
        provider=result.metadata.provider_id,
        status="invalid",
        evidence=EXTERNAL_MAPPING_EVIDENCE,
        confidence=float(result.mapper_confidence or 0.0),
        warnings=result.warnings
        + ((result.error,) if result.error is not None else ()),
    )


def _merge_evidence_candidates(
    *candidate_groups: Sequence[ReactionEvidenceCandidate],
) -> tuple[ReactionEvidenceCandidate, ...]:
    values = []
    seen = set()
    for candidate in (
        candidate for group in candidate_groups for candidate in group
    ):
        key = (
            candidate.provider,
            candidate.status,
            candidate.evidence,
            tuple(
                (
                    edit.edit_type,
                    edit.atom_1.local_environment_id,
                    (
                        edit.atom_2.local_environment_id
                        if edit.atom_2 is not None
                        else "H"
                    ),
                    edit.old_order,
                    edit.new_order,
                )
                for edit in candidate.edits
            ),
        )
        if key in seen:
            continue
        seen.add(key)
        values.append(candidate)
    return tuple(values)


def _mapped_core_proposal(
    *,
    base: ReactionAnalysis,
    mapped: ParsedReaction,
    mapping: ExternalAtomMappingResult,
    evidence: str,
    confidence: float,
) -> tuple[Optional[ReactionCoreProjection], bool]:
    """Build a mapped core on base coordinates without a SMILES round-trip."""
    from .reaction_core import build_reaction_core_projection

    atom_map_overrides = _mapping_in_original_coordinates(base, mapped)
    if atom_map_overrides is None or mapping.normalization is None:
        return None, False
    core = build_reaction_core_projection(
        reactants=base.reactants,
        products=base.products,
        edits=mapping.normalization.edits,
        stereo_changes=mapping.normalization.stereo_changes,
        evidence=evidence,
        confidence=confidence,
        atom_map_overrides=atom_map_overrides,
    )
    return core, True


def analyze_reaction_with_external_mapping(
    reaction_smiles: str,
    provider: AtomMappingProvider,
    *,
    base_analysis: Optional[ReactionAnalysis] = None,
    force_resolved_shadow: bool = False,
) -> ExternalMappingAssessment:
    """Use an external mapper only as configured, preserving conflict evidence."""
    from .reaction_api import featurize_reaction

    base = base_analysis or featurize_reaction(reaction_smiles)
    if not base.valid:
        return ExternalMappingAssessment(
            input_reaction_smiles=reaction_smiles,
            status="not_requested_invalid_reaction",
            analysis=base,
            provider_metadata=provider.metadata,
        )
    if any(
        component.atom_mapped
        for component in (*base.reactants, *base.products)
    ):
        return ExternalMappingAssessment(
            input_reaction_smiles=reaction_smiles,
            status="not_requested_supplied_mapping",
            analysis=base,
            provider_metadata=provider.metadata,
        )
    if base.reaction_signature is not None and not force_resolved_shadow:
        return ExternalMappingAssessment(
            input_reaction_smiles=reaction_smiles,
            status="not_requested_resolved_internal_evidence",
            analysis=base,
            provider_metadata=provider.metadata,
        )
    mapping = provider.map_reactions((reaction_smiles,))[0]
    external_candidate = _external_candidate(mapping)
    if not mapping.valid or mapping.normalization is None:
        warnings = tuple(
            sorted(set(base.warnings).union(mapping.warnings, ("EXTERNAL_MAPPING_FAILED",)))
        )
        return ExternalMappingAssessment(
            input_reaction_smiles=reaction_smiles,
            status="external_mapping_failed",
            analysis=replace(
                base,
                evidence_candidates=_merge_evidence_candidates(
                    base.evidence_candidates,
                    (external_candidate,),
                ),
                warnings=warnings,
            ),
            provider_metadata=provider.metadata,
            mapping_result=mapping,
            warnings=("EXTERNAL_MAPPING_FAILED",),
        )

    mapped_parsed = parse_reaction_smiles(str(mapping.mapped_reaction_smiles))
    external_profile = normalized_edit_profile(
        mapping.normalization.edits,
        mapped_parsed.reactants,
        mapped_parsed.products,
    )
    signature_matches = bool(
        base.reaction_signature is not None
        and normalized_edit_profile(
            base.reaction_signature.edits,
            base.reactants,
            base.products,
        )
        == external_profile
    )
    matches = tuple(
        hypothesis
        for hypothesis in base.edit_hypotheses
        if normalized_edit_profile(
            hypothesis.edits,
            base.reactants,
            base.products,
        )
        == external_profile
    )
    matched_ids = tuple(hypothesis.hypothesis_id for hypothesis in matches)
    if force_resolved_shadow and base.reaction_signature is not None:
        if not signature_matches:
            warning = "EXTERNAL_MAPPING_SIGNATURE_CONFLICT"
            reaction_core, _ = _mapped_core_proposal(
                base=base,
                mapped=mapped_parsed,
                mapping=mapping,
                evidence="conflicting_edit_evidence",
                confidence=mapping.normalization.confidence,
            )
            proposal_warning = "REACTION_CORE_CONFLICTING_EVIDENCE_PROPOSAL"
            if reaction_core is not None:
                reaction_core = replace(
                    reaction_core,
                    warnings=tuple(
                        sorted(
                            set(reaction_core.warnings).union(
                                (proposal_warning,)
                            )
                        )
                    ),
                )
            return ExternalMappingAssessment(
                input_reaction_smiles=reaction_smiles,
                status="external_mapping_signature_conflict",
                analysis=replace(
                    base,
                    evidence_candidates=_merge_evidence_candidates(
                        base.evidence_candidates,
                        (external_candidate,),
                    ),
                    reaction_core=reaction_core,
                    warnings=tuple(
                        sorted(
                            set(base.warnings).union(
                                (warning,)
                                + (
                                    (proposal_warning,)
                                    if reaction_core is not None
                                    else ()
                                )
                            )
                        )
                    ),
                ),
                provider_metadata=provider.metadata,
                mapping_result=mapping,
                warnings=(warning, "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW"),
            )

        # Forced mapping of an already resolved reaction is a shadow operation:
        # it supplies atom correspondence for the minimized graphic without
        # replacing the internally resolved identity or interpretation. The
        # display label is rerendered so compatible core context is not lost.
        confidence = min(
            mapping.normalization.confidence,
            min(
                (edit.confidence for edit in base.reaction_signature.edits),
                default=1.0,
            ),
        )
        evidence = "external_mapping_internal_consensus"
        assessment_warnings = (
            "EXTERNAL_MAPPING_INTERNAL_CONSENSUS",
            "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW",
        )
        reaction_core, coordinates_projected = _mapped_core_proposal(
            base=base,
            mapped=mapped_parsed,
            mapping=mapping,
            evidence=evidence,
            confidence=confidence,
        )
        if not coordinates_projected:
            warning = "EXTERNAL_MAPPING_CORE_COORDINATE_PROJECTION_FAILED"
            return ExternalMappingAssessment(
                input_reaction_smiles=reaction_smiles,
                status="external_mapping_signature_unavailable",
                analysis=replace(
                    base,
                    evidence_candidates=_merge_evidence_candidates(
                        base.evidence_candidates,
                        (external_candidate,),
                    ),
                    warnings=tuple(
                        sorted(set(base.warnings).union((warning,)))
                    ),
                ),
                provider_metadata=provider.metadata,
                mapping_result=mapping,
                warnings=(warning, "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW"),
            )
        if reaction_core is None:
            warning = "EXTERNAL_MAPPING_SIGNATURE_UNAVAILABLE"
            return ExternalMappingAssessment(
                input_reaction_smiles=reaction_smiles,
                status="external_mapping_signature_unavailable",
                analysis=replace(
                    base,
                    evidence_candidates=_merge_evidence_candidates(
                        base.evidence_candidates,
                        (external_candidate,),
                    ),
                    warnings=tuple(
                        sorted(set(base.warnings).union((warning,)))
                    ),
                ),
                provider_metadata=provider.metadata,
                mapping_result=mapping,
                warnings=(warning, "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW"),
            )
        effective = _attach_core_and_rerender(
            replace(
                base,
                evidence_candidates=_merge_evidence_candidates(
                    base.evidence_candidates,
                    (external_candidate,),
                ),
            ),
            reaction_core,
            warnings=assessment_warnings,
        )
        return ExternalMappingAssessment(
            input_reaction_smiles=reaction_smiles,
            status="external_mapping_internal_consensus",
            analysis=effective,
            provider_metadata=provider.metadata,
            mapping_result=mapping,
            warnings=assessment_warnings,
        )

    if base.edit_hypotheses and len(matches) != 1:
        status: ExternalMappingAssessmentStatus = (
            "external_mapping_ambiguous_hypothesis_match"
            if matches
            else "external_mapping_hypothesis_conflict"
        )
        warning = (
            "EXTERNAL_MAPPING_AMBIGUOUS_HYPOTHESIS_MATCH"
            if matches
            else "EXTERNAL_MAPPING_HYPOTHESIS_CONFLICT"
        )
        reaction_core, _ = _mapped_core_proposal(
            base=base,
            mapped=mapped_parsed,
            mapping=mapping,
            evidence="conflicting_edit_evidence",
            confidence=mapping.normalization.confidence,
        )
        proposal_warning = "REACTION_CORE_CONFLICTING_EVIDENCE_PROPOSAL"
        if reaction_core is not None:
            reaction_core = replace(
                reaction_core,
                warnings=tuple(
                    sorted(
                        set(reaction_core.warnings).union((proposal_warning,))
                    )
                ),
            )
        return ExternalMappingAssessment(
            input_reaction_smiles=reaction_smiles,
            status=status,
            analysis=replace(
                base,
                evidence_candidates=_merge_evidence_candidates(
                    base.evidence_candidates,
                    (external_candidate,),
                ),
                reaction_core=reaction_core,
                warnings=tuple(
                    sorted(
                        set(base.warnings).union(
                            (warning,)
                            + (
                                (proposal_warning,)
                                if reaction_core is not None
                                else ()
                            )
                        )
                    )
                ),
            ),
            provider_metadata=provider.metadata,
            mapping_result=mapping,
            matched_hypothesis_ids=matched_ids,
            warnings=(warning, "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW"),
        )

    if matches:
        status = "external_mapping_internal_consensus"
        evidence = "external_mapping_internal_consensus"
        confidence = min(
            mapping.normalization.confidence,
            matches[0].confidence,
        )
        assessment_warnings = (
            "EXTERNAL_MAPPING_INTERNAL_CONSENSUS",
            "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW",
        )
    else:
        status = "external_mapping_only"
        evidence = EXTERNAL_MAPPING_EVIDENCE
        confidence = mapping.normalization.confidence
        assessment_warnings = (
            "EXTERNAL_MAPPING_WITHOUT_INTERNAL_CONSENSUS",
            "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW",
        )
    override = replace(
        mapping.normalization,
        evidence=evidence,
        confidence=confidence,
        warnings=tuple(
            sorted(
                set(mapping.normalization.warnings).union(assessment_warnings)
            )
        ),
    )
    mapped_analysis = featurize_reaction(
        str(mapping.mapped_reaction_smiles),
        _mapped_edit_override=override,
        _mapped_provider=provider.metadata.provider_id,
    )
    if mapped_analysis.reaction_signature is None:
        warning = "EXTERNAL_MAPPING_SIGNATURE_UNAVAILABLE"
        return ExternalMappingAssessment(
            input_reaction_smiles=reaction_smiles,
            status="external_mapping_signature_unavailable",
            analysis=(
                _attach_core_and_rerender(
                    replace(
                        base,
                        evidence_candidates=_merge_evidence_candidates(
                            base.evidence_candidates,
                            (external_candidate,),
                        ),
                    ),
                    mapped_analysis.reaction_core,
                    warnings=(
                        warning,
                        "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW",
                    ),
                )
                if mapped_analysis.reaction_core is not None
                else replace(
                    base,
                    evidence_candidates=_merge_evidence_candidates(
                        base.evidence_candidates,
                        (external_candidate,),
                    ),
                    warnings=tuple(
                        sorted(set(base.warnings).union((warning,)))
                    ),
                )
            ),
            provider_metadata=provider.metadata,
            mapping_result=mapping,
            matched_hypothesis_ids=matched_ids,
            warnings=(warning, "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW"),
        )
    effective_warnings = tuple(
        sorted(set(mapped_analysis.warnings).union(assessment_warnings))
    )
    effective_signature = replace(
        mapped_analysis.reaction_signature,
        warnings=tuple(
            sorted(
                set(mapped_analysis.reaction_signature.warnings).union(
                    assessment_warnings
                )
            )
        ),
    )
    effective = replace(
        mapped_analysis,
        input_reaction_smiles=reaction_smiles,
        evidence_candidates=_merge_evidence_candidates(
            base.evidence_candidates,
            mapped_analysis.evidence_candidates,
            (external_candidate,),
        ),
        edit_hypotheses=base.edit_hypotheses,
        reaction_signature=effective_signature,
        warnings=effective_warnings,
    )
    return ExternalMappingAssessment(
        input_reaction_smiles=reaction_smiles,
        status=status,
        analysis=effective,
        provider_metadata=provider.metadata,
        mapping_result=mapping,
        matched_hypothesis_ids=matched_ids,
        warnings=assessment_warnings,
    )


class RxnMapperProvider:
    """Optional, instance-scoped RXNMapper adapter for offline conversion."""

    _MODEL_RELATIVE_PATH = (
        "models/transformers/albert_heads_8_uspto_all_1310k/"
        "pytorch_model.bin"
    )

    def __init__(
        self,
        *,
        batch_size: int = 16,
        mapper: Any = None,
    ) -> None:
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        self._batch_size = batch_size
        self._mapper = mapper
        try:
            provider_version = metadata.version("rxnmapper")
        except metadata.PackageNotFoundError:
            provider_version = "unavailable"
        model_sha256 = None
        if provider_version != "unavailable":
            try:
                model_path = resources.files("rxnmapper").joinpath(
                    self._MODEL_RELATIVE_PATH
                )
                model_sha256 = _sha256_file(model_path)
            except (FileNotFoundError, ModuleNotFoundError, OSError):
                model_sha256 = None
        self._metadata = AtomMappingProviderMetadata(
            provider_id="rxnmapper",
            provider_version=provider_version,
            model_id="albert_heads_8_uspto_all_1310k",
            model_sha256=model_sha256,
        )

    @staticmethod
    def is_available() -> bool:
        """Return whether the optional RXNMapper package can be imported."""
        return importlib.util.find_spec("rxnmapper") is not None

    @property
    def metadata(self) -> AtomMappingProviderMetadata:
        return self._metadata

    def _get_mapper(self) -> Any:
        if self._mapper is None:
            if not self.is_available():
                raise RuntimeError(
                    "RXNMapper is not installed; install requirements-mapping.txt"
                )
            from rxnmapper import RXNMapper

            self._mapper = RXNMapper()
        return self._mapper

    def _failure(
        self,
        reaction_smiles: str,
        error: str,
    ) -> ExternalAtomMappingResult:
        return ExternalAtomMappingResult(
            input_reaction_smiles=reaction_smiles,
            mapped_reaction_smiles=None,
            metadata=self.metadata,
            mapper_confidence=None,
            structure_preserved=False,
            reactant_mapping_coverage=0.0,
            product_mapping_coverage=0.0,
            shared_mapped_heavy_atom_count=0,
            normalization=None,
            valid=False,
            warnings=("EXTERNAL_MAPPER_FAILED",),
            error=error,
        )

    def _validate_raw(self, reaction: str, raw: Any) -> ExternalAtomMappingResult:
        if isinstance(raw, dict):
            mapped = raw.get("mapped_rxn") or raw.get("mapped_smiles")
            confidence = raw.get("confidence")
        else:
            mapped = str(raw) if raw is not None else None
            confidence = None
        if not mapped:
            return self._failure(reaction, "EMPTY_MAPPER_RESULT")
        numeric_confidence = (
            float(confidence) if isinstance(confidence, (int, float)) else None
        )
        return validate_external_atom_mapping(
            reaction,
            str(mapped),
            provider_metadata=self.metadata,
            mapper_confidence=numeric_confidence,
        )

    def map_reactions(
        self, reaction_smiles: Iterable[str]
    ) -> Tuple[ExternalAtomMappingResult, ...]:
        """Map batches, isolating a failed batch down to individual reactions."""
        reactions = tuple(str(reaction).strip() for reaction in reaction_smiles)
        if not reactions:
            return ()
        try:
            mapper = self._get_mapper()
        except Exception as exc:
            return tuple(self._failure(reaction, str(exc)) for reaction in reactions)
        outputs: list[ExternalAtomMappingResult] = []
        for offset in range(0, len(reactions), self._batch_size):
            batch = reactions[offset : offset + self._batch_size]
            try:
                raw_results = tuple(mapper.get_attention_guided_atom_maps(list(batch)))
                if len(raw_results) != len(batch):
                    raise RuntimeError("RXNMapper returned an unexpected batch size")
                outputs.extend(
                    self._validate_raw(reaction, raw)
                    for reaction, raw in zip(batch, raw_results)
                )
            except Exception as batch_exc:
                for reaction in batch:
                    try:
                        raw = mapper.get_attention_guided_atom_maps([reaction])[0]
                        outputs.append(self._validate_raw(reaction, raw))
                    except Exception as exc:
                        outputs.append(
                            self._failure(
                                reaction,
                                f"{type(batch_exc).__name__}: {exc}",
                            )
                        )
        return tuple(outputs)


__all__ = [
    "analyze_reaction_with_external_mapping",
    "AtomMappingProvider",
    "AtomMappingProviderMetadata",
    "EXTERNAL_ATOM_MAPPING_SCHEMA_VERSION",
    "EXTERNAL_MAPPING_EVIDENCE",
    "ExternalAtomMappingResult",
    "ExternalMappingAssessment",
    "ExternalMappingAssessmentStatus",
    "normalized_edit_profile",
    "RxnMapperProvider",
    "validate_external_atom_mapping",
]
