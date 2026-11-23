"""Compare retained functional groups on edit-aligned precursor partners."""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import asdict, dataclass, field
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Iterator, get_args

from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator

from .api import component_annotation_atom_order
from .chemistry.rdkit_utils import parse_smiles
from .chemistry.smarts_cache import compile_smarts
from .molecular_motifs import load_molecular_motif_definitions
from .models import SiteType
from .reaction_api import featurize_reaction
from .reaction_models import ReactionAnalysis, ReactionComponent
from .reaction_signatures import reaction_signature_definition_versions


_DEFINITIONS = Path(__file__).with_name("definitions")


@dataclass(frozen=True)
class _AdditionalContext:
    feature: str
    smarts: str
    atom_positions: tuple[int, ...]
    supersedes: tuple[str, ...]


@dataclass(frozen=True)
class _Rules:
    version: str
    additions: tuple[_AdditionalContext, ...]
    departures: tuple[tuple[str, str], ...]
    element_motifs: tuple[tuple[str, tuple[str, ...]], ...]
    competing_site_types: tuple[str, ...]
    fingerprint_radius: int
    fingerprint_bits: int
    use_chirality: bool
    max_exact_heavy_bonds: int
    positional_exclusions: tuple[str, ...]


@lru_cache(maxsize=1)
def _rules() -> _Rules:
    path = _DEFINITIONS / "step_precedent_context.v1.json"
    data = json.loads(path.read_bytes())
    if (
        data["definition_id"] != "step_precedent_context.v1"
        or data["schema_version"] != "1.0"
        or not data["definition_version"]
        or data["motif_source"] != "molecular_motifs.v1.json"
    ):
        raise ValueError("Invalid step context definition")
    motifs = {d["id"] for d in load_molecular_motif_definitions()}
    ids = set(motifs)
    additions = []
    for row in data["additional_contexts"]:
        feature = row["id"]
        query = compile_smarts(row["smarts"], validate=True)
        maps = {
            a.GetAtomMapNum(): a.GetIdx() for a in query.GetAtoms() if a.GetAtomMapNum()
        }
        atoms = row["context_atom_maps"]
        if (
            not feature
            or feature in ids
            or not atoms
            or len(set(atoms)) != len(atoms)
            or not set(atoms) <= maps.keys()
            or not set(row["supersedes_at_same_atoms"]) <= motifs
        ):
            raise ValueError("Invalid additional step context")
        ids.add(feature)
        additions.append(
            _AdditionalContext(
                feature,
                row["smarts"],
                tuple(maps[i] for i in atoms),
                tuple(row["supersedes_at_same_atoms"]),
            )
        )
    elements = data["element_resolved_motifs"]
    if not elements.keys() <= motifs or not all(elements.values()):
        raise ValueError("Unknown or empty element-specific motif")
    departures = [
        (el, group)
        for group, members in data["departure_equivalences"].items()
        for el in members
    ]
    if len(dict(departures)) != len(departures):
        raise ValueError("Overlapping departing element classes")
    sites = tuple(data["competing_site_types"])
    if (
        not sites
        or len(sites) != len(set(sites))
        or not set(sites) <= set(get_args(SiteType))
    ):
        raise ValueError("Invalid competing site type")
    fp = data["fingerprint"]
    if (
        type(fp["radius"]) is not int
        or fp["radius"] < 1
        or type(fp["bits"]) is not int
        or fp["bits"] < 1
        or type(fp["use_chirality"]) is not bool
    ):
        raise ValueError("Invalid step context fingerprint definition")
    positions = data["positional_context"]
    radius = positions["max_exact_heavy_bonds"]
    excluded = tuple(positions["exclude_composite_features"])
    if type(radius) is not int or not 1 <= radius <= 6 or not set(excluded) <= ids:
        raise ValueError("Invalid step context positional definition")
    basis = json.dumps(
        [reaction_signature_definition_versions(), rdBase.rdkitVersion], sort_keys=True
    ).encode()
    digest = hashlib.sha256(
        path.read_bytes() + (_DEFINITIONS / data["motif_source"]).read_bytes() + basis
    ).hexdigest()
    return _Rules(
        f"{data['definition_id']}@{data['definition_version']}:{digest[:16]}",
        tuple(additions),
        tuple(departures),
        tuple((k, tuple(v)) for k, v in elements.items()),
        sites,
        fp["radius"],
        fp["bits"],
        fp["use_chirality"],
        radius,
        excluded,
    )


@dataclass(frozen=True, order=True)
class StepContextPartner:
    """A graph-contributing reactant and its retained or uncertain tracked features."""

    edit_anchor: str
    precursor_smiles: str
    retained_features: tuple[str, ...]
    uncertain_features: tuple[str, ...]


@dataclass(frozen=True)
class StepPrecedentContext:
    """A versioned, precursor-order-invariant inventory on the reacting partners."""

    status: str
    partners: tuple[StepContextPartner, ...]
    requirements: tuple[str, ...]
    warnings: tuple[str, ...]
    definition_version: str
    has_specified_stereo: bool = False
    schema_version: str = "step_precedent_context.v1"
    product_element_excess: tuple[tuple[str, int], ...] = ()

    def to_dict(self) -> dict:
        """Serialize the molecular evidence and version."""
        return asdict(self)


@dataclass(frozen=True)
class StepContextComparison:
    """Query coverage and extra source chemistry on the aligned precursor partners."""

    status: str
    covered: tuple[str, ...]
    missing: tuple[str, ...]
    unresolved: tuple[str, ...]
    minimum_partner_similarity: float | None
    warnings: tuple[str, ...]
    definition_version: str
    reverse_status: str
    source_only: tuple[str, ...]
    source_unresolved: tuple[str, ...]
    schema_version: str = "step_context_comparison.v1"

    def to_dict(self) -> dict:
        """Serialize the alignment result and unsatisfied requirements."""
        return asdict(self)


def _component_anchors(analysis: ReactionAnalysis) -> dict[int, dict[int, str]]:
    retained = {
        (t.before_state.component_index, t.before_state.atom_index)
        for t in analysis.reaction_core.atom_transitions
        if t.before_state and t.after_state
    }
    departures = dict(_rules().departures)
    tokens: dict[tuple[int, int], list[str]] = defaultdict(list)
    labels = {}
    for edit in analysis.observation.edits:
        for atom, other in ((edit.atom_1, edit.atom_2), (edit.atom_2, edit.atom_1)):
            if atom is None or atom.side != "reactant":
                continue
            key = (atom.component_index, atom.atom_index)
            if key not in retained:
                continue
            labels[key] = ("aromatic-" if atom.aromatic else "") + atom.element
            other_element = other.element if other else "H"
            if (
                edit.edit_type == "broken"
                and other
                and (other.component_index, other.atom_index) not in retained
            ):
                other_element = departures.get(other_element, other_element)
            tokens[key].append(
                f"{edit.edit_type}({other_element},{edit.old_order or '-'}>{edit.new_order or '-'})"
            )
    components: dict[int, dict[int, str]] = defaultdict(dict)
    for (component, atom), changes in tokens.items():
        components[component][atom] = (
            f"{labels[component, atom]}[{','.join(sorted(changes))}]"
        )
    return components


def _feature_occurrences(
    component: ReactionComponent, mol: Chem.Mol
) -> set[tuple[str, frozenset[int]]]:
    # A noncanonical SMILES traversal can still reorder branches. Molecular
    # annotations refer to that traversal; edits and continuity use the input mol.
    annotation_to_input = component_annotation_atom_order(mol)
    observed = set()
    element_motifs = dict(_rules().element_motifs)
    for motif in component.molecule_analysis.interpretation.motifs:
        feature = motif.motif_id
        atoms = frozenset(annotation_to_input[i] for i in motif.atom_indices)
        allowed = element_motifs.get(feature)
        if allowed:
            actual = sorted(
                mol.GetAtomWithIdx(i).GetSymbol()
                for i in atoms
                if mol.GetAtomWithIdx(i).GetSymbol() in allowed
            )
            feature = f"{feature}({','.join(actual)})"
        observed.add((feature, atoms))
    motif_atoms = {atoms for _, atoms in observed}
    for site in component.molecule_analysis.reactive_site_hypotheses:
        if site.site_type in _rules().competing_site_types:
            atoms = frozenset(annotation_to_input[i] for i in site.atom_indices)
            if atoms not in motif_atoms:
                observed.add((site.canonical_signature, atoms))
    for row in _rules().additions:
        pattern = compile_smarts(row.smarts, validate=False)
        for match in mol.GetSubstructMatches(pattern):
            atoms = frozenset(match[i] for i in row.atom_positions)
            for masked in row.supersedes:
                observed.discard((masked, atoms))
            observed.add((row.feature, atoms))
    return observed


def _nearby_distances(
    mol: Chem.Mol, origins: frozenset[int], radius: int
) -> dict[int, int]:
    """Minimum heavy-bond distance from the feature, stopping at the configured radius."""
    distances = {i: 0 for i in origins}
    pending = deque(origins)
    while pending:
        i = pending.popleft()
        distance = distances[i]
        if distance == radius:
            continue
        for atom in mol.GetAtomWithIdx(i).GetNeighbors():
            j = atom.GetIdx()
            if j not in distances and atom.GetAtomicNum() > 1:
                distances[j] = distance + 1
                pending.append(j)
    return distances


def _positional_features(
    mol: Chem.Mol,
    anchors: dict[int, str],
    occurrences: list[tuple[str, frozenset[int], bool]],
) -> Iterator[tuple[str, bool]]:
    """Compare reactant positions of surviving motifs, including their multiplicities."""
    rules = _rules()
    radius = rules.max_exact_heavy_bonds
    tracked = [o for o in occurrences if o[0] not in rules.positional_exclusions]
    at_atom: dict[int, list[int]] = defaultdict(list)
    for i, (_, atoms, _) in enumerate(tracked):
        for atom in atoms:
            at_atom[atom].append(i)
    for i, (feature, atoms, confirmed) in enumerate(tracked):
        distances = _nearby_distances(mol, atoms, radius)
        near = [
            (distance, anchors[a]) for a, distance in distances.items() if a in anchors
        ]
        if near:
            closest = min(distance for distance, _ in near)
            labels = sorted(label for distance, label in near if distance == closest)
            target = "edit" if len(anchors) == 1 else f"edit({'+'.join(labels)})"
            unit = "bond" if closest == 1 else "bonds"
            yield f"{feature} to {target} = {closest} {unit}", confirmed
        else:
            yield f"{feature} to edit > {radius} bonds", confirmed
        nearby: dict[int, int] = {}
        for atom, distance in distances.items():
            for j in at_atom.get(atom, ()):
                if j > i:
                    nearby[j] = min(distance, nearby.get(j, radius + 1))
        for j, distance in nearby.items():
            if distance:
                other, _, other_confirmed = tracked[j]
                a, b = sorted((feature, other))
                unit = "bond" if distance == 1 else "bonds"
                yield f"{a} to {b} = {distance} {unit}", confirmed and other_confirmed


def _partner(
    analysis: ReactionAnalysis, component: ReactionComponent, anchors: dict[int, str]
) -> StepContextPartner:
    mol = parse_smiles(component.input_smiles)
    centers = {
        t.before_state.atom_index
        for t in analysis.reaction_core.atom_transitions
        if t.before_state
        and t.before_state.component_index == component.component_index
    }
    remote = [
        r
        for r in analysis.reaction_core.remote_subgraphs
        if r.side == "reactant" and r.component_index == component.component_index
    ]
    retained = {i for r in remote if r.continuity == "retained" for i in r.atom_indices}
    departing = {
        i for r in remote if r.continuity == "departing" for i in r.atom_indices
    }
    occurrences = []
    for feature, atoms in _feature_occurrences(component, mol):
        if centers & atoms or atoms <= departing:
            continue
        occurrences.append((feature, atoms, atoms <= retained))
    feature_states = [(feature, known) for feature, _, known in occurrences]
    feature_states.extend(_positional_features(mol, anchors, occurrences))
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return StepContextPartner(
        "+".join(sorted(anchors.values())),
        Chem.MolToSmiles(mol),
        tuple(sorted(feature for feature, known in feature_states if known)),
        tuple(sorted(feature for feature, known in feature_states if not known)),
    )


def _numbered_requirements(
    partners: tuple[StepContextPartner, ...],
) -> list[tuple[StepContextPartner, str]]:
    counts = Counter(p.edit_anchor for p in partners)
    seen = Counter()
    result = []
    for partner in partners:
        anchor = partner.edit_anchor
        seen[anchor] += 1
        label = f"{anchor}#{seen[anchor]}" if counts[anchor] > 1 else anchor
        result.append((partner, label))
    return result


def _features(
    partner: StepContextPartner, label: str
) -> Iterator[tuple[str, str, int, bool]]:
    confirmed = Counter(partner.retained_features)
    all_features = confirmed + Counter(partner.uncertain_features)
    for feature, count in sorted(all_features.items()):
        for i in range(1, count + 1):
            requirement = f"{label} / {feature}" + (f" #{i}" if count > 1 else "")
            yield requirement, feature, i, i <= confirmed[feature]


def build_step_context(reaction_smiles: str) -> StepPrecedentContext:
    """Observe registered retained motifs on precursors contributing to the edited product.

    Group identities, counts, and short graph distances are partner-specific.
    Positional requirements summarize the reactant graph of retained motifs.
    They do not establish full scaffold or three-dimensional alignment.
    """
    return context_from_reaction_analysis(featurize_reaction(reaction_smiles))


def context_from_reaction_analysis(analysis: ReactionAnalysis) -> StepPrecedentContext:
    """Reuse a previously computed reaction analysis without remapping its inputs."""
    version = _rules().version
    if (
        not analysis.valid
        or not analysis.reaction_signature
        or not analysis.reaction_core
        or not analysis.observation
        or not analysis.reaction_completeness
        or analysis.reaction_completeness.status != "verified"
    ):
        completeness = analysis.reaction_completeness
        warnings = {"VERIFIED_EDIT_CONTEXT_UNAVAILABLE"}
        if analysis.error:
            warnings.add(analysis.error)
        if completeness:
            warnings.update(completeness.warnings)
        return StepPrecedentContext(
            "unresolved", (), (), tuple(sorted(warnings)), version,
            product_element_excess=(
                tuple(sorted(completeness.product_element_excess.items()))
                if completeness else ()
            ),
        )
    anchors = _component_anchors(analysis)
    partners = tuple(
        sorted(
            _partner(analysis, c, anchors[c.component_index])
            for c in analysis.reactants
            if c.component_index in anchors
        )
    )
    if not partners:
        return StepPrecedentContext(
            "unresolved", (), (), ("RETAINED_EDIT_ANCHORS_UNAVAILABLE",), version
        )
    stereo = any(
        _has_specified_stereo(c.input_smiles)
        for c in (*analysis.reactants, *analysis.products)
    )
    requirements = tuple(
        sorted(
            key
            for p, label in _numbered_requirements(partners)
            for key, _, _, _ in _features(p, label)
        )
    )
    warnings = (
        ("UNCERTAIN_FEATURE_RETENTION",)
        if any(p.uncertain_features for p in partners)
        else ()
    )
    return StepPrecedentContext(
        "resolved", partners, requirements, warnings, version, stereo
    )


def _has_specified_stereo(smiles: str) -> bool:
    mol = parse_smiles(smiles)
    return any(
        a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()
    ) or any(
        b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
        for b in mol.GetBonds()
    )


@lru_cache(maxsize=1024)
def _fingerprint(smiles: str) -> frozenset[int]:
    rules = _rules()
    fp = rdFingerprintGenerator.GetMorganGenerator(
        radius=rules.fingerprint_radius,
        fpSize=rules.fingerprint_bits,
        includeChirality=rules.use_chirality,
    )
    return frozenset(fp.GetFingerprint(parse_smiles(smiles)).GetOnBits())


@dataclass
class _Coverage:
    covered: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    unresolved: list[str] = field(default_factory=list)

    def compare(
        self, partner: StepContextPartner, label: str, counterpart: StepContextPartner
    ) -> None:
        confirmed = Counter(counterpart.retained_features)
        possible = confirmed + Counter(counterpart.uncertain_features)
        for key, feature, number, known in _features(partner, label):
            if not known or confirmed[feature] < number <= possible[feature]:
                self.unresolved.append(key)
            elif number <= confirmed[feature]:
                self.covered.append(key)
            else:
                self.missing.append(key)

    def status(self, ambiguous: bool) -> str:
        if ambiguous:
            return "unresolved"
        if self.missing:
            return "missing"
        return "unresolved" if self.unresolved else "covered"


def _groups(
    context: StepPrecedentContext,
) -> dict[str, list[tuple[StepContextPartner, str]]]:
    groups = defaultdict(list)
    for partner, label in _numbered_requirements(context.partners):
        groups[partner.edit_anchor].append((partner, label))
    return groups


def compare_step_contexts(
    query: StepPrecedentContext, source: StepPrecedentContext
) -> StepContextComparison:
    """Check each feature on its corresponding partner and retain ambiguous assignments.

    Alignment requires matching retained edit-anchor signatures. Identical signatures
    across multiple partners remain ambiguous unless the complete groups of partners
    have identical molecular identities and feature inventories.
    """
    version = _rules().version
    warnings = set(query.warnings + source.warnings)
    query_groups, source_groups = _groups(query), _groups(source)
    unaligned = None
    if (
        query.status != "resolved"
        or source.status != "resolved"
        or query.definition_version != version
        or source.definition_version != version
    ):
        unaligned = "CONTEXT_COMPARISON_UNAVAILABLE"
    elif query_groups.keys() != source_groups.keys():
        unaligned = "EDIT_ANCHOR_ROLES_DIFFER"
    if unaligned:
        return StepContextComparison(
            status="unresolved",
            covered=(),
            missing=(),
            unresolved=query.requirements,
            minimum_partner_similarity=None,
            warnings=tuple(sorted(warnings | {unaligned})),
            definition_version=version,
            reverse_status="unresolved",
            source_only=(),
            source_unresolved=source.requirements,
        )
    forward, reverse = _Coverage(), _Coverage()
    similarities = []
    ambiguous = False
    for role, queries in query_groups.items():
        sources = source_groups[role]
        if len(queries) != len(sources) or (len(queries) > 1 and queries != sources):
            ambiguous = True
            forward.unresolved.extend(
                key for p, label in queries for key, _, _, _ in _features(p, label)
            )
            reverse.unresolved.extend(
                key for s, label in sources for key, _, _, _ in _features(s, label)
            )
            continue
        for (p, p_label), (s, s_label) in zip(queries, sources):
            left, right = (
                _fingerprint(p.precursor_smiles),
                _fingerprint(s.precursor_smiles),
            )
            similarities.append(len(left & right) / len(left | right))
            forward.compare(p, p_label, s)
            reverse.compare(s, s_label, p)
    if ambiguous:
        warnings.add("AMBIGUOUS_PARTNER_ALIGNMENT")
    if not query.requirements:
        warnings.add("NO_TRACKED_SPECTATOR_REQUIREMENTS")
    if reverse.missing:
        warnings.add("SOURCE_ONLY_RETAINED_CONTEXT")
    if reverse.unresolved:
        warnings.add("SOURCE_CONTEXT_COMPARISON_UNRESOLVED")
    return StepContextComparison(
        status=forward.status(ambiguous),
        covered=tuple(sorted(forward.covered)),
        missing=tuple(sorted(forward.missing)),
        unresolved=tuple(sorted(forward.unresolved)),
        minimum_partner_similarity=None if ambiguous else min(similarities),
        warnings=tuple(sorted(warnings)),
        definition_version=version,
        reverse_status=reverse.status(ambiguous),
        source_only=tuple(sorted(reverse.missing)),
        source_unresolved=tuple(sorted(reverse.unresolved)),
    )
