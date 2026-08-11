"""Dataset-driven condition-aware selectivity proof of concept.

The POC treats selectivity as a choice among structurally valid alternatives.
It currently supports the common single-connection replacement topology: one
heavy-atom bond is broken and one heavy-atom bond is formed through a shared
electrophilic endpoint.  Alternative connection endpoints are obtained from
``reactive_taxonomy`` and validated by applying the same graph edit.

No named reaction family or functional-group compatibility rule participates
in candidate generation or scoring.  A small hashed listwise linear model
learns site and condition preferences from observed choice sets.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache
from typing import Any, Dict, Iterable, Mapping, Sequence, Tuple

from condition_registry import CONDITION_RECIPE_COMPONENT_BUCKETS
from reactive_taxonomy import analyze_molecule, featurize_reaction
from reactive_taxonomy.reaction_graph_editing import bond_type, set_total_hydrogens


SELECTIVITY_POC_SCHEMA_VERSION = "1.0"
SELECTIVITY_POC_MODEL_ID = "conditional_edit_choice_poc.v1"
FUNCTIONAL_GROUP_COMPETITION_AUDIT_ID = (
    "structural_endpoint_competition_audit.v1"
)


def _digest(value: Any, *, length: int = 20) -> str:
    serialized = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:length]


def _reaction_sides(reaction_smiles: str) -> tuple[str, str]:
    parts = reaction_smiles.split(">")
    if len(parts) != 3 or not parts[0] or not parts[2]:
        raise ValueError("reaction SMILES must contain reactants>agents>products")
    return parts[0], parts[2]


def _condition_value_tokens(prefix: str, value: Any) -> list[str]:
    if value is None:
        return [f"{prefix}=<missing>"]
    if isinstance(value, Mapping):
        values = []
        for key in sorted(value, key=lambda item: str(item)):
            nested = f"{prefix}.{key}" if prefix else str(key)
            values.extend(_condition_value_tokens(nested, value[key]))
        return values
    if isinstance(value, (tuple, list, set, frozenset)):
        values = []
        for item in sorted(value, key=lambda member: str(member)):
            values.extend(_condition_value_tokens(prefix, item))
        return values
    if isinstance(value, float):
        normalized = format(value, ".8g")
    elif isinstance(value, bool):
        normalized = str(value).lower()
    else:
        normalized = str(value).strip().casefold()
    return [f"{prefix}={normalized}"]


def condition_tokens_from_mapping(conditions: Mapping[str, Any]) -> Tuple[str, ...]:
    """Flatten an arbitrary canonical recipe mapping into stable feature tokens."""

    return tuple(sorted(set(_condition_value_tokens("condition", conditions))))


def condition_tokens_from_recipe(recipe: Mapping[str, Any]) -> Tuple[str, ...]:
    """Project a resolved recipe to chemistry fields, excluding provenance IDs."""

    tokens = set()
    for bucket in CONDITION_RECIPE_COMPONENT_BUCKETS:
        for component in recipe.get(bucket) or ():
            identity = (
                component.get("substance_id")
                or component.get("cas")
                or component.get("canonical_name")
                or component.get("raw_identifier")
                or "<unresolved>"
            )
            normalized_identity = str(identity).strip().casefold()
            tokens.add(f"condition.{bucket}.identity={normalized_identity}")
            role = component.get("primary_role")
            if role:
                tokens.add(
                    f"condition.{bucket}.role={str(role).strip().casefold()}"
                )
                tokens.add(
                    f"condition.{bucket}.identity_role="
                    f"{normalized_identity}|{str(role).strip().casefold()}"
                )
    for field in (
        "temperature_c",
        "time_h",
        "concentration_m",
        "atmosphere",
    ):
        tokens.update(_condition_value_tokens(f"condition.{field}", recipe.get(field)))
    for stage in recipe.get("stages") or ():
        stage_index = int(stage.get("stage_index") or 0)
        for field in ("temperature_c", "time_h", "atmosphere"):
            tokens.update(
                _condition_value_tokens(
                    f"condition.stage{stage_index}.{field}",
                    stage.get(field),
                )
            )
    return tuple(sorted(tokens))


@dataclass(frozen=True)
class ReactionOutcomeCandidate:
    """One structurally validated endpoint choice for a fixed edit topology."""

    candidate_id: str
    component_index: int
    atom_index: int
    element: str
    site_type: str
    site_signature: str
    availability: str
    product_smiles: str
    feature_tokens: Tuple[str, ...]
    is_selected: bool
    schema_version: str = SELECTIVITY_POC_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-compatible candidate representation."""

        return asdict(self)


@dataclass(frozen=True)
class FunctionalGroupCompetitionOutcome:
    """Review-facing summary of one structurally valid endpoint outcome."""

    candidate_id: str
    component_index: int
    atom_index: int
    element: str
    site_type: str
    site_signature: str
    availability: str
    product_smiles: str


@dataclass(frozen=True)
class FunctionalGroupCompetitionWarning:
    """Non-ranking warning that the intended edit has alternative endpoints."""

    code: str
    message: str
    selected_outcome: FunctionalGroupCompetitionOutcome
    competing_outcomes: Tuple[FunctionalGroupCompetitionOutcome, ...]
    assessment_mode: str = "structural_endpoint_enumeration"
    ranking_impact: str = "none"
    conditions_evaluated: bool = False
    audit_id: str = FUNCTIONAL_GROUP_COMPETITION_AUDIT_ID
    schema_version: str = SELECTIVITY_POC_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-compatible warning representation."""

        return asdict(self)


@dataclass(frozen=True)
class ReactionChoiceSet:
    """Observed selection and all structurally available competing endpoints."""

    choice_set_id: str
    reaction_smiles: str
    precursor_smiles: str
    selected_product_smiles: str
    electrophile_component_index: int
    electrophile_atom_index: int
    leaving_group_component_index: int
    leaving_group_atom_index: int
    selected_candidate_id: str
    candidates: Tuple[ReactionOutcomeCandidate, ...]
    condition_tokens: Tuple[str, ...]
    reference_id: str
    label_strength: float
    warnings: Tuple[str, ...] = ()
    schema_version: str = SELECTIVITY_POC_SCHEMA_VERSION

    def __post_init__(self) -> None:
        candidate_ids = {candidate.candidate_id for candidate in self.candidates}
        if len(candidate_ids) != len(self.candidates):
            raise ValueError("choice-set candidate IDs must be unique")
        if self.selected_candidate_id not in candidate_ids:
            raise ValueError("selected candidate must belong to the choice set")
        if sum(candidate.is_selected for candidate in self.candidates) != 1:
            raise ValueError("choice set must contain exactly one selected candidate")
        if not 0.0 < self.label_strength <= 1.0:
            raise ValueError("label strength must be in (0, 1]")

    @property
    def selected_candidate(self) -> ReactionOutcomeCandidate:
        """Return the observed or intended candidate."""

        return next(
            candidate
            for candidate in self.candidates
            if candidate.candidate_id == self.selected_candidate_id
        )

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-compatible choice-set representation."""

        return asdict(self)


@dataclass(frozen=True)
class RankedReactionOutcome:
    """One model-ranked outcome from a selectivity choice set."""

    candidate_id: str
    site_signature: str
    element: str
    product_smiles: str
    probability: float
    is_selected: bool


@dataclass(frozen=True)
class SelectivityAssessment:
    """Condition-aware ranking, margin, uncertainty, and direct data support."""

    choice_set_id: str
    desired_candidate_id: str
    desired_probability: float
    best_competitor_probability: float
    probability_margin: float
    normalized_entropy: float
    exact_condition_reference_support: int
    evidence_status: str
    ranked_outcomes: Tuple[RankedReactionOutcome, ...]
    model_id: str = SELECTIVITY_POC_MODEL_ID
    schema_version: str = SELECTIVITY_POC_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-compatible assessment."""

        return asdict(self)


@dataclass(frozen=True)
class ChoiceModelTrainingReport:
    """Deterministic optimization summary for one fitted POC model."""

    choice_set_count: int
    candidate_count: int
    epoch_count: int
    initial_loss: float
    final_loss: float
    model_id: str = SELECTIVITY_POC_MODEL_ID


def _atom_key(reference: Any) -> tuple[int, int]:
    return int(reference.component_index), int(reference.atom_index)


def _single_connection_replacement(analysis: Any) -> tuple[Any, Any, Any]:
    observation = analysis.observation
    if observation is None or not observation.valid:
        raise ValueError("reaction does not have a valid graph observation")
    formed = tuple(
        edit
        for edit in observation.edits
        if edit.edit_type == "formed" and edit.atom_1 and edit.atom_2
    )
    broken = tuple(
        edit
        for edit in observation.edits
        if edit.edit_type == "broken" and edit.atom_1 and edit.atom_2
    )
    if len(formed) != 1 or len(broken) != 1:
        raise ValueError(
            "selectivity POC requires one formed and one broken heavy-atom bond"
        )
    formed_atoms = {_atom_key(formed[0].atom_1), _atom_key(formed[0].atom_2)}
    broken_atoms = {_atom_key(broken[0].atom_1), _atom_key(broken[0].atom_2)}
    shared = formed_atoms.intersection(broken_atoms)
    if len(shared) != 1:
        raise ValueError("formed and broken bonds must share one endpoint")
    electrophile_key = next(iter(shared))
    electrophile = next(
        atom
        for atom in (formed[0].atom_1, formed[0].atom_2)
        if _atom_key(atom) == electrophile_key
    )
    selected = next(
        atom
        for atom in (formed[0].atom_1, formed[0].atom_2)
        if _atom_key(atom) != electrophile_key
    )
    leaving_group = next(
        atom
        for atom in (broken[0].atom_1, broken[0].atom_2)
        if _atom_key(atom) != electrophile_key
    )
    return electrophile, selected, leaving_group


def _component_offsets(molecules: Sequence[Any]) -> tuple[int, ...]:
    offsets = []
    offset = 0
    for molecule in molecules:
        offsets.append(offset)
        offset += molecule.GetNumAtoms()
    return tuple(offsets)


def _candidate_product(
    component_smiles: Sequence[str],
    *,
    electrophile: Any,
    leaving_group: Any,
    candidate_component_index: int,
    candidate_endpoint: Any,
    formed_order: str,
) -> str | None:
    from rdkit import Chem

    if candidate_endpoint.required_bond_capacity_decrement:
        return None
    molecules = tuple(Chem.MolFromSmiles(smiles) for smiles in component_smiles)
    if any(molecule is None for molecule in molecules):
        return None
    typed_molecules = tuple(molecule for molecule in molecules if molecule is not None)
    offsets = _component_offsets(typed_molecules)
    combined = typed_molecules[0]
    for molecule in typed_molecules[1:]:
        combined = Chem.CombineMols(combined, molecule)
    editable = Chem.RWMol(combined)
    electrophile_index = (
        offsets[int(electrophile.component_index)] + int(electrophile.atom_index)
    )
    leaving_index = (
        offsets[int(leaving_group.component_index)] + int(leaving_group.atom_index)
    )
    candidate_index = (
        offsets[candidate_component_index]
        + int(candidate_endpoint.endpoint.atom_index)
    )
    if candidate_index in {electrophile_index, leaving_index}:
        return None
    if editable.GetBondBetweenAtoms(electrophile_index, leaving_index) is None:
        return None
    editable.RemoveBond(electrophile_index, leaving_index)
    fragments = Chem.GetMolFrags(
        editable.GetMol(),
        asMols=False,
        sanitizeFrags=False,
    )
    leaving_fragment = next(
        (set(fragment) for fragment in fragments if leaving_index in fragment),
        set(),
    )
    if not leaving_fragment or electrophile_index in leaving_fragment:
        return None
    if candidate_index in leaving_fragment:
        return None
    if not set_total_hydrogens(
        combined,
        editable,
        candidate_index,
        int(candidate_endpoint.required_hydrogen_delta),
    ):
        return None
    candidate_atom = editable.GetAtomWithIdx(candidate_index)
    candidate_atom.SetFormalCharge(
        int(candidate_atom.GetFormalCharge())
        + int(candidate_endpoint.required_formal_charge_delta)
    )
    if editable.GetBondBetweenAtoms(electrophile_index, candidate_index) is not None:
        return None
    editable.AddBond(
        electrophile_index,
        candidate_index,
        bond_type(formed_order),
    )
    for atom_index in sorted(leaving_fragment, reverse=True):
        editable.RemoveAtom(int(atom_index))
    product = editable.GetMol()
    try:
        Chem.SanitizeMol(product)
    except (ValueError, RuntimeError):
        return None
    for atom in product.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)


def _endpoint_feature_tokens(endpoint: Any) -> Tuple[str, ...]:
    atom = endpoint.endpoint
    values = {
        f"site_type={endpoint.source_site_type}",
        f"source_signature={endpoint.source_signature}",
        f"element={atom.element}",
        f"formal_charge={atom.formal_charge}",
        f"aromatic={atom.aromatic}",
        f"hybridization={atom.hybridization}",
        f"availability={endpoint.availability}",
        f"hydrogen_delta={endpoint.required_hydrogen_delta}",
        f"charge_delta={endpoint.required_formal_charge_delta}",
    }
    values.update(f"environment={value}" for value in atom.local_environment)
    values.update(f"context={value}" for value in atom.context_tokens)
    values.update(f"annotation={value}" for value in endpoint.annotation_tokens)
    return tuple(sorted(values))


def build_reaction_choice_set(
    reaction_smiles: str,
    conditions: Mapping[str, Any] | Iterable[str],
    *,
    reference_id: str = "",
    label_strength: float = 0.5,
) -> ReactionChoiceSet:
    """Enumerate validated endpoint alternatives for an observed reaction.

    The POC holds the reaction partners and single-connection replacement
    topology fixed.  It varies endpoints on the component containing the
    observed connection partner, which captures intramolecular functional-group
    competition without inventing additional reactant equivalents.
    """

    if not 0.0 < label_strength <= 1.0:
        raise ValueError("label strength must be in (0, 1]")
    reactant_side, product_side = _reaction_sides(reaction_smiles)
    component_smiles = tuple(reactant_side.split("."))
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid or analysis.error:
        raise ValueError(analysis.error or "reaction analysis failed")
    electrophile, selected, leaving_group = _single_connection_replacement(analysis)
    formed_edit = next(
        edit for edit in analysis.observation.edits if edit.edit_type == "formed"
    )
    observed_component_index = int(selected.component_index)
    molecule_analysis = analyze_molecule(component_smiles[observed_component_index])
    endpoint_by_atom: dict[int, list[Any]] = {}
    for interfaces in molecule_analysis.connectivity_hypotheses:
        for endpoint in interfaces.connection_endpoints:
            atom_index = endpoint.endpoint.atom_index
            if atom_index is None:
                continue
            endpoint_by_atom.setdefault(int(atom_index), []).append(endpoint)

    provisional = []
    for atom_index, endpoints in sorted(endpoint_by_atom.items()):
        endpoint = sorted(
            endpoints,
            key=lambda value: (
                value.source_site_type,
                value.source_signature,
                value.site_id,
            ),
        )[0]
        product_smiles = _candidate_product(
            component_smiles,
            electrophile=electrophile,
            leaving_group=leaving_group,
            candidate_component_index=observed_component_index,
            candidate_endpoint=endpoint,
            formed_order=str(formed_edit.new_order or "SINGLE"),
        )
        if product_smiles is None:
            continue
        combined_features = set()
        for alternative in endpoints:
            combined_features.update(_endpoint_feature_tokens(alternative))
        site_signature = "|".join(
            (
                str(endpoint.endpoint.element),
                str(endpoint.source_site_type),
                str(endpoint.source_signature),
                str(endpoint.endpoint.hybridization),
                str(endpoint.endpoint.formal_charge),
            )
        )
        candidate_id = "ROC1:" + _digest(
            (
                observed_component_index,
                atom_index,
                site_signature,
                product_smiles,
            )
        )
        is_selected = (
            observed_component_index == int(selected.component_index)
            and atom_index == int(selected.atom_index)
        )
        provisional.append(
            ReactionOutcomeCandidate(
                candidate_id=candidate_id,
                component_index=observed_component_index,
                atom_index=atom_index,
                element=str(endpoint.endpoint.element),
                site_type=str(endpoint.source_site_type),
                site_signature=site_signature,
                availability=str(endpoint.availability),
                product_smiles=product_smiles,
                feature_tokens=tuple(sorted(combined_features)),
                is_selected=is_selected,
            )
        )
    candidates = tuple(
        sorted(
            provisional,
            key=lambda candidate: (
                candidate.component_index,
                candidate.atom_index,
                candidate.candidate_id,
            ),
        )
    )
    selected_candidates = tuple(
        candidate for candidate in candidates if candidate.is_selected
    )
    if len(selected_candidates) != 1:
        raise ValueError("observed connection endpoint was not uniquely enumerated")
    if isinstance(conditions, Mapping):
        condition_tokens = condition_tokens_from_mapping(conditions)
    else:
        condition_tokens = tuple(sorted(set(str(value) for value in conditions)))
    selected_product = product_side
    warnings = tuple(analysis.warnings)
    identity = {
        "reaction": reaction_smiles,
        "conditions": condition_tokens,
        "reference": reference_id,
    }
    return ReactionChoiceSet(
        choice_set_id="RCS1:" + _digest(identity),
        reaction_smiles=reaction_smiles,
        precursor_smiles=reactant_side,
        selected_product_smiles=selected_product,
        electrophile_component_index=int(electrophile.component_index),
        electrophile_atom_index=int(electrophile.atom_index),
        leaving_group_component_index=int(leaving_group.component_index),
        leaving_group_atom_index=int(leaving_group.atom_index),
        selected_candidate_id=selected_candidates[0].candidate_id,
        candidates=candidates,
        condition_tokens=condition_tokens,
        reference_id=reference_id,
        label_strength=float(label_strength),
        warnings=warnings,
    )


def build_reaction_choice_set_from_record(
    record: Mapping[str, Any],
    *,
    label_strength: float = 0.5,
) -> ReactionChoiceSet:
    """Build one choice set from the canonical converted-record contract."""

    reaction_smiles = str(record.get("reaction_smiles") or "")
    if not reaction_smiles:
        raise ValueError("converted record does not contain reaction_smiles")
    recipe = record.get("resolved_recipe") or {}
    if not isinstance(recipe, Mapping):
        raise ValueError("converted record resolved_recipe must be a mapping")
    return build_reaction_choice_set(
        reaction_smiles,
        condition_tokens_from_recipe(recipe),
        reference_id=str(record.get("reference_id") or ""),
        label_strength=label_strength,
    )


def _competition_outcome(
    candidate: ReactionOutcomeCandidate,
) -> FunctionalGroupCompetitionOutcome:
    return FunctionalGroupCompetitionOutcome(
        candidate_id=candidate.candidate_id,
        component_index=candidate.component_index,
        atom_index=candidate.atom_index,
        element=candidate.element,
        site_type=candidate.site_type,
        site_signature=candidate.site_signature,
        availability=candidate.availability,
        product_smiles=candidate.product_smiles,
    )


@lru_cache(maxsize=20_000)
def detect_functional_group_competition(
    reaction_smiles: str,
) -> FunctionalGroupCompetitionWarning | None:
    """Return a review-only warning for distinct alternative endpoint products.

    Unsupported edit topologies do not produce a warning. This audit uses only
    graph-derived endpoint availability; it deliberately does not estimate
    condition-dependent selectivity and cannot affect retrosynthesis ranking.
    """

    try:
        choice_set = build_reaction_choice_set(reaction_smiles, ())
    except ValueError:
        return None
    selected = choice_set.selected_candidate
    competitors = tuple(
        candidate
        for candidate in choice_set.candidates
        if not candidate.is_selected
        and candidate.product_smiles != selected.product_smiles
    )
    if not competitors:
        return None
    competitor_elements = ", ".join(
        sorted({candidate.element for candidate in competitors})
    )
    count = len(competitors)
    message = (
        "Possible functional-group competition: the selected "
        f"{selected.element} endpoint has {count} structurally valid "
        f"alternative endpoint{'s' if count != 1 else ''} "
        f"({competitor_elements}) on the same reactant component. "
        "This is a screening warning only; reaction conditions and relative "
        "outcome probabilities have not been evaluated."
    )
    return FunctionalGroupCompetitionWarning(
        code="POSSIBLE_FUNCTIONAL_GROUP_COMPETITION",
        message=message,
        selected_outcome=_competition_outcome(selected),
        competing_outcomes=tuple(
            _competition_outcome(candidate) for candidate in competitors
        ),
    )


def _hashed_index(token: str, feature_dimension: int) -> tuple[int, float]:
    digest = hashlib.blake2b(token.encode("utf-8"), digest_size=8).digest()
    value = int.from_bytes(digest, byteorder="big", signed=False)
    index = value % feature_dimension
    sign = 1.0 if (value >> 63) == 0 else -1.0
    return index, sign


class ConditionalEditChoiceModel:
    """Small deterministic listwise model over generic edit choices."""

    def __init__(
        self,
        *,
        feature_dimension: int = 512,
        l2_penalty: float = 1e-4,
    ) -> None:
        if feature_dimension < 32:
            raise ValueError("feature dimension must be at least 32")
        if l2_penalty < 0.0:
            raise ValueError("L2 penalty cannot be negative")
        self.feature_dimension = int(feature_dimension)
        self.l2_penalty = float(l2_penalty)
        self.weights = [0.0] * self.feature_dimension
        self.training_choice_sets: Tuple[ReactionChoiceSet, ...] = ()

    def _features(
        self,
        choice_set: ReactionChoiceSet,
        candidate: ReactionOutcomeCandidate,
    ) -> Dict[int, float]:
        candidate_tokens = set(candidate.feature_tokens)
        candidate_tokens.add(f"site_signature={candidate.site_signature}")
        for competitor in choice_set.candidates:
            if competitor.candidate_id == candidate.candidate_id:
                continue
            candidate_tokens.add(
                "competition="
                f"{candidate.site_signature}>{competitor.site_signature}"
            )
        tokens = set(candidate_tokens)
        for condition in choice_set.condition_tokens:
            for candidate_token in candidate_tokens:
                tokens.add(f"{condition}|{candidate_token}")
        features: Dict[int, float] = {}
        for token in sorted(tokens):
            index, sign = _hashed_index(token, self.feature_dimension)
            features[index] = features.get(index, 0.0) + sign
        return features

    def _score(self, features: Mapping[int, float]) -> float:
        return sum(self.weights[index] * value for index, value in features.items())

    @staticmethod
    def _softmax(scores: Sequence[float]) -> tuple[float, ...]:
        maximum = max(scores)
        exponentials = tuple(math.exp(score - maximum) for score in scores)
        denominator = sum(exponentials)
        return tuple(value / denominator for value in exponentials)

    def predict_probabilities(
        self, choice_set: ReactionChoiceSet
    ) -> Dict[str, float]:
        """Return normalized probabilities for every available edit choice."""

        features = tuple(
            self._features(choice_set, candidate)
            for candidate in choice_set.candidates
        )
        probabilities = self._softmax(tuple(self._score(value) for value in features))
        return {
            candidate.candidate_id: probability
            for candidate, probability in zip(
                choice_set.candidates,
                probabilities,
                strict=True,
            )
        }

    def _loss(self, choice_sets: Sequence[ReactionChoiceSet]) -> float:
        total_weight = sum(choice_set.label_strength for choice_set in choice_sets)
        if total_weight <= 0.0:
            return 0.0
        value = 0.0
        for choice_set in choice_sets:
            probability = self.predict_probabilities(choice_set)[
                choice_set.selected_candidate_id
            ]
            value -= choice_set.label_strength * math.log(max(probability, 1e-15))
        value /= total_weight
        value += 0.5 * self.l2_penalty * sum(
            weight * weight for weight in self.weights
        )
        return value

    def fit(
        self,
        choice_sets: Iterable[ReactionChoiceSet],
        *,
        epochs: int = 250,
        learning_rate: float = 0.15,
    ) -> ChoiceModelTrainingReport:
        """Fit a deterministic full-batch weighted listwise softmax model."""

        values = tuple(sorted(choice_sets, key=lambda value: value.choice_set_id))
        if not values:
            raise ValueError("at least one choice set is required")
        if any(len(value.candidates) < 2 for value in values):
            raise ValueError("training choice sets require at least two candidates")
        if epochs < 1 or learning_rate <= 0.0:
            raise ValueError("epochs and learning rate must be positive")
        self.weights = [0.0] * self.feature_dimension
        initial_loss = self._loss(values)
        total_weight = sum(value.label_strength for value in values)
        feature_cache = {
            (choice_set.choice_set_id, candidate.candidate_id): self._features(
                choice_set, candidate
            )
            for choice_set in values
            for candidate in choice_set.candidates
        }
        for _ in range(epochs):
            gradient = [0.0] * self.feature_dimension
            for choice_set in values:
                candidate_features = tuple(
                    feature_cache[(choice_set.choice_set_id, candidate.candidate_id)]
                    for candidate in choice_set.candidates
                )
                probabilities = self._softmax(
                    tuple(self._score(features) for features in candidate_features)
                )
                for candidate, probability, features in zip(
                    choice_set.candidates,
                    probabilities,
                    candidate_features,
                    strict=True,
                ):
                    target = float(
                        candidate.candidate_id == choice_set.selected_candidate_id
                    )
                    factor = choice_set.label_strength * (probability - target)
                    for index, feature_value in features.items():
                        gradient[index] += factor * feature_value
            for index, weight in enumerate(self.weights):
                normalized = gradient[index] / total_weight
                normalized += self.l2_penalty * weight
                self.weights[index] -= learning_rate * normalized
        self.training_choice_sets = values
        return ChoiceModelTrainingReport(
            choice_set_count=len(values),
            candidate_count=sum(len(value.candidates) for value in values),
            epoch_count=epochs,
            initial_loss=round(initial_loss, 8),
            final_loss=round(self._loss(values), 8),
        )

    def assess(self, choice_set: ReactionChoiceSet) -> SelectivityAssessment:
        """Assess the intended outcome against all enumerated competitors."""

        probabilities = self.predict_probabilities(choice_set)
        ranked = tuple(
            RankedReactionOutcome(
                candidate_id=candidate.candidate_id,
                site_signature=candidate.site_signature,
                element=candidate.element,
                product_smiles=candidate.product_smiles,
                probability=round(probabilities[candidate.candidate_id], 8),
                is_selected=candidate.is_selected,
            )
            for candidate in sorted(
                choice_set.candidates,
                key=lambda value: (
                    -probabilities[value.candidate_id],
                    value.candidate_id,
                ),
            )
        )
        desired_probability = probabilities[choice_set.selected_candidate_id]
        competitor_probabilities = tuple(
            probability
            for candidate_id, probability in probabilities.items()
            if candidate_id != choice_set.selected_candidate_id
        )
        best_competitor = max(competitor_probabilities, default=0.0)
        entropy = -sum(
            probability * math.log(max(probability, 1e-15))
            for probability in probabilities.values()
        )
        normalized_entropy = (
            entropy / math.log(len(probabilities)) if len(probabilities) > 1 else 0.0
        )
        selected_signature = choice_set.selected_candidate.site_signature
        references = {
            training.reference_id
            for training in self.training_choice_sets
            if training.reference_id
            and training.condition_tokens == choice_set.condition_tokens
            and training.selected_candidate.site_signature == selected_signature
        }
        if references:
            evidence_status = "direct_condition_support"
        elif self.training_choice_sets:
            evidence_status = "model_only"
        else:
            evidence_status = "unfitted"
        return SelectivityAssessment(
            choice_set_id=choice_set.choice_set_id,
            desired_candidate_id=choice_set.selected_candidate_id,
            desired_probability=round(desired_probability, 8),
            best_competitor_probability=round(best_competitor, 8),
            probability_margin=round(desired_probability - best_competitor, 8),
            normalized_entropy=round(normalized_entropy, 8),
            exact_condition_reference_support=len(references),
            evidence_status=evidence_status,
            ranked_outcomes=ranked,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Serialize learned parameters without retaining training records."""

        return {
            "model_id": SELECTIVITY_POC_MODEL_ID,
            "schema_version": SELECTIVITY_POC_SCHEMA_VERSION,
            "feature_dimension": self.feature_dimension,
            "l2_penalty": self.l2_penalty,
            "weights": tuple(self.weights),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ConditionalEditChoiceModel":
        """Load a serialized POC model with strict version validation."""

        if value.get("model_id") != SELECTIVITY_POC_MODEL_ID:
            raise ValueError("unsupported selectivity model ID")
        if value.get("schema_version") != SELECTIVITY_POC_SCHEMA_VERSION:
            raise ValueError("unsupported selectivity model schema")
        model = cls(
            feature_dimension=int(value["feature_dimension"]),
            l2_penalty=float(value["l2_penalty"]),
        )
        weights = tuple(float(weight) for weight in value.get("weights") or ())
        if len(weights) != model.feature_dimension:
            raise ValueError("serialized weight dimension does not match the model")
        model.weights = list(weights)
        return model


__all__ = [
    "ChoiceModelTrainingReport",
    "ConditionalEditChoiceModel",
    "RankedReactionOutcome",
    "ReactionChoiceSet",
    "ReactionOutcomeCandidate",
    "SELECTIVITY_POC_MODEL_ID",
    "SELECTIVITY_POC_SCHEMA_VERSION",
    "SelectivityAssessment",
    "build_reaction_choice_set",
    "build_reaction_choice_set_from_record",
    "condition_tokens_from_mapping",
    "condition_tokens_from_recipe",
]
