"""Data-driven route-action policy learned from validated replay choices."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from functools import lru_cache
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping, Optional, Sequence

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from .chemistry import digest
from .generic_models import GenericDisconnectionCandidate
from .route_action_conversion import iter_route_action_evaluations
from .route_action_evaluation import (
    RouteActionCandidate,
    RouteActionEvaluation,
)


ROUTE_ACTION_POLICY_SCHEMA_VERSION = "1.1"
ROUTE_POLICY_EXAMPLE_SCHEMA_VERSION = "1.0"
_DEFINITION_PATH = (
    Path(__file__).resolve().parent
    / "definitions"
    / "route_action_policy.v1.json"
)


@dataclass(frozen=True)
class RouteActionPolicyDefinition:
    """Validated training and planner-integration configuration."""

    definition_id: str
    feature_dimension: int
    l2_penalty: float
    epochs: int
    learning_rate: float
    strategy_equivalent_label_strength: float
    planner_probability_deficit_weight: float
    baseline_rank_logit_weight: float
    validation_residual_scales: tuple[float, ...]
    minimum_validation_examples_for_activation: int
    product_fingerprint_radius: int
    product_fingerprint_size: int
    feature_groups: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.definition_id.startswith("route_action_policy.v1@"):
            raise ValueError("unexpected route-action policy definition")
        if self.feature_dimension < 128:
            raise ValueError("route-action feature dimension is too small")
        if self.l2_penalty < 0.0:
            raise ValueError("route-action L2 penalty cannot be negative")
        if self.epochs < 1 or self.learning_rate <= 0.0:
            raise ValueError("route-action training schedule must be positive")
        if not 0.0 < self.strategy_equivalent_label_strength <= 1.0:
            raise ValueError("strategy-equivalent label strength must be in (0, 1]")
        if self.planner_probability_deficit_weight < 0.0:
            raise ValueError("planner probability-deficit weight cannot be negative")
        if self.baseline_rank_logit_weight < 0.0:
            raise ValueError("baseline rank-logit weight cannot be negative")
        if (
            not self.validation_residual_scales
            or any(value < 0.0 for value in self.validation_residual_scales)
            or tuple(sorted(set(self.validation_residual_scales)))
            != self.validation_residual_scales
        ):
            raise ValueError("validation residual scales must be unique and ordered")
        if self.minimum_validation_examples_for_activation < 1:
            raise ValueError("minimum validation examples must be positive")
        if self.product_fingerprint_radius < 1:
            raise ValueError("product fingerprint radius must be positive")
        if self.product_fingerprint_size < 64:
            raise ValueError("product fingerprint size is too small")
        if not self.feature_groups:
            raise ValueError("route-action feature groups cannot be empty")


@lru_cache(maxsize=1)
def load_route_action_policy_definition() -> RouteActionPolicyDefinition:
    """Load the versioned route-action learning definition."""

    value = json.loads(_DEFINITION_PATH.read_text(encoding="utf-8"))
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported route-action policy definition schema")
    return RouteActionPolicyDefinition(
        definition_id=str(value["definition_id"]),
        feature_dimension=int(value["feature_dimension"]),
        l2_penalty=float(value["l2_penalty"]),
        epochs=int(value["epochs"]),
        learning_rate=float(value["learning_rate"]),
        strategy_equivalent_label_strength=float(
            value["strategy_equivalent_label_strength"]
        ),
        planner_probability_deficit_weight=float(
            value["planner_probability_deficit_weight"]
        ),
        baseline_rank_logit_weight=float(value["baseline_rank_logit_weight"]),
        validation_residual_scales=tuple(
            float(item) for item in value["validation_residual_scales"]
        ),
        minimum_validation_examples_for_activation=int(
            value["minimum_validation_examples_for_activation"]
        ),
        product_fingerprint_radius=int(value["product_fingerprint_radius"]),
        product_fingerprint_size=int(value["product_fingerprint_size"]),
        feature_groups=tuple(str(item) for item in value["feature_groups"]),
    )


@dataclass(frozen=True)
class RoutePolicyCandidate:
    """Training projection of one validated single-step candidate."""

    candidate_id: str
    candidate_rank: int
    precursor_smiles: str
    abstraction_level: str
    template_id: str
    operator_id: str
    synthon_signature: str
    score: float
    structural_score_band: int
    strategic_complexity_score: float
    strategic_class: str
    independent_reference_support: int
    precursor_compatibility_disposition: str
    supervision_label: str
    source_patent_precedent_overlap: bool

    @classmethod
    def from_replay_candidate(
        cls, candidate: RouteActionCandidate
    ) -> "RoutePolicyCandidate":
        """Project a replay candidate without adding inferred chemistry."""

        return cls(
            candidate_id=digest(
                "RPCA1",
                candidate.template_id,
                candidate.precursor_smiles,
                candidate.strategy_id,
            ),
            candidate_rank=candidate.candidate_rank,
            precursor_smiles=candidate.precursor_smiles,
            abstraction_level=candidate.abstraction_level,
            template_id=candidate.template_id,
            operator_id=candidate.operator_id,
            synthon_signature=candidate.synthon_signature,
            score=candidate.score,
            structural_score_band=candidate.structural_score_band,
            strategic_complexity_score=candidate.strategic_complexity_score,
            strategic_class=candidate.strategic_class,
            independent_reference_support=candidate.independent_reference_support,
            precursor_compatibility_disposition=(
                candidate.precursor_compatibility_disposition
            ),
            supervision_label=candidate.supervision_label,
            source_patent_precedent_overlap=(
                candidate.source_patent_precedent_overlap
            ),
        )


@dataclass(frozen=True)
class RoutePolicyExample:
    """One observed route decision among validated single-step alternatives."""

    example_id: str
    route_id: str
    patent_id: Optional[str]
    split: str
    target_product_smiles: str
    retrosynthetic_depth: int
    observed_remaining_steps: int
    maximum_route_depth: int
    selected_candidate_ids: tuple[str, ...]
    label_source: str
    label_strength: float
    candidates: tuple[RoutePolicyCandidate, ...]
    schema_version: str = ROUTE_POLICY_EXAMPLE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != ROUTE_POLICY_EXAMPLE_SCHEMA_VERSION:
            raise ValueError("unsupported route-policy example schema")
        if len(self.candidates) < 2:
            raise ValueError("route-policy examples require alternatives")
        candidate_ids = {candidate.candidate_id for candidate in self.candidates}
        if len(candidate_ids) != len(self.candidates):
            raise ValueError("route-policy candidate IDs must be unique")
        if not self.selected_candidate_ids or not set(
            self.selected_candidate_ids
        ).issubset(candidate_ids):
            raise ValueError("route-policy selected candidates are inconsistent")
        if self.label_source not in {"observed_exact", "strategy_equivalent"}:
            raise ValueError("unsupported route-policy label source")
        if not 0.0 < self.label_strength <= 1.0:
            raise ValueError("route-policy label strength must be in (0, 1]")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible training example."""

        value = asdict(self)
        value["selected_candidate_ids"] = list(self.selected_candidate_ids)
        value["candidates"] = [asdict(candidate) for candidate in self.candidates]
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RoutePolicyExample":
        """Reconstruct a route-policy example."""

        fields = dict(value)
        fields["selected_candidate_ids"] = tuple(
            str(item) for item in fields.get("selected_candidate_ids") or ()
        )
        fields["candidates"] = tuple(
            RoutePolicyCandidate(**item) for item in fields.get("candidates") or ()
        )
        return cls(**fields)


def build_route_policy_examples(
    evaluations: Iterable[RouteActionEvaluation],
    *,
    definition: Optional[RouteActionPolicyDefinition] = None,
) -> tuple[RoutePolicyExample, ...]:
    """Build listwise route choices only where replay recovered the observation."""

    policy = definition or load_route_action_policy_definition()
    examples = []
    for route in evaluations:
        route_id = route.source_route_id or route.tree_id
        for step in route.steps:
            if step.search_status != "searched" or len(step.candidates) < 2:
                continue
            candidates = tuple(
                RoutePolicyCandidate.from_replay_candidate(candidate)
                for candidate in step.candidates
            )
            exact = tuple(
                candidate.candidate_id
                for candidate in candidates
                if candidate.supervision_label == "observed_exact"
            )
            strategies = tuple(
                candidate.candidate_id
                for candidate in candidates
                if candidate.supervision_label == "strategy_equivalent"
            )
            if exact:
                selected = exact
                source = "observed_exact"
                strength = 1.0
            elif strategies:
                selected = strategies
                source = "strategy_equivalent"
                strength = policy.strategy_equivalent_label_strength
            else:
                continue
            examples.append(
                RoutePolicyExample(
                    example_id=digest(
                        "RPE1",
                        step.evaluation_id,
                        source,
                        *selected,
                    ),
                    route_id=route_id,
                    patent_id=route.patent_id,
                    split=route.split or "unknown",
                    target_product_smiles=(
                        step.observed_action.target_smiles or ""
                    ),
                    retrosynthetic_depth=step.retrosynthetic_depth,
                    observed_remaining_steps=step.observed_remaining_steps,
                    maximum_route_depth=route.maximum_depth,
                    selected_candidate_ids=selected,
                    label_source=source,
                    label_strength=strength,
                    candidates=candidates,
                )
            )
    return tuple(sorted(examples, key=lambda item: item.example_id))


def _hashed_index(token: str, dimension: int) -> tuple[int, float]:
    token_digest = hashlib.blake2b(token.encode("utf-8"), digest_size=8).digest()
    value = int.from_bytes(token_digest, byteorder="big", signed=False)
    return value % dimension, 1.0 if (value >> 63) == 0 else -1.0


@lru_cache(maxsize=100_000)
def _product_fingerprint_bits(
    smiles: str,
    radius: int,
    size: int,
) -> tuple[int, ...]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return ()
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=size,
    )
    return tuple(int(value) for value in generator.GetFingerprint(molecule).GetOnBits())


@dataclass(frozen=True)
class RoutePolicyCandidateAssessment:
    """Learned preference for one already validated action."""

    candidate: GenericDisconnectionCandidate
    probability: float
    policy_rank: int
    original_rank: int


@dataclass(frozen=True)
class RouteActionPolicyTrainingReport:
    """Deterministic training and held-out ranking metrics."""

    training_example_count: int
    training_candidate_count: int
    epoch_count: int
    initial_loss: float
    final_loss: float
    selected_residual_scale: float
    policy_active: bool
    activation_reason: str
    baseline_metrics_by_split: dict[str, dict[str, Optional[float] | int]]
    metrics_by_split: dict[str, dict[str, Optional[float] | int]]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible report."""

        return asdict(self)


class RouteActionPolicyModel:
    """Deterministic hashed listwise policy over validated disconnections."""

    def __init__(
        self,
        *,
        definition: Optional[RouteActionPolicyDefinition] = None,
        weights: Optional[Sequence[float]] = None,
        residual_scale: float = 1.0,
    ) -> None:
        self.definition = definition or load_route_action_policy_definition()
        if weights is None:
            self.weights = [0.0] * self.definition.feature_dimension
        else:
            if len(weights) != self.definition.feature_dimension:
                raise ValueError("route-action policy weight dimension mismatch")
            self.weights = [float(value) for value in weights]
        if residual_scale < 0.0:
            raise ValueError("route-action residual scale cannot be negative")
        self.residual_scale = float(residual_scale)

    @property
    def model_id(self) -> str:
        """Return a content-derived identity for the learned policy."""

        normalized_weights = [
            0.0 if abs(weight) < 5e-13 else round(weight, 12)
            for weight in self.weights
        ]
        payload = json.dumps(
            normalized_weights,
            separators=(",", ":"),
        )
        return digest(
            "RAPM1",
            self.definition.definition_id,
            f"residual_scale={self.residual_scale:.8f}",
            payload,
        )

    @property
    def planner_probability_deficit_weight(self) -> float:
        """Return the versioned planner integration weight."""

        if self.residual_scale <= 0.0:
            return 0.0
        return self.definition.planner_probability_deficit_weight

    def _features(
        self,
        product_smiles: str,
        retrosynthetic_depth: int,
        candidate: RoutePolicyCandidate | GenericDisconnectionCandidate,
        candidate_rank: int,
    ) -> dict[int, float]:
        score_bin = max(0, min(20, int(round(float(candidate.score) * 20))))
        complexity_bin = max(
            -20,
            min(20, int(round(float(candidate.strategic_complexity_score) * 10))),
        )
        support_bin = min(
            12,
            int(math.log2(max(1, int(candidate.independent_reference_support)))),
        )
        groups = set(self.definition.feature_groups)
        tokens = {"bias"}
        if "existing_candidate_rank" in groups:
            tokens.add(f"candidate_rank={min(candidate_rank, 20)}")
        if "structural_score_band" in groups:
            tokens.add(f"structural_band={int(candidate.structural_score_band)}")
        if "candidate_score_bin" in groups:
            tokens.add(f"score_bin={score_bin}")
        if "strategic_complexity_bin" in groups:
            tokens.add(f"complexity_bin={complexity_bin}")
        if "strategic_class" in groups:
            tokens.add(f"strategic_class={candidate.strategic_class}")
        if "abstraction_level" in groups:
            tokens.add(f"level={candidate.abstraction_level}")
        if "precursor_compatibility" in groups:
            tokens.add(
                "compatibility="
                f"{candidate.precursor_compatibility_disposition}"
            )
        if "independent_reference_support_bin" in groups:
            tokens.add(f"support_bin={support_bin}")
        if "retrosynthetic_depth" in groups:
            tokens.add(f"depth={min(max(0, retrosynthetic_depth), 8)}")
        if "operator_identity" in groups:
            tokens.add(f"operator={candidate.operator_id}")
        if "synthon_identity" in groups:
            tokens.add(f"synthon={candidate.synthon_signature}")
        if {"product_morgan_bits", "product_bit_by_operator"}.intersection(groups):
            product_bits = _product_fingerprint_bits(
                product_smiles,
                self.definition.product_fingerprint_radius,
                self.definition.product_fingerprint_size,
            )
            for bit in product_bits:
                if "product_morgan_bits" in groups:
                    tokens.add(f"product_bit={bit}")
                if "product_bit_by_operator" in groups:
                    tokens.add(
                        f"product_bit={bit}|operator={candidate.operator_id}"
                    )
        features: dict[int, float] = {}
        for token in sorted(tokens):
            index, sign = _hashed_index(token, self.definition.feature_dimension)
            features[index] = features.get(index, 0.0) + sign
        return features

    def _score(self, features: Mapping[int, float]) -> float:
        return sum(self.weights[index] * value for index, value in features.items())

    def _candidate_scores(
        self,
        features: Sequence[Mapping[int, float]],
        candidate_ranks: Sequence[int],
    ) -> tuple[float, ...]:
        return tuple(
            -self.definition.baseline_rank_logit_weight * max(0, rank - 1)
            + self.residual_scale * self._score(candidate_features)
            for candidate_features, rank in zip(
                features,
                candidate_ranks,
                strict=True,
            )
        )

    @staticmethod
    def _softmax(scores: Sequence[float]) -> tuple[float, ...]:
        if not scores:
            return ()
        maximum = max(scores)
        exponentials = tuple(math.exp(score - maximum) for score in scores)
        denominator = sum(exponentials)
        return tuple(value / denominator for value in exponentials)

    def predict_example_probabilities(
        self, example: RoutePolicyExample
    ) -> dict[str, float]:
        """Return normalized preferences for one replay choice set."""

        features = tuple(
            self._features(
                example.target_product_smiles,
                example.retrosynthetic_depth,
                candidate,
                candidate.candidate_rank,
            )
            for candidate in example.candidates
        )
        probabilities = self._softmax(
            self._candidate_scores(
                features,
                tuple(candidate.candidate_rank for candidate in example.candidates),
            )
        )
        return {
            candidate.candidate_id: probability
            for candidate, probability in zip(
                example.candidates,
                probabilities,
                strict=True,
            )
        }

    def _loss(self, examples: Sequence[RoutePolicyExample]) -> float:
        total_weight = sum(example.label_strength for example in examples)
        if total_weight <= 0.0:
            return 0.0
        loss = 0.0
        for example in examples:
            probabilities = self.predict_example_probabilities(example)
            target_count = len(example.selected_candidate_ids)
            selected_loss = sum(
                -math.log(max(probabilities[candidate_id], 1e-15))
                for candidate_id in example.selected_candidate_ids
            ) / target_count
            loss += example.label_strength * selected_loss
        loss /= total_weight
        loss += 0.5 * self.definition.l2_penalty * sum(
            weight * weight for weight in self.weights
        )
        return loss

    def fit(
        self,
        examples: Iterable[RoutePolicyExample],
        *,
        evaluation_examples: Iterable[RoutePolicyExample] = (),
    ) -> RouteActionPolicyTrainingReport:
        """Fit a deterministic full-batch weighted listwise softmax model."""

        training = tuple(sorted(examples, key=lambda item: item.example_id))
        if not training:
            raise ValueError("route-action policy training requires examples")
        self.weights = [0.0] * self.definition.feature_dimension
        self.residual_scale = 1.0
        feature_cache = {
            (example.example_id, candidate.candidate_id): self._features(
                example.target_product_smiles,
                example.retrosynthetic_depth,
                candidate,
                candidate.candidate_rank,
            )
            for example in training
            for candidate in example.candidates
        }
        initial_loss = self._loss(training)
        total_weight = sum(example.label_strength for example in training)
        for _ in range(self.definition.epochs):
            gradient = [0.0] * self.definition.feature_dimension
            for example in training:
                candidate_features = tuple(
                    feature_cache[(example.example_id, candidate.candidate_id)]
                    for candidate in example.candidates
                )
                probabilities = self._softmax(
                    self._candidate_scores(
                        candidate_features,
                        tuple(
                            candidate.candidate_rank
                            for candidate in example.candidates
                        ),
                    )
                )
                target_probability = 1.0 / len(example.selected_candidate_ids)
                selected = set(example.selected_candidate_ids)
                for candidate, probability, features in zip(
                    example.candidates,
                    probabilities,
                    candidate_features,
                    strict=True,
                ):
                    target = target_probability if candidate.candidate_id in selected else 0.0
                    factor = (
                        self.residual_scale
                        * example.label_strength
                        * (probability - target)
                    )
                    for index, feature_value in features.items():
                        gradient[index] += factor * feature_value
            for index, weight in enumerate(self.weights):
                normalized = gradient[index] / total_weight
                normalized += self.definition.l2_penalty * weight
                self.weights[index] -= self.definition.learning_rate * normalized
        trained_loss = self._loss(training)
        evaluation = tuple(evaluation_examples)
        validation = tuple(item for item in evaluation if item.split == "validation")
        minimum_validation = (
            self.definition.minimum_validation_examples_for_activation
        )
        if len(validation) >= minimum_validation:
            best_scale = 0.0
            best_key = (-1.0, -1.0, 0.0)
            for scale in self.definition.validation_residual_scales:
                self.residual_scale = scale
                metrics = evaluate_route_action_policy(self, validation)["validation"]
                key = (
                    float(metrics["mean_reciprocal_rank"] or 0.0),
                    float(metrics["top1"] or 0.0),
                    -scale,
                )
                if key > best_key:
                    best_key = key
                    best_scale = scale
            self.residual_scale = best_scale
            activation_reason = (
                "heldout_validation_selected_residual"
                if best_scale > 0.0
                else "heldout_validation_selected_baseline"
            )
        else:
            self.residual_scale = 0.0
            activation_reason = (
                "insufficient_validation_examples:"
                f"{len(validation)}<{minimum_validation}"
            )
        evaluated = (*training, *evaluation)
        return RouteActionPolicyTrainingReport(
            training_example_count=len(training),
            training_candidate_count=sum(len(item.candidates) for item in training),
            epoch_count=self.definition.epochs,
            initial_loss=round(initial_loss, 8),
            final_loss=round(trained_loss, 8),
            selected_residual_scale=self.residual_scale,
            policy_active=self.residual_scale > 0.0,
            activation_reason=activation_reason,
            baseline_metrics_by_split=evaluate_route_action_baseline(evaluated),
            metrics_by_split=evaluate_route_action_policy(self, evaluated),
        )

    def rank_candidates(
        self,
        product_smiles: str,
        candidates: Sequence[GenericDisconnectionCandidate],
        *,
        retrosynthetic_depth: int,
    ) -> tuple[RoutePolicyCandidateAssessment, ...]:
        """Rerank already validated candidates without changing admission."""

        values = tuple(candidates)
        if not values:
            return ()
        features = tuple(
            self._features(product_smiles, retrosynthetic_depth, candidate, rank)
            for rank, candidate in enumerate(values, 1)
        )
        probabilities = self._softmax(
            self._candidate_scores(features, tuple(range(1, len(values) + 1)))
        )
        ordered = sorted(
            zip(values, probabilities, range(1, len(values) + 1), strict=True),
            key=lambda item: (-item[1], item[2], item[0].template_id),
        )
        return tuple(
            RoutePolicyCandidateAssessment(
                candidate=candidate,
                probability=round(probability, 8),
                policy_rank=rank,
                original_rank=original_rank,
            )
            for rank, (candidate, probability, original_rank) in enumerate(ordered, 1)
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a compact JSON-compatible learned model."""

        sparse_weights = []
        for index, weight in enumerate(self.weights):
            normalized = 0.0 if abs(weight) < 5e-13 else round(weight, 12)
            if normalized:
                sparse_weights.append([index, normalized])
        return {
            "schema_version": ROUTE_ACTION_POLICY_SCHEMA_VERSION,
            "definition": asdict(self.definition),
            "model_id": self.model_id,
            "residual_scale": self.residual_scale,
            "sparse_weights": sparse_weights,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RouteActionPolicyModel":
        """Reconstruct and validate a serialized route-action model."""

        if value.get("schema_version") != ROUTE_ACTION_POLICY_SCHEMA_VERSION:
            raise ValueError("unsupported route-action policy model schema")
        raw_definition = dict(value.get("definition") or {})
        raw_definition["feature_groups"] = tuple(
            str(item) for item in raw_definition.get("feature_groups") or ()
        )
        raw_definition["validation_residual_scales"] = tuple(
            float(item)
            for item in raw_definition.get("validation_residual_scales") or ()
        )
        definition = RouteActionPolicyDefinition(**raw_definition)
        weights = [0.0] * definition.feature_dimension
        for raw_index, raw_weight in value.get("sparse_weights") or ():
            index = int(raw_index)
            if index < 0 or index >= len(weights):
                raise ValueError("route-action sparse weight index is invalid")
            weights[index] = float(raw_weight)
        model = cls(
            definition=definition,
            weights=weights,
            residual_scale=float(value.get("residual_scale") or 0.0),
        )
        if value.get("model_id") != model.model_id:
            raise ValueError("route-action policy model ID mismatch")
        return model


def evaluate_route_action_policy(
    model: RouteActionPolicyModel,
    examples: Iterable[RoutePolicyExample],
) -> dict[str, dict[str, Optional[float] | int]]:
    """Report top-1 and reciprocal-rank recovery by source split."""

    grouped: dict[str, list[RoutePolicyExample]] = {}
    for example in examples:
        grouped.setdefault(example.split or "unknown", []).append(example)
    report: dict[str, dict[str, Optional[float] | int]] = {}
    for split, values in sorted(grouped.items()):
        top1 = 0
        reciprocal_rank = 0.0
        for example in values:
            probabilities = model.predict_example_probabilities(example)
            ranked_ids = sorted(
                probabilities,
                key=lambda candidate_id: (-probabilities[candidate_id], candidate_id),
            )
            selected = set(example.selected_candidate_ids)
            rank = next(
                index
                for index, candidate_id in enumerate(ranked_ids, 1)
                if candidate_id in selected
            )
            top1 += int(rank == 1)
            reciprocal_rank += 1.0 / rank
        count = len(values)
        report[split] = {
            "example_count": count,
            "top1": top1 / count if count else None,
            "mean_reciprocal_rank": reciprocal_rank / count if count else None,
        }
    return report


def evaluate_route_action_baseline(
    examples: Iterable[RoutePolicyExample],
) -> dict[str, dict[str, Optional[float] | int]]:
    """Report recovery from the original validated single-step order."""

    grouped: dict[str, list[RoutePolicyExample]] = {}
    for example in examples:
        grouped.setdefault(example.split or "unknown", []).append(example)
    report: dict[str, dict[str, Optional[float] | int]] = {}
    for split, values in sorted(grouped.items()):
        top1 = 0
        reciprocal_rank = 0.0
        for example in values:
            ranked = sorted(
                example.candidates,
                key=lambda candidate: (
                    candidate.candidate_rank,
                    candidate.candidate_id,
                ),
            )
            selected = set(example.selected_candidate_ids)
            rank = next(
                index
                for index, candidate in enumerate(ranked, 1)
                if candidate.candidate_id in selected
            )
            top1 += int(rank == 1)
            reciprocal_rank += 1.0 / rank
        count = len(values)
        report[split] = {
            "example_count": count,
            "top1": top1 / count if count else None,
            "mean_reciprocal_rank": reciprocal_rank / count if count else None,
        }
    return report


def _write_deterministic_gzip(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw_handle, mtime=0
            ) as gzip_handle:
                with io.TextIOWrapper(gzip_handle, encoding="utf-8") as handle:
                    json.dump(value, handle, sort_keys=True, separators=(",", ":"))
        os.replace(temporary, path)
    except BaseException:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise


def save_route_action_policy(
    model: RouteActionPolicyModel,
    destination: str | Path,
) -> None:
    """Write a deterministic compressed route-action policy."""

    _write_deterministic_gzip(Path(destination), model.to_dict())


def load_route_action_policy(source: str | Path) -> RouteActionPolicyModel:
    """Load and validate a compressed route-action policy."""

    with gzip.open(Path(source), "rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("route-action policy artifact must contain an object")
    return RouteActionPolicyModel.from_dict(value)


def train_route_action_policy_from_replay(
    source_replay: str | Path,
    output_model: str | Path,
    *,
    report_path: Optional[str | Path] = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Train on train-split replay choices and audit held-out route splits."""

    source = Path(source_replay)
    output = Path(output_model)
    report_output = (
        Path(report_path) if report_path is not None else Path(f"{output}.report.json")
    )
    if not source.is_file():
        raise FileNotFoundError(source)
    if not overwrite:
        for path in (output, report_output):
            if path.exists():
                raise FileExistsError(path)
    evaluations = tuple(iter_route_action_evaluations(source))
    examples = build_route_policy_examples(evaluations)
    training = tuple(example for example in examples if example.split == "train")
    held_out = tuple(example for example in examples if example.split != "train")
    model = RouteActionPolicyModel()
    training_report = model.fit(training, evaluation_examples=held_out)
    save_route_action_policy(model, output)
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    output_hash = hashlib.sha256(output.read_bytes()).hexdigest()
    counts = Counter(example.label_source for example in examples)
    report = {
        "route_action_policy_schema_version": ROUTE_ACTION_POLICY_SCHEMA_VERSION,
        "definition_id": model.definition.definition_id,
        "model_id": model.model_id,
        "source": {"path": str(source.resolve()), "sha256": source_hash},
        "output": {"path": str(output.resolve()), "sha256": output_hash},
        "example_count": len(examples),
        "example_counts_by_label_source": dict(sorted(counts.items())),
        "training": training_report.to_dict(),
        "warnings": [
            "The model ranks only candidates admitted by deterministic chemistry.",
            "Unchosen alternatives provide relative preference, not hard negatives.",
            "Patent and close-analogue leakage must be controlled before release evaluation.",
        ]
        + (
            []
            if training_report.policy_active
            else [
                "The learned residual is inactive and has zero planner influence: "
                f"{training_report.activation_reason}."
            ]
        ),
    }
    report_output.parent.mkdir(parents=True, exist_ok=True)
    report_output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return report


__all__ = [
    "ROUTE_ACTION_POLICY_SCHEMA_VERSION",
    "ROUTE_POLICY_EXAMPLE_SCHEMA_VERSION",
    "RouteActionPolicyDefinition",
    "RouteActionPolicyModel",
    "RouteActionPolicyTrainingReport",
    "RoutePolicyCandidate",
    "RoutePolicyCandidateAssessment",
    "RoutePolicyExample",
    "build_route_policy_examples",
    "evaluate_route_action_baseline",
    "evaluate_route_action_policy",
    "load_route_action_policy",
    "load_route_action_policy_definition",
    "save_route_action_policy",
    "train_route_action_policy_from_replay",
]
