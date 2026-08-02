"""Chemist-oriented ML grouping for weak-label condition recipes.

The proof of concept separates the decisive condition core from protocol
context.  Catalyst/ligand, activation reagent, acid/base, and redox components
define the learned core.  Solvents, additives, and operating conditions remain
attached as searchable cross-references and never disappear from the evidence.

This module is isolated from production retrieval and condition identity.
Generated group identifiers are snapshot-bound review proposals.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import statistics
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple


CONDITION_GROUPING_POC_SCHEMA_VERSION = "2.0"
CONDITION_GROUPING_POC_DEFINITION_VERSION = "condition_grouping_poc.v2.1"
CONDITION_MATERIAL_ID_NAMESPACE = "CMAT1"
CONDITION_CORE_ID_NAMESPACE = "CCOREPOC1"
CONDITION_GROUP_ID_NAMESPACE = "CGPOC2"

# These roles describe what chemists ordinarily name as the decisive system.
# The policy is explicit and reviewable; source column names alone do not assign
# roles, because the weak-label input already supplies contextual role labels.
CORE_COMPONENT_ROLES = frozenset(
    {
        "acid",
        "base",
        "catalyst",
        "condensation_agent",
        "coupling_reagent",
        "ligand",
        "oxidant",
        "reductant",
    }
)
MEDIUM_COMPONENT_ROLES = frozenset({"solvent"})
MODIFIER_COMPONENT_ROLES = frozenset({"additive"})

SUPPORTED_SIMILARITY_THRESHOLD = 0.9
SUPPORTED_MARGIN_THRESHOLD = 0.05
PROVISIONAL_SIMILARITY_THRESHOLD = 0.8
PROVISIONAL_MARGIN_THRESHOLD = 0.02

_COMPONENT_PATTERN = re.compile(r"^\s*(.*?)\s*\[([^\[\]]+)\]\s*$")
_REQUIRED_COLUMNS = (
    "source_reaction_type",
    "condition_recipe_id",
    "condition_display",
    "yield_pct",
    "temperature_c",
    "time_h",
    "condition_identity_uncertainty",
)


@dataclass(frozen=True)
class ParsedConditionComponent:
    """One normalized component and its source-provided contextual role."""

    name: str
    normalized_name: str
    role: str


@dataclass(frozen=True)
class ParsedConditionDisplay:
    """Order-invariant material system parsed from a weak-label display."""

    components: Tuple[ParsedConditionComponent, ...]
    declared_absences: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()

    @property
    def identity_payload(self) -> Dict[str, Any]:
        """Return the complete material payload used for exact identity."""

        return {
            "components": tuple(
                (component.role, component.normalized_name)
                for component in self.components
            ),
            "declared_absences": self.declared_absences,
            "schema_version": CONDITION_GROUPING_POC_SCHEMA_VERSION,
        }

    @property
    def material_id(self) -> str:
        """Return a stable identity for the complete normalized material set."""

        return _digest(CONDITION_MATERIAL_ID_NAMESPACE, self.identity_payload)

    @property
    def canonical_display(self) -> str:
        """Render the complete role-aware material system for review."""

        return _render_components(self.components, self.declared_absences)

    @property
    def core_components(self) -> Tuple[ParsedConditionComponent, ...]:
        """Return components that define the chemist-facing condition core."""

        return tuple(
            component
            for component in self.components
            if component.role in CORE_COMPONENT_ROLES
        )

    @property
    def context_components(self) -> Tuple[ParsedConditionComponent, ...]:
        """Return medium, modifier, and unresolved contextual components."""

        return tuple(
            component
            for component in self.components
            if component.role not in CORE_COMPONENT_ROLES
        )

    @property
    def core_declared_absences(self) -> Tuple[str, ...]:
        """Retain only absences that alter an observed decisive subsystem."""

        roles = {component.role for component in self.core_components}
        values = set()
        if "catalyst" in roles and "ligand" not in roles:
            if "ligand" in self.declared_absences:
                values.add("ligand")
        return tuple(sorted(values))

    @property
    def core_id(self) -> Optional[str]:
        """Return the exact decisive-core identity, or ``None`` if unresolved."""

        if not self.core_components:
            return None
        return _digest(
            CONDITION_CORE_ID_NAMESPACE,
            {
                "components": tuple(
                    (component.role, component.normalized_name)
                    for component in self.core_components
                ),
                "declared_absences": self.core_declared_absences,
                "core_role_policy": sorted(CORE_COMPONENT_ROLES),
                "schema_version": CONDITION_GROUPING_POC_SCHEMA_VERSION,
            },
        )

    @property
    def core_display(self) -> str:
        """Render only the decisive system used for grouping."""

        return _render_components(
            self.core_components,
            self.core_declared_absences,
        )

    @property
    def core_feature_tokens(self) -> Tuple[str, ...]:
        """Return interpretable decisive-core tokens for latent grouping."""

        tokens = []
        by_role: Dict[str, list[str]] = defaultdict(list)
        for component in self.core_components:
            by_role[component.role].append(component.normalized_name)
            tokens.append(f"core_role:{component.role}")
            tokens.append(
                f"core_component:{component.role}={component.normalized_name}"
            )
        for role in self.core_declared_absences:
            tokens.append(f"core_absence:{role}")
        topology = "|".join(
            f"{role}:{len(names)}" for role, names in sorted(by_role.items())
        )
        if self.core_declared_absences:
            topology += "|" + "|".join(
                f"absence:{role}" for role in self.core_declared_absences
            )
        tokens.append(f"core_topology:{topology}")
        for catalyst in sorted(by_role.get("catalyst", ())):
            for ligand in sorted(by_role.get("ligand", ())):
                tokens.append(f"core_pair:catalyst_ligand={catalyst}|{ligand}")
            for base in sorted(by_role.get("base", ())):
                tokens.append(f"core_pair:catalyst_base={catalyst}|{base}")
        activation_roles = ("condensation_agent", "coupling_reagent")
        for activation_role in activation_roles:
            for reagent in sorted(by_role.get(activation_role, ())):
                for base in sorted(by_role.get("base", ())):
                    tokens.append(
                        f"core_pair:{activation_role}_base={reagent}|{base}"
                    )
        return tuple(sorted(tokens))

    def role_names(self, role: str) -> Tuple[str, ...]:
        """Return normalized component names assigned to one contextual role."""

        return tuple(
            component.normalized_name
            for component in self.components
            if component.role == role
        )


@dataclass
class MaterialAggregate:
    """Dataset evidence for one complete normalized material system."""

    parsed: ParsedConditionDisplay
    observation_count: int = 0
    recipe_ids: set[str] = field(default_factory=set)
    reaction_types: Counter[str] = field(default_factory=Counter)
    yield_sum: float = 0.0
    yield_count: int = 0
    uncertain_count: int = 0
    temperatures: Counter[float] = field(default_factory=Counter)
    times: Counter[float] = field(default_factory=Counter)
    source_displays: Counter[str] = field(default_factory=Counter)

    @property
    def material_id(self) -> str:
        return self.parsed.material_id


@dataclass
class ConditionCoreAggregate:
    """Evidence and material variants belonging to one exact decisive core."""

    core_id: str
    core_display: str
    feature_tokens: Tuple[str, ...]
    materials: list[MaterialAggregate] = field(default_factory=list)

    @property
    def observation_count(self) -> int:
        return sum(material.observation_count for material in self.materials)

    @property
    def recipe_ids(self) -> set[str]:
        return set().union(*(material.recipe_ids for material in self.materials))


@dataclass(frozen=True)
class ConditionGroupingRun:
    """Paths and report returned by one completed POC run."""

    report: Mapping[str, Any]
    report_json_path: Path
    report_markdown_path: Path
    groups_path: Path
    assignments_path: Path
    model_path: Path


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    )


def _identity_analyzer(value: Iterable[str]) -> Iterable[str]:
    """Return pre-tokenized features for scikit-learn serialization."""

    return value


def _digest(namespace: str, payload: Any) -> str:
    value = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
    return f"{namespace}:{value}"


def _normalized_text(value: str) -> str:
    return " ".join(unicodedata.normalize("NFKC", value).strip().split()).casefold()


def _display_text(value: str) -> str:
    return " ".join(unicodedata.normalize("NFKC", value).strip().split())


def _render_components(
    components: Iterable[ParsedConditionComponent],
    declared_absences: Iterable[str] = (),
) -> str:
    values = [f"{component.name} [{component.role}]" for component in components]
    values.extend(f"No {role}" for role in declared_absences)
    return "; ".join(values)


def parse_condition_display(value: str) -> ParsedConditionDisplay:
    """Parse ``name [role]`` segments without inventing missing chemistry."""

    components = []
    declared_absences = set()
    warnings = set()
    for raw_segment in str(value or "").split(";"):
        segment = _display_text(raw_segment)
        if not segment:
            continue
        absence = re.fullmatch(r"no\s+([A-Za-z_ -]+)", segment, flags=re.I)
        if absence:
            declared_absences.add(
                _normalized_text(absence.group(1)).replace(" ", "_")
            )
            continue
        match = _COMPONENT_PATTERN.fullmatch(segment)
        if match is None:
            warnings.add(f"UNPARSED_SEGMENT:{segment}")
            continue
        name = _display_text(match.group(1))
        role = _normalized_text(match.group(2)).replace(" ", "_")
        if not name or not role:
            warnings.add(f"EMPTY_COMPONENT_OR_ROLE:{segment}")
            continue
        components.append(
            ParsedConditionComponent(
                name=name,
                normalized_name=_normalized_text(name),
                role=role,
            )
        )
    ordered = tuple(
        sorted(
            components,
            key=lambda item: (item.role, item.normalized_name, item.name),
        )
    )
    if not ordered and not declared_absences:
        warnings.add("NO_PARSED_CONDITION_COMPONENTS")
    if "ligand" in declared_absences and any(
        component.role == "ligand" for component in ordered
    ):
        warnings.add("CONTRADICTORY_LIGAND_PRESENCE_AND_ABSENCE")
    return ParsedConditionDisplay(
        components=ordered,
        declared_absences=tuple(sorted(declared_absences)),
        warnings=tuple(sorted(warnings)),
    )


def _optional_float(value: Any) -> Optional[float]:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def _truthy(value: Any) -> bool:
    return str(value or "").strip().casefold() in {"1", "true", "yes", "y"}


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_weak_label_materials(
    dataset_path: str | Path,
) -> Tuple[Tuple[MaterialAggregate, ...], Dict[str, Any]]:
    """Aggregate the weak-label CSV by complete normalized material system."""

    source = Path(dataset_path)
    aggregates: Dict[str, MaterialAggregate] = {}
    row_count = 0
    parse_warning_counts: Counter[str] = Counter()
    unique_recipe_ids = set()
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = sorted(set(_REQUIRED_COLUMNS).difference(reader.fieldnames or ()))
        if missing:
            raise ValueError(
                "weak-label condition grouping dataset is missing columns: "
                + ", ".join(missing)
            )
        for row in reader:
            row_count += 1
            parsed = parse_condition_display(row["condition_display"])
            parse_warning_counts.update(parsed.warnings)
            aggregate = aggregates.get(parsed.material_id)
            if aggregate is None:
                aggregate = MaterialAggregate(parsed=parsed)
                aggregates[parsed.material_id] = aggregate
            aggregate.observation_count += 1
            recipe_id = str(row["condition_recipe_id"] or "").strip()
            if recipe_id:
                aggregate.recipe_ids.add(recipe_id)
                unique_recipe_ids.add(recipe_id)
            reaction_type = str(row["source_reaction_type"] or "").strip()
            if reaction_type:
                aggregate.reaction_types[reaction_type] += 1
            outcome = _optional_float(row["yield_pct"])
            if outcome is not None:
                aggregate.yield_sum += outcome
                aggregate.yield_count += 1
            aggregate.uncertain_count += int(
                _truthy(row["condition_identity_uncertainty"])
            )
            temperature = _optional_float(row["temperature_c"])
            if temperature is not None:
                aggregate.temperatures[temperature] += 1
            time_h = _optional_float(row["time_h"])
            if time_h is not None:
                aggregate.times[time_h] += 1
            source_display = str(row["condition_display"] or "").strip()
            if source_display:
                aggregate.source_displays[source_display] += 1
    ordered = tuple(sorted(aggregates.values(), key=lambda item: item.material_id))
    audit = {
        "row_count": row_count,
        "unique_material_count": len(ordered),
        "unique_recipe_id_count": len(unique_recipe_ids),
        "parse_warning_counts": dict(sorted(parse_warning_counts.items())),
    }
    return ordered, audit


def build_condition_cores(
    materials: Iterable[MaterialAggregate],
) -> Tuple[Tuple[ConditionCoreAggregate, ...], Tuple[MaterialAggregate, ...]]:
    """Factor exact materials into decisive cores and unresolved core records."""

    cores: Dict[str, ConditionCoreAggregate] = {}
    unresolved = []
    for material in materials:
        core_id = material.parsed.core_id
        if core_id is None:
            unresolved.append(material)
            continue
        core = cores.get(core_id)
        if core is None:
            core = ConditionCoreAggregate(
                core_id=core_id,
                core_display=material.parsed.core_display,
                feature_tokens=material.parsed.core_feature_tokens,
            )
            cores[core_id] = core
        core.materials.append(material)
    return (
        tuple(sorted(cores.values(), key=lambda item: item.core_id)),
        tuple(sorted(unresolved, key=lambda item: item.material_id)),
    )


def _percentile(values: Sequence[float], fraction: float) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(float(value) for value in values)
    if len(ordered) == 1:
        return ordered[0]
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    remainder = position - lower
    return ordered[lower] * (1.0 - remainder) + ordered[upper] * remainder


def _top_items(counter: Counter[str], limit: int = 10) -> Tuple[Dict[str, Any], ...]:
    return tuple(
        {"value": value, "count": count}
        for value, count in sorted(counter.items(), key=lambda item: (-item[1], item[0]))[
            :limit
        ]
    )


def _role_component_counts(
    materials: Iterable[MaterialAggregate],
) -> Dict[str, Counter[str]]:
    values: Dict[str, Counter[str]] = defaultdict(Counter)
    for material in materials:
        for component in material.parsed.components:
            values[component.role][component.name] += material.observation_count
    return values


def _role_profile_counts(
    materials: Iterable[MaterialAggregate], role: str
) -> Counter[str]:
    values: Counter[str] = Counter()
    for material in materials:
        names = material.parsed.role_names(role)
        if names:
            values[" + ".join(names)] += material.observation_count
    return values


def _counter_sum(
    materials: Iterable[MaterialAggregate], attribute: str
) -> Counter[Any]:
    result: Counter[Any] = Counter()
    for material in materials:
        result.update(getattr(material, attribute))
    return result


def _numeric_cross_reference(
    values: Counter[float], observation_count: int
) -> Dict[str, Any]:
    count = sum(values.values())
    if not count:
        return {
            "observed_count": 0,
            "coverage_rate": 0.0,
            "minimum": None,
            "median": None,
            "maximum": None,
            "mean": None,
        }
    ordered = sorted(values.items())
    midpoint = (count - 1) / 2
    cumulative = 0
    median = ordered[-1][0]
    for value, frequency in ordered:
        cumulative += frequency
        if cumulative > midpoint:
            median = value
            break
    return {
        "observed_count": count,
        "coverage_rate": round(count / observation_count, 6),
        "minimum": ordered[0][0],
        "median": median,
        "maximum": ordered[-1][0],
        "mean": round(
            sum(value * frequency for value, frequency in ordered) / count,
            6,
        ),
    }


def _topology_key(core: ConditionCoreAggregate) -> str:
    roles = []
    if core.materials:
        parsed = core.materials[0].parsed
        counts = Counter(component.role for component in parsed.core_components)
        roles.extend(f"{role}:{count}" for role, count in sorted(counts.items()))
        roles.extend(
            f"absence:{role}" for role in parsed.core_declared_absences
        )
    return "|".join(roles)


def _safe_round(value: Optional[float], digits: int = 6) -> Optional[float]:
    return round(float(value), digits) if value is not None else None


def _assignment_status(similarity: float, margin: float) -> str:
    """Apply an explicit, uncalibrated POC review policy."""

    if (
        similarity >= SUPPORTED_SIMILARITY_THRESHOLD
        and margin >= SUPPORTED_MARGIN_THRESHOLD
    ):
        return "supported"
    if (
        similarity >= PROVISIONAL_SIMILARITY_THRESHOLD
        and margin >= PROVISIONAL_MARGIN_THRESHOLD
    ):
        return "provisional"
    return "review"


def _render_report(report: Mapping[str, Any]) -> str:
    dataset = report["dataset"]
    model = report["model"]
    quality = report["quality"]
    lines = [
        "# Chemist-oriented weak-label condition grouping POC",
        "",
        "The learned identity covers decisive condition cores. Solvents, additives,",
        "temperature, and time remain attached as cross-referenced protocol context.",
        "Reaction labels and outcomes were excluded from clustering.",
        "",
        "## Dataset",
        "",
        f"- Rows: {dataset['row_count']:,}",
        f"- Exact recipe IDs: {dataset['unique_recipe_id_count']:,}",
        f"- Complete material systems: {dataset['unique_material_count']:,}",
        f"- Exact decisive cores: {dataset['unique_core_count']:,}",
        f"- Material systems without a decisive core: {dataset['unresolved_core_material_count']:,}",
        f"- Input SHA-256: `{dataset['sha256']}`",
        "",
        "## Model",
        "",
        f"- Requested groups: {model['requested_cluster_count']:,}",
        f"- Populated groups: {model['populated_cluster_count']:,}",
        f"- TF-IDF features: {model['feature_count']:,}",
        f"- Latent dimensions: {model['latent_dimension_count']:,}",
        f"- Seed: {model['seed']}",
        "",
        "## Diagnostics",
        "",
        f"- Sampled silhouette score: {quality['sampled_silhouette_score']}",
        f"- Mean centroid similarity: {quality['mean_centroid_similarity']}",
        f"- Median assignment margin: {quality['median_assignment_margin']}",
        f"- Median exact cores per group: {quality['median_unique_cores_per_group']}",
        f"- Median material variants per group: {quality['median_unique_materials_per_group']}",
        f"- Median observations per group: {quality['median_observations_per_group']}",
        f"- Supported core assignments: {quality['core_assignment_status_counts']['supported']:,}",
        f"- Provisional core assignments: {quality['core_assignment_status_counts']['provisional']:,}",
        f"- Review core assignments: {quality['core_assignment_status_counts']['review']:,}",
        "",
        "## Caveats",
        "",
        "- `CGPOC2` identifiers are bound to this dataset snapshot and model settings.",
        "- Groups are statistical proposals, not reagent-interchangeability claims.",
        "- Additives are contextual in this POC even though some are mechanistically critical.",
        "- Production identities require chemistry review and frozen definitions.",
        "",
        "## Largest groups",
        "",
        "| Group | Rows | Cores | Material variants | Prototype core |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for group in report["largest_groups"]:
        prototype = str(group["prevalence_default_core_display"]).replace(
            "|", "\\|"
        )
        lines.append(
            f"| `{group['condition_group_id']}` | {group['observation_count']:,} | "
            f"{group['unique_core_count']:,} | {group['unique_material_count']:,} | "
            f"{prototype} |"
        )
    lines.append("")
    return "\n".join(lines)


def run_condition_grouping_poc(
    dataset_path: str | Path,
    output_dir: str | Path,
    *,
    cluster_count: int = 256,
    latent_dimensions: int = 48,
    seed: int = 17,
    silhouette_sample_size: int = 3000,
) -> ConditionGroupingRun:
    """Group decisive cores and write context-preserving review artifacts."""

    if cluster_count < 2:
        raise ValueError("cluster_count must be at least two")
    if latent_dimensions < 2:
        raise ValueError("latent_dimensions must be at least two")
    try:
        import joblib
        import numpy as np
        from sklearn.cluster import MiniBatchKMeans
        from sklearn.decomposition import TruncatedSVD
        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn.metrics import silhouette_score
        from sklearn.preprocessing import normalize
    except ImportError as exc:  # pragma: no cover - optional POC environment
        raise RuntimeError(
            "condition grouping POC requires scikit-learn, numpy, and joblib"
        ) from exc

    source = Path(dataset_path)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    materials, audit = load_weak_label_materials(source)
    cores, unresolved_materials = build_condition_cores(materials)
    if len(cores) <= cluster_count:
        raise ValueError(
            "cluster_count must be smaller than the number of exact cores"
        )

    vectorizer = TfidfVectorizer(
        analyzer=_identity_analyzer,
        lowercase=False,
        token_pattern=None,
        binary=True,
        norm="l2",
        min_df=2,
        dtype=np.float64,
    )
    matrix = vectorizer.fit_transform([core.feature_tokens for core in cores])
    dimension_count = min(latent_dimensions, matrix.shape[0] - 1, matrix.shape[1] - 1)
    if dimension_count < 2:
        raise ValueError("condition core features are insufficient for grouping")
    svd = TruncatedSVD(n_components=dimension_count, random_state=seed)
    latent = normalize(svd.fit_transform(matrix), norm="l2")
    model = MiniBatchKMeans(
        n_clusters=cluster_count,
        random_state=seed,
        batch_size=min(1024, len(cores)),
        n_init=10,
        max_iter=300,
        reassignment_ratio=0.0,
    )
    sample_weights = np.sqrt(
        np.asarray([core.observation_count for core in cores], dtype=float)
    )
    labels = model.fit_predict(latent, sample_weight=sample_weights)
    normalized_centers = normalize(model.cluster_centers_, norm="l2")
    similarities = latent @ normalized_centers.T
    assigned_similarity = similarities[np.arange(len(cores)), labels]
    sorted_similarity = np.sort(similarities, axis=1)
    margins = sorted_similarity[:, -1] - sorted_similarity[:, -2]

    by_label: Dict[int, list[int]] = defaultdict(list)
    for index, label in enumerate(labels.tolist()):
        by_label[int(label)].append(index)
    input_sha256 = _sha256_file(source)
    group_id_by_label: Dict[int, str] = {}
    group_records = []
    assignment_rows = []
    core_status_counts: Counter[str] = Counter()
    material_status_counts: Counter[str] = Counter()

    for label, indices in sorted(by_label.items()):
        medoid_index = sorted(
            indices,
            key=lambda index: (
                -float(assigned_similarity[index]),
                -cores[index].observation_count,
                cores[index].core_id,
            ),
        )[0]
        prototype = cores[medoid_index]
        default_index = sorted(
            indices,
            key=lambda index: (
                -cores[index].observation_count,
                -float(assigned_similarity[index]),
                cores[index].core_id,
            ),
        )[0]
        prevalence_default = cores[default_index]
        group_id = _digest(
            CONDITION_GROUP_ID_NAMESPACE,
            {
                "definition_version": CONDITION_GROUPING_POC_DEFINITION_VERSION,
                "input_sha256": input_sha256,
                "cluster_count": cluster_count,
                "latent_dimensions": dimension_count,
                "seed": seed,
                "prototype_core_id": prototype.core_id,
            },
        )
        group_id_by_label[label] = group_id
        member_cores = [cores[index] for index in indices]
        member_materials = [
            material for core in member_cores for material in core.materials
        ]
        observation_count = sum(
            material.observation_count for material in member_materials
        )
        recipe_ids = set().union(
            *(material.recipe_ids for material in member_materials)
        )
        reactions: Counter[str] = Counter()
        topologies: Counter[str] = Counter()
        yield_sum = yield_count = uncertain_count = 0
        for core in member_cores:
            topologies[_topology_key(core)] += core.observation_count
        for material in member_materials:
            reactions.update(material.reaction_types)
            yield_sum += material.yield_sum
            yield_count += material.yield_count
            uncertain_count += material.uncertain_count
        role_components = _role_component_counts(member_materials)
        temperatures = _counter_sum(member_materials, "temperatures")
        times = _counter_sum(member_materials, "times")
        member_similarities = [float(assigned_similarity[index]) for index in indices]
        member_margins = [float(margins[index]) for index in indices]
        representatives = tuple(
            {
                "condition_core_id": cores[index].core_id,
                "core_display": cores[index].core_display,
                "observation_count": cores[index].observation_count,
                "material_variant_count": len(cores[index].materials),
                "centroid_similarity": round(float(assigned_similarity[index]), 6),
                "assignment_margin": round(float(margins[index]), 6),
            }
            for index in sorted(
                indices,
                key=lambda index: (
                    -float(assigned_similarity[index]),
                    -cores[index].observation_count,
                    cores[index].core_id,
                ),
            )[:10]
        )
        core_role_components = {
            role: _top_items(role_components[role])
            for role in sorted(CORE_COMPONENT_ROLES)
            if role_components.get(role)
        }
        core_role_consensus = {
            role: {
                "value": ranked[0]["value"],
                "observation_count": ranked[0]["count"],
                "coverage_rate": round(
                    ranked[0]["count"] / observation_count, 6
                ),
            }
            for role, ranked in core_role_components.items()
            if ranked
        }
        contextual_roles = {
            role: _top_items(counter)
            for role, counter in sorted(role_components.items())
            if role not in CORE_COMPONENT_ROLES
        }
        group_records.append(
            {
                "schema_version": CONDITION_GROUPING_POC_SCHEMA_VERSION,
                "definition_version": CONDITION_GROUPING_POC_DEFINITION_VERSION,
                "condition_group_id": group_id,
                "snapshot_cluster_label": label,
                "statistical_prototype_core_id": prototype.core_id,
                "statistical_prototype_core_display": prototype.core_display,
                "prevalence_default_core_id": prevalence_default.core_id,
                "prevalence_default_core_display": prevalence_default.core_display,
                "prevalence_default_observation_count": (
                    prevalence_default.observation_count
                ),
                "unique_core_count": len(member_cores),
                "unique_material_count": len(member_materials),
                "exact_recipe_id_count": len(recipe_ids),
                "observation_count": observation_count,
                "mean_yield_pct": round(yield_sum / yield_count, 4)
                if yield_count
                else None,
                "yield_observation_count": yield_count,
                "condition_identity_uncertainty_rate": round(
                    uncertain_count / observation_count, 6
                ),
                "mean_centroid_similarity": round(
                    statistics.fmean(member_similarities), 6
                ),
                "p10_centroid_similarity": _safe_round(
                    _percentile(member_similarities, 0.1)
                ),
                "median_centroid_similarity": _safe_round(
                    _percentile(member_similarities, 0.5)
                ),
                "mean_assignment_margin": round(
                    statistics.fmean(member_margins), 6
                ),
                "core_topology_purity": round(
                    max(topologies.values(), default=0) / observation_count, 6
                ),
                "top_core_components": core_role_components,
                "core_role_consensus": core_role_consensus,
                "context_cross_references": {
                    "solvent_components": contextual_roles.get("solvent", ()),
                    "solvent_systems": _top_items(
                        _role_profile_counts(member_materials, "solvent")
                    ),
                    "additive_components": contextual_roles.get("additive", ()),
                    "additive_systems": _top_items(
                        _role_profile_counts(member_materials, "additive")
                    ),
                    "other_context_components": {
                        role: values
                        for role, values in contextual_roles.items()
                        if role not in MEDIUM_COMPONENT_ROLES
                        and role not in MODIFIER_COMPONENT_ROLES
                    },
                    "temperature_c": _numeric_cross_reference(
                        temperatures, observation_count
                    ),
                    "time_h": _numeric_cross_reference(times, observation_count),
                },
                "top_reaction_types": _top_items(reactions),
                "representative_cores": representatives,
                "most_common_cores": tuple(
                    {
                        "condition_core_id": cores[index].core_id,
                        "core_display": cores[index].core_display,
                        "observation_count": cores[index].observation_count,
                        "material_variant_count": len(cores[index].materials),
                    }
                    for index in sorted(
                        indices,
                        key=lambda index: (
                            -cores[index].observation_count,
                            cores[index].core_id,
                        ),
                    )[:10]
                ),
                "warnings": (
                    "EXPERIMENTAL_SNAPSHOT_BOUND_GROUP_ID",
                    "STATISTICAL_GROUP_NOT_INTERCHANGEABILITY_CLAIM",
                ),
            }
        )
        for index in indices:
            core = cores[index]
            status = _assignment_status(
                float(assigned_similarity[index]), float(margins[index])
            )
            core_status_counts[status] += 1
            for material in core.materials:
                material_status_counts[status] += 1
                other_context = tuple(
                    component.name
                    for component in material.parsed.context_components
                    if component.role not in MEDIUM_COMPONENT_ROLES
                    and component.role not in MODIFIER_COMPONENT_ROLES
                )
                assignment_rows.append(
                    {
                        "material_id": material.material_id,
                        "condition_core_id": core.core_id,
                        "condition_group_id": group_id,
                        "centroid_similarity": round(
                            float(assigned_similarity[index]), 6
                        ),
                        "assignment_margin": round(float(margins[index]), 6),
                        "assignment_status": status,
                        "observation_count": material.observation_count,
                        "exact_recipe_id_count": len(material.recipe_ids),
                        "core_display": core.core_display,
                        "solvent_cross_reference": " + ".join(
                            material.parsed.role_names("solvent")
                        ),
                        "additive_cross_reference": " + ".join(
                            material.parsed.role_names("additive")
                        ),
                        "other_context_cross_reference": " + ".join(other_context),
                        "complete_condition_display": (
                            material.parsed.canonical_display
                        ),
                        "assignment_reason": "LEARNED_DECISIVE_CORE_GROUP",
                        "parse_warnings": "|".join(material.parsed.warnings),
                    }
                )

    for material in unresolved_materials:
        material_status_counts["review"] += 1
        assignment_rows.append(
            {
                "material_id": material.material_id,
                "condition_core_id": "",
                "condition_group_id": "",
                "centroid_similarity": "",
                "assignment_margin": "",
                "assignment_status": "review",
                "observation_count": material.observation_count,
                "exact_recipe_id_count": len(material.recipe_ids),
                "core_display": "",
                "solvent_cross_reference": " + ".join(
                    material.parsed.role_names("solvent")
                ),
                "additive_cross_reference": " + ".join(
                    material.parsed.role_names("additive")
                ),
                "other_context_cross_reference": "",
                "complete_condition_display": material.parsed.canonical_display,
                "assignment_reason": "NO_DECISIVE_CORE_COMPONENT",
                "parse_warnings": "|".join(material.parsed.warnings),
            }
        )

    group_records.sort(
        key=lambda item: (-item["observation_count"], item["condition_group_id"])
    )
    assignment_rows.sort(key=lambda item: item["material_id"])
    sampled_size = min(max(0, silhouette_sample_size), len(cores))
    sampled_silhouette = None
    if sampled_size >= max(2, len(by_label)):
        sampled_silhouette = float(
            silhouette_score(
                latent,
                labels,
                metric="cosine",
                sample_size=sampled_size if sampled_size < len(cores) else None,
                random_state=seed,
            )
        )
    core_counts = [record["unique_core_count"] for record in group_records]
    material_counts = [record["unique_material_count"] for record in group_records]
    observation_counts = [record["observation_count"] for record in group_records]
    unresolved_observation_count = sum(
        material.observation_count for material in unresolved_materials
    )
    report: Dict[str, Any] = {
        "artifact_type": "chemist_oriented_condition_grouping_poc",
        "schema_version": CONDITION_GROUPING_POC_SCHEMA_VERSION,
        "definition_version": CONDITION_GROUPING_POC_DEFINITION_VERSION,
        "dataset": {
            "path": str(source.resolve()),
            "sha256": input_sha256,
            **audit,
            "unique_core_count": len(cores),
            "unresolved_core_material_count": len(unresolved_materials),
            "unresolved_core_observation_count": unresolved_observation_count,
        },
        "identity_factorization": {
            "material_identity": "all reported components and absences",
            "decisive_core_identity": {
                "roles": tuple(sorted(CORE_COMPONENT_ROLES)),
                "conditional_absence": "No ligand only when catalyst is present",
            },
            "medium_cross_reference_roles": tuple(sorted(MEDIUM_COMPONENT_ROLES)),
            "modifier_cross_reference_roles": tuple(
                sorted(MODIFIER_COMPONENT_ROLES)
            ),
            "operating_cross_references": ("temperature_c", "time_h"),
        },
        "feature_policy": {
            "included": (
                "decisive_core_component_names",
                "decisive_core_roles",
                "decisive_core_topology",
                "catalyst_ligand_pairs",
                "catalyst_base_pairs",
                "activation_reagent_base_pairs",
                "catalyst_conditioned_ligand_absence",
            ),
            "cross_referenced_not_clustered": (
                "solvents",
                "additives",
                "other_context_components",
                "temperature_c",
                "time_h",
            ),
            "excluded": (
                "source_reaction_type",
                "reactive_site_labels_and_signatures",
                "yield_pct",
                "z_score",
                "procedure_text",
            ),
        },
        "model": {
            "kind": "core_tfidf_svd_minibatch_kmeans",
            "requested_cluster_count": cluster_count,
            "populated_cluster_count": len(group_records),
            "feature_count": int(matrix.shape[1]),
            "latent_dimension_count": dimension_count,
            "seed": seed,
            "sample_weight": "sqrt_core_observation_count",
            "scikit_learn_version": __import__("sklearn").__version__,
        },
        "quality": {
            "sampled_silhouette_size": sampled_size,
            "sampled_silhouette_score": _safe_round(sampled_silhouette),
            "mean_centroid_similarity": round(
                float(np.mean(assigned_similarity)), 6
            ),
            "p10_centroid_similarity": _safe_round(
                _percentile(assigned_similarity.tolist(), 0.1)
            ),
            "median_assignment_margin": _safe_round(
                _percentile(margins.tolist(), 0.5)
            ),
            "p10_assignment_margin": _safe_round(
                _percentile(margins.tolist(), 0.1)
            ),
            "median_unique_cores_per_group": _safe_round(
                _percentile(core_counts, 0.5)
            ),
            "median_unique_materials_per_group": _safe_round(
                _percentile(material_counts, 0.5)
            ),
            "median_observations_per_group": _safe_round(
                _percentile(observation_counts, 0.5)
            ),
            "singleton_core_group_count": sum(value == 1 for value in core_counts),
            "core_assignment_status_counts": {
                status: core_status_counts[status]
                for status in ("supported", "provisional", "review")
            },
            "material_assignment_status_counts": {
                status: material_status_counts[status]
                for status in ("supported", "provisional", "review")
            },
            "assignment_status_policy": {
                "supported": {
                    "minimum_similarity": SUPPORTED_SIMILARITY_THRESHOLD,
                    "minimum_margin": SUPPORTED_MARGIN_THRESHOLD,
                },
                "provisional": {
                    "minimum_similarity": PROVISIONAL_SIMILARITY_THRESHOLD,
                    "minimum_margin": PROVISIONAL_MARGIN_THRESHOLD,
                },
                "review": "otherwise or no decisive core",
                "calibrated": False,
            },
        },
        "largest_groups": tuple(group_records[:20]),
        "warnings": (
            "POC_ONLY_NOT_PRODUCTION_CONDITION_IDENTITY",
            "GROUP_IDS_CHANGE_WHEN_DATA_OR_MODEL_SETTINGS_CHANGE",
            "WEAK_LABEL_SOURCE_COMPONENT_ROLES_NOT_INDEPENDENTLY_VALIDATED",
            "CONTEXT_EXCLUSION_DOES_NOT_IMPLY_CONTEXT_IS_CHEMICALLY_IRRELEVANT",
        ),
    }

    groups_path = destination / "condition_groups.jsonl"
    with groups_path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in group_records:
            handle.write(_canonical_json(record) + "\n")
    assignments_path = destination / "condition_group_assignments.csv"
    with assignments_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(assignment_rows[0]))
        writer.writeheader()
        writer.writerows(assignment_rows)
    model_path = destination / "condition_grouping_model.joblib"
    joblib.dump(
        {
            "schema_version": CONDITION_GROUPING_POC_SCHEMA_VERSION,
            "definition_version": CONDITION_GROUPING_POC_DEFINITION_VERSION,
            "input_sha256": input_sha256,
            "core_roles": tuple(sorted(CORE_COMPONENT_ROLES)),
            "vectorizer": vectorizer,
            "svd": svd,
            "cluster_model": model,
            "group_id_by_label": group_id_by_label,
        },
        model_path,
        compress=3,
    )
    report["artifacts"] = {
        "groups": groups_path.name,
        "assignments": assignments_path.name,
        "model": model_path.name,
        "model_sha256": _sha256_file(model_path),
    }
    report_json_path = destination / "condition_grouping_report.json"
    report_json_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report_markdown_path = destination / "condition_grouping_report.md"
    report_markdown_path.write_text(_render_report(report), encoding="utf-8")
    return ConditionGroupingRun(
        report=report,
        report_json_path=report_json_path,
        report_markdown_path=report_markdown_path,
        groups_path=groups_path,
        assignments_path=assignments_path,
        model_path=model_path,
    )


__all__ = [
    "CONDITION_GROUPING_POC_DEFINITION_VERSION",
    "CONDITION_GROUPING_POC_SCHEMA_VERSION",
    "CORE_COMPONENT_ROLES",
    "ConditionCoreAggregate",
    "ConditionGroupingRun",
    "MaterialAggregate",
    "ParsedConditionComponent",
    "ParsedConditionDisplay",
    "build_condition_cores",
    "load_weak_label_materials",
    "parse_condition_display",
    "run_condition_grouping_poc",
]
