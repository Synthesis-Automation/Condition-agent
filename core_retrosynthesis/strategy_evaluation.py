"""Frozen, evidence-grouped evaluation of single-step strategy search.

Condition-index observations supply the corpus and shared leakage grouping.
Only the training partition supplies executable operators to either engine.
"""

from __future__ import annotations

import hashlib
import json
import random
from collections import defaultdict
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable

from condition_recommender.evaluation import grouped_holdout_split, _leakage_summary
from condition_recommender.evaluation_panel import _tokens
from condition_recommender.generic_indexing import load_generic_index
from .chemistry import canonical_smiles
from .generic_compiler import analyze_generic_reaction, compile_generic_templates
from .generic_library import build_generic_library, save_generic_library
from .generic_search import disconnect_operator_ladder_detailed
from .strategy_search import disconnect_strategies_detailed
from .strategy_identity import build_strategy_id


STRATEGY_EVALUATION_VERSION = "single_step_strategy_evaluation.v2"


def freeze_strategy_panel(
    source_path: str | Path,
    output: str | Path,
    *,
    size: int = 1500,
    query_limit: int = 60,
    seed: int = 2026091401,
    exclusion_paths: Iterable[str | Path] = (),
) -> Path:
    """Freeze raw inputs and grouped splits without examining search outcomes."""

    if size < 2 or query_limit < 1:
        raise ValueError("positive query limit and at least two observations required")
    exclusion_paths = tuple(exclusion_paths)
    destination = Path(output)
    path = destination / "panel.json"
    if path.exists():
        raise FileExistsError("refusing to overwrite frozen strategy panel")
    excluded = set()
    for previous_path in exclusion_paths:
        for row in load_generic_index(previous_path).rows:
            excluded.update(_tokens(row))
    index = load_generic_index(source_path)
    if size > len(index.rows):
        raise ValueError("panel exceeds source corpus")
    positions = random.Random(seed).sample(range(len(index.rows)), len(index.rows))
    selected = []
    for offset in range(0, len(positions), 500):
        for row in index.select(positions[offset : offset + 500]):
            if not (_tokens(row) & excluded):
                selected.append(row)
            if len(selected) == size:
                break
        if len(selected) == size:
            break
    if len(selected) != size:
        raise ValueError("insufficient observations after exclusions")
    split = grouped_holdout_split(selected, seed=seed, test_fraction=0.2)
    # Evaluate diverse structural annotations without using them as routing keys.
    strata = defaultdict(list)
    for row in split.test_rows:
        strata[row.transformation_class or "unannotated"].append(row)
    for values in strata.values():
        values.sort(
            key=lambda row: hashlib.sha256(
                f"{seed}:{row.observation_id}".encode()
            ).digest()
        )
    queries = []
    while any(strata.values()) and len(queries) < query_limit:
        for key in sorted(strata):
            if strata[key] and len(queries) < query_limit:
                queries.append(strata[key].pop(0).observation_id)

    def record(row: Any) -> dict[str, Any]:
        return {
            "reaction_id": row.reaction_id,
            "observation_id": row.observation_id,
            "reaction_smiles": row.reaction_smiles,
            "reference_id": row.reference_id,
            "transformation_class": row.transformation_class,
            "source_dataset": row.source_dataset,
        }

    payload = {
        "definition_id": STRATEGY_EVALUATION_VERSION,
        "seed": seed,
        "source_path": str(source_path),
        "source_count": len(index.rows),
        "exclusion_paths": [str(value) for value in exclusion_paths],
        "train": [record(row) for row in split.train_rows],
        "test": [record(row) for row in split.test_rows],
        "query_ids": queries,
        "leakage": _leakage_summary(split),
        "protocol": {
            "top_k": 5,
            "max_templates_per_level": 100,
            "max_validations_per_level": 30,
            "engines": ["flat", "strategy"],
            "no_tuning_after_outcomes": True,
        },
    }
    destination.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def evaluate_strategy_panel(panel_path: str | Path) -> dict[str, Any]:
    """Search valid products and report reference availability separately."""

    panel_path = Path(panel_path)
    output = panel_path.parent
    execution = output / "execution.json"
    if execution.exists():
        raise FileExistsError("frozen strategy evaluation already started")
    panel = json.loads(panel_path.read_text(encoding="utf-8"))
    root = Path(__file__).resolve().parents[1]
    hashes = {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for package in ("core_retrosynthesis", "reactive_taxonomy")
        for path in (root / package).rglob("*")
        if path.suffix in {".py", ".json"}
    }
    execution.write_text(
        json.dumps(
            {
                "panel_sha256": hashlib.sha256(panel_path.read_bytes()).hexdigest(),
                "code_sha256": hashes,
                "evaluation_definition_id": STRATEGY_EVALUATION_VERSION,
                "query_policy": "search_all_valid_single_products",
                "recovery_policy": "report_available_reference_denominators",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    library = build_generic_library(
        panel["train"], levels=("L2", "L1", "L0"), admission_mode="data_driven"
    )
    save_generic_library(library, output / "training_operators.json.gz")
    selected = set(panel["query_ids"])
    cases = []
    for row in panel["test"]:
        if row["observation_id"] not in selected:
            continue
        compilation = compile_generic_templates(
            row, levels=("L2", "L1", "L0"), admission_mode="data_driven"
        )
        case = {
            "observation_id": row["observation_id"],
            "reaction_smiles": row["reaction_smiles"],
            "transformation_class": row["transformation_class"],
            "source_rejection": compilation.rejection_reason,
            "engines": {},
        }
        identity = None
        expected_strategy = ""
        expected_precursors = None
        sides = row["reaction_smiles"].split(">")
        target = canonical_smiles(sides[2]) if len(sides) == 3 else None
        if not target or "." in target:
            case["target_rejection"] = "single_valid_product_required"
            target = None
        if compilation.templates:
            template = compilation.templates[0]
            precedent = template.precedents[0]
            identity = analyze_generic_reaction(precedent.mapped_reaction_smiles)
            expected_strategy = (
                build_strategy_id(
                    template.operator_id,
                    identity.disconnection_site_key,
                    identity.synthon_signature,
                )
                if identity
                and identity.disconnection_site_key
                and identity.synthon_signature
                else ""
            )
            expected_precursors = precedent.precursor_smiles
        case["target_smiles"] = target or ""
        case["expected_strategy_id"] = expected_strategy
        case["source_reference_available"] = bool(compilation.templates)
        if target:
            for engine in ("flat", "strategy"):
                start = perf_counter()
                options = dict(
                    max_templates_to_apply=100, max_candidates_to_validate=30
                )
                if engine == "flat":
                    candidates, diagnostics = disconnect_operator_ladder_detailed(
                        target, library, top_k=5, **options
                    )
                    serialized = {
                        "candidates": [c.to_dict() for c in candidates],
                        "search_diagnostics": diagnostics.to_dict(),
                    }
                else:
                    result = disconnect_strategies_detailed(
                        target, library, top_k_strategies=5, **options
                    )
                    candidates = tuple(
                        c
                        for strategy in result.strategies
                        for c in strategy.realizations
                    )
                    diagnostics = result.diagnostics
                    serialized = result.to_dict()
                case["engines"][engine] = {
                    "seconds": perf_counter() - start,
                    "strategy_count": len(
                        {c.strategy_id for c in candidates if c.strategy_id}
                    ),
                    "realization_count": len(candidates),
                    "strategy_recovered": (
                        any(c.strategy_id == expected_strategy for c in candidates)
                        if expected_strategy else None
                    ),
                    "site_recovered": (
                        any(
                            c.disconnection_site_key == identity.disconnection_site_key
                            for c in candidates
                        )
                        if identity and identity.disconnection_site_key else None
                    ),
                    "exact_precursors_recovered": (
                        any(c.precursor_smiles == expected_precursors for c in candidates)
                        if expected_precursors is not None else None
                    ),
                    "invalid_returned_count": sum(
                        c.forward_validation_status != "verified_signature"
                        for c in candidates
                    ),
                    "validation_attempts": diagnostics.validation_attempt_count,
                    "result": serialized,
                }
        cases.append(case)
        (output / "progress.json").write_text(
            json.dumps(
                {"completed_queries": len(cases), "total_queries": len(selected)}
            ),
            encoding="utf-8",
        )
    summary = {}
    for engine in ("flat", "strategy"):
        values = [case["engines"].get(engine, {}) for case in cases]
        denominator = max(1, len(values))
        summary[engine] = {
            "coverage": sum(bool(v.get("realization_count")) for v in values)
            / denominator,
            **{
                key: sum(v.get(key, 0) for v in values) / denominator
                for key in (
                    "strategy_count",
                    "validation_attempts",
                    "seconds",
                )
            },
            "invalid_returned_count": sum(
                v.get("invalid_returned_count", 0) for v in values
            ),
        }
        recovery_keys = (
            "strategy_recovered", "site_recovered", "exact_precursors_recovered",
        )
        recovery_denominators = {
            key: sum(v.get(key) is not None for v in values)
            for key in recovery_keys
        }
        summary[engine]["recovery_denominators"] = recovery_denominators
        for key, count in recovery_denominators.items():
            summary[engine][key] = (
                sum(bool(v.get(key)) for v in values) / count if count else None
            )
    report = {
        "definition_id": STRATEGY_EVALUATION_VERSION,
        "panel_definition_id": panel["definition_id"],
        "query_count": len(cases),
        "source_compilation_failures": sum(
            bool(case["source_rejection"]) for case in cases
        ),
        "target_query_failures": sum(not case["target_smiles"] for case in cases),
        "training_source_rows": library.source_row_count,
        "training_accepted_rows": library.accepted_observation_count,
        "leakage": panel["leakage"],
        "summary": summary,
        "cases": cases,
    }
    (output / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
