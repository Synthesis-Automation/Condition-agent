"""Reference-disjoint POC for grammar-independent Fischer edit retrieval.

The Fischer source label selects the evaluation cohort only. It is never
passed to reaction featurization, signature construction, or retrieval.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable, Mapping


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from condition_recommender.edit_prototypes import (  # noqa: E402
    anonymous_edit_prototype,
    anonymous_edit_prototype_from_hypothesis,
    anonymous_edit_similarity,
)
from reactive_taxonomy import featurize_reaction  # noqa: E402
from reactive_taxonomy.chemistry.rdkit_utils import parse_smiles  # noqa: E402


POC_SCHEMA_VERSION = "1.1"
DEFAULT_SOURCE = (
    Path(__file__).parents[1]
    / "data-processor"
    / "reaction_dataset"
    / "Fischer_indole_synthesis.csv"
)


def _heavy_atom_count(smiles: str) -> int:
    molecule = parse_smiles(smiles)
    return int(molecule.GetNumHeavyAtoms()) if molecule is not None else 0


def _precursor_mode(reaction_smiles: str) -> str:
    left = reaction_smiles.split(">", 1)[0]
    substantial = sum(
        _heavy_atom_count(smiles) >= 2
        for smiles in left.split(".")
        if smiles
    )
    return "single_substrate" if substantial == 1 else "multi_component"


def _reference_split(reference: str) -> str:
    digest = hashlib.sha256(reference.strip().encode("utf-8")).digest()
    return "holdout" if int.from_bytes(digest[:2], "big") % 5 == 0 else "train"


def _read_unique_rows(path: Path) -> tuple[dict[str, str], ...]:
    rows = []
    seen = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            reaction_smiles = str(row.get("reaction_smiles") or "").strip()
            if not reaction_smiles or reaction_smiles in seen:
                continue
            seen.add(reaction_smiles)
            rows.append(dict(row))
    return tuple(rows)


def _source_row_count(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def _counter(values: Iterable[str]) -> dict[str, int]:
    return dict(sorted(Counter(values).items()))


def evaluate_fischer_edit_poc(
    source: str | Path = DEFAULT_SOURCE,
    *,
    max_unique_reactions: int | None = None,
    edit_graph_threshold: float = 0.6,
) -> dict[str, Any]:
    """Evaluate correspondence and name-free retrieval on one source cohort."""
    path = Path(source)
    rows = _read_unique_rows(path)
    if max_unique_reactions is not None:
        rows = rows[:max_unique_reactions]
    started = time.perf_counter()
    observations = []
    for row in rows:
        reaction_smiles = row["reaction_smiles"]
        analysis = featurize_reaction(reaction_smiles)
        signature = (
            asdict(analysis.reaction_signature)
            if analysis.reaction_signature is not None
            else None
        )
        prototype = (
            anonymous_edit_prototype(signature)
            if isinstance(signature, Mapping)
            else None
        )
        hypotheses = tuple(
            asdict(hypothesis)
            for hypothesis in analysis.edit_hypotheses
        )
        hypothesis_prototypes = tuple(
            prototype
            for hypothesis in hypotheses
            for prototype in (
                anonymous_edit_prototype_from_hypothesis(hypothesis),
            )
            if prototype is not None
        )
        reference = str(row.get("reference") or "").strip()
        observations.append(
            {
                "reaction_id": str(row.get("reaction_id") or ""),
                "reference": reference,
                "split": _reference_split(reference),
                "precursor_mode": _precursor_mode(reaction_smiles),
                "evidence_quality": analysis.evidence_quality,
                "signature": signature,
                "prototype": prototype,
                "edit_hypotheses": hypotheses,
                "hypothesis_prototypes": hypothesis_prototypes,
                "warnings": tuple(analysis.warnings),
            }
        )

    train = tuple(
        observation
        for observation in observations
        if observation["split"] == "train"
        and observation["prototype"] is not None
    )
    holdout = tuple(
        observation
        for observation in observations
        if observation["split"] == "holdout"
        and observation["prototype"] is not None
    )
    retrieval = Counter()
    retrieval_examples: dict[str, list[str]] = {
        "exact_l3": [],
        "edit_graph": [],
        "cross_mode_edit_graph": [],
        "no_match": [],
    }
    for query in holdout:
        exact = tuple(
            candidate
            for candidate in train
            if candidate["signature"]["bond_edit_signature_key"]
            == query["signature"]["bond_edit_signature_key"]
        )
        approximate = tuple(
            candidate
            for candidate in train
            if anonymous_edit_similarity(
                query["prototype"],
                candidate["prototype"],
            )
            >= edit_graph_threshold
        )
        cross_mode = tuple(
            candidate
            for candidate in approximate
            if candidate["precursor_mode"] != query["precursor_mode"]
        )
        if exact:
            retrieval["exact_l3"] += 1
            retrieval_examples["exact_l3"].append(query["reaction_id"])
        if approximate:
            retrieval["edit_graph"] += 1
            retrieval_examples["edit_graph"].append(query["reaction_id"])
        if cross_mode:
            retrieval["cross_mode_edit_graph"] += 1
            retrieval_examples["cross_mode_edit_graph"].append(
                query["reaction_id"]
            )
        if not approximate:
            retrieval["no_match"] += 1
            retrieval_examples["no_match"].append(query["reaction_id"])

    signatures = tuple(
        observation for observation in observations if observation["signature"]
    )
    reference_disjoint_retrieval = Counter()
    for query in signatures:
        independent = tuple(
            candidate
            for candidate in signatures
            if candidate["reference"] != query["reference"]
        )
        exact = tuple(
            candidate
            for candidate in independent
            if candidate["signature"]["bond_edit_signature_key"]
            == query["signature"]["bond_edit_signature_key"]
        )
        approximate = tuple(
            candidate
            for candidate in independent
            if anonymous_edit_similarity(
                query["prototype"],
                candidate["prototype"],
            )
            >= edit_graph_threshold
        )
        cross_mode = tuple(
            candidate
            for candidate in approximate
            if candidate["precursor_mode"] != query["precursor_mode"]
        )
        if exact:
            reference_disjoint_retrieval["exact_l3"] += 1
        if approximate:
            reference_disjoint_retrieval["edit_graph"] += 1
        if cross_mode:
            reference_disjoint_retrieval["cross_mode_edit_graph"] += 1
        if not approximate:
            reference_disjoint_retrieval["no_match"] += 1
    named_families = tuple(
        str(observation["signature"].get("named_family") or "")
        for observation in signatures
        if observation["signature"].get("named_family")
    )
    warning_counts = Counter(
        warning
        for observation in observations
        for warning in observation["warnings"]
    )
    ambiguous_queries = tuple(
        observation
        for observation in observations
        if observation["edit_hypotheses"]
    )
    hypothesis_retrieval = Counter()
    for query in ambiguous_queries:
        prototypes = query["hypothesis_prototypes"]
        if len(prototypes) != len(query["edit_hypotheses"]):
            hypothesis_retrieval["prototype_incomplete"] += 1
            continue
        independent = tuple(
            candidate
            for candidate in signatures
            if candidate["reference"] != query["reference"]
        )
        robust = tuple(
            candidate
            for candidate in independent
            if all(
                anonymous_edit_similarity(
                    prototype,
                    candidate["prototype"],
                )
                >= edit_graph_threshold
                for prototype in prototypes
            )
        )
        cross_mode = tuple(
            candidate
            for candidate in robust
            if candidate["precursor_mode"] != query["precursor_mode"]
        )
        if robust:
            hypothesis_retrieval["consensus_edit_graph"] += 1
        else:
            hypothesis_retrieval["no_consensus_match"] += 1
        if cross_mode:
            hypothesis_retrieval["cross_mode_consensus_edit_graph"] += 1
    return {
        "schema_version": POC_SCHEMA_VERSION,
        "artifact_type": "grammar_independent_edit_poc",
        "source": str(path),
        "source_label_role": "evaluation_cohort_only",
        "source_label_used_for_routing": False,
        "reaction_name_used_in_signature_identity": False,
        "row_count": _source_row_count(path),
        "unique_reaction_count": len(observations),
        "unique_reference_count": len(
            {observation["reference"] for observation in observations}
        ),
        "split_counts": _counter(
            observation["split"] for observation in observations
        ),
        "precursor_mode_counts": _counter(
            observation["precursor_mode"] for observation in observations
        ),
        "evidence_quality_counts": _counter(
            observation["evidence_quality"] for observation in observations
        ),
        "signature_count": len(signatures),
        "signature_coverage": round(
            len(signatures) / max(len(observations), 1),
            6,
        ),
        "fragmented_signature_count": sum(
            observation["evidence_quality"]
            == "fragmented_scaffold_correspondence"
            for observation in signatures
        ),
        "ambiguous_correspondence_count": sum(
            observation["evidence_quality"]
            == "ambiguous_atom_correspondence"
            for observation in observations
        ),
        "edit_hypothesis_query_count": len(ambiguous_queries),
        "fully_prototyped_edit_hypothesis_query_count": sum(
            len(observation["hypothesis_prototypes"])
            == len(observation["edit_hypotheses"])
            for observation in ambiguous_queries
        ),
        "shared_prototype_edit_hypothesis_query_count": sum(
            len(
                {
                    prototype.prototype_id
                    for prototype in observation["hypothesis_prototypes"]
                }
            )
            == 1
            for observation in ambiguous_queries
            if observation["hypothesis_prototypes"]
        ),
        "reference_disjoint_hypothesis_retrieval_counts": dict(
            sorted(hypothesis_retrieval.items())
        ),
        "named_family_count": len(named_families),
        "anonymous_prototype_count": len(
            {
                observation["prototype"].prototype_id
                for observation in signatures
                if observation["prototype"] is not None
            }
        ),
        "holdout_signed_query_count": len(holdout),
        "retrieval_counts": dict(sorted(retrieval.items())),
        "retrieval_examples": {
            key: values[:10] for key, values in retrieval_examples.items()
        },
        "reference_disjoint_signed_query_count": len(signatures),
        "reference_disjoint_retrieval_counts": dict(
            sorted(reference_disjoint_retrieval.items())
        ),
        "top_warnings": dict(warning_counts.most_common(15)),
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-unique-reactions", type=int)
    args = parser.parse_args()
    report = evaluate_fischer_edit_poc(
        args.source,
        max_unique_reactions=args.max_unique_reactions,
    )
    payload = json.dumps(report, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    print(payload)


if __name__ == "__main__":
    main()
