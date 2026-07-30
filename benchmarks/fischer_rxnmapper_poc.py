"""Evaluate RXNMapper as untrusted atom-correspondence evidence on Fischer data.

The Fischer source label selects the evaluation cohort only.  It is not passed
to RXNMapper, reaction featurization, edit normalization, or comparison.
Mapper output remains review evidence throughout this POC.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import statistics
import sys
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy import (  # noqa: E402
    ExternalAtomMappingResult,
    ReactionAtomReference,
    ReactionComponent,
    ReactionEdit,
    RxnMapperProvider,
    featurize_reaction,
)
from reactive_taxonomy.chemistry.rdkit_utils import parse_smiles  # noqa: E402
from reactive_taxonomy.reaction_parser import parse_reaction_smiles  # noqa: E402


POC_SCHEMA_VERSION = "1.0"
DEFAULT_SOURCE = (
    PROJECT_ROOT
    / "data-processor"
    / "reaction_dataset"
    / "Fischer_indole_synthesis.csv"
)
DEFAULT_OUTPUT_DIRECTORY = PROJECT_ROOT / "results" / "fischer_rxnmapper_poc"


def _read_rows(path: Path) -> tuple[dict[str, str], ...]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return tuple(dict(row) for row in csv.DictReader(handle))


def _unique_rows(
    rows: Sequence[dict[str, str]],
) -> tuple[dict[str, str], ...]:
    unique = {}
    for row in rows:
        reaction = str(row.get("reaction_smiles") or "").strip()
        if reaction:
            unique.setdefault(reaction, row)
    return tuple(unique.values())


def _counter(values: Iterable[str]) -> dict[str, int]:
    return dict(sorted(Counter(values).items()))


def _canonical_without_maps(smiles: str) -> str:
    from rdkit import Chem

    molecule = parse_smiles(smiles)
    if molecule is None:
        return smiles
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _component(
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


def _endpoint_token(
    atom: Optional[ReactionAtomReference],
    reactants: Sequence[ReactionComponent],
    products: Sequence[ReactionComponent],
) -> tuple[Any, ...]:
    if atom is None:
        return ("H",)
    components = products if atom.side == "product" else reactants
    component = _component(components, atom.component_index)
    return (
        atom.side,
        _canonical_without_maps(component.input_smiles) if component else "?",
        atom.element,
        atom.formal_charge,
        atom.aromatic,
        atom.hybridization,
        atom.local_environment_id,
    )


def _detailed_edit_profile(
    edits: Sequence[ReactionEdit],
    reactants: Sequence[ReactionComponent],
    products: Sequence[ReactionComponent],
) -> tuple[tuple[Any, ...], ...]:
    tokens = []
    for edit in edits:
        endpoints = tuple(
            sorted(
                (
                    _endpoint_token(edit.atom_1, reactants, products),
                    _endpoint_token(edit.atom_2, reactants, products),
                ),
                key=lambda value: json.dumps(value, sort_keys=True),
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
    return tuple(
        sorted(tokens, key=lambda value: json.dumps(value, sort_keys=True))
    )


def _bond_type(edit: ReactionEdit, order: Optional[str]) -> str:
    if edit.atom_2 is None:
        return f"{edit.atom_1.element}-H:{order or 'NONE'}"
    elements = sorted((edit.atom_1.element, edit.atom_2.element))
    return f"{elements[0]}-{elements[1]}:{order or 'NONE'}"


def _coarse_edit_profile(edits: Sequence[ReactionEdit]) -> dict[str, tuple[str, ...]]:
    return {
        "formed_bond_types": tuple(
            sorted(
                _bond_type(edit, edit.new_order)
                for edit in edits
                if edit.edit_type == "formed"
            )
        ),
        "broken_bond_types": tuple(
            sorted(
                _bond_type(edit, edit.old_order)
                for edit in edits
                if edit.edit_type == "broken"
            )
        ),
        "order_changes": tuple(
            sorted(
                f"{_bond_type(edit, edit.old_order)}>{edit.new_order or 'NONE'}"
                for edit in edits
                if edit.edit_type == "order_changed"
            )
        ),
        "hydrogen_changes": tuple(
            sorted(
                f"{_bond_type(edit, edit.old_order)}>{edit.new_order or 'NONE'}"
                for edit in edits
                if edit.edit_type == "hydrogen_change"
            )
        ),
    }


def _signature_coarse_profile(signature: Any) -> dict[str, tuple[str, ...]]:
    return {
        "formed_bond_types": tuple(signature.formed_bond_types),
        "broken_bond_types": tuple(signature.broken_bond_types),
        "order_changes": tuple(signature.order_changes),
        "hydrogen_changes": tuple(signature.hydrogen_changes),
    }


def _profile_id(profile: Any) -> str:
    canonical = json.dumps(profile, sort_keys=True, separators=(",", ":"))
    return "EMP1:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:24]


def _confidence_summary(values: Sequence[float]) -> dict[str, Optional[float]]:
    if not values:
        return {
            "count": 0,
            "minimum": None,
            "p10": None,
            "median": None,
            "p90": None,
            "maximum": None,
            "mean": None,
        }
    ordered = sorted(values)

    def percentile(fraction: float) -> float:
        index = round((len(ordered) - 1) * fraction)
        return round(ordered[index], 6)

    return {
        "count": len(ordered),
        "minimum": round(ordered[0], 6),
        "p10": percentile(0.10),
        "median": round(statistics.median(ordered), 6),
        "p90": percentile(0.90),
        "maximum": round(ordered[-1], 6),
        "mean": round(statistics.fmean(ordered), 6),
    }


def _reverse_reactant_order(reaction: str) -> str:
    if ">>" in reaction:
        left, right = reaction.split(">>", 1)
        return f"{'.'.join(reversed(left.split('.')))}>>{right}"
    left, middle, right = reaction.split(">")
    return f"{'.'.join(reversed(left.split('.')))}>{middle}>{right}"


def _alternate_component_serialization(reaction: str) -> str:
    from rdkit import Chem

    def alternate_side(side: str) -> str:
        output = []
        for token in side.split("."):
            molecule = parse_smiles(token)
            if molecule is None or molecule.GetNumAtoms() < 2:
                output.append(token)
                continue
            reversed_molecule = Chem.RenumberAtoms(
                molecule,
                tuple(reversed(range(molecule.GetNumAtoms()))),
            )
            output.append(
                Chem.MolToSmiles(
                    reversed_molecule,
                    canonical=False,
                    isomericSmiles=True,
                )
            )
        return ".".join(output)

    if ">>" in reaction:
        left, right = reaction.split(">>", 1)
        return f"{alternate_side(left)}>>{alternate_side(right)}"
    left, middle, right = reaction.split(">")
    return (
        f"{alternate_side(left)}>{alternate_side(middle)}>"
        f"{alternate_side(right)}"
    )


def _map_with_progress(
    provider: RxnMapperProvider,
    reactions: Sequence[str],
    *,
    chunk_size: int,
    progress: bool,
) -> tuple[ExternalAtomMappingResult, ...]:
    results = []
    for offset in range(0, len(reactions), chunk_size):
        batch = reactions[offset : offset + chunk_size]
        results.extend(provider.map_reactions(batch))
        if progress:
            print(
                f"mapped {min(offset + len(batch), len(reactions))}/"
                f"{len(reactions)}",
                file=sys.stderr,
                flush=True,
            )
    return tuple(results)


def evaluate_fischer_rxnmapper_poc(
    source: str | Path = DEFAULT_SOURCE,
    *,
    batch_size: int = 16,
    stability_sample_size: int = 24,
    max_unique_reactions: Optional[int] = None,
    progress: bool = False,
) -> tuple[dict[str, Any], tuple[dict[str, Any], ...]]:
    """Return the POC summary and one auditable record per unique reaction."""
    source_path = Path(source)
    rows = _read_rows(source_path)
    all_unique_rows = _unique_rows(rows)
    unique_rows = all_unique_rows
    if max_unique_reactions is not None:
        unique_rows = unique_rows[:max_unique_reactions]
    reactions = tuple(row["reaction_smiles"].strip() for row in unique_rows)
    duplicate_counts = Counter(
        str(row.get("reaction_smiles") or "").strip() for row in rows
    )
    started = time.perf_counter()
    analyses = []
    for index, row in enumerate(unique_rows, start=1):
        analyses.append(featurize_reaction(row["reaction_smiles"]))
        if progress and index % 25 == 0:
            print(
                f"analyzed {index}/{len(unique_rows)} existing outcomes",
                file=sys.stderr,
                flush=True,
            )
    provider = RxnMapperProvider(batch_size=batch_size)
    mapping_started = time.perf_counter()
    mapping_results = _map_with_progress(
        provider,
        reactions,
        chunk_size=batch_size,
        progress=progress,
    )
    mapping_elapsed = time.perf_counter() - mapping_started

    records = []
    agreement_counts: Counter[str] = Counter()
    disposition_counts: Counter[str] = Counter()
    valid_confidences = []
    for row, analysis, mapping in zip(unique_rows, analyses, mapping_results):
        external_coarse = None
        external_detailed = None
        mapped_parsed = None
        if mapping.valid and mapping.normalization is not None:
            valid_confidences.append(float(mapping.mapper_confidence or 0.0))
            mapped_parsed = parse_reaction_smiles(
                str(mapping.mapped_reaction_smiles)
            )
            external_coarse = _coarse_edit_profile(
                mapping.normalization.edits
            )
            external_detailed = _detailed_edit_profile(
                mapping.normalization.edits,
                mapped_parsed.reactants,
                mapped_parsed.products,
            )

        signature_coarse_agreement = None
        signature_detailed_agreement = None
        if analysis.reaction_signature is not None and external_coarse is not None:
            agreement_counts["existing_signature_evaluable"] += 1
            signature_coarse_agreement = (
                external_coarse
                == _signature_coarse_profile(analysis.reaction_signature)
            )
            signature_detailed_agreement = (
                external_detailed
                == _detailed_edit_profile(
                    analysis.reaction_signature.edits,
                    analysis.reactants,
                    analysis.products,
                )
            )
            agreement_counts[
                "existing_signature_coarse_agreement"
                if signature_coarse_agreement
                else "existing_signature_coarse_disagreement"
            ] += 1
            agreement_counts[
                "existing_signature_detailed_agreement"
                if signature_detailed_agreement
                else "existing_signature_detailed_disagreement"
            ] += 1

        hypothesis_match_count = None
        if analysis.edit_hypotheses and external_detailed is not None:
            agreement_counts["ambiguous_hypothesis_evaluable"] += 1
            hypothesis_profiles = tuple(
                _detailed_edit_profile(
                    hypothesis.edits,
                    analysis.reactants,
                    analysis.products,
                )
                for hypothesis in analysis.edit_hypotheses
            )
            hypothesis_match_count = sum(
                profile == external_detailed for profile in hypothesis_profiles
            )
            if hypothesis_match_count == 1:
                agreement_counts["unique_hypothesis_match"] += 1
            elif hypothesis_match_count > 1:
                agreement_counts["multiple_hypothesis_matches"] += 1
            else:
                agreement_counts["no_hypothesis_match"] += 1

        if not mapping.valid:
            disposition = "invalid_mapper_result"
        elif signature_detailed_agreement:
            disposition = "agrees_with_existing_signature"
        elif hypothesis_match_count == 1:
            disposition = "selects_one_hypothesis_review_required"
        elif analysis.edit_hypotheses:
            disposition = "ambiguous_or_novel_mapper_edits_review_required"
        else:
            disposition = "mapper_only_review_required"
        disposition_counts[disposition] += 1

        records.append(
            {
                "reaction_id": str(row.get("reaction_id") or ""),
                "reference": str(row.get("reference") or ""),
                "reaction_smiles": row["reaction_smiles"],
                "duplicate_row_count": duplicate_counts[row["reaction_smiles"]],
                "existing_evidence_quality": analysis.evidence_quality,
                "existing_signature_id": (
                    analysis.reaction_signature.signature_id
                    if analysis.reaction_signature is not None
                    else None
                ),
                "existing_edit_hypothesis_count": len(
                    analysis.edit_hypotheses
                ),
                "mapper_valid": mapping.valid,
                "mapper_confidence": mapping.mapper_confidence,
                "structure_preserved": mapping.structure_preserved,
                "reactant_mapping_coverage": mapping.reactant_mapping_coverage,
                "product_mapping_coverage": mapping.product_mapping_coverage,
                "shared_mapped_heavy_atom_count": (
                    mapping.shared_mapped_heavy_atom_count
                ),
                "mapped_reaction_smiles": mapping.mapped_reaction_smiles,
                "external_edit_profile": external_coarse,
                "external_edit_profile_id": (
                    _profile_id(external_coarse)
                    if external_coarse is not None
                    else None
                ),
                "normalized_edits": (
                    [asdict(edit) for edit in mapping.normalization.edits]
                    if mapping.normalization is not None
                    else []
                ),
                "signature_coarse_agreement": signature_coarse_agreement,
                "signature_detailed_agreement": (
                    signature_detailed_agreement
                ),
                "hypothesis_match_count": hypothesis_match_count,
                "review_disposition": disposition,
                "warnings": list(mapping.warnings),
                "error": mapping.error,
            }
        )

    stability_rows = tuple(
        sorted(
            unique_rows,
            key=lambda row: hashlib.sha256(
                row["reaction_smiles"].encode("utf-8")
            ).hexdigest(),
        )[:stability_sample_size]
    )
    base_by_reaction = {
        record["reaction_smiles"]: record for record in records
    }
    order_variants = tuple(
        _reverse_reactant_order(row["reaction_smiles"])
        for row in stability_rows
    )
    serialization_variants = tuple(
        _alternate_component_serialization(row["reaction_smiles"])
        for row in stability_rows
    )
    variant_results = _map_with_progress(
        provider,
        order_variants + serialization_variants,
        chunk_size=batch_size,
        progress=progress,
    )
    order_results = variant_results[: len(stability_rows)]
    serialization_results = variant_results[len(stability_rows) :]

    def stability_count(
        variant_values: Sequence[ExternalAtomMappingResult],
    ) -> tuple[int, int]:
        valid_count = 0
        stable_count = 0
        for row, variant in zip(stability_rows, variant_values):
            base = base_by_reaction[row["reaction_smiles"]]
            if not variant.valid or variant.normalization is None:
                continue
            valid_count += 1
            if (
                _coarse_edit_profile(variant.normalization.edits)
                == base["external_edit_profile"]
            ):
                stable_count += 1
        return valid_count, stable_count

    order_valid, order_stable = stability_count(order_results)
    serialization_valid, serialization_stable = stability_count(
        serialization_results
    )
    valid_mappings = tuple(result for result in mapping_results if result.valid)
    full_product_mapping_count = sum(
        result.product_mapping_coverage == 1.0 for result in mapping_results
    )
    structure_preserved_count = sum(
        result.structure_preserved for result in mapping_results
    )
    weighted_valid_rows = sum(
        duplicate_counts[result.input_reaction_smiles]
        for result in valid_mappings
    )
    report = {
        "schema_version": POC_SCHEMA_VERSION,
        "artifact_type": "fischer_rxnmapper_correspondence_poc",
        "source": str(source_path),
        "source_label_role": "evaluation_cohort_only",
        "source_label_used_for_mapping_or_edit_identity": False,
        "mapper_output_used_for_admission": False,
        "row_count": len(rows),
        "source_unique_reaction_count": len(all_unique_rows),
        "analyzed_unique_reaction_count": len(unique_rows),
        "provider": asdict(provider.metadata),
        "existing_system": {
            "unique_evidence_quality_counts": _counter(
                analysis.evidence_quality for analysis in analyses
            ),
            "unique_signature_count": sum(
                analysis.reaction_signature is not None
                for analysis in analyses
            ),
            "unique_ambiguous_correspondence_count": sum(
                analysis.evidence_quality == "ambiguous_atom_correspondence"
                for analysis in analyses
            ),
        },
        "external_mapping": {
            "valid_unique_reaction_count": len(valid_mappings),
            "valid_row_count": weighted_valid_rows,
            "structure_preserved_unique_reaction_count": (
                structure_preserved_count
            ),
            "full_product_mapping_unique_reaction_count": (
                full_product_mapping_count
            ),
            "failure_count": len(mapping_results) - len(valid_mappings),
            "confidence": _confidence_summary(valid_confidences),
            "confidence_threshold_counts": {
                "gte_0_5": sum(value >= 0.5 for value in valid_confidences),
                "gte_0_7": sum(value >= 0.7 for value in valid_confidences),
                "gte_0_9": sum(value >= 0.9 for value in valid_confidences),
            },
        },
        "agreement": dict(sorted(agreement_counts.items())),
        "review_disposition_counts": dict(sorted(disposition_counts.items())),
        "stability": {
            "sample_size": len(stability_rows),
            "reactant_order_valid_count": order_valid,
            "reactant_order_stable_edit_profile_count": order_stable,
            "alternate_serialization_valid_count": serialization_valid,
            "alternate_serialization_stable_edit_profile_count": (
                serialization_stable
            ),
        },
        "representative_reaction": next(
            (
                {
                    key: record[key]
                    for key in (
                        "reaction_id",
                        "reaction_smiles",
                        "existing_evidence_quality",
                        "mapper_confidence",
                        "product_mapping_coverage",
                        "external_edit_profile",
                        "hypothesis_match_count",
                        "review_disposition",
                        "mapped_reaction_smiles",
                    )
                }
                for record in records
                if record["reaction_smiles"]
                == (
                    "O=C1CCCCC1.Cl.NNc1ccc(F)cc1"
                    ">>Fc1ccc2[nH]c3c(c2c1)CCCC3"
                )
            ),
            None,
        ),
        "interpretation": {
            "external_mapping_is_observation_candidate": True,
            "external_mapping_is_ground_truth": False,
            "single_provider_confidence_is_admission_authority": False,
            "newly_promoted_record_count": 0,
        },
        "timing_seconds": {
            "mapping": round(mapping_elapsed, 3),
            "total": round(time.perf_counter() - started, 3),
        },
    }
    return report, tuple(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=DEFAULT_OUTPUT_DIRECTORY,
    )
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--stability-sample-size", type=int, default=24)
    parser.add_argument("--max-unique-reactions", type=int)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be positive")
    if args.stability_sample_size < 0:
        raise SystemExit("--stability-sample-size cannot be negative")
    if args.max_unique_reactions is not None and args.max_unique_reactions < 1:
        raise SystemExit("--max-unique-reactions must be positive")
    report, records = evaluate_fischer_rxnmapper_poc(
        args.source,
        batch_size=args.batch_size,
        stability_sample_size=args.stability_sample_size,
        max_unique_reactions=args.max_unique_reactions,
        progress=args.progress,
    )
    args.output_directory.mkdir(parents=True, exist_ok=True)
    report_path = args.output_directory / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    records_path = args.output_directory / "mapped_reactions.jsonl.gz"
    with gzip.open(records_path, "wt", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
