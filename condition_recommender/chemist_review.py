"""Blind, chemistry-readable review packets for generic recommendations."""

from __future__ import annotations

import csv
import hashlib
import html
import json
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

from rdkit import Chem
from rdkit.Chem import Draw

from .compatibility import assess_recipe_compatibility
from .evaluation import grouped_holdout_split
from .generic_api import recommend_indexed_signature
from .generic_indexing import (
    GenericIndexedReaction,
    build_generic_index_from_rows,
    load_generic_index,
)

CHEMIST_REVIEW_PACKET_SCHEMA_VERSION = "1.1"


def _digest(*values: str) -> str:
    return hashlib.sha256("\0".join(values).encode("utf-8")).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _reactant_components(reaction_smiles: str) -> tuple[str, ...]:
    left = (
        reaction_smiles.split(">>", maxsplit=1)[0]
        if ">>" in reaction_smiles
        else reaction_smiles.split(">", maxsplit=1)[0]
    )
    return tuple(value for value in left.split(".") if value)


def _reaction_center_atoms(
    signature: Mapping[str, Any],
) -> Dict[int, list[int]]:
    values: Dict[int, set[int]] = {}
    for edit in signature.get("edits") or ():
        if not isinstance(edit, Mapping):
            continue
        for key in ("atom_1", "atom_2"):
            atom = edit.get(key)
            if not isinstance(atom, Mapping) or atom.get("side") != "reactant":
                continue
            if atom.get("component_index") is None or atom.get("atom_index") is None:
                continue
            values.setdefault(int(atom["component_index"]), set()).add(
                int(atom["atom_index"])
            )
    return {key: sorted(atoms) for key, atoms in values.items()}


def _write_structure_svg(row: GenericIndexedReaction, path: Path) -> bool:
    components = _reactant_components(row.reaction_smiles)
    molecules = []
    highlights = []
    legends = []
    center_atoms = _reaction_center_atoms(row.signature)
    for index, smiles in enumerate(components):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            continue
        molecules.append(molecule)
        highlights.append(center_atoms.get(index, []))
        legends.append(f"Reactant {index + 1}")
    if not molecules:
        return False
    svg = Draw.MolsToGridImage(
        molecules,
        molsPerRow=2,
        subImgSize=(360, 260),
        legends=legends,
        highlightAtomLists=highlights,
        useSVG=True,
    )
    path.write_text(str(svg), encoding="utf-8")
    return True


def _negative_control(
    query: GenericIndexedReaction,
    train_rows: Sequence[GenericIndexedReaction],
    excluded_core_ids: set[str],
) -> tuple[Optional[GenericIndexedReaction], str]:
    ordered = sorted(
        train_rows,
        key=lambda row: (
            _digest(query.observation_id, row.observation_id),
            row.observation_id,
        ),
    )
    structural_mismatch = None
    query_bond_key = str(query.signature.get("bond_edit_signature_key") or "")
    for row in ordered:
        if row.recipe_core_id in excluded_core_ids:
            continue
        assessment = assess_recipe_compatibility(
            query.signature,
            row.resolved_recipe,
        )
        if not assessment.compatible:
            return row, "hard_recipe_compatibility_conflict"
        if (
            structural_mismatch is None
            and str(row.signature.get("bond_edit_signature_key") or "")
            != query_bond_key
        ):
            structural_mismatch = row
    if structural_mismatch is not None:
        return structural_mismatch, "incompatible_bond_edit_precedent"
    return None, ""


def _html_list(values: Sequence[Any], *, empty: str = "None reported") -> str:
    if not values:
        return f"<p class=\"muted\">{html.escape(empty)}</p>"
    return "<ul>" + "".join(
        f"<li>{html.escape(str(value))}</li>" for value in values
    ) + "</ul>"


def _review_packet_html(packet_rows: Sequence[Mapping[str, Any]]) -> str:
    cases = []
    for case in packet_rows:
        candidates = []
        for candidate in case.get("candidates") or ():
            recipe = json.dumps(
                candidate.get("resolved_recipe") or {},
                ensure_ascii=False,
                indent=2,
                sort_keys=True,
            )
            candidates.append(
                "<article class=\"candidate\">"
                f"<h3>{html.escape(str(candidate['candidate_id']))}</h3>"
                "<h4>Canonical condition recipe</h4>"
                f"<pre>{html.escape(recipe)}</pre>"
                "<h4>Matching chemistry</h4>"
                f"{_html_list(candidate.get('matching_chemistry') or ())}"
                "<h4>Cautions</h4>"
                f"{_html_list(candidate.get('cautions') or ())}"
                "</article>"
            )
        transformation = json.dumps(
            case.get("observed_transformation") or {},
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        structure = case.get("structure_svg")
        image = (
            f'<img src="{html.escape(str(structure))}" '
            'alt="Highlighted reactant structures">'
            if structure
            else '<p class="muted">Structure rendering unavailable.</p>'
        )
        cases.append(
            "<section class=\"case\">"
            f"<h2>{html.escape(str(case['case_id']))}</h2>"
            f"{image}"
            "<h3>Reaction SMILES</h3>"
            f"<code>{html.escape(str(case['query_reaction_smiles']))}</code>"
            "<h3>Observed structural transformation</h3>"
            f"<pre>{html.escape(transformation)}</pre>"
            f"<div class=\"candidate-grid\">{''.join(candidates)}</div>"
            "</section>"
        )
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Blind generic-condition chemist review</title>
<style>
body { font-family: system-ui, sans-serif; margin: 2rem auto; max-width: 1200px;
       color: #17212b; line-height: 1.45; }
.notice { background: #eef6ff; border-left: 5px solid #356da8; padding: 1rem; }
.case { border-top: 3px solid #27384a; margin-top: 2.5rem; padding-top: 1rem; }
.candidate-grid { display: grid; gap: 1rem;
                  grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); }
.candidate { border: 1px solid #b7c2cc; border-radius: 6px; padding: 1rem; }
img { max-width: 100%; background: white; }
pre, code { background: #f4f6f8; overflow-wrap: anywhere; white-space: pre-wrap; }
pre { padding: .75rem; }
.muted { color: #65717d; }
@media print { .case { break-before: page; } }
</style>
</head>
<body>
<h1>Blind generic-condition chemist review</h1>
<div class="notice">
Assess every candidate independently and record its ID in <b>review_form.csv</b>.
Do not open <b>answer_key.jsonl</b> until the form is complete. Allowed decisions:
compatible, compatible_with_caution, incompatible, or uncertain.
</div>
""" + "".join(cases) + """
</body>
</html>
"""


def generate_chemist_review_packet(
    index_path: str | Path,
    output_dir: str | Path,
    *,
    max_cases: int = 300,
    seed: int = 17,
    top_k: int = 3,
    minimum_pool_size: int = 1,
) -> Dict[str, Any]:
    """Generate randomized review candidates and a separate answer key."""
    if max_cases < 1 or top_k < 1:
        raise ValueError("review case and candidate limits must be positive")
    index = load_generic_index(index_path)
    split = grouped_holdout_split(index.rows, test_fraction=0.2, seed=seed)
    train_index = build_generic_index_from_rows(split.train_rows)
    selected_queries = sorted(
        split.test_rows,
        key=lambda row: (
            _digest(str(seed), row.canonical_reaction_id, row.observation_id),
            row.observation_id,
        ),
    )[:max_cases]
    destination = Path(output_dir)
    structures = destination / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    packet_rows = []
    answer_rows = []
    form_rows = []
    control_count = 0
    for case_number, row in enumerate(selected_queries, start=1):
        case_id = f"REVIEW-{case_number:04d}"
        structure_path = structures / f"{case_id}.svg"
        has_structure = _write_structure_svg(row, structure_path)
        result = recommend_indexed_signature(
            row.signature,
            train_index,
            query_reaction_smiles=row.reaction_smiles,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
        )
        candidates = []
        answers = {}
        recommended_cores = {
            recommendation.recipe_core_id for recommendation in result.recommendations
        }
        for recommendation in result.recommendations:
            candidate_id = "CAND1:" + _digest(
                case_id,
                recommendation.recipe_core_id,
            )[:20]
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "recipe_core_id": recommendation.recipe_core_id,
                    "recipe_variant_ids": recommendation.recipe_variant_ids,
                    "resolved_recipe": recommendation.resolved_recipe,
                    "matching_chemistry": recommendation.explanation,
                    "cautions": recommendation.cautions,
                    "precedent_reaction_ids": (
                        recommendation.precedent_reaction_ids
                    ),
                    "precedent_reference_ids": (
                        recommendation.precedent_reference_ids
                    ),
                }
            )
            answers[candidate_id] = {
                "source_kind": "recommendation",
                "model_rank": recommendation.rank,
                "reason": recommendation.retrieval_level,
            }
        negative, negative_reason = _negative_control(
            row,
            split.train_rows,
            recommended_cores,
        )
        if negative is not None:
            candidate_id = "CAND1:" + _digest(
                case_id,
                negative.recipe_core_id,
                "negative",
            )[:20]
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "recipe_core_id": negative.recipe_core_id,
                    "recipe_variant_ids": (negative.recipe_id,),
                    "resolved_recipe": negative.resolved_recipe,
                    "matching_chemistry": (),
                    "cautions": (
                        "Assess compatibility with the query transformation",
                    ),
                    "precedent_reaction_ids": (negative.reaction_id,),
                    "precedent_reference_ids": (
                        (negative.reference_id,) if negative.reference_id else ()
                    ),
                }
            )
            answers[candidate_id] = {
                "source_kind": "negative_control",
                "model_rank": None,
                "reason": negative_reason,
            }
            control_count += 1
        candidates.sort(
            key=lambda candidate: (
                _digest(str(seed), case_id, candidate["candidate_id"]),
                candidate["candidate_id"],
            )
        )
        packet_rows.append(
            {
                "case_id": case_id,
                "query_reaction_smiles": row.reaction_smiles,
                "structure_svg": (
                    structure_path.relative_to(destination).as_posix()
                    if has_structure
                    else None
                ),
                "observed_transformation": {
                    "transformation_class": row.transformation_class,
                    "named_family": row.named_family or None,
                    "formed_bond_types": row.signature.get("formed_bond_types") or (),
                    "broken_bond_types": row.signature.get("broken_bond_types") or (),
                    "order_changes": row.signature.get("order_changes") or (),
                    "hydrogen_changes": row.signature.get("hydrogen_changes") or (),
                    "reaction_scope": (
                        (row.signature.get("topology") or {}).get(
                            "reaction_scope"
                        )
                    ),
                },
                "retrieval_level": result.retrieval_level,
                "warnings": result.warnings,
                "candidates": candidates,
                "reviewer_decision": "",
                "reviewer_comments": "",
            }
        )
        for candidate in candidates:
            answer_rows.append(
                {
                    "case_id": case_id,
                    "candidate_id": candidate["candidate_id"],
                    **answers[candidate["candidate_id"]],
                }
            )
            form_rows.append(
                {
                    "case_id": case_id,
                    "candidate_id": candidate["candidate_id"],
                    "structure_svg": (
                        structure_path.relative_to(destination).as_posix()
                        if has_structure
                        else ""
                    ),
                    "chemist_decision": "",
                    "chemist_confidence": "",
                    "chemist_comments": "",
                }
            )
    (destination / "review_packet.jsonl").write_text(
        "".join(
            json.dumps(value, ensure_ascii=False, sort_keys=True) + "\n"
            for value in packet_rows
        ),
        encoding="utf-8",
    )
    (destination / "review_packet.html").write_text(
        _review_packet_html(packet_rows),
        encoding="utf-8",
    )
    (destination / "answer_key.jsonl").write_text(
        "".join(
            json.dumps(value, ensure_ascii=False, sort_keys=True) + "\n"
            for value in answer_rows
        ),
        encoding="utf-8",
    )
    with (destination / "review_form.csv").open(
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        fields = (
            "case_id",
            "candidate_id",
            "structure_svg",
            "chemist_decision",
            "chemist_confidence",
            "chemist_comments",
        )
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(form_rows)
    report = {
        "schema_version": CHEMIST_REVIEW_PACKET_SCHEMA_VERSION,
        "artifact_type": "generic_chemist_review_packet",
        "index_path": str(Path(index_path)),
        "index_sha256": _file_sha256(Path(index_path)),
        "seed": seed,
        "max_cases": max_cases,
        "top_k": top_k,
        "minimum_pool_size": minimum_pool_size,
        "case_count": len(packet_rows),
        "candidate_count": len(answer_rows),
        "negative_control_count": control_count,
        "split": {
            "train_row_count": len(split.train_rows),
            "test_row_count": len(split.test_rows),
            "leakage_group_count": len(
                set(split.train_group_ids) & set(split.test_group_ids)
            ),
        },
        "files": {
            "packet": "review_packet.jsonl",
            "html_packet": "review_packet.html",
            "answer_key": "answer_key.jsonl",
            "review_form": "review_form.csv",
            "structures": "structures/",
        },
    }
    (destination / "review_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


__all__ = [
    "CHEMIST_REVIEW_PACKET_SCHEMA_VERSION",
    "generate_chemist_review_packet",
]
