"""Build and persist structurally diverse generic template libraries."""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Literal

from .chemistry import digest

try:
    import orjson
except ImportError:  # pragma: no cover - compatibility outside the web runtime
    orjson = None

from .generic_compiler import compile_generic_templates
from .generic_models import (
    GenericCoreTemplate,
    GenericGraphOperator,
    GenericHandleCompletionGroup,
    GenericTemplateLibrary,
    GenericTemplatePrecedent,
)
from .retrieval_index import build_generic_retrieval_index


GENERIC_LIBRARY_DEFINITION = {
    # Preserve the serialized v3 identity across the Python-package rename.
    "definition_id": "generic_core_retrosynthesis_poc.v3",
    "routing_source": "normalized_graph_edits",
    "named_families_used_for_routing": False,
    "source_round_trip_required": True,
    "operator_identity": "handle_independent_normalized_graph_edits",
    "supported_transformations": [
        "acyl_substitution",
        "c_c_coupling",
        "carbonyl_oxidation",
        "carbonyl_reduction",
        "conjugate_addition",
        "carbonyl_condensation",
        "ring_formation",
    ],
}


AdmissionCallback = Callable[[Dict[str, Any]], None]


def _row_reference_id(row: Dict[str, Any]) -> str:
    direct = str(row.get("reference_id") or "").strip()
    if direct:
        return direct
    return str((row.get("reference_identity") or {}).get("reference_id") or "")


def _row_provenance(row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "source_shard": str(
            row.get("_build_source_shard")
            or row.get("_sampling_shard")
            or row.get("source_path")
            or ""
        ),
        "source_row_number": row.get("_build_source_row_number")
        or row.get("source_row_number"),
    }


def _precedent_identity(precedent: GenericTemplatePrecedent) -> str:
    return precedent.reaction_id or precedent.mapped_reaction_smiles


def _context_key(precedent: GenericTemplatePrecedent) -> str:
    return json.dumps(
        asdict(precedent.context),
        sort_keys=True,
        separators=(",", ":"),
    )


def select_context_representatives(
    precedents: Iterable[GenericTemplatePrecedent],
    limit: int,
) -> tuple[GenericTemplatePrecedent, ...]:
    """Select deterministic representatives across distinct context bins."""

    unique = {
        _precedent_identity(precedent): precedent for precedent in precedents
    }
    by_context: dict[str, list[GenericTemplatePrecedent]] = {}
    for precedent in unique.values():
        by_context.setdefault(_context_key(precedent), []).append(precedent)
    representatives = []
    for context_key, values in by_context.items():
        representatives.append(
            (
                hashlib.sha256(context_key.encode()).hexdigest(),
                min(values, key=_precedent_identity),
            )
        )
    selected = [
        precedent
        for _, precedent in sorted(
            representatives,
            key=lambda item: (item[0], _precedent_identity(item[1])),
        )[:limit]
    ]
    if len(selected) < limit:
        selected_ids = {_precedent_identity(value) for value in selected}
        remaining = sorted(
            (
                value
                for value in unique.values()
                if _precedent_identity(value) not in selected_ids
            ),
            key=lambda value: hashlib.sha256(
                _precedent_identity(value).encode()
            ).hexdigest(),
        )
        selected.extend(remaining[: limit - len(selected)])
    return tuple(sorted(selected, key=_precedent_identity))


def build_generic_library(
    rows: Iterable[Dict[str, Any]],
    *,
    engine: Literal["reaction_core", "rdchiral"] = "reaction_core",
    levels: Iterable[Literal["L0", "L1", "L2"]] = ("L1", "L2"),
    admission_mode: Literal["supported", "data_driven"] = "supported",
    max_precedents_per_template: int = 8,
    admission_callback: AdmissionCallback | None = None,
) -> GenericTemplateLibrary:
    """Compile and aggregate a generic source-round-tripped library."""

    if max_precedents_per_template < 1:
        raise ValueError("max precedents must be positive")
    source_count = 0
    accepted_count = 0
    annotated_accepted_count = 0
    unannotated_accepted_count = 0
    rejections: Counter[str] = Counter()
    templates: dict[str, GenericCoreTemplate] = {}
    support_units: dict[str, set[str]] = {}
    operator_support_units: dict[str, set[str]] = {}
    operator_observations: dict[str, set[str]] = {}
    completion_support_units: dict[str, set[str]] = {}
    completion_observations: dict[str, set[str]] = {}
    for row in rows:
        source_count += 1
        result = compile_generic_templates(
            row,
            engine=engine,
            levels=levels,
            admission_mode=admission_mode,
        )
        if not result.templates:
            rejections[str(result.rejection_reason or "unknown_rejection")] += 1
            if admission_callback is not None:
                admission_callback(
                    {
                        **_row_provenance(row),
                        "status": "rejected",
                        "reaction_id": str(row.get("reaction_id") or ""),
                        "reference_id": _row_reference_id(row),
                        "reason": str(
                            result.rejection_reason or "unknown_rejection"
                        ),
                        "stage": result.rejection_stage or "unknown",
                        "diagnostics": dict(result.diagnostics or {}),
                    }
                )
            continue
        accepted_count += 1
        if result.templates[0].named_annotations:
            annotated_accepted_count += 1
        else:
            unannotated_accepted_count += 1
        first_precedent = result.templates[0].precedents[0]
        observation_key = _precedent_identity(first_precedent)
        support_key = first_precedent.reference_id or observation_key
        accepted_template_ids = []
        accepted_operator_ids = set()
        accepted_completion_ids = set()
        for compiled in result.templates:
            current = templates.get(compiled.template_id)
            support_units.setdefault(compiled.template_id, set()).add(support_key)
            accepted_template_ids.append(compiled.template_id)
            accepted_operator_ids.add(compiled.operator_id)
            completion_id = digest(
                "COMP2",
                compiled.operator_id,
                compiled.handle_signature,
            )
            accepted_completion_ids.add(completion_id)
            operator_support_units.setdefault(compiled.operator_id, set()).add(
                support_key
            )
            operator_observations.setdefault(compiled.operator_id, set()).add(
                observation_key
            )
            completion_support_units.setdefault(completion_id, set()).add(
                support_key
            )
            completion_observations.setdefault(completion_id, set()).add(
                observation_key
            )
            if current is None:
                templates[compiled.template_id] = compiled
                continue
            precedents = select_context_representatives(
                (*current.precedents, compiled.precedents[0]),
                max_precedents_per_template,
            )
            templates[compiled.template_id] = replace(
                current,
                observation_support=current.observation_support + 1,
                precedents=tuple(precedents),
            )
        if admission_callback is not None:
            admission_callback(
                {
                    **_row_provenance(row),
                    "status": "accepted",
                    "reaction_id": first_precedent.reaction_id,
                    "reference_id": first_precedent.reference_id,
                    "observation_key": observation_key,
                    "support_key": support_key,
                    "template_ids": sorted(set(accepted_template_ids)),
                    "operator_ids": sorted(accepted_operator_ids),
                    "completion_group_ids": sorted(accepted_completion_ids),
                    "named_annotations": list(
                        result.templates[0].named_annotations
                    ),
                    "diagnostics": dict(result.diagnostics or {}),
                }
            )
    finalized = tuple(
        replace(
            template,
            independent_reference_support=len(support_units[template.template_id]),
        )
        for template in sorted(templates.values(), key=lambda item: item.template_id)
    )
    by_operator: dict[str, list[GenericCoreTemplate]] = {}
    for template in finalized:
        by_operator.setdefault(template.operator_id, []).append(template)
    operators = []
    for operator_id, members in sorted(by_operator.items()):
        operators.append(
            GenericGraphOperator(
                operator_id=operator_id,
                operator_signature=members[0].operator_signature,
                edit_tokens=tuple(
                    sorted(
                        {
                            token
                            for member in members
                            for token in member.edit_tokens
                        }
                    )
                ),
                realization_ids=tuple(
                    sorted(
                        {
                            member.realization_id
                            for member in members
                            if member.realization_id
                        }
                    )
                ),
                abstraction_levels=tuple(
                    sorted({member.abstraction_level for member in members})
                ),
                observation_support=len(operator_observations[operator_id]),
                independent_reference_support=len(
                    operator_support_units[operator_id]
                ),
                named_annotations=tuple(
                    sorted(
                        {
                            annotation
                            for member in members
                            for annotation in member.named_annotations
                        }
                    )
                ),
            )
        )
    by_completion: dict[str, list[GenericCoreTemplate]] = {}
    for template in finalized:
        completion_id = digest(
            "COMP2",
            template.operator_id,
            template.handle_signature,
        )
        by_completion.setdefault(completion_id, []).append(template)
    completion_groups = tuple(
        GenericHandleCompletionGroup(
            completion_group_id=completion_id,
            operator_id=members[0].operator_id,
            completion_signature=members[0].handle_signature,
            synthon_signatures=tuple(
                sorted({member.synthon_signature for member in members})
            ),
            realization_ids=tuple(
                sorted({member.realization_id for member in members})
            ),
            template_ids=tuple(sorted(member.template_id for member in members)),
            handle_signatures=tuple(
                sorted({member.handle_signature for member in members})
            ),
            observation_support=len(completion_observations[completion_id]),
            independent_reference_support=len(
                completion_support_units[completion_id]
            ),
        )
        for completion_id, members in sorted(by_completion.items())
    )
    library = GenericTemplateLibrary(
        templates=finalized,
        source_row_count=source_count,
        accepted_observation_count=accepted_count,
        rejection_counts=dict(sorted(rejections.items())),
        definition={
            **GENERIC_LIBRARY_DEFINITION,
            "compiler_engine": engine,
            "admission_mode": admission_mode,
            "annotated_accepted_observation_count": annotated_accepted_count,
            "unannotated_accepted_observation_count": unannotated_accepted_count,
        },
        operators=tuple(operators),
        completion_groups=completion_groups,
    )
    return replace(
        library,
        retrieval_index=build_generic_retrieval_index(library),
    )


def save_generic_library(
    library: GenericTemplateLibrary,
    destination: str | Path,
) -> None:
    """Write canonical gzip JSON."""

    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        library.to_dict(),
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as handle:
            handle.write(payload)


def load_generic_library(source: str | Path) -> GenericTemplateLibrary:
    """Load a generic gzip JSON library."""

    with gzip.open(Path(source), "rb") as handle:
        payload = handle.read()
    value = (
        orjson.loads(payload)
        if orjson is not None
        else json.loads(payload.decode("utf-8"))
    )
    if not isinstance(value, dict):
        raise ValueError("generic library must contain an object")
    return GenericTemplateLibrary.from_dict(value)


__all__ = [
    "GENERIC_LIBRARY_DEFINITION",
    "build_generic_library",
    "load_generic_library",
    "save_generic_library",
]
