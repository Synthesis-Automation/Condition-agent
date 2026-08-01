"""Shared conversion for sp² C–O and C–S coupling datasets."""

from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List

from reactive_taxonomy import featurize_reaction

from ..condition_normalization import normalize_conditions, optional_float
from ..models import AdmissionTier, RecommendationRecord
from .profile_fields import profile_classes
from .signature_serialization import flattened_signature_fields, signature_record_fields

_CONFIG = {
    "O": {
        "interpretation": "sp2_c_o_substitution",
        "family": "c_o_coupling",
        "connection": "C_O",
    },
    "S": {
        "interpretation": "sp2_c_s_substitution",
        "family": "c_s_coupling",
        "connection": "C_S",
    },
}


def convert_row(
    row: Dict[str, Any], source_row_number: int, *, element: str
) -> RecommendationRecord:
    config = _CONFIG[element.upper()]
    analysis = featurize_reaction(str(row.get("reaction_smiles") or ""))
    conditions, yield_pct = (
        normalize_conditions(row),
        optional_float(row.get("yield_pct")),
    )
    reasons: List[str] = []
    if not analysis.valid:
        tier, reasons = AdmissionTier.REJECTED, ["invalid_reaction_smiles"]
    elif yield_pct is None or not 0 <= yield_pct <= 100:
        tier, reasons = AdmissionTier.REJECTED, ["missing_or_invalid_yield"]
    elif analysis.evidence_quality == "unresolved":
        tier, reasons = AdmissionTier.REJECTED, ["unresolved_transformation"]
    else:
        selected = analysis.selected_candidate
        if analysis.evidence_quality != "exact_product_reconstruction":
            reasons.append("reaction_not_exactly_verified")
        if (
            selected is None
            or selected.annotation_id != config["interpretation"]
        ):
            reasons.append(
                f"not_verified_as_{config['interpretation']}"
            )
        if (
            analysis.family_environment is None
            or analysis.family_environment.family_id != config["family"]
        ):
            reasons.append(f"missing_{config['family']}_environment")
        if (
            analysis.product_connection is None
            or analysis.product_connection.connection_type != config["connection"]
        ):
            reasons.append(f"missing_{config['connection'].lower()}_product_connection")
        if not conditions.complete:
            reasons.append("incomplete_condition_identity")
        if "MULTIPLE_PRODUCTS" in analysis.warnings:
            reasons.append("multiple_products")
        tier = AdmissionTier.REVIEW if reasons else AdmissionTier.VERIFIED
        if not reasons:
            reasons = [
                f"exact_single_event_{config['connection'].lower()}_with_conditions"
            ]
    return RecommendationRecord(
        reaction_id=str(row.get("reaction_id") or f"row-{source_row_number}"),
        source_row_number=source_row_number,
        reaction_smiles=str(row.get("reaction_smiles") or ""),
        admission_tier=tier,
        admission_reasons=tuple(sorted(set(reasons))),
        evidence_quality=analysis.evidence_quality,
        named_family=analysis.named_family,
        reaction_label=(
            asdict(analysis.reaction_label)
            if analysis.reaction_label is not None
            else None
        ),
        yield_pct=yield_pct,
        temperature_c=optional_float(row.get("temperature_c")),
        time_h=optional_float(row.get("time_h")),
        conditions=conditions,
        family_environment=asdict(analysis.family_environment)
        if analysis.family_environment
        else None,
        product_connection=asdict(analysis.product_connection)
        if analysis.product_connection
        else None,
        spectator_groups=tuple(asdict(group) for group in analysis.spectator_groups),
        **signature_record_fields(analysis),
        source={
            "reaction_type": row.get("reaction_type", ""),
            "reference": row.get("reference", ""),
            "notes": row.get("notes", ""),
        },
    )


def flatten_record(record: RecommendationRecord) -> Dict[str, Any]:
    environment = record.family_environment or {}
    partners = {item.get("role"): item for item in environment.get("partners") or []}
    electrophile, nucleophile = (
        partners.get("electrophile", {}),
        partners.get("nucleophile", {}),
    )
    features = nucleophile.get("features") or {}
    electrophile_steric, electrophile_electronic = profile_classes(electrophile)
    def joined(values: Iterable[Any]) -> str:
        return "|".join(str(value) for value in values)
    return {
        "schema_version": record.schema_version,
        "reaction_id": record.reaction_id,
        "source_row_number": record.source_row_number,
        "admission_tier": record.admission_tier.value,
        "admission_reasons": joined(record.admission_reasons),
        "reaction_smiles": record.reaction_smiles,
        "evidence_quality": record.evidence_quality,
        "reaction_display_label": (record.reaction_label or {}).get("concise", ""),
        "reaction_display_label_detailed": (record.reaction_label or {}).get(
            "detailed", ""
        ),
        "reaction_display_status": (record.reaction_label or {}).get("status", ""),
        "yield_pct": record.yield_pct if record.yield_pct is not None else "",
        "condition_recipe_id": record.conditions.recipe_id,
        "conditions_complete": record.conditions.complete,
        "catalyst_cas": joined(record.conditions.catalyst_cas),
        "reagent_cas": joined(record.conditions.reagent_cas),
        "solvent_cas": joined(record.conditions.solvent_cas),
        "electrophile_label": electrophile.get("chemist_label", ""),
        "electrophile_context": electrophile.get("anchor_context", ""),
        "electrophile_handle": electrophile.get("handle_token", ""),
        "electrophile_steric_class": electrophile_steric,
        "electrophile_electronic_class": electrophile_electronic,
        "electrophile_flags": joined(electrophile.get("flags") or []),
        "nucleophile_label": nucleophile.get("chemist_label", ""),
        "nucleophile_family": features.get("derived_family", ""),
        "nucleophile_contexts": joined(features.get("retained_contexts") or []),
        "nucleophile_availability": features.get("availability", ""),
        "nucleophile_flags": joined(nucleophile.get("flags") or []),
        "product_connection_label": (record.product_connection or {}).get(
            "concise_label", ""
        ),
        "product_connection_json": json.dumps(
            record.product_connection, ensure_ascii=False, sort_keys=True
        )
        if record.product_connection
        else "",
        "spectator_group_ids": joined(
            sorted({item["group_id"] for item in record.spectator_groups})
        ),
        "family_environment_json": json.dumps(
            record.family_environment, ensure_ascii=False, sort_keys=True
        )
        if record.family_environment
        else "",
        "transformation_class": record.transformation_class or "",
        "transformation_confidence": record.transformation_confidence,
        "family_confidence": record.family_confidence,
        **flattened_signature_fields(record.reaction_signature),
        "reference": record.source.get("reference", ""),
        "notes": record.source.get("notes", ""),
    }


def convert_file(
    input_path: str | Path, output_dir: str | Path, *, element: str
) -> Dict[str, Any]:
    element = element.upper()
    config = _CONFIG[element]
    input_path, output_dir = Path(input_path), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    records = [
        convert_row(row, index, element=element)
        for index, row in enumerate(rows, start=2)
    ]
    by_tier = {
        tier: [record for record in records if record.admission_tier == tier]
        for tier in AdmissionTier
    }
    for tier, values in by_tier.items():
        flat = [flatten_record(record) for record in values]
        fields = (
            list(flat[0])
            if flat
            else list(flatten_record(convert_row({}, 0, element=element)))
        )
        with (output_dir / f"{tier.value}.csv").open(
            "w", encoding="utf-8-sig", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(flat)
    report = {
        "schema_version": "1.0",
        "element": element,
        "transformation": config["connection"],
        "total_rows": len(records),
        "tier_counts": {tier.value: len(by_tier[tier]) for tier in AdmissionTier},
        "reason_counts": dict(
            sorted(
                Counter(
                    reason for record in records for reason in record.admission_reasons
                ).items()
            )
        ),
        "label_coverage": {
            tier.value: sum(bool(record.reaction_label) for record in by_tier[tier])
            for tier in AdmissionTier
        },
        "product_connection_coverage": {
            tier.value: sum(bool(record.product_connection) for record in by_tier[tier])
            for tier in AdmissionTier
        },
    }
    (output_dir / "conversion_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return report


__all__ = ["convert_file", "convert_row", "flatten_record"]
