"""Convert literature-style Suzuki rows into versioned recommendation records."""

from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

from reactive_taxonomy import featurize_reaction

from ..condition_normalization import normalize_conditions, optional_float
from ..models import AdmissionTier, RecommendationRecord


def _admission(
    *, valid: bool, evidence: str, family: str | None,
    yield_pct: float | None, conditions_complete: bool,
    has_environment: bool, warnings: Iterable[str],
) -> Tuple[AdmissionTier, Tuple[str, ...]]:
    reasons: List[str] = []
    warning_set = set(warnings)
    if not valid:
        return AdmissionTier.REJECTED, ("invalid_reaction_smiles",)
    if yield_pct is None or not 0.0 <= yield_pct <= 100.0:
        return AdmissionTier.REJECTED, ("missing_or_invalid_yield",)
    if evidence == "unresolved":
        return AdmissionTier.REJECTED, ("unresolved_transformation",)
    if evidence == "ambiguous":
        reasons.append("ambiguous_participating_sites")
    if evidence != "exact_product_reconstruction":
        reasons.append("reaction_not_exactly_verified")
    if family != "suzuki_miyaura":
        reasons.append("not_verified_as_suzuki_miyaura")
    if not has_environment:
        reasons.append("missing_suzuki_family_environment")
    if not conditions_complete:
        reasons.append("incomplete_condition_identity")
    if "MULTIPLE_PRODUCTS" in warning_set:
        reasons.append("multiple_products")
    if reasons:
        return AdmissionTier.REVIEW, tuple(sorted(set(reasons)))
    return AdmissionTier.VERIFIED, ("exact_single_event_suzuki_with_conditions",)


def _source_family_label(analysis: Any, source_reaction_type: str) -> Tuple[str | None, str]:
    """Use declared dataset family only to narrow partial review labels."""
    if analysis.selected_candidate is not None:
        return analysis.reaction_label, analysis.reaction_label_status
    normalized_type = source_reaction_type.strip().lower().replace("-", "_")
    if normalized_type not in {"suzuki", "suzuki_miyaura"}:
        return analysis.reaction_label, analysis.reaction_label_status
    labels = sorted({
        str(candidate.reaction_label).split(" → ", 1)[0]
        for candidate in analysis.candidates
        if "suzuki_miyaura" in candidate.compatible_named_families and candidate.reaction_label
    })
    if len(labels) == 1:
        return f"{labels[0]} →", "source_family_reactant_only"
    if labels:
        return " OR ".join(f"({label})" for label in labels) + " →", "source_family_ambiguous_reactants"
    return analysis.reaction_label, analysis.reaction_label_status


def convert_row(row: Dict[str, Any], source_row_number: int) -> RecommendationRecord:
    reaction_smiles = str(row.get("reaction_smiles") or "").strip()
    analysis = featurize_reaction(reaction_smiles)
    yield_pct = optional_float(row.get("yield_pct"))
    conditions = normalize_conditions(row)
    family_environment = asdict(analysis.family_environment) if analysis.family_environment else None
    reaction_label, reaction_label_status = _source_family_label(
        analysis, str(row.get("reaction_type") or "")
    )
    tier, reasons = _admission(
        valid=analysis.valid,
        evidence=analysis.evidence_quality,
        family=analysis.named_family,
        yield_pct=yield_pct,
        conditions_complete=conditions.complete,
        has_environment=family_environment is not None,
        warnings=analysis.warnings,
    )
    return RecommendationRecord(
        reaction_id=str(row.get("reaction_id") or f"row-{source_row_number}"),
        source_row_number=source_row_number,
        reaction_smiles=reaction_smiles,
        admission_tier=tier,
        admission_reasons=reasons,
        evidence_quality=analysis.evidence_quality,
        named_family=analysis.named_family,
        reaction_label=reaction_label,
        reaction_label_status=reaction_label_status,
        yield_pct=yield_pct,
        temperature_c=optional_float(row.get("temperature_c")),
        time_h=optional_float(row.get("time_h")),
        conditions=conditions,
        family_environment=family_environment,
        spectator_groups=tuple(asdict(group) for group in analysis.spectator_groups),
        source={
            "reaction_type": str(row.get("reaction_type") or ""),
            "reference": str(row.get("reference") or ""),
            "reactant_cas": str(row.get("reactant_cas") or ""),
            "product_cas": str(row.get("product_cas") or ""),
            "notes": str(row.get("notes") or ""),
            "stages": str(row.get("stages") or ""),
            "steps": str(row.get("steps") or ""),
        },
    )


def _partner(record: RecommendationRecord, role: str) -> Dict[str, Any]:
    environment = record.family_environment or {}
    return next((item for item in environment.get("partners", []) if item.get("role") == role), {})


def flatten_record(record: RecommendationRecord) -> Dict[str, Any]:
    electrophile = _partner(record, "electrophile")
    transfer = _partner(record, "transfer_partner")
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
        "named_family": record.named_family or "",
        "reaction_label": record.reaction_label or "",
        "reaction_label_status": record.reaction_label_status,
        "yield_pct": record.yield_pct if record.yield_pct is not None else "",
        "temperature_c": record.temperature_c if record.temperature_c is not None else "",
        "time_h": record.time_h if record.time_h is not None else "",
        "condition_recipe_id": record.conditions.recipe_id,
        "conditions_complete": record.conditions.complete,
        "catalyst_cas": joined(record.conditions.catalyst_cas),
        "reagent_cas": joined(record.conditions.reagent_cas),
        "solvent_cas": joined(record.conditions.solvent_cas),
        "electrophile_label": electrophile.get("chemist_label", ""),
        "electrophile_context": electrophile.get("anchor_context", ""),
        "electrophile_handle": electrophile.get("handle_token", ""),
        "electrophile_steric_class": (electrophile.get("steric") or {}).get("class", ""),
        "electrophile_electronic_class": (electrophile.get("electronic") or {}).get("class", ""),
        "electrophile_flags": joined(electrophile.get("flags") or []),
        "transfer_label": transfer.get("chemist_label", ""),
        "transfer_context": transfer.get("anchor_context", ""),
        "transfer_handle": transfer.get("handle_token", ""),
        "transfer_steric_class": (transfer.get("steric") or {}).get("class", ""),
        "transfer_electronic_class": (transfer.get("electronic") or {}).get("class", ""),
        "transfer_flags": joined(transfer.get("flags") or []),
        "spectator_group_ids": joined(sorted({item["group_id"] for item in record.spectator_groups})),
        "family_environment_json": json.dumps(record.family_environment, ensure_ascii=False, sort_keys=True) if record.family_environment else "",
        "spectator_groups_json": json.dumps(record.spectator_groups, ensure_ascii=False, sort_keys=True),
        "reference": record.source.get("reference", ""),
        "notes": record.source.get("notes", ""),
    }


def _write_csv(path: Path, records: List[RecommendationRecord]) -> None:
    rows = [flatten_record(record) for record in records]
    fieldnames = list(flatten_record(records[0]).keys()) if records else list(flatten_record(RecommendationRecord(
        reaction_id="", source_row_number=0, reaction_smiles="", admission_tier=AdmissionTier.REJECTED,
        admission_reasons=(), evidence_quality="", named_family=None, reaction_label=None,
        reaction_label_status="unavailable",
        yield_pct=None, temperature_c=None, time_h=None,
        conditions=normalize_conditions({}),
    )).keys())
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def convert_file(input_path: str | Path, output_dir: str | Path) -> Dict[str, Any]:
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        source_rows = list(csv.DictReader(handle))
    records = [convert_row(row, index) for index, row in enumerate(source_rows, start=2)]
    by_tier = {
        tier: [record for record in records if record.admission_tier == tier]
        for tier in AdmissionTier
    }
    for tier, tier_records in by_tier.items():
        _write_csv(output_dir / f"{tier.value}.csv", tier_records)
    reason_counts = Counter(reason for record in records for reason in record.admission_reasons)
    report = {
        "schema_version": "1.0",
        "input_file": str(input_path),
        "total_rows": len(records),
        "tier_counts": {tier.value: len(by_tier[tier]) for tier in AdmissionTier},
        "reason_counts": dict(sorted(reason_counts.items())),
        "verified_yield_mean": round(sum(record.yield_pct or 0.0 for record in by_tier[AdmissionTier.VERIFIED]) / len(by_tier[AdmissionTier.VERIFIED]), 3) if by_tier[AdmissionTier.VERIFIED] else None,
        "temperature_coverage": sum(record.temperature_c is not None for record in records),
        "time_coverage": sum(record.time_h is not None for record in records),
        "reaction_label_coverage": {
            tier.value: sum(bool(record.reaction_label) for record in by_tier[tier])
            for tier in AdmissionTier
        },
        "reaction_label_status_counts": dict(sorted(Counter(
            record.reaction_label_status for record in records
        ).items())),
    }
    with (output_dir / "conversion_report.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    markdown = f"""# Suzuki recommendation pilot conversion

- Source rows: {report['total_rows']}
- Verified: {report['tier_counts']['verified']}
- Review: {report['tier_counts']['review']}
- Rejected: {report['tier_counts']['rejected']}
- Mean verified yield: {report['verified_yield_mean']}%
- Temperature coverage: {report['temperature_coverage']}/{report['total_rows']}
- Time coverage: {report['time_coverage']}/{report['total_rows']}
- Review-label coverage: {report['reaction_label_coverage']['review']}/{report['tier_counts']['review']}

## Admission policy

Verified records require exact product reconstruction, a unique Suzuki–Miyaura
family environment, a yield from 0–100%, and non-empty catalyst, reagent, and
solvent CAS identity groups. Review records retain usable but incomplete or
unverified observations. Invalid, yield-invalid, and unresolved reactions are
rejected.

## Reason counts

""" + "\n".join(
        f"- `{reason}`: {count}" for reason, count in report["reason_counts"].items()
    ) + "\n"
    (output_dir / "conversion_report.md").write_text(markdown, encoding="utf-8")
    return report


__all__ = ["convert_file", "convert_row", "flatten_record"]
