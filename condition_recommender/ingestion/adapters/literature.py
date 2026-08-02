"""Adapter for the literature reaction-dataset CSV contract."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator, List, Tuple

from ..models import (
    CanonicalSourceObservation,
    ConditionComponentClaim,
    ConditionInput,
    ConditionStageInput,
    OutcomeInput,
    ReactionEvidenceInput,
    SourceIdentifier,
)
from .base import (
    clean_text,
    observation_id,
    optional_float,
    optional_int,
    raw_fields,
    source_provenance,
    split_identifiers,
    supplied_mapping_status,
    validate_headers,
)


class LiteratureCsvAdapter:
    """Normalize the established literature CSV without interpreting chemistry."""

    adapter_id = "literature_csv.v1"
    adapter_version = "1.0"
    corpus_id = "literature_reaction_dataset"
    required_columns = (
        "reaction_id",
        "reaction_type",
        "yield_pct",
        "reaction_smiles",
        "reference",
        "reagent_cas",
        "catalyst_cas",
        "solvent_cas",
    )

    @staticmethod
    def _components(row: dict[str, str]) -> Tuple[ConditionComponentClaim, ...]:
        components: List[ConditionComponentClaim] = []
        for field_name, role in (
            ("catalyst_cas", "catalyst"),
            ("reagent_cas", "reagent"),
            ("solvent_cas", "solvent"),
        ):
            for position, value in enumerate(
                split_identifiers(row.get(field_name)), start=1
            ):
                components.append(
                    ConditionComponentClaim(
                        component_key=f"{field_name}_{position}",
                        source_slot=field_name,
                        source_role_hint=role,
                        identifiers=(SourceIdentifier("cas", value, field_name),),
                        provenance={"source_position": position},
                    )
                )
        return tuple(components)

    def iter_observations(
        self, path: Path, *, source_sha256: str
    ) -> Iterator[CanonicalSourceObservation]:
        """Stream one normalized observation for every literature row."""
        validate_headers(path, self.required_columns)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row_number, row in enumerate(csv.DictReader(handle), start=2):
                warnings: List[str] = []
                reaction_smiles = clean_text(row.get("reaction_smiles"))
                record_id = clean_text(row.get("reaction_id")) or (
                    f"{path.stem}:row-{row_number}"
                )
                temperature, warning = optional_float(row.get("temperature_c"))
                if warning:
                    warnings.append(f"{warning}:temperature_c")
                time_h, warning = optional_float(row.get("time_h"))
                if warning:
                    warnings.append(f"{warning}:time_h")
                declared_stages, warning = optional_int(row.get("stages"))
                if warning:
                    warnings.append(f"{warning}:stages")
                components = self._components(row)
                stage = ConditionStageInput(
                    stage_index=1,
                    component_keys=tuple(
                        component.component_key for component in components
                    ),
                    temperature_c=temperature,
                    time_h=time_h,
                    provenance={"source": "top_level_columns"},
                )
                if declared_stages is not None and declared_stages > 1:
                    warnings.append("UNSTRUCTURED_MULTISTAGE_CONDITIONS")
                outcome_value, warning = optional_float(row.get("yield_pct"))
                if warning:
                    warnings.append(f"{warning}:yield_pct")
                outcomes: List[OutcomeInput] = []
                raw_yield = clean_text(row.get("yield_pct"))
                if raw_yield:
                    outcomes.append(
                        OutcomeInput(
                            outcome_type="reported_yield_pct",
                            value=outcome_value,
                            unit="percent",
                            raw_value=raw_yield,
                            source_field="yield_pct",
                            metadata={"measurement_basis": "source_unspecified"},
                        )
                    )
                for product_index in range(1, 8):
                    field_name = f"product_yield_{product_index}"
                    raw_value = clean_text(row.get(field_name))
                    if not raw_value:
                        continue
                    value, warning = optional_float(raw_value)
                    if warning:
                        warnings.append(f"{warning}:{field_name}")
                    outcomes.append(
                        OutcomeInput(
                            outcome_type="reported_product_yield_pct",
                            value=value,
                            unit="percent",
                            raw_value=raw_value,
                            source_field=field_name,
                            metadata={"product_position": product_index},
                        )
                    )
                if not reaction_smiles:
                    warnings.append("MISSING_REACTION_SMILES")
                yield CanonicalSourceObservation(
                    observation_id=observation_id(
                        adapter_id=self.adapter_id,
                        source_sha256=source_sha256,
                        row_number=row_number,
                        record_id=record_id,
                    ),
                    observation_kind="structure_backed",
                    source=source_provenance(
                        adapter=self,
                        path=path,
                        source_sha256=source_sha256,
                        row_number=row_number,
                        record_id=record_id,
                        source_groups={
                            "source_declared_family": clean_text(
                                row.get("reaction_type")
                            ),
                            "stage_count": clean_text(row.get("stages")),
                            "step_count": clean_text(row.get("steps")),
                        },
                        reference=clean_text(row.get("reference")),
                    ),
                    reaction=ReactionEvidenceInput(
                        evidence_kind="source_structure",
                        reaction_smiles=reaction_smiles or None,
                        supplied_mapping_status=supplied_mapping_status(
                            reaction_smiles
                        ),
                        source_reaction_type=clean_text(row.get("reaction_type")),
                        structure_available=bool(reaction_smiles),
                    ),
                    conditions=ConditionInput(
                        components=components,
                        stages=(stage,),
                        declared_stage_count=declared_stages,
                        procedure_text=clean_text(row.get("experimental_procedure")),
                        warnings=tuple(
                            value
                            for value in warnings
                            if "CONDITION" in value or "stage" in value.lower()
                        ),
                    ),
                    outcomes=tuple(outcomes),
                    ingestion_status="accepted" if reaction_smiles else "review",
                    warnings=tuple(sorted(set(warnings))),
                    raw_fields=raw_fields(row),
                )


__all__ = ["LiteratureCsvAdapter"]
