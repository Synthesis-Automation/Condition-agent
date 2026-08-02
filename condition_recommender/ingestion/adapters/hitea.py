"""Adapter for the HiTEA approved full-dataset CSV contract."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

from ..models import (
    CanonicalSourceObservation,
    ConditionAmountInput,
    ConditionComponentClaim,
    ConditionComponentGroup,
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
    raw_fields,
    source_provenance,
    supplied_mapping_status,
    validate_headers,
)


_ABSENT_VALUES = {"", "none", "null", "n/a", "na"}


def _present(value: object) -> bool:
    return clean_text(value).casefold() not in _ABSENT_VALUES


def _amount(
    value: object, field_name: str
) -> tuple[Optional[ConditionAmountInput], Optional[str]]:
    raw = clean_text(value)
    if not raw:
        return None, None
    number, warning = optional_float(raw)
    return ConditionAmountInput(number, "equivalent", raw), (
        f"{warning}:{field_name}" if warning else None
    )


class HiTeaCsvAdapter:
    """Preserve mapped reactions, screening provenance, and condition claims."""

    adapter_id = "hitea_approved_csv.v1"
    adapter_version = "1.0"
    corpus_id = "hitea"
    required_columns = (
        "Dataset entry number",
        "ReactionClass",
        "SCREEN_ID",
        "NOTEBOOK_ID",
        "REACTION_ID",
        "RXN_SMILES",
        "Product_Yield_PCT_Area_UV",
        "Solvent_1_Name",
        "Reaction_T",
        "Reaction_Time_hrs",
    )

    @staticmethod
    def _condition_claims(
        row: dict[str, str], warnings: List[str]
    ) -> tuple[
        Tuple[ConditionComponentClaim, ...], Tuple[ConditionComponentGroup, ...]
    ]:
        components: List[ConditionComponentClaim] = []
        groups: List[ConditionComponentGroup] = []
        for slot in (1, 2):
            group_key = f"catalyst_{slot}"
            amount, warning = _amount(
                row.get(f"Catalyst_{slot}_eq"), f"Catalyst_{slot}_eq"
            )
            if warning:
                warnings.append(warning)
            shorthand = clean_text(row.get(f"Catalyst_{slot}_Short_Hand"))
            group_has_components = False
            for position in (1, 2):
                id_field = f"Catalyst_{slot}_ID[{position}]"
                smiles_field = f"catalyst_{slot}_ID_{position}_SMILES"
                identifiers = []
                if _present(row.get(id_field)):
                    identifiers.append(
                        SourceIdentifier("mfcd", clean_text(row[id_field]), id_field)
                    )
                if _present(row.get(smiles_field)):
                    identifiers.append(
                        SourceIdentifier(
                            "smiles", clean_text(row[smiles_field]), smiles_field
                        )
                    )
                if not identifiers:
                    continue
                group_has_components = True
                components.append(
                    ConditionComponentClaim(
                        component_key=f"{group_key}_component_{position}",
                        source_slot=id_field,
                        source_role_hint="catalyst",
                        identifiers=tuple(identifiers),
                        group_key=group_key,
                        provenance={"source_position": position},
                    )
                )
            if shorthand and not group_has_components:
                components.append(
                    ConditionComponentClaim(
                        component_key=f"{group_key}_shorthand",
                        source_slot=f"Catalyst_{slot}_Short_Hand",
                        source_role_hint="catalyst",
                        identifiers=(
                            SourceIdentifier(
                                "name",
                                shorthand,
                                f"Catalyst_{slot}_Short_Hand",
                            ),
                        ),
                        amount=amount,
                        group_key=group_key,
                    )
                )
            if group_has_components or shorthand:
                groups.append(
                    ConditionComponentGroup(
                        group_key=group_key,
                        source_slot=f"Catalyst_{slot}",
                        display_name=shorthand,
                        amount=amount,
                    )
                )
        for slot in (1, 2):
            id_field = f"Reagent_{slot}_ID"
            name_field = f"Reagent_{slot}_Short_Hand"
            identifiers = []
            if _present(row.get(id_field)):
                identifiers.append(
                    SourceIdentifier("mfcd", clean_text(row[id_field]), id_field)
                )
            if _present(row.get(name_field)):
                identifiers.append(
                    SourceIdentifier("name", clean_text(row[name_field]), name_field)
                )
            if not identifiers:
                continue
            amount, warning = _amount(
                row.get(f"Reagent_{slot}_eq"), f"Reagent_{slot}_eq"
            )
            if warning:
                warnings.append(warning)
            components.append(
                ConditionComponentClaim(
                    component_key=f"reagent_{slot}",
                    source_slot=id_field,
                    source_role_hint="reagent",
                    identifiers=tuple(identifiers),
                    amount=amount,
                )
            )
        solvent = clean_text(row.get("Solvent_1_Name"))
        if solvent:
            components.append(
                ConditionComponentClaim(
                    component_key="solvent_1",
                    source_slot="Solvent_1_Name",
                    source_role_hint="solvent",
                    identifiers=(SourceIdentifier("name", solvent, "Solvent_1_Name"),),
                )
            )
        return tuple(components), tuple(groups)

    def iter_observations(
        self, path: Path, *, source_sha256: str
    ) -> Iterator[CanonicalSourceObservation]:
        """Stream normalized HiTEA screening observations."""
        validate_headers(path, self.required_columns)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row_number, row in enumerate(csv.DictReader(handle), start=2):
                warnings: List[str] = []
                record_id = (
                    clean_text(row.get("REACTION_ID"))
                    or clean_text(row.get("Dataset entry number"))
                    or f"{path.stem}:row-{row_number}"
                )
                reaction_smiles = clean_text(row.get("RXN_SMILES"))
                temperature, warning = optional_float(row.get("Reaction_T"))
                if warning:
                    warnings.append(f"{warning}:Reaction_T")
                time_h, warning = optional_float(row.get("Reaction_Time_hrs"))
                if warning:
                    warnings.append(f"{warning}:Reaction_Time_hrs")
                components, groups = self._condition_claims(row, warnings)
                stage = ConditionStageInput(
                    stage_index=1,
                    component_keys=tuple(
                        component.component_key for component in components
                    ),
                    temperature_c=temperature,
                    time_h=time_h,
                    provenance={"source": "screening_columns"},
                )
                yield_value, warning = optional_float(
                    row.get("Product_Yield_PCT_Area_UV")
                )
                if warning:
                    warnings.append(f"{warning}:Product_Yield_PCT_Area_UV")
                metadata = {}
                for field_name, output_name in (
                    ("Product_UV_Analysis_Wavelength_nm", "wavelength_nm"),
                    ("Product_Yield_Mass_Ion_Count", "mass_ion_count"),
                    ("Product_Selectivity", "selectivity"),
                ):
                    value, numeric_warning = optional_float(row.get(field_name))
                    if numeric_warning:
                        warnings.append(f"{numeric_warning}:{field_name}")
                    if clean_text(row.get(field_name)):
                        metadata[output_name] = value
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
                            "screen_id": clean_text(row.get("SCREEN_ID")),
                            "notebook_id": clean_text(row.get("NOTEBOOK_ID")),
                            "reaction_class": clean_text(row.get("ReactionClass")),
                            "keyword": clean_text(row.get("KeyWord_STD")),
                            "reaction_group": clean_text(row.get("ReactionGroup")),
                        },
                    ),
                    reaction=ReactionEvidenceInput(
                        evidence_kind="source_structure",
                        reaction_smiles=reaction_smiles or None,
                        supplied_mapping_status=supplied_mapping_status(
                            reaction_smiles
                        ),
                        source_reaction_type=clean_text(row.get("ReactionClass")),
                        source_labels={
                            "keyword": clean_text(row.get("KeyWord_STD")),
                            "reaction_group": clean_text(row.get("ReactionGroup")),
                        },
                        structure_available=bool(reaction_smiles),
                    ),
                    conditions=ConditionInput(
                        components=components,
                        component_groups=groups,
                        stages=(stage,),
                        declared_stage_count=1,
                    ),
                    outcomes=(
                        OutcomeInput(
                            outcome_type="uv_area_yield_pct",
                            value=yield_value,
                            unit="percent",
                            raw_value=clean_text(row.get("Product_Yield_PCT_Area_UV")),
                            source_field="Product_Yield_PCT_Area_UV",
                            metadata=metadata,
                        ),
                    ),
                    ingestion_status="accepted" if reaction_smiles else "review",
                    warnings=tuple(sorted(set(warnings))),
                    raw_fields=raw_fields(row),
                )


__all__ = ["HiTeaCsvAdapter"]
