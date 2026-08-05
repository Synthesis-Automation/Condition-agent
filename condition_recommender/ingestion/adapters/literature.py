"""Adapter for the literature reaction-dataset CSV contract."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterator, List, Mapping, Sequence, Tuple

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
    """Normalize the structured literature CSV without interpreting chemistry."""

    adapter_id = "literature_csv.v2"
    adapter_version = "2.0"
    corpus_id = "literature_reaction_dataset"
    required_columns = (
        "reaction_id",
        "reaction_type",
        "yield_pct",
        "reactant_smiles",
        "product_smiles",
        "reaction_smiles",
        "reference",
        "reagent_cas",
        "catalyst_cas",
        "solvent_cas",
        "product_yields_json",
        "reactants_json",
        "products_json",
        "reagents_json",
        "catalysts_json",
        "solvents_json",
    )

    @staticmethod
    def _json_value(
        row: Mapping[str, str],
        field_name: str,
        expected_type: type,
        warnings: List[str],
    ) -> tuple[Any, bool]:
        """Decode one structured source field and retain invalid rows for review."""
        raw_value = clean_text(row.get(field_name))
        if not raw_value:
            warnings.append(f"MISSING_STRUCTURED_FIELD:{field_name}")
            return expected_type(), False
        try:
            value = json.loads(raw_value)
        except json.JSONDecodeError:
            warnings.append(f"INVALID_JSON:{field_name}")
            return expected_type(), False
        if not isinstance(value, expected_type):
            warnings.append(f"INVALID_JSON_TYPE:{field_name}")
            return expected_type(), False
        return value, True

    @classmethod
    def _component_records(
        cls,
        row: Mapping[str, str],
        field_name: str,
        warnings: List[str],
    ) -> tuple[Tuple[Mapping[str, Any], ...], bool]:
        value, valid = cls._json_value(row, field_name, list, warnings)
        records: List[Mapping[str, Any]] = []
        for position, item in enumerate(value, start=1):
            if not isinstance(item, Mapping):
                warnings.append(f"INVALID_COMPONENT_RECORD:{field_name}:{position}")
                valid = False
                continue
            records.append(item)
        return tuple(records), valid

    @staticmethod
    def _values(value: Any) -> Tuple[str, ...]:
        """Normalize a scalar or list-valued structured source attribute."""
        source_values: Sequence[Any] = value if isinstance(value, list) else (value,)
        return tuple(cleaned for item in source_values if (cleaned := clean_text(item)))

    @classmethod
    def _components(
        cls,
        row: Mapping[str, str],
        structured: Mapping[str, Tuple[Mapping[str, Any], ...]],
        structured_validity: Mapping[str, bool],
        warnings: List[str],
    ) -> Tuple[ConditionComponentClaim, ...]:
        components: List[ConditionComponentClaim] = []
        for json_field, flat_field, role in (
            ("catalysts_json", "catalyst_cas", "catalyst"),
            ("reagents_json", "reagent_cas", "reagent"),
            ("solvents_json", "solvent_cas", "solvent"),
        ):
            records = structured[json_field]
            structured_cas: List[str] = []
            role_components: List[ConditionComponentClaim] = []
            for source_position, record in enumerate(records, start=1):
                cas_values = cls._values(record.get("cas_rn"))
                structured_cas.extend(cas_values)
                if not cas_values:
                    warnings.append(
                        f"MISSING_COMPONENT_IDENTIFIER:{json_field}:{source_position}"
                    )
                    continue
                source_index, index_warning = optional_int(record.get("index"))
                if index_warning or source_index is None or source_index < 1:
                    warnings.append(
                        f"INVALID_COMPONENT_INDEX:{json_field}:{source_position}"
                    )
                    source_index = source_position
                identifier_field = f"{json_field}[{source_index}].cas_rn"
                provenance: Dict[str, Any] = {
                    "source_position": source_position,
                    "source_component_index": source_index,
                }
                amd_values = cls._values(record.get("amd"))
                if amd_values:
                    provenance["additional_material_description"] = amd_values
                role_components.append(
                    ConditionComponentClaim(
                        component_key=f"{role}_{source_index}",
                        source_slot=json_field,
                        source_role_hint=role,
                        identifiers=tuple(
                            SourceIdentifier("cas", value, identifier_field)
                            for value in cas_values
                        ),
                        provenance=provenance,
                    )
                )

            flat_cas = list(split_identifiers(row.get(flat_field)))
            if structured_cas != flat_cas:
                warnings.append(f"CONFLICTING_STRUCTURED_FIELD:{json_field}:{flat_field}")
            if role_components or (structured_validity[json_field] and not flat_cas):
                components.extend(role_components)
                continue

            warnings.append(f"STRUCTURED_CONDITIONS_FALLBACK:{json_field}")
            for position, value in enumerate(flat_cas, start=1):
                components.append(
                    ConditionComponentClaim(
                        component_key=f"{role}_{position}",
                        source_slot=flat_field,
                        source_role_hint=role,
                        identifiers=(SourceIdentifier("cas", value, flat_field),),
                        provenance={
                            "source_position": position,
                            "structured_field": json_field,
                            "structured_fallback": True,
                        },
                    )
                )
        return tuple(components)

    @classmethod
    def _outcomes(
        cls,
        row: Mapping[str, str],
        products: Tuple[Mapping[str, Any], ...],
        warnings: List[str],
    ) -> Tuple[OutcomeInput, ...]:
        outcomes: List[OutcomeInput] = []
        raw_yield = clean_text(row.get("yield_pct"))
        if raw_yield:
            value, warning = optional_float(raw_yield)
            if warning:
                warnings.append(f"{warning}:yield_pct")
            outcomes.append(
                OutcomeInput(
                    outcome_type="reported_yield_pct",
                    value=value,
                    unit="percent",
                    raw_value=raw_yield,
                    source_field="yield_pct",
                    metadata={"measurement_basis": "source_unspecified"},
                )
            )

        yield_mapping, yields_valid = cls._json_value(
            row, "product_yields_json", dict, warnings
        )
        products_by_index = {
            clean_text(product.get("index")): product
            for product in products
            if clean_text(product.get("index"))
        }
        product_yields: List[tuple[str, Any]] = []
        if yields_valid:
            product_yields = sorted(
                ((clean_text(key), value) for key, value in yield_mapping.items()),
                key=lambda item: (
                    0 if item[0].isdigit() else 1,
                    int(item[0]) if item[0].isdigit() else item[0],
                ),
            )
        else:
            warnings.append("STRUCTURED_OUTCOMES_FALLBACK:product_yields_json")
            product_yields = [
                (str(index), row.get(f"product_yield_{index}"))
                for index in range(1, 8)
                if clean_text(row.get(f"product_yield_{index}"))
            ]

        for product_index, source_value in product_yields:
            raw_value = clean_text(source_value)
            if not raw_value:
                continue
            value, warning = optional_float(raw_value)
            if warning:
                warnings.append(
                    f"{warning}:product_yields_json[{product_index}]"
                )
            product = products_by_index.get(product_index, {})
            product_raw_yield = clean_text(product.get("yield_pct"))
            if product_raw_yield and product_raw_yield != raw_value:
                warnings.append(
                    f"CONFLICTING_PRODUCT_YIELD:products_json:{product_index}"
                )
            flat_field = (
                f"product_yield_{product_index}"
                if product_index.isdigit() and 1 <= int(product_index) <= 7
                else ""
            )
            flat_value = clean_text(row.get(flat_field)) if flat_field else ""
            if flat_value and flat_value != raw_value:
                warnings.append(
                    f"CONFLICTING_PRODUCT_YIELD:{flat_field}:{product_index}"
                )
            metadata: Dict[str, Any] = {"product_position": product_index}
            product_cas = cls._values(product.get("cas_rn"))
            if product_cas:
                metadata["product_cas"] = product_cas
            product_smiles = clean_text(product.get("smiles"))
            if product_smiles:
                metadata["product_smiles"] = product_smiles
            outcomes.append(
                OutcomeInput(
                    outcome_type="reported_product_yield_pct",
                    value=value,
                    unit="percent",
                    raw_value=raw_value,
                    source_field=f"product_yields_json[{product_index}]",
                    metadata=metadata,
                )
            )
        return tuple(outcomes)

    @staticmethod
    def _structured_side_smiles(
        records: Tuple[Mapping[str, Any], ...],
    ) -> str:
        smiles = tuple(clean_text(record.get("smiles")) for record in records)
        return ".".join(smiles) if smiles and all(smiles) else ""

    def iter_observations(
        self, path: Path, *, source_sha256: str
    ) -> Iterator[CanonicalSourceObservation]:
        """Stream one normalized observation for every literature row."""
        validate_headers(path, self.required_columns)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row_number, row in enumerate(csv.DictReader(handle), start=2):
                warnings: List[str] = []
                structured: Dict[str, Tuple[Mapping[str, Any], ...]] = {}
                structured_validity: Dict[str, bool] = {}
                for field_name in (
                    "reactants_json",
                    "products_json",
                    "reagents_json",
                    "catalysts_json",
                    "solvents_json",
                ):
                    records, valid = self._component_records(
                        row, field_name, warnings
                    )
                    structured[field_name] = records
                    structured_validity[field_name] = valid

                reaction_smiles = clean_text(row.get("reaction_smiles"))
                reactant_smiles = clean_text(row.get("reactant_smiles"))
                product_smiles = clean_text(row.get("product_smiles"))
                structured_reactants = self._structured_side_smiles(
                    structured["reactants_json"]
                )
                structured_products = self._structured_side_smiles(
                    structured["products_json"]
                )
                if structured_reactants and reactant_smiles != structured_reactants:
                    warnings.append(
                        "CONFLICTING_STRUCTURED_FIELD:reactants_json:reactant_smiles"
                    )
                if structured_products and product_smiles != structured_products:
                    warnings.append(
                        "CONFLICTING_STRUCTURED_FIELD:products_json:product_smiles"
                    )
                side_reaction = (
                    f"{reactant_smiles}>>{product_smiles}"
                    if reactant_smiles and product_smiles
                    else ""
                )
                if not reaction_smiles and side_reaction:
                    reaction_smiles = side_reaction
                    warnings.append("REACTION_SMILES_RECONSTRUCTED_FROM_SIDES")
                elif reaction_smiles and side_reaction and reaction_smiles != side_reaction:
                    warnings.append("CONFLICTING_REACTION_SMILES_SIDES")
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
                components = self._components(
                    row, structured, structured_validity, warnings
                )
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
                outcomes = self._outcomes(
                    row, structured["products_json"], warnings
                )
                source_structure_warnings = tuple(
                    item.strip()
                    for item in clean_text(row.get("structure_warnings")).split("|")
                    if item.strip()
                )
                warnings.extend(
                    f"SOURCE_STRUCTURE_WARNING:{item}"
                    for item in source_structure_warnings
                )
                if not reaction_smiles:
                    warnings.append("MISSING_REACTION_SMILES")
                requires_review = any(
                    item == "MISSING_REACTION_SMILES"
                    or item == "CONFLICTING_REACTION_SMILES_SIDES"
                    or item.startswith("CONFLICTING_STRUCTURED_FIELD:")
                    or item.startswith("INVALID_JSON")
                    or item.startswith("SOURCE_STRUCTURE_WARNING:")
                    for item in warnings
                )
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
                            "duplicate_reaction_ids": clean_text(
                                row.get("duplicate_reaction_ids")
                            ),
                            "source_files": clean_text(row.get("source_files")),
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
                        source_labels={
                            "title": clean_text(row.get("title")),
                            "authors": clean_text(row.get("authors")),
                            "citation": clean_text(row.get("citation")),
                            "reactants": tuple(
                                dict(item) for item in structured["reactants_json"]
                            ),
                            "products": tuple(
                                dict(item) for item in structured["products_json"]
                            ),
                        },
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
                    outcomes=outcomes,
                    ingestion_status="review" if requires_review else "accepted",
                    warnings=tuple(sorted(set(warnings))),
                    raw_fields=raw_fields(row),
                )


__all__ = ["LiteratureCsvAdapter"]
