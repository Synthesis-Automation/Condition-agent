"""Adapter for the cleaned USPTO reaction-condition CSV contract."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator, List, Mapping, Tuple

from ..models import (
    CanonicalSourceObservation,
    ConditionComponentClaim,
    ConditionInput,
    ConditionStageInput,
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


_ABSENT_VALUES = {"", "none", "null", "n/a", "na", "nan"}


def _present(value: object) -> bool:
    return clean_text(value).casefold() not in _ABSENT_VALUES


def _complete_reaction_smiles(value: str) -> bool:
    """Check only that source text supplies nonempty reactant and product sides."""
    if ">>" in value:
        reactants, products = value.split(">>", 1)
        return bool(reactants.strip() and products.strip())
    parts = value.split(">")
    return len(parts) == 3 and bool(parts[0].strip() and parts[2].strip())


class UsptoConditionCsvAdapter:
    """Preserve USPTO mapped structures and SMILES condition claims."""

    adapter_id = "uspto_condition_csv.v1"
    adapter_version = "1.0"
    corpus_id = "uspto_condition_reactions_cleaned"
    required_columns = (
        "source",
        "canonical_rxn",
        "catalyst1",
        "solvent1",
        "solvent2",
        "reagent1",
        "reagent2",
        "remapped_rxn",
        "confidence",
    )

    @staticmethod
    def _condition_claims(
        row: Mapping[str, str],
    ) -> Tuple[ConditionComponentClaim, ...]:
        components: List[ConditionComponentClaim] = []
        for role, fields in (
            ("catalyst", ("catalyst1",)),
            ("solvent", ("solvent1", "solvent2")),
            ("reagent", ("reagent1", "reagent2")),
        ):
            for position, field_name in enumerate(fields, start=1):
                value = clean_text(row.get(field_name))
                if not _present(value):
                    continue
                components.append(
                    ConditionComponentClaim(
                        component_key=f"{role}_{position}",
                        source_slot=field_name,
                        source_role_hint=role,
                        identifiers=(
                            SourceIdentifier("smiles", value, field_name),
                        ),
                        provenance={
                            "source_position": position,
                            "source_declared_role": role,
                        },
                    )
                )
        return tuple(components)

    @staticmethod
    def _select_reaction_smiles(
        row: Mapping[str, str], warnings: List[str]
    ) -> tuple[str, str]:
        mapped = clean_text(row.get("remapped_rxn"))
        canonical = clean_text(row.get("canonical_rxn"))
        mapped_status = supplied_mapping_status(mapped)
        if (
            mapped_status == "supplied_unvalidated"
            and _complete_reaction_smiles(mapped)
        ):
            return mapped, "remapped_rxn"

        if not mapped:
            warnings.append("MISSING_REMAPPED_REACTION_SMILES")
        elif mapped_status != "supplied_unvalidated":
            warnings.append("REMAPPED_REACTION_HAS_NO_ATOM_MAPS")
        else:
            warnings.append("INCOMPLETE_REMAPPED_REACTION_SMILES")

        if canonical and _complete_reaction_smiles(canonical):
            warnings.append("CANONICAL_REACTION_SMILES_FALLBACK")
            return canonical, "canonical_rxn"
        if canonical:
            warnings.append("INCOMPLETE_CANONICAL_REACTION_SMILES")
        else:
            warnings.append("MISSING_CANONICAL_REACTION_SMILES")
        warnings.append("MISSING_REACTION_SMILES")
        return "", ""

    def iter_observations(
        self, path: Path, *, source_sha256: str
    ) -> Iterator[CanonicalSourceObservation]:
        """Stream normalized observations from the cleaned USPTO dataset."""
        validate_headers(path, self.required_columns)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row_number, row in enumerate(csv.DictReader(handle), start=2):
                warnings: List[str] = []
                patent_id = clean_text(row.get("source"))
                record_id = (
                    f"{patent_id}:row-{row_number}"
                    if patent_id
                    else f"{path.stem}:row-{row_number}"
                )
                reaction_smiles, selected_field = self._select_reaction_smiles(
                    row, warnings
                )
                confidence, warning = optional_float(row.get("confidence"))
                if warning:
                    warnings.append(f"{warning}:confidence")
                elif confidence is not None and not 0.0 <= confidence <= 1.0:
                    warnings.append("MAPPING_CONFIDENCE_OUT_OF_RANGE")

                components = self._condition_claims(row)
                stage = ConditionStageInput(
                    stage_index=1,
                    component_keys=tuple(
                        component.component_key for component in components
                    ),
                    provenance={"source": "condition_smiles_columns"},
                )
                requires_review = not reaction_smiles or any(
                    item.startswith("INCOMPLETE_")
                    or item.startswith("INVALID_")
                    or item == "MAPPING_CONFIDENCE_OUT_OF_RANGE"
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
                            "patent_id": patent_id,
                            "mapping_confidence": clean_text(
                                row.get("confidence")
                            ),
                        },
                        reference=patent_id,
                    ),
                    reaction=ReactionEvidenceInput(
                        evidence_kind="source_structure",
                        reaction_smiles=reaction_smiles or None,
                        supplied_mapping_status=supplied_mapping_status(
                            reaction_smiles
                        ),
                        source_labels={
                            "reaction_smiles_source_field": selected_field,
                            "mapping_confidence": confidence,
                            "canonical_reaction_smiles": clean_text(
                                row.get("canonical_rxn")
                            ),
                        },
                        structure_available=bool(reaction_smiles),
                    ),
                    conditions=ConditionInput(
                        components=components,
                        stages=(stage,),
                    ),
                    ingestion_status="review" if requires_review else "accepted",
                    warnings=tuple(sorted(set(warnings))),
                    raw_fields=raw_fields(row),
                )


__all__ = ["UsptoConditionCsvAdapter"]
