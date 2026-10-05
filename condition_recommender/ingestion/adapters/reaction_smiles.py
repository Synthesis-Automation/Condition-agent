"""Source adapters for the higher-level release's two reaction CSV exports."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator

from ..models import CanonicalSourceObservation, ConditionInput, ReactionEvidenceInput
from .base import (
    clean_text,
    observation_id,
    raw_fields,
    source_provenance,
    supplied_mapping_status,
    validate_headers,
)


class ReactionSmilesCsvAdapter:
    """Preserve explicit reaction strings without interpreting chemical validity."""

    adapter_id = "reaction_smiles_csv.v1"
    adapter_version = "1.0"
    corpus_id = "reaction_smiles_export"
    required_columns = ("id", "rxn_smiles")

    def iter_observations(
        self, path: Path, *, source_sha256: str
    ) -> Iterator[CanonicalSourceObservation]:
        """Keep every source row, including missing structures for review."""
        validate_headers(path, self.required_columns)
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for row_number, row in enumerate(csv.DictReader(handle), start=2):
                record_id = clean_text(row.get("id")) or f"row-{row_number}"
                reaction = clean_text(row.get("rxn_smiles"))
                abstraction = isinstance(self, HigherLevelAbstractionCsvAdapter)
                warnings = (
                    ("ALGORITHMIC_ABSTRACTION_NOT_OBSERVED_REACTION",)
                    if abstraction
                    else (("MISSING_REACTION_SMILES",) if not reaction else ())
                )
                yield CanonicalSourceObservation(
                    observation_id=observation_id(
                        adapter_id=self.adapter_id,
                        source_sha256=source_sha256,
                        row_number=row_number,
                        record_id=record_id,
                    ),
                    observation_kind="algorithmic_abstraction"
                    if abstraction
                    else "structure_backed",
                    source=source_provenance(
                        adapter=self,
                        path=path,
                        source_sha256=source_sha256,
                        row_number=row_number,
                        record_id=record_id,
                    ),
                    reaction=ReactionEvidenceInput(
                        evidence_kind="source_abstraction"
                        if abstraction
                        else "source_structure",
                        # Abstract isotope tags encode the author's abstractions;
                        # they must never become an observed molecular reaction.
                        reaction_smiles=None if abstraction else reaction or None,
                        structure_available=bool(reaction) and not abstraction,
                        supplied_mapping_status=(
                            "not_applicable"
                            if abstraction
                            else supplied_mapping_status(reaction)
                        ),
                        source_labels={"abstracted_reaction_smiles": reaction}
                        if abstraction
                        else {},
                    ),
                    conditions=ConditionInput(
                        warnings=("CONDITION_DETAILS_UNAVAILABLE",)
                    ),
                    ingestion_status="review" if warnings else "accepted",
                    warnings=warnings,
                    raw_fields=raw_fields(row),
                )


class HigherLevelAbstractionCsvAdapter(ReactionSmilesCsvAdapter):
    """Archive algorithmic abstractions separately from physical reaction evidence."""

    adapter_id = "higher_level_abstraction_csv.v1"
    corpus_id = "higher_level_reaction_abstractions"
    artifact_suffix = ".abstractions.jsonl.gz"
