"""Conservative extraction of C-N, C-O, and C-S retrosynthetic templates."""

from __future__ import annotations

import io
from contextlib import redirect_stdout
from dataclasses import dataclass
from typing import Any, Dict, Optional

from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun
from rdchiral.template_extractor import extract_from_reaction

from reactive_taxonomy import featurize_reaction

from ...chemistry import (
    canonical_smiles,
    contributing_reactants,
    digest,
    split_reaction_smiles,
)
from .event_normalization import normalize_single_cx_event
from .models import CxBondKind, TemplatePrecedent
from ...mapping import materialize_atom_mapping


ALLOWED_CORE_EVIDENCE = frozenset({"verified", "inferred"})
ALLOWED_HETEROATOMS = frozenset({"N", "O", "S"})


@dataclass(frozen=True)
class ExtractedTemplate:
    """Internal successful extraction result prior to aggregation."""

    template_id: str
    bond_kind: CxBondKind
    reaction_smarts: str
    product_smarts: str
    precursor_smarts: str
    intra_only: bool
    dimer_only: bool
    precedent: TemplatePrecedent


@dataclass(frozen=True)
class ExtractionResult:
    """One accepted template or one deterministic rejection reason."""

    template: Optional[ExtractedTemplate]
    rejection_reason: Optional[str]


def _observation(row: Dict[str, Any]) -> Dict[str, Any]:
    return dict(row.get("reaction_observation") or row.get("observation") or {})


def _reaction_smiles(row: Dict[str, Any], observation: Dict[str, Any]) -> str:
    return str(
        row.get("reaction_smiles")
        or observation.get("input_reaction_smiles")
        or row.get("input_reaction_smiles")
        or ""
    )


def _reference_id(row: Dict[str, Any]) -> str:
    direct = str(row.get("reference_id") or "").strip()
    if direct:
        return direct
    identity = row.get("reference_identity") or {}
    return str(identity.get("reference_id") or "").strip()


def _formed_bond_kind(observation: Dict[str, Any]) -> Optional[CxBondKind]:
    heavy_formed = []
    for edit in observation.get("edits") or ():
        if edit.get("edit_type") != "formed" or edit.get("atom_2") is None:
            continue
        heavy_formed.append(edit)
    if len(heavy_formed) != 1:
        return None
    edit = heavy_formed[0]
    if str(edit.get("new_order") or "").upper() != "SINGLE":
        return None
    elements = {
        str((edit.get("atom_1") or {}).get("element") or ""),
        str((edit.get("atom_2") or {}).get("element") or ""),
    }
    if "C" not in elements:
        return None
    heteroatoms = elements.intersection(ALLOWED_HETEROATOMS)
    if len(heteroatoms) != 1:
        return None
    return f"C-{next(iter(heteroatoms))}"  # type: ignore[return-value]


def _round_trip(
    reaction_smarts: str,
    product_smiles: str,
    expected_precursors: str,
) -> bool:
    canonical_product = canonical_smiles(product_smiles)
    if canonical_product is None:
        return False
    try:
        with redirect_stdout(io.StringIO()):
            outcomes = rdchiralRun(
                rdchiralReaction(reaction_smarts),
                rdchiralReactants(canonical_product),
            )
    except Exception:
        return False
    return expected_precursors in {
        canonical
        for outcome in outcomes
        if (canonical := canonical_smiles(str(outcome))) is not None
    }


def extract_cx_template(row: Dict[str, Any]) -> ExtractionResult:
    """Extract one C-X rule from a trusted, fully observed reaction row."""

    core = dict(row.get("reaction_core") or {})
    observation = _observation(row)
    if not core or not observation:
        return ExtractionResult(None, "missing_core_or_observation")
    quality_status = (core.get("quality") or {}).get("status")
    if quality_status not in {"pass", "review"}:
        return ExtractionResult(None, "core_quality_blocked")
    if core.get("evidence_status") not in ALLOWED_CORE_EVIDENCE:
        return ExtractionResult(None, "core_evidence_not_allowed")
    if normalize_single_cx_event(core, observation) is None:
        return ExtractionResult(None, "not_single_event")
    completeness = observation.get("completeness") or row.get(
        "reaction_completeness"
    ) or {}
    if completeness.get("status") != "verified":
        return ExtractionResult(None, "product_completeness_not_verified")
    products = observation.get("products") or ()
    if len(products) != 1:
        return ExtractionResult(None, "not_single_product")
    if observation.get("stereo_changes"):
        return ExtractionResult(None, "stereo_change_out_of_scope")
    if any(
        edit.get("edit_type") == "order_changed"
        for edit in observation.get("edits") or ()
    ):
        return ExtractionResult(None, "bond_order_change_out_of_scope")
    bond_kind = _formed_bond_kind(observation)
    if bond_kind is None:
        return ExtractionResult(None, "not_single_cx_bond_formation")

    source_reaction = _reaction_smiles(row, observation)
    materialized = materialize_atom_mapping(source_reaction)
    if materialized is None:
        return ExtractionResult(None, "atom_mapping_unavailable")
    mapped_reaction = materialized.reaction_smiles
    split = split_reaction_smiles(mapped_reaction)
    if split is None:
        return ExtractionResult(None, "invalid_reaction_smiles")
    reactants, product = split
    expected_precursors = contributing_reactants(reactants, product)
    canonical_product = canonical_smiles(product)
    if expected_precursors is None or canonical_product is None:
        return ExtractionResult(None, "participant_canonicalization_failed")
    mapped_analysis = featurize_reaction(mapped_reaction)
    mapped_core = mapped_analysis.reaction_core
    if (
        not mapped_analysis.valid
        or mapped_core is None
        or mapped_core.quality.status != "pass"
    ):
        return ExtractionResult(None, "materialized_mapping_not_verified")
    if mapped_core.center_transition_key != core.get("center_transition_key"):
        return ExtractionResult(None, "materialized_mapping_core_conflict")

    try:
        with redirect_stdout(io.StringIO()):
            raw = extract_from_reaction(
                {"_id": "cx_poc", "reactants": reactants, "products": product}
            )
    except Exception:
        return ExtractionResult(None, "template_extraction_failed")
    if not raw or not raw.get("reaction_smarts"):
        return ExtractionResult(None, "template_extraction_failed")
    reaction_smarts = str(raw["reaction_smarts"])
    if not _round_trip(reaction_smarts, product, expected_precursors):
        return ExtractionResult(None, "source_round_trip_failed")

    template_id = digest("CXT1", bond_kind, reaction_smarts)
    reference_id = _reference_id(row)
    support_identity = str(
        core.get("mapping_equivalence_key")
        or core.get("core_id")
        or row.get("reaction_id")
        or row.get("observation_id")
        or template_id
    )
    precedent = TemplatePrecedent(
        reaction_id=str(row.get("reaction_id") or ""),
        observation_id=str(row.get("observation_id") or ""),
        reference_id=reference_id,
        support_unit_id=(
            f"reference:{reference_id}"
            if reference_id
            else f"mapping_equivalence:{support_identity}"
        ),
        core_id=str(core.get("core_id") or ""),
        mapping_evidence=materialized.evidence,
        mapping_confidence=materialized.confidence,
        product_smiles=canonical_product,
        precursor_smiles=expected_precursors,
        mapped_reaction_smiles=mapped_reaction,
    )
    return ExtractionResult(
        ExtractedTemplate(
            template_id=template_id,
            bond_kind=bond_kind,
            reaction_smarts=reaction_smarts,
            product_smarts=str(raw["products"]),
            precursor_smarts=str(raw["reactants"]),
            intra_only=bool(raw.get("intra_only")),
            dimer_only=bool(raw.get("dimer_only")),
            precedent=precedent,
        ),
        None,
    )


__all__ = ["ExtractionResult", "ExtractedTemplate", "extract_cx_template"]
