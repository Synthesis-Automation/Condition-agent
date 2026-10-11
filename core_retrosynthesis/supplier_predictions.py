"""Structure-normalized external proposals, without implied reaction validation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from typing import Any, Mapping, Sequence

from rdkit import Chem


def _audit_value(value: Any) -> Any:
    """Retain malformed numeric tokens without emitting nonstandard JSON."""
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    if isinstance(value, dict):
        return {key: _audit_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_audit_value(item) for item in value]
    return value


def canonical_supplier_smiles(smiles: str) -> str:
    """Normalize a parsed structure, preserving stereochemistry and multiplicity."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("A nonempty SMILES is required")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not mol.GetNumAtoms():
        raise ValueError(f"Invalid molecular structure: {smiles}")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


@dataclass(frozen=True)
class SupplierPrediction:
    """A parseable hypothesis; scores do not establish feasibility or evidence."""

    candidate_id: str
    supplier: str
    supplier_rank: int
    target_smiles: str
    precursor_smiles: tuple[str, ...]
    proposed_reaction_smiles: str
    model_score: float | None
    model_probability: float | None
    component_ranks: tuple[tuple[str, int | None], ...]
    evidence_kind: str = "predicted"
    validation_status: str = "structures_parsed_reaction_unverified"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible proposal with explicit model-only scores."""
        value = asdict(self)
        value["precursor_smiles"] = list(self.precursor_smiles)
        value["component_ranks"] = dict(self.component_ranks)
        return value


def normalize_supplier_predictions(
    target_smiles: str, rows: Sequence[Mapping[str, Any]], *, supplier: str = "retrochimera",
) -> dict[str, Any]:
    """Deduplicate precursor multisets, retaining rejected rows and original ranks.

    Parsing is the only chemical gate here. No atom correspondence, feasibility,
    literature support, or taxonomy operator is inferred from a model score.
    """
    target = canonical_supplier_smiles(target_smiles)
    candidates, rejected, duplicates = [], [], []
    seen: dict[tuple[str, ...], str] = {}
    for rank, row in enumerate(rows, 1):
        try:
            if canonical_supplier_smiles(row["product"]) != target:
                raise ValueError("Supplier product differs from requested target")
            raw = row["reactants"]
            if not isinstance(raw, list) or not raw:
                raise ValueError("Supplier must return a nonempty reactant list")
            precursors = tuple(sorted(
                part for item in raw for part in canonical_supplier_smiles(item).split(".")
            ))
            if target in precursors:
                raise ValueError("Proposal contains unchanged target as a precursor")
            scores = []
            for key in ("score", "probability"):
                value = row.get(key)
                if value is not None:
                    if isinstance(value, bool) or not isinstance(value, (int, float)):
                        raise ValueError(f"Invalid model {key}")
                    value = float(value)
                    if not math.isfinite(value) or (key == "probability" and not 0 <= value <= 1):
                        raise ValueError(f"Invalid model {key}")
                scores.append(value)
            ranks = row.get("individual_ranks", {})
            if not isinstance(ranks, dict) or any(
                not isinstance(k, str) or (v is not None and (type(v) is not int or v < 0))
                for k, v in ranks.items()
            ):
                raise ValueError("Invalid zero-based component ranks")
            identity = json.dumps(["supplier_prediction.v1", supplier, target, precursors])
            candidate_id = hashlib.sha256(identity.encode()).hexdigest()
            if precursors in seen:
                duplicates.append({"supplier_rank": rank, "candidate_id": seen[precursors], "raw": dict(row)})
                continue
            seen[precursors] = candidate_id
            candidates.append(SupplierPrediction(
                candidate_id, supplier, rank, target, precursors,
                ".".join(precursors) + ">>" + target, scores[0], scores[1], tuple(sorted(ranks.items())),
            ).to_dict())
        except (ValueError, KeyError, TypeError) as exc:
            rejected.append({"supplier_rank": rank, "reason": str(exc), "raw": _audit_value(dict(row))})
    return {
        "schema_version": "supplier_predictions.v1", "supplier": supplier,
        "target_smiles": target, "candidates": candidates, "rejected": rejected,
        "duplicates": duplicates, "raw_candidate_count": len(rows),
        "limitations": [
            "Predicted reactions; structure parsing does not verify the transformation.",
            "Model probability is shortlist-normalized, not reaction success or yield.",
            "Component ranks are zero-based; supplier ranks are one-based.",
            "Precedent support, conditions, selectivity and starting-material access remain unassessed.",
        ],
    }
