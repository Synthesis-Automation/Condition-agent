"""Conditional mass/amount checks for explicitly assigned source structures.

Only immediate, complete parenthetical mass/amount pairs are extracted. Literal
labels locate a passage; they never resolve identities or verify assignments.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Literal

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors

_NUMBER = r"[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?"


@lru_cache(maxsize=1)
def quantity_consistency_rules() -> dict[str, Any]:
    """Load validated, versioned units and discrepancy tolerance."""
    path = Path(__file__).with_name("definitions") / "quantity_consistency.v1.json"
    rules = json.loads(path.read_text(encoding="utf-8"))
    if (rules.get("schema_version") != "quantity_consistency.v1"
            or rules.get("policy") != "retain_reports_without_selecting_or_correcting_values"
            or not rules.get("definition_version")):
        raise ValueError("Unsupported quantity consistency definition")
    tolerance = rules.get("relative_tolerance")
    if type(tolerance) not in (int, float) or not math.isfinite(tolerance) or not 0 < tolerance < 1:
        raise ValueError("Quantity consistency tolerance must be between zero and one")
    for key in ("mass_units_g", "amount_units_mol"):
        units = rules.get(key)
        if not isinstance(units, dict) or not units:
            raise ValueError("Quantity consistency units must be nonempty mappings")
        for unit, factor in units.items():
            if (not isinstance(unit, str) or not unit or type(factor) not in (int, float)
                    or not math.isfinite(factor) or factor <= 0):
                raise ValueError("Invalid quantity consistency conversion")
    return rules


@dataclass(frozen=True)
class SourceQuantityCheck:
    """Reported quantities and conditional graph arithmetic, never corrected data."""

    status: Literal["consistent", "conflicting", "not_assessed"]
    matched_identifier: str
    source_text: str
    source_start: int
    source_end: int
    reported_mass: float | None
    mass_unit: str
    reported_amount: float | None
    amount_unit: str
    molecular_weight: float | None
    expected_mass_g: float | None
    relative_error: float | None
    reason: str
    definition_version: str
    schema_version: str = "source_quantity_check.v1"


def audit_source_quantities(
    smiles: str, passage: str, identifiers: tuple[str, ...],
) -> tuple[SourceQuantityCheck, ...]:
    """Check exact adjacent reports against the supplied complete molecular graph.

    Handles mass-first and amount-first pairs, including SI prefixes. Unsupported
    units, solution/concentration annotations and nonadjacent values are not
    guessed. Empty results mean no eligible pair was found, not that the source
    is consistent. Salts/isotopes remain in the graph; mixtures and invalid or
    wildcard structures retain reports with not_assessed status.
    """
    rules = quantity_consistency_rules()
    digest = hashlib.sha256(json.dumps(rules, sort_keys=True).encode()).hexdigest()[:16]
    version = f"{rules['definition_version']}@sha256:{digest}"
    units = sorted({*rules["mass_units_g"], *rules["amount_units_mol"]}, key=lambda x: (-len(x), x))
    unit_pattern = "|".join(re.escape(unit).replace("u", "[uµμ]") for unit in units)
    pair = rf"\s*\(\s*({_NUMBER})\s*({unit_pattern})\s*,\s*({_NUMBER})\s*({unit_pattern})\s*\)"
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles)
    valid = molecule is not None and all(
        atom.GetAtomicNum() > 0 and not atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()
    )
    # Multiple carbon-containing fragments may be a mixture/solvate with unknown
    # composition. Do not assume one mole of the concatenated graph.
    mixture = valid and sum(
        any(atom.GetAtomicNum() == 6 for atom in fragment.GetAtoms())
        for fragment in Chem.GetMolFrags(molecule, asMols=True)
    ) > 1
    mw = Descriptors.MolWt(molecule) if valid and not mixture else None
    checks: dict[tuple[int, int], SourceQuantityCheck] = {}
    for identifier in dict.fromkeys(identifiers):
        if not identifier.strip():
            continue
        pattern = r"(?<!\w)" + re.escape(identifier) + r"(?!\w)" + pair
        for match in re.finditer(pattern, passage, re.IGNORECASE):
            first, first_unit, second, second_unit = match.groups()
            first_unit = first_unit.casefold().replace("µ", "u").replace("μ", "u")
            second_unit = second_unit.casefold().replace("µ", "u").replace("μ", "u")
            if first_unit in rules["amount_units_mol"]:
                first, second, first_unit, second_unit = second, first, second_unit, first_unit
            if first_unit not in rules["mass_units_g"] or second_unit not in rules["amount_units_mol"]:
                continue
            mass, amount = float(first), float(second)
            eligible = mw is not None and all(math.isfinite(value) and value > 0 for value in (mass, amount))
            expected = amount * rules["amount_units_mol"][second_unit] * mw if eligible else None
            if expected is not None and (not math.isfinite(expected) or expected <= 0):
                eligible, expected = False, None
            error = abs(mass * rules["mass_units_g"][first_unit] - expected) / expected if eligible else None
            if error is not None and not math.isfinite(error):
                eligible, error = False, None
            status = ("conflicting" if error > rules["relative_tolerance"] else "consistent") if eligible else "not_assessed"
            key = (match.start(1), match.end())
            if key in checks:
                continue
            checks[key] = SourceQuantityCheck(
                status, identifier, match.group(), match.start(), match.end(),
                mass if math.isfinite(mass) else None, first_unit,
                amount if math.isfinite(amount) else None, second_unit, mw, expected, error,
                "Conditional on supplied graph/material assignment; both reports are retained."
                if eligible else "Invalid/ambiguous graph, mixture or invalid quantity; no conversion accepted.",
                version,
            )
    return tuple(checks[key] for key in sorted(checks))
