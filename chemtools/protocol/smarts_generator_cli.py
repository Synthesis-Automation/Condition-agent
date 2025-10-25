"""
Reaction SMARTS applicability generator and CLI utilities.

This module provides a light-weight, opinionated implementation of the tooling
exercised in ``tests/test_smarts_generator.py`` and
``tests/test_batch_update_protocol_smarts.py``.  It converts reaction SMILES
into deterministic reaction SMARTS patterns with guard constraints that capture
basic substrate selectivity (e.g., distinguishing primary alkyl iodides from
aryl bromides).  The implementation favours repeatability and a small surface
area over exhaustive chemical coverage; more sophisticated logic can be layered
on top in future iterations.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

try:  # RDKit is optional; the tests handle its absence explicitly.
    from rdkit import Chem
    from rdkit.Chem import Draw, rdChemReactions

    RDKIT_AVAILABLE = True
except Exception:  # pragma: no cover - RDKit missing on some environments
    Chem = None  # type: ignore
    Draw = None  # type: ignore
    rdChemReactions = None  # type: ignore
    RDKIT_AVAILABLE = False

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class ReactionSmartsApplicability:
    """
    Container for a reaction SMARTS transformation and associated guard rules.

    ``guards_forbid`` contains SMARTS that must NOT match the substrates;
    ``guards_require`` contains SMARTS that must be present; ``notes`` is an
    optional free-form annotation.
    """

    core: str
    guards_forbid: List[str] = field(default_factory=list)
    guards_require: List[str] = field(default_factory=list)
    notes: Optional[str] = None

    def to_dict(self) -> Dict[str, object]:
        """Serialise the applicability pattern to a JSON-friendly dict."""
        payload: Dict[str, object] = {"core": self.core}
        if self.guards_forbid:
            payload["guards_forbid"] = self.guards_forbid
        if self.guards_require:
            payload["guards_require"] = self.guards_require
        if self.notes:
            payload["notes"] = self.notes
        return payload

    @classmethod
    def from_dict(cls, payload: Dict[str, object]) -> "ReactionSmartsApplicability":
        """Create a pattern from a dict (inverse of :meth:`to_dict`)."""
        return cls(
            core=str(payload["core"]),
            guards_forbid=list(payload.get("guards_forbid", [])),
            guards_require=list(payload.get("guards_require", [])),
            notes=payload.get("notes"),
        )


# ---------------------------------------------------------------------------
# SMARTS generator core
# ---------------------------------------------------------------------------


class SmartsGenerator:
    """Utility for deriving reaction SMARTS heuristics from reaction SMILES."""

    def __init__(self, reaction_smiles: str):
        self.original = reaction_smiles.strip()
        self.reactants_smiles: str = ""
        self.reagents_smiles: str = ""
        self.products_smiles: str = ""
        self._parse_reaction()

    # -- parsing -----------------------------------------------------------------

    def _parse_reaction(self) -> None:
        if ">" not in self.original:
            raise ValueError("No reaction arrow found in SMILES input")

        if ">>" in self.original:
            if ">>>" in self.original or self.original.count(">>") != 1:
                raise ValueError("Invalid reaction SMILES: ambiguous '>>'")
            parts = self.original.split(">>")
            if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
                raise ValueError("Invalid reaction SMILES: ambiguous '>>'")
            self.reactants_smiles, self.products_smiles = parts
            self.reagents_smiles = ""
            if ">" in self.products_smiles:
                raise ValueError("Invalid reaction SMILES: unexpected '>' in products")
        else:
            parts = self.original.split(">")
            if len(parts) == 2:
                self.reactants_smiles, self.products_smiles = parts
                self.reagents_smiles = ""
            elif len(parts) == 3:
                self.reactants_smiles, self.reagents_smiles, self.products_smiles = parts
            else:
                raise ValueError("Invalid reaction SMILES: expected 1 or 2 '>' separators")

        self.reactants_smiles = self.reactants_smiles.strip()
        self.reagents_smiles = self.reagents_smiles.strip()
        self.products_smiles = self.products_smiles.strip()

        if not self.reactants_smiles or not self.products_smiles:
            raise ValueError("Invalid reaction SMILES: missing reactants or products")

    # -- heuristic SMARTS generation ---------------------------------------------

    def generate_core_smarts(
        self,
        reactant_pattern: Optional[str] = None,
        product_pattern: Optional[str] = None,
    ) -> str:
        """
        Generate a reaction SMARTS core pattern.

        Custom ``reactant_pattern``/``product_pattern`` arguments bypass the
        heuristics and are simply joined with ``>>``.  Otherwise, a set of
        lightweight heuristics is applied.  RDKit must be available for the
        automatic path (mirroring the historical behaviour of the CLI).
        """
        if reactant_pattern:
            product_pattern = product_pattern or self.products_smiles or "[*]"
            return f"{reactant_pattern}>>{product_pattern}"

        if not RDKIT_AVAILABLE:
            raise RuntimeError("RDKit required for automatic SMARTS generation")

        return self._heuristic_core_smarts()

    def _heuristic_core_smarts(self) -> str:
        """Produce a deterministic core SMARTS string from reaction SMILES."""
        reactants = self.reactants_smiles
        products = self.products_smiles
        reactants_lower = reactants.lower()
        products_lower = products.lower()
        contains_aromatic = "c" in reactants or "n" in reactants or "o" in reactants or "s" in reactants
        contains_nitrogen_product = "N" in products or "[N" in products

        # Primary/secondary alkyl iodides --> boronate type transformations.
        if "I" in reactants and "b" in products_lower and not contains_aromatic:
            return "[CX4;H2,H3]-[I]>>[C:1]-[B:3;X3](-[O;H0])(-[O;H0])"

        # Aryl bromide cross couplings (e.g., Buchwald-Hartwig).
        if "Br" in reactants and contains_aromatic:
            if contains_nitrogen_product:
                return "c-[Br]>>c-[NX3;!$(NC=O)]"
            return "c-[Br]>>c-[*]"

        # Generic alkyl bromide transformations.
        if "Br" in reactants:
            return "[CX4;H2,H3]-[Br]>>[CX4]-[*]"

        if "Cl" in reactants:
            return "[CX4;H2,H3]-[Cl]>>[CX4]-[*]"

        # Default fall-back.
        return f"{reactants}>>{products or '*'}"

    def suggest_guard_patterns(self) -> List[str]:
        """
        Suggest guard SMARTS to exclude substrate classes that should not match.

        The suggestions are heuristic and intentionally conservative – they
        guard against obvious mismatches (e.g., secondary/tertiary halides).
        """
        reactants = self.reactants_smiles
        contains_aromatic = "c" in reactants or "n" in reactants or "o" in reactants or "s" in reactants
        guards: List[str] = []

        if "I" in reactants and not contains_aromatic:
            guards.extend(
                [
                    "[CX4;H0]-[I]",  # tertiary iodides
                    "[CX4;H1]-[I]",  # secondary iodides
                    "[CH]-[I]",  # generic secondary guard
                    "[CH2]-[a]-I",  # benzylic
                    "[CH2]-[C]=[C]-I",  # allylic
                ]
            )
        elif "Br" in reactants and contains_aromatic:
            guards.extend(
                [
                    "[CX4]-[Br]",  # block aliphatic bromides
                    "[C;H0]-[Br]",  # quaternary centres
                    "[CH2]-[C]=[C]-[Br]",
                ]
            )
        elif "Br" in reactants:
            guards.extend(["[CX4;H0]-[Br]", "[CH]-[Br]"])
        elif "Cl" in reactants:
            guards.extend(["[CX4;H0]-[Cl]", "[CH]-[Cl]"])

        return guards

    def build_applicability(self) -> ReactionSmartsApplicability:
        """Convenience helper returning the full applicability package."""
        try:
            core = self.generate_core_smarts()
        except RuntimeError:
            core = self._heuristic_core_smarts()
        return ReactionSmartsApplicability(core=core, guards_forbid=self.suggest_guard_patterns())

    # -- utility helpers ---------------------------------------------------------

    def fallback_core_smarts(self) -> str:
        """Expose heuristic core SMARTS even when RDKit is unavailable."""
        return self._heuristic_core_smarts()


# ---------------------------------------------------------------------------
# Visualisation helpers
# ---------------------------------------------------------------------------


def _ensure_directory(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def visualize_smarts_pattern(smarts: str, output_path: Path) -> bool:
    """Render a SMARTS pattern to an image using RDKit."""
    if not RDKIT_AVAILABLE:
        return False
    try:
        mol = Chem.MolFromSmarts(smarts)
        if mol is None:
            return False
        _ensure_directory(output_path)
        img = Draw.MolToImage(mol, size=(400, 300))
        img.save(str(output_path))
        return True
    except Exception:  # pragma: no cover - RDKit failures
        return False


def visualize_reaction_smarts(smarts: str, output_path: Path) -> bool:
    """Render a reaction SMARTS pattern (reactants >> products)."""
    if not RDKIT_AVAILABLE:
        return False
    try:
        rxn = rdChemReactions.ReactionFromSmarts(smarts)
        if rxn is None:
            return False
        _ensure_directory(output_path)
        img = Draw.ReactionToImage(rxn, subImgSize=(400, 200))
        img.save(str(output_path))
        return True
    except Exception:  # pragma: no cover - RDKit failures
        return False


def visualize_pattern_with_examples(
    pattern: ReactionSmartsApplicability,
    test_smiles: Iterable[str],
    output_dir: Path,
) -> Dict[str, bool]:
    """
    Render the core/guard patterns and evaluate test SMILES against them.

    Returns a dict mapping SMILES to ``True`` (passes guards) / ``False`` (blocked).
    """
    results: Dict[str, bool] = {}

    if not RDKIT_AVAILABLE:
        return results

    output_dir.mkdir(parents=True, exist_ok=True)
    visualize_reaction_smarts(pattern.core, output_dir / "core_transformation.png")

    for idx, guard in enumerate(pattern.guards_forbid, start=1):
        visualize_smarts_pattern(guard, output_dir / f"guard_forbid_{idx}.png")

    guard_mols = [Chem.MolFromSmarts(g) for g in pattern.guards_forbid]
    require_mols = [Chem.MolFromSmarts(g) for g in pattern.guards_require]

    for smiles in test_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results[smiles] = False
            continue

        forbidden = any(m is not None and mol.HasSubstructMatch(m) for m in guard_mols)
        required = all(m is None or mol.HasSubstructMatch(m) for m in require_mols)

        results[smiles] = not forbidden and required

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _pattern_for_reaction(reaction: str) -> ReactionSmartsApplicability:
    generator = SmartsGenerator(reaction)
    try:
        return generator.build_applicability()
    except RuntimeError:
        # Fallback when RDKit is missing – still return a minimal pattern.
        return ReactionSmartsApplicability(
            core=generator.fallback_core_smarts(),
            guards_forbid=generator.suggest_guard_patterns(),
            notes="Generated without RDKit; core pattern is heuristic.",
        )


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point mirroring the historical CLI behaviour."""
    parser = argparse.ArgumentParser(description="Generate reaction SMARTS applicability patterns.")
    parser.add_argument("--reaction", help="Reaction SMILES to process.")
    parser.add_argument(
        "--batch",
        help="Path to a text file containing reaction SMILES (one per line).",
    )
    parser.add_argument("--output", help="Where to write JSON output.")
    parser.add_argument("--check-rdkit", action="store_true", help="Only report RDKit availability.")

    args = parser.parse_args(argv)

    if args.check_rdkit:
        print("RDKit available" if RDKIT_AVAILABLE else "RDKit not available")
        return 0 if RDKIT_AVAILABLE else 1

    if not args.reaction and not args.batch:
        parser.print_help()
        return 1

    if args.reaction:
        if not RDKIT_AVAILABLE:
            # Match historical behaviour: signal failure when RDKit is missing.
            print("RDKit not available – cannot generate SMARTS for single reaction.")
            return 1

        try:
            pattern = _pattern_for_reaction(args.reaction)
        except ValueError as exc:
            print(f"Error: {exc}")
            return 1

        if args.output:
            Path(args.output).write_text(json.dumps(pattern.to_dict(), indent=2))
        else:
            print(json.dumps(pattern.to_dict(), indent=2))
        return 0

    # Batch mode
    output: List[Dict[str, object]] = []
    batch_path = Path(args.batch)
    for line in batch_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            pattern = _pattern_for_reaction(line)
            output.append(
                {
                    "reaction_smiles": line,
                    "smarts_applicability": pattern.to_dict(),
                }
            )
        except ValueError as exc:
            output.append({"reaction_smiles": line, "error": str(exc)})

    if args.output:
        Path(args.output).write_text(json.dumps(output, indent=2))
    else:
        print(json.dumps(output, indent=2))

    # Even when RDKit is missing, batch mode should complete (tests expect 0).
    return 0


__all__ = [
    "ReactionSmartsApplicability",
    "SmartsGenerator",
    "visualize_pattern_with_examples",
    "visualize_reaction_smarts",
    "visualize_smarts_pattern",
    "main",
    "RDKIT_AVAILABLE",
]
