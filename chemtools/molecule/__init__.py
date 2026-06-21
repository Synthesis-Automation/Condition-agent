"""Molecule-level feature extraction and motif detection."""

from .featurize import featurize_molecule, print_motif_rank_breakdown
from .models import MoleculeFeatures, MoleculeParseResult
from .motifs import (
    CompoundPattern,
    build_compound_registry,
    calculate_smarts_complexity,
    choose_best_compound_hit,
    classify_compound_batch,
    classify_compound_smiles,
    detect_motifs,
    discover_undocumented_motifs,
)
from .parse import parse_molecule

__all__ = [
    "CompoundPattern",
    "MoleculeFeatures",
    "MoleculeParseResult",
    "build_compound_registry",
    "calculate_smarts_complexity",
    "choose_best_compound_hit",
    "classify_compound_batch",
    "classify_compound_smiles",
    "detect_motifs",
    "discover_undocumented_motifs",
    "featurize_molecule",
    "parse_molecule",
    "print_motif_rank_breakdown",
]
