"""
Retrosynthetic Disconnector — Core RDKit Engine.

Provides the algorithmic layer for:
  1. Retron detection   : match SMARTS retrons against a target molecule
  2. Transform application : generate precursor SMILES via retro-transforms
  3. Scoring            : rank disconnections by chemical quality
  4. Complexity scoring : BertzCT complexity delta (does this step simplify?)

Public API:
    from chemtools.retro.disconnector import rank_disconnections, find_retrons

Usage:
    results = rank_disconnections("c1ccc(-c2ccccc2)cc1")
    for r in results:
        print(r.description, r.precursor_1, r.precursor_2, r.difficulty)
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class RetronMatch:
    """A retron that matched in the target molecule."""
    retron_name: str
    reaction_name: str
    difficulty: float
    description: str
    notes: str
    precursor_hints: List[str]
    match_count: int               # How many times the SMARTS matched
    product_smarts: str


@dataclass
class DisconnectionResult:
    """A single retrosynthetic disconnection (one step back in a route)."""
    target_smiles: str
    precursor_1: str               # First precursor SMILES (electrophile synthon → reagent)
    precursor_2: str               # Second precursor SMILES (nucleophile synthon → reagent)
    reaction_name: str             # Taxonomy ID (e.g., "suzuki_miyaura")
    retron_name: str               # Which retron triggered this
    difficulty: float              # 0.0–1.0 (lower = easier)
    complexity_delta: float        # BertzCT(target) - max(BertzCT(p1), BertzCT(p2)); positive = simplified
    fragment_balance: float        # 0.0–1.0 (1.0 = perfectly equal halves)
    overall_score: float           # Combined ranking score (higher = better)
    description: str
    notes: str
    precursor_hints: List[str]     # Human-readable precursor type names


# ---------------------------------------------------------------------------
# Retron matching
# ---------------------------------------------------------------------------

def find_retrons(mol_or_smiles: Any, max_difficulty: float = 1.0) -> List[RetronMatch]:
    """
    Match all retron SMARTS patterns against a target molecule.

    Args:
        mol_or_smiles: RDKit Mol object or SMILES string of the target.
        max_difficulty: Filter out retrons harder than this threshold.

    Returns:
        List of RetronMatch objects sorted by difficulty ascending (easiest first).
    """
    try:
        from rdkit import Chem
    except ImportError:
        logger.error("RDKit is required for retron matching")
        return []

    # Resolve mol
    if isinstance(mol_or_smiles, str):
        mol = Chem.MolFromSmiles(mol_or_smiles)
        if mol is None:
            logger.warning(f"Could not parse SMILES for retron matching: {mol_or_smiles!r}")
            return []
    else:
        mol = mol_or_smiles

    from chemtools.util.smarts_cache import compile_smarts
    from .retron_patterns import RETRONS

    matches: List[RetronMatch] = []

    for retron in RETRONS:
        if retron["difficulty"] > max_difficulty:
            continue

        pattern = compile_smarts(retron["product_smarts"])
        if pattern is None:
            continue

        try:
            submatches = mol.GetSubstructMatches(pattern)
        except Exception:
            continue

        if submatches:
            matches.append(RetronMatch(
                retron_name=retron["name"],
                reaction_name=retron["reaction_name"],
                difficulty=retron["difficulty"],
                description=retron["description"],
                notes=retron["notes"],
                precursor_hints=retron.get("precursor_hints", []),
                match_count=len(submatches),
                product_smarts=retron["product_smarts"],
            ))

    # Sort: easiest first, then by number of matches (more common bonds = more useful)
    matches.sort(key=lambda m: (m.difficulty, -m.match_count))
    return matches


# ---------------------------------------------------------------------------
# Molecular complexity
# ---------------------------------------------------------------------------

def _bertz_complexity(mol: Any) -> float:
    """Compute BertzCT complexity. Returns 0.0 if RDKit unavailable."""
    try:
        from rdkit.Chem.GraphDescriptors import BertzCT
        return float(BertzCT(mol))
    except Exception:
        return 0.0


def _fragment_balance(mol1: Any, mol2: Any) -> float:
    """
    Compute fragment balance: 1.0 if perfectly equal halves, 0.0 if wildly imbalanced.
    Uses heavy atom count as a simple proxy.
    """
    try:
        n1 = mol1.GetNumHeavyAtoms() if mol1 else 1
        n2 = mol2.GetNumHeavyAtoms() if mol2 else 1
        smaller = min(n1, n2)
        larger = max(n1, n2)
        return smaller / larger if larger > 0 else 0.0
    except Exception:
        return 0.5


# ---------------------------------------------------------------------------
# Retrosynthetic transform application
# ---------------------------------------------------------------------------

def _apply_retro_transforms(
    mol: Any, retron: Dict[str, Any]
) -> List[Tuple[str, str]]:
    """
    Apply retrosynthetic transforms for a given retron to produce precursor pairs.

    Strategy:
      1. Find substructure matches in the target.
      2. For each matched bond (the bond that would be FORMED in forward reaction),
         fragment the molecule at that bond.
      3. Add appropriate leaving-group or functional-group "caps" to each fragment
         to produce realistic precursor SMILES.
      4. Validate precursor SMILES and return unique, canonicalized pairs.

    Returns:
        List of (precursor_1_smiles, precursor_2_smiles) tuples.
        Empty list if transform fails.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdmolops
    except ImportError:
        return []

    results: List[Tuple[str, str]] = []
    retron_name = retron["name"]

    # Use dedicated transform logic per retron class
    _transform_fn = _RETRON_TRANSFORMS.get(retron_name) or _RETRON_TRANSFORMS.get(
        retron.get("reaction_name", "")
    )

    if _transform_fn:
        try:
            pairs = _transform_fn(mol, retron)
            results.extend(pairs)
        except Exception as exc:
            logger.debug(f"Custom transform for {retron_name!r} failed: {exc}")

    # Generic fallback: bond fragmentation by retron match atoms
    if not results:
        pairs = _generic_bond_fragment(mol, retron)
        results.extend(pairs)

    # Deduplicate and validate — skip any pair where either SMILES is invalid
    from rdkit import rdBase
    seen: set = set()
    clean: List[Tuple[str, str]] = []
    for p1, p2 in results:
        if not p1 or not p2:
            continue
        try:
            with rdBase.BlockLogs():
                mol1 = Chem.MolFromSmiles(p1)
                mol2 = Chem.MolFromSmiles(p2)
            if mol1 is None or mol2 is None:
                continue
            canon_1 = Chem.MolToSmiles(mol1)
            canon_2 = Chem.MolToSmiles(mol2)
        except Exception:
            continue

        if canon_1 and canon_2 and (canon_1, canon_2) not in seen:
            seen.add((canon_1, canon_2))
            clean.append((canon_1, canon_2))

    return clean[:3]  # Return at most 3 distinct precursor pairs per retron


def _generic_bond_fragment(mol: Any, retron: Dict[str, Any]) -> List[Tuple[str, str]]:
    """
    Generic fragmentation: find the bond indicated by the retron SMARTS,
    break it, and return the two fragments. Each fragment gets a placeholder
    functional group based on the reaction class.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdmolops
        from chemtools.util.smarts_cache import compile_smarts
    except ImportError:
        return []

    pattern = compile_smarts(retron["product_smarts"])
    if pattern is None:
        return []

    try:
        matches = mol.GetSubstructMatches(pattern)
    except Exception:
        return []

    if not matches:
        return []

    results = []
    # Use first match (most common case)
    match = matches[0]
    if len(match) < 2:
        return []

    # Try to find a single bond connecting mapped atoms 1 and 2 in the match
    atom1_idx = match[0]
    atom2_idx = match[1]

    bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
    if bond is None:
        # Try adjacent atoms if direct bond not found
        return []

    bond_idx = bond.GetIdx()

    # Fragment at this bond
    try:
        frag_mol = rdmolops.FragmentOnBonds(
            mol, [bond_idx], addDummies=True, dummyLabels=[(0, 0)]
        )
        frags = rdmolops.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    except Exception as exc:
        logger.debug(f"Bond fragmentation failed for {retron['name']}: {exc}")
        return []

    if len(frags) != 2:
        return []

    # Convert fragments to SMILES, replacing dummy atoms ([*]) with functional groups
    reaction_class = _get_reaction_class(retron)
    caps = _get_dummy_caps(reaction_class)

    p1 = _replace_dummy(frags[0], caps[0])
    p2 = _replace_dummy(frags[1], caps[1])

    if p1 and p2:
        results.append((p1, p2))

    return results


def _replace_dummy(frag_mol: Any, cap: str) -> str:
    """
    Replace [*] dummy atoms in a fragment with a real functional group cap.
    Returns canonical SMILES or empty string on failure.
    """
    try:
        from rdkit import Chem
        from rdkit import rdBase

        # Use SMILES-level replacement: convert to SMILES, replace dummy notation
        smi = Chem.MolToSmiles(frag_mol)
        # RDKit dummy atom in SMILES is [*] or [1*] etc.
        import re
        smi_capped = re.sub(r"\[\d*\*\]|\[Xe\]", cap, smi)

        # Guard: a leading '=' (e.g. '=OC...') is never valid SMILES.
        # This happens when [*] was the first atom and cap starts with '='.
        if smi_capped.startswith("="):
            return ""

        # Validate — suppress RDKit's parse-error console noise
        with rdBase.BlockLogs():
            mol_check = Chem.MolFromSmiles(smi_capped)
        if mol_check is None:
            # Try without cap (leave as-is but remove dummy)
            smi_plain = re.sub(r"\[\d*\*\]|\[Xe\]", "", smi)
            with rdBase.BlockLogs():
                mol_check = Chem.MolFromSmiles(smi_plain)
            if mol_check:
                return Chem.MolToSmiles(mol_check)
            return ""

        return Chem.MolToSmiles(mol_check)
    except Exception:
        return ""


def _get_reaction_class(retron: Dict[str, Any]) -> str:
    """Derive a simple reaction class from retron metadata."""
    name = retron.get("name", "")
    rxn = retron.get("reaction_name", "")

    if "suzuki" in rxn or "negishi" in rxn or "biaryl" in name:
        return "c_c_coupling"
    if "buchwald" in rxn or "aryl_amine" in name:
        return "c_n_coupling"
    if "williamson" in rxn or "ether" in name:
        return "c_o_coupling"
    if "wittig" in rxn or "heck" in rxn:
        return "alkene_forming"
    if "aldol" in rxn:
        return "c_c_addition"
    if "reduction" in rxn:
        return "reduction"
    if "oxidation" in rxn or "oxidation" in name:
        return "oxidation"
    return "generic"


def _get_dummy_caps(reaction_class: str) -> Tuple[str, str]:
    """Return (cap1, cap2) SMILES atoms/groups for each fragment dummy position."""
    caps = {
        # aryl-Br + aryl-B(OH)2  (Suzuki)
        # OB(O) cap: B bonds directly to ring C. B(O)O would put O between B and ring.
        "c_c_coupling": ("Br", "OB(O)"),
        # aryl-Br + amine: amine fragment dummy removed → implicit H on N remains
        "c_n_coupling": ("Br", ""),
        # alkyl-Br + alcohol: alcohol dummy removed
        "c_o_coupling": ("Br", ""),
        # Wittig: each fragment loses one C of the alkene; each becomes a carbonyl
        "alkene_forming": ("C=O", "C=O"),
        # aldol: one fragment gets a carbonyl appendage
        "c_c_addition": ("C=O", ""),
        # parent carbonyl (for alcohol retrons)
        "reduction": ("C=O", ""),
        # parent alcohol (for oxidation retrons)
        "oxidation": ("O", ""),
        "generic": ("", ""),
    }
    return caps.get(reaction_class, ("", ""))


# ---------------------------------------------------------------------------
# Dedicated transform functions (higher accuracy for key reaction classes)
# ---------------------------------------------------------------------------

def _transform_biaryl_suzuki(mol: Any, retron: Dict[str, Any]) -> List[Tuple[str, str]]:
    """For biaryl bonds: precursors are aryl bromide + phenylboronic acid."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdmolops
        from chemtools.util.smarts_cache import compile_smarts
    except ImportError:
        return []

    pattern = compile_smarts("[c:1]-[c:2]")
    if pattern is None:
        return []

    matches = mol.GetSubstructMatches(pattern)
    results = []

    for match in matches[:2]:  # Try up to 2 biaryl bonds
        atom1_idx, atom2_idx = match[0], match[1]
        bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
        if bond is None:
            continue

        bond_idx = bond.GetIdx()
        try:
            frag_mol = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=True)
            frags = rdmolops.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
        except Exception:
            continue

        if len(frags) != 2:
            continue

        # Fragment 1 → aryl halide (replace dummy with Br)
        p1 = _replace_dummy(frags[0], "Br")
        # Fragment 2 → aryl boronic acid (replace dummy with OB(O), NOT B(O)O)
        # OB(O) keeps B as the chain atom that bonds directly to the ring carbon.
        # B(O)O would insert O between B and the ring, producing a borate ester
        # (B-O-Ar) instead of a boronic acid (Ar-B(OH)2). This gives OBOc1ccoc1
        # (B has no C neighbor) vs OB(O)c1ccoc1 (B directly bonded to C). ✓
        p2 = _replace_dummy(frags[1], "OB(O)")

        if p1 and p2:
            results.append((p1, p2))

    return results


def _transform_aryl_amine_buchwald(mol: Any, retron: Dict[str, Any]) -> List[Tuple[str, str]]:
    """For aryl amine bonds: aryl bromide + NH amine."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdmolops
        from chemtools.util.smarts_cache import compile_smarts
    except ImportError:
        return []

    # Match Ar-N bond
    pattern = compile_smarts("[c:1]-[NX3:2]")
    if pattern is None:
        return []

    matches = mol.GetSubstructMatches(pattern)
    results = []

    def _fragment_role(frag: Any) -> str:
        """Classify the cleavage-site fragment as aryl-electrophile or amine partner."""
        for atom in frag.GetAtoms():
            if atom.GetAtomicNum() != 0:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 7:
                    return "amine"
                if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
                    return "aryl"
        return "other"

    for match in matches[:2]:
        atom_c, atom_n = match[0], match[1]
        bond = mol.GetBondBetweenAtoms(atom_c, atom_n)
        if bond is None:
            continue

        bond_idx = bond.GetIdx()
        try:
            frag_mol = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=True)
            frags = rdmolops.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
        except Exception:
            continue

        if len(frags) != 2:
            continue

        role0 = _fragment_role(frags[0])
        role1 = _fragment_role(frags[1])

        if role0 == "aryl" and role1 == "amine":
            aryl_frag, amine_frag = frags[0], frags[1]
        elif role0 == "amine" and role1 == "aryl":
            aryl_frag, amine_frag = frags[1], frags[0]
        else:
            # Fallback for ambiguous fragments: prefer the fragment whose
            # cleavage atom is an aromatic carbon as the electrophile.
            aryl_frag, amine_frag = frags[0], frags[1]

        p_arene = _replace_dummy(aryl_frag, "Br")
        p_amine = _replace_dummy(amine_frag, "")  # remove dummy; N retains implicit H

        if not p_arene or not p_amine:
            continue

        p_amine_mol = Chem.MolFromSmiles(p_amine)
        if p_amine_mol and any(a.GetAtomicNum() == 7 for a in p_amine_mol.GetAtoms()):
            results.append((p_arene, p_amine))

    return results


def _transform_alkene_wittig(mol: Any, retron: Dict[str, Any]) -> List[Tuple[str, str]]:
    """For C=C bonds: aldehyde/ketone fragment + phosphonium ylide."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdmolops
        from chemtools.util.smarts_cache import compile_smarts
    except ImportError:
        return []

    pattern = compile_smarts("[CX3:1]=[CX3:2]")
    if pattern is None:
        return []

    matches = mol.GetSubstructMatches(pattern)
    results = []

    for match in matches[:1]:
        atom1, atom2 = match[0], match[1]
        bond = mol.GetBondBetweenAtoms(atom1, atom2)
        if bond is None:
            continue

        bond_idx = bond.GetIdx()
        try:
            frag_mol = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=True)
            frags = rdmolops.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
        except Exception:
            continue

        if len(frags) != 2:
            continue

        # Each fragment becomes a carbonyl precursor (replace dummy with =O)
        p1 = _replace_dummy(frags[0], "=O")
        p2 = _replace_dummy(frags[1], "=O")

        if p1 and p2:
            # Describe as aldehyde/ketone + aldehyde/ketone (HWE precursors)
            results.append((p1, p2))

    return results


def _transform_amide(mol: Any, retron: Dict[str, Any]) -> List[Tuple[str, str]]:
    """For amide bonds: carboxylic acid + amine."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdmolops
        from chemtools.util.smarts_cache import compile_smarts
    except ImportError:
        return []

    # C(=O)N bond
    pattern = compile_smarts("[CX3:1](=[O])-[NX3:2]")
    if pattern is None:
        return []

    matches = mol.GetSubstructMatches(pattern)
    results = []

    for match in matches[:2]:
        atom_c = match[0]
        atom_n = match[1]
        bond = mol.GetBondBetweenAtoms(atom_c, atom_n)
        if bond is None:
            continue

        bond_idx = bond.GetIdx()
        try:
            frag_mol = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=True)
            frags = rdmolops.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
        except Exception:
            continue

        if len(frags) != 2:
            continue

        # Carbonyl fragment → carboxylic acid (replace dummy with O for OH)
        p_acid = _replace_dummy(frags[0], "O")
        # Amine fragment → remove dummy; N retains implicit H
        p_amine = _replace_dummy(frags[1], "")

        if p_acid and p_amine:
            results.append((p_acid, p_amine))

    return results


# Registry of dedicated transform functions
_RETRON_TRANSFORMS: Dict[str, Any] = {
    "biaryl_suzuki": _transform_biaryl_suzuki,
    "suzuki_miyaura": _transform_biaryl_suzuki,
    "aryl_amine_buchwald": _transform_aryl_amine_buchwald,
    "buchwald_hartwig_amination": _transform_aryl_amine_buchwald,
    "alkene_wittig": _transform_alkene_wittig,
    "wittig_reaction": _transform_alkene_wittig,
    "amide_direct": _transform_amide,
    "amide_coupling": _transform_amide,
}


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def _score_disconnection(
    target_mol: Any,
    prec1_smiles: str,
    prec2_smiles: str,
    retron_match: RetronMatch,
) -> Tuple[float, float, float]:
    """
    Score a disconnection on three axes:
      - complexity_delta: BertzCT(target) - max(BertzCT(p1), BertzCT(p2))
        Positive = simplification. Normalized to 0–1 by dividing by target complexity.
      - fragment_balance: min(n_atoms) / max(n_atoms), 0–1.
      - overall_score: weighted combination (higher = better).

    Returns:
        (complexity_delta, fragment_balance, overall_score)
    """
    try:
        from rdkit import Chem
    except ImportError:
        return (0.0, 0.5, 0.5 - retron_match.difficulty)

    target_ct = _bertz_complexity(target_mol)

    p1_mol = Chem.MolFromSmiles(prec1_smiles) if prec1_smiles else None
    p2_mol = Chem.MolFromSmiles(prec2_smiles) if prec2_smiles else None

    p1_ct = _bertz_complexity(p1_mol) if p1_mol else 0.0
    p2_ct = _bertz_complexity(p2_mol) if p2_mol else 0.0
    max_prec_ct = max(p1_ct, p2_ct)

    # Complexity delta: positive = step reduces complexity (good)
    complexity_delta = target_ct - max_prec_ct

    # Normalize to -1..+1 scale
    norm_delta = complexity_delta / max(target_ct, 1.0)

    # Fragment balance
    balance = _fragment_balance(p1_mol, p2_mol)

    # Overall score: favour easy, complex-reducing, balanced disconnections
    overall = (
        0.40 * (1.0 - retron_match.difficulty)      # ease of reaction
        + 0.35 * max(0.0, min(1.0, norm_delta))      # complexity reduction
        + 0.25 * balance                              # convergence
    )

    return (complexity_delta, balance, overall)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def rank_disconnections(
    target_smiles: str,
    top_k: int = 5,
    max_difficulty: float = 0.8,
) -> List[DisconnectionResult]:
    """
    Main retrosynthesis function: find and rank retrosynthetic disconnections.

    For the given target molecule, this function:
      1. Parses and validates the SMILES.
      2. Matches all applicable retron SMARTS patterns.
      3. Applies retrosynthetic transforms to generate precursor pairs.
      4. Scores each disconnection (complexity delta + balance + difficulty).
      5. Returns the top-k disconnections ranked by overall score.

    Args:
        target_smiles: SMILES of the target molecule.
        top_k: Maximum number of disconnections to return.
        max_difficulty: Filter out retrons harder than this (0=trivial, 1=heroic).

    Returns:
        List of DisconnectionResult, sorted best-first (highest overall_score).
        Returns empty list if molecule cannot be parsed.
    """
    try:
        from rdkit import Chem
    except ImportError:
        logger.error("RDKit required for rank_disconnections")
        return []

    mol = Chem.MolFromSmiles(target_smiles)
    if mol is None:
        logger.warning(f"Cannot parse SMILES for disconnection: {target_smiles!r}")
        return []

    # Step 1: Find matching retrons
    retron_matches = find_retrons(mol, max_difficulty=max_difficulty)
    if not retron_matches:
        logger.info(f"No retrons matched for: {target_smiles!r}")
        return []

    from .retron_patterns import get_retron_by_name

    concrete_results: List[DisconnectionResult] = []
    conceptual_results: List[DisconnectionResult] = []

    for match in retron_matches:
        retron_def = get_retron_by_name(match.retron_name)
        if not retron_def:
            continue

        # Step 2: Apply retrosynthetic transforms
        precursor_pairs = _apply_retro_transforms(mol, retron_def)

        if not precursor_pairs:
            # Keep conceptual matches only as a low-priority fallback so they do
            # not outrank concrete precursor pairs elsewhere in the route.
            conceptual_results.append(DisconnectionResult(
                target_smiles=target_smiles,
                precursor_1="",
                precursor_2="",
                reaction_name=match.reaction_name,
                retron_name=match.retron_name,
                difficulty=match.difficulty,
                complexity_delta=0.0,
                fragment_balance=0.0,
                overall_score=0.0,
                description=match.description,
                notes=match.notes,
                precursor_hints=match.precursor_hints,
            ))
            continue

        for p1_smi, p2_smi in precursor_pairs:
            # Step 3: Score
            complexity_delta, balance, overall = _score_disconnection(
                mol, p1_smi, p2_smi, match
            )

            concrete_results.append(DisconnectionResult(
                target_smiles=target_smiles,
                precursor_1=p1_smi,
                precursor_2=p2_smi,
                reaction_name=match.reaction_name,
                retron_name=match.retron_name,
                difficulty=match.difficulty,
                complexity_delta=complexity_delta,
                fragment_balance=balance,
                overall_score=overall,
                description=match.description,
                notes=match.notes,
                precursor_hints=match.precursor_hints,
            ))

    # Sort concrete results by chemical score; conceptual-only matches stay last.
    concrete_results.sort(key=lambda r: r.overall_score, reverse=True)
    conceptual_results.sort(key=lambda r: r.difficulty)
    results = concrete_results + conceptual_results

    # Deduplicate by reaction_name (keep top score per reaction class)
    seen_rxn: set = set()
    deduped: List[DisconnectionResult] = []
    for r in results:
        if r.reaction_name not in seen_rxn:
            seen_rxn.add(r.reaction_name)
            deduped.append(r)
        elif len(deduped) < top_k:
            deduped.append(r)

    return deduped[:top_k]
