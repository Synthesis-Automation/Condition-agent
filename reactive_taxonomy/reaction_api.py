"""Public reaction featurization API."""

from __future__ import annotations

from typing import Any, Dict, List

from .chemistry.rdkit_utils import parse_smiles

from .labels import available_styles
from .reaction_bond_changes import supplied_map_bond_changes
from .reaction_candidates import enumerate_reaction_candidates
from .reaction_labels import render_reactant_label, render_reaction_label
from .reaction_environments import build_reaction_family_environment
from .reaction_models import ReactionAnalysis, ReactionCandidate
from .reaction_operators import apply_operator
from .reaction_parser import parse_reaction_smiles
from .reaction_spectators import derive_spectator_groups


def _canonical_without_maps(smiles: str) -> str | None:
    from rdkit import Chem
    mol = parse_smiles(smiles)
    if mol is None: return None
    for atom in mol.GetAtoms(): atom.SetAtomMapNum(0)
    try: return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception: return None


def featurize_reaction(
    reaction_smiles: str,
    *,
    label_style: str = "unicode",
    max_candidates: int = 500,
) -> ReactionAnalysis:
    """Analyze a reaction with site grammars and product reconstruction."""
    parsed = parse_reaction_smiles(reaction_smiles)
    if label_style not in available_styles():
        return ReactionAnalysis(reaction_smiles, False, error=f"UNKNOWN_LABEL_STYLE:{label_style}")
    if not parsed.valid:
        return ReactionAnalysis(reaction_smiles, False, parsed.reactants, parsed.agents, parsed.products, warnings=parsed.warnings, error=parsed.error)
    observed_products = {_canonical_without_maps(component.input_smiles) for component in parsed.products}
    observed_products.discard(None)
    raw = enumerate_reaction_candidates(parsed.reactants)
    warnings = list(parsed.warnings)
    if len(raw) > max_candidates:
        raw = raw[:max_candidates]
        warnings.append("CANDIDATE_LIMIT_REACHED")
    candidates: List[ReactionCandidate] = []
    exact: List[ReactionCandidate] = []
    for grammar, assignment in raw:
        predicted, changes = apply_operator(grammar, assignment, parsed.reactants)
        predicted_canonical = _canonical_without_maps(predicted) if predicted else None
        verification = "exact_product_reconstruction" if predicted_canonical in observed_products else ("product_mismatch" if predicted else "construction_failed")
        label = render_reaction_label(grammar, assignment, style=label_style)
        candidate = ReactionCandidate(
            grammar_id=grammar["id"], transformation_class=grammar["transformation_class"],
            role_assignments=assignment, predicted_bond_changes=changes,
            predicted_product_smiles=predicted_canonical, verification=verification,
            reaction_label=label,
            compatible_named_families=tuple(grammar.get("compatible_named_families") or []),
        )
        candidates.append(candidate)
        if verification == "exact_product_reconstruction": exact.append(candidate)
    selected = None
    evidence = "reactant_grammar_only" if candidates else "unresolved"
    if len(exact) == 1:
        selected, evidence = exact[0], "exact_product_reconstruction"
    elif len(exact) > 1:
        signatures = {
            (candidate.grammar_id, tuple(sorted(site.canonical_signature for site in candidate.role_assignments.values())))
            for candidate in exact
        }
        if len(signatures) == 1:
            selected, evidence = exact[0], "exact_product_reconstruction"
            warnings.append("SYMMETRY_EQUIVALENT_ASSIGNMENTS")
        else:
            evidence = "ambiguous"
            warnings.append("AMBIGUOUS_PARTICIPATING_SITES")
    mapped_changes = tuple(supplied_map_bond_changes(reaction_smiles))
    spectators = derive_spectator_groups(parsed.reactants, selected, evidence)
    family_environment = build_reaction_family_environment(parsed.reactants, selected, spectators, evidence)
    reaction_label = selected.reaction_label if selected else None
    reaction_label_status = "exact_product" if selected else "unavailable"
    if selected is None and candidates:
        reactant_labels = sorted({
            render_reactant_label(assignment, style=label_style)
            for _, assignment in raw
        })
        if len(reactant_labels) == 1:
            reaction_label = f"{reactant_labels[0]} →"
            reaction_label_status = "reactant_only"
        elif reactant_labels:
            reaction_label = " OR ".join(f"({label})" for label in reactant_labels) + " →"
            reaction_label_status = "ambiguous_reactants"
    return ReactionAnalysis(
        input_reaction_smiles=reaction_smiles, valid=True,
        reactants=parsed.reactants, agents=parsed.agents, products=parsed.products,
        candidates=tuple(candidates), selected_candidate=selected,
        transformation_class=selected.transformation_class if selected else None,
        compatible_named_families=selected.compatible_named_families if selected else (),
        named_family=selected.compatible_named_families[0] if selected and len(selected.compatible_named_families) == 1 else None,
        reaction_label=reaction_label,
        reaction_label_status=reaction_label_status,
        evidence_quality=evidence, mapped_bond_changes=mapped_changes,
        spectator_groups=spectators,
        family_environment=family_environment,
        warnings=tuple(warnings),
    )


__all__ = ["featurize_reaction"]
