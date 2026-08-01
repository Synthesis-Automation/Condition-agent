"""Reaction-SMILES parsing without automatic atom mapping."""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from typing import List, Tuple

from .api import interpret_molecular_reactivity, observe_molecular_structure
from .models import MoleculeAnalysis, MolecularInterpretation
from .reaction_models import ReactionComponent


_MAP_RE = re.compile(r":\d+\]")


@dataclass(frozen=True)
class ParsedReaction:
    valid: bool
    reactants: Tuple[ReactionComponent, ...] = ()
    agents: Tuple[ReactionComponent, ...] = ()
    products: Tuple[ReactionComponent, ...] = ()
    warnings: Tuple[str, ...] = ()
    error: str | None = None


def _components(text: str, side: str) -> Tuple[List[ReactionComponent], List[str]]:
    components: List[ReactionComponent] = []
    warnings: List[str] = []
    for index, smiles in enumerate(token for token in text.split(".") if token):
        structure = observe_molecular_structure(smiles)
        analysis = MoleculeAnalysis(
            structure=structure,
            interpretation=MolecularInterpretation(),
        )
        if not analysis.structure.valid:
            warnings.append(f"INVALID_{side.upper()}_COMPONENT:{index}")
        components.append(ReactionComponent(
            side=side, component_index=index, input_smiles=smiles,
            canonical_smiles=analysis.structure.canonical_smiles or "",
            atom_mapped=bool(_MAP_RE.search(smiles)), molecule_analysis=analysis,
        ))
    return components, warnings


def parse_reaction_smiles(
    reaction_smiles: str,
    *,
    include_molecular_interpretation: bool = True,
) -> ParsedReaction:
    text = str(reaction_smiles or "").strip()
    if not text:
        return ParsedReaction(False, error="EMPTY_REACTION_SMILES")
    if ">>" in text:
        left, right = text.split(">>", 1)
        middle = ""
    else:
        parts = text.split(">")
        if len(parts) != 3:
            return ParsedReaction(False, error="INVALID_REACTION_FORMAT")
        left, middle, right = parts
    if not left.strip(): return ParsedReaction(False, error="EMPTY_REACTANT_SIDE")
    if not right.strip(): return ParsedReaction(False, error="EMPTY_PRODUCT_SIDE")
    reactants, rw = _components(left, "reactant")
    agents, aw = _components(middle, "agent")
    products, pw = _components(right, "product")
    warnings = rw + aw + pw
    if len(products) > 1: warnings.append("MULTIPLE_PRODUCTS")
    valid = all(
        component.molecule_analysis.structure.valid
        for component in (*reactants, *products)
    )
    parsed = ParsedReaction(valid, tuple(reactants), tuple(agents), tuple(products), tuple(warnings), None if valid else "INVALID_COMPONENT")
    return (
        interpret_parsed_molecules(parsed)
        if include_molecular_interpretation and parsed.valid
        else parsed
    )


def interpret_parsed_molecules(parsed: ParsedReaction) -> ParsedReaction:
    """Attach optional molecular annotations after structural observation."""

    def annotate(component: ReactionComponent) -> ReactionComponent:
        structure = component.molecule_analysis.structure
        return replace(
            component,
            molecule_analysis=MoleculeAnalysis(
                structure=structure,
                interpretation=interpret_molecular_reactivity(structure),
            ),
        )

    return replace(
        parsed,
        reactants=tuple(annotate(component) for component in parsed.reactants),
        agents=tuple(annotate(component) for component in parsed.agents),
        products=tuple(annotate(component) for component in parsed.products),
    )


__all__ = [
    "ParsedReaction",
    "interpret_parsed_molecules",
    "parse_reaction_smiles",
]
