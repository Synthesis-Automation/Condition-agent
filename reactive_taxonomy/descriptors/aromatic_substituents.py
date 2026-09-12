"""Position-aware, explicitly heuristic aromatic substituent contributions."""

from __future__ import annotations

from typing import Any, Tuple

from ..chemistry.smarts_cache import compile_smarts
from .models import ElectronicContribution
from .registry import descriptor_rules


def aromatic_substituent_contributions(
    mol: Any, center: int,
) -> Tuple[ElectronicContribution, ...]:
    """Separate induction and resonance on the same six-membered aromatic ring.

    Graph attachment and ring position are observations; contribution magnitudes
    are versioned chemistry priors, not measured electronic properties. Other
    ring sizes and remote fused-ring pathways remain outside this rule set.
    """
    rings = tuple(
        tuple(ring) for ring in mol.GetRingInfo().AtomRings()
        if len(ring) == 6 and center in ring
        and all(mol.GetAtomWithIdx(index).GetIsAromatic() for index in ring)
    )
    rules = descriptor_rules()["electronic"]["aromatic_substituents"]
    values = []
    seen = set()
    for rule in rules["rules"]:
        query = compile_smarts(rule["smarts"], validate=True)
        positions = {atom.GetAtomMapNum(): atom.GetIdx() for atom in query.GetAtoms()}
        for match in mol.GetSubstructMatches(query):
            attachment, ipso = match[positions[1]], match[positions[2]]
            distances = []
            for ring in rings:
                if ipso not in ring or ipso == center:
                    continue
                delta = abs(ring.index(center) - ring.index(ipso))
                distances.append(min(delta, len(ring) - delta))
            if not distances:
                continue
            relation = rules["relations"][str(min(distances))]
            key = (rule["id"], attachment, ipso)
            if key in seen:
                continue
            seen.add(key)
            for pathway in ("inductive", "resonance"):
                contribution = float(rule[pathway][relation])
                if not contribution:
                    continue
                values.append(ElectronicContribution(
                    source_id=f"aromatic_substituent:{rule['id']}",
                    effect="withdrawing" if contribution > 0 else "donating",
                    pathway=pathway,
                    positional_relation=relation,
                    contribution=contribution,
                    atom_indices=tuple(sorted(match)),
                ))
    return tuple(sorted(values, key=lambda item: (
        item.source_id, item.pathway, item.positional_relation, item.atom_indices,
    )))
