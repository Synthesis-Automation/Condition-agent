"""Recover documented atom correspondence from identical indexed reactions."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.chemistry.reaction_smiles import split_reaction_fields
from reactive_taxonomy.chemistry.source_mapping import assess_source_mapping, source_material_graphs
from reactive_taxonomy.reaction_parser import parse_reaction_smiles
from reactive_taxonomy.shared_reaction_core import whole_reaction_graph_key

from .generic_indexing import GenericReactionIndex
from .shared_core_index import SharedCoreIndex


@dataclass(frozen=True)
class MappingSource:
    """An indexed observation inspected for graph identity and atom correspondence."""

    reaction_id: str
    observation_id: str
    reference_id: str
    mapping_equivalence_key: str | None


@dataclass(frozen=True)
class ExactReactionMapping:
    """A bounded lookup of recorded atom maps, independent of conditions and yield."""

    status: str
    index_source_identity: str
    candidates_examined: int = 0
    mapped_reaction_smiles: str | None = None
    sources: tuple[MappingSource, ...] = ()
    schema_version: str = "exact_reaction_mapping.v1"

    def to_dict(self) -> dict[str, Any]:
        """Serialize the outcome and observation-level evidence."""
        return asdict(self)


def _checked_mapping_key(original: str, mapped: str) -> str | None:
    """Require agreement of all graphs, observed product origins and verified edits."""
    if assess_source_mapping(original, mapped).issues:
        return None
    analysis = featurize_reaction(mapped)
    core = analysis.reaction_core
    if (not core or core.evidence_status != "verified" or core.quality.status == "blocked"
            or not analysis.reaction_completeness or analysis.reaction_completeness.status != "verified"):
        return None
    return core.mapping_equivalence_key


def recover_exact_reaction_mapping(
    reaction_smiles: str, index: GenericReactionIndex, shared_index: SharedCoreIndex, *,
    candidate_limit: int = 128,
) -> ExactReactionMapping:
    """Find a validated source map for an unresolved, unmapped reaction.

    Both full sides, stereo and component multiplicity must agree. Query agents
    remain as submitted. Conflicting edit graphs or a capped exact lookup leave
    correspondence unresolved. Only this structural evidence is transferred.
    """
    if type(candidate_limit) is not int or not 1 <= candidate_limit <= 512:
        raise ValueError("candidate_limit must be an integer between 1 and 512")
    shared_index.validate_binding(index)
    result = ExactReactionMapping("not_found", shared_index.source_identity)
    parsed = parse_reaction_smiles(reaction_smiles, include_molecular_interpretation=False)
    if not parsed.valid or any(not c.molecule_analysis.structure.valid for c in parsed.agents):
        return replace(result, status="invalid_query")
    if any(c.atom_mapped for c in (*parsed.reactants, *parsed.products, *parsed.agents)):
        return replace(result, status="supplied_mapping_preserved")
    left, agents, right = split_reaction_fields(reaction_smiles)
    inputs, outputs = source_material_graphs(left), source_material_graphs(right)
    key = whole_reaction_graph_key(inputs, outputs)
    positions, truncated = shared_index.lookup("whole_reaction", key, candidate_limit)
    checked: dict[str, str | None] = {}
    sources = []
    for row in index.select(positions):
        before, _, after = split_reaction_fields(row.reaction_smiles)
        if source_material_graphs(before) != inputs or source_material_graphs(after) != outputs:
            continue
        candidate = f"{before}>{agents}>{after}"
        # A capped lookup cannot establish agreement. Preserve the graph-identical
        # source IDs but skip costly validation and leave their edit keys unset.
        if not truncated and candidate not in checked:
            checked[candidate] = _checked_mapping_key(reaction_smiles, candidate)
        sources.append(MappingSource(row.reaction_id, row.observation_id, row.reference_id, checked.get(candidate)))
    result = replace(result, candidates_examined=len(positions), sources=tuple(sources))
    if truncated:
        return replace(result, status="truncated")
    keys = set(checked.values()) - {None}
    if len(keys) > 1:
        return replace(result, status="conflicting_mappings")
    if not keys or None in checked.values():
        return replace(result, status="no_validated_mapping" if sources else "not_found")
    selected = next(candidate for candidate, mapping_key in checked.items() if mapping_key)
    return replace(result, status="recovered", mapped_reaction_smiles=selected)
