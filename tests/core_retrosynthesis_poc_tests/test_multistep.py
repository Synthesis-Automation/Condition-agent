"""Bounded deterministic multistep retrosynthesis regressions."""

from __future__ import annotations

import csv
from dataclasses import replace

from core_retrosynthesis_poc import (
    GenericDisconnectionCandidate,
    plan_multistep_routes,
    route_tree_distance,
)
from cas_tools import (
    CanonicalMoleculeIndex,
    MoleculeIndexMatch,
    build_canonical_molecule_index,
)


class _LiteratureIndex:
    def __init__(
        self,
        known: tuple[str, ...] = (),
        *,
        source_role: str | None = "reactant",
    ) -> None:
        from cas_tools import molecule_identity

        self.known = {molecule_identity(value).canonical_smiles for value in known}
        self.source_role = source_role

    def lookup(self, identity, *, provenance_limit: int = 5):
        if identity.canonical_smiles not in self.known:
            return None
        source_record = {"reaction_id": "known"}
        if self.source_role is not None:
            source_record["source_role"] = self.source_role
        return MoleculeIndexMatch(
            canonical_smiles=identity.canonical_smiles,
            inchi_key=identity.inchi_key,
            occurrence_count=1,
            source_records=(source_record,),
        )


def _candidate(product: str, precursors: str, score: float = 0.9):
    return GenericDisconnectionCandidate(
        target_smiles=product,
        precursor_smiles=precursors,
        proposed_reaction_smiles=f"{precursors}>>{product}",
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="reaction_core",
        template_id=f"template:{product}:{precursors}",
        score=score,
        context_similarity=0.0,
        product_similarity=score,
        precursor_similarity=score,
        template_specificity=1.0,
        independent_reference_support=1,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key="site",
        precedent_reaction_ids=("reaction",),
        operator_id="operator",
        realization_id="realization",
        operator_signature="signature",
        synthon_signature="synthon",
    )


def _expander(mapping):
    def expand(product: str, top_k: int):
        return tuple(mapping.get(product, ()))[:top_k]

    return expand


def test_two_depth_route_requires_every_leaf_to_be_terminal() -> None:
    expansions = {
        "CCCCCCCC": (_candidate("CCCCCCCC", "C.CCCCCC"),),
        "CCCCCC": (_candidate("CCCCCC", "CCN.CCO"),),
    }

    result = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=50.0,
        expander=_expander(expansions),
    )

    assert len(result.routes) == 1
    route = result.routes[0]
    assert route.solved is True
    assert route.reaction_count == 2
    assert route.maximum_depth == 2
    assert {leaf.canonical_smiles for leaf in route.leaves} == {
        "C",
        "CCN",
        "CCO",
    }
    assert all(leaf.terminal for leaf in route.leaves)
    assert route.route_tree.root.smiles == "CCCCCCCC"
    assert route.route_tree.reaction_count == route.reaction_count
    assert route.route_tree.root.reaction is not None
    assert len(route.route_tree.root.reaction.children) == 2
    expanded_child = next(
        child
        for child in route.route_tree.root.reaction.children
        if child.smiles == "CCCCCC"
    )
    assert expanded_child.reaction is not None
    assert route_tree_distance(route.route_tree, route.route_tree) == 0.0
    assert result.diagnostics.proposed_actions == 2
    assert result.diagnostics.validation_attempts == 2
    assert result.diagnostics.valid_actions == 2
    assert result.diagnostics.first_solution_expansion == 2

    serialized = result.to_dict()
    assert serialized["routes"][0]["route_tree"]["tree_id"]
    assert serialized["route_postprocessing"]["distance_matrix"] == [[0.0]]


def test_root_never_terminates_without_one_disconnection() -> None:
    result = plan_multistep_routes(
        "CC",
        object(),
        _LiteratureIndex(("CC",)),
        max_depth=2,
        molecular_weight_threshold=150.0,
        expander=_expander({"CC": (_candidate("CC", "C.O"),)}),
    )

    assert result.routes[0].reaction_count == 1
    assert result.routes[0].solved is True


def test_high_weight_literature_match_is_terminal() -> None:
    known = "CCCCCCCCCC"
    result = plan_multistep_routes(
        "CCCCCCCCCCCC",
        object(),
        _LiteratureIndex((known,)),
        max_depth=2,
        molecular_weight_threshold=20.0,
        expander=_expander(
            {"CCCCCCCCCCCC": (_candidate("CCCCCCCCCCCC", f"C.{known}"),)}
        ),
    )

    matched = next(
        leaf for leaf in result.routes[0].leaves if leaf.canonical_smiles == known
    )
    assert matched.molecular_weight > 20.0
    assert matched.terminal_reasons == ("literature_reactant_match",)
    assert matched.terminal_evidence == "role_supported_literature"
    assert matched.catalog_role_status == "reactant_supported"
    assert matched.literature_match is not None


def test_role_provenance_round_trips_through_canonical_index(tmp_path) -> None:
    known = "CCCCCCCCCC"
    source = tmp_path / "starting-materials.csv"
    output = tmp_path / "starting-materials.sqlite"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("compound_smiles", "reaction_id", "source_role"),
        )
        writer.writeheader()
        writer.writerow(
            {
                "compound_smiles": known,
                "reaction_id": "RXN-1",
                "source_role": "reactant",
            }
        )
    build_canonical_molecule_index(
        source,
        output,
        provenance_columns=("reaction_id", "source_role"),
    )

    with CanonicalMoleculeIndex(output) as index:
        result = plan_multistep_routes(
            "CCCCCCCCCCCC",
            object(),
            index,
            max_depth=2,
            molecular_weight_threshold=20.0,
            expander=_expander(
                {"CCCCCCCCCCCC": (_candidate("CCCCCCCCCCCC", f"C.{known}"),)}
            ),
        )

    assert result.routes
    matched = next(
        leaf for leaf in result.routes[0].leaves if leaf.canonical_smiles == known
    )
    assert matched.catalog_role_status == "reactant_supported"


def test_depth_limit_returns_partial_route_and_cycle_is_rejected() -> None:
    expansions = {
        "CCCCCCCC": (_candidate("CCCCCCCC", "CCCCCC"),),
        "CCCCCC": (
            _candidate("CCCCCC", "CCCCCCCC", score=0.95),
            _candidate("CCCCCC", "CCCCC", score=0.90),
        ),
    }

    result = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=20.0,
        expander=_expander(expansions),
    )

    assert not result.routes
    assert result.partial_routes
    assert result.partial_routes[0].leaves[0].unresolved_reason == "maximum_depth"
    assert result.diagnostics.rejected_cycles == 1


def test_nonreactant_and_untyped_catalog_matches_are_not_strong_terminals() -> None:
    known = "CCCCCCCCCC"
    expansion = _expander({"CCCCCCCCCCCC": (_candidate("CCCCCCCCCCCC", f"C.{known}"),)})

    catalyst_result = plan_multistep_routes(
        "CCCCCCCCCCCC",
        object(),
        _LiteratureIndex((known,), source_role="catalyst"),
        max_depth=2,
        molecular_weight_threshold=20.0,
        expander=expansion,
    )
    untyped_result = plan_multistep_routes(
        "CCCCCCCCCCCC",
        object(),
        _LiteratureIndex((known,), source_role=None),
        max_depth=2,
        molecular_weight_threshold=20.0,
        expander=expansion,
    )

    assert not catalyst_result.routes
    assert not untyped_result.routes
    assert catalyst_result.partial_routes[0].leaves[-1].catalog_role_status in {
        "nonreactant_only",
        "no_match",
    }


def test_legacy_untyped_terminal_requires_explicit_compatibility_switch() -> None:
    known = "CCCCCCCCCC"
    result = plan_multistep_routes(
        "CCCCCCCCCCCC",
        object(),
        _LiteratureIndex((known,), source_role=None),
        max_depth=2,
        molecular_weight_threshold=20.0,
        allow_untyped_literature_terminals=True,
        expander=_expander(
            {"CCCCCCCCCCCC": (_candidate("CCCCCCCCCCCC", f"C.{known}"),)}
        ),
    )

    assert result.routes
    matched = next(
        leaf for leaf in result.routes[0].leaves if leaf.canonical_smiles == known
    )
    assert matched.terminal_evidence == "untyped_literature"


def test_search_does_not_stop_at_first_more_expensive_solved_route() -> None:
    direct = _candidate("CCCCCCCC", "C.C", score=0.0)
    direct = replace(direct, abstraction_level="L0")
    result = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=50.0,
        top_k_routes=1,
        expander=_expander(
            {
                "CCCCCCCC": (
                    direct,
                    _candidate("CCCCCCCC", "CCCCCCC", score=1.0),
                ),
                "CCCCCCC": (_candidate("CCCCCCC", "CC.CCC", score=1.0),),
            }
        ),
    )

    assert result.routes[0].reaction_count == 2
    assert result.routes[0].route_cost < 2.55
    assert result.ranking_policy_definition_id == "multistep_ranking.v1"


def test_search_retains_distinct_paths_to_the_same_leaf_state() -> None:
    result = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=50.0,
        top_k_routes=2,
        expander=_expander(
            {
                "CCCCCCCC": (
                    _candidate("CCCCCCCC", "CCCCCCC", score=0.95),
                    _candidate("CCCCCCCC", "CCCCCCN", score=0.94),
                ),
                "CCCCCCC": (_candidate("CCCCCCC", "CC.CCC", score=0.95),),
                "CCCCCCN": (_candidate("CCCCCCN", "CC.CCC", score=0.95),),
            }
        ),
    )

    assert len(result.routes) == 2
    assert result.diagnostics.retained_alternative_paths >= 1
    assert route_tree_distance(
        result.routes[0].route_tree,
        result.routes[1].route_tree,
    ) > 0.0
    matrix = result.to_dict()["route_postprocessing"]["distance_matrix"]
    assert matrix[0][1] == matrix[1][0]
