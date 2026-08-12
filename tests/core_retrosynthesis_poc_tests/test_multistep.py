"""Bounded deterministic multistep retrosynthesis regressions."""

from __future__ import annotations

from core_retrosynthesis_poc import (
    GenericDisconnectionCandidate,
    plan_multistep_routes,
)
from cas_tools import MoleculeIndexMatch


class _LiteratureIndex:
    def __init__(self, known: tuple[str, ...] = ()) -> None:
        from cas_tools import molecule_identity

        self.known = {
            molecule_identity(value).canonical_smiles for value in known
        }

    def lookup(self, identity, *, provenance_limit: int = 5):
        if identity.canonical_smiles not in self.known:
            return None
        return MoleculeIndexMatch(
            canonical_smiles=identity.canonical_smiles,
            inchi_key=identity.inchi_key,
            occurrence_count=1,
            source_records=({"reaction_id": "known"},),
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


def test_root_never_terminates_without_one_disconnection() -> None:
    result = plan_multistep_routes(
        "CC",
        object(),
        _LiteratureIndex(("CC",)),
        max_depth=2,
        molecular_weight_threshold=150.0,
        expander=_expander(
            {"CC": (_candidate("CC", "C.O"),)}
        ),
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
            {
                "CCCCCCCCCCCC": (
                    _candidate("CCCCCCCCCCCC", f"C.{known}"),
                )
            }
        ),
    )

    matched = next(
        leaf for leaf in result.routes[0].leaves if leaf.canonical_smiles == known
    )
    assert matched.molecular_weight > 20.0
    assert matched.terminal_reasons == ("literature_catalog_match",)
    assert matched.literature_match is not None


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

