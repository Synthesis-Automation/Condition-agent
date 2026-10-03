"""Composite actions retain chemistry, dependency uncertainty and physical cost."""

from dataclasses import replace
import json

import pytest

from core_retrosynthesis import build_generic_library
from core_retrosynthesis.chemistry import canonical_smiles, digest
from core_retrosynthesis.composite_actions import (
    CompositeStrategyDefinition,
    assess_composite_dependency,
    build_composite_strategy_catalog,
    load_composite_action_policy,
    load_composite_strategy_catalog,
    save_composite_strategy_catalog,
    search_composite_actions,
)
from core_retrosynthesis.generic_compiler import analyze_generic_reaction
from core_retrosynthesis.generic_models import (
    GenericDisconnectionCandidate,
    GenericSearchDiagnostics,
)
from tests.core_retrosynthesis_tests.test_generic_diverse_retrosynthesis import (
    _row_from_reaction,
)


# Independently mapped physical steps; map numbers need not agree across steps.
ACTIVATION = (
    "[CH3:1][CH2:2][OH:3].[CH3:4][S:5](=[O:6])(=[O:7])[Cl:8]>>"
    "[CH3:1][CH2:2][O:3][S:5]([CH3:4])(=[O:6])=[O:7]"
)
SUBSTITUTION = (
    "[CH3:11][CH2:12][O:13][S:15]([CH3:14])(=[O:16])=[O:17].[NH3:18]>>"
    "[CH3:11][CH2:12][NH2:18]"
)
NITRATION = (
    "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1.[O:7]=[N+:8]([O-:9])[OH:10]>>"
    "[c:1]1([N+:8](=[O:7])[O-:9])[cH:2][cH:3][cH:4][cH:5][cH:6]1"
)
REDUCTION = (
    "[c:11]1([N+:18](=[O:17])[O-:19])[cH:12][cH:13][cH:14][cH:15][cH:16]1>>"
    "[c:11]1([NH2:18])[cH:12][cH:13][cH:14][cH:15][cH:16]1"
)


def candidate(reaction: str, operator: str = "") -> GenericDisconnectionCandidate:
    """Build a graph-backed candidate without inventing correspondence."""
    precursor, _, target = reaction.split(">")
    identity = analyze_generic_reaction(reaction)
    assert identity is not None
    return GenericDisconnectionCandidate(
        target_smiles=canonical_smiles(target),
        precursor_smiles=canonical_smiles(precursor),
        proposed_reaction_smiles=f"{canonical_smiles(precursor)}>>{canonical_smiles(target)}",
        condition_query_reaction_smiles=reaction,
        transformation_kind=None,
        abstraction_level="L1",
        compiler_engine="test",
        template_id="test:template",
        score=0.8,
        context_similarity=0.0,
        product_similarity=1.0,
        precursor_similarity=1.0,
        template_specificity=1.0,
        independent_reference_support=2,
        forward_validation_status="verified_signature",
        center_transition_key="test:center",
        disconnection_site_key="test:site",
        precedent_reaction_ids=("test:precedent",),
        operator_id=operator or digest("OP1", identity.operator_signature),
    )


def strategy(
    first, second, relationship="handle_progression"
) -> CompositeStrategyDefinition:
    return CompositeStrategyDefinition(
        digest("CRV1OP1", relationship, first.operator_id, second.operator_id),
        relationship,
        first.operator_id,
        second.operator_id,
        ("TRAIN1", "TRAIN2"),
        2,
        (),
    )


@pytest.mark.parametrize(
    "reactions,expected",
    [
        ((ACTIVATION, SUBSTITUTION), "created_handle_consumed"),
        ((NITRATION, REDUCTION), "created_handle_consumed"),
    ],
)
def test_common_sequences_recheck_target_site(reactions, expected):
    first, second = map(candidate, reactions)
    assessment, tree = assess_composite_dependency(
        first, second, first.target_smiles, "handle_progression"
    )
    assert assessment.admitted, assessment
    assert assessment.dependency_class == expected
    assert assessment.evidence["ambiguity_invariant"]
    assert tree.route_kind == "planned" and tree.reaction_count == 2
    assert tree.root.reaction.evidence.evidence_kind == "predicted"
    assert all(not node.terminal for node in tree.root.reaction.children)


def test_independent_sites_are_rejected_on_same_intermediate():
    first = candidate(
        "[CH3:1][CH2:2][CH2:3][CH2:4][CH2:5][OH:6]>>[CH3:1][CH2:2][CH2:3][CH2:4][CH:5]=[O:6]"
    )
    second = candidate(
        "[CH3:11][CH2:12][CH2:13][CH2:14][CH:15]=[O:16].[Br:17][Br:18]>>[CH2:11]([Br:17])[CH2:12][CH2:13][CH2:14][CH:15]=[O:16]"
    )
    result, _ = assess_composite_dependency(
        first, second, first.target_smiles, "same_site_coupled"
    )
    assert not result.admitted and result.status == "independent_sites"


def test_symmetry_dependent_site_overlap_remains_unresolved():
    first = candidate(
        "[OH:1][CH2:2][CH2:3][CH2:4][CH:5]=[O:6]>>[OH:1][CH2:2][CH2:3][CH2:4][CH2:5][OH:6]"
    )
    second = candidate(
        "[OH:11][CH2:12][CH2:13][CH2:14][CH2:15][OH:16]>>[O:11]=[CH:12][CH2:13][CH2:14][CH2:15][OH:16]"
    )
    result, _ = assess_composite_dependency(
        first, second, first.target_smiles, "same_site_coupled"
    )
    assert not result.admitted and result.relationship_class == "lineage_ambiguous"
    assert not result.evidence["ambiguity_invariant"]


def test_conflicting_declared_relationship_and_structures_are_retained():
    first, second = candidate(ACTIVATION), candidate(SUBSTITUTION)
    result, _ = assess_composite_dependency(
        first, second, first.target_smiles, "same_site_coupled"
    )
    assert not result.admitted and result.status == "conflicting"
    altered = replace(second, condition_query_reaction_smiles=REDUCTION)
    result, _ = assess_composite_dependency(
        first, altered, first.target_smiles, "handle_progression"
    )
    assert not result.admitted and result.status == "conflicting"


@pytest.fixture(scope="module")
def composite_library():
    rows = tuple(
        _row_from_reaction(r, reaction_id=f"r:{i}", reference_id=f"patent:{i}")
        for i, r in enumerate((ACTIVATION, SUBSTITUTION, NITRATION, REDUCTION))
    )
    return build_generic_library(
        rows,
        levels=("L0", "L1", "L2"),
        admission_mode="data_driven",
        core_admission_policy="validated_departures",
    )


@pytest.mark.parametrize(
    "reactions", [(ACTIVATION, SUBSTITUTION), (NITRATION, REDUCTION)]
)
def test_real_operator_composition_roundtrips_and_preserves_partners(
    composite_library, reactions
):
    first, second = map(candidate, reactions)
    result = search_composite_actions(
        second.target_smiles, composite_library, (strategy(first, second),)
    )
    assert result.actions, result.to_dict()
    action = result.actions[0].to_dict()
    assert action["physical_step_count"] == action["physical_step_cost"] == 2
    assert action["logical_action_count"] == 1
    assert action["dependency"]["admitted"]
    assert len(action["physical_steps"]) == 2
    assert all(step["precedent_reaction_ids"] for step in action["physical_steps"])
    assert (
        action["condition_compatibility_status"]
        == action["one_pot_status"]
        == "not_assessed"
    )
    if reactions[0] == ACTIVATION:
        assert "N" in action["terminal_precursor_smiles"].split(".")
    assert action["route_tree"]["maximum_depth"] == 2
    reordered = search_composite_actions(
        "NCC" if reactions[0] == ACTIVATION else "c1ccccc1N",
        composite_library,
        (strategy(first, second),),
    )
    assert reordered.to_dict() == result.to_dict()


def test_invalid_atom_maps_do_not_get_catalogue_precedence():
    first, second = candidate(ACTIVATION), candidate(SUBSTITUTION)
    second = replace(
        second,
        condition_query_reaction_smiles=SUBSTITUTION.replace("[CH2:12]", "[CH2:11]"),
    )
    assessment, _ = assess_composite_dependency(
        first, second, first.target_smiles, "handle_progression"
    )
    assert not assessment.admitted
    assert assessment.warnings or assessment.evidence


def test_catalog_identity_and_validation(tmp_path):
    first, second = candidate(ACTIVATION), candidate(SUBSTITUTION)
    definition = strategy(first, second)
    catalog = build_composite_strategy_catalog((definition,))
    shuffled = build_composite_strategy_catalog(
        (replace(definition, training_patent_ids=("TRAIN2", "TRAIN1")),)
    )
    assert shuffled.catalog_id == catalog.catalog_id
    path = tmp_path / "catalog.json"
    save_composite_strategy_catalog(catalog, path)
    assert load_composite_strategy_catalog(path) == catalog
    raw = json.loads(path.read_text())
    assert "cases" not in raw
    raw["strategies"][0]["first_operator_id"] = "changed"
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="identity mismatch"):
        load_composite_strategy_catalog(path)
    with pytest.raises(ValueError, match="duplicate"):
        build_composite_strategy_catalog((definition, definition))
    assert load_composite_action_policy()["physical_step_cost"] == 2


def test_rejected_dependency_does_not_erase_single_step_fallback(composite_library):
    first, second = candidate(ACTIVATION), candidate(SUBSTITUTION)
    result = search_composite_actions(
        second.target_smiles,
        composite_library,
        (strategy(first, second, "same_site_coupled"),),
    )
    assert not result.actions
    assert result.one_step_fallbacks
    assert result.diagnostics.dependency_rejected_count > 0
    assert result.dependency_reviews[0]["dependency"]["status"] == "conflicting"


def test_failed_physical_compatibility_is_filtered_before_scoring(composite_library):
    first, second = candidate(ACTIVATION), candidate(SUBSTITUTION)
    first = replace(first, reaction_compatibility_disposition="reject")

    def searcher(target, _library, **_kwargs):
        values = (
            (second,)
            if target == second.target_smiles
            else (first,)
            if target == first.target_smiles
            else ()
        )
        return values, GenericSearchDiagnostics(validation_attempt_count=len(values))

    result = search_composite_actions(
        second.target_smiles,
        composite_library,
        (strategy(first, second),),
        searcher=searcher,
    )
    assert not result.actions
    assert result.dependency_reviews[0]["dependency"]["status"] == "blocked"
