"""Architecture regressions for observation, interpretation, and rendering."""

from dataclasses import asdict
import json

import reactive_taxonomy.reaction_interpretation as interpretation_module
from reactive_taxonomy import (
    build_observation_signature,
    featurize_reaction,
    interpret_reaction,
    load_reaction_operators,
    observe_reaction,
    render_reaction,
)


MAPPED_SUZUKI = (
    "[c:1]1([Br:2])[cH:3][cH:4][cH:5][cH:6][cH:7]1."
    "[c:8]1([B:9]([OH:10])[OH:11])[cH:12][cH:13][cH:14][cH:15][cH:16]1"
    ">>[c:1]1(-[c:8]2[cH:12][cH:13][cH:14][cH:15][cH:16]2)"
    "[cH:3][cH:4][cH:5][cH:6][cH:7]1"
)
UNMAPPED_SUZUKI = (
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
)


def test_observation_is_grammar_free_and_builds_generic_products() -> None:
    observation = observe_reaction(MAPPED_SUZUKI)

    assert observation.valid
    assert observation.edits
    assert observation.topology is not None
    assert observation.core is not None
    signature = build_observation_signature(observation)
    assert signature is not None
    assert signature.named_family is None
    assert signature.compatible_named_families == ()
    payload = json.dumps(asdict(observation), sort_keys=True)
    assert '"grammar_id"' not in payload
    assert '"named_family"' not in payload


def test_interpretation_adds_optional_family_semantics() -> None:
    observation = observe_reaction(UNMAPPED_SUZUKI)
    interpretation = interpret_reaction(observation)

    assert interpretation.selected_candidate is not None
    assert interpretation.selected_candidate.grammar_id == (
        "boron_transfer_coupling"
    )
    assert interpretation.named_family == "suzuki_miyaura"
    assert {partner.role for partner in interpretation.partners} == {
        "electrophile",
        "transfer_partner",
    }
    assert all(partner.role is not None for partner in interpretation.partners)
    signature = build_observation_signature(observation)
    assert signature is not None
    assert all(partner.role is None for partner in signature.partners)
    rendered = render_reaction(observation, interpretation)
    assert rendered.source == "verified_grammar"


def test_signature_and_core_do_not_depend_on_loaded_grammars(monkeypatch) -> None:
    reaction = "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]"
    interpreted = featurize_reaction(reaction)
    monkeypatch.setattr(
        interpretation_module,
        "enumerate_reaction_candidates",
        lambda _reactants: (),
    )
    generic = featurize_reaction(reaction)

    assert interpreted.reaction_signature is not None
    assert generic.reaction_signature is not None
    assert interpreted.reaction_signature.signature_id == (
        generic.reaction_signature.signature_id
    )
    assert interpreted.reaction_core is not None
    assert generic.reaction_core is not None
    assert interpreted.reaction_core.core_id == generic.reaction_core.core_id
    assert interpreted.selected_candidate is not None
    assert generic.named_family is None
    assert interpreted.reaction_label is not None
    assert generic.reaction_label is not None
    assert interpreted.reaction_label.source == "verified_grammar"
    assert generic.reaction_label.source != "verified_grammar"


def test_public_operator_registry_excludes_grammar_metadata() -> None:
    operators = load_reaction_operators()

    assert operators
    assert len({operator.operator_id for operator in operators}) == len(operators)
    assert all(not hasattr(operator, "grammar_ids") for operator in operators)


def test_cycloaddition_uses_generic_topology_renderer() -> None:
    reaction = "CC#C.CN=[N+]=[N-]>>Cc1nnn(C)c1"
    result = featurize_reaction(reaction)

    assert result.reaction_label is not None
    assert result.reaction_label.source == "generic_topology"
    assert result.reaction_label.concise == (
        "C≡C + N=N=N → aromatic 5-membered C₂N₃ ring"
    )


def test_analysis_serializes_one_canonical_reaction_label_contract() -> None:
    result = featurize_reaction(UNMAPPED_SUZUKI)
    payload = asdict(result)

    assert set(payload) & {
        "display_label",
        "reaction_display_label",
        "reaction_label_status",
    } == set()
    assert payload["reaction_label"]["concise"]
    assert payload["reaction_label"]["detailed"]
    assert payload["reaction_label"]["schema_version"] == "2.0"
    assert "reaction_label" not in payload["candidates"][0]
    assert payload["candidates"][0]["grammar_label"]
