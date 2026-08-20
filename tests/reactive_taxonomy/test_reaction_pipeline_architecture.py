"""Architecture regressions for observation, interpretation, and rendering."""

from dataclasses import asdict
import json

import reactive_taxonomy.reaction_patterns as pattern_module
import reactive_taxonomy.reaction_api as reaction_api_module
from reactive_taxonomy import (
    build_observation_signature,
    featurize_reaction,
    interpret_reaction,
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


def test_observation_is_interpretation_free_and_builds_generic_products() -> None:
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
    assert '"annotation_id"' not in payload
    assert '"named_family"' not in payload
    assert '"transformation_class"' not in payload
    assert '"interpretation_label"' not in payload
    assert all(
        not hasattr(component, "molecule_analysis")
        for component in observation.reactants + observation.products
    )


def test_interpretation_adds_optional_family_semantics() -> None:
    observation = observe_reaction(UNMAPPED_SUZUKI)
    interpretation = interpret_reaction(observation)

    assert interpretation.primary_pattern_id == "organoboron_c_c_coupling_like"
    assert interpretation.named_family == "suzuki_miyaura"
    assert interpretation.partners == ()
    signature = build_observation_signature(observation)
    assert signature is not None
    assert all(partner.role is None for partner in signature.partners)
    from reactive_taxonomy import build_reaction_render_context

    rendered = render_reaction(build_reaction_render_context(
        observation=observation,
        interpretation=interpretation,
    ))
    assert rendered.basis == "literal_edits"


def test_signature_and_core_do_not_depend_on_optional_annotations(monkeypatch) -> None:
    reaction = "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]"
    interpreted = featurize_reaction(reaction)
    monkeypatch.setattr(
        pattern_module,
        "load_reaction_pattern_definitions",
        lambda: (),
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
    assert generic.named_family is None
    assert interpreted.reaction_label is not None
    assert generic.reaction_label is not None
    assert interpreted.reaction_label.text == generic.reaction_label.text


def test_signature_is_finalized_before_optional_molecular_interpretation(
    monkeypatch,
) -> None:
    """Keep graph identity construction ahead of optional annotations."""
    calls = []
    original_observation = reaction_api_module.build_reaction_observation
    original_signature = reaction_api_module.build_observation_signature
    original_interpretation = reaction_api_module.interpret_parsed_molecules

    def record_observation(*args, **kwargs):
        calls.append("observation")
        return original_observation(*args, **kwargs)

    def record_signature(*args, **kwargs):
        calls.append("signature")
        return original_signature(*args, **kwargs)

    def record_interpretation(*args, **kwargs):
        calls.append("molecular_interpretation")
        return original_interpretation(*args, **kwargs)

    monkeypatch.setattr(
        reaction_api_module,
        "build_reaction_observation",
        record_observation,
    )
    monkeypatch.setattr(
        reaction_api_module,
        "build_observation_signature",
        record_signature,
    )
    monkeypatch.setattr(
        reaction_api_module,
        "interpret_parsed_molecules",
        record_interpretation,
    )

    result = featurize_reaction(MAPPED_SUZUKI)

    assert result.reaction_signature is not None
    assert calls[:3] == [
        "observation",
        "signature",
        "molecular_interpretation",
    ]


def test_molecular_reactivity_hypotheses_cannot_change_observation(monkeypatch) -> None:
    reaction = UNMAPPED_SUZUKI
    baseline = featurize_reaction(reaction)
    monkeypatch.setattr(
        reaction_api_module,
        "interpret_parsed_molecules",
        lambda parsed, **_: parsed,
    )
    without_annotations = featurize_reaction(reaction)

    assert asdict(baseline.observation) == asdict(without_annotations.observation)
    assert baseline.reaction_core is not None
    assert without_annotations.reaction_core is not None
    assert baseline.reaction_core.core_id == without_annotations.reaction_core.core_id
    assert baseline.reaction_signature is not None
    assert without_annotations.reaction_signature is not None
    assert baseline.reaction_signature.signature_id == (
        without_annotations.reaction_signature.signature_id
    )


def test_cycloaddition_uses_generic_topology_renderer() -> None:
    reaction = "CC#C.CN=[N+]=[N-]>>Cc1nnn(C)c1"
    result = featurize_reaction(reaction)

    assert result.reaction_label is not None
    assert result.reaction_label.basis == "ring_topology"
    assert result.reaction_label.text == (
        "R¹–C≡C–H + R–N₃ → aromatic 5-membered C₂N₃ ring"
    )


def test_analysis_serializes_one_canonical_reaction_label_contract() -> None:
    result = featurize_reaction(UNMAPPED_SUZUKI)
    payload = asdict(result)

    assert set(payload) & {
        "display_label",
        "reaction_display_label",
        "reaction_label_status",
    } == set()
    assert payload["reaction_label"]["text"]
    assert payload["reaction_label"]["basis"]
    assert payload["reaction_label"]["schema_version"] == "1.1"
    assert payload["interpretation"]["pattern_matches"]
