from __future__ import annotations

from chemtools.featurizers import reaction_path


def test_analyze_reaction_featurization_shares_cleanup_and_options(monkeypatch) -> None:
    seen: dict[str, object] = {}

    def fake_featurize_reaction(reaction_smiles: str, options=None):
        seen["reaction_smiles"] = reaction_smiles
        seen["options"] = dict(options or {})
        return {
            "reaction_type": "C_N_Coupling",
            "reaction_key": "|Ar-Cl -> AromN-H | bond_formed: C-N | bond_broken: C-Cl",
        }

    monkeypatch.setattr(reaction_path, "featurize_reaction", fake_featurize_reaction)
    reaction_path.cached_featurize_reaction.cache_clear()
    try:
        result = reaction_path.analyze_reaction_featurization(
            "Cl.Clc1ccncc1.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1.Cl",
            base_options={"skip_bond_analysis": True},
            llm_assist_options={
                "enabled": True,
                "provider": "openai",
                "model": "gpt-test",
            },
            cleanup=True,
        )
    finally:
        reaction_path.cached_featurize_reaction.cache_clear()

    assert seen["reaction_smiles"] == "Clc1ccncc1>>CN(C)c1ccncc1"
    options = seen["options"]
    assert isinstance(options, dict)
    assert options["include_reaction_type"] is True
    assert options["detailed"] is True
    assert options["motif_site_filter"] == "substituent"
    assert options["skip_bond_analysis"] is True
    assert options["llm_assist"]["model"] == "gpt-test"
    assert result.cleanup_stats["cleanup_applied"] == 1
    assert result.detected_reaction_type == "C_N_Coupling"
