from reactive_taxonomy import featurize_molecule, featurize_reaction, validate_taxonomy


def test_functional_groups_and_site_environments_are_public() -> None:
    result = featurize_molecule("Brc1ccccc1C#N")
    assert result.valid
    assert {group.group_id for group in result.functional_groups} >= {"aryl_halide", "nitrile"}
    assert result.site_environments
    leaving_environment = next(
        environment for environment in result.site_environments
        if environment.site_id == next(site.site_id for site in result.sites if site.site_type == "leaving_group")
    )
    assert leaving_environment.steric["method"] == "graph_local_v1"
    assert any(group["group_id"] == "nitrile" for group in leaving_environment.nearby_groups)


def test_reaction_spectators_exclude_consumed_handle_and_keep_nitrile() -> None:
    reaction = "Brc1ccc(C#N)cc1.OB(O)c1ccccc1>>N#Cc1ccc(-c2ccccc2)cc1"
    result = featurize_reaction(reaction)
    assert result.evidence_quality == "exact_product_reconstruction"
    ids = {group.group_id for group in result.spectator_groups}
    assert "nitrile" in ids
    assert "aryl_halide" not in ids
    nitrile = next(group for group in result.spectator_groups if group.group_id == "nitrile")
    assert nitrile.graph_distance is not None
    assert nitrile.unchanged_evidence == "exact_product_reconstruction_event_exclusion"


def test_functional_group_registry_is_validated() -> None:
    assert validate_taxonomy() == []
