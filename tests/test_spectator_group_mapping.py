from app import A_convert_to_hte_format as hte_convert
from chemtools.featurizers.formatters.utils import group_id_from_motif_id
from chemtools.featurizers.unified import featurize_reaction
from chemtools.recommend import recommender as hte
from chemtools.util import rdkit_helpers


def test_formatter_group_id_preserves_hetero_scaffold_for_c_h() -> None:
    assert group_id_from_motif_id("HeteroAr-H") == "HeteroAr"
    assert group_id_from_motif_id("AromN-H") == "AromN"
    assert group_id_from_motif_id("Ar-H") == "H"


def test_recommender_spectator_groups_include_hetero_scaffold() -> None:
    groups = hte._spectator_groups_from_motifs({"HeteroAr-H", "Ar-CF3", "Ar-H"})
    assert "HeteroAr" in groups
    assert "CF3" in groups
    assert "H" not in groups


def test_converter_group_id_preserves_hetero_scaffold_for_c_h() -> None:
    assert hte_convert._group_id_from_motif_id("HeteroAr-H") == "HeteroAr"
    assert hte_convert._group_id_from_motif_id("AromN-H") == "AromN"


def test_reaction_aggregates_include_named_scaffold_spectator_group() -> None:
    if not rdkit_helpers.rdkit_available():
        return

    rxn = "Brc1ccc(cc1)c2ccncc2.COc1ccc(B(O)O)cc1>>COc1ccc(-c1ccc(cc1)c2ccncc2)cc1"
    payload = featurize_reaction(rxn, options={"detailed": True, "motif_site_filter": "substituent"})
    aggregates = payload.get("aggregates") or {}
    groups = set(aggregates.get("spectator_groups_ranked") or [])
    reaction_key = str(payload.get("reaction_key") or "")

    assert "Pyridine" in reaction_key
    assert "Pyridine" in groups


def test_pyridine_context_is_not_misclassified_as_reacted() -> None:
    if not rdkit_helpers.rdkit_available():
        return

    rxn = "NS(=O)(=O)c1cccc(F)n1.COc1ccc(B(O)O)cc1>>COc1ccc(NS(=O)(=O)c2cccc(F)n2)cc1"
    payload = featurize_reaction(rxn, options={"detailed": True, "motif_site_filter": "substituent"})
    aggregates = payload.get("aggregates") or {}

    reacted = set(aggregates.get("reacted_motifs") or [])
    spectators = set(aggregates.get("spectator_motifs") or [])
    groups = set(aggregates.get("spectator_groups_ranked") or [])

    assert "Pyridine" not in reacted
    assert "Pyridine" in spectators
    assert "Pyridine" in groups
