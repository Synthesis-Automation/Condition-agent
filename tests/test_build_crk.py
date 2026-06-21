"""Tests for build_crk() and CrkResult in chemtools.reaction.featurize."""
import pytest

from chemtools.reaction.featurize import (
    CrkResult,
    build_crk,
    featurize_reaction,
    get_crk_options,
)
from chemtools.reaction.typing import extract_reaction_key

# Standard Suzuki coupling: bromobenzene + phenylboronic acid → biphenyl
_SUZUKI = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"


class TestCrkResultShape:
    def test_fields_present(self):
        result = build_crk(_SUZUKI)
        assert hasattr(result, "reacted_motifs")
        assert hasattr(result, "formed_motifs")
        assert hasattr(result, "spectator_motifs")
        assert hasattr(result, "reaction_key")
        assert hasattr(result, "aggregates")

    def test_field_types(self):
        result = build_crk(_SUZUKI)
        assert isinstance(result.reacted_motifs, list)
        assert isinstance(result.formed_motifs, list)
        assert isinstance(result.spectator_motifs, list)
        assert isinstance(result.reaction_key, str)
        assert isinstance(result.aggregates, dict)

    def test_reaction_key_nonempty_for_valid_reaction(self):
        result = build_crk(_SUZUKI)
        assert len(result.reaction_key) > 0


class TestCrkResultParity:
    """build_crk() must produce the same motif lists as featurize_reaction(options=get_crk_options())."""

    def _full(self, smiles: str) -> dict:
        return featurize_reaction(smiles, options=get_crk_options())

    def test_reacted_motifs_parity(self):
        crk = build_crk(_SUZUKI)
        full_agg = self._full(_SUZUKI).get("aggregates", {})
        assert crk.reacted_motifs == full_agg.get("reacted_motifs", [])

    def test_formed_motifs_parity(self):
        crk = build_crk(_SUZUKI)
        full_agg = self._full(_SUZUKI).get("aggregates", {})
        assert crk.formed_motifs == full_agg.get("formed_motifs", [])

    def test_spectator_motifs_parity(self):
        crk = build_crk(_SUZUKI)
        full_agg = self._full(_SUZUKI).get("aggregates", {})
        assert crk.spectator_motifs == full_agg.get("spectator_motifs", [])

    def test_reaction_key_core_parity(self):
        """The reaction key (before any | events: suffix) must match."""
        crk_key = build_crk(_SUZUKI).reaction_key
        full_key = self._full(_SUZUKI).get("reaction_key", "")
        # build_crk skips event-signature annotation; strip that suffix if present
        crk_core = crk_key.split("| events:")[0].strip()
        full_core = full_key.split("| events:")[0].strip()
        assert crk_core == full_core


class TestExtractReactionKeyIntegration:
    def test_returns_four_tuple(self):
        result = extract_reaction_key(_SUZUKI)
        assert isinstance(result, tuple)
        assert len(result) == 4

    def test_tuple_types(self):
        reacted, formed, spectator, key = extract_reaction_key(_SUZUKI)
        assert isinstance(reacted, list)
        assert isinstance(formed, list)
        assert isinstance(spectator, list)
        assert isinstance(key, str)

    def test_key_nonempty(self):
        _, _, _, key = extract_reaction_key(_SUZUKI)
        assert len(key) > 0


class TestBuildCrkRobustness:
    def test_invalid_smiles_does_not_raise(self):
        """build_crk should return a valid CrkResult even for malformed input."""
        try:
            result = build_crk("not>>valid")
            assert isinstance(result, CrkResult)
        except Exception:
            pytest.fail("build_crk() raised unexpectedly on invalid SMILES")

    def test_empty_reaction_key_for_trivial_smiles(self):
        result = build_crk("not>>valid")
        assert isinstance(result.reaction_key, str)
        assert isinstance(result.reacted_motifs, list)
