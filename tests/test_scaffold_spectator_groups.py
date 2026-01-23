import csv
import json
from pathlib import Path

import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util import rdkit_helpers


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_scaffold_motifs_can_appear_in_spectator_groups() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    sample_path = repo_root / "examples" / "sample_reactions.csv"
    scaffold_path = repo_root / "chemtools" / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"

    assert sample_path.exists(), f"Missing sample reactions CSV: {sample_path}"
    assert scaffold_path.exists(), f"Missing scaffold motifs file: {scaffold_path}"

    scaffold_ids = set()
    payload = json.loads(scaffold_path.read_text(encoding="utf-8"))
    for entry in payload.get("compounds", []) or []:
        if isinstance(entry, dict) and entry.get("id"):
            scaffold_ids.add(str(entry["id"]).strip())

    assert scaffold_ids, "No scaffold motif IDs loaded"

    found = None
    with sample_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rxn = (row.get("rxn_smiles_clean") or "").strip()
            if not rxn:
                continue
            result = featurize_reaction(
                rxn,
                options={
                    "include_agent_roles": False,
                    "include_roles": False,
                    "include_rdkit_props": False,
                },
            )
            reaction = result.get("reaction") if isinstance(result, dict) else None
            if not isinstance(reaction, dict):
                continue
            aggregates = reaction.get("aggregates") or {}
            spectators = aggregates.get("spectator_groups_combined") or []
            hit = sorted(set(spectators) & scaffold_ids)
            if hit:
                found = {"reaction": rxn, "spectators": hit}
                break

    assert found is not None, "No scaffold motifs found in spectator groups for sample reactions"
