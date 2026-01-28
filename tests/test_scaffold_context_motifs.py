import csv
import json
from pathlib import Path

import pytest

from chemtools.featurizers.unified import featurize_molecule
from chemtools.util import rdkit_helpers


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_scaffold_motifs_are_context_only() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    sample_path = repo_root / "examples" / "sample_compounds.csv"
    scaffold_path = repo_root / "chemtools" / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"

    assert sample_path.exists(), f"Missing sample compounds CSV: {sample_path}"
    assert scaffold_path.exists(), f"Missing scaffold motifs file: {scaffold_path}"

    scaffold_ids = set()
    payload = json.loads(scaffold_path.read_text(encoding="utf-8"))
    for entry in payload.get("compounds", []) or []:
        if isinstance(entry, dict) and entry.get("id"):
            scaffold_ids.add(str(entry["id"]).strip())

    assert scaffold_ids, "No scaffold motif IDs loaded"

    with sample_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    assert rows, "sample_compounds.csv is empty"

    for row in rows:
        smiles = (row.get("smiles") or "").strip()
        if not smiles:
            continue
        result = featurize_molecule(smiles)
        motif_ids = {
            str(m.get("compound_id")).strip()
            for m in (result.get("motifs") or [])
            if m.get("compound_id")
        }
        context_ids = {
            str(m.get("compound_id")).strip()
            for m in (result.get("context_motifs") or [])
            if m.get("compound_id")
        }

        assert motif_ids.isdisjoint(scaffold_ids), f"Scaffold ID leaked into motifs: {motif_ids & scaffold_ids}"
        assert context_ids.issubset(scaffold_ids), f"Non-scaffold ID found in context_motifs: {context_ids - scaffold_ids}"
