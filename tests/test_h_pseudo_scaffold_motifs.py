import json
from pathlib import Path

from chemtools.featurizers.calculable import classify_reactant_smiles, get_reactant_type_features
from chemtools.featurizers.molecule import featurize_molecule
from chemtools.taxonomy import loader as taxonomy_loader


def _taxonomy_data_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "chemtools" / "taxonomy" / "data"


def test_h_scaffold_exists_for_allowlisted_pseudo_motifs() -> None:
    path = _taxonomy_data_dir() / "organic_groups.v1.3.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    groups = [entry for entry in (payload.get("groups") or []) if isinstance(entry, dict)]
    h_entries = [entry for entry in groups if str(entry.get("id") or "") == "H"]
    assert len(h_entries) == 1
    assert str(h_entries[0].get("kind") or "") == "scaffold"


def test_only_allowlisted_h_compounds_exist() -> None:
    payload = taxonomy_loader.load_organic_compounds()
    compounds = [entry for entry in (payload.get("compounds") or []) if isinstance(entry, dict)]

    found = set()
    for entry in compounds:
        a_ref = str(entry.get("A") or "")
        if a_ref != "H":
            continue
        b_ref = str(entry.get("B") or "")
        compound_id = str(entry.get("id") or "")
        if not compound_id:
            compound_id = f"{a_ref}{b_ref}" if b_ref.startswith("-") else f"{a_ref}-{b_ref}"
        found.add(compound_id)
        assert bool(entry.get("smarts") or entry.get("smarts_any"))

    assert found == {"H-NH2", "H-CONH2"}


def test_ammonia_maps_to_h_nh2() -> None:
    features = get_reactant_type_features("N")
    assert "H-NH2" in set(features.get("member_types") or [])
    best = classify_reactant_smiles("N")
    assert isinstance(best, dict)
    assert str(best.get("category") or "") == "H-NH2"


def test_formamide_maps_to_h_conh2() -> None:
    bundle = featurize_molecule("NC=O")
    motif_ids = {str(m.get("compound_id") or "") for m in (bundle.get("motifs") or [])}
    context_ids = {str(m.get("compound_id") or "") for m in (bundle.get("context_motifs") or [])}
    assert "H-CONH2" in (motif_ids | context_ids)


def test_h_pseudo_motifs_do_not_overmatch() -> None:
    methylamine = featurize_molecule("CN")
    methylamine_ids = {str(m.get("compound_id") or "") for m in (methylamine.get("motifs") or [])}
    assert "H-NH2" not in methylamine_ids

    acetamide = featurize_molecule("CC(=O)N")
    acetamide_ids = {str(m.get("compound_id") or "") for m in (acetamide.get("motifs") or [])}
    assert "H-CONH2" not in acetamide_ids
