from chemtools.reaction.featurize import featurize_reaction
from chemtools.molecule.featurize import featurize_molecule


def test_boc_amine_detects_alkyl_oconh2_motif() -> None:
    bundle = featurize_molecule("CC(C)(C)OC(N)=O")
    motif_ids = {str(m.get("compound_id") or "") for m in (bundle.get("motifs") or [])}
    context_ids = {str(m.get("compound_id") or "") for m in (bundle.get("context_motifs") or [])}
    assert "Alkyl-OCONH2" in (motif_ids | context_ids)


def test_c_n_substitution_with_boc_amine_keeps_nucleophile_in_roles() -> None:
    reaction = "CC(C)(C)OC(N)=O.Brc1ccc2c(ccc3ccccc32)c1>>CC(C)(C)OC(=O)Nc1ccc2c(ccc3ccccc32)c1"
    result = featurize_reaction(reaction, options={"include_reaction_type": True, "detailed": True})
    role_summary = ((result.get("extended") or {}).get("role_classification") or {}).get("reactants") or {}
    reactants = role_summary.get("reactants") or []
    categories = {str(row.get("category") or "") for row in reactants if isinstance(row, dict)}
    assert "Ar-Br" in categories
    assert "Alkyl-OCONH2" in categories
