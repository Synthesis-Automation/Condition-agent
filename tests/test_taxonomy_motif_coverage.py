from __future__ import annotations

from chemtools.taxonomy.reaction_catalog import load_reaction_catalog


def test_priority_motifs_are_covered_by_reaction_slots() -> None:
    definitions, _ = load_reaction_catalog()
    covered = set()
    for defn in definitions.values():
        for slot in defn.reactants.values():
            covered.update(slot.allowed)
        for slot in defn.products.values():
            covered.update(slot.allowed)

    assert "Alkyl-N(R)CO2R" in covered
    assert "Ar-C=N" in covered
    assert "Alkyl-H" in covered
    assert "Alkyl-CF3" in covered


def test_reductive_amination_and_trifluoromethylation_include_new_motifs() -> None:
    definitions, _ = load_reaction_catalog()

    reductive = definitions["Reductive_amination"]
    reductive_reactants = reductive.reactants["electrophile"].allowed
    reductive_products = reductive.products["product"].allowed
    assert "Ar-C=N" in reductive_reactants
    assert "Alkyl-N(R)CO2R" in reductive_products

    cf3 = definitions["Trifluoromethylation"]
    cf3_reactants = cf3.reactants["substrate"].allowed
    cf3_products = cf3.products["product"].allowed
    assert "Alkyl-H" in cf3_reactants
    assert "Alkyl-CF3" in cf3_products
