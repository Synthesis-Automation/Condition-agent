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

    assert "Alkyl-N(R)COOR" not in covered
    assert "Ar-C=N" in covered
    assert "Alkyl-H" in covered
    assert "Alkyl-CF3" in covered


def test_reductive_amination_and_trifluoromethylation_include_new_motifs() -> None:
    definitions, _ = load_reaction_catalog()

    reductive = definitions["Reductive_amination"]
    reductive_reactants = reductive.reactants["electrophile"].allowed
    reductive_products = reductive.products["product"].allowed
    assert "Ar-C=N" in reductive_reactants
    assert "Alkyl-N(R)COOR" not in reductive_products
    assert "Ar-NHR" not in reductive_products
    assert "Bn-NHR" in reductive_products

    cf3 = definitions["Trifluoromethylation"]
    cf3_reactants = cf3.reactants["substrate"].allowed
    cf3_products = cf3.products["product"].allowed
    assert "Alkyl-H" in cf3_reactants
    assert "Alkyl-CF3" in cf3_products


def test_alkyl_substitution_includes_phosphoryl_composite_motifs() -> None:
    definitions, _ = load_reaction_catalog()
    sn2 = definitions["Alkyl_Nucleophilic_Substitution"]
    nucleophile_allowed = sn2.reactants["nucleophile"].allowed
    product_allowed = sn2.products["product"].allowed

    assert "Ar-PO2OH" in nucleophile_allowed
    assert "Alkyl-PO2OH" in nucleophile_allowed
    assert "RCH2-PO2NH2" in nucleophile_allowed
    assert "Ar-PO2OR" in product_allowed
    assert "Alkyl-PO2NHR" in product_allowed
    assert "R3C-PO2NR2" in product_allowed

