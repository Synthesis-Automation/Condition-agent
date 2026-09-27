"""Display label cleanup must preserve explicit names and never invent chemistry."""

import pytest

from visualization.scheme_labels import compact_condition_labels, compact_yield_label


@pytest.mark.parametrize(("texts", "reactants", "expected"), [
    (["Core 26.3 g (76.8 mmol); stannane 36.0 g (99.7 mmol); Pd[P(t-Bu)3]2 1.96 g (3.84 mmol)."],
     (), ("Pd[P(t-Bu)3]2",)),
    (["NMP 400 mL; 45 °C, 48 h."], (), ("NMP",)),
    (["C14 14.3 g (42.8 mmol) in THF 200 mL; add 1 M HCl, 156 mL."], (), ("THF", "HCl")),
    (["Initial screen: P6 0.10 mmol, aniline 0.20 mmol, DIPEA 0.20 mmol, DMSO 1.0 mL."],
     (), ("aniline", "DIPEA", "DMSO")),
    (["Initial screen: P6 0.10 mmol, aniline 0.20 mmol, DIPEA 0.20 mmol, DMSO 1.0 mL."],
     ("Aniline", "P6"), ("DIPEA", "DMSO")),
    (["NaHCO3, dry THF, 16 h"], (), ("NaHCO3", "THF")),
    (["Pd(PPh3)4 (5 mol%); Pd[P(t-Bu)3]2 1.96 g; THF"],
     (), ("Pd(PPh3)4", "Pd[P(t-Bu)3]2", "THF")),
    (["N,N-dimethylformamide,1,4-dioxane,2,6-lutidine"],
     (), ("N,N-dimethylformamide", "1,4-dioxane", "2,6-lutidine")),
    (["(1,2-bis(diphenylphosphino)ethane)palladium(II) dichloride (0.1 mmol)"],
     (), ("(1,2-bis(diphenylphosphino)ethane)palladium(II) dichloride",)),
    (["Solvent: ethyl acetate", "THF; THF", "dry THF"], (), ("ethyl acetate", "THF")),
    (["H2; N2; O2; Br2; Fe(CO)5"], (), ("H2", "N2", "O2", "Br2", "Fe(CO)5")),
])
def test_explicit_ingredient_labels_and_formula_digits_are_preserved(texts, reactants, expected) -> None:
    assert compact_condition_labels(texts, reactants) == expected


@pytest.mark.parametrize("text", [
    "Ice/brine quench; ethyl acetate extraction; MgSO4 drying; silica, 5–9% ethyl acetate/petroleum ether.",
    "Conditions not supplied; yield unreported; unknown; room temperature; overnight",
    "The mixture was stirred with an appropriate catalyst until the reaction was complete.",
    "See the patent example for reaction details and proposed catalyst choices.",
    "No Pd(PPh3)4; avoid THF; without DIPEA",
    "No NaHCO3 in THF",
    "The mixture was heated in THF",
    "1 M; 156 mL; 45 °C; 48 h; 5 mol%",
    "Core; stannane; C14; P6; substrate 12; compound A1",
])
def test_workup_unknown_conditions_and_narrative_are_omitted(text: str) -> None:
    assert compact_condition_labels([text]) == ()


def test_exact_reactant_exclusion_does_not_guess_aliases_or_other_names() -> None:
    assert compact_condition_labels(["aniline; anisidine; PhNH2"], reactant_names=["aniline"]) == (
        "anisidine", "PhNH2",
    )


@pytest.mark.parametrize(("text", "expected"), [
    ("90% (23 g)", "90%"), ("28 mg, 4%", "4%"), ("Isolated yield: 90 %", "90%"),
    ("80–90%", "80–90%"), ("<5%", "<5%"), ("90% and 90%", "90%"),
    (None, None), ("", None), ("unreported", None), ("quantitative", None),
    ("28 mg", None), ("Not reported; 95% purity", None), ("95% purity", None),
    ("40% conversion", None), ("40% or 60%", None),
])
def test_yield_requires_an_unambiguous_explicit_percentage(text: str | None, expected: str | None) -> None:
    assert compact_yield_label(text) == expected
