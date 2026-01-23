import pandas as pd
import pytest

from chemtools.util import rdkit_helpers

from app.A_convert_to_hte_format import process_reaction_dataset


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_reacted_motifs_joined_without_alt_ids(tmp_path) -> None:
    smiles = "Sc1ccccc1Br.Brc1ccccc1C1=NCCN1>>c1ccc2c(c1)Sc1ccccc1N1CCN=C21"

    input_path = tmp_path / "input.csv"
    output_path = tmp_path / "output.csv"

    pd.DataFrame(
        [
            {
                "reaction_smiles": smiles,
                "yield_pct": 50,
            }
        ]
    ).to_csv(input_path, index=False)

    process_reaction_dataset(
        str(input_path),
        str(output_path),
        drop_no_catalyst=False,
    )

    out = pd.read_csv(output_path)
    assert len(out) == 1

    row = out.iloc[0]
    reactant_1 = row["reactant_1"]
    reactant_2 = row["reactant_2"]

    assert set(reactant_1.split("|")) == {"Ar-Br", "Ar-SH"}
    assert set(reactant_2.split("|")) == {"Ar-Br", "RCH2-NHR"}
    assert "Alkyl-NHR" not in reactant_2
    assert pd.isna(row["reactant_3"]) or row["reactant_3"] == ""


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_aromn_scaffold_override_in_converter(tmp_path) -> None:
    smiles = "c1cn[nH]c1.Clc1ccc(I)cc1>>Clc1ccc(-n2cccn2)cc1"

    input_path = tmp_path / "input.csv"
    output_path = tmp_path / "output.csv"

    pd.DataFrame(
        [
            {
                "reaction_smiles": smiles,
                "yield_pct": 50,
            }
        ]
    ).to_csv(input_path, index=False)

    process_reaction_dataset(
        str(input_path),
        str(output_path),
        drop_no_catalyst=False,
    )

    out = pd.read_csv(output_path)
    assert len(out) == 1

    row = out.iloc[0]
    assert "Pyrazole" in str(row["reactant_1"])
    assert "AromN-H" not in str(row["reactant_1"])
    assert "AromN-H" not in str(row["reactant_2"])
