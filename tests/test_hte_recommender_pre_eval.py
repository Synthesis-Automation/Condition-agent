import pandas as pd

from chemtools.HTE.recommender import HTERecommender
import chemtools.HTE.recommender as hte_mod


def test_recommend_uses_product_to_identify_spectators(tmp_path, monkeypatch):
    def fake_featurize(smiles: str):
        mapping = {
            "Sc1ccccc1": {
                "motifs": [
                    {"compound_id": "Ar-SH", "alt_compound_ids": []}
                ]
            },
            "Clc1ccc(I)cc1": {
                "motifs": [
                    {"compound_id": "Ar-Cl", "alt_compound_ids": []},
                    {"compound_id": "Ar-I", "alt_compound_ids": []},
                ]
            },
            "Oc1ccccc1Sc1ccc(Cl)cc1": {
                "motifs": [
                    {"compound_id": "Ar-Cl", "alt_compound_ids": []},
                    {"compound_id": "Ar-SAr", "alt_compound_ids": []},
                ]
            },
        }
        return mapping.get(smiles, {"motifs": []})

    monkeypatch.setattr(hte_mod, "featurize_molecule", fake_featurize)

    df = pd.DataFrame(
        [
            {
                "reaction_type": "Ar-I|Ar-SH -> Ar-SAr || Ar-Cl",
                "yield": 80.0,
                "z_score": 1.0,
                "reactant_1": "Ar-SH",
                "reactant_2": "Ar-I",
                "catalyst": "Pd",
                "ligand": "L1",
                "base": "K2CO3",
                "solvent": "DMF",
                "additive": "",
            },
            {
                "reaction_type": "Ar-Cl|Ar-SH -> Ar-SAr || None",
                "yield": 70.0,
                "z_score": 0.5,
                "reactant_1": "Ar-SH",
                "reactant_2": "Ar-Cl",
                "catalyst": "Ni",
                "ligand": "L2",
                "base": "K3PO4",
                "solvent": "MeCN",
                "additive": "",
            },
        ]
    )
    csv_path = tmp_path / "hte.csv"
    df.to_csv(csv_path, index=False)

    rec = HTERecommender(str(csv_path))
    result = rec.recommend(
        reactant_a_smiles="Sc1ccccc1",
        reactant_b_smiles="Clc1ccc(I)cc1",
        product_smiles="Oc1ccccc1Sc1ccc(Cl)cc1",
        top_k=1,
        min_experiments=1,
    )

    assert result.total_matching_experiments == 1
    assert result.recommendations
    assert result.recommendations[0].catalyst == "Pd"
