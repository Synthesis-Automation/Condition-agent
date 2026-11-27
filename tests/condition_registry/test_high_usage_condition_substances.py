"""Identity regressions for frequently referenced salts, complexes, and formulations."""
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from condition_registry import ConditionRegistry, get_registry


@pytest.mark.parametrize(("identifier", "cas", "formula", "old_fragment"), [
    ("LiBH4", "16949-15-8", "H4BLi", "[BH4-]"),
    ("NaBH4", "16940-66-2", "H4BNa", "B"),
    ("LiAlH4", "16853-85-3", "H4AlLi", "[AlH3]"),
    ("KNO2", "7758-09-0", "KNO2", "O=N[O-]"),
    ("NaNO2", "7632-00-0", "NNaO2", "O=NO"),
    ("Ir[dF(CF3)ppy]2(dtbbpy)PF6", "870987-63-6", "C42H34F16IrN4P", "F[P-](F)(F)(F)(F)F"),
    ("Pd2dba3", "51364-51-3", "C51H42O3Pd2", "P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    ("PdCl2(PPh3)2", "13965-03-2", "C36H30Cl2P2Pd", "C1CCC2=NCCCN2CC1"),
    ("Cu(OTf)2", "34946-82-2", "C2CuF6O6S2", "O=S(=O)([O-])C(F)(F)F"),
    ("Sodium thiosulfate pentahydrate", "10102-17-7", "H10Na2O8S2", "O=S([O-])([O-])=S"),
    ("Lithium triethylborohydride", "22560-16-3", "C6H16BLi", "CC[B-](CC)CC"),
    ("lithium tri-sec-butylborohydride", "38721-52-7", "C12H28BLi", "CCC(C)[B-](C(C)CC)C(C)CC"),
])
def test_high_usage_materials_require_the_whole_composition(
    identifier: str, cas: str, formula: str, old_fragment: str,
) -> None:
    registry: ConditionRegistry = get_registry()
    named = registry.resolve_identifier(identifier, identifier_type="name")
    assert named.status == "resolved" and named.substance.cas == cas
    mol = Chem.MolFromSmiles(named.substance.smiles)
    assert rdMolDescriptors.CalcMolFormula(mol) == formula
    assert Chem.GetFormalCharge(mol) == 0
    partial = registry.resolve_identifier(old_fragment, identifier_type="smiles")
    partial_ids = {partial.substance.substance_id} if partial.substance else set(partial.candidates)
    assert named.substance.substance_id not in partial_ids


def test_superoxide_includes_the_radical_anion_and_potassium() -> None:
    identity = get_registry().resolve(cas="12030-88-5")
    mol = Chem.MolFromSmiles(identity.substance.smiles)
    assert rdMolDescriptors.CalcMolFormula(mol) == "KO2"
    assert Chem.GetFormalCharge(mol) == 0
    assert sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()) == 1


def test_recovered_citrate_dihydrate_keeps_the_original_registry_id() -> None:
    registry = get_registry()
    by_cas = registry.resolve(cas="6132-04-3")
    by_name = registry.resolve(name="Trisodium citrate dihydrate")
    assert by_cas.status == by_name.status == "resolved"
    assert by_cas.substance == by_name.substance
    assert by_cas.substance.substance_id == "registry-row:26345"
    mol = Chem.MolFromSmiles(by_cas.substance.smiles)
    assert rdMolDescriptors.CalcMolFormula(mol) == "C6H9Na3O9"
    assert Chem.GetFormalCharge(mol) == 0
    assert len(Chem.GetMolFrags(mol)) == 6


@pytest.mark.parametrize(("smiles", "cas"), [
    ("CS(C)=O", "67-68-5"),
    ("C[S+](C)[O-]", "67-68-5"),
    ("[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]", "2206-27-1"),
])
def test_verified_dmso_graph_alias_preserves_ordinary_and_isotope_identities(
    smiles: str, cas: str,
) -> None:
    result = get_registry().resolve_identifier(smiles, identifier_type="smiles")
    assert result.status == "resolved" and result.substance.cas == cas


@pytest.mark.parametrize("cas", ["917-54-4", "97-93-8", "1148148-01-9"])
def test_uncertain_graphs_remain_resolvable_by_cas(cas: str) -> None:
    identity = get_registry().resolve(cas=cas)
    assert identity.status == "resolved"
    assert identity.substance.smiles is None
