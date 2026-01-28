"""
Test the CLI programmatically by simulating input.
"""
import sys
from pathlib import Path
from io import StringIO

REPO_ROOT = Path(__file__).resolve().parents[0]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Import after adding to path
from app.Cpd_rxn_featurization_cli import _print_readable
from chemtools.featurizers.unified import featurize_molecule, featurize_reaction

print("=" * 80)
print("TEST 1: Molecule with Extended Output")
print("=" * 80)

mol_payload = featurize_molecule("c1ccccc1Br", options={"detailed": True})
_print_readable(mol_payload, show_rdkit=False, show_extended=True)

print("\n" + "=" * 80)
print("TEST 2: Suzuki Reaction with Extended Output")
print("=" * 80)

rxn_payload = featurize_reaction(
    "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1",
    options={"detailed": True, "confirm_coupling_products": True}
)
_print_readable(rxn_payload, show_roles=True, show_rdkit=False, show_extended=True)

print("\n✅ CLI display test completed!")
