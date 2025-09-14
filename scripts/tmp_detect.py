import os, sys
sys.path.insert(0, os.path.abspath('.'))
from chemtools import router
from chemtools.util.rdkit_helpers import rdkit_available
print('rdkit available?', rdkit_available())
reactants = ['Cc1ccc(C)c(Br)c1', 'OB(O)c1ccccc1']
print('reactants:', reactants)
print(router.detect_family(reactants))
