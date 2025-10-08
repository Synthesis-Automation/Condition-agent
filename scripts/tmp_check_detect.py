import os, sys
sys.path.insert(0, os.path.abspath('.'))
from chemtools.reaction_type_detector import detect_reaction_type, is_available
s = "Cc1ccc(C)c(Br)c1.OB(O)c1ccccc1>>Cc1ccc(C)c(-c2ccccc2)c1"
print('rxn_insight available?', is_available())
print(detect_reaction_type(s))
