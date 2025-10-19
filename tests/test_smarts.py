from rdkit import Chem
from rdkit.Chem import AllChem, Draw

rxn = AllChem.ReactionFromSmarts('[#6H2,#6H3;D1,D2;!$([#6]-[a]);!$([#6]-[#6]=[#6]):1]-[I:2]>>[#6:1]-[B:3](-[O]-[*])(-[O]-[*])')
img = Draw.ReactionToImage(rxn, subImgSize=(300,200))
img.save("applicability.png")