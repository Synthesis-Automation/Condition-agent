
from chemtools import detect_reaction
import json

rxn = "Brc1ccccc1.CCNCC>>CCN(CC)c1ccccc1"
result = detect_reaction(rxn)
print(json.dumps(result, indent=2))
