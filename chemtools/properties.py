from typing import Dict, Any
DATA = {
 "7778-53-2":{"role":"BASE","token":"K3PO4","pKa_DMSO":30.0},
 "1310-58-3":{"role":"BASE","token":"KOH","pKa_water":15.7},
 "7732-18-5":{"role":"SOLVENT","token":"Water","KT":{"alpha":1.17,"beta":0.47,"pi*":1.09}},
 "7681-65-4":{"role":"CATALYST","token":"CuI"},
 "72-52-8":{"role":"LIGAND","token":"Phenanthroline"}
}
def lookup(query: str) -> Dict[str, Any]:
    q=(query or '').strip().lower()
    for uid, rec in DATA.items():
        if q==uid.lower() or q==rec['token'].lower():
            out={"uid":uid, **rec}; return {"found":True, "record":out}
    return {"found":False, "record":None}
