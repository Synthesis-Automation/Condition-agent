from typing import Dict, Any, List
import os, json
DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "reactions_sample.jsonl")
def _load() -> List[Dict[str, Any]]:
    rows=[]; 
    if os.path.exists(DATA_PATH):
        with open(DATA_PATH,"r",encoding="utf-8") as f:
            for line in f:
                line=line.strip()
                if not line: continue
                try: rows.append(json.loads(line))
                except: pass
    return rows
def _sim(a: Dict[str, Any], b: Dict[str, Any]) -> float:
    if a.get("bin")==b.get("bin"): return 1.0
    an=(a.get("nuc_class") or '').lower(); bn=(b.get("nuc_class") or '').lower()
    return 0.6 if an and bn and an==bn else 0.0
def knn(family: str, features: Dict[str, Any], k: int=50, relax: Dict[str, Any]|None=None)->Dict[str,Any]:
    rows=_load(); scored=[]
    for r in rows:
        if r.get("rxn_type")!="Ullmann C–N": continue
        s=_sim(features, r.get("features", {}))
        if s>0.0: scored.append((s,r))
    scored.sort(key=lambda x:(-x[0], -(x[1].get("yield_value") or 0)))
    top=[r for s,r in scored[:k]]; support=len(top)
    proto=f"proto_{support}_{abs(hash(features.get('bin','NA'))) % 10000}"
    precedents=[{"reaction_id":r.get("reaction_id"),"yield":r.get("yield_value"),
                 "core":r.get("condition_core"),"base_uid":r.get("base_uid"),
                 "solvent_uid":r.get("solvent_uid"),"T_C":r.get("T_C"),"time_h":r.get("time_h")} for r in top[:10]]
    return {"prototype_id": proto, "support": support, "precedents": precedents}
