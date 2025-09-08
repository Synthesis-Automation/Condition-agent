# -*- coding: utf-8 -*-
import numpy as np
import chemtools.precedent as prec
import chemtools.reaction_similarity as rs

rows = [
    {"reaction_id":"R1","rxn_type":"Ullmann C–N","yield_value":50.0,
     "features":{"LG":"Br","nuc_class":"aniline","bin":"LG:Br|NUC:aniline"},
     "reaction_smiles":"LOW_RS"},
    {"reaction_id":"R2","rxn_type":"Ullmann C–N","yield_value":50.0,
     "features":{"LG":"Cl","nuc_class":"phenol","bin":"LG:Cl|NUC:phenol"},
     "reaction_smiles":"HIGH_RS"},
]

prec._knn_cached.cache_clear()
prec._load = lambda: rows
rs.drfp_available = lambda : True
rs.encode_drfp_cached = lambda s, n_bits=4096, radius=3: (np.ones(8, dtype='uint8') if 'HIGH' in s else np.zeros(8, dtype='uint8'))

features = {"LG":"Br","nuc_class":"aniline","bin":"LG:Br|NUC:aniline"}
relax = {"use_drfp":True, "reaction_smiles":"HIGH_QUERY", "drfp_weight":0.9, "strict_bin":False, "min_candidates":2}
print(prec.knn("Ullmann_CN", features, k=2, relax=relax))
