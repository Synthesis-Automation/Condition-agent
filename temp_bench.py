import sys, time
sys.path.insert(0, ".")
from chemtools.recommend.api import recommend
from chemtools.recommend.models import RecommendationRequest

REACTION = 'Brc1ccc(C(F)(F)F)cc1.c1ccc(B(O)O)cc1>>FC(F)(F)c1ccc(-c2ccccc2)cc1'

for strategy in ['literature', 'similarity']:
    t0 = time.perf_counter()
    req = RecommendationRequest(reaction_smiles=REACTION, strategy=strategy, top_k=5)
    result = recommend(req)
    elapsed = time.perf_counter() - t0
    rec_obj = result.recommendation
    recs = (rec_obj.recommendations or []) if rec_obj else []
    by_source = (rec_obj.recommendations_by_source or {}) if rec_obj else {}
    total = getattr(rec_obj, 'total_matching_experiments', 0) if rec_obj else 0
    ptype = getattr(rec_obj, 'predicted_reaction_type', '') if rec_obj else ''
    print(f'=== {strategy.upper()} ({elapsed:.2f}s) ===')
    print(f'  predicted={ptype} matches={total} recs={len(recs)} keys={list(by_source.keys())}')
    for i, r in enumerate(recs[:5], 1):
        print(f'  {i}. {r.catalyst} / {r.ligand} / {r.base}')
    print()
