
print("Script start")
import sys
import traceback

try:
    from chemtools.recommend.recommender import HTERecommender, PROJECT_ROOT
    from pathlib import Path
    
    db_path = PROJECT_ROOT / "data" / "protocol_db_v2"
    recommender = HTERecommender(db_path)
    
    rxn = "O=C1CCCC1.Brc2ccc(C(C)=O)cc2>>O=C(C)c3ccc(C4C(CCC4)=O)cc3"
    print(f"Query: {rxn}")
    
    # Passing default min_experiments (2) but with protocols source
    results = recommender.recommend(
        reactant_a_smiles="O=C1CCCC1",
        reactant_b_smiles="Brc2ccc(C(C)=O)cc2",
        product_smiles="O=C(C)c3ccc(C4C(CCC4)=O)cc3",
        top_k=5,
        source_group="protocols"
    )
    
    print(f"Recommendations Found: {len(results.recommendations)}")
    
    for i, rec in enumerate(results.recommendations):
        print(f"[{i+1}] ID: {rec.reaction_id} | Type: {rec.reaction_type}")

except Exception:
    traceback.print_exc()
print("Script end")
