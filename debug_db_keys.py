
from chemtools.recommend.recommender import _load_hte_database_cached, PROJECT_ROOT
import pandas as pd

db_path = PROJECT_ROOT / "data" / "protocol_db_v2"
df, indexed, patterns, transformation_indices = _load_hte_database_cached(str(db_path))

print(f"Total transformation keys: {len(transformation_indices)}")
for key in transformation_indices.keys():
    if "Ar-Br" in key:
        print(f"Key: '{key}'")

