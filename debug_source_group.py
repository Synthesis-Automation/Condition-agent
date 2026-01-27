
from chemtools.recommend.recommender import _load_hte_database_cached, PROJECT_ROOT
import pandas as pd

db_path = PROJECT_ROOT / "data" / "protocol_db_v2"
df, _, _, _ = _load_hte_database_cached(str(db_path))

print(f"Columns: {df.columns.tolist()}")
if "Source_Group" in df.columns:
    print(f"Unique groups: {df['Source_Group'].unique()}")
    print(f"Sample values: {df['Source_Group'].head().tolist()}")
else:
    print("Source_Group column MISSING")
