import os
os.environ['CHEMTOOLS_DISABLE_RDKIT']='1'
from fastapi.testclient import TestClient
from app.main import app
client = TestClient(app)
resp = client.post('/api/v1/condition-core/validate-dataset', json={'path':'data/reaction_dataset/Ullman-C-N.jsonl','limit': 25})
print(resp.status_code)
print(resp.json())
