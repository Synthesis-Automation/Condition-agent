import requests
import json

# Test Fusion endpoint
response = requests.post(
    'http://localhost:8000/api/v1/recommend/fusion',
    json={
        'reaction': 'Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1',
        'k': 50,
        'max_variants': 3
    }
)

print("Status:", response.status_code)
print("\nJSON Output:")
print(json.dumps(response.json(), indent=2))
