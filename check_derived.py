import json

with open('chemtools/featurizers/calculable_features.json') as f:
    data = json.load(f)

shortcuts = data.get('derived_shortcuts', [])
print(f"Total derived_shortcuts: {len(shortcuts)}\n")

print("First 5 shortcuts:")
for i, s in enumerate(shortcuts[:5], 1):
    token = s.get('token', '')
    derive = s.get('derive', '')[:100]
    print(f"{i}. {token}")
    print(f"   Derive: {derive}...")
    print()
