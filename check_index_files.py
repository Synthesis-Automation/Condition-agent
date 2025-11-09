import json

with open('build/unified_index_complete/index.json') as f:
    data = json.load(f)

print("Rule files in index vs actual files:")
print("=" * 60)

for rule in data['rules']:
    source_file = rule['source_file']
    id = rule['id']
    
    from pathlib import Path
    exists = Path(source_file).exists()
    status = "✅" if exists else "❌"
    
    print(f"{status} {id}")
    print(f"   Index: {source_file}")
    if not exists:
        # Try to find similar file
        dir_path = Path(source_file).parent
        if dir_path.exists():
            similar = list(dir_path.glob(f"*{id.split('_')[0]}*.json"))
            if similar:
                print(f"   Found: {similar[0]}")
    print()
