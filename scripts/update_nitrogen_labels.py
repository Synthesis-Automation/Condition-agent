import json
fpath = r'chemtools/taxonomy/data/organic_compounds.v1.3.json'
with open(fpath, 'r', encoding='utf-8') as f:
    d = json.load(f)

mapping = {
    'Hydrazine': 'NHN<',
    'Hydrazide': 'CONHNH<',
    'Urea': 'NHCON<',
    'Thiourea': 'NHCSN<',
    'Sulfamide': 'NHSO2N<'
}

for c in d['compounds']:
    for k, v in mapping.items():
        if k in c['name']:
            # Replace -K with -V
            c['name'] = c['name'].replace('-' + k, '-' + v)
            # Just in case some don't have the dash
            c['name'] = c['name'].replace(k, v)
    
    # Clean up any potential double dashes formed by accident
    c['name'] = c['name'].replace('--', '-')

with open(fpath, 'w', encoding='utf-8') as f:
    json.dump(d, f, indent=2)

print("Batch update completed.")
