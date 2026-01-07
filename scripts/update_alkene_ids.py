import json
fpath = r'chemtools/taxonomy/data/organic_compounds.v1.3.json'
with open(fpath, 'r', encoding='utf-8') as f:
    d = json.load(f)

for c in d['compounds']:
    # Update Substituent ID references
    if c.get('B') == '-CH=CH2':
        c['B'] = '-C*=CH2'
    if c.get('B') == '-CH=CR2':
        c['B'] = '-C*=C'
    
    # Update Names
    c['name'] = c['name'].replace('-CH=CH2', '-C*=CH2').replace('-CH=CR2', '-C*=C')
    
    # Update Descriptions
    if 'description' in c:
        c['description'] = c['description'].replace('-CH=CH2', '-C*=CH2').replace('-CH=CR2', '-C*=C')
        if '-C*=C' in c['description']:
             c['description'] = c['description'].replace('internal/substituted', 'general')

with open(fpath, 'w', encoding='utf-8') as f:
    json.dump(d, f, indent=2)

print("Update completed.")
