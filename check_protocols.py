from chemtools.protocol.indexer import ProtocolIndexer

idx = ProtocolIndexer.load()
print(f'Total protocols: {len(idx.records)}')
print('\nProtocol families:')
families = {}
for fname, rec in idx.records.items():
    if rec.family not in families:
        families[rec.family] = []
    families[rec.family].append(fname)

for family, files in sorted(families.items()):
    print(f'  {family}: {len(files)} protocol(s)')
    for f in files[:2]:  # Show first 2
        print(f'    - {f}')
