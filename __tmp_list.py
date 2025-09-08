from chemtools.registry import search
print('BASE sample:', len(search(role='BASE', limit=5)), search(role='BASE', limit=5))
print('LIGAND sample:', len(search(role='LIGAND', limit=5)), search(role='LIGAND', limit=5))
