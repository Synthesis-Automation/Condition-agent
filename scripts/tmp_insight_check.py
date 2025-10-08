import importlib
try:
    m = importlib.import_module('rxn_insight')
    print('MODULE', m.__name__)
    attrs = sorted([a for a in dir(m) if not a.startswith('_')])
    print('ATTRS', attrs)
    for name in ['analyze','analyze_reaction','classify','classify_reaction','analyze_smiles','classify_smiles']:
        fn = getattr(m, name, None)
        print(name, 'callable?', callable(fn))
    cls = getattr(m, 'Reaction', None)
    print('Has Reaction class?', cls is not None)
except Exception as e:
    print('IMPORT_ERROR', repr(e))
