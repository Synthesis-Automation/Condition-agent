from rxn_insight import Reaction
import inspect
sig = None
try:
    print('Reaction:', Reaction)
    r = Reaction('CCBr.B(O)O>>')
    print('dir:', [a for a in dir(r) if not a.startswith('_')])
    for name in ['analyze','classify','analyze_reaction','classification','classify_reaction','analyze_smiles','classify_smiles']:
        attr = getattr(r, name, None)
        print(name, type(attr), 'callable?', callable(attr))
        if callable(attr):
            try:
                res = attr()
                print(name, 'result type:', type(res))
                print('result keys/attrs:', getattr(res,'__dict__', None) and list(getattr(res,'__dict__').keys()))
            except Exception as e:
                print(name, 'call error:', repr(e))
    # Try attribute access for classification
    try:
        print('classification value:', getattr(r, 'classification', None))
    except Exception as e:
        print('classification get error:', repr(e))
except Exception as e:
    print('ERROR constructing Reaction:', repr(e))
