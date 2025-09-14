import rxn_insight, importlib, inspect
m = importlib.import_module('rxn_insight.classification')
print('classification module attrs:', [a for a in dir(m) if not a.startswith('_')])
for name in ['analyze', 'classify', 'classify_reaction', 'analyze_reaction', 'predict']:
    obj = getattr(m, name, None)
    print(name, 'callable?', callable(obj), 'type', type(obj))
    if callable(obj):
        try:
            import json
            s = "Cc1ccc(C)c(Br)c1.OB(O)c1ccccc1>>Cc1ccc(C)c(-c2ccccc2)c1"
            out = obj(s)
            print(name, 'returned type', type(out))
            print(name, 'returned', out)
        except Exception as e:
            print(name, 'call error:', repr(e))
