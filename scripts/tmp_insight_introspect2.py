from rxn_insight import Reaction
r = Reaction(r'')
print('constructed ok')
print('dir:', [a for a in dir(r) if not a.startswith('_')])
print('classification attr:', getattr(r,'classification', None))
try:
    print('analyze:', callable(getattr(r,'analyze', None)))
except Exception as e:
    print('analyze err', e)
