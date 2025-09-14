from rxn_insight import Reaction
s = "Cc1ccc(C)c(Br)c1.OB(O)c1ccccc1>>Cc1ccc(C)c(-c2ccccc2)c1"
r = Reaction(s)
print('constructed ok')
print('classification:', getattr(r,'classification', None))
print('rxn_class:', getattr(r, 'classification', None))
try:
    print('has analyze?', callable(getattr(r,'analyze', None)))
    print('has classify?', callable(getattr(r,'classify', None)))
except Exception as e:
    print('err', e)
