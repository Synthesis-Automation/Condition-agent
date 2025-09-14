from rxn_insight.classification import ReactionClassifier
s = "Cc1ccc(C)c(Br)c1.OB(O)c1ccccc1>>Cc1ccc(C)c(-c2ccccc2)c1"
rc = ReactionClassifier(s)
print('methods:', [m for m in dir(rc) if not m.startswith('_')])
try:
    res = rc.tag_reaction(s)
    print('tag_reaction type:', type(res))
    print('result:', res)
except Exception as e:
    print('error:', repr(e))
