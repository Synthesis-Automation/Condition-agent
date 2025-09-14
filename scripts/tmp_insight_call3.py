from rxn_insight.classification import ReactionClassifier
import inspect
rc = ReactionClassifier()
print('ReactionClassifier methods:', [m for m in dir(rc) if not m.startswith('_')])
try:
    s = "Cc1ccc(C)c(Br)c1.OB(O)c1ccccc1>>Cc1ccc(C)c(-c2ccccc2)c1"
    res = rc.tag_reaction(s)
    print('tag_reaction OK type:', type(res))
    print('res:', res)
except Exception as e:
    print('tag_reaction error:', repr(e))
