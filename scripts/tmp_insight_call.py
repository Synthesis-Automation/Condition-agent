import rxn_insight
s = "Cc1ccc(C)c(Br)c1.OB(O)c1ccccc1>>Cc1ccc(C)c(-c2ccccc2)c1"
print('module has classification:', hasattr(rxn_insight, 'classification'))
print('callable?', callable(getattr(rxn_insight,'classification', None)))
try:
    res = rxn_insight.classification(s)
    print('classification result type:', type(res))
    print('result:', res)
except Exception as e:
    print('classification call error:', repr(e))
