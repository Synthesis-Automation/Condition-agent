from rxn_insight.classification import ReactionClassifier
s = "Cc1ccc(C)c(Br)c1.OB(O)c1ccccc1>>Cc1ccc(C)c(-c2ccccc2)c1"
rc = ReactionClassifier(s)
try:
    rc.classify_reaction()
    print('rxn_class:', rc.reaction)
    print('name_reaction:', rc.name_reaction())
    # attributes after classify
    print('reaction:', rc.reaction)
    print('is_cc_coupling:', rc.is_cc_coupling())
except Exception as e:
    print('error:', repr(e))
