from rxn_insight.classification import ReactionClassifier
from rxn_insight.database import Database
s = "Cc1ccc(C)c(Br)c1.OB(O)c1ccccc1>>Cc1ccc(C)c(-c2ccccc2)c1"
rc = ReactionClassifier(s)
rc.classify_reaction()
db = Database()
try:
    name = rc.name_reaction(db.smirks_db)
    print('name:', name)
except Exception as e:
    print('name error:', repr(e))
print('is_cc_coupling:', rc.is_cc_coupling())
print('is_heteroatom_alkylation:', rc.is_heteroatom_alkylation())
