import importlib
mdb = importlib.import_module('rxn_insight.database')
attrs = [a for a in dir(mdb) if not a.startswith('_')]
print('database attrs:', attrs)
print('has DB classes?', [x for x in attrs if 'Database' in x or 'db' in x.lower()])
try:
    from rxn_insight.database import Database
    print('Database class:', Database)
except Exception as e:
    print('Database import error:', repr(e))
