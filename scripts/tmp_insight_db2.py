from rxn_insight.database import Database
import inspect
print(Database)
print([m for m in dir(Database) if not m.startswith('_')])
db = Database()
print([m for m in dir(db) if not m.startswith('_')])
for name in ['smirks_db','get_reaction_template','analyze_reactions','curate_smirks']:
    print(name, 'present?', hasattr(db,name), getattr(db,name, None))
