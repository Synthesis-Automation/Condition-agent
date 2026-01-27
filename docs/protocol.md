# the folder /protocol contain our old system of protocol-based recommendation. it used the json format database. I want to move to using of similar uniformed csv-based dataformat like the HTE recommender. is it feasible?

Yes, it’s feasible, but it isn’t a simple file‑extension swap. The protocol DB is hierarchical (reaction_setup, chemicals, conditions, workup/purification, etc.), while the HTE recommender CSVs are flat rows. You’ll need either a normalized CSV schema or a hybrid.

Why it’s non‑trivial

The protocol system reads JSON directly in indexer.py, recommend.py, validate_protocols.py, and matcher.py.
The recommender pulls full protocol details (e.g., reaction_setup, chemicals, conditions) from JSON to build results. A single flat CSV can’t represent those nested lists without losing structure.
Feasible paths

Hybrid (lowest risk): Keep JSON for full protocol details, add a protocols.csv for metadata + reaction fields (reaction_smiles, SMARTS, family, tags, source). Update the indexer to read CSV but still load JSON for full details. This keeps the recommendation output intact while giving you a uniform CSV for indexing.
Fully CSV (more work): Use multiple CSVs (normalized tables) and reconstruct protocol dicts in code.
protocols.csv: protocol_id, reaction_smiles, reaction_SMARTS (semicolon), family, tags, notes, source fields, etc.
protocol_chemicals.csv: protocol_id, setup_index, chem_index, name, role, amounts, etc.
protocol_conditions.csv: protocol_id, setup_index, cond_index, temperature_C, time_h, atmosphere, details
protocol_workup.csv, protocol_purification.csv, optional protocol_notes.csv
Then write a loader to rebuild the in‑memory JSON‑like structure for the existing API output.
If you want to proceed, tell me:

Do you need to preserve full procedure details (reaction_setup/workup/purification), or is metadata + core conditions enough?
Do you want hybrid (CSV for index, JSON for details) or fully CSV (normalized tables)?
Once you choose, I can draft the CSV schema and a converter, and update the loaders accordingly.

# my answer

1) in one csv, but could support multiple csv
2) original_procedure, reaction_setup_json, source (in one column)

# csv_utils.py is the unfinished work
