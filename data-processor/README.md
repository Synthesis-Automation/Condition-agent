# SciFinder RDF Processor

`data-processor/` has been reduced to the active pipeline only:

- `Scifinder_rdf_processer.py` (Qt GUI entrypoint)
- `process_reactions.py` (core RDF/TXT parsing helpers used by the GUI)

## Install

```powershell
pip install -r requirements.txt
```

Optional packages:

```powershell
pip install PyQt6
pip install rdkit
```

## Run GUI

```powershell
python .\data-processor\Scifinder_rdf_processer.py
```

## Run CLI Core (direct)

```powershell
python .\data-processor\process_reactions.py --rdf .\Reaction_2024-1.rdf --txt .\Reaction_2024-1.txt --out .\reactions_2024-1.csv
```

## Notes

- `Scifinder_rdf_processer.py` imports `parse_rdf` and `assemble_rows` from `process_reactions.py`.
- RDKit is optional. If unavailable, SMILES-derived fields stay blank.
- `xxhash` is optional. If unavailable, hash fields (e.g. `CondSig`, `FamSig`) are blank.
