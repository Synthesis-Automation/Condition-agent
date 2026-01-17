# Chem Assistant Quickstart (Featurization/Analysis)

## Install

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## CLI

```powershell
python -m chem_assistant.chemtools_cli
```

Commands:

- help
- tools
- clear
- verbose on/off
- quit

## Agent API

```python
from chem_assistant.chemtools_agent import ChemToolsAgent

agent = ChemToolsAgent()
print(agent.run("Featurize molecule: c1ccccc1O"))
```

## Tool families

- Analysis: normalize SMILES/reactions, taxonomy classification, reaction summaries
- Detection: reaction type detection and motif ids
- Featurizers: unified bundles, motif-based features, reaction-pair features
- Calculable: feature tokens for reactant types

## Notes

- RDKit is required for SMARTS-driven detection and RDKit properties.
