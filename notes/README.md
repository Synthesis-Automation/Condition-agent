# Notes Folder — Chemistry Knowledge Base

This folder contains reaction-type-specific chemistry notes extracted from
literature sources. Think of it as the chemistry equivalent of `CLAUDE.md` —
guidelines and caveats that an expert chemist would write in the margins of
their notebook.

## How It Works

Notes are organized per reaction type:

```
notes/
  suzuki_miyaura.md        ← solvent warnings, catalyst notes, side reactions
  buchwald_hartwig.md
  general.md               ← cross-reaction insights
```

Each file is appended to automatically by the `intake` command and read at
query time by the `read_reaction_notes` tool.

## Adding Notes

### From a URL or PDF (recommended)

```bash
# Fetch a webpage and extract generalizable chemistry knowledge
python chem_coworker/cli.py intake https://www.orgsyn.org/demo.aspx?prep=v102p0086

# From a PDF
python chem_coworker/cli.py intake paper.pdf --reaction-type suzuki_miyaura

# From a local text file
python chem_coworker/cli.py intake my_lab_notes.txt
```

### Manually

Just edit or create a `.md` file directly. The agent reads markdown and
understands bullet-point chemistry notes without any indexing step.

## Note Format

The LLM extracts notes in this structure (you can also write manually):

```markdown
## Source: Org. Syn. 2024  ·  2025-01-15

### Solvent Notes
- ✓ THF/H₂O 3:1 — excellent for Pd couplings with K₂CO₃
- ✗ DMF — causes proto-deboronation of arylboronic acids under basic conditions

### Reagent and Catalyst Notes
- XPhos Pd G3 preferred for hindered or electron-rich aryl chlorides
- Avoid Pd(OAc)₂ without phosphine — rapid decomposition to Pd black

### Side Reactions
- Homocoupling of ArB(OH)₂ when O₂ not excluded or Pd > 3 mol%

### Critical Conditions
- Add base LAST — prevents premature boronate hydrolysis
- Degas rigorously — O₂ kills catalyst in <5 min
```

## Configuration

To use a different folder:

```
CHEM_NOTES_PATH=/path/to/your/notes
```

## Tips

- Short, principle-based bullets search better than verbatim procedures
- Include the source name at the end of each bullet: `(Molander, OrgSyn 2024)`
- The `general.md` file is a catch-all for cross-reaction insights
