# Knowledge Base — Chemistry Knowledge Base

This folder is the two-layer knowledge store for ChemCoworker. Think of it as the
chemistry equivalent of `CLAUDE.md` — curated principles an expert chemist would write
in the margins of their notebook, backed by the original source documents for full-detail
lookups when needed.

## Structure

```
knowledge_base/
  notes/
    reactions/       — per-reaction-type notes (conditions, warnings, HTE data)
    mechanisms/      — mechanistic principle notes
    substrates/      — substrate class notes
    protocols/       — practical technique / procedure notes
    _index.json      — auto-generated; run build_index.py to rebuild
    build_index.py   — scans subfolders, parses front matter, writes _index.json
  sources/           — raw source files (HTML→text, PDFs, text files)
```

## How It Works

**Layer 1 — Notes** (`notes/`): Distilled, principle-based knowledge extracted by the
LLM. Fast to search, always loaded at query time. Full experimental procedures are
intentionally omitted here — only generalizable insights are kept.

**Layer 2 — Sources** (`sources/`): The original raw text saved during intake. Web pages
are saved as `.txt`; PDFs are saved as-is. Each note's front matter includes a
`source_file:` field pointing to the corresponding source file.

The agent can retrieve full source text on demand using the `read_literature_source` tool,
which looks up by filename or the original URL via `.meta.json` sidecars.

## Adding Notes

### Single URL or PDF

```bash
# Fetch a webpage and extract chemistry knowledge (reactions subfolder, default)
python chem_coworker/cli.py intake https://www.orgsyn.org/demo.aspx?prep=v92p0195

# Specify a note type
python chem_coworker/cli.py intake paper.pdf --note-type mechanisms

# Use a cheaper model for extraction
python chem_coworker/cli.py intake paper.pdf --extract-model deepseek-v3.2 --extract-provider aliyun
```

### Batch intake from a CSV file

`knowledge_base/intake_from_csv.py` reads URLs from a CSV file and processes each one.

**CSV formats supported:**

Plain list — one URL per line (e.g. `orgsyn_v102.csv`):

```text
https://orgsyn.org/demo.aspx?prep=v102p0001
https://orgsyn.org/demo.aspx?prep=v102p0019
```

CSV with header — URL column is auto-detected (`url`, `link`, `href`, etc.):

```csv
url,title
https://example.com/paper1,Suzuki Review
```

**Examples:**

```bash
# Basic run — processes all URLs in the CSV
python knowledge_base/intake_from_csv.py knowledge_base/orgsyn_v102.csv

# Use a specific model and provider
python knowledge_base/intake_from_csv.py knowledge_base/orgsyn_v102.csv \
    --model deepseek-v3.2 --provider aliyun

# Preview URLs without fetching anything
python knowledge_base/intake_from_csv.py knowledge_base/orgsyn_v102.csv --dry-run

# Resume a previously interrupted run (skips already-succeeded URLs)
python knowledge_base/intake_from_csv.py knowledge_base/orgsyn_v102.csv --resume

# Slower, more polite crawl
python knowledge_base/intake_from_csv.py knowledge_base/orgsyn_v102.csv --delay 5

# Extract as mechanism notes instead of reaction notes
python knowledge_base/intake_from_csv.py knowledge_base/orgsyn_v102.csv --note-type mechanisms

# CSV where the URL is in a non-standard column
python knowledge_base/intake_from_csv.py my_papers.csv --url-column doi_url
```

A log file is written automatically alongside the CSV as `{name}_log.jsonl`, recording
success/failure, detected reaction types, and timing for each URL. Use `--resume` to
skip URLs that already appear in the log as successful.

### Manually

Edit or create a `.md` file directly in the appropriate subfolder. See `notes/SCHEMA.md`
for the canonical front matter format.

## Note Format

Each note file has YAML front matter followed by source sections:

```markdown
---
id: suzuki_miyaura
type: reaction
title: Suzuki-Miyaura Cross-Coupling
bond_formed: [C-C]
metal: [palladium]
tags: [pd, boronic_acid, cross_coupling]
source_file: demo.aspx.txt
url: https://www.orgsyn.org/demo.aspx?prep=v102p0086
date: 2025-01-15
---

## Source: Org. Syn. 2024  ·  2025-01-15

### Solvent Notes
- THF/H₂O 3:1 — excellent for Pd couplings with K₂CO₃
- DMF — causes proto-deboronation of arylboronic acids under basic conditions

### Reagent and Catalyst Notes
- XPhos Pd G3 preferred for hindered or electron-rich aryl chlorides

### Side Reactions
- Homocoupling of ArB(OH)₂ when O₂ not excluded or Pd > 3 mol%

Full procedure available in source file.
```

## Index

`notes/_index.json` enables fast faceted search by type, bond, metal, and tags without
reading every file. Rebuild it after adding notes manually:

```bash
python knowledge_base/notes/build_index.py
```

The `intake` command rebuilds the index automatically after each extraction.

## Configuration

Override the default folders via environment variables:

```env
CHEM_NOTES_PATH=/path/to/your/notes      # overrides knowledge_base/notes/
CHEM_DOCS_PATH=/path/to/your/sources     # overrides knowledge_base/sources/
```
