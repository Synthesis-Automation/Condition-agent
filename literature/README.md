# Literature Folder

Drop chemistry papers, procedures, and lab notes here. ChemCoworker will
search them automatically via the `search_literature` tool.

## Supported Formats

| Format | Notes |
|--------|-------|
| `.txt` | Plain text — fastest to load |
| `.md`  | Markdown — headings help chunking |
| `.pdf` | Requires `pypdf` or `pdfminer.six`: `pip install pypdf` |

## Optional Metadata

Add a sidecar `{filename}.meta.json` alongside any document to provide
citation details that appear in search results:

```json
{
  "title": "Palladium-Catalyzed Suzuki Coupling of Hindered Aryl Chlorides",
  "doi": "10.1021/ja123456",
  "year": 2021
}
```

Without a `.meta.json`, the filename is used as the title.

## Example Files

```
literature/
  suzuki_hindered_substrates.txt          ← paste procedure from a paper
  suzuki_hindered_substrates.txt.meta.json
  buchwald_hartwig_amination_review.pdf
  my_lab_procedures.md                    ← private lab protocols
  snar_scope_2023.txt
```

## How It Works

- Documents are indexed automatically on first search (no setup step)
- Searches are keyword-based: query terms are matched against document chunks
- Results include an excerpt (~600 chars) and source citation
- The LLM combines literature findings with HTE database recommendations

## Configuration

To use a different folder, set the environment variable:
```
CHEM_DOCS_PATH=/path/to/your/docs
```

## Tips for Best Results

- Paste the **Experimental Section** of papers (procedures are dense with conditions)
- Include the **Supporting Information** — it often has substrate scope tables
- For lab procedures, use clear headings: `## Suzuki Coupling (General Procedure A)`
- Shorter, focused documents search better than entire textbooks
