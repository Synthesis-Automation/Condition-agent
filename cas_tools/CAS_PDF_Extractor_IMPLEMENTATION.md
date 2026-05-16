CAS Number Extractor

Files added under /cas_tools:

- cas_number_extractor.py: core extraction logic for text, CSV/TSV, PDF, Word (.docx/.docm), and Excel (.xlsx/.xlsm/.xls)
- cas_number_extractor_gui.py: PyQt folder-based GUI matching the Scifinder RDF processor workflow

Behavior:

- Select one folder and scan all supported files recursively.
- Generic text files are accepted even when the extension is unusual, as long as the content decodes as text.
- Output is written to <selected-folder>/<folder_name>_cas_numbers.csv.
- CSV columns: source_file, relative_path, file_type, location, cas_number, context

Run:

```powershell
python .\cas_tools\cas_number_extractor_gui.py
```

Notes:

- PDF support requires pypdf or pdfminer.six.
- .docx/.docm use python-docx when available, otherwise a zip/xml fallback.
- .xlsx/.xlsm use openpyxl when available, otherwise a zip/xml fallback.
- Legacy .xls relies on pandas and a compatible Excel engine in the environment.
