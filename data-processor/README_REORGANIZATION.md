# RDF File Reorganization

## Overview

This directory contains tools to reorganize RDF files from a complex nested structure into a flat, organized structure by reaction type.

## Problem

The `original_dataset/` folder had a complex structure:

- Multiple levels of nesting (reaction types → subfolders → year folders → RDF files)
- Year-based subfolders (e.g., `2020-2024/`, `2023-2025/`)
- Duplicate filenames across different year folders (e.g., `1000.rdf` in both `2020-2024/` and `2023-2025/`)
- Inconsistent folder organization across reaction types

## Solution

The `reorganize_rdf_files.py` script flattens the structure to:

```
reorganized/
├── Amide_formation/
│   ├── 1.rdf
│   ├── 2.rdf
│   └── ...
├── C_N_Coupling/
│   ├── 1.rdf
│   ├── 2020-2024_1000.rdf    # Renamed to avoid conflict
│   ├── 2023-2025_1000.rdf    # Renamed to avoid conflict
│   └── ...
├── Suzuki/
│   └── ...
└── ...
```

## Features

✅ **Flat Structure**: All RDF files moved to `/reaction_type/*.rdf`
✅ **Smart Renaming**: Handles duplicate filenames by adding year/folder prefixes
✅ **Preservation**: Original files kept intact (copies, not moves)
✅ **Tracking**: Creates JSON mapping file for full traceability
✅ **Safe**: Dry-run mode to preview changes before execution

## Usage

### 1. Preview Changes (Dry Run)

```bash
cd data-processor
python reorganize_rdf_files.py --dry-run
```

This shows what will be done without making any changes.

### 2. Execute Reorganization

```bash
cd data-processor
python reorganize_rdf_files.py
```

You'll be prompted to confirm before files are copied.

### 3. Custom Paths

```bash
python reorganize_rdf_files.py --source my_data --dest output_folder
```

### 4. Skip Backup Mapping

```bash
python reorganize_rdf_files.py --no-backup
```

## Renaming Strategy

When duplicate filenames are detected, the script uses this priority:

1. **Original name** - If unique, keep as-is
2. **Year prefix** - Add year folder name (e.g., `2020-2024_1000.rdf`)
3. **Folder prefix** - Add parent subfolder name (e.g., `Chan-Lam_1000.rdf`)
4. **Combined prefix** - Year + folder (e.g., `2020_Chan-Lam_1000.rdf`)
5. **Sequential number** - Add counter (e.g., `1000_2.rdf`)

## Output Files

### Reorganized RDF Files

Location: `data-processor/reorganized/`

Flat structure with all RDF files organized by reaction type.

### Mapping File

Location: `data-processor/reorganization_mapping_YYYYMMDD_HHMMSS.json`

Contains complete tracking information:

- Timestamp of reorganization
- Source and destination paths
- Statistics per reaction type
- Full mapping of original → new paths
- Metadata (year folders, parent folders, relative paths)

Example mapping entry:

```json
{
  "C:\\...\\C_N_Coupling\\2020-2024\\1000.rdf": {
    "new_path": "C:\\...\\reorganized\\C_N_Coupling\\2020-2024_1000.rdf",
    "original_name": "1000.rdf",
    "new_name": "2020-2024_1000.rdf",
    "year_folder": "2020-2024",
    "parent_folder": null,
    "relative_path": "2020-2024\\1000.rdf"
  }
}
```

## Statistics (Current Dataset)

From the latest analysis:

- **Total Reaction Types**: 17
- **Total RDF Files**: 628
- **Files Renamed**: 94 (15%)
- **Files Preserved**: 534 (85%)

Top reaction types by file count:

1. Suzuki: 113 files (0 renamed)
2. Amide_formation: 86 files (0 renamed)
3. C_N_Coupling: 81 files (40 renamed)
4. C–C bond formation & cross-couplings: 68 files (28 renamed)
5. Multicomponent & heterocycle-forming: 35 files (0 renamed)

## Important Notes

⚠️ **Original Files Preserved**: The script creates **copies**, not moves. Original files remain untouched.

⚠️ **Check Before Deleting**: Review the reorganized folder and mapping file before deleting originals.

⚠️ **Year Folders**: Files from year subfolders are automatically prefixed with the year to avoid conflicts.

⚠️ **Special Characters**: Folder names are sanitized for safe filename usage (spaces → underscores, etc.).

## Troubleshooting

### "Source directory not found"

Ensure you're running from `data-processor/` or specify `--source` with correct path.

### Files not showing up

Check if files are in deeper subfolders - the script recursively scans all levels.

### Duplicate name conflicts

The script automatically handles this - check the mapping file to see how files were renamed.

## Example Workflow

```bash
# 1. Preview what will happen
python reorganize_rdf_files.py --dry-run | tee preview.txt

# 2. Review the preview
cat preview.txt

# 3. Execute if everything looks good
python reorganize_rdf_files.py

# 4. Verify results
ls -R reorganized/

# 5. Check mapping file
cat reorganization_mapping_*.json

# 6. After verification, optionally archive originals
# (DO NOT DELETE until you've confirmed everything works!)
```

## Future Improvements

Potential enhancements:

- [ ] Add option to move instead of copy
- [ ] Support for excluding specific folders
- [ ] Progress bar for large datasets
- [ ] Validation of RDF file integrity
- [ ] Automatic duplicate detection in reorganized folder

## Questions?

Contact the maintainer or check the script's `--help`:

```bash
python reorganize_rdf_files.py --help
```
