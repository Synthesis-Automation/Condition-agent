# RDF Processor - Recursive Subfolder Support

## Changes Made

Updated the SciFinder RDF Processor to **recursively scan all subfolders** for RDF files and improved the output file naming convention.

### 1. Recursive RDF File Discovery

**Previous behavior:**
- Only scanned RDF files in the selected folder
- Ignored any subfolders

**New behavior:**
- Scans the selected folder **and all subfolders** recursively
- Uses `os.walk()` to traverse the entire directory tree
- Example: Selecting `C_N_Coupling/` will process RDF files in:
  - `C_N_Coupling/2020-2022/`
  - `C_N_Coupling/2023-2025/`
  - `C_N_Coupling/2020-2024/`
  - Any other subfolders

### 2. Improved File Naming Convention

**Output file naming:**
- JSONL: `data/reaction_dataset/{category}.jsonl`
- Markdown: `{selected_folder}/{category}.md`

**Examples:**
- Select `original_dataset/C_N_Coupling/` → `C_N_Coupling.jsonl`
- Select `original_dataset/Suzuki/` → `Suzuki.jsonl`
- Select `original_dataset/C_O_Coupling/` → `C_O_Coupling.jsonl`

The files automatically include data from all subfolders (recursive scanning).

### 3. Enhanced UI Display

**File list improvements:**
- Shows relative paths from the selected folder
- Displays subfolder structure clearly
- Example display:
  ```
  2020-2022/reaction_001.rdf
  2020-2022/reaction_002.rdf
  2023-2025/reaction_003.rdf
  2023-2025/reaction_004.rdf
  ```

**Updated labels:**
- Added hint: "(includes all subfolders)" next to folder selection
- Updated file count: "Found X RDF files (including subfolders)"
- Improved info notes showing the naming pattern

### 4. Code Changes

**Modified functions:**

1. **`_find_rdf_files()`** - RDFWorker class
   ```python
   # Old: Only scanned current folder
   for file in os.listdir(self.folder_path):
       if file.lower().endswith('.rdf'):
           ...
   
   # New: Recursively scans all subfolders
   for root, dirs, files in os.walk(self.folder_path):
       for file in files:
           if file.lower().endswith('.rdf'):
               ...
   ```

2. **`_update_file_list()`** - RDFProcessorWindow class
   ```python
   # Updated to use os.walk() and show relative paths
   for root, dirs, files in os.walk(folder_path):
       for file in files:
           if file.lower().endswith('.rdf'):
               rel_path = os.path.relpath(full_path, folder_path)
               self.file_list.addItem(rel_path)
   ```

3. **Output path calculation** - `run_processing()` method
   ```python
   # Improved naming logic:
   # - Detects category from folder hierarchy
   # - Uses simple category name (e.g., C_N_Coupling)
   # - Handles paths like original_dataset/C_N_Coupling/2020-2022/
   # - Output files automatically replace old versions
   ```

## Usage

### Typical Workflow

1. **Select a category folder** (e.g., `C_N_Coupling/`, `Suzuki/`, `C_O_Coupling/`)
2. **All RDF files in subfolders are automatically detected**
3. **Click "Process RDF Files"**
4. **Output files are created:**
   - `data/reaction_dataset/C_N_Coupling.jsonl` (for chemtools)
   - `original_dataset/C_N_Coupling/C_N_Coupling.md` (for records)

### Example Directory Structure

```
data-processor/original_dataset/
├── C_N_Coupling/
│   ├── 2020-2022/
│   │   ├── reaction_001.rdf
│   │   └── reaction_002.rdf
│   ├── 2023-2025/
│   │   ├── reaction_003.rdf
│   │   └── reaction_004.rdf
│   └── 2020-2024/
│       └── reaction_005.rdf
├── Suzuki/
│   ├── 2018-2020/
│   └── 2021-2023/
└── C_O_Coupling/
    └── 2019-2024/
```

**Selecting `C_N_Coupling/` will process:**
- All 5 RDF files from all three subfolders
- Output: `C_N_Coupling.jsonl` with all reactions

## Benefits

1. **Simplified Workflow**: No need to process each subfolder separately
2. **Chemtools Compatible**: Output files use standard naming convention (no `_combined` suffix)
3. **Better Organization**: Category-level files instead of per-subfolder files
4. **Transparency**: UI shows exactly which files from which subfolders are included
5. **Backward Compatible**: Still works with folders that have no subfolders
6. **Automatic Replacement**: New processing automatically replaces old dataset files

## Testing Recommendations

1. Test with a category folder containing multiple subfolders
2. Verify all RDF files are found and displayed
3. Check output file naming follows the pattern
4. Confirm JSONL is saved to `data/reaction_dataset/`
5. Confirm Markdown is saved to the selected folder

---

**Status**: Implementation complete and ready for testing
**Files Modified**: `data-processor/Scifinder_rdf_processer.py`
