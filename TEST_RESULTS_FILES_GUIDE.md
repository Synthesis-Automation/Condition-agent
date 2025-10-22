# Test Results Output Files Summary

## 📁 Generated Files Overview

After running the comprehensive cross-family recommendation test on all 320 sample reactions, the following output files have been generated:

---

## 1. **test_all_sample_reactions_results.json** 
*Raw data (machine-readable)*

- **Format:** JSON
- **Content:** Complete raw test results including all reaction details, recommendations, and metadata
- **Use case:** For programmatic analysis, further processing, or integration with other tools
- **Size:** Contains full nested data structures with all recommendation details

---

## 2. **test_results_concise.txt** ⭐ *RECOMMENDED FOR QUICK SCANNING*
*Best for quick reference*

- **Format:** Plain text, highly readable
- **Content:** One-line summary per reaction showing:
  - Reaction number and description
  - SMILES
  - **Top recommended condition** (Catalyst/Ligand/Base/Solvent/Temp/Time)
  - Precedent count in brackets [N precedents]
- **Example:**
  ```
  1. Suzuki - Simple Ph-Ph
     Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
     → Base: Potassium carbonate | Solvent: Ethanol [2 precedent(s)]
  ```
- **Use case:** Quick lookup of conditions for specific reactions

---

## 3. **test_results_detailed.txt**
*Comprehensive view with all options*

- **Format:** Plain text with structured sections
- **Content:** Complete details for each reaction including:
  - Full SMILES
  - Detected reaction family
  - **ALL recommended conditions** (not just top one)
  - Each condition broken down as: `Catalyst | Ligand | Base | Solvent | Temp | Time [N precedents]`
- **Example:**
  ```
  [1] Suzuki - Simple Ph-Ph
  ────────────────────────────────
  SMILES: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
  Detected Family: Unknown

  RECOMMENDED CONDITIONS (7 options):

    Option 1: Base: Potassium carbonate | Solvent: Ethanol [2 precedent(s)]
    Option 2: Base: Sodium hydrogen carbonate | Solvent: Water [1 precedent(s)]
    Option 3: Base: Potassium carbonate | Solvent: Water [1 precedent(s)]
    ...
  ```
- **Use case:** When you need all alternative conditions for a specific reaction

---

## 4. **test_results_readable.txt**
*Alternative readable format*

- **Format:** Plain text with different layout
- **Content:** Similar to detailed but with slightly different formatting
- **Use case:** Alternative view if you prefer a different layout style

---

## 5. **test_results_summary.csv** 📊
*Spreadsheet-compatible*

- **Format:** CSV (comma-separated values)
- **Columns:**
  - Index
  - Description
  - SMILES
  - Reaction Type (C-C Coupling, C-N Coupling, etc.)
  - Reaction Family (detected)
  - Number of Recommendations
  - Top Condition (simplified string)
- **Example:**
  ```csv
  Index,Description,SMILES,Reaction Type,Reaction Family,Num Recommendations,Top Condition
  1,"Suzuki - Simple Ph-Ph","Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",C-C Coupling,Unknown,7,"Potassium carbonate / Ethanol"
  ```
- **Use case:** Import into Excel, Google Sheets, or data analysis tools for filtering, sorting, and pivot tables

---

## 6. **test_results_condensed.md**
*Markdown format with top reactions*

- **Format:** Markdown
- **Content:** Top 5 reactions from top 10 reaction types with top 3 conditions each
- **Use case:** GitHub README, documentation, or markdown-compatible viewers

---

## 📊 Quick Stats from Test Results

- **Total Reactions:** 320
- **Success Rate:** 100%
- **Total Recommendations:** 1,550
- **Average Recommendations per Reaction:** 4.8
- **Execution Time:** ~61 minutes (~11.5s per reaction)

---

## 🎯 Recommended Workflow

### For Quick Lookup:
1. **Open:** `test_results_concise.txt`
2. **Search:** Ctrl+F for your reaction type or description
3. **Read:** Top condition directly shown

### For Detailed Analysis:
1. **Open:** `test_results_detailed.txt`
2. **Find your reaction**
3. **Review all recommended condition options**
4. **Choose based on your laboratory constraints**

### For Data Analysis:
1. **Open:** `test_results_summary.csv` in Excel/Sheets
2. **Filter/Sort** by reaction type, number of recommendations
3. **Create pivot tables** for statistical analysis

---

## 💡 Understanding the Condition Format

Conditions are shown in this format:
```
Catalyst: [name] | Ligand: [name] | Base: [name] | Solvent: [name] | Temp: [value] | Time: [value] [N precedent(s)]
```

**Only present fields are shown.** For example:
- `Base: K2CO3 | Solvent: THF [5 precedents]` - Only base and solvent specified
- `Catalyst: Pd(OAc)2 | Ligand: XPhos | Base: Cs2CO3 | Solvent: Toluene | Temp: 110°C` - Full condition set

**Precedent count** `[N precedent(s)]` indicates:
- How many similar reactions in the database used this exact condition combination
- Higher count = more validated/reliable conditions
- Conditions without precedent count are extrapolated from cross-family matching

---

## 🔧 Regenerating Results

To regenerate these files from the JSON:

```bash
# Generate all readable formats
python generate_detailed_readable.py

# Generate CSV and markdown
python generate_readable_results.py
```

To rerun the entire test:
```bash
python test_all_sample_reactions.py
```

---

## 📖 Example Use Cases

### Use Case 1: "I need conditions for Buchwald-Hartwig coupling"
1. Open `test_results_concise.txt`
2. Search for "Buchwald-Hartwig" or "C-N"
3. Find your specific substrate match
4. Use the top recommended condition

### Use Case 2: "I want to compare all C-C coupling conditions"
1. Open `test_results_summary.csv` in Excel
2. Filter "Reaction Type" column for "C-C Coupling"
3. Sort by "Num Recommendations" (descending) to see best-covered reactions
4. Analyze patterns in "Top Condition" column

### Use Case 3: "I need alternatives if my first choice doesn't work"
1. Open `test_results_detailed.txt`
2. Find your exact reaction
3. Try Option 1, if it fails, try Option 2, etc.
4. Options are ranked by precedent count and similarity

---

## ✅ Files Generated

- ✅ `test_all_sample_reactions_results.json` - Raw JSON data
- ✅ `test_results_concise.txt` - **Quick reference (RECOMMENDED)**
- ✅ `test_results_detailed.txt` - Complete conditions
- ✅ `test_results_readable.txt` - Alternative format
- ✅ `test_results_summary.csv` - Spreadsheet format
- ✅ `test_results_condensed.md` - Markdown summary

All files are in the root directory: `c:\Git-softwares\Condition-agent\`
