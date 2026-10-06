# RxnOptBench — reconstructed screening records

**Completed data, not just a reconstruction toolkit.** This package was generated from the user-uploaded `multiple_choice_all_varying.jsonl`. No outside dataset, original PDF, reference subset, DOI resolver or Crossref response was used.

## Start here

- **Human review:** open `RxnOptBench_screening_review.xlsx`.
- **Condition-recommendation / data-processing input:** use `rxnoptbench_screening_rows.jsonl` for native nested objects, or `rxnoptbench_screening_rows.csv` for a flat table with JSON-encoded object columns.
- **Quality and provenance:** read `QUALITY_REPORT.md` and `DATA_DICTIONARY.md` before filtering or training.

All text outputs use UTF-8. CSV files include a UTF-8 byte-order mark for Excel compatibility. JSONL and Markdown use UTF-8 without a byte-order mark. XLSX is a ZIP-based XML workbook, not a plain-text encoding.

## Reconstructed coverage

| Quantity | Result |
|---|---:|
| Retained screening rows | 4,768 |
| Released table records | 552 |
| Exact source identifier strings | 236 |
| Case-insensitive source reference groups | 232 |
| Sum of source-entry counts in table metadata | 5,451 |
| Source entries with no retained option | 683 |
| Fraction represented by retained options | 87.47% |
| Source condition keys preserved | 26 |
| Distinct canonical reaction groups | 358 |

The 5,451 figure is a **sum of metadata counts**, not an independently verified count of unique experiments. Repeated content and a case-colliding table identifier exist in the release. The 4,768 output rows represent all options in the supplied file, with none fabricated or silently removed.

## Files

| File | Contents |
|---|---|
| `RxnOptBench_screening_review.xlsx` | Six compact review sheets: Read me, Screening rows (17 essential fields for all 4,768 rows), Tables, References, QC issues, SMILES audit (the 58 invalid condition strings). Use CSV/JSONL for all fields and the complete SMILES audit. |
| `rxnoptbench_screening_rows.csv` / `.jsonl` | 4,768 screening rows; 139 fields including raw/canonical reaction SMILES, original conditions, outcomes, row pointers and flags. |
| `rxnoptbench_screening_tables.csv` / `.jsonl` | 552 table records. JSONL additionally embeds the complete corresponding original benchmark record. |
| `rxnoptbench_references.csv` / `.jsonl` | 232 case-normalized reference groups, all 236 raw identifier variants, inferred DOI links, per-reference counts. |
| `rxnoptbench_condition_items.csv` | 33,527 linked condition-item records, one source item per row, including original key, constant/varying origin, name and SMILES. Non-molecular settings such as temperature are included. |
| `rxnoptbench_reaction_molecules.csv` | 1,441 reaction-side molecule occurrences across the 552 table records; raw names/SMILES and validation results. |
| `rxnoptbench_smiles_validation.csv` | Audit of 1,772 distinct non-empty SMILES strings occurring across reaction-side and condition items. |
| `rxnoptbench_quality_issues.csv` | 724 information/warning records; not 724 distinct erroneous experiments. |
| `rxnoptbench_quality_report.json` | Machine-readable counts, source hash and checks performed. |
| `rxnoptbench_condition_key_map.json` | Original condition keys → flattened column suffixes. |
| `raw/multiple_choice_all_varying.jsonl` | Unmodified source file. |
| `reconstruct_rxnoptbench.py` | Offline script that regenerates the core CSV/JSONL/JSON outputs and source copy. |
| `DATA_DICTIONARY.md`, `QUALITY_REPORT.md`, `README.md` | Interpretation and schema documentation. |
| `manifest.json`, `SHA256SUMS.txt` | Package provenance and integrity hashes. |

## Reconstruction method

For each benchmark record, the script expands every observed `options[i]` and merges it with `input.constant_conditions`. It pairs this option only with the **same-index** `meta.scores[i]`, `meta.yields[i]`, `meta.effective_scores[i]`, `meta.option_relative_scores[i]` and `meta.option_source_entries[i]`. There is no Cartesian-product expansion and no copying of the best yield to other options.

Rows are presented in source-record order, then in increasing source `row_index`, rather than shuffled option order. `option_index` and a stable compact row ID remain available for exact backtracking. The full source record remains in table-level JSONL.

### What is derived, and what is preserved

**Preserved:** compound names and source SMILES; constant/varying condition dictionaries; yields/selectivity numbers; benchmark best-option labels; source IDs, row pointers and metadata; all 26 condition keys. Empty lists, empty SMILES and absent keys are distinguishable in the nested representations.

**Derived and labelled:** concatenated `reactants>>products`; RDKit-canonicalized SMILES; exact-text temperature/time conversions; canonical reaction grouping; DOI string decoding; case-insensitive reference grouping; duplicated-content signatures; completeness and validation flags.

**Not inferred:** missing base/additive roles, catalyst loading, reaction scale, concentrations, reaction class, original footnotes, analytical yield method, source inequality notation, article title/authors/date, reaction centers or atom mapping. Some quantities may be embedded in original condition text; the text is preserved, but no invented standalone quantity is supplied.

## Important interpretation

**SMILES:** all reaction-side strings parse with RDKit 2025.09.4. This does not establish chemical identity, atom balance, stereochemical correctness or correspondence to the original paper. There are 58 distinct non-empty condition strings that fail RDKit parsing, affecting 436 retained rows. Organometallic representation conventions can also cause parsing failures. Nothing is automatically repaired or deleted.

**Conditions:** the source key `reagents` is not recategorized into base/additive/oxidant/reductant. `Reaction_Type`, where present, is preserved as `reaction_type_source`; it is not an organic transformation-class assignment. A missing key, an empty list and an explicitly named absence are not equated.

**Temperature/time:** only a whole-string single numeric Celsius value or a whole-string single time duration is converted. `rt`/ambient remain numerically unspecified; ranges, inequalities, approximate values and staged protocols are not averaged, summed or silently reduced. The full text is always retained.

**Outcomes:** `yield`, `ee`, `dr`, `er` and `rr` are the released numeric values. No attempt is made to recover original formatting such as `trace`, `<5%`, `>20:1`, method qualifiers or uncertainty. `effective_score` and `relative_score` remain separate from yield. `is_benchmark_best` and `is_max_yield_retained` are different fields and need not mean the same thing.

**Reference information:** DOI strings are decoded using explicit publisher-specific patterns. In particular, Nature-style suffixes use hyphens rather than the periods used in several other publisher patterns. The result is marked `doi_inferred`, with `doi_resolver_verified=false`. All 232 reference-group strings are decoded, but links have not been externally verified. Titles/authors/journal/publication dates are not in this upload and have not been guessed. Internal page/table indices are not presented as printed page or table numbers. A missing `si`/`main` token is labelled `unspecified`, not assumed to be main text.

**Case variants and repeats:** 236 exact source IDs collapse to 232 under case-insensitive grouping. Two records (`T0464`, `T0495`) have parent-table IDs differing only in case but different contents (6 versus 11 retained options); both are retained and flagged. Thirty groups contain identical released reaction/condition/outcome content within a normalized source reference, accounting for 32 additional appearances. These remain in the exports and are not assumed to be independent measurements or extraction errors.

**Coverage gap:** 683 metadata-counted source entries have no retained option in this file, across 254 table records. Their conditions/outcomes are not reconstructed. The uploaded file alone does not establish why each individual entry was removed. Per-table absent source row indices are supplied only as provenance gaps.

## Validation

The script checked array alignment, unique row identifiers, unique in-range source row pointers, absence of constant/varying key conflicts, and full option-key coverage. It compared all 4,768 reconstructed records with the supplied structures, condition objects, score objects, source pointers, relative/effective scores and best-option labels. It also checked CSV and JSONL re-reading. The 552 embedded original table records match the input exactly as JSON objects.

This is **round-trip validation against the uploaded benchmark**, not independent validation against original publications or the authors' separate reference subset.

## Reproduce core data offline

Python 3.10 or newer is sufficient for the script. For the same structure-validation results, use the RDKit version recorded in `requirements_optional.txt`.

```bash
python reconstruct_rxnoptbench.py --input raw/multiple_choice_all_varying.jsonl --out regenerated
```

The script still extracts records without RDKit, but canonical structures and parsing results are then unavailable. It does not regenerate the optional Excel review workbook or the explanatory Markdown files. No internet connection or Hugging Face account is needed for the reconstruction itself.

## Sources and snapshot

Input: `multiple_choice_all_varying.jsonl`, 2,637,997 bytes.

SHA-256: `252cd2bca895fe0b0d2d1a1fd47da2ca69d478da53774cd8977f6a2d28853d58`

User-provided source context (not fetched for this reconstruction):

- https://huggingface.co/datasets/songjhPKU/RxnOptBench/blob/main/data/jsonl/multiple_choice_all_varying.jsonl
- https://arxiv.org/html/2610.02242v1

The raw data remain subject to the original source's applicable terms. This reconstruction does not assign a new license to third-party dataset content.
