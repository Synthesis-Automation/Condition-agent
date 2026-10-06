# Quality report — reconstructed RxnOptBench screening data

## Executive findings

**4,768 retained options** were reconstructed from **552 table records**, with no dropped or fabricated options. The input provides **5,451** total source-entry counts, leaving **683** entries that are not represented by a retained option. This covers **87.47%** of the metadata sum, not a proven percentage of distinct original experiments.

| Check | Result |
|---|---:|
| Structure/conditions/scores/provenance round-trip comparisons | 4,768 / 4,768 passed |
| Complete original table records retained | 552 / 552 |
| Tables with valid reaction-side SMILES syntax | 552 / 552 |
| Reaction-side source SMILES strings, distinct | 786 |
| Condition SMILES strings, non-empty and distinct | 1,025 |
| Condition SMILES strings failing RDKit parsing, distinct | 58 |
| Rows affected by a condition SMILES parsing failure | 436 |
| Table records affected by a condition SMILES parsing failure | 93 |
| Condition-item appearances containing a wildcard | 2 |
| Source case-variant groups | 4 |
| Case-colliding parent-table identifier groups | 1 (two different records) |
| Identical released-content groups retained | 30 |
| Additional appearances of identical released content | 32 |
| Table records with unretained source entries | 254 |

Parsing results use RDKit **2025.09.4**. A parsing failure is a representation warning, not proof that a reagent is chemically impossible. A parsing success is not proof that the source structure is the correct compound.

## Outcome availability

| Released metric | Non-null retained rows |
|---|---:|
| yield | 4,768 |
| ee | 149 |
| dr | 337 |
| er | 159 |
| rr | 152 |

Metric counts overlap: an entry can have more than one selectivity metric. There are **828** rows with a released numeric yield of zero and **1,479** rows with yield at or below 10. No zeros were imputed.

## Conservative physical-value parsing

A single explicit numeric Celsius temperature was extracted for **2,548** rows. A single explicit duration was converted to hours for **4,114** rows. The remaining values stay in their original textual form. Room temperature was not assigned an assumed numerical value; staged protocols were not collapsed.

## Reference and table identity

There are **236** literal `meta.doi` strings and **232** case-insensitive groups. The four sets of variants are:

- `10_1039_D4GC06293K` / `10_1039_d4gc06293k`
- `10_1039_D5GC03495G` / `10_1039_d5gc03495g`
- `10_1039_D5SC01073J` / `10_1039_d5sc01073j`
- `10_1039_D5SC05212B` / `10_1039_d5sc05212b`

The parent-table case collision is:

| Export record | Parent-table ID | Retained options | Action |
|---|---|---:|---|
| T0464 | `10_1039_D5GC03495G_si_13_table_0` | 6 | Preserved, flagged |
| T0495 | `10_1039_d5gc03495g_si_13_table_0` | 11 | Preserved, flagged |

Do not merge those table records solely by lowercasing their parent-table IDs. Their input contents differ. Internal indices are preserved without asserting printed page/table numbering.

Every DOI link is inferred from a source identifier, with publisher-specific suffix punctuation, and is explicitly unverified. No bibliographic title, author list, journal name or publication date was added.

## What the checks do not establish

The supplied source is already a processed benchmark. Original table footnotes, the original formatting of limits/ratios, analytical methods, complete experimental quantities and all discarded source rows may be absent. None can be recovered merely by unpacking the supplied JSONL.

No original paper/SI, external reference-subset file or DOI/Crossref record was used. No atom mapping, atom-balance verification, identity correction, catalyst-role inference or reaction classification was performed. Chemical quality therefore still requires source-level review before treating the output as a fully curated condition-recommendation training set.

## Audit files

`rxnoptbench_quality_issues.csv` contains per-table gaps, per-row condition parsing warnings, the two case-colliding records and repeated-content groups. It has 724 issue records; these are not 724 unique problematic experiments.

`rxnoptbench_smiles_validation.csv` includes each distinct non-empty SMILES string, associated source names/roles, canonicalization result, parse status and affected record IDs. `rxnoptbench_condition_items.csv` retains exact item-level placement, including missing strings.

`rxnoptbench_quality_report.json` contains all counts, validation scope and the raw-input hash.
