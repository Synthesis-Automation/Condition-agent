# Data dictionary

## Conventions

`*_json` columns are native arrays/objects in JSONL and JSON-encoded strings in CSV. Empty CSV numeric cells represent JSON null; the JSONL is authoritative for distinctions between null, empty strings, absent keys and empty arrays. All raw condition objects are preserved.

The primary row unit is **one retained benchmark option**, not necessarily one unique original experiment. Best-option labels are preserved separately from maximum yield. No original source entries absent from the options are invented.

## Screening-row fields

| Field | Meaning / provenance |
|---|---|
| `screening_row_id` | Derived stable identifier: input table ordinal plus original zero-based benchmark option index, e.g. T0001_O0004. |
| `table_record_id` | Derived table-record identifier, T0001–T0552, in uploaded record order. |
| `parent_table_id` | Unmodified meta.parent_table_id; preserve case. Not a printed table label. |
| `benchmark_question_id` | Unmodified benchmark id; links the reconstructed row to the original record. |
| `source_record_number` | One-based ordinal of the nonblank input record. |
| `reference_id` | Case-folded meta.doi source identifier; used only for reference grouping. |
| `doi_source_id` | Literal source meta.doi identifier, including its original case and underscore separators. |
| `doi_inferred` | DOI string decoded using a labelled publisher-specific rule; not independently verified. |
| `doi_url_inferred` | Resolver URL formed from doi_inferred; availability has not been checked. |
| `doi_decode_rule` | Named decoding rule used for the source DOI identifier. |
| `doi_resolver_verified` | False throughout: external DOI resolution was not performed. |
| `source_file` | Unmodified meta.source_file pointer to an upstream JSON file; that file is not included unless it is the uploaded input. |
| `source_scope_token` | si/main token parsed from parent_table_id, otherwise unspecified. Does not assert publication section when absent. |
| `source_page_index` | Integer decoded from the parent-table identifier; not asserted to be a printed page number. |
| `source_table_index` | Integer decoded from the parent-table identifier; not asserted to be a printed table number. |
| `reaction_smiles` | Reactant source SMILES joined by dots, then >>, then product source SMILES; no conditions inserted and no atom mapping added. |
| `reaction_smiles_canonical` | RDKit-canonicalized reaction-side SMILES, preserving molecule-list order and available stereochemistry; no atom mapping added. |
| `reactant_smiles` | Reactant source SMILES joined by dots in input order. |
| `product_smiles` | Product source SMILES joined by dots in input order. |
| `reactant_smiles_canonical` | Canonicalized individual reactant SMILES joined by dots in input order. |
| `product_smiles_canonical` | Canonicalized individual product SMILES joined by dots in input order. |
| `reaction_group_id` | Truncated SHA-256 of sorted canonical reactant/product molecule lists; an exact representation grouping, not a reaction class. |
| `reactant_names` | Source reactant content/name strings joined with &#124; for display. |
| `product_names` | Source product content/name strings joined with &#124; for display. |
| `reaction_smiles_valid` | True only if all supplied non-empty reaction-side strings parse as molecules with the recorded RDKit version; not chemical correctness. |
| `reaction_has_dummy_atoms` | Whether any reaction-side atom is a wildcard/dummy atom. |
| `main_product_smiles_source` | Unmodified meta.main_product_smiles. |
| `option_index` | Original zero-based index in the shuffled benchmark options array. |
| `row_index` | Unmodified zero-based source row_index from option_source_entries; output is sorted by this within a table. |
| `entry_idx` | Unmodified original entry label as a string, not converted to an integer. |
| `rxn_id` | Unmodified source-entry rxn_id; all values in this upload are the string 0. |
| `catalyst` | Convenience display of content strings from the exact source key `catalysts`, joined with &#124;; full objects remain in conditions_json. |
| `catalyst_smiles` | Non-empty original SMILES strings from `catalysts`, joined with &#124;; no molecular validation or repair is implied by this column. |
| `ligand` | Convenience display of content strings from the exact source key `ligand`, joined with &#124;; full objects remain in conditions_json. |
| `ligand_smiles` | Non-empty original SMILES strings from `ligand`, joined with &#124;; no molecular validation or repair is implied by this column. |
| `reagents` | Convenience display of content strings from the exact source key `reagents`, joined with &#124;; full objects remain in conditions_json. |
| `reagents_smiles` | Non-empty original SMILES strings from `reagents`, joined with &#124;; no molecular validation or repair is implied by this column. |
| `solvents` | Convenience display of content strings from the exact source key `solvents`, joined with &#124;; full objects remain in conditions_json. |
| `solvents_smiles` | Non-empty original SMILES strings from `solvents`, joined with &#124;; no molecular validation or repair is implied by this column. |
| `temperature` | Convenience display of content strings from the exact source key `reaction_temperature`, joined with &#124;; full objects remain in conditions_json. |
| `time` | Convenience display of content strings from the exact source key `reaction_time`, joined with &#124;; full objects remain in conditions_json. |
| `atmosphere` | Convenience display of content strings from the exact source key `atmosphere`, joined with &#124;; full objects remain in conditions_json. |
| `light_condition` | Convenience display of content strings from the exact source key `light_condition`, joined with &#124;; full objects remain in conditions_json. |
| `electrolyte` | Convenience display of content strings from the exact source key `Electrolyte`, joined with &#124;; full objects remain in conditions_json. |
| `anode` | Convenience display of content strings from the exact source key `Anode(+)`, joined with &#124;; full objects remain in conditions_json. |
| `cathode` | Convenience display of content strings from the exact source key `Cathode(-)`, joined with &#124;; full objects remain in conditions_json. |
| `current` | Convenience display of content strings from the exact source key `Constant_Current`, joined with &#124;; full objects remain in conditions_json. |
| `current_density` | Convenience display of content strings from the exact source key `Current_Density`, joined with &#124;; full objects remain in conditions_json. |
| `cell_voltage` | Convenience display of content strings from the exact source key `Constant_Cell_Voltage`, joined with &#124;; full objects remain in conditions_json. |
| `potential` | Convenience display of content strings from the exact source key `Constant_Potential`, joined with &#124;; full objects remain in conditions_json. |
| `charge_quantity` | Convenience display of content strings from the exact source key `Constant_Quantity`, joined with &#124;; full objects remain in conditions_json. |
| `mode` | Convenience display of content strings from the exact source key `Main_Mode`, joined with &#124;; full objects remain in conditions_json. |
| `reaction_type_source` | Convenience display of content strings from the exact source key `Reaction_Type`, joined with &#124;; full objects remain in conditions_json. |
| `stirring_speed` | Convenience display of content strings from the exact source key `speed`, joined with &#124;; full objects remain in conditions_json. |
| `pressure` | Convenience display of content strings from the exact source key `pressure`, joined with &#124;; full objects remain in conditions_json. |
| `pH` | Convenience display of content strings from the exact source key `PH`, joined with &#124;; full objects remain in conditions_json. |
| `temperature_C` | Conservative derived numeric Celsius temperature only for a whole-string single explicit °C value. No rt/ambient assumption. |
| `temperature_parse_status` | Whether the original text was a single explicit Celsius value, ambient without a numeric value, missing, or not safely reducible. |
| `time_h` | Conservative derived duration in hours for a single whole-string numeric duration; multi-stage durations are not summed. |
| `time_parse_status` | Whether a single duration was converted, the value was missing/empty, or the text was not safely reducible. |
| `yield` | Unmodified released numeric yield. Source values are in the 0–100 range; workbook displays percentage points without rescaling. |
| `ee` | Unmodified released numeric ee; no conversion or sign correction. |
| `dr` | Unmodified released numeric dr; original ratio/inequality notation is not reconstructed. |
| `er` | Unmodified released numeric er; original ratio/inequality notation is not reconstructed. |
| `rr` | Unmodified released numeric rr; original ratio/inequality notation is not reconstructed. |
| `effective_score` | Source meta.effective_scores[option_index]; not substituted for yield and not recalculated. |
| `relative_score` | Source meta.option_relative_scores[option_index]; not a universal yield or probability. |
| `is_benchmark_best` | True when option_index occurs in the original answer array; ties retained. |
| `is_max_yield_retained` | Derived exact equality to the maximum released numeric yield among this table record’s retained options. |
| `is_max_effective_score_retained` | Derived exact equality to the maximum released effective score among retained options. |
| `is_global_best_entry_source` | Whether entry_idx matches an original meta.global_best_entry_indices value. |
| `varying_keys` | Unmodified list of condition keys varied by the table’s benchmark question. |
| `condition_keys` | Sorted keys present in merged conditions; includes explicit empty-list keys. |
| `empty_condition_keys` | Present condition keys whose source value is an explicit empty list. |
| `missing_condition_keys` | Global condition keys absent from this row; absence is not treated as an explicit negative. |
| `invalid_condition_smiles_count` | Number of non-empty condition SMILES items failing RDKit parsing in this row. |
| `empty_condition_smiles_item_count` | Number of condition items with empty SMILES; physical settings may legitimately have no molecular SMILES. |
| `quality_flags` | Machine-readable warnings/information; rows are retained regardless of flag. |
| `conditions_json` | Native JSON object in JSONL, JSON-encoded cell in CSV. Exact constant + varying condition merge. |
| `constant_conditions_json` | Exact input.constant_conditions object. |
| `varying_conditions_json` | Exact options[option_index] object. |
| `reactants_json` | Exact input.reactants array, preserving name/SMILES alignment. |
| `products_json` | Exact input.products array, preserving name/SMILES alignment. |
| `scores_json` | Exact meta.scores[option_index] object. |
| `source_entry_json` | Exact meta.option_source_entries[option_index] object. |
| `cond__alternating_constant_forward_v` | Display content from exact source condition key `Alternating_Constant_Forward_V`. Original lists/objects remain in conditions_json. |
| `cond_smiles__alternating_constant_forward_v` | Original non-empty SMILES from exact source condition key `Alternating_Constant_Forward_V`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__alternating_constant_frequency` | Display content from exact source condition key `Alternating_Constant_Frequency`. Original lists/objects remain in conditions_json. |
| `cond_smiles__alternating_constant_frequency` | Original non-empty SMILES from exact source condition key `Alternating_Constant_Frequency`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__anode_plus` | Display content from exact source condition key `Anode(+)`. Original lists/objects remain in conditions_json. |
| `cond_smiles__anode_plus` | Original non-empty SMILES from exact source condition key `Anode(+)`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__cathode_minus` | Display content from exact source condition key `Cathode(-)`. Original lists/objects remain in conditions_json. |
| `cond_smiles__cathode_minus` | Original non-empty SMILES from exact source condition key `Cathode(-)`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__constant_cell_voltage` | Display content from exact source condition key `Constant_Cell_Voltage`. Original lists/objects remain in conditions_json. |
| `cond_smiles__constant_cell_voltage` | Original non-empty SMILES from exact source condition key `Constant_Cell_Voltage`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__constant_current` | Display content from exact source condition key `Constant_Current`. Original lists/objects remain in conditions_json. |
| `cond_smiles__constant_current` | Original non-empty SMILES from exact source condition key `Constant_Current`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__constant_potential` | Display content from exact source condition key `Constant_Potential`. Original lists/objects remain in conditions_json. |
| `cond_smiles__constant_potential` | Original non-empty SMILES from exact source condition key `Constant_Potential`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__constant_quantity` | Display content from exact source condition key `Constant_Quantity`. Original lists/objects remain in conditions_json. |
| `cond_smiles__constant_quantity` | Original non-empty SMILES from exact source condition key `Constant_Quantity`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__current_density` | Display content from exact source condition key `Current_Density`. Original lists/objects remain in conditions_json. |
| `cond_smiles__current_density` | Original non-empty SMILES from exact source condition key `Current_Density`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__electrolyte` | Display content from exact source condition key `Electrolyte`. Original lists/objects remain in conditions_json. |
| `cond_smiles__electrolyte` | Original non-empty SMILES from exact source condition key `Electrolyte`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__main_mode` | Display content from exact source condition key `Main_Mode`. Original lists/objects remain in conditions_json. |
| `cond_smiles__main_mode` | Original non-empty SMILES from exact source condition key `Main_Mode`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__ph` | Display content from exact source condition key `PH`. Original lists/objects remain in conditions_json. |
| `cond_smiles__ph` | Original non-empty SMILES from exact source condition key `PH`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__reaction_type` | Display content from exact source condition key `Reaction_Type`. Original lists/objects remain in conditions_json. |
| `cond_smiles__reaction_type` | Original non-empty SMILES from exact source condition key `Reaction_Type`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__atmosphere` | Display content from exact source condition key `atmosphere`. Original lists/objects remain in conditions_json. |
| `cond_smiles__atmosphere` | Original non-empty SMILES from exact source condition key `atmosphere`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__catalysts` | Display content from exact source condition key `catalysts`. Original lists/objects remain in conditions_json. |
| `cond_smiles__catalysts` | Original non-empty SMILES from exact source condition key `catalysts`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__cooling_heating_condition` | Display content from exact source condition key `cooling_heating_condition`. Original lists/objects remain in conditions_json. |
| `cond_smiles__cooling_heating_condition` | Original non-empty SMILES from exact source condition key `cooling_heating_condition`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__is_variable_current` | Display content from exact source condition key `is_variable_current`. Original lists/objects remain in conditions_json. |
| `cond_smiles__is_variable_current` | Original non-empty SMILES from exact source condition key `is_variable_current`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__ligand` | Display content from exact source condition key `ligand`. Original lists/objects remain in conditions_json. |
| `cond_smiles__ligand` | Original non-empty SMILES from exact source condition key `ligand`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__light_condition` | Display content from exact source condition key `light_condition`. Original lists/objects remain in conditions_json. |
| `cond_smiles__light_condition` | Original non-empty SMILES from exact source condition key `light_condition`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__pressure` | Display content from exact source condition key `pressure`. Original lists/objects remain in conditions_json. |
| `cond_smiles__pressure` | Original non-empty SMILES from exact source condition key `pressure`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__reaction_temperature` | Display content from exact source condition key `reaction_temperature`. Original lists/objects remain in conditions_json. |
| `cond_smiles__reaction_temperature` | Original non-empty SMILES from exact source condition key `reaction_temperature`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__reaction_time` | Display content from exact source condition key `reaction_time`. Original lists/objects remain in conditions_json. |
| `cond_smiles__reaction_time` | Original non-empty SMILES from exact source condition key `reaction_time`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__reagents` | Display content from exact source condition key `reagents`. Original lists/objects remain in conditions_json. |
| `cond_smiles__reagents` | Original non-empty SMILES from exact source condition key `reagents`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__solvents` | Display content from exact source condition key `solvents`. Original lists/objects remain in conditions_json. |
| `cond_smiles__solvents` | Original non-empty SMILES from exact source condition key `solvents`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__speed` | Display content from exact source condition key `speed`. Original lists/objects remain in conditions_json. |
| `cond_smiles__speed` | Original non-empty SMILES from exact source condition key `speed`. Item-level alignment is preserved in JSON and condition_items. |
| `cond__vacuum_condition` | Display content from exact source condition key `vacuum_condition`. Original lists/objects remain in conditions_json. |
| `cond_smiles__vacuum_condition` | Original non-empty SMILES from exact source condition key `vacuum_condition`. Item-level alignment is preserved in JSON and condition_items. |
| `content_signature` | Full SHA-256 of normalized reference ID plus original reactants/products/merged conditions/score object; provenance excluded. |
| `identical_content_group_size` | Number of released rows with the same content_signature. No deduplication performed. |

## Additional table-level fields

Fields shared with screening rows have the same meaning. The complete original record is retained only in table JSONL.

| Field | Meaning / provenance |
|---|---|
| `n_reconstructed_options` | Number of options retained in this uploaded table record. |
| `total_entries_in_source_table` | Unmodified meta.total_entries_in_source_table; not verified against original PDFs. |
| `n_source_entries_not_retained` | Source-entry count minus number of retained options. |
| `missing_source_row_indices` | Source row-index positions in the metadata range with no retained option; no data values imputed. |
| `constant_condition_keys` | Sorted source constant-condition keys. |
| `all_condition_keys` | Sorted union of constant and varying condition keys for this table record. |
| `yield_min_retained` | Derived minimum of source yields among retained options. |
| `yield_max_retained` | Derived maximum of source yields among retained options. |
| `yield_mean_retained` | Derived arithmetic mean of source yields among retained options. |
| `best_yield_source` | Unmodified meta.best_yield; distinct from an independently computed maximum. |
| `best_effective_score_source` | Unmodified meta.best_effective_score. |
| `best_option_indices_source` | Unmodified answer array of benchmark best-option indices. |
| `best_entry_indices_source` | Unmodified meta.best_entry_indices. |
| `global_best_yield_source` | Unmodified meta.global_best_yield. |
| `global_best_effective_score_source` | Unmodified meta.global_best_effective_score. |
| `global_best_entry_indices_source` | Unmodified meta.global_best_entry_indices. |
| `n_rows_with_invalid_condition_smiles` | Count of this table record’s retained rows with at least one condition parsing failure. |
| `source_record_json` | Complete original benchmark record; included in table JSONL, omitted from table summary CSV. |

## Reference fields

| Field | Meaning / provenance |
|---|---|
| `reference_id` | Case-folded meta.doi source identifier; used only for reference grouping. |
| `doi_inferred` | DOI string decoded using a labelled publisher-specific rule; not independently verified. |
| `doi_url_inferred` | Resolver URL formed from doi_inferred; availability has not been checked. |
| `doi_decode_rule` | Named decoding rule used for the source DOI identifier. |
| `doi_resolver_verified` | False throughout: external DOI resolution was not performed. |
| `source_id_variants` | Original case-sensitive source identifiers retained in this reference group. |
| `n_source_id_variants` | Number of original source identifier variants. |
| `n_table_records` | Number of released table records in the indicated reference group or structure usage set. |
| `n_reconstructed_rows` | Number of retained screening rows in this reference group. |
| `total_entries_in_source_table_sum` | Sum of source-entry metadata counts in this reference group. |
| `n_source_entries_not_retained` | Source-entry count minus number of retained options. |
| `table_record_ids` | Compact exported table IDs linked to this reference group or structure usage. |
| `parent_table_ids` | Unmodified parent-table IDs linked to this reference group. |
| `bibliographic_metadata_status` | Explicit note that article titles, authors, journal and dates are absent from the upload. |

## Condition-item fields

| Field | Meaning / provenance |
|---|---|
| `screening_row_id` | Derived stable identifier: input table ordinal plus original zero-based benchmark option index, e.g. T0001_O0004. |
| `table_record_id` | Derived table-record identifier, T0001–T0552, in uploaded record order. |
| `condition_key` | Original condition key, not an inferred chemical role. |
| `condition_origin` | Whether this specific condition item originated in constant_conditions or in the varying option. |
| `item_index` | Zero-based item index within the source condition-key list. |
| `content_source` | Exact source condition item content/name string. |
| `smiles_source` | Exact source SMILES string; not repaired. |
| `smiles_canonical` | RDKit-generated canonical SMILES when parsing succeeded. |
| `smiles_status` | valid, invalid, missing or not_checked; valid refers only to parsing. |
| `has_dummy_atom` | Whether a parsed structure contains an atomic-number-zero dummy/wildcard atom. |

## Reaction-molecule fields

| Field | Meaning / provenance |
|---|---|
| `table_record_id` | Derived table-record identifier, T0001–T0552, in uploaded record order. |
| `parent_table_id` | Unmodified meta.parent_table_id; preserve case. Not a printed table label. |
| `reference_id` | Case-folded meta.doi source identifier; used only for reference grouping. |
| `side` | Original reaction-side array: reactants or products. |
| `molecule_index` | Zero-based molecule index within the original reaction-side array. |
| `name` | Exact source reaction-side content/name string. |
| `smiles_source` | Exact source SMILES string; not repaired. |
| `smiles_canonical` | RDKit-generated canonical SMILES when parsing succeeded. |
| `smiles_status` | valid, invalid, missing or not_checked; valid refers only to parsing. |
| `has_dummy_atom` | Whether a parsed structure contains an atomic-number-zero dummy/wildcard atom. |

## Distinct-SMILES audit fields

| Field | Meaning / provenance |
|---|---|
| `smiles_source` | Exact source SMILES string; not repaired. |
| `smiles_canonical` | RDKit-generated canonical SMILES when parsing succeeded. |
| `smiles_status` | valid, invalid, missing or not_checked; valid refers only to parsing. |
| `has_dummy_atom` | Whether a parsed structure contains an atomic-number-zero dummy/wildcard atom. |
| `source_names` | Distinct original names associated with this SMILES, without reconciliation. |
| `source_roles` | Original reaction sides or condition keys associated with this SMILES. |
| `n_table_records` | Number of released table records in the indicated reference group or structure usage set. |
| `n_screening_rows_as_condition` | Number of retained screening rows using this SMILES as a condition item. |
| `n_occurrences_in_exported_molecule_and_condition_item_tables` | Occurrence count across reaction_molecules (once per table occurrence) and condition_items (once per retained row occurrence). |
| `table_record_ids` | Compact exported table IDs linked to this reference group or structure usage. |
| `screening_row_ids_as_condition` | All retained row IDs using this structure as a condition item. |

## Quality-issue fields

| Field | Meaning / provenance |
|---|---|
| `severity` | Information/warning classification, not a calibrated confidence score. |
| `scope` | Issue applies to a row, table or repeated-content group. |
| `record_id` | ID corresponding to scope; repeated-content groups use the full content signature. |
| `issue_code` | Machine-readable issue category. |
| `detail` | Issue details as JSON-encoded text or source-preserving explanatory text. |
