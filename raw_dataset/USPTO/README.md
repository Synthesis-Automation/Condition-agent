# USPTO condition dataset cleanup

`USPTO_condition_pred_category_maped_rxn.csv` is retained as the immutable
source artifact. `USPTO_condition_reactions_cleaned.csv` is the database-ready
subset produced by:

```powershell
python -m raw_dataset.USPTO.clean_uspto_dataset
```

The cleanup keeps a row only when:

- at least one catalyst or reagent is present;
- every supplied condition is valid SMILES;
- the canonical and mapped reactions contain parseable structures;
- the product-side heavy atoms have complete, internally consistent mapping;
- every product map has a reactant-side source with the same element;
- the mapping describes at least one bond, hydrogen, or charge change; and
- RXNMapper confidence is at least `0.5`.

Exact duplicates are removed after projection to the output schema. The
source-split and reaction-label columns `dataset`, `rxn_category`, and
`rxn_class_name` are intentionally omitted. The patent identifier, raw
condition SMILES, mapped and unmapped reactions, and mapping confidence are
preserved. Counts and rejection reasons are recorded in
`USPTO_condition_reactions_cleanup.json`.
