# Equivalents Calculation from Rule Files

## Overview

The rule-based recommendation system now automatically calculates chemical equivalents from the `amount` field in the standardized rule format. This feature parses percentage-based amounts and converts them to equivalents suitable for experimental protocols.

## Implementation

### Parsing Logic

The `_parse_amount_to_equivalents()` function in `chemtools/output_formatter.py` handles various amount formats:

#### 1. Explicit Equivalents
```json
{
  "amount": "2.0 eq"
}
```
→ `equivalents: 2.0`

#### 2. Percentage Ranges (Catalysts/Ligands)
```json
{
  "name": "Pd2(dba)3",
  "role": "metal_source",
  "amount": "0.5–1.5%"
}
```
→ `equivalents: 0.01` (midpoint: (0.5 + 1.5) / 2 / 100 = 0.01)

For roles like `metal_source`, `metal_precursor`, `catalyst`, and `ligand`, percentages are treated as **mol%** and divided by 100.

#### 3. Percentage (Bases/Reagents)
```json
{
  "name": "K2CO3",
  "role": "base",
  "amount": "200%"
}
```
→ `equivalents: 2.0` (200 / 100 = 2.0)

For roles like `base`, `additive`, and `partner`, percentages represent direct equivalents.

#### 4. No Amount Specified
```json
{
  "name": "THF/H2O",
  "role": "solvent",
  "amount": ""
}
```
→ `equivalents: null`

### Role-Specific Conversion

The conversion logic differs by chemical role:

**Mol% Roles** (divide by 100):
- `metal_source`
- `metal_precursor`
- `catalyst`
- `ligand`

**Direct Equivalent Roles** (divide by 100):
- `base`
- `additive`
- `partner`
- `boron_partner`

**No Conversion**:
- `solvent` (equivalents remain null)
- `reagent` (default fallback)

### Range Handling

When an amount is specified as a range (e.g., "0.5–1.5%", "1.0-3.0%"), the parser:
1. Extracts both bounds
2. Calculates the midpoint
3. Applies role-specific conversion

Supported separators: `–` (en dash), `-` (hyphen), `~` (tilde)

## Example Output

### Rule File Entry
```json
{
  "id": "SCDB-SUZ-ARBRI-GENERAL-SPhos",
  "reagents": [
    {
      "name": "Pd2(dba)3",
      "role": "metal_source",
      "amount": "0.5–1.5%"
    },
    {
      "name": "SPhos",
      "role": "ligand",
      "amount": "1.0–3.0%"
    },
    {
      "name": "K2CO3",
      "role": "base",
      "amount": "200%"
    },
    {
      "name": "THF/H2O (4:1)",
      "role": "solvent",
      "amount": ""
    }
  ]
}
```

### Recommendation Output
```json
{
  "chemicals": [
    {
      "name": "Tris(dibenzylideneacetone)dipalladium(0)",
      "abbreviation": "Pd2(dba)3",
      "role": "metal_precursor",
      "equivalents": 0.01
    },
    {
      "name": "SPhos",
      "role": "ligand",
      "equivalents": 0.02
    },
    {
      "name": "Potassium carbonate",
      "abbreviation": "K2CO3",
      "role": "base",
      "equivalents": 2.0
    },
    {
      "name": "THF/H2O",
      "role": "solvent",
      "equivalents": null
    }
  ]
}
```

## Code Flow

1. **Rule Matching**: `chem.rules.match()` finds matching rule entry
2. **Extract Reagents**: `_convert_rule_match_to_recommendations()` iterates through `reagents` array
3. **Parse Amount**: For each reagent, `_normalize_rule_string_value()` calls `_parse_amount_to_equivalents(amount, role)`
4. **Calculate**: Based on role and amount format, equivalents are calculated
5. **Enrich**: Database lookup adds full name, CAS, SMILES
6. **Output**: Final JSON includes equivalents for each chemical

## Testing

Run the test script to verify equivalents calculation:

```bash
python test_equivalents.py
```

Expected output:
```
=== Recommendation 1 ===
Chemicals:
  - Pd2(dba)3                      [metal_precursor] equivalents: 0.010
  - SPhos                          [ligand         ] equivalents: 0.020
  - K2CO3                          [base           ] equivalents: 2.000
  - THF/H2O                        [solvent        ] equivalents: None
```

## Backward Compatibility

The implementation maintains backward compatibility:

### Legacy Format Support
If no `reagents` array exists, the system falls back to extracting from legacy keys:
```json
{
  "pd_source": "Pd2(dba)3 (1 mol%)",
  "ligand": "SPhos (2 mol%)",
  "base": "K2CO3 (2 eq)"
}
```

The parser attempts to extract "1 eq" or "2 mol%" patterns from the value strings.

### Fallback Logic
```python
# First: Try parsing from amount parameter (new format)
if amount:
    equivalents = _parse_amount_to_equivalents(amount, role)

# Fallback: Extract from value string (legacy format)
if equivalents is None:
    eq_match = re.search(r"([+-]?\d+(?:\.\d+)?)\s*eq", value)
```

## Benefits

1. **Experimental Accuracy**: Chemists get precise loading values for experimental protocols
2. **Consistency**: All recommendations include standardized equivalents
3. **Automation**: No manual calculation needed from percentage values
4. **Flexibility**: Supports ranges, explicit eq, and percentage formats
5. **Role-Aware**: Correctly interprets mol% for catalysts vs. equivalents for bases

## Future Enhancements

Potential improvements:
- Support for weight % (requires molecular weight data)
- Volume-based amounts for solvents (mL, L)
- Concentration-based specifications (M, mM)
- Custom conversion rules per reaction family
