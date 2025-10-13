# Quick Start: Interactive CLI

## 1. Set up your API key

```bash
# For Aliyun (default)
export DASHSCOPE_API_KEY=your_api_key_here

# Or for OpenAI
export OPENAI_API_KEY=your_api_key_here
```

## 2. Run the CLI

```bash
python app/cli_recommend.py
```

## 3. Minimal Example

Input a valid reaction SMILES when prompted:

```
Reaction SMILES: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Requirements: (press Enter to skip)
```

That's it! The system only requires valid SMILES - all constraints are optional.

## 4. Example with Constraints

```
Reaction SMILES: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Requirements: no strong base, room temperature, avoid expensive catalysts
```

The LLM will parse this and extract:
- base_strength: moderate
- temperature.max: 30
- cost_level: low

## 5. What Makes a Valid SMILES?

✅ Valid format:
- `reactants>>products`
- Must have `>>` separator
- Non-empty reactants side
- Multiple reactants separated by `.`

❌ Invalid format:
- `invalid`
- `C=C` (no products)
- `>>products` (no reactants)
- `broken>>` (empty products)

## 6. If SMILES is Invalid

The CLI will ask you to fix it:

```
❌ VALIDATION ISSUES (must fix):
1. Invalid reaction SMILES format

Please fix the validation issues above:
> Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
```

## 7. Optional Clarifications

If your requirements are ambiguous, the LLM may ask for clarification:

```
⚠️  CLARIFICATIONS NEEDED (optional):
1. You mentioned "room temperature" - should this be strict (max 25°C) or flexible?

Please provide clarifications (or press Enter to skip):
> Flexible up to 40C is fine
```

**Note**: These are optional - you can press Enter to skip if you don't want to provide more detail.

## 8. Submission

Once valid, you'll see a preview and confirmation:

```
✅ Request is ready for submission!

Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Type: Suzuki_Miyaura
Constraints: 2 specified

Submit this request? (yes/no): yes
```

## Common Use Cases

### Case 1: Just Want Recommendations (No Constraints)

```
Reaction: CCBr.O=C(O)C>>CCC(=O)O
Requirements: (Enter)
```

### Case 2: Temperature Limit

```
Reaction: CCBr.O=C(O)C>>CCC(=O)O
Requirements: no high temperature
```

### Case 3: Multiple Constraints

```
Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Requirements: room temp, no strong base, cheap catalysts, prefer DMF
```

### Case 4: Prefer Specific Metal Catalyst

```
Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Requirements: use copper catalyst
```

### Case 5: Exclude Specific Reagent

```
Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Requirements: no Pd(PPh3)4, prefer other Pd source
```

### Case 6: Mixed Positive and Negative Constraints

```
Reaction: Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Requirements: use copper catalyst, no strong base, room temperature, avoid DMF
```

## Tips

1. **SMILES First**: Make sure your reaction SMILES is correct - this is the only required field
2. **Be Natural**: Describe requirements in plain English - the LLM will understand
3. **Skip if Unsure**: All constraints are optional - press Enter to skip if you don't know
4. **Iterate**: If the LLM misunderstands, you can clarify in the refinement loop
5. **Save Draft**: If interrupted, your progress is saved to `draft_request.json`

## Troubleshooting

**"Failed to initialize LLM client"**
→ Set your API key: `export DASHSCOPE_API_KEY=your_key`

**"Invalid reaction SMILES format"**
→ Make sure format is `reactants>>products` with `>>` separator

**LLM asks too many questions**
→ Press Enter to skip optional clarifications

**Want to exit?**
→ Type `quit` or press Ctrl+C

## Help

```bash
python app/cli_recommend.py --help
```
