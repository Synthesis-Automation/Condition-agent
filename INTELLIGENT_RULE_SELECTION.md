# Intelligent Rule-Based Condition Selection

## Overview

The split mode now features **intelligent substrate-aware condition selection** that automatically identifies the most suitable reaction conditions from rule files based on detected substrate features.

## How It Works

### 1. Feature Detection
When you query a reaction, the system:
- Analyzes the reaction SMILES using `FeatureAnalyzer`
- Detects chemical features from reactants (halide types, functional groups, etc.)
- Examples:
  - `sp2_chloride_present` - Aryl chloride substrate
  - `sp2_bromide_present` - Aryl bromide substrate
  - `primary_amine_present` - Primary amine nucleophile
  - `aniline_present` - Aniline derivative
  - `ortho_substitution_present` - Ortho-substituted substrate

### 2. Rule Matching with Specificity Scoring
The system evaluates all base_rules in the rule file and scores them by:

**Priority 1: Feature Specificity** (Higher = Better)
- Specific halides (`sp2_chloride_present`, `sp2_bromide_present`): Score 10
- Challenging features (`ortho_substitution_present`, `electron_poor_boronate_present`): Score 8
- General features (`sp2_halide_present`, `sp2_boron_present`): Score 3

**Priority 2: Number of Matched Features**
- More matched features = better understanding of substrate

**Priority 3: Rule Focus**
- Rules with fewer total conditions are more focused

### 3. Best Match Selection
- ✓ Marks the best matched rule based on substrate features
- Shows which features were matched
- Displays full conditions for that specific rule
- Lists applicable alternatives with ✓ indicators

## Examples

### Example 1: Aryl Chloride Suzuki Coupling

**Query:**
```bash
python app/unified_interactive_cli.py --rxn "Clc1ccc(C(F)(F)F)cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccc(C(F)(F)F)cc2)cc1" --split
```

**Result:**
```
✓ Best Matched Conditions: CataCXium A–Pd G3 / Cs2CO3
  Matched features: sp2_chloride_present
  (Frequently near the top; complementary to XPhos/P(tBu)3 choices.)
  • Precatalyst: CataCXium A–Pd G3
  • Ligand: CataCXium A
  • Temperature: 90–110 °C
```

**Why this was selected:**
- Detected `sp2_chloride_present` (specificity score: 10)
- CataCXium A rule specifically targets aryl chlorides
- Better than general dtbpf rule for challenging substrates

### Example 2: Primary Aniline C-N Coupling

**Query:**
```bash
python app/unified_interactive_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --split
```

**Result:**
```
✓ Best Matched Conditions: Primary amines (aliphatic or aniline)
  Matched features: primary_amine_present, aniline_present
  (Primary NH2 nucleophiles (aliphatic or aryl). First shot: t‑AmOH/Cs2CO3.)
  • Ligand: tBuBrettPhos or BippyPhos
  • Solvent: MeCN, t‑AmOH, dioxane, or toluene
  • Base: K3PO4 (aq) or Cs2CO3
```

**Why this was selected:**
- Detected both `primary_amine_present` and `aniline_present`
- 2 matched features with high combined specificity
- More targeted than generic amine coupling conditions

### Example 3: Regular Aryl Bromide

**Query:**
```bash
python app/unified_interactive_cli.py --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" --split
```

**Result:**
```
✓ Best Matched Conditions: CataCXium A–Pd G3 / Cs2CO3
  Matched features: sp2_bromide_present
  • Precatalyst: CataCXium A–Pd G3
  • Temperature: 90–110 °C

Alternative Conditions:
  1. dtbpf / K3PO4 (most frequent hit) ✓
     General, scalable; good for broad aryl bromide/iodide...
```

## Feature Specificity Hierarchy

The system uses a built-in specificity hierarchy:

| Feature | Specificity Score | Description |
|---------|-------------------|-------------|
| `sp2_chloride_present` | 10 | Specific halide type (challenging) |
| `sp2_bromide_present` | 10 | Specific halide type (moderate) |
| `sp2_iodide_present` | 10 | Specific halide type (easy) |
| `ortho_substitution_present` | 8 | Sterically hindered |
| `electron_poor_boronate_present` | 8 | Less reactive nucleophile |
| `five_member_heteroaryl_halide_present` | 8 | Special substrate class |
| `base_sensitive_functionality_present` | 7 | Requires special handling |
| `acid_labile_group_present` | 7 | Requires special handling |
| `sp2_halide_present` | 3 | General feature |
| `sp2_boron_present` | 3 | General feature |

## Benefits

1. **Automatic Optimization**: No manual rule selection needed
2. **Substrate-Aware**: Matches conditions to your specific substrate
3. **Challenging Substrates**: Prioritizes specialized conditions for difficult cases
4. **Transparent Reasoning**: Shows which features were matched
5. **Alternative Options**: Still provides other applicable conditions
6. **Research-Ready**: Conditions are immediately actionable in the lab

## Technical Details

### Scoring Algorithm

For each base_rule in the rule file:

1. **Check feature matching:**
   - `all`: All listed features must be present
   - `any`: At least one listed feature must be present

2. **Calculate score:**
   ```python
   score = (
       specificity_score,    # Sum (for 'all') or max (for 'any') of matched features
       matched_count,        # Number of matched features
       -total_conditions     # Negative (prefer focused rules)
   )
   ```

3. **Select best match:**
   - Sort rules by score (descending)
   - Pick highest-scored rule
   - Show matched features to user

### Fallback Behavior

If feature detection fails or no rules match:
- Falls back to showing default_rule conditions
- Labeled as "Default Conditions" instead of "Best Matched Conditions"
- Still provides all alternative conditions

## Usage Tips

1. **Trust the Match**: The ✓ Best Matched Conditions are optimized for your substrate
2. **Check Alternatives**: Other ✓-marked conditions are also applicable
3. **Challenging Substrates**: System automatically selects stronger conditions when needed
4. **Multiple Features**: Rules matching multiple features are preferred
5. **Specificity Wins**: Specific features (like `sp2_chloride`) outrank general ones (like `sp2_halide`)

## Future Enhancements

Potential improvements:
- Machine learning-based feature importance
- Historical success rate tracking
- User feedback integration
- Custom specificity scoring per reaction family
- Quantitative yield predictions
