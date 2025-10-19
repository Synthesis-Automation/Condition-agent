# Enhanced Visualization Feature

## What's New

The SMARTS Generator now creates **comprehensive visualization images** that show:

### 1. Full Reaction View
The complete transformation in one glance

### 2. Individual Pattern Components
Side-by-side view of:
- **Reactant Pattern** (labeled in blue) - Shows what substrates match
- **Product Pattern** (labeled in green) - Shows what products form

### 3. SMARTS Text
Each pattern component includes the actual SMARTS string for reference

## Example Output

When you run:
```powershell
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCI.B2pin2>>CCCBpin" \
  --visualize
```

The `core_transformation.png` image contains:

```
┌─────────────────────────────────────────────────────────────┐
│  Reaction SMARTS Pattern Visualization                      │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Complete Reaction Transformation                           │
│  [Reactant Structure] ──→ [Product Structure]               │
│                                                              │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Individual Pattern Components:                             │
│                                                              │
│  Reactant Pattern              Product Pattern              │
│  [Structure diagram]           [Structure diagram]          │
│  [C:1]-[I:2]                   [C:1]-[B:3]                  │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

## Benefits

### ✅ See Both Views at Once
- Top: Overall transformation
- Bottom: Detailed components

### ✅ Understand Pattern Scope
- Blue reactant pattern shows what substrates match
- Green product pattern shows expected products

### ✅ SMARTS Reference
- Pattern text included for copy/paste
- Easy to verify against your protocol

### ✅ Better Documentation
- Single image shows complete picture
- Include in presentations and reports

## Comparison: Before vs After

### Before (Old Version)
- Single reaction drawing only
- No individual pattern views
- No SMARTS text display

### After (New Enhanced Version)
- ✅ Full reaction transformation
- ✅ Individual reactant pattern (blue label)
- ✅ Individual product pattern (green label)
- ✅ SMARTS text for each component
- ✅ Clear section labels and organization

## Usage Example

### Generate Enhanced Visualization

```powershell
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCCCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCCCCCB1OC(C)(C)C(C)(C)O1" \
  --visualize \
  --viz-dir my_protocol_viz
```

### What You Get

**File:** `my_protocol_viz/core_transformation.png`

Contains:
1. **Header** - "Reaction SMARTS Pattern Visualization"
2. **Section 1** - "Complete Reaction Transformation"
   - Full reaction showing reactant → product
3. **Section 2** - "Individual Pattern Components"
   - Left side: Reactant pattern with blue label and SMARTS text
   - Right side: Product pattern with green label and SMARTS text

### Guard Patterns

**Files:** `guard_forbid_1.png`, `guard_forbid_2.png`, etc.

Each shows:
- Visual structure of the forbidden pattern
- Helpful for understanding what substrates are excluded

## Tips

### 1. Use for Protocol Documentation
Include the `core_transformation.png` in your protocol documentation to show:
- What the protocol does (full reaction)
- What substrates work (reactant pattern)
- What products form (product pattern)

### 2. Verify Pattern Correctness
Look at the individual patterns to ensure:
- Reactant pattern matches intended substrates
- Product pattern shows expected transformation
- SMARTS text is correct

### 3. Share with Team
Single comprehensive image is perfect for:
- Team meetings
- Protocol reviews
- Training materials
- Literature documentation

## Technical Details

### Image Layout

```
Total image size: ~1200 x 500+ pixels (auto-sized)

Layout:
├─ Title (top)
├─ Full reaction section (~1000 x 300 px)
│  └─ Centered reaction drawing
├─ Component label
└─ Individual patterns section
   ├─ Reactant pattern (left, ~500 x 300 px)
   │  ├─ Blue label
   │  ├─ Structure
   │  └─ SMARTS text
   └─ Product pattern (right, ~500 x 300 px)
      ├─ Green label
      ├─ Structure
      └─ SMARTS text
```

### Dependencies

- **RDKit**: For molecular structure drawing
- **Pillow (PIL)**: For image composition and labels
- **Fallback**: If Pillow not available, creates simple reaction drawing

### Installation

If you see "PIL/Pillow not available" warning:
```powershell
pip install Pillow
```

## Examples

### Example 1: Alkyl Iodide Borylation

**Command:**
```powershell
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCI.B2pin2>>CCCBpin" \
  --visualize
```

**Result:** `smarts_visualizations/core_transformation.png` shows:
- **Top:** C-I + B2pin2 → C-Bpin transformation
- **Bottom Left:** Alkyl iodide pattern structure
- **Bottom Right:** Alkyl boronic ester pattern structure

### Example 2: With Testing

**Command:**
```powershell
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCI.B2pin2>>CCCBpin" \
  --visualize \
  --test-smiles "CCCCI,CC(C)I,CC(C)(C)I"
```

**Result:**
1. Enhanced visualization images
2. Test results showing which molecules match
3. Easy to correlate structure to match results

## Summary

The enhanced visualization provides:
- ✅ **Complete context** - Full reaction + individual patterns
- ✅ **Better understanding** - See both overall and detailed views
- ✅ **Easy verification** - SMARTS text included for checking
- ✅ **Professional output** - Suitable for documentation and presentations
- ✅ **Single comprehensive image** - Everything in one file

**Try it now:**
```powershell
python -m chemtools.protocol.smarts_generator_cli --interactive
```

When prompted for visualization, answer "Y" to see the enhanced images!
