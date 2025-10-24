# Batch Reagent Registry Addition - Usage Example

## Test Run Results

Successfully tested the batch reagent processor on the mapped markdown file.

### Test 1: High-Frequency Additives

```powershell
python batch_add_reagents_to_registry.py --dry-run --max-reagents 3 --min-occurrences 1000
```

**Results:**
- ✅ water (2,696 occurrences) - CAS found: 7732-18-5 → Would add
- ✅ Potassium 2-ethylhexanoate (1,722 occurrences) - CAS found: 3164-85-0 → Would add
- ⚠️ "3" (1,606 occurrences) - CAS not found (generic placeholder name)

**Summary:** 2/3 reagents ready to add (66% success rate)

### Test 2: Condensation Agents (Coupling Reagents)

```powershell
python batch_add_reagents_to_registry.py --dry-run --max-reagents 5 --min-occurrences 500 --skip-roles additive base ligand metal_precursor preformed_metal_catalyst solvent
```

**Results:**
- ✓ DIC (1,718 occurrences) - Already exists in registry
- ✓ T3P (862 occurrences) - Already exists in registry
- ✓ BOP (604 occurrences) - Already exists in registry
- ✅ PyBrOP (524 occurrences) - CAS found: 132705-51-2 → Would add

**Summary:** 3/4 already exist, 1 new reagent ready to add

**✅ Validation:** Coupling reagents correctly mapped to `condensation_agent` role

## Recommended Production Workflow

### Step 1: Preview High-Priority Reagents

Focus on condensation_agents (coupling reagents), bases, and catalysts with high occurrence counts:

```powershell
# Condensation agents (coupling reagents)
python batch_add_reagents_to_registry.py --dry-run --min-occurrences 100 --skip-roles additive ligand solvent --output preview_condensation.json

# Bases
python batch_add_reagents_to_registry.py --dry-run --min-occurrences 100 --skip-roles additive ligand solvent condensation_agent metal_precursor preformed_metal_catalyst --output preview_bases.json
```

### Step 2: Add High-Priority Reagents

```powershell
# Add coupling reagents with ≥100 occurrences
python batch_add_reagents_to_registry.py --min-occurrences 100 --delay 1.5 --skip-roles additive base ligand metal_precursor preformed_metal_catalyst solvent --output add_condensation_agents.json

# Add bases with ≥100 occurrences
python batch_add_reagents_to_registry.py --min-occurrences 100 --delay 1.5 --skip-roles additive ligand solvent condensation_agent metal_precursor preformed_metal_catalyst --output add_bases.json

# Add metal precursors with ≥100 occurrences
python batch_add_reagents_to_registry.py --min-occurrences 100 --delay 1.5 --skip-roles additive base ligand solvent condensation_agent preformed_metal_catalyst --output add_metal_precursors.json
```

### Step 3: Add Medium-Priority Reagents

```powershell
# Lower threshold to 50 occurrences
python batch_add_reagents_to_registry.py --min-occurrences 50 --delay 2.0 --skip-roles additive ligand solvent --output add_medium_priority.json
```

### Step 4: Review and Add Remaining

```powershell
# Add all remaining reagents (may have lower quality CAS matches)
python batch_add_reagents_to_registry.py --min-occurrences 10 --delay 2.0 --output add_all_remaining.json
```

### Step 5: Review Results

```powershell
# View summary statistics
cat add_condensation_agents.json | ConvertFrom-Json | Select-Object total_processed, already_exists, cas_not_found, added_successfully, failed

# View which reagents were added
cat add_condensation_agents.json | ConvertFrom-Json | Select-Object -ExpandProperty details | Where-Object { $_.status -eq "added" } | Select-Object name, role, cas, family
```

## Expected Results by Role

Based on test runs, here's what to expect:

### Condensation Agent (36 reagents)

**Likely already exist:** DIC, T3P, BOP, EDC, CDI (common coupling reagents)

**Expected new additions:** 5-10 specialized coupling reagents (PyBrOP, COMU, TCFH, etc.)

### Base (63 reagents)

**Likely already exist:** K₃PO₄, Cs₂CO₃, K₂CO₃, Et₃N (common bases)

**Expected new additions:** 20-30 specialized bases

### Metal Precursor (114 reagents)

**Likely already exist:** CuI, PdCl₂, Pd(OAc)₂ (common catalysts)

**Expected new additions:** 40-60 specialized metal salts

### Preformed Metal Catalyst (27 reagents)

**Likely already exist:** Common Buchwald precatalysts

**Expected new additions:** 10-15 specialized precatalysts

### Ligand (146 reagents)

**Likely already exist:** XPhos, BrettPhos, SPhos, PPh₃ (common ligands)

**Expected new additions:** 50-80 specialized ligands

### Additive (97 reagents)

**Expected new additions:** 30-50 (many are generic or have unclear names)

### Solvent (69 reagents)

**Likely already exist:** Most common solvents (DMF, THF, PhMe, etc.)

**Expected new additions:** 10-20 specialized solvents

## Troubleshooting Common Issues

### Issue: Generic Names Not Found

**Example:** Reagent named "3" or "N" or "1"

**Solution:** These are placeholder names in the dataset. Options:
1. Skip them (they're likely already in registry with proper names)
2. Manually identify and add via UI
3. Add to skip list in script

### Issue: Abbreviations Not Resolved

**Example:** "tBuXPhos" vs "2-Di-tert-butylphosphino-2',4',6'-triisopropylbiphenyl"

**Solution:** PubChem search works better with:
- Full chemical names
- Common abbreviations (like "DMF")
- Trade names

For specialized ligands/catalysts, may need manual addition.

### Issue: Already Exists

**Status:** This is expected and good! It means your registry already has good coverage.

**Action:** Review the "already_exists" count to see coverage gaps.

## Quality Control Checklist

After batch processing:

- [ ] Review `added_successfully` count - reasonable?
- [ ] Check `cas_not_found` list - are these real reagents or placeholders?
- [ ] Verify `failed` entries - can they be added manually?
- [ ] Spot-check 5-10 added entries in registry JSON files
- [ ] Confirm families are appropriate for the role
- [ ] Check abbreviations and aliases look correct
- [ ] Validate SMILES strings if present
- [ ] Commit changes to git

## Performance Expectations

With default settings (`--delay 1.0`):

- **10 reagents:** ~15-20 seconds
- **50 reagents:** ~1-2 minutes
- **100 reagents:** ~2-3 minutes
- **500 reagents:** ~10-15 minutes (full dataset)

Rate limiting delay prevents PubChem API throttling.

## Success Metrics

Target metrics for a successful batch run:

- ✅ **CAS found rate:** >70% (some generic names expected to fail)
- ✅ **Already exists rate:** 30-50% (shows good existing coverage)
- ✅ **Error rate:** <5% (network issues, API errors)
- ✅ **Added successfully:** 100-200 reagents (out of 552 total)

## Next Steps After Batch Addition

1. **Validate entries** using reagent_taxonomy_ui.py
2. **Assign correct families** (may need refinement from defaults)
3. **Add role-specific properties:**
   - Ligands: donors, denticity
   - Bases: basicity, nucleophilicity
   - Catalysts: metal, oxidation states
4. **Generate embeddings** for similarity search
5. **Test recommendations** with new reagents included

---

**Ready to proceed with production batch processing!** 🚀
