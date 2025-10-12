# Sample CAS Numbers for Pure LLM Workflow Testing

## Test Suite Overview

This document provides a comprehensive list of diverse CAS numbers to test and improve the Pure LLM workflow. The samples are organized by reagent category to ensure broad coverage.

## Test Samples (16 total)

### 1. Bases (3 samples)

| CAS | Name | Expected Role | Expected Family | Notes |
|-----|------|---------------|-----------------|-------|
| 121-44-8 | Triethylamine | base | tertiary_amines_aliphatic | Common organic base, moderate basicity |
| 1633-05-2 | Potassium tert-butoxide | base | alkoxides | Strong base, very nucleophilic |
| 7664-38-2 | Phosphoric acid | acid | inorganic_acids | Weak acid, commonly used |

### 2. Catalysts (3 samples)

| CAS | Name | Expected Role | Expected Family | Notes |
|-----|------|---------------|-----------------|-------|
| 14221-01-3 | Tetrakis(triphenylphosphine)palladium(0) | catalyst | palladium_catalysts | Pd(0) for cross-coupling reactions |
| 31274-51-8 | Bis(diphenylphosphino)ferrocene]dichloropalladium(II) | catalyst | palladium_catalysts | Pd(II) with dppf ligand |
| 14808-79-8 | Grubbs Catalyst 1st Gen | catalyst | ruthenium_catalysts | Olefin metathesis catalyst |

### 3. Solvents (3 samples)

| CAS | Name | Expected Role | Expected Family | Notes |
|-----|------|---------------|-----------------|-------|
| 67-56-1 | Methanol | solvent | alcohols | Protic polar solvent |
| 109-99-9 | Tetrahydrofuran | solvent | ethers | Aprotic polar solvent |
| 127-19-5 | N,N-Dimethylacetamide | solvent | amides | Polar aprotic solvent |

### 4. Ligands (2 samples)

| CAS | Name | Expected Role | Expected Family | Notes |
|-----|------|---------------|-----------------|-------|
| 603-35-0 | Triphenylphosphine | ligand | phosphine_ligands | Monodentate phosphine |
| 154-23-4 | BINAP | ligand | phosphine_ligands | Chiral bidentate phosphine |

### 5. Reducing Agents (1 sample)

| CAS | Name | Expected Role | Expected Family | Notes |
|-----|------|---------------|-----------------|-------|
| 16940-66-2 | Sodium borohydride | reducing_agent | metal_hydrides | Mild reducing agent |

### 6. Lewis Acids (1 sample)

| CAS | Name | Expected Role | Expected Family | Notes |
|-----|------|---------------|-----------------|-------|
| 7446-70-0 | Aluminum chloride | lewis_acid | metal_halides | Strong Lewis acid catalyst |

### 7. Oxidizing Agents (1 sample)

| CAS | Name | Expected Role | Expected Family | Notes |
|-----|------|---------------|-----------------|-------|
| 7722-84-1 | Hydrogen peroxide | oxidizing_agent | peroxides | Common mild oxidizer |

### 8. Coupling Reagents (1 sample)

| CAS | Name | Expected Role | Expected Family | Notes |
|-----|------|---------------|-----------------|-------|
| 2524-03-0 | Dicyclohexylcarbodiimide (DCC) | coupling_reagent | carbodiimides | Peptide coupling reagent |

### 9. Additives (1 sample)

| CAS | Name | Expected Role | Expected Family | Notes |
|-----|------|---------------|-----------------|-------|
| 7757-93-9 | Calcium phosphate | additive | inorganic_salts | Inorganic additive |

## Additional Challenging Cases

### Edge Cases to Test

| CAS | Name | Expected Role | Challenge |
|-----|------|---------------|-----------|
| 7647-01-0 | Hydrochloric acid | acid | Both reagent and solvent in some contexts |
| 71-43-2 | Benzene | solvent | Potentially misclassified as starting material |
| 7732-18-5 | Water | solvent | Multi-functional (solvent, reagent, product) |
| 64-17-5 | Ethanol | solvent | Can be both solvent and protecting group removal |
| 7758-98-7 | Copper(II) sulfate | oxidizing_agent / lewis_acid | Dual role |

### Complex Catalysts

| CAS | Name | Expected Role | Notes |
|-----|------|---------------|-------|
| 308135-20-0 | XPhos Pd G3 | catalyst | Modern precatalyst, complex structure |
| 1314132-12-3 | Ruphos Pd G3 | catalyst | Air-stable precatalyst |
| 244187-81-3 | [(1,5-COD)IrCl]2 | catalyst | Iridium dimer catalyst |

## Testing Strategy

### Phase 1: Core Functionality (8 samples)
Test one sample from each major category:
1. Triethylamine (base)
2. Pd(PPh3)4 (catalyst)
3. Methanol (solvent)
4. Triphenylphosphine (ligand)
5. Sodium borohydride (reducing agent)
6. Aluminum chloride (Lewis acid)
7. Hydrogen peroxide (oxidizer)
8. DCC (coupling reagent)

**Expected Results**:
- Role classification accuracy: >90%
- Field assignment completeness: >80%
- Verification pass rate: >70%

### Phase 2: Accuracy Refinement (16 samples)
Test all samples in the main list to identify patterns in misclassifications.

**Focus Areas**:
- Confidence scores for each classification
- Common errors in family selection
- Field completeness by role type
- Verification issues by category

### Phase 3: Edge Cases (5 samples)
Test challenging multi-functional reagents.

**Goals**:
- Handle ambiguous cases gracefully
- Provide reasoning for difficult classifications
- Suggest review when uncertainty is high

## Running the Tests

###Quick Test (Component Testing)
```bash
cd data-processor
python quick_llm_test.py
```

This tests:
- Step 2: Role classification (8 samples)
- Step 3: Field assignment (based on detected roles)
- Step 4: Verification (based on assigned fields)

**Output**: Console summary with accuracy stats

### Full Test (End-to-End Workflow)
```bash
cd data-processor
python test_pure_llm_samples.py
```

This tests:
- Complete workflow (identity → role → fields → verification)
- All 16 samples
- Detailed JSON output

**Output**: `test_pure_llm_results.json` with full results

### Manual Testing (UI)
1. Launch UI: `python reagent_taxonomy_qt.py`
2. Select "Pure LLM Workflow"
3. For each sample:
   - Enter CAS number
   - Select "Auto-detect (LLM)" for role
   - Click Generate
   - Review output

## Success Criteria

### Quantitative Metrics
- **Role Classification Accuracy**: ≥90% match with expected roles
- **Average Confidence**: ≥85% for correct classifications
- **Field Assignment Success**: ≥80% complete and valid
- **Verification Approval**: ≥70% entries pass quality check
- **Performance**: <10 seconds per reagent (all 4 steps)

### Qualitative Assessment
- **Reasoning Quality**: LLM provides chemically sound explanations
- **Error Handling**: Graceful degradation for missing/invalid data
- **Consistency**: Same CAS number yields same classification across runs
- **Completeness**: All required fields populated with valid values

## Improvement Areas

Based on testing, focus improvements on:

1. **Role Ambiguity**: Multi-functional reagents (e.g., Cu salts as both oxidizer and Lewis acid)
2. **Family Selection**: Ensure families match role requirements
3. **Field Validation**: Check allowed values are used correctly
4. **Missing Data**: Handle incomplete PubChem records gracefully
5. **Verification Rigor**: Balance between strictness and flexibility

## Expected Challenges

### 1. Multifunctional Reagents
- **Example**: HCl (acid AND solvent)
- **Solution**: Prompt engineering to pick primary role in organic synthesis context

### 2. Complex Catalysts
- **Example**: Pd(II) complexes with multiple ligands
- **Solution**: Enhanced identity resolution from PubChem

### 3. Rare Reagents
- **Example**: Novel organocatalysts
- **Solution**: Fallback to family definitions when examples are sparse

### 4. Inorganic Reagents
- **Example**: Metal salts without clear organic chemistry role
- **Solution**: Use stoichiometry and reaction context (when available)

## Next Steps After Testing

1. **Analyze Results**: Identify patterns in errors/failures
2. **Refine Prompts**: Adjust chemistry-specific instructions
3. **Enhance Examples**: Add more diverse reagent examples to family definitions
4. **Update Families**: Create missing families for common reagent types
5. **Documentation**: Document best practices and limitations
6. **User Feedback**: Collect feedback from chemists using the system

## Test Execution Checklist

- [ ] Configure LLM client (Aliyun/DeepSeek credentials)
- [ ] Run quick_llm_test.py (component testing)
- [ ] Record accuracy metrics
- [ ] Identify top 3 failure modes
- [ ] Test 5 edge cases manually
- [ ] Document unexpected behaviors
- [ ] Propose prompt improvements
- [ ] Re-test with updated prompts
- [ ] Measure performance improvement
- [ ] Update workflow documentation

## Notes

- All tests use `deepseek-v3.2-exp` model (default)
- Temperature set to 0.0 for deterministic results
- Registry path: `../data/reagents` (relative to data-processor/)
- Test scripts include rate limiting (0.5s delays) to avoid API throttling
