# ML Recommendation Error Handling Guide

## Overview

The ML recommendation system now provides comprehensive error handling with helpful, actionable guidance when:
- Reaction type cannot be detected
- Detected type is not supported by ML models
- No precedents are found in the database
- System errors occur during processing

## Error Scenarios

### 1. Unsupported Reaction Type

**When it happens:**
- Auto-detection identifies a reaction type not in the ML training dataset
- User manually selects an unsupported type

**Error Message:**
```
**Cannot Proceed with ML Recommendations**

🔍 **Detection Result:**
**Auto-detected (rxn-insight):** Esterification
  Class: Acylation
  Confidence: 75.00%

❌ **No ML model available** for: `Esterification`

**Available ML reaction types:**
  - C-N Coupling (Cu) (`C_N_Coupling_Cu`)
  - C-N Coupling (Pd) (`C_N_Coupling_Pd`)
  - C-N Coupling (Ni) (`C_N_Coupling_Ni`)
  - Amide Formation (`Amide_Coupling`)
  - Suzuki Coupling (`Suzuki_CC`)

**What to do:**
1. ✅ **Try rule-based recommendations** instead
2. 🔄 **Manually select** a supported reaction type from the dropdown
3. 📝 **Verify your SMILES** represents a supported reaction
4. 📖 **Check if your reaction** is similar to a supported type
```

**Resolution Steps:**
1. Click "Get Rule Recommendations" button instead
2. If reaction is similar to a supported type, manually select it
3. Verify SMILES format is correct
4. Consider if your reaction can be modeled as a supported type

---

### 2. No Precedents Found

**When it happens:**
- Reaction type is supported but no similar precedents exist in database
- Substrates are too unusual or complex
- Functional groups don't match training data

**Error Message (Supported Type):**
```
**No ML Recommendations Found**

🔍 **Auto-Detection Result:**
**Auto-detected (rxn-insight):** Buchwald_CN (Pd-catalyzed)
  Class: Heteroatom Alkylation and Arylation
  Catalysts: Pd
  Confidence: 85.00%

**Reaction Type:** Buchwald_CN (Pd-catalyzed)
**ML Family:** C_N_Coupling_Pd
**Detection Confidence:** 85.0%

**Why this happened:**

✅ Reaction type `C_N_Coupling_Pd` is supported, but no precedents were found.

**Possible reasons:**
1. 📊 **No similar reactions** in the precedent database
2. 🔬 **Unusual substrates** or functional groups
3. 💾 **Dataset not loaded** (first-run issue)
4. 🎯 **Detection error** - Wrong family detected

**What to do:**
1. ✅ **Try rule-based recommendations** (click 'Get Rule Recommendations')
2. 🔄 **Manually select** the correct reaction type if auto-detect was wrong
3. 🧪 **Simplify substrates** - Try a simpler model reaction first
4. 📖 **Check literature** for similar transformations

💡 **Note:** 12 precedents were found but none matched your substrates closely enough.
   Try adjusting the reaction SMILES or selecting a different reaction type.
```

**Resolution Steps:**
1. **Try rule-based recommendations** - May still provide useful conditions
2. **Simplify substrates** - Test with simpler model substrates first
3. **Check detection** - Verify auto-detected type is correct
4. **Manual selection** - Try related reaction types
5. **Literature search** - Look for similar published reactions

---

### 3. Low Detection Confidence Warning

**When it happens:**
- Auto-detection confidence < 50%
- Reaction pattern is ambiguous
- Multiple possible reaction types

**Console Warning:**
```
⚠️ WARNING: Low detection confidence (42.0%)
   Consider manually selecting the reaction type for better results.
```

**What it means:**
- Detection may be incorrect
- Results may not be reliable
- Manual selection recommended

**Resolution Steps:**
1. Review your SMILES for errors
2. Manually select the reaction type from dropdown
3. Verify functional groups match expected reaction
4. Consider if reaction fits into standard categories

---

### 4. System Errors

**When it happens:**
- Dataset loading fails
- SMILES parsing error
- Memory/resource issues
- Unexpected exceptions

**Error Message:**
```
**ML Recommendation System Error**

❌ **Error Type:** `FileNotFoundError`
**Message:** Dataset file not found: data/precedents/cn_coupling_pd.jsonl

🔍 **Detection Info:**
**Auto-detected (rxn-insight):** Buchwald_CN (Pd-catalyzed)
  Class: Heteroatom Alkylation and Arylation
  Catalysts: Pd

**Reaction Type:** C_N_Coupling_Pd
**SMILES:** Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1

**Common Issues & Solutions:**

1. 💾 **Dataset Loading Issue**
   - The ML precedent database may not be loaded
   - Try running a simple test reaction first to warm up the system
   - Check that dataset files exist in `data/` directory

**What to try:**
1. ✅ **Use rule-based recommendations** instead (click 'Get Rule Recommendations')
2. 🔄 **Manually select** the reaction type instead of auto-detect
3. 📝 **Simplify your SMILES** or try a simpler model reaction
4. 🔁 **Restart the application** if this persists

**Technical Details:**
```
Traceback (most recent call last):
  File "chemtools/recommend.py", line 123, in recommend_from_reaction
    dataset = load_dataset(family)
FileNotFoundError: Dataset file not found: data/precedents/cn_coupling_pd.jsonl
```
```

**Resolution Steps:**
1. **Check dataset files** - Ensure data files are present
2. **Restart application** - May resolve loading issues
3. **Use rule-based** - Fallback to rule-based recommendations
4. **Report bug** - If error persists, report with details

---

## Error Types by Category

### Detection Failures

| Error | Cause | Solution |
|-------|-------|----------|
| Cannot detect type | Unusual reaction | Manual selection |
| Low confidence | Ambiguous pattern | Manual selection |
| Wrong type detected | Similar patterns | Manual override |

### Data Availability

| Error | Cause | Solution |
|-------|-------|----------|
| Unsupported family | Not in training data | Try rule-based |
| No precedents | Unusual substrates | Simplify or manual search |
| Dataset missing | Installation issue | Check data/ directory |

### System Errors

| Error | Cause | Solution |
|-------|-------|----------|
| SMILES parsing | Invalid syntax | Verify SMILES format |
| Memory error | Large dataset | Reduce top_k parameter |
| File not found | Missing data | Check installation |

---

## Best Practices

### For Users

1. **Start Simple**
   - Test with known, simple reactions first
   - Gradually increase complexity

2. **Verify SMILES**
   - Use standard reaction SMILES format: `reactants>>products`
   - Include catalysts when relevant: `reactants>catalyst>products`

3. **Check Detection**
   - Review auto-detected type before proceeding
   - Manually select if confidence is low (<60%)

4. **Use Fallbacks**
   - Try rule-based recommendations when ML fails
   - Consult literature for unusual transformations

5. **Provide Context**
   - Include relevant functional groups in SMILES
   - Use realistic substrates

### For Developers

1. **Error Messages Should:**
   - Explain what happened clearly
   - Provide actionable next steps
   - Offer alternatives (rule-based)
   - Include technical details when relevant

2. **Always Include:**
   - Detection confidence scores
   - Available alternatives
   - Troubleshooting steps
   - Fallback options

3. **Handle Gracefully:**
   - Never show raw stack traces to users
   - Provide context-specific guidance
   - Suggest rule-based as fallback
   - Log errors for debugging

---

## Testing Error Handling

Run the comprehensive error handling test:

```powershell
$env:PYTHONIOENCODING='utf-8'
python scripts/test_ml_error_handling.py
```

This tests:
- Unsupported reaction types
- No precedents scenarios
- Valid reactions (control)
- Error message formatting

---

## Configuration

### Adjusting Behavior

**Disable auto-detection:**
```python
# Always require manual selection
ML_FAMILY_MAP = {
    # Remove "Auto-detect": None
    "C-N Coupling (Pd)": "C_N_Coupling_Pd",
    ...
}
```

**Set confidence threshold:**
```python
# In detect_and_map_reaction_type()
MIN_CONFIDENCE = 0.6  # Require 60% confidence
if detection_confidence < MIN_CONFIDENCE:
    return {
        "success": False,
        "message": f"Detection confidence too low ({detection_confidence:.1%})"
    }
```

**Customize error messages:**
Edit the error message templates in `get_ml_recommendations()` function.

---

## Future Enhancements

Planned improvements:

1. **Confidence Thresholds**
   - Block recommendations below confidence threshold
   - Suggest manual selection for borderline cases

2. **Substrate Similarity**
   - Show closest matching precedents even when not recommending
   - Suggest modifications to improve matching

3. **Alternative Types**
   - Suggest related reaction types to try
   - Show similarity scores across families

4. **User Feedback**
   - Allow users to report incorrect detections
   - Learn from corrections over time

5. **Batch Processing**
   - Handle multiple failures gracefully
   - Provide summary of successes/failures

---

## Summary

The ML error handling system provides:

✅ **Clear explanations** - Users understand what went wrong  
✅ **Actionable guidance** - Specific steps to resolve issues  
✅ **Context awareness** - Error messages adapt to situation  
✅ **Fallback options** - Always suggest rule-based alternative  
✅ **Technical details** - Developers can debug issues  
✅ **User-friendly** - No raw stack traces or jargon  

This comprehensive error handling ensures users can always make progress, even when ML recommendations fail, by providing clear paths forward through rule-based recommendations, manual selection, or troubleshooting steps.
