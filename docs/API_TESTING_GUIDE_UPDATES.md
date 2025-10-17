# API Testing Guide - Updates Summary

## What's New

The API Testing Guide has been comprehensively updated with:

### 1. **Catalyst Filtering Support** ✅
- Added catalyst filter examples for all applicable endpoints
- Demonstrates `relax.catalyst_class` parameter usage
- Covers Pd, Cu, Ni, Ru, Co, and other catalyst types

### 2. **Auto-Detection Testing** ✅
- Dedicated section on reaction type auto-detection
- Examples of omitting `reaction_type` for automatic detection
- Combined auto-detection + catalyst filtering examples
- Multi-reaction testing script to verify detection accuracy

### 3. **Enhanced Examples**
- **Example 1**: Compare different catalysts (Pd vs Cu vs Ni) for same reaction
- **Example 2**: Rule-based matching with catalyst filter
- **Example 3**: Auto-detection with catalyst preference
- **New Section 6**: Comprehensive auto-detection workflows

### 4. **Quick Reference Cards**
- Auto-detection usage patterns (3 methods)
- Catalyst filter syntax reference
- Common catalyst values table
- Reaction type values with auto-detect option
- Auto-detection behavior priority table
- Rerank strategies reference

### 5. **Common Mistakes Section**
- Wrong field names (reaction vs reaction_smiles)
- Missing `-Depth` parameter for nested objects
- Database/catalyst mismatch warnings

## Key Sections Added/Updated

### Section 2: Rule-Based Recommendation
**Before:** Basic usage only  
**After:** 
- Basic usage (no filter)
- With catalyst filter (NEW!)
- Multiple catalyst examples (Cu, Pd)

### Section 3: ML-Based Recommendation
**Before:** Basic usage only  
**After:**
- Basic usage (no filter)
- With catalyst filter (NEW!)
- Multiple catalyst examples (Pd, Cu)
- Combined filters (catalyst + rerank strategy)

### Section 5: Reaction Type Detection
**Before:** Simple detection endpoint  
**After:** Detailed explanation of detection response fields

### Section 6: Auto-Detection with Recommendations (NEW!)
- ML recommendations with auto-detection
- Auto-detection + catalyst filter
- Multi-reaction type testing script

### Catalyst Filtering Examples Section (NEW!)
- Example 1: Compare catalysts for C-N coupling
- Example 2: Rule-based with catalyst filter
- Example 3: Auto-detection + catalyst filter

## Testing Scenarios Now Covered

### ✅ Basic Testing
- Health check
- Rule-based matching
- ML-based recommendations
- Fusion (deprecated)
- Reaction type detection
- Protocol recommendations (CLI)

### ✅ Catalyst Filtering
- Filter by Pd, Cu, Ni, Ru, Co
- Combined with rule-based matching
- Combined with ML recommendations
- Multi-catalyst comparison

### ✅ Auto-Detection
- Omit reaction_type completely
- Set reaction_type to null
- Auto-detect with catalyst filter
- Verify detection accuracy across reaction types

### ✅ Advanced Workflows
- Auto-detect + catalyst + rerank
- Compare multiple catalysts
- Database selection based on catalyst
- Filter validation and metadata

## PowerShell Examples Added

### 1. Single Catalyst Filter
```powershell
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    relax = @{ catalyst_class = "Cu" }
} | ConvertTo-Json -Depth 5
```

### 2. Auto-Detection
```powershell
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    # No reaction_type - will auto-detect
} | ConvertTo-Json
```

### 3. Auto-Detection + Catalyst
```powershell
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    relax = @{ catalyst_class = "Pd" }
} | ConvertTo-Json -Depth 5
```

### 4. Multi-Catalyst Comparison
```powershell
# Test Pd, Cu, Ni sequentially
foreach ($cat in @("Pd", "Cu", "Ni")) {
    $body = @{
        reaction = "..."
        relax = @{ catalyst_class = $cat }
    } | ConvertTo-Json -Depth 5
    
    $result = Invoke-RestMethod ...
    Write-Host "Results for $cat`: $($result.recommendations.Count)"
}
```

### 5. Detection Verification
```powershell
# Test multiple reaction types
$reactions = @{
    "C-N Coupling" = "Brc1ccccc1.Nc1ccccc1>>..."
    "Suzuki" = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>..."
    "Amide" = "CC(=O)O.NCc1ccccc1>>..."
}

foreach ($type in $reactions.Keys) {
    $body = @{ reaction = $reactions[$type] } | ConvertTo-Json
    $result = Invoke-RestMethod ...
    Write-Host "$type detected as: $($result.meta.reaction_type)"
}
```

## Quick Reference Tables Added

### Catalyst Values Table
| Value | Description | Example Reactions |
|-------|-------------|-------------------|
| "Pd" | Palladium | Suzuki, Buchwald-Hartwig |
| "Cu" | Copper | Ullmann, Chan-Lam |
| "Ni" | Nickel | Negishi, Kumada |

### Reaction Types Table
| Value | Description | Auto-Detect |
|-------|-------------|-------------|
| "Suzuki" | C-C coupling | ❌ Manual |
| "C_N_Coupling" | C-N coupling | ❌ Manual |
| null / omitted | Auto-detect | ✅ Yes |

### Auto-Detection Priority Table
1. User-specified reaction_type (if provided)
2. ML-based (rxn-insight, if available)
3. Rule-based (chemtools.router)
4. Fallback to "Unknown"

## API Endpoint Summary Table Updated

| Endpoint | Catalyst Filter | Auto-Detect |
|----------|-----------------|-------------|
| `/match` | ✅ Yes | ❌ No |
| `/api/v1/recommend/conditions` | ✅ Yes | ✅ Yes |
| `/api/v1/recommend` | ✅ Yes | ✅ Yes |
| `/api/v1/reaction/detect-type` | N/A | ✅ Yes (purpose) |

## How to Use the Updated Guide

### For Quick Testing:
1. Go to "Test Commands (Copy & Paste)" section
2. Copy the relevant example
3. Adjust reaction SMILES if needed
4. Run in PowerShell

### For Catalyst Exploration:
1. Check "Catalyst Filtering Examples" section
2. Use Example 1 to compare all catalyst types
3. Pick the best catalyst based on results

### For Auto-Detection:
1. Check "Section 6: Auto-Detection with Recommendations"
2. Use the multi-reaction testing script
3. Verify detection accuracy for your use case

### For Advanced Workflows:
1. Combine auto-detection + catalyst filtering
2. Use rerank strategies for better results
3. Check metadata in response for detection details

## Testing Checklist

- [ ] Health check passes
- [ ] Rule-based matching works
- [ ] ML recommendations work
- [ ] Catalyst filtering works (Pd, Cu, Ni)
- [ ] Auto-detection works (omit reaction_type)
- [ ] Auto-detection + catalyst filter works
- [ ] Detection endpoint returns expected type
- [ ] Multiple catalysts can be compared
- [ ] Reranking strategies work
- [ ] OpenAPI docs accessible at /docs

## Files Updated

1. **docs/API_TESTING_GUIDE.md**
   - Added catalyst filtering sections
   - Added auto-detection sections
   - Added quick reference cards
   - Added common mistakes section
   - Updated all examples with new features

## Next Steps

### For Users:
1. Try the new catalyst filtering examples
2. Test auto-detection with your reactions
3. Compare different catalysts for your use case
4. Report any issues or unexpected behavior

### For Developers:
1. Ensure all endpoints support catalyst filtering
2. Verify auto-detection accuracy
3. Add more catalyst types if needed
4. Improve detection confidence reporting

## Related Documentation

- `CATALYST_ENDPOINT_SUPPORT.md` - Technical details on catalyst filtering
- `WEB_CLI_SYNC_COMPLETE.md` - CLI updates for catalyst support
- `app/main.py` - Endpoint implementations
- `chemtools/contracts.py` - Request/response models

## Support

If you encounter issues:
1. Check "Common Mistakes" section
2. Verify database files exist
3. Check server logs for errors
4. Test with `/docs` interactive API
5. Report bugs with example request/response

---

**Status:** ✅ Complete  
**Version:** Updated for catalyst filtering + auto-detection  
**Date:** October 17, 2025
