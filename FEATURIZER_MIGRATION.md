# Featurizer Migration: ullmann → molecular

## Summary

The reaction-specific `ullmann.py` featurizer has been removed and consolidated into the general `molecular.py` featurizer. This change removes naming that was tied to a specific reaction type (Ullmann) when the features are actually general C-N coupling features applicable to both Ullmann and Buchwald-Hartwig reactions.

## Changes Made

### 1. **Removed Files**
- ❌ `chemtools/featurizers/ullmann.py` - deleted

### 2. **Updated Files**
- ✅ `chemtools/featurizers/molecular.py` - now contains full implementation
- ✅ `chemtools/featurizers/__init__.py` - removed ullmann import
- ✅ `chemtools/context.py` - changed `ullmann()` → `molecular()`
- ✅ `AGENTS.md` - updated documentation
- ✅ `docs/README.md` - updated references

### 3. **API Endpoints**
The API endpoint `/api/v1/featurize/ullmann` remains as a **deprecated alias** that redirects to `/api/v1/featurize/molecular` for backward compatibility. See `app/main.py` line 116.

## Migration Guide

### For Code

**Before:**
```python
from chemtools.featurizers.ullmann import featurize
# or
from chemtools.featurizers import ullmann

result = featurize("Brc1ccccc1", "Nc1ccccc1")
# or
result = ullmann.featurize("Brc1ccccc1", "Nc1ccccc1")
```

**After:**
```python
from chemtools.featurizers.molecular import featurize
# or
from chemtools.featurizers import molecular

result = featurize("Brc1ccccc1", "Nc1ccccc1")
# or
result = molecular.featurize("Brc1ccccc1", "Nc1ccccc1")
```

### For API Calls

**Before:**
```bash
POST /api/v1/featurize/ullmann
{
  "electrophile": "Brc1ccccc1",
  "nucleophile": "Nc1ccccc1"
}
```

**After (Recommended):**
```bash
POST /api/v1/featurize/molecular
{
  "electrophile": "Brc1ccccc1",
  "nucleophile": "Nc1ccccc1"
}
```

**Note:** The old endpoint still works but returns a deprecation header.

### For Context API

**Before:**
```python
from chemtools import chem

featurizer = chem.featurizers.ullmann()
result = featurizer.featurize("Brc1ccccc1", "Nc1ccccc1")
```

**After:**
```python
from chemtools import chem

featurizer = chem.featurizers.molecular()
result = featurizer.featurize("Brc1ccccc1", "Nc1ccccc1")
```

## Rationale

### Why Remove Reaction-Specific Featurizers?

1. **Misleading Name**: The features extracted (LG type, electrophile class, nucleophile basicity, etc.) are general C-N coupling features, not Ullmann-specific.

2. **Applicable to Multiple Reactions**: The same features work for:
   - Ullmann reactions (Cu-catalyzed C-N coupling)
   - Buchwald-Hartwig reactions (Pd-catalyzed C-N coupling)
   - Nickel-catalyzed C-N coupling
   - Any aryl halide + nucleophile coupling

3. **Reduces Confusion**: Having both `ullmann.py` and `molecular.py` with the same functionality was confusing.

4. **Better Abstraction**: `molecular` is a more appropriate name for substrate featurization.

## Feature Set

The `molecular.featurize()` function extracts:

### Electrophile Features
- **LG** (Leaving Group): Br, Cl, I, OTf, UNK
- **elec_class**: aryl, vinyl, alkyl
- **ortho_count**: "0", "1", "2+"
- **para_EWG**: boolean (para electron-withdrawing group)
- **heteroaryl**: boolean

### Nucleophile Features
- **nuc_class**: aniline, indole, amine_primary, amine_secondary, phenol, amide_deactivated
- **n_basicity**: aromatic_primary, aliphatic_primary, secondary, deactivated
- **steric_alpha**: low, med, high

### Combined
- **bin**: Coarse bin key like "LG:Br|NUC:aniline"

### Optional (when CHEMTOOLS_ATTACH_ROLE_AWARE=1)
- **role_aware**: Role-specific feature vectors for electrophile and nucleophile

## Testing

```bash
# Test the import
python -c "from chemtools.featurizers import molecular; print('✓ Import works')"

# Test featurization
python -c "
from chemtools.featurizers import molecular
result = molecular.featurize('Brc1ccccc1', 'Nc1ccccc1')
print('LG:', result['LG'])
print('bin:', result['bin'])
print('✓ Featurization works')
"

# Test API endpoint
curl -X POST http://localhost:8000/api/v1/featurize/molecular \
  -H "Content-Type: application/json" \
  -d '{"electrophile": "Brc1ccccc1", "nucleophile": "Nc1ccccc1"}'
```

## Backward Compatibility

### API Endpoint
The `/api/v1/featurize/ullmann` endpoint is maintained as a **deprecated alias** with:
- Warning logged to server logs
- `X-Deprecated: true` header
- `Link: </api/v1/featurize/molecular>; rel="successor-version"` header

### Python Code
No backward compatibility shim is provided for direct Python imports. Code must be updated to use `molecular` instead of `ullmann`.

## Timeline

- **October 15, 2025**: Migration completed, `ullmann.py` removed
- **Deprecation Period**: API endpoint `/api/v1/featurize/ullmann` will be maintained for 6 months
- **End of Life**: April 2026 - deprecated endpoint may be removed

## Questions?

If you have questions about this migration or need assistance updating your code, please:
1. Check the examples above
2. Review `chemtools/featurizers/molecular.py` for implementation details
3. Refer to API documentation at `/docs`

## Related Files

- `chemtools/featurizers/molecular.py` - Main implementation
- `chemtools/featurizers/__init__.py` - Module exports
- `app/main.py` - API endpoints (lines 116-129)
- `chemtools/context.py` - Context API (line 549)
- `chemtools/integrations/mcp/tools/featurize.py` - MCP integration
