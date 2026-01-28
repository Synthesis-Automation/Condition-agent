# Logical Improvement Plan for Reaction Detection

## Problem Summary

Current detection runs **before** calculating reacted_motifs, so it can't use the most reliable information (actual motif consumption). This causes systematic misclassifications for cross-coupling reactions.

## Solution: Two-Pass Detection with Pattern Validation

### Phase 1: Slot-Based Detection (Existing)
- Fast, works without products
- Uses motif presence in reactants
- Can have false positives (e.g., Ar-H in boronic acids)

### Phase 2: Pattern Validation (New)
- Uses reacted_motifs + formed_motifs
- High-confidence pattern matching
- Corrects known false positives
- Provides audit trail

## Key Principles

1. **Non-Breaking**: Slot-based detection still runs first
2. **Logical Patterns**: Use proven reactant→product transformations
3. **High Confidence**: Only override when pattern is unambiguous
4. **Transparent**: Track original detection and correction reason
5. **Extensible**: Easy to add new patterns

## Implementation Logic

### Pattern Matching Rules

Each pattern follows this logic:
```
IF (reacted_motifs contains [Group A] AND [Group B])
   AND (formed_motifs contains [Expected Product])
   AND (initial_detection != Expected_Type)
THEN correct to Expected_Type with high confidence
```

### Common Patterns

#### 1. Suzuki-Miyaura
```
Reacted: (Ar-B(OH)2 OR Ar-B(OR)2) AND (Ar-Br OR Ar-Cl OR Ar-I)
Formed: Ar-Ar OR Ar-Alkenyl
→ Suzuki_miyaura (confidence: 0.95)
```

#### 2. Buchwald-Hartwig (C-N Coupling)
```
Reacted: (Ar-Br OR Ar-Cl OR Ar-I) AND (Ar-NH2 OR RNH2 OR Ar-NHR)
Formed: Ar-NHR OR Ar-NR2
→ C_N_Coupling (confidence: 0.95)
```

#### 3. Stille Coupling
```
Reacted: (Ar-Br OR Ar-Cl OR Ar-I) AND (Ar-SnR3 OR Alkenyl-SnR3)
Formed: Ar-Ar
→ Stille (confidence: 0.95)
```

#### 4. Negishi Coupling
```
Reacted: (Ar-Br OR Ar-Cl OR Ar-I) AND (Ar-ZnX OR Alkyl-ZnX)
Formed: Ar-Ar
→ Negishi (confidence: 0.95)
```

#### 5. Heck Reaction
```
Reacted: (Ar-Br OR Ar-Cl OR Ar-I) AND Alkenyl
Formed: Ar-Alkenyl
→ Heck (confidence: 0.95)
```

#### 6. Sonogashira Coupling
```
Reacted: (Ar-Br OR Ar-Cl OR Ar-I) AND (Alkynyl OR Ar-Alkynyl)
Formed: Ar-Alkynyl
→ Sonogashira (confidence: 0.95)
```

### Exclusion Rules

#### Arylation_Ar_H False Positive
```
IF initial_detection == "Arylation_Ar_H"
   AND reacted_motifs contains organometallic nucleophile
   (Ar-B(OH)2, Ar-SnR3, Ar-ZnX, etc.)
THEN mark as "Unknown" and let role classification handle it
```

This prevents C-H activation from being called when there's clearly a cross-coupling nucleophile present.

## Integration Points

### File: `chemtools/featurizers/formatters/reaction.py`

**Location: After aggregates calculation (around line 260)**

```python
# Current code
aggregates = aggregate_reaction_features(
    reactant_bundles,
    product_motif_ids=product_motif_ids,
)

# NEW: Validate detection with reacted motifs
from .detection_validation import validate_detection_with_reacted_motifs

validated = validate_detection_with_reacted_motifs(
    initial_detection=reaction_type,
    initial_confidence=detection.confidence,
    reacted_motifs=aggregates.get('reacted_motifs', []),
    formed_motifs=aggregates.get('formed_motifs', []),
    spectator_motifs=aggregates.get('spectator_motifs', []),
)

# Update reaction type if corrected
if validated.get('corrected_from'):
    reaction_type = validated['reaction_type']
    # Store validation metadata
    detection_payload['validation'] = {
        'original_detection': validated['corrected_from'],
        'validated_detection': validated['reaction_type'],
        'validation_method': validated['validation_method'],
        'validation_reason': validated['reason'],
        'validation_confidence': validated['confidence'],
    }
```

## Benefits

### Accuracy Improvements
- **Suzuki**: Fixes false "Arylation_Ar_H" (common)
- **Buchwald-Hartwig**: Fixes false "Arylation_Ar_H" (common)
- **Other cross-couplings**: Systematic correction

### Transparency
- Original detection preserved
- Correction reason documented
- Confidence scores maintained
- Audit trail for debugging

### Maintainability
- Pattern rules are clear and explicit
- Easy to add new patterns
- No changes to core detection system
- Backward compatible

## Testing Strategy

### Unit Tests
```python
def test_suzuki_correction():
    result = validate_detection_with_reacted_motifs(
        initial_detection="Arylation_Ar_H",
        initial_confidence=1.0,
        reacted_motifs=["Ar-B(OH)2", "Ar-Br"],
        formed_motifs=["Ar-Ar"],
        spectator_motifs=[],
    )
    assert result['reaction_type'] == "Suzuki_miyaura"
    assert result['corrected_from'] == "Arylation_Ar_H"
```

### Integration Tests
```python
def test_end_to_end_suzuki():
    payload = featurize_reaction("c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1")
    assert payload['reaction_type'] == "Suzuki_miyaura"
    # Should have validation metadata
    assert 'validation' in payload.get('detection', {})
```

### Regression Tests
- Ensure legitimate Arylation_Ar_H reactions still work
- Verify no false corrections
- Check confidence scores are appropriate

## Rollout Plan

### Phase 1: Implementation (1-2 days)
1. Create `detection_validation.py` module
2. Implement pattern matching logic
3. Add unit tests

### Phase 2: Integration (1 day)
1. Integrate into `featurize_reaction()`
2. Update CLI to show validation info
3. Add integration tests

### Phase 3: Validation (1-2 days)
1. Run on test dataset
2. Analyze corrections
3. Tune patterns if needed

### Phase 4: Documentation (1 day)
1. Update API docs
2. Add usage examples
3. Document pattern rules

## Expected Impact

### Quantitative
- Estimated 10-20% reduction in cross-coupling misclassifications
- 95%+ accuracy for major cross-coupling families
- Zero regression on existing correct detections

### Qualitative
- Conditions recommendation more accurate
- Better user trust in system
- Easier debugging (validation metadata)
- Foundation for future improvements

## Future Enhancements

### Short-term
- Add more cross-coupling patterns (Kumada, Negishi variants)
- Support C-O coupling patterns
- Add C-S coupling patterns

### Long-term
- ML-based pattern learning
- Confidence calibration using historical data
- Multi-step reaction detection
- Regioselectivity prediction

## Conclusion

This two-pass approach is **logically sound** because:

1. **Uses the right information**: Reacted motifs are the ground truth
2. **Non-destructive**: Original detection preserved for analysis
3. **Transparent**: Clear reasoning for corrections
4. **Extensible**: Easy to add patterns as we learn more
5. **Practical**: Solves real user pain points immediately

The key insight is that we don't need to rewrite detection from scratch—we just add a smart validation layer that catches systematic errors using proven chemical logic.
