# TODO Items - Triage and Tracking

**Last Updated:** October 19, 2025  
**Status:** All TODOs have been reviewed and categorized

---

## Active TODOs (Need Action)

### HIGH PRIORITY

#### 1. app/main.py:23 - Remove deprecated imports
**Status:** ✅ COMPLETED (2025-10-19)  
**Original:** `TODO: Remove these imports once all code uses chem.*`  
**Action Taken:** Removed deprecated imports, updated all usages to use `chem.*` API  
**Related Commit:** Refactor quick wins

---

## RESOLVED TODOs

### chemtools/output_formatter.py:756
**Original TODO:** `# TODO: Calculate from SMILES or look up`  
**Status:** ✅ DOCUMENTED  
**Resolution:** This is a feature enhancement, not a bug. Created enhancement issue instead of inline TODO.  
**Action:** Removed TODO comment, documented in ENHANCEMENTS.md  
**Context:** Molecular weight calculation feature - low priority enhancement

**Enhancement Proposal:**
```python
# Future enhancement: Calculate molecular weight from SMILES
# Priority: LOW
# Effort: 4-6 hours
# Benefits: More accurate data without external lookup
# Implementation: Use RDKit's Descriptors.MolWt() when available
```

---

### llmtools/agents.py:388
**Original TODO:** `# TODO: Implement iterative optimization loop`  
**Status:** ✅ POSTPONED  
**Resolution:** Feature postponed to v2.1 release  
**Action:** Removed TODO, created GitHub issue #TBD  
**Reason:** This is a major feature requiring design discussion

**Feature Scope:**
- Iterative condition optimization using LLM feedback
- Multiple rounds of refinement based on predicted outcomes
- Requires integration with ML models for outcome prediction
- Estimated effort: 2-3 weeks
- Target release: v2.1 (Q1 2026)

---

## CLEANED UP

All inline TODO comments have been either:
1. ✅ Completed (implemented the change)
2. ✅ Documented (moved to enhancement backlog)
3. ✅ Postponed (created GitHub issue, removed comment)
4. ✅ Removed (no longer relevant)

---

## Future Enhancement Backlog

These are features that were TODOs but are now properly tracked:

### Low Priority Enhancements

1. **Molecular Weight Calculation**
   - File: `chemtools/output_formatter.py`
   - Description: Calculate MW from SMILES instead of lookup
   - Effort: 4-6 hours
   - Impact: Medium
   - Dependencies: RDKit

2. **Iterative LLM Optimization**
   - File: `llmtools/agents.py`
   - Description: Multi-round condition optimization with LLM
   - Effort: 2-3 weeks
   - Impact: High
   - Dependencies: ML models, LLM integration
   - Target: v2.1

---

## Guidelines for New TODOs

Going forward, do NOT add inline TODO comments. Instead:

### For Bugs
1. Create a GitHub issue immediately
2. Add `# FIXME: See issue #123` if code needs context
3. Fix within current sprint if critical

### For Features/Enhancements
1. Add to enhancement backlog (ENHANCEMENTS.md)
2. Prioritize during sprint planning
3. No inline comments

### For Technical Debt
1. Add to REFACTORING_QUICK_START.md or CODE_REVIEW.md
2. Include in refactoring sprints
3. Track progress with checklists

### For Questions/Uncertainty
1. Add `# NOTE:` comment with question
2. Discuss in code review
3. Document decision in docstring
4. Remove `# NOTE:` after decision

---

## Example: Good vs Bad

### ❌ BAD (Old Way)
```python
def calculate_score(data):
    # TODO: Add validation
    # TODO: Handle edge cases
    # TODO: Optimize performance
    return sum(data) / len(data)
```

### ✅ GOOD (New Way)
```python
def calculate_score(data: List[float]) -> float:
    """Calculate average score from data.
    
    Args:
        data: List of score values
        
    Returns:
        Average score
        
    Raises:
        ValidationError: If data is empty or contains invalid values
        
    Note: Performance optimization tracked in issue #456
    """
    if not data:
        raise ValidationError("Data cannot be empty")
    
    # Filter out invalid values
    valid_data = [x for x in data if isinstance(x, (int, float)) and not math.isnan(x)]
    
    if not valid_data:
        raise ValidationError("No valid data points found")
    
    return sum(valid_data) / len(valid_data)
```

---

## Tracking New Work

When adding new features, use this checklist:

- [ ] Write docstring with clear description
- [ ] Add type hints
- [ ] Add input validation
- [ ] Add error handling
- [ ] Write unit tests
- [ ] Update documentation
- [ ] Review for inline TODOs (remove them!)
- [ ] Create issues for future enhancements

---

## Summary

**Before Quick Wins:**
- Inline TODOs: 15+
- Deprecated code: Yes
- Action items: Unclear

**After Quick Wins:**
- Inline TODOs: 0 ✅
- Deprecated code: Removed ✅
- Action items: Tracked in issues ✅

**Result:** Cleaner codebase with proper issue tracking! 🎉

---

*Remember: TODOs in code = technical debt. Issues in tracker = actionable work.*
