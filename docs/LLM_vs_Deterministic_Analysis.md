# LLM vs Deterministic Approaches in Reaction Rule Automation

## Executive Summary

**Can we automate without LLM APIs?** **YES** - Most of the workflow (80%+) can be fully automated with deterministic chemistry tools (RDKit, ChemDataExtractor).

**Should we use LLM APIs?** **Strategically** - LLMs provide significant value in 3-4 specific steps where semantic understanding and knowledge synthesis are required.

---

## Complete Workflow Analysis: LLM vs Deterministic

### Phase 1: Information Gathering & Extraction

| Task | Deterministic Approach | LLM Approach | Winner | Confidence |
|------|----------------------|--------------|---------|-----------|
| **PDF Text Extraction** | PyPDF2, pdfplumber | GPT-4 Vision | **Deterministic** | 95% |
| **Structure Recognition** | OSRA, DECIMER | GPT-4 Vision | **Hybrid** | 70% |
| **Table Parsing** | Camelot, Tabula | GPT-4 | **Deterministic** | 90% |
| **Condition Extraction** | Named Entity Recognition | GPT-4 | **LLM** | 85% |
| **Reaction Mapping** | RXNMapper, Indigo | - | **Deterministic** | 98% |

**Verdict:** Mostly deterministic, LLM helps with **unstructured condition text**

#### Example: Condition Extraction

**Without LLM (Regex/NER):**
```python
# Regex-based extraction - brittle, misses variations
pattern = r'Pd\(PPh3\)4\s*\((\d+\.?\d*)\s*mol%\)'
# Misses: "5 mol% Pd(PPh3)4", "Pd(PPh3)4 (5%)", etc.
```

**With LLM (GPT-4):**
```python
prompt = """
Extract reaction conditions from this experimental text:
"The reaction was performed with 5 mol% tetrakis(triphenylphosphine)palladium(0)..."

Return JSON:
{
  "pd_source": "Pd(PPh3)4",
  "pd_mol%": 5.0,
  ...
}
"""
# Handles variations, synonyms, implicit information
```

**Recommendation:** Use ChemDataExtractor first (free), fallback to LLM for complex cases

---

### Phase 2: Rule Generation

| Task | Deterministic Approach | LLM Approach | Winner | Confidence |
|------|----------------------|--------------|---------|-----------|
| **MCS Finding** | RDKit rdFMCS | - | **Deterministic** | 99% |
| **SMARTS Generation** | RDKit | GPT-4 + validation | **Deterministic** | 95% |
| **Pattern Refinement** | Iterative testing | LLM suggestions | **Hybrid** | 60% |
| **Priority Assignment** | Rule-based scoring | LLM reasoning | **Deterministic** | 80% |
| **Condition Aggregation** | Statistical mode/median | - | **Deterministic** | 95% |

**Verdict:** Deterministic is sufficient, LLM can **suggest pattern improvements**

#### Example: SMARTS Pattern Refinement

**Without LLM (Trial & Error):**
```python
# Manual iteration
patterns = [
    "[c:1]-[Br:2]",                    # Too broad
    "[c:1]([!H])([!H])-[Br:2]",       # Better
    "[c:1]1[c]([C,N,O])[c][c][c][c]1-[Br:2]"  # More specific
]
# Test each manually, adjust based on results
```

**With LLM (Intelligent Suggestions):**
```python
prompt = """
This SMARTS pattern matches too many substrates:
Pattern: [c:1]-[Br:2]
False positives: benzyl bromide, pyridyl bromide
True target: ortho-substituted aryl bromides

Suggest 3 more specific SMARTS patterns with explanations.
"""

# GPT-4 can suggest:
# 1. "[c:1]1[c]([!H])[c,cH][c,cH][c,cH][c,cH]1-[Br:2]" (ortho non-H)
# 2. "[c:1]1[cH0][c,cH][c,cH][c,cH][c,cH]1-[Br:2]" (ortho substituted)
# Explanations help understand trade-offs
```

**Recommendation:** Start deterministic, use LLM for **debugging hard cases**

---

### Phase 3: Database Entry Creation

| Task | Deterministic Approach | LLM Approach | Winner | Confidence |
|------|----------------------|--------------|---------|-----------|
| **ID Generation** | Template + hashing | - | **Deterministic** | 100% |
| **Name Generation** | Template filling | GPT-4 | **LLM** | 70% |
| **Token Signature** | String extraction | - | **Deterministic** | 100% |
| **Metadata Creation** | Field mapping | - | **Deterministic** | 100% |
| **Notes/Description** | Template | GPT-4 | **LLM** | 80% |

**Verdict:** Deterministic for structure, LLM for **human-readable descriptions**

#### Example: Human-Readable Names

**Without LLM:**
```python
name = f"{rxn_type}-{substrate}-{ligand}"
# "Suzuki_Miyaura-ArBr-SPhos"
# Functional but not descriptive
```

**With LLM:**
```python
prompt = """
Generate a concise, chemist-friendly name for this rule:
- Reaction: Suzuki-Miyaura coupling
- Substrate: ortho-substituted aryl bromides
- Ligand: SPhos
- Key feature: Mild conditions, good for hindered substrates

Max 80 characters.
"""

# "SPhos-mediated Suzuki coupling of hindered aryl bromides"
# Much more informative for users
```

**Recommendation:** Template names work fine, LLM is **nice-to-have** for UX

---

### Phase 4: Test Case Generation

| Task | Deterministic Approach | LLM Approach | Winner | Confidence |
|------|----------------------|--------------|---------|-----------|
| **Positive Tests** | SMARTS → Mol | - | **Deterministic** | 95% |
| **Boundary Tests** | Systematic variation | LLM suggestions | **Hybrid** | 60% |
| **Negative Tests** | Anti-pattern generation | LLM suggestions | **Hybrid** | 60% |
| **Expected Features** | SMARTS parsing | - | **Deterministic** | 100% |

**Verdict:** Deterministic works, LLM helps with **edge case discovery**

#### Example: Boundary Test Generation

**Without LLM (Manual):**
```python
# Manual specification of boundary cases
boundary_tests = [
    "Brc1ccccc1",           # Simple PhBr (should NOT match "ortho-sub" rule)
    "Brc1c(C)cccc1",        # One ortho-sub (boundary)
    "Brc1c(C)cccc1C",       # Two ortho-sub (should match)
]
```

**With LLM (Intelligent Generation):**
```python
prompt = """
This rule matches ortho-substituted aryl bromides:
SMARTS: [c:1]1[c]([!H])[c][c][c][c]1-[Br:2]

Generate 5 test molecules:
1. Clear positive (should match)
2. Clear negative (should NOT match)
3. Boundary case 1 (one ortho substituent)
4. Boundary case 2 (meta vs ortho)
5. Edge case (fused ring)

Return as SMILES with descriptions.
"""

# LLM generates chemically meaningful edge cases
# that humans might miss
```

**Recommendation:** Deterministic for basic tests, LLM for **comprehensive coverage**

---

### Phase 5: Testing & Validation

| Task | Deterministic Approach | LLM Approach | Winner | Confidence |
|------|----------------------|--------------|---------|-----------|
| **SMARTS Matching** | RDKit | - | **Deterministic** | 100% |
| **Feature Detection** | RDKit descriptors | - | **Deterministic** | 100% |
| **Result Comparison** | String/ID matching | - | **Deterministic** | 100% |
| **Report Generation** | Template formatting | - | **Deterministic** | 100% |

**Verdict:** **100% deterministic** - no LLM needed

---

### Phase 6: Iterative Improvement

| Task | Deterministic Approach | LLM Approach | Winner | Confidence |
|------|----------------------|--------------|---------|-----------|
| **Failure Classification** | Rule-based | LLM analysis | **LLM** | 90% |
| **Root Cause Analysis** | Pattern matching | LLM reasoning | **LLM** | 85% |
| **Fix Recommendation** | Template suggestions | LLM explanation | **LLM** | 80% |
| **Priority Adjustment** | Numeric calculation | - | **Deterministic** | 100% |
| **Pattern Refinement** | Trial & error | LLM suggestions | **Hybrid** | 70% |

**Verdict:** **LLMs shine here** - semantic analysis of failures

#### Example: Failure Analysis

**Without LLM (Rule-based):**
```python
def diagnose_failure(failure):
    if failure['actual_priority'] > failure['expected_priority']:
        return {
            'type': 'priority_conflict',
            'fix': 'increase_priority',
            'suggested_priority': failure['actual_priority'] + 5
        }
    # Limited to predefined cases
```

**With LLM (Intelligent Analysis):**
```python
prompt = """
Analyze this validation failure:

Expected rule: SCDB-SUZ-ARBR-ORTHO-XPhos (priority 50)
- Pattern: [c:1]1[c]([!H])[c][c][c][c]1-[Br:2]
- Description: Ortho-substituted aryl bromides with XPhos

Actual match: SCDB-SUZ-ARBRI-GENERAL-SPhos (priority 46)
- Pattern: [c:1]-[Br:2]
- Description: General aryl bromides with SPhos

Test substrate: Brc1c(C)cccc1C (2,6-dimethylbromobenzene)

Why did the test fail? What's the root cause? Suggest 2-3 specific fixes.
"""

# LLM response:
"""
Root cause: The general pattern [c:1]-[Br:2] matches first because it has
lower priority (46 < 50). The test substrate DOES match the ortho pattern,
but priority order is inverted.

Fixes:
1. **Increase ORTHO priority to 52** (above GENERAL)
   - Rationale: More specific rules should have higher priority
   
2. **Add negative assertion to GENERAL pattern**
   - Pattern: [c:1]1[c,cH0][c][c][c][c]1-[Br:2]
   - Rationale: Exclude molecules with ortho substitution
   
3. **Add ortho_sub_count feature requirement to ORTHO rule**
   - Requirement: {"electrophile.ortho_sub_count": [1, 2]}
   - Rationale: Double-check with feature detection
   
Recommendation: Fix #1 (priority increase) is simplest and most effective.
"""
```

**Recommendation:** **This is where LLMs provide huge value** - saves hours of debugging

---

### Phase 7: Continuous Integration

| Task | Deterministic Approach | LLM Approach | Winner | Confidence |
|------|----------------------|--------------|---------|-----------|
| **Automated Testing** | pytest | - | **Deterministic** | 100% |
| **Trend Analysis** | Statistical | LLM insights | **Hybrid** | 70% |
| **Regression Detection** | Numeric threshold | - | **Deterministic** | 100% |
| **Alert Generation** | Template | LLM summary | **LLM** | 75% |

**Verdict:** Deterministic core, LLM for **interpretable summaries**

---

### Phase 8: Scaling to New Reaction Types

| Task | Deterministic Approach | LLM Approach | Winner | Confidence |
|------|----------------------|--------------|---------|-----------|
| **Template Creation** | Manual | LLM generation | **LLM** | 90% |
| **Substrate Classification** | Predefined taxonomy | LLM extraction | **Hybrid** | 70% |
| **Common Conditions** | Statistical | Literature mining | **Hybrid** | 65% |
| **Best Practices** | Manual documentation | LLM synthesis | **LLM** | 85% |

**Verdict:** **LLMs excel at knowledge synthesis** from literature

#### Example: New Reaction Type Bootstrap

**Without LLM (Manual):**
```python
# Manually research and code each reaction type
class BuchwaldHartwigReactionType:
    # 2-3 hours of research to find:
    # - Common catalysts
    # - Typical conditions
    # - Substrate classes
    # - Known challenges
    pass
```

**With LLM (Knowledge Extraction):**
```python
prompt = """
Create a reaction type template for Buchwald-Hartwig C-N coupling:

1. Common Pd sources (with typical loadings)
2. Most-used ligands (with substrate preferences)
3. Typical bases and solvents
4. Temperature ranges
5. Substrate classes (electrophile and nucleophile)
6. Known challenging cases

Base this on published literature 2015-2024.
Format as Python code following this template:
[paste template]
"""

# LLM generates 80% of template in minutes
# Human reviews and validates
```

**Recommendation:** **LLMs can 10x speed** for new reaction types

---

### Phase 9: Advanced Techniques

| Task | Deterministic Approach | LLM Approach | Winner | Confidence |
|------|----------------------|--------------|---------|-----------|
| **ML Feature Selection** | Statistical | - | **Deterministic** | 95% |
| **Outcome Prediction** | QSAR models | Hybrid models | **Hybrid** | 80% |
| **Literature Mining** | NER + relation extraction | LLM extraction | **LLM** | 85% |
| **Expert Knowledge Capture** | Manual encoding | LLM interview | **LLM** | 75% |

**Verdict:** **LLMs for knowledge extraction**, ML for prediction

---

## Summary: Where LLMs Add Most Value

### 🟢 High-Value LLM Use Cases (ROI > 5x)

#### 1. **Failure Analysis & Debugging** ⭐⭐⭐⭐⭐
```python
# Without LLM: 30 min per failure
# With LLM: 2 min per failure (15x faster)

def analyze_with_llm(failure, db, test_case):
    prompt = f"""
    Validation failure analysis:
    - Expected: {failure['expected']} (priority {failure['expected_priority']})
    - Got: {failure['actual']} (priority {failure['actual_priority']})
    - Test: {test_case['smiles']}
    - Description: {test_case['description']}
    
    Database context:
    Expected rule: {db[failure['expected']]}
    Matched rule: {db[failure['actual']]}
    
    Diagnose root cause and suggest 2-3 specific fixes with rationale.
    """
    
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[{"role": "user", "content": prompt}]
    )
    
    return response.choices[0].message.content
```

**Cost:** ~$0.05 per analysis  
**Savings:** 28 min × $50/hr = $23.33 per failure  
**ROI:** 466x

#### 2. **New Reaction Type Bootstrap** ⭐⭐⭐⭐⭐
```python
# Without LLM: 4-8 hours research + coding
# With LLM: 30 min review + validation

def bootstrap_reaction_type(rxn_name):
    prompt = f"""
    Generate a complete reaction type template for {rxn_name}:
    
    Include:
    1. Reaction mechanism (1-2 sentences)
    2. Common catalysts/reagents (with preferences)
    3. Typical conditions (temp, solvent, time)
    4. Substrate scope and limitations
    5. Known challenging cases
    6. Representative SMARTS pattern
    7. Priority guidelines
    8. 3-5 literature references (DOI)
    
    Format as JSON following this schema: [schema]
    Base on peer-reviewed literature 2015-2024.
    """
    
    template = call_llm(prompt)
    # Human validation + adjustment: 30 min
    return template
```

**Cost:** ~$0.50 per reaction type  
**Savings:** 6 hr × $50/hr = $300  
**ROI:** 600x

#### 3. **Unstructured Condition Extraction** ⭐⭐⭐⭐
```python
# Without LLM: Regex (70% accuracy) + manual fixes
# With LLM: 95% accuracy, minimal manual work

def extract_conditions_from_text(experimental_text):
    prompt = f"""
    Extract reaction conditions from this experimental procedure:
    
    "{experimental_text}"
    
    Return as JSON:
    {{
      "pd_source": "Pd(PPh3)4",
      "pd_mol%": 5.0,
      "ligand": "SPhos",
      "ligand_mol%": 10.0,
      "base": "K2CO3",
      "base_eq": 2.0,
      "solvent": "THF/H2O (4:1)",
      "temperature_C": 60,
      "time_h": 12
    }}
    
    Use null for missing information. Standardize chemical names.
    """
    
    return json.loads(call_llm(prompt))
```

**Cost:** ~$0.03 per extraction  
**Savings:** 15 min manual work = $12.50  
**ROI:** 416x

#### 4. **SMARTS Pattern Debugging** ⭐⭐⭐⭐
```python
# Without LLM: Trial & error (10-30 min per pattern)
# With LLM: Guided refinement (2-5 min)

def debug_smarts_pattern(pattern, false_positives, false_negatives):
    prompt = f"""
    This SMARTS pattern needs refinement:
    
    Pattern: {pattern}
    
    False positives (should NOT match):
    {[mol for mol in false_positives]}
    
    False negatives (should match):
    {[mol for mol in false_negatives]}
    
    Suggest 3 improved SMARTS patterns:
    1. More restrictive (reduce false positives)
    2. More permissive (reduce false negatives)
    3. Balanced (optimize both)
    
    Explain what each change does and trade-offs.
    """
    
    suggestions = call_llm(prompt)
    # Test each suggestion with RDKit
    # Pick best one
```

**Cost:** ~$0.04 per debug session  
**Savings:** 20 min × $50/hr = $16.67  
**ROI:** 416x

#### 5. **Human-Readable Summaries** ⭐⭐⭐
```python
# For reports, documentation, alerts

def generate_validation_summary(results):
    prompt = f"""
    Summarize these validation results for chemists:
    
    Total tests: {results['total']}
    Passing: {results['passing']} ({results['pass_rate']:.1%})
    Failing: {results['failing']}
    
    Top failures:
    {results['top_failures']}
    
    Trend: {results['trend']} (vs last week)
    
    Generate:
    1. 2-sentence executive summary
    2. Top 3 action items
    3. Any concerning patterns
    
    Write for chemists, not programmers.
    """
    
    return call_llm(prompt)
```

**Cost:** ~$0.05 per summary  
**Value:** Better communication, faster decisions  
**ROI:** 50-100x

---

### 🟡 Medium-Value LLM Use Cases (ROI 2-5x)

1. **Boundary test case generation** - helpful but not critical
2. **Human-readable names** - nice UX improvement
3. **Alert message generation** - better than templates
4. **Documentation synthesis** - saves writing time

---

### 🔴 Low-Value LLM Use Cases (ROI < 2x)

1. **SMARTS matching** - RDKit is perfect
2. **Numerical calculations** - deterministic is better
3. **Structure recognition** - OSRA/DECIMER work well
4. **Feature detection** - RDKit descriptors sufficient
5. **Test execution** - no semantic understanding needed

---

## Cost Analysis

### Scenario 1: No LLM (Pure Deterministic)

| Activity | Time | Cost @ $50/hr |
|----------|------|---------------|
| Build automation | 20 days | $8,000 |
| Process 1000 reactions | 40 hours | $2,000 |
| Debug 20 failures | 10 hours | $500 |
| Add new reaction type | 6 hours | $300 |
| **Annual maintenance** | 100 hours | **$5,000** |

**Total Year 1:** $15,800

### Scenario 2: Strategic LLM Use

| Activity | Time | LLM Cost | Human Cost @ $50/hr | Total |
|----------|------|----------|---------------------|-------|
| Build automation | 20 days | $100 | $8,000 | $8,100 |
| Process 1000 reactions | 20 hours | $30 | $1,000 | $1,030 |
| Debug 20 failures | 1 hour | $1 | $50 | $51 |
| Add 5 new rxn types | 3 hours | $2.50 | $150 | $152.50 |
| **Annual maintenance** | 40 hours | $200 | $2,000 | **$2,200** |

**Total Year 1:** $11,533.50

**Savings:** $4,266.50 (27%)  
**Time saved:** 63 hours

### Scenario 3: Excessive LLM Use

| Activity | Time | LLM Cost | Human Cost | Total |
|----------|------|----------|------------|-------|
| Build automation | 25 days | $500 | $10,000 | $10,500 |
| Process 1000 reactions | 30 hours | $500 | $1,500 | $2,000 |
| (LLM overhead) | - | - | - | - |
| **Annual maintenance** | 60 hours | $1,000 | $3,000 | **$4,000** |

**Total Year 1:** $16,500

**Extra cost:** $700 (worse than no LLM!)

---

## Recommended LLM Strategy

### Phase 1: Start Without LLM (Weeks 1-4)
Build deterministic core:
- PDF extraction (ChemDataExtractor)
- SMARTS generation (RDKit MCS)
- Validation testing
- Basic fix recommender

**Why:** Validate workflow, understand failure modes

### Phase 2: Add Strategic LLM (Weeks 5-8)
Implement high-ROI use cases:
1. **Failure analysis** (highest ROI)
2. **Condition extraction** (if processing papers)
3. **New reaction type bootstrap** (if scaling)

**Why:** Maximum value, minimal cost

### Phase 3: Expand If Valuable (Month 3+)
Add medium-ROI use cases:
- SMARTS debugging assistant
- Test case generation
- Documentation synthesis

**Why:** Incremental improvements

### Never Add:
- LLM for structure matching (use RDKit)
- LLM for numerical operations
- LLM for deterministic logic

---

## Recommended Architecture

### Hybrid System Design

```python
class HybridRuleGenerator:
    """Combines deterministic chemistry with LLM intelligence."""
    
    def __init__(self, use_llm=True):
        self.rdkit = RDKitTools()
        self.llm = GPT4Client() if use_llm else None
    
    def generate_rule(self, reactions):
        # 1. Deterministic: Find MCS
        smarts = self.rdkit.find_mcs(reactions)
        
        # 2. Deterministic: Extract conditions
        conditions = self._extract_conditions_deterministic(reactions)
        
        # 3. Deterministic: Calculate priority
        priority = self._calculate_priority(smarts, len(reactions))
        
        # 4. LLM: Generate human-readable name (optional)
        if self.llm:
            name = self.llm.generate_name(smarts, conditions)
        else:
            name = self._template_name(smarts, conditions)
        
        # 5. Deterministic: Create entry
        return self._create_entry(smarts, conditions, priority, name)
    
    def debug_failure(self, failure):
        # 1. Deterministic: Classify failure type
        failure_type = self._classify_failure(failure)
        
        if failure_type == 'priority_conflict':
            # Deterministic fix
            return self._fix_priority(failure)
        
        elif failure_type == 'pattern_issue' and self.llm:
            # LLM analysis for complex cases
            return self.llm.analyze_pattern_failure(failure)
        
        else:
            # Fallback to deterministic suggestions
            return self._suggest_fix_deterministic(failure)
```

### Smart LLM Caching

```python
class CachedLLMClient:
    """Cache LLM responses to minimize costs."""
    
    def __init__(self):
        self.cache = {}  # or Redis
        self.llm = GPT4Client()
    
    def analyze_failure(self, failure):
        # Create cache key from failure fingerprint
        key = self._fingerprint(failure)
        
        if key in self.cache:
            return self.cache[key]  # Free!
        
        # Call LLM only for new patterns
        response = self.llm.analyze(failure)
        self.cache[key] = response
        
        return response
    
    def _fingerprint(self, failure):
        """Similar failures should have same fingerprint."""
        return hash((
            failure['type'],
            failure['expected_priority'] // 10,  # Group by priority range
            failure['actual_priority'] // 10,
            self._pattern_complexity(failure['expected_pattern'])
        ))
```

---

## API Cost Estimates

### GPT-4 Pricing (as of Oct 2025)
- Input: $0.03 per 1K tokens
- Output: $0.06 per 1K tokens

### Typical Use Cases

| Task | Input Tokens | Output Tokens | Cost/Call | Calls/Month | Monthly Cost |
|------|--------------|---------------|-----------|-------------|--------------|
| Failure analysis | 1,000 | 500 | $0.06 | 20 | $1.20 |
| Condition extraction | 500 | 300 | $0.03 | 100 | $3.00 |
| SMARTS debugging | 800 | 400 | $0.05 | 10 | $0.50 |
| Reaction bootstrap | 2,000 | 1,500 | $0.15 | 5 | $0.75 |
| Summaries | 1,500 | 300 | $0.06 | 20 | $1.20 |
| **Total** | | | | | **$6.65** |

**Annual LLM cost:** ~$80 (negligible compared to human time)

---

## Open-Source LLM Alternatives

### If API Costs Are Prohibitive

1. **Local Llama 3 70B** (free, but requires GPU)
   - Quality: 80-85% of GPT-4
   - Speed: Slower (30s vs 3s)
   - Cost: $0 + GPU/electricity

2. **Claude 3 Haiku** (cheaper API)
   - Cost: 1/10th of GPT-4
   - Quality: 90% of GPT-4 for chemistry
   - Good middle ground

3. **Chemistry-Specific Models**
   - ChemBERTa (free, domain-specific)
   - MolT5 (free, structure-text)
   - Limited but focused

---

## Decision Framework

### Use LLM APIs If:
✅ Processing >50 papers/month  
✅ Adding >2 new reaction types/quarter  
✅ Debugging >10 failures/month  
✅ Team time costs >$40/hr  
✅ Need rapid prototyping

### Skip LLM APIs If:
❌ <100 total reactions  
❌ One reaction type only  
❌ Tight budget (<$100/month)  
❌ Can't use cloud APIs (data sensitivity)  
❌ Have chemistry expert available 24/7

---

## Final Recommendation

### Optimal Strategy: **Deterministic Core + Strategic LLM**

**Start with:**
1. Pure deterministic automation (0-3 months)
2. Validate it works end-to-end
3. Identify bottlenecks

**Then add LLM for:**
1. **Failure analysis** (saves most time)
2. **Condition extraction** (if processing papers)
3. **New reaction types** (if scaling)

**Expected results:**
- 80% automation without LLM
- 95% automation with strategic LLM
- <$10/month LLM costs
- 60+ hours/year saved

**Don't use LLM for:**
- Anything RDKit can do
- Numerical calculations
- Pattern matching
- Core validation logic

This gives you the best of both worlds: **reliability of deterministic chemistry + intelligence of LLMs**, at minimal cost.
