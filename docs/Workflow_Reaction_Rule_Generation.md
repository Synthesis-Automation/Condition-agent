# Systematic Workflow for Reaction Rule Generation, Testing, and Improvement

## Overview

A comprehensive workflow to build reaction condition databases from literature, datasets, and web sources, with continuous validation and improvement cycles.

---

## Phase 1: Information Gathering & Extraction

### 1.1 Source Identification

**Input Sources:**

- Published articles (PDF)
- Reaction datasets (CSV, JSON, RDF, SDF)
- Web searches (SciFinder, Reaxys, patent databases)
- Experimental notebooks

**Tools & Methods:**

```python
# Example: Multi-source data collection
sources = {
    "pdf_articles": ["doi_10.1021/...", "path/to/paper.pdf"],
    "datasets": ["USPTO_grants.csv", "CAS_reactions.json"],
    "web_searches": ["SciFinder query: 'Suzuki coupling heteroaryl'"],
    "internal_data": ["lab_notebook_2024.xlsx"]
}
```

### 1.2 Content Extraction

**For PDFs:**

1. **Text extraction** (PyPDF2, pdfplumber)
   - Extract experimental sections
   - Parse reaction schemes
   - Identify conditions tables

2. **Chemical structure recognition** (OSRA, ChemDataExtractor)
   - Convert images to SMILES
   - Extract reaction arrows and mapping

3. **Table parsing**
   - Identify condition tables
   - Extract substrate scope data

**For Datasets:**

1. **Format standardization**
   - Convert all to common format (JSONL recommended)
   - Map field names to standard schema
   - Validate SMILES/reaction strings

2. **Quality filtering**
   - Remove duplicates
   - Filter by yield thresholds
   - Validate chemical structures

**Tool Recommendations:**

```bash
# PDF processing
- ChemDataExtractor (chemistry-aware NLP)
- RXNMapper (reaction mapping)
- DECIMER (optical structure recognition)

# Dataset processing
- RDKit (structure validation)
- pandas (data manipulation)
- reaction-miner (pattern extraction)
```

---

## Phase 2: Rule Generation

### 2.1 SMARTS Pattern Creation (Primary Method)

**Approach:** Use SMARTS for substrate-specific matching

**Workflow:**

1. **Identify substrate classes**

   ```
   Examples:
   - Aryl halides: [c:1]-[Br,I,Cl:2]
   - Heteroaryl: [c:1]1[n,o,s][c,n]...[Br:2]
   - Vinyl halides: [C:1]=[C:2]-[Br:3]
   - Bulky substrates: [c:1]([!H])([!H])
   ```

2. **Create hierarchical patterns (general → specific)**

   ```
   Priority 0:  DEFAULT (any ArBr) → Pd(PPh3)4
   Priority 40: GENERAL-ArBr → Pd2(dba)3/SPhos
   Priority 50: ORTHO-ArBr → Pd2(dba)3/XPhos
   Priority 60: HETERO-ArBr → Pd(OAc)2/SPhos
   ```

3. **Test pattern specificity**

   ```python
   from rdkit import Chem
   
   pattern = Chem.MolFromSmarts("[c:1]-[Br:2].[B:3](-[O])(-[O])-[c:4]([!H])([!H])")
   
   # Test against known substrates
   test_cases = {
       "simple": "Brc1ccccc1.OB(O)c1ccccc1",  # Should NOT match
       "bulky": "Brc1ccccc1.OB(O)c1c(C)cccc1C"  # Should match
   }
   ```

### 2.2 Feature Requirements (Secondary/Complementary)

**Use for properties hard to express in SMARTS:**

```json
{
  "feature_requirements": {
    "electrophile.electronics": ["electron-poor"],
    "electrophile.ring_hetero_count": [1, 2, 3],
    "nucleophile.ortho_sub_count": [2]
  }
}
```

**When to use:**

- Electronic properties (electron-rich/poor)
- Computed properties (logP, HOMO-LUMO)
- Complex geometric features

### 2.3 Priority Assignment Strategy

```
0-10:    DEFAULT rules (catch-all, lowest priority)
20-40:   GENERAL methods (broad applicability)
45-60:   SPECIFIC substrates (particular functional groups)
65-80:   VERY SPECIFIC (multiple criteria, challenging cases)
85-100:  SCHEME-BASED (exact literature matches)
```

**Rule of thumb:**

- Higher specificity = Higher priority
- Literature-validated > General knowledge
- Scheme matches > Database matches

---

## Phase 3: Database Entry Creation

### 3.1 Entry Template

```json
{
  "id": "SCDB-{RXNTYPE}-{SUBSTRATE}-{LIGAND}",
  "reaction_type": "Suzuki_Miyaura",
  "name": "Human-readable description",
  "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
  "token_signature": ["Suzuki", "ArBr", "SPhos", "K2CO3"],
  "conditions": {
    "pd_source": ["Pd2(dba)3"],
    "ligands": ["SPhos"],
    "boron_partner": ["aryl-B(OH)2 1.2-1.3 eq"],
    "base": ["K2CO3 2.0 eq"],
    "solvent": ["THF/H2O (4:1)"],
    "temperature_C": [45, 60],
    "time_h": [8, 12],
    "loadings": {
      "Pd_mol%": [0.5, 1.5],
      "ligand_mol%": [1.0, 3.0],
      "base_eq": [2.0, 2.5]
    }
  },
  "env": {
    "features_from_smiles": {
      "electrophile.lg_class": "parse from SMARTS"
    },
    "feature_requirements": {}
  },
  "evidence": {
    "provenance": "DOI:10.1021/... or dataset:USPTO_2024",
    "last_updated": "2025-10-05T00:00:00Z",
    "confidence": 0.95
  },
  "notes": [
    "Validated on 20+ substrates",
    "Preferred for electron-rich aryl halides"
  ],
  "priority": 50
}
```

### 3.2 Automated Entry Generation Script

```python
def generate_rule_from_literature(paper_data):
    """
    Extract and format rule from parsed paper data.
    """
    rule = {
        "id": f"SCDB-{paper_data['reaction_type']}-{generate_substrate_key(paper_data)}",
        "reaction_type": paper_data['reaction_type'],
        "name": paper_data['title'],
        "rxn_smiles_min": generate_smarts_from_examples(paper_data['examples']),
        "conditions": extract_conditions(paper_data),
        "priority": assign_priority_from_scope(paper_data),
        "evidence": {
            "provenance": paper_data['doi'],
            "last_updated": datetime.now().isoformat(),
            "confidence": calculate_confidence(paper_data)
        }
    }
    return rule

def generate_smarts_from_examples(examples):
    """
    Analyze multiple reaction examples to generate generalized SMARTS.
    """
    # Find common substructures
    common_pattern = find_mcs(examples)
    # Generate SMARTS with atom mapping
    return create_reaction_smarts(common_pattern)
```

---

## Phase 4: Test Case Generation

### 4.1 Create Representative Test Cases

**For each rule, generate 1-3 test reactions:**

```python
# tests/sample_reactions.py
REACTION_DB_TEST_REACTIONS = {
    "SCDB-SUZ-BULKY-NUC-XPHOS": {
        "smiles": "Brc1ccccc1.OB(O)c1c(C)cccc1C>>Cc1cccc(C)c1-c1ccccc1",
        "description": "PhBr + 2,6-dimethylphenylboronic acid (hindered nucleophile)",
        "expected_features": {
            "lg_class": "Br",
            "nucleophile_ortho_sub": 2
        },
        "literature_yield": 0.85,
        "literature_ref": "DOI:10.1021/..."
    }
}
```

**Test case selection criteria:**

1. **Positive match** - Should match this specific rule
2. **Boundary case** - Near decision boundary with other rules
3. **Negative case** - Should NOT match (for specificity testing)

### 4.2 Validation Script Template

```python
def validate_single_rule(rule_id, test_data, db):
    """
    Test if a rule can be reached with its test case.
    """
    result = match(db, rxn_smiles=test_data['smiles'])
    
    # Check if correct rule matched
    if result.entry_id == rule_id:
        return {"status": "PASS", "result": result}
    else:
        return {
            "status": "FAIL",
            "expected": rule_id,
            "got": result.entry_id,
            "priority_expected": get_priority(db, rule_id),
            "priority_got": result.priority,
            "reason": diagnose_mismatch(rule_id, result, test_data)
        }
```

---

## Phase 5: Testing & Validation

### 5.1 Automated Testing Pipeline

```python
# scripts/validate_all_rules.py
def run_validation_suite(database_path, test_cases_path):
    """
    Run comprehensive validation and generate reports.
    """
    db = load_db(database_path)
    tests = load_test_cases(test_cases_path)
    
    results = []
    for rule_id, test_data in tests.items():
        result = validate_single_rule(rule_id, test_data, db)
        results.append(result)
    
    # Generate reports
    generate_markdown_report(results, "validation_report.md")
    generate_conditions_report(results, "conditions_report.md")
    generate_json_report(results, "validation_results.json")
    
    # Calculate metrics
    pass_rate = sum(1 for r in results if r['status'] == 'PASS') / len(results)
    
    return {
        "pass_rate": pass_rate,
        "results": results,
        "recommendations": generate_fix_recommendations(results)
    }
```

### 5.2 Metrics to Track

```python
metrics = {
    "coverage": {
        "total_rules": 23,
        "tested_rules": 23,
        "passing_rules": 10,
        "pass_rate": 0.435
    },
    "specificity": {
        "true_positives": 10,   # Correct matches
        "false_positives": 0,    # Wrong matches
        "false_negatives": 13    # Missed matches
    },
    "priority_conflicts": {
        "expected_behavior": 5,  # Lower priority beaten by higher
        "unexpected": 8          # Need investigation
    }
}
```

### 5.3 Quality Checks

1. **SMARTS syntax validation**

   ```python
   from rdkit import Chem
   rxn = Chem.ReactionFromSmarts(smarts_pattern)
   if not rxn:
       report_error("Invalid SMARTS pattern")
   ```

2. **Priority consistency check**

   ```python
   # Ensure more specific rules have higher priority
   if is_subset(rule_A_pattern, rule_B_pattern):
       assert rule_A.priority > rule_B.priority
   ```

3. **Condition completeness**

   ```python
   required_fields = ["pd_source", "base", "solvent", "temperature_C"]
   for field in required_fields:
       assert field in entry['conditions']
   ```

---

## Phase 6: Iterative Improvement

### 6.1 Analyze Failures

**Categorize failure types:**

```python
def diagnose_failure(failure_result):
    """
    Determine why a rule didn't match.
    """
    categories = {
        "smarts_too_broad": check_if_matches_too_many(failure_result),
        "smarts_too_specific": check_if_matches_none(failure_result),
        "priority_too_low": check_priority_conflict(failure_result),
        "feature_not_detected": check_feature_requirements(failure_result),
        "pattern_syntax_error": check_smarts_validity(failure_result)
    }
    
    return [cat for cat, is_issue in categories.items() if is_issue]
```

### 6.2 Fix Strategies

**1. SMARTS Too Broad → Make More Specific**

```python
# Before: Matches everything
"[c:1]-[Br:2]"

# After: Requires heteroatom in ring
"[c:1]1[n,o,s][c,n][c,n][c,n][c,n]1-[Br:2]"
```

**2. SMARTS Too Specific → Relax Constraints**

```python
# Before: Too rigid, doesn't match variations
"[c:1]1[cH][cH][cH][cH][cH]1-[Br:2]"

# After: Allow substitution
"[c:1]1[c,cH][c,cH][c,cH][c,cH][c,cH]1-[Br:2]"
```

**3. Priority Conflict → Adjust Hierarchy**

```python
# If rule A (priority 50) should beat rule B (priority 60):
rule_A.priority = 65  # Increase above B
```

**4. Feature Detection → Use SMARTS Instead**

```python
# Before: Rely on feature detection
"feature_requirements": {"nucleophile.ortho_sub_count": [2]}

# After: Encode in SMARTS
"rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])(-[O])-[c:4]([!H])([!H])>>[c:1]-[c:4]"
```

### 6.3 Automated Fix Suggestions

```python
def suggest_fixes(failure_result):
    """
    Automatically suggest fixes based on failure type.
    """
    suggestions = []
    
    if failure_result['reason'] == 'priority_too_low':
        conflicting_rule = failure_result['got']
        suggestions.append({
            "action": "increase_priority",
            "rule": failure_result['expected'],
            "new_priority": conflicting_rule.priority + 5
        })
    
    if failure_result['reason'] == 'smarts_too_broad':
        suggestions.append({
            "action": "add_specificity",
            "pattern": failure_result['pattern'],
            "recommendation": "Add constraints for specific substituents"
        })
    
    return suggestions
```

---

## Phase 7: Continuous Integration

### 7.1 CI/CD Pipeline

```yaml
# .github/workflows/validate_rules.yml
name: Reaction Rules Validation

on: [push, pull_request]

jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.10'
      - name: Install dependencies
        run: |
          pip install -r requirements.txt
      - name: Run validation suite
        run: |
          python scripts/validate_all_rules.py
      - name: Check pass rate threshold
        run: |
          python scripts/check_metrics.py --min-pass-rate 0.80
      - name: Upload reports
        uses: actions/upload-artifact@v2
        with:
          name: validation-reports
          path: docs/*_report.md
```

### 7.2 Version Control Strategy

```bash
# Tag each validated version
git tag -a v1.0-suzuki -m "Suzuki rules: 43.5% pass rate, 23 entries"
git push origin v1.0-suzuki

# Track metrics over time
metrics/
├── v1.0_metrics.json
├── v1.1_metrics.json
└── metrics_history.csv
```

---

## Phase 8: Scaling to New Reaction Types

### 8.1 Reaction Type Template

```python
# chemtools/reaction_types/{reaction_name}.py

class BuchwaldHartwigReactionType:
    """
    Template for C-N coupling reactions.
    """
    
    @staticmethod
    def get_base_smarts():
        return "[c:1]-[X:2].[N:3]>>[c:1]-[N:3]"
    
    @staticmethod
    def get_substrate_classes():
        return {
            "electrophile": ["ArBr", "ArCl", "ArI", "ArOTf"],
            "nucleophile": ["primary_amine", "secondary_amine", "aniline"]
        }
    
    @staticmethod
    def get_common_conditions():
        return {
            "pd_sources": ["Pd2(dba)3", "Pd(OAc)2"],
            "ligands": ["BINAP", "XPhos", "tBuXPhos", "BrettPhos"],
            "bases": ["NaOtBu", "Cs2CO3", "K3PO4"],
            "solvents": ["toluene", "dioxane", "THF"]
        }
```

### 8.2 Multi-Reaction Workflow

```python
reaction_types = [
    "Suzuki_Miyaura",
    "Buchwald_Hartwig",
    "Negishi",
    "Heck",
    "Sonogashira",
    "Ullmann_CN"
]

for rxn_type in reaction_types:
    # 1. Gather literature
    papers = search_literature(rxn_type, year_range="2020-2025")
    
    # 2. Extract rules
    rules = extract_rules_from_papers(papers)
    
    # 3. Generate test cases
    tests = generate_test_cases(rules)
    
    # 4. Validate
    results = validate_rules(rules, tests)
    
    # 5. Improve
    improved_rules = improve_based_on_results(rules, results)
    
    # 6. Save
    save_to_database(f"{rxn_type.lower()}_db.json", improved_rules)
```

---

## Phase 9: Advanced Techniques

### 9.1 Machine Learning Integration

```python
from sklearn.ensemble import RandomForestClassifier

def ml_assisted_rule_generation(reaction_data):
    """
    Use ML to identify important features and suggest SMARTS patterns.
    """
    # Extract features from successful reactions
    features = extract_molecular_features(reaction_data)
    X = features[['descriptors']]
    y = features[['success']]
    
    # Train model
    model = RandomForestClassifier()
    model.fit(X, y)
    
    # Get feature importance
    important_features = get_top_features(model, n=10)
    
    # Suggest SMARTS based on important features
    suggested_patterns = generate_smarts_from_features(important_features)
    
    return suggested_patterns
```

### 9.2 Crowdsourced Validation

```python
def create_validation_task(rule, test_case):
    """
    Generate human-reviewable validation tasks.
    """
    return {
        "rule_id": rule['id'],
        "substrate": render_structure(test_case['smiles']),
        "predicted_conditions": rule['conditions'],
        "question": "Would these conditions work for this substrate?",
        "options": ["Yes", "No", "Needs modification"],
        "literature_context": rule['evidence']['provenance']
    }
```

### 9.3 Reaction Outcome Prediction

```python
def predict_reaction_outcome(substrate, conditions, model):
    """
    Predict yield and selectivity for given conditions.
    """
    features = encode_reaction(substrate, conditions)
    prediction = model.predict(features)
    
    return {
        "predicted_yield": prediction['yield'],
        "confidence": prediction['confidence'],
        "alternative_conditions": suggest_alternatives(substrate, model)
    }
```

---

## Complete Workflow Summary

```
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 1: Information Gathering                                      │
│  └─ PDF extraction, Dataset parsing, Web scraping                   │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 2: Rule Generation                                            │
│  └─ SMARTS patterns, Feature requirements, Priority assignment      │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 3: Database Entry Creation                                    │
│  └─ JSON schema, Condition extraction, Metadata annotation          │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 4: Test Case Generation                                       │
│  └─ Positive/negative cases, Boundary testing                       │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 5: Testing & Validation                                       │
│  └─ Automated testing, Metrics tracking, Quality checks             │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 6: Iterative Improvement                                      │
│  └─ Failure analysis, Pattern refinement, Priority adjustment       │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 7: Continuous Integration                                     │
│  └─ CI/CD pipeline, Version control, Regression testing             │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 8: Scale to New Reaction Types                                │
│  └─ Template reuse, Multi-reaction databases                        │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 9: Advanced Techniques                                        │
│  └─ ML integration, Crowdsourcing, Outcome prediction               │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Key Success Factors

### ✅ Best Practices

1. **SMARTS-first approach** - More reliable than complex feature detection
2. **Hierarchical priorities** - Clear precedence rules
3. **Comprehensive testing** - Every rule has test cases
4. **Iterative refinement** - Continuous improvement based on validation
5. **Version control** - Track changes and metrics over time
6. **Documentation** - Clear provenance and notes for each rule

### ⚠️ Common Pitfalls to Avoid

1. Over-reliance on feature detection (hard to maintain)
2. Vague SMARTS patterns that match too broadly
3. No test cases (can't validate if rules are reachable)
4. Ignoring priority conflicts
5. Not tracking metrics over time
6. Missing literature provenance

### 📊 Success Metrics

- **Coverage:** % of known substrates covered by rules
- **Pass Rate:** % of test cases correctly matched
- **Precision:** % of matches that are correct
- **Recall:** % of expected matches that are found
- **F1 Score:** Harmonic mean of precision and recall

---

## Tools & Resources

### Recommended Stack

- **Chemistry:** RDKit, OpenBabel
- **NLP:** ChemDataExtractor, Grobid
- **Data:** pandas, numpy
- **Testing:** pytest, hypothesis
- **CI/CD:** GitHub Actions, GitLab CI
- **Visualization:** RDKit, matplotlib, seaborn

### Community Resources

- RDKit mailing list
- CompChem forums
- ACS Division of Computers in Chemistry
- GitHub reaction databases (USPTO, ORD)

---

*End of Workflow Document*
