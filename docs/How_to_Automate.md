# How to Automate the Reaction Rule Workflow

## Executive Summary

This document provides a practical guide for automating the reaction rule generation workflow based on our experience with Suzuki coupling rules.

## Three Levels of Automation

### Level 1: Semi-Automated (Ready Now) ✅

**What's Automated:**
- Validation testing
- Report generation  
- Status monitoring
- Basic fix recommendations

**Usage:**
```bash
# Check system status
python -m automation.orchestrator status

# Run validation
python -m automation.orchestrator validate

# Analyze failures (basic)
python -m automation.orchestrator improve --auto-fix
```

**Best For:**
- Current workflow
- Getting started with automation
- Manual review of changes

### Level 2: Batch Automated (Requires Implementation) 🟡

**What Would Be Automated:**
- Dataset batch processing
- Rule generation from clusters
- Test case generation
- Continuous validation
- Automated git commits

**Implementation Needed:**
- Dataset clustering (`automation/ingestion/dataset_processor.py`)
- MCS-based SMARTS extraction (`automation/rule_generator.py`)
- Continuous validator service (`automation/validation/validator_service.py`)

**Best For:**
- Processing large reaction datasets
- Periodic batch updates
- Scheduled improvements

### Level 3: Fully Automated (Advanced) ⚪

**What Would Be Automated:**
- PDF monitoring and extraction
- Real-time rule updates
- ML-assisted pattern generation
- Crowdsourced validation
- Production deployment

**Implementation Needed:**
- PDF extraction pipeline (ChemDataExtractor, OSRA)
- Airflow/Prefect orchestration
- Monitoring dashboard
- ML models for rule suggestion

**Best For:**
- Production systems
- Large organizations
- Continuous learning

---

## Quick Start: Current Automation (Level 1)

### 1. Install Dependencies

```bash
# Core dependencies (already installed)
pip install -r requirements.txt

# Optional automation extras
pip install -r requirements-automation.txt
```

### 2. Check System Status

```bash
python -m automation.orchestrator status
```

Output shows:
- Number of database entries
- Priority distribution
- Test case count
- Validation pass rate
- Available automation components

### 3. Run Validation

```bash
python -m automation.orchestrator validate
```

This runs the existing `validate_all_suzuki_rules.py` script and displays results.

### 4. Full Pipeline (Current State)

```bash
python -m automation.orchestrator run-full-pipeline
```

Executes all implemented steps in sequence.

---

## Incremental Implementation Plan

### Phase 1: Enhance Fix Recommender (1-2 days)

**Goal:** Automatically diagnose and fix common issues

**Implementation:**
1. Enhance `automation/validation/fix_recommender.py`:
   - Add SMARTS pattern analysis
   - Detect patterns that are too broad/specific
   - Suggest pattern refinements
   - Auto-fix priority conflicts

2. Integration:
   ```bash
   python -m automation.orchestrator improve --auto-fix
   ```

**Value:** Reduces manual debugging time by 50%+

### Phase 2: Dataset Batch Processing (3-5 days)

**Goal:** Process reaction datasets automatically

**Implementation:**
1. Implement `automation/ingestion/dataset_processor.py`:
   - Load and validate datasets
   - Filter by reaction type and yield
   - Cluster similar reactions using DRFP

2. Implement `automation/rule_generator.py`:
   - Extract MCS from reaction clusters
   - Generate SMARTS patterns
   - Calculate priorities automatically
   - Create database entries

3. Usage:
   ```bash
   python -m automation.orchestrator process-dataset \
       --dataset data/uspto_grants.csv \
       --reaction-type Suzuki \
       --min-yield 0.70 \
       --output automation/output/suzuki_rules
   ```

**Value:** Scale to thousands of reactions automatically

### Phase 3: Continuous Validation (2-3 days)

**Goal:** Monitor database health continuously

**Implementation:**
1. Implement `automation/validation/validator_service.py`:
   - Schedule periodic validation
   - Track metrics over time
   - Generate trend reports
   - Alert on regressions

2. Usage:
   ```bash
   # Run validation every 6 hours
   python -m automation.orchestrator continuous-validation --interval-hours 6
   ```

**Value:** Catch issues early, track improvement over time

### Phase 4: PDF Extraction (1-2 weeks)

**Goal:** Extract rules from literature PDFs

**Implementation:**
1. Set up dependencies:
   - ChemDataExtractor (chemistry NLP)
   - OSRA (optical structure recognition)
   - Grobid (PDF parsing)

2. Implement `automation/ingestion/pdf_monitor.py`:
   - Watch directory for new PDFs
   - Extract reaction schemes
   - Parse experimental conditions
   - Generate rules automatically

3. Usage:
   ```bash
   # Start monitoring
   python -m automation.orchestrator watch-pdfs --pdf-dir ./incoming_papers
   ```

**Value:** Extract rules from literature automatically

### Phase 5: Orchestration & Monitoring (1 week)

**Goal:** Production-ready automation

**Implementation:**
1. Set up Airflow/Prefect DAGs
2. Create monitoring dashboard (Streamlit)
3. Add alerting (Slack, email)
4. Containerize with Docker
5. CI/CD pipeline

**Value:** Reliable, monitored production system

---

## Architecture Decision Guide

### When to Use What

| Scenario | Recommended Approach |
|----------|---------------------|
| Small dataset (<100 reactions) | Manual + Level 1 validation |
| Medium dataset (100-1000) | Level 2 batch processing |
| Large dataset (1000+) | Level 2 + parallel processing |
| Periodic updates | Level 3 continuous validation |
| Literature monitoring | Level 3 PDF extraction |
| Production deployment | Level 3 full automation |

### Technology Choices

**Workflow Orchestration:**
- **Simple:** CLI (`automation/orchestrator.py`)
- **Medium:** Schedule library (Python cron)
- **Advanced:** Airflow/Prefect

**Data Processing:**
- **Simple:** Pandas + RDKit
- **Medium:** + scikit-learn clustering
- **Advanced:** + Spark for large datasets

**Monitoring:**
- **Simple:** Terminal output + markdown reports
- **Medium:** Streamlit dashboard
- **Advanced:** Grafana + Prometheus

---

## Implementation Checklist

### Ready to Use Now ✅
- [x] CLI orchestrator
- [x] Validation runner
- [x] Status monitoring
- [x] Report generation
- [x] Basic fix recommender

### Next Steps (Prioritized)
- [ ] Enhance fix recommender with SMARTS analysis
- [ ] Implement continuous validation service
- [ ] Add dataset clustering
- [ ] Implement rule generation from clusters
- [ ] Add automated test case generation
- [ ] Set up monitoring dashboard
- [ ] Implement PDF extraction pipeline
- [ ] Add Airflow orchestration
- [ ] Containerize with Docker
- [ ] CI/CD integration

---

## Example Workflows

### Workflow 1: Weekly Dataset Update

```bash
# 1. Process new dataset
python -m automation.orchestrator process-dataset \
    --dataset data/new_reactions.csv \
    --reaction-type Suzuki \
    --min-yield 0.70

# 2. Validate new rules
python -m automation.orchestrator validate

# 3. Fix issues automatically
python -m automation.orchestrator improve --auto-fix

# 4. Re-validate
python -m automation.orchestrator validate

# 5. Generate reports
python scripts/generate_conditions_report.py > docs/conditions_report.md

# 6. Commit if improved
git add data/conditionDB/*.json
git commit -m "Update: Added rules from new dataset"
```

### Workflow 2: Continuous Monitoring

```bash
# Start continuous validation (runs in background)
python -m automation.orchestrator continuous-validation --interval-hours 6 &

# Monitor dashboard
streamlit run automation/monitoring/dashboard.py
```

### Workflow 3: Literature Update

```bash
# 1. Drop PDF in incoming_papers/
# 2. Monitor extracts rules automatically
python -m automation.orchestrator watch-pdfs --pdf-dir ./incoming_papers &

# 3. Validation runs automatically
# 4. Alerts sent if issues detected
```

---

## Cost-Benefit Analysis

### Time Investment

| Phase | Implementation Time | Maintenance | Annual Savings |
|-------|-------------------|-------------|----------------|
| Level 1 (Current) | 0 days (done) | 1 hr/month | 10 hrs/year |
| Fix Recommender | 2 days | 0.5 hr/month | 40 hrs/year |
| Batch Processing | 5 days | 1 hr/month | 100 hrs/year |
| Continuous Validation | 3 days | 0.5 hr/month | 20 hrs/year |
| PDF Extraction | 10 days | 2 hrs/month | 150 hrs/year |
| Full Automation | 20 days | 4 hrs/month | 300 hrs/year |

### When to Automate

**Automate if:**
- Processing >100 reactions/month
- Multiple reaction types to manage
- Team size >3 people
- Need reproducibility
- Frequent updates from literature

**Stay Manual if:**
- <50 reactions total
- One-time project
- Rapid prototyping phase
- No development resources

---

## Best Practices

### Start Small
1. Use Level 1 automation first
2. Validate manually until confident
3. Automate incrementally
4. Keep manual override capability

### Version Control
```bash
# Tag each validated version
git tag -a v1.0-suzuki -m "Suzuki rules: 43.5% pass rate"

# Track metrics
echo "$(date),suzuki,0.435,23" >> metrics/history.csv
```

### Testing
- Test automation on small datasets first
- Validate automation output manually
- Keep regression tests
- Monitor for drift

### Documentation
- Document all automation decisions
- Keep examples of manual vs automated output
- Track what was automated and why
- Record failure modes

---

## Troubleshooting

### Common Issues

**Issue:** Automation produces invalid SMARTS
- **Solution:** Add validation step, fallback to manual review
- **Prevention:** Test on known examples first

**Issue:** Too many false positives
- **Solution:** Increase specificity thresholds
- **Prevention:** Start with high yield cutoffs

**Issue:** Rules conflict after automation
- **Solution:** Run priority conflict detector
- **Prevention:** Use hierarchical priority assignment

### Debugging

```bash
# Run with verbose output
python -m automation.orchestrator validate --verbose

# Check logs
tail -f automation/logs/validation.log

# Inspect intermediate outputs
ls automation/output/
```

---

## Resources

### Documentation
- `docs/Workflow_Reaction_Rule_Generation.md` - Complete manual workflow
- `docs/Automation_Framework.md` - Full automation implementation
- `automation/README.md` - Component status and roadmap

### Code
- `automation/orchestrator.py` - Main CLI entrypoint
- `scripts/validate_all_suzuki_rules.py` - Validation script
- `automation/validation/fix_recommender.py` - Fix suggestions

### Data
- `data/conditionDB/` - Reaction databases
- `tests/sample_reactions.py` - Test cases
- `docs/*_report.md` - Validation reports

---

## Next Actions

1. **Immediate (This Week):**
   - Run `python -m automation.orchestrator status`
   - Review current validation results
   - Identify top 3 failing rules to fix manually

2. **Short-term (This Month):**
   - Enhance fix recommender
   - Implement continuous validation
   - Process one new dataset with batch automation

3. **Long-term (This Quarter):**
   - Full batch processing pipeline
   - Monitoring dashboard
   - Consider PDF extraction if needed

---

**Remember:** Automation is a journey, not a destination. Start with what provides the most value for the least effort, and build incrementally based on real needs.
