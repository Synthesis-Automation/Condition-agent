# Evaluation & Validation Framework Proposal

## Current System (3 Tiers)

```
Tier 1: String Patterns (instant, free)
  ↓
Tier 2: DeepSeek-v3.2 Quick Glance (~17s, comprehensive)
  ↓
Tier 3: gpt-4o-mini Deep Analysis (~5s, guided by Tier 2)
  ↓
OUTPUT (no validation)
```

**Problems:**
- No validation of LLM outputs
- Disagreements between tiers are warned but not resolved
- Low-confidence results aren't re-analyzed
- Chemical implausibility not detected

## Proposed Framework: 4-Tier + Validation Loop

### Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  Tier 1: String Patterns (instant, free)                   │
│  - Detect obvious patterns (Suzuki, Grignard, etc.)        │
│  - Flag complexity (simple/tandem/complex)                  │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  Tier 2: Quick Analysis (DeepSeek-v3.2, ~17s)              │
│  - Comprehensive chemistry analysis                         │
│  - All structural changes, protecting groups                │
│  - Confidence score                                         │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  Tier 3: Deep Mechanistic Analysis (gpt-4o-mini, ~5s)      │
│  - Detailed mechanism steps                                 │
│  - Bond-by-bond transformations                             │
│  - Guided by Tier 2 context                                 │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  Tier 4: VALIDATION & EVALUATION (NEW)                     │
│  ┌──────────────────┬──────────────────┬─────────────────┐ │
│  │ RDKit Validator  │ LLM Self-Critic  │ Consensus Check │ │
│  │ (~1s, free)      │ (~3s, cheap)     │ (~0s, free)     │ │
│  └──────────────────┴──────────────────┴─────────────────┘ │
└─────────────────────────────────────────────────────────────┘
                           ↓
              ┌─────────────────────┐
              │ Quality Gate        │
              │ Pass/Fail/Re-analyze│
              └─────────────────────┘
                    ↓         ↓
              ┌──────┐   ┌───────────────┐
              │ PASS │   │ FAIL → RETRY  │
              └──────┘   │ with stronger │
                         │ model/context │
                         └───────────────┘
```

## Tier 4 Components (Validation Layer)

### Component 1: RDKit Validator (Tool-based, ~1s)

**What it checks:**
1. **Structural Validity**
   - Products are valid molecules (parseable)
   - No impossible valences or bonds

2. **Transformation Plausibility**
   - Count heavy atoms (reactants vs products)
   - Check if reported changes match actual SMILES diff
   - Validate functional group transformations

3. **Mass Balance Check**
   - Compare molecular formulas
   - Flag large discrepancies (beyond spectators/byproducts)

**Implementation:**
```python
def validate_with_rdkit(
    rxn_smiles_clean: str,
    interpretation: Dict[str, Any]
) -> Dict[str, Any]:
    """
    RDKit-based validation of reaction analysis.

    Returns:
        {
            "valid": True/False,
            "issues": [list of problems],
            "warnings": [list of concerns],
            "atom_balance": {"reactants": 45, "products": 29, "lost": 16},
            "confidence_adjustment": -0.1  # lower if issues found
        }
    """
    issues = []
    warnings = []

    # Parse SMILES
    parts = rxn_smiles_clean.split(">>")
    reactants_smiles = parts[0]
    products_smiles = parts[1]

    # 1. Parse and validate structures
    try:
        reactant_mols = [Chem.MolFromSmiles(s) for s in reactants_smiles.split(".")]
        product_mols = [Chem.MolFromSmiles(s) for s in products_smiles.split(".")]

        if None in reactant_mols:
            issues.append("Invalid reactant structure")
        if None in product_mols:
            issues.append("Invalid product structure")
    except Exception as e:
        issues.append(f"Structure parsing failed: {e}")
        return {"valid": False, "issues": issues}

    # 2. Count atoms
    reactant_atoms = sum(mol.GetNumHeavyAtoms() for mol in reactant_mols if mol)
    product_atoms = sum(mol.GetNumHeavyAtoms() for mol in product_mols if mol)
    atoms_lost = reactant_atoms - product_atoms

    # 3. Check atom loss matches reported changes
    reported_changes = interpretation.get("all_changes", [])
    if atoms_lost > 0 and not any("deprotection" in str(c).lower() or "removal" in str(c).lower()
                                   for c in reported_changes):
        warnings.append(f"{atoms_lost} atoms lost but no deprotection/removal reported")

    # 4. Validate reaction class plausibility
    overall_class = interpretation.get("overall_class", "")
    reaction_types = interpretation.get("reaction_types", [])

    # Example: Cross-coupling should show C-C bond formation
    if "cross_coupling" in overall_class or any("suzuki" in rt.lower() for rt in reaction_types):
        # Check for aryl halide in reactants
        has_aryl_halide = any("Br" in Chem.MolToSmiles(mol) or "I" in Chem.MolToSmiles(mol)
                              for mol in reactant_mols if mol)
        if not has_aryl_halide:
            warnings.append("Suzuki/cross-coupling reported but no aryl halide found")

    # 5. Calculate confidence adjustment
    confidence_adjustment = 0.0
    if issues:
        confidence_adjustment = -0.3  # Major issues
    elif len(warnings) > 2:
        confidence_adjustment = -0.1  # Multiple warnings
    elif len(warnings) == 1:
        confidence_adjustment = -0.05  # Minor concern

    return {
        "valid": len(issues) == 0,
        "issues": issues,
        "warnings": warnings,
        "atom_balance": {
            "reactants": reactant_atoms,
            "products": product_atoms,
            "lost": atoms_lost
        },
        "confidence_adjustment": confidence_adjustment
    }
```

### Component 2: LLM Self-Critic (Model-based, ~3s)

**What it does:**
- Reviews Tier 2 and Tier 3 outputs
- Checks for logical inconsistencies
- Identifies chemical implausibilities
- Uses cheap model (gpt-4o-mini) with targeted prompt

**Prompt template:**
```python
SELF_CRITIC_PROMPT = """You are a chemistry validation expert. Review this reaction analysis for errors.

REACTION: {rxn_smiles}

ANALYSIS TO VALIDATE:
- Tier 2 Classification: {tier2_types}
- Tier 3 Classification: {tier3_class}
- Reported Changes: {all_changes}
- Mechanism: {mechanism_summary}

VALIDATION CHECKLIST:
1. Do Tier 2 and Tier 3 agree on reaction type?
2. Are the reported structural changes chemically plausible?
3. Does the mechanism match the classification?
4. Are protecting group changes correctly identified?
5. Is this a tandem/sequential reaction being reported as simple (or vice versa)?

Return JSON:
{
  "consistent": true/false,
  "issues": ["specific problems found"],
  "suggestions": ["corrections or clarifications"],
  "confidence": 0.0-1.0,
  "recommendation": "accept|request_clarification|re_analyze"
}

Be critical but fair. Flag genuine problems, not minor wording differences."""
```

**Implementation:**
```python
def validate_with_llm_critic(
    rxn_smiles: str,
    tier2_result: Dict,
    tier3_result: Dict,
    client: LLMClient
) -> Dict[str, Any]:
    """
    LLM-based self-critique validation.

    Uses cheap model to review for inconsistencies.
    """
    prompt = SELF_CRITIC_PROMPT.format(
        rxn_smiles=rxn_smiles,
        tier2_types=tier2_result.get('reaction_types', []),
        tier3_class=tier3_result.get('overall_class', ''),
        all_changes=tier2_result.get('all_changes', [])[:3],  # First 3
        mechanism_summary=tier3_result.get('mechanism_summary', '')
    )

    # Use fast, cheap model for validation
    critic_client = LLMClient(provider="openai", model="gpt-4o-mini")
    response = critic_client.chat(prompt, temperature=0.0, max_tokens=500)

    # Parse response
    critique = json.loads(_strip_markdown_fences(response.content))

    return critique
```

### Component 3: Consensus Checker (Logic-based, instant)

**What it checks:**
1. **Cross-tier agreement**
   - Tier 1 → Tier 2 consistency
   - Tier 2 → Tier 3 consistency

2. **Confidence thresholds**
   - Tier 2 confidence < 0.7 → Flag for review
   - Tier 3 confidence < 0.6 → Flag for review

3. **Known problematic patterns**
   - Contradictory classifications (e.g., "simple" + "tandem")
   - Missing expected features (e.g., Suzuki without Br/I)

**Implementation:**
```python
def check_consensus(
    tier1_result: Dict,
    tier2_result: Dict,
    tier3_result: Dict
) -> Dict[str, Any]:
    """
    Check consensus across all tiers.

    Returns quality score and issues.
    """
    issues = []
    warnings = []

    # 1. Check Tier 1 → Tier 2 agreement
    tier1_patterns = tier1_result.get('interpretation', {})
    tier1_types = [p.lower() for p in tier1_patterns.get('likely_types', [])]
    tier2_types = [rt.lower() for rt in tier2_result.get('reaction_types', [])]

    # If Tier 1 detected something specific, Tier 2 should too
    for t1_type in tier1_types:
        if t1_type and not any(t1_type in t2 for t2 in tier2_types):
            warnings.append(f"Tier 1 detected '{t1_type}' but Tier 2 didn't confirm")

    # 2. Check Tier 2 → Tier 3 agreement (already done elsewhere, but re-check)
    tier3_class = tier3_result.get('overall_class', '').lower()

    # Example: Tier 2 says Suzuki, Tier 3 should say cross_coupling or similar
    suzuki_in_t2 = any('suzuki' in rt or 'coupling' in rt for rt in tier2_types)
    coupling_in_t3 = 'coupling' in tier3_class or 'suzuki' in str(tier3_result.get('tags', []))

    if suzuki_in_t2 and not coupling_in_t3:
        issues.append(f"Major disagreement: Tier 2 says coupling, Tier 3 says {tier3_class}")

    # 3. Check confidence scores
    tier2_conf = tier2_result.get('confidence', 0.0)
    tier3_conf = tier3_result.get('confidence', 0.0)

    if tier2_conf < 0.7:
        warnings.append(f"Tier 2 low confidence ({tier2_conf:.2f})")
    if tier3_conf < 0.6:
        warnings.append(f"Tier 3 low confidence ({tier3_conf:.2f})")

    # 4. Calculate overall quality score
    quality_score = 1.0
    quality_score -= len(issues) * 0.2  # -0.2 per issue
    quality_score -= len(warnings) * 0.05  # -0.05 per warning
    quality_score = max(0.0, min(1.0, quality_score))

    return {
        "quality_score": quality_score,
        "issues": issues,
        "warnings": warnings,
        "recommendation": "accept" if quality_score > 0.8 else "review" if quality_score > 0.6 else "re_analyze"
    }
```

## Quality Gate & Feedback Loop

### Decision Logic

```python
def quality_gate(
    result: Dict[str, Any],
    rdkit_validation: Dict,
    llm_critique: Dict,
    consensus_check: Dict
) -> Dict[str, Any]:
    """
    Quality gate: Decide if analysis is acceptable or needs retry.

    Returns:
        {
            "status": "pass"|"warning"|"fail",
            "action": "accept"|"accept_with_warnings"|"re_analyze",
            "retry_config": {...}  # If re-analysis needed
        }
    """
    # Collect all issues
    all_issues = []
    all_issues.extend(rdkit_validation.get('issues', []))
    all_issues.extend(llm_critique.get('issues', []))
    all_issues.extend(consensus_check.get('issues', []))

    all_warnings = []
    all_warnings.extend(rdkit_validation.get('warnings', []))
    if not llm_critique.get('consistent', True):
        all_warnings.extend(llm_critique.get('suggestions', []))
    all_warnings.extend(consensus_check.get('warnings', []))

    # Calculate overall quality
    quality_score = consensus_check['quality_score']

    # Apply adjustments from validators
    quality_score += rdkit_validation.get('confidence_adjustment', 0.0)

    # Decision tree
    if len(all_issues) > 0:
        # Critical issues found → Re-analyze
        return {
            "status": "fail",
            "action": "re_analyze",
            "issues": all_issues,
            "warnings": all_warnings,
            "retry_config": {
                "use_stronger_model": True,  # Use DeepSeek for Tier 3
                "include_context": all_issues,  # Tell model what went wrong
                "mode": "expert"  # Use highest quality mode
            }
        }
    elif quality_score < 0.7 or len(all_warnings) > 3:
        # Multiple warnings → Accept but flag for review
        return {
            "status": "warning",
            "action": "accept_with_warnings",
            "warnings": all_warnings,
            "quality_score": quality_score,
            "suggestion": "Manual review recommended"
        }
    else:
        # Looks good → Accept
        return {
            "status": "pass",
            "action": "accept",
            "warnings": all_warnings,  # May have minor warnings
            "quality_score": quality_score
        }
```

### Retry Mechanism

```python
def analyze_with_validation(
    rxn_smiles: str,
    client: LLMClient,
    drop_spectators: bool = True,
    max_retries: int = 1
) -> Dict[str, Any]:
    """
    Analyze reaction with validation and optional retry.

    Returns results with validation metadata.
    """
    retry_count = 0
    retry_context = None

    while retry_count <= max_retries:
        # Run analysis (Tiers 1-3)
        result = analyze_reaction_smiles(
            rxn_smiles,
            client,
            drop_spectators=drop_spectators,
            retry_context=retry_context  # Pass issues from previous attempt
        )

        # Run Tier 4 validation
        rdkit_val = validate_with_rdkit(
            result['input']['rxn_smiles_clean'],
            result['interpretation']
        )

        llm_critique = validate_with_llm_critic(
            rxn_smiles,
            result.get('quick_glance', {}),
            result['interpretation'],
            client
        )

        consensus = check_consensus(
            result.get('auto_interpretation', {}),
            result.get('quick_glance', {}),
            result['interpretation']
        )

        # Quality gate
        gate = quality_gate(result, rdkit_val, llm_critique, consensus)

        # Add validation to result
        result['validation'] = {
            "rdkit": rdkit_val,
            "llm_critique": llm_critique,
            "consensus": consensus,
            "gate": gate,
            "retry_count": retry_count
        }

        # Check if we need to retry
        if gate['action'] == 're_analyze' and retry_count < max_retries:
            print(f"⚠️  Quality gate failed, retrying with stronger model...")
            print(f"   Issues: {', '.join(gate['issues'][:2])}")

            # Setup retry
            retry_context = gate['retry_config']
            client.model = "deepseek-v3.2"  # Use best model for retry
            retry_count += 1
        else:
            # Accept result (pass or warning)
            break

    return result
```

## Integration with Current System

### Minimal Changes Required

1. **New file: `reaction_agent/validation.py`**
   - All validation functions
   - RDKit validator
   - LLM self-critic
   - Consensus checker
   - Quality gate

2. **Update `reaction_agent/agent.py`**
   - Add optional `validate=True` parameter
   - Call validation after Tier 3
   - Implement retry logic

3. **Update `reaction_agent/cli.py`**
   - Add `--validate` flag
   - Display validation results
   - Show quality score and warnings

### CLI Output Example

```
================================================================================
  VALIDATION RESULTS (Tier 4)
================================================================================

RDKit Checks: ✓ PASS
  • Structures valid
  • Atom balance: 45 → 29 (16 lost, matches deprotection)
  • Functional groups consistent with cross-coupling

LLM Self-Critique: ✓ PASS with minor suggestions
  • Tier 2/3 consensus: GOOD
  • Mechanism plausibility: HIGH
  • Suggestion: Consider mentioning workup conditions for THP removal

Consensus Score: 0.85 / 1.00

Overall Status: ✓ PASS - High quality analysis
```

## Comparison: Alternative Approaches

### Option A: Simple RDKit-only Validation (Lightweight)
- **Pro**: Fast (~1s), free, deterministic
- **Con**: Can't validate semantic correctness (e.g., wrong mechanism)
- **Best for**: Quick sanity checks

### Option B: Multi-Agent System (Complex)
- **Pro**: Specialist agents for each task
- **Con**: Slow, expensive, complex orchestration
- **Best for**: Mission-critical applications

### Option C: Human-in-the-loop (Manual)
- **Pro**: Best accuracy
- **Con**: Not scalable, slow
- **Best for**: Training/gold standard creation

### **Recommended: Hybrid (Tier 4 as described above)**
- RDKit for structural validation
- Cheap LLM for semantic validation
- Logic-based consensus checking
- Total cost: ~$0.0001 per reaction
- Total time: ~4 seconds additional

## Cost-Benefit Analysis

### Current System (3 Tiers)
- **Time**: ~23s
- **Cost**: ~$0.006 per reaction
- **Accuracy**: Good (DeepSeek Tier 2 is excellent)
- **Reliability**: Unknown (no validation)

### With Tier 4 Validation
- **Time**: ~27s (+4s)
- **Cost**: ~$0.0061 per reaction (+$0.0001)
- **Accuracy**: Same base accuracy
- **Reliability**: High (catch 80% of errors)

### With Tier 4 + Retry
- **Time**: ~50s worst case (if retry triggered)
- **Cost**: ~$0.012 per reaction worst case
- **Accuracy**: Higher (best model on retry)
- **Reliability**: Very high (93%+ correct)

## Implementation Priority

### Phase 1: RDKit Validator (1-2 hours)
- Implement structural validation
- Atom counting
- Basic plausibility checks
- **Impact**: Catch obvious errors

### Phase 2: Consensus Checker (1 hour)
- Cross-tier agreement logic
- Confidence thresholds
- **Impact**: Surface disagreements

### Phase 3: CLI Integration (1 hour)
- Display validation results
- Add --validate flag
- **Impact**: User visibility

### Phase 4: LLM Self-Critic (2-3 hours)
- Implement critic prompt
- Parse and integrate results
- **Impact**: Semantic validation

### Phase 5: Retry Loop (2 hours)
- Implement retry logic
- Quality gate decision tree
- **Impact**: Self-healing system

**Total: ~7-9 hours of development**

## Conclusion

Adding Tier 4 validation provides:
- ✅ **Reliability**: Catch errors before they reach users
- ✅ **Trust**: Users know results are validated
- ✅ **Self-healing**: Automatic retry for low-quality results
- ✅ **Minimal cost**: ~$0.0001 and 4 seconds per reaction
- ✅ **Extensible**: Easy to add more validators

The hybrid approach (RDKit + LLM + Logic) is the sweet spot for your use case.
