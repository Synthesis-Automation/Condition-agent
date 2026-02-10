"""
Chemistry-Specific Prompt Templates
====================================

Reusable prompt templates for common chemistry tasks:
- Reaction condition recommendation
- Retrosynthesis planning
- Mechanism explanation
- Literature search and summarization
- Reagent selection and optimization
- Safety and hazard assessment

Templates use placeholder substitution and can be customized.
"""

from typing import Any, Dict, Optional


class PromptTemplate:
    """Base class for prompt templates with variable substitution."""
    
    def __init__(self, template: str, **defaults):
        """
        Initialize template.
        
        Args:
            template: Template string with {variable} placeholders
            **defaults: Default values for variables
        """
        self.template = template
        self.defaults = defaults
    
    def format(self, **kwargs) -> str:
        """
        Format template with provided variables.
        
        Args:
            **kwargs: Variable values (overrides defaults)
            
        Returns:
            Formatted prompt string
        """
        # Merge defaults with provided kwargs
        values = {**self.defaults, **kwargs}
        return self.template.format(**values)
    
    def __call__(self, **kwargs) -> str:
        """Shortcut for format()."""
        return self.format(**kwargs)


# =============================================================================
# CONDITION RECOMMENDATION
# =============================================================================

CONDITION_RECOMMENDATION = PromptTemplate(
    template="""You are an expert synthetic chemist specializing in reaction optimization.

**Task**: Recommend optimal reaction conditions for the following transformation.

**Reaction SMILES**: {reaction_smiles}
**Reaction Type**: {reaction_type}

**Requirements**:
1. Suggest appropriate catalyst, ligand (if needed), base, and solvent
2. Recommend temperature and reaction time
3. Provide rationale based on literature precedents
4. Consider functional group compatibility
5. Highlight potential side reactions or challenges

**Additional Context**: {context}

Please provide your recommendations in a structured format:
- Catalyst:
- Ligand:
- Base:
- Solvent:
- Temperature:
- Time:
- Rationale:
- Potential Issues:
""",
    context="None provided"
)


# =============================================================================
# RETROSYNTHESIS PLANNING
# =============================================================================

RETROSYNTHESIS = PromptTemplate(
    template="""You are an expert synthetic chemist skilled in retrosynthetic analysis.

**Task**: Propose a retrosynthetic route for the target molecule.

**Target SMILES**: {target_smiles}
**Target Name**: {target_name}

**Constraints**:
- Maximum {max_steps} steps
- Use {availability} starting materials
- Consider {complexity} synthetic complexity

Provide:
1. Full retrosynthetic tree
2. Forward synthetic route
3. Key transformations and reagents
4. Estimated overall yield
5. Potential bottlenecks

Format your response as a step-by-step plan.
""",
    target_name="Unknown",
    max_steps=5,
    availability="commercially available",
    complexity="moderate"
)


# =============================================================================
# MECHANISM EXPLANATION
# =============================================================================

MECHANISM_EXPLANATION = PromptTemplate(
    template="""You are an expert in reaction mechanisms and physical organic chemistry.

**Task**: Explain the mechanism of the following reaction.

**Reaction SMILES**: {reaction_smiles}
**Reaction Type**: {reaction_type}
**Catalyst/Reagents**: {reagents}

Provide:
1. Step-by-step mechanism with electron pushing
2. Key intermediates and transition states
3. Rate-determining step
4. Stereochemical considerations (if applicable)
5. Literature references supporting the mechanism

**Detail Level**: {detail_level}
""",
    reagents="Not specified",
    detail_level="Detailed with all intermediates"
)


# =============================================================================
# LITERATURE SEARCH
# =============================================================================

LITERATURE_SEARCH = PromptTemplate(
    template="""You are a chemistry research assistant specializing in literature analysis.

**Task**: Find and summarize relevant literature for this reaction.

**Reaction SMILES**: {reaction_smiles}
**Reaction Type**: {reaction_type}
**Focus**: {focus}

Provide:
1. Key publications (with DOI if possible)
2. Common reaction conditions used
3. Typical yields and selectivities
4. Recent advances or improvements
5. Industrial applications (if any)

**Search Scope**: {scope}
**Year Range**: {year_range}
""",
    focus="General conditions and scope",
    scope="Academic literature",
    year_range="Last 10 years"
)


# =============================================================================
# REAGENT SELECTION
# =============================================================================

REAGENT_SELECTION = PromptTemplate(
    template="""You are an expert in chemical reagent selection and optimization.

**Task**: Recommend the best {reagent_type} for this transformation.

**Reaction SMILES**: {reaction_smiles}
**Reaction Type**: {reaction_type}
**Current Conditions**: {current_conditions}

**Selection Criteria**:
- Selectivity: {selectivity_priority}
- Cost: {cost_priority}
- Safety: {safety_priority}
- Scalability: {scale_priority}

Compare at least {num_options} options with pros/cons for each.
Provide final recommendation with justification.
""",
    reagent_type="base",
    current_conditions="Not specified",
    selectivity_priority="High",
    cost_priority="Medium",
    safety_priority="High",
    scale_priority="Medium",
    num_options=3
)


# =============================================================================
# REAGENT REGISTRY REVIEW
# =============================================================================

REAGENT_REGISTRY_REVIEW = PromptTemplate(
    template="""You are an expert chemical data curator helping vet entries for a reagent taxonomy database.

Review the proposed registry entry and respond in JSON only. Do not include any additional commentary.

## Reagent Identity
- Name: {name}
- CAS: {cas}
- Synonyms: {synonyms}
- Abbreviations: {abbreviations}

## Deterministic Resolution
- Resolver source: {resolver_source}
- Resolver name: {resolver_name}
- Resolver SMILES: {resolver_smiles}
- Token hits for family: {family_tokens}
- Default family used: {used_default_family}
- Debug notes: {debug_log}

## Proposed Assignment
- Role: {role}
- Family ID: {family_id}
- Family definition: {family_definition}
- Family notes: {family_notes}
- Family keywords: {family_keywords}
- Family examples: {family_examples}

## Existing Conflicts
- Already in registry (same role): {existing_same_role}
- Already in registry (other roles): {existing_other_roles}

## Task
Evaluate the proposed assignment carefully. Pay special attention to:
1. **Role accuracy**: If current role is "other_reagent", STRONGLY recommend a more specific role if applicable
2. **Family precision**: Suggest the most specific family that matches the reagent's chemistry
3. **Field reliability**: Identify critical missing fields (e.g., metal, oxidation_states, basicity, etc.)
4. **Synonym completeness**: Add any well-known synonyms or trade names

Respond with strict JSON using this schema:
{{
  "status": "confirm" | "needs_review" | "reject",
  "proposed_role": "string",
  "proposed_family": "string",
  "confidence": 0-1 float,
  "justification": "short explanation",
  "alerts": ["list of warnings or required actions"],
  "suggested_synonyms": ["list of additional synonyms (if any)"],
  "field_suggestions": {{
    "field_name": "suggested_value",
    "field_name2": "suggested_value2"
  }}
}}

If you are uncertain or data is inconsistent, set status to "needs_review" and include alerts explaining why.
If the deterministic assignment appears wrong, set status to "reject" and provide the corrected role/family you recommend.
If current role is "other_reagent", ALWAYS propose a more specific role unless truly uncertain.
""",
    synonyms="None provided",
    abbreviations="[]",
    resolver_source="Unknown",
    resolver_name="Unknown",
    resolver_smiles="Unknown",
    family_tokens="[]",
    used_default_family=False,
    debug_log="[]",
    role="Unknown",
    family_id="Unknown",
    family_definition="Definition not available",
    family_notes="Notes not available",
    family_keywords="[]",
    family_examples="[]",
    existing_same_role="None",
    existing_other_roles="None",
)


# =============================================================================
# SAFETY ASSESSMENT
# =============================================================================

SAFETY_ASSESSMENT = PromptTemplate(
    template="""You are a chemistry safety expert with expertise in hazard assessment.

**Task**: Assess safety considerations for this reaction.

**Reaction SMILES**: {reaction_smiles}
**Reagents**: {reagents}
**Solvents**: {solvents}
**Scale**: {scale}

Evaluate:
1. Chemical hazards (toxicity, reactivity, flammability)
2. Process hazards (exotherms, gas evolution, pressure)
3. Required PPE and engineering controls
4. Safe handling and storage
5. Waste disposal considerations
6. Emergency response procedures

**Risk Tolerance**: {risk_level}
**Setting**: {setting}
""",
    reagents="Not specified",
    solvents="Not specified",
    scale="Laboratory (mg-g scale)",
    risk_level="Standard academic lab",
    setting="Academic research laboratory"
)


# =============================================================================
# REACTION TROUBLESHOOTING
# =============================================================================

TROUBLESHOOTING = PromptTemplate(
    template="""You are an expert synthetic chemist specializing in reaction troubleshooting.

**Problem**: {problem_description}

**Reaction Details**:
- SMILES: {reaction_smiles}
- Type: {reaction_type}
- Current Conditions: {current_conditions}
- Observed Result: {observed_result}
- Expected Result: {expected_result}

**Analytical Data Available**: {analytical_data}

Provide:
1. Likely causes of the problem
2. Diagnostic experiments to identify the issue
3. Suggested modifications to reaction conditions
4. Alternative approaches if needed
5. Step-by-step troubleshooting plan

Priority: {priority}
""",
    problem_description="Reaction not proceeding as expected",
    current_conditions="Not specified",
    observed_result="No conversion",
    expected_result="High yield of product",
    analytical_data="TLC, NMR",
    priority="High - need solution quickly"
)


# =============================================================================
# SPECTROSCOPY INTERPRETATION
# =============================================================================

SPECTROSCOPY_INTERPRETATION = PromptTemplate(
    template="""You are an expert in spectroscopic analysis and structure elucidation.

**Task**: Interpret spectroscopic data to confirm product structure.

**Expected Product SMILES**: {product_smiles}
**Reaction**: {reaction_smiles}

**Spectroscopic Data**:
{spectroscopy_data}

Analyze:
1. Is the data consistent with the expected structure?
2. Key diagnostic peaks/signals
3. Alternative structures if data doesn't match
4. Additional experiments needed for confirmation
5. Confidence level in structure assignment

**Available Techniques**: {available_techniques}
""",
    reaction_smiles="Not provided",
    spectroscopy_data="Paste NMR, MS, IR data here",
    available_techniques="NMR (1H, 13C), MS, IR"
)


# =============================================================================
# PROTOCOL GENERATION
# =============================================================================

PROTOCOL_GENERATION = PromptTemplate(
    template="""You are an expert synthetic chemist writing experimental procedures.

**Task**: Write a detailed experimental protocol for this reaction.

**Reaction SMILES**: {reaction_smiles}
**Reaction Type**: {reaction_type}
**Recommended Conditions**: {conditions}
**Scale**: {scale}

Generate a complete protocol including:
1. Materials and equipment needed
2. Step-by-step procedure
3. Reaction setup and monitoring
4. Workup procedure
5. Purification method
6. Characterization requirements
7. Safety precautions
8. Expected yield and appearance

**Format**: {format_style}
**Detail Level**: {detail_level}
""",
    conditions="Use recommended conditions from previous analysis",
    scale="1 mmol",
    format_style="Journal-style experimental section",
    detail_level="Detailed enough for reproduction by trained chemist"
)


# =============================================================================
# REACTION FEATURIZATION REVIEW
# =============================================================================

REACTION_FEATURIZATION_REVIEW = PromptTemplate(
    template="""You are a chemistry quality-control assistant for a deterministic reaction featurization pipeline.

Task: Review the deterministic output and suggest a corrected reaction type only when confidence is low or edge-case signals are present.

Reaction SMILES: {reaction_smiles}
Normalized reaction: {normalized_reaction}

Deterministic reaction type: {deterministic_reaction_type}
Deterministic confidence: {deterministic_confidence}
Mapping warning: {mapping_warning}

CRK raw key: {reaction_key_raw}
CRK final key: {reaction_key}

Reacted motifs: {reacted_motifs}
Formed motifs: {formed_motifs}
Spectator motifs: {spectator_motifs}
Product broad tags: {product_broad_tags}
Product reactive motifs: {product_motifs_reactive}

Allowed taxonomy reaction types:
{reaction_type_candidates}

Rules:
1. Suggest ONLY one reaction type from the allowed taxonomy list, or "Unknown".
2. Prefer deterministic result when evidence is weak.
3. Do not invent new reaction labels.
4. Keep rationale short and chemistry-based.

Respond with ONLY valid JSON:
{{
  "suggested_reaction_type": "reaction_type_id_or_Unknown",
  "confidence": 0.0,
  "rationale": "short explanation",
  "requires_human_review": false,
  "uncertainty_flags": ["optional_flag"]
}}
""",
    reaction_smiles="",
    normalized_reaction="",
    deterministic_reaction_type="Unknown",
    deterministic_confidence="0.0",
    mapping_warning="None",
    reaction_key_raw="None",
    reaction_key="None",
    reacted_motifs="[]",
    formed_motifs="[]",
    spectator_motifs="[]",
    product_broad_tags="[]",
    product_motifs_reactive="[]",
    reaction_type_candidates="Unknown",
)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_template(name: str) -> PromptTemplate:
    """
    Get a prompt template by name.
    
    Args:
        name: Template name (case-insensitive)
        
    Returns:
        PromptTemplate instance
        
    Raises:
        KeyError: If template name not found
    """
    templates = {
        "condition_recommendation": CONDITION_RECOMMENDATION,
        "retrosynthesis": RETROSYNTHESIS,
        "mechanism": MECHANISM_EXPLANATION,
        "literature": LITERATURE_SEARCH,
        "reagent_selection": REAGENT_SELECTION,
        "reagent_registry_review": REAGENT_REGISTRY_REVIEW,
        "safety": SAFETY_ASSESSMENT,
        "troubleshooting": TROUBLESHOOTING,
        "spectroscopy": SPECTROSCOPY_INTERPRETATION,
        "protocol": PROTOCOL_GENERATION,
        "reagent_role_classification": REAGENT_ROLE_CLASSIFICATION,
        "reagent_field_assignment": REAGENT_FIELD_ASSIGNMENT,
        "reagent_entry_verification": REAGENT_ENTRY_VERIFICATION,
        "reaction_featurization_review": REACTION_FEATURIZATION_REVIEW,
    }
    
    key = name.lower().replace(" ", "_").replace("-", "_")
    if key not in templates:
        available = ", ".join(templates.keys())
        raise KeyError(f"Template '{name}' not found. Available: {available}")
    
    return templates[key]


def list_templates() -> Dict[str, str]:
    """
    List all available templates with descriptions.
    
    Returns:
        Dict mapping template names to descriptions
    """
    return {
        "condition_recommendation": "Recommend optimal reaction conditions",
        "retrosynthesis": "Plan retrosynthetic routes",
        "mechanism": "Explain reaction mechanisms",
        "literature": "Search and summarize literature",
        "reagent_selection": "Select optimal reagents",
        "safety": "Assess safety and hazards",
        "troubleshooting": "Debug problematic reactions",
        "spectroscopy": "Interpret spectroscopic data",
        "protocol": "Generate experimental protocols",
        "reagent_role_classification": "Classify reagent role from chemistry",
        "reagent_field_assignment": "Assign family and fields for reagent role",
        "reagent_entry_verification": "Verify reagent entry for errors",
        "reaction_featurization_review": "Review uncertain reaction featurization output",
    }


# =============================================================================
# REAGENT TAXONOMY CLASSIFICATION (New LLM Workflow)
# =============================================================================

REAGENT_ROLE_CLASSIFICATION = PromptTemplate(
    template="""You are an expert chemical reagent classifier. Classify this reagent into the correct role category.

## Reagent Information
- Name: {name}
- CAS: {cas}
- SMILES: {smiles}
- Molecular Formula: {molecular_formula}
- Synonyms: {synonyms}

## Available Roles

**ligand**: Phosphines, NHCs, diimines, bidentate/multidentate donor ligands for metal coordination
**metal_catalyst**: Metal salts or pre-ligated complexes that provide the catalytic metal source (e.g., Pd(OAc)2, Pd(PPh3)4)
**base**: Br酶nsted or Lewis bases - amides, alkoxides, carbonates, phosphazenes, superbases
**acid**: Mineral acids, sulfonic acids, Lewis acids used as activators or promoters
**condensation_agent**: Carbodiimides, uronium salts, phosphonium activators for amide/ester formation
**oxidant**: Terminal oxidants and co-oxidants (peroxides, hypervalent iodine, Oxone, O2)
**reductant**: Hydrides, silanes, metal powders, organic electron donors
**additive**: Phase-transfer catalysts, halide scavengers, fluoride sources, reaction modulators
**solvent**: Reaction media categorized by polarity, coordination ability, proticity
**organo_catalyst**: Small-molecule organocatalysts (cinchona, phosphoric acids, NHCs, thioureas)
**enzyme**: Biocatalysts (isolated enzymes or whole-cell systems)
**other_reagent**: Generic reagents that don't fit above categories (use as LAST RESORT)

## Task

Analyze the chemical structure and properties to determine the most appropriate role.

Respond with ONLY valid JSON (no markdown, no explanation):
{{
  "role": "base",
  "confidence": 0.95,
  "reasoning": "Tertiary aliphatic amine with strong Br酶nsted basicity"
}}

IMPORTANT:
- Choose the MOST SPECIFIC role that fits
- Use "other_reagent" ONLY if truly uncertain
- Confidence should be 0.0-1.0
- Reasoning should be chemistry-focused (1-2 sentences)
""",
    name="Unknown",
    cas="Unknown",
    smiles="Unknown",
    molecular_formula="Unknown",
    synonyms="None",
)


REAGENT_FIELD_ASSIGNMENT = PromptTemplate(
    template="""You are an expert chemical database curator. Assign this reagent to the correct family and populate all required fields.

## Reagent Information
- Name: {name}
- CAS: {cas}
- SMILES: {smiles}
- Molecular Formula: {molecular_formula}
- Role: {role}

## Available Families for Role '{role}'

{families_description}

## Required Fields for Role '{role}'

{fields_schema}

## Examples from Database

{examples}

## Task

1. Select the most specific family that matches this reagent's chemistry
2. Assign values to ALL required fields based on chemical properties
3. Suggest abbreviations and additional synonyms if commonly known

Respond with ONLY valid JSON (no markdown, no explanation):
{{
  "family": "family_id",
  "fields": {{
    "field1": "value1",
    "field2": "value2"
  }},
  "abbreviations": ["ABBR1", "ABBR2"],
  "additional_synonyms": ["synonym1"],
  "confidence": 0.92,
  "reasoning": "Brief chemistry-based justification"
}}

IMPORTANT:
- All field values must match allowed options in schema
- Include ALL required fields (no null/missing values)
- Abbreviations should be well-known (not made up)
- Confidence should reflect certainty of classification
""",
    name="Unknown",
    cas="Unknown",
    smiles="Unknown",
    molecular_formula="Unknown",
    role="Unknown",
    families_description="No families available",
    fields_schema="No fields required",
    examples="No examples available",
)


REAGENT_ENTRY_VERIFICATION = PromptTemplate(
    template="""You are a quality control reviewer for a chemical reagent database. Check this proposed entry for errors and inconsistencies.

## Proposed Entry

{entry_json}

## Verification Checklist

1. **Chemical accuracy**: Does SMILES match the described role and properties?
2. **Field consistency**: Are field values logically consistent? (e.g., "superbase" + "weak nucleophilicity" is suspicious)
3. **Oxidation states**: For metals, are oxidation states chemically reasonable?
4. **Missing information**: Are there obvious gaps in critical fields?
5. **Obvious mistakes**: Wrong metal element, impossible coordination, etc.

## Task

Review the entry and identify any errors or warnings.

Respond with ONLY valid JSON (no markdown, no explanation):
{{
  "approved": true,
  "issues": [
    {{
      "severity": "error",
      "field": "oxidation_states",
      "message": "Pd cannot have oxidation state +5"
    }}
  ],
  "suggestions": [
    "Consider adding pKa value to notes",
    "Verify SMILES stereochemistry"
  ]
}}

IMPORTANT:
- Set approved=false if ANY "error" severity issues exist
- Use "warning" severity for minor issues that don't prevent saving
- Keep messages concise and actionable
- Focus on chemistry correctness, not formatting
""",
    entry_json="{}",
)


MULTI_SOURCE_SYNTHESIS = PromptTemplate(
    template="""You are an expert synthetic chemist analyzing multiple recommendation sources for reaction conditions.

## Reaction

{reaction_smiles}

## ML-based Recommendations (DRFP Similarity with Precedent Database)

Top 3 conditions from structural similarity search:
{ml_conditions}

## Rule-based Recommendations (SCDB Pattern Matching)

Conditions from curated reaction database:
{rule_conditions}

## Protocol-based Recommendations (Literature Procedures)

Conditions from published protocols:
{protocol_conditions}

## User Constraints

{constraints}

## Your Task

Analyze all three sources and synthesize a final recommendation:

1. **Identify consensus**: Where do sources agree? This indicates high confidence.
2. **Explain discrepancies**: Why do sources disagree? Which is more appropriate for THIS substrate?
3. **Assess confidence**: High (all agree), Medium (2/3 agree), Low (all differ)
4. **Consider constraints**: User requirements like scale, cost, air sensitivity
5. **Provide ONE best recommendation** with clear rationale citing specific sources
6. **Suggest 1-2 backup conditions** for if the main recommendation fails

## Response Format

Respond with ONLY valid JSON (no markdown fences):
{{
  "consensus_analysis": {{
    "catalyst": {{
      "agreement": "high|medium|low",
      "consensus_value": "Pd(PPh3)4 or null if no consensus",
      "notes": "Why sources agree/disagree"
    }},
    "solvent": {{"agreement": "...", "consensus_value": "...", "notes": "..."}},
    "temperature": {{"agreement": "...", "consensus_value": "...", "notes": "..."}},
    "base": {{"agreement": "...", "consensus_value": "...", "notes": "..."}}
  }},
  "confidence_level": "high|medium|low",
  "confidence_reasoning": "Why this confidence level based on source agreement",
  "recommended_condition": {{
    "catalyst": "Best catalyst choice",
    "ligand": "Ligand if applicable",
    "solvent": "Best solvent",
    "temperature": "Best temperature",
    "base": "Best base",
    "additive": "Any additives needed",
    "rationale": "2-3 sentence explanation citing which sources support this and why it's best given the substrate and constraints"
  }},
  "backup_conditions": [
    {{
      "catalyst": "...",
      "solvent": "...",
      "temperature": "...",
      "base": "...",
      "when_to_use": "Specific scenario where this is better than main recommendation"
    }}
  ],
  "warnings": [
    "Any red flags, functional group incompatibilities, or practical concerns"
  ],
  "source_comparison": {{
    "ml_contribution": "What valuable insight did ML precedents provide?",
    "rule_contribution": "What did rule-based database add?",
    "protocol_contribution": "What did literature protocols add?"
  }}
}}

CRITICAL:
- Cite specific sources when explaining choices (e.g., "ML precedents show..." or "SCDB rule indicates...")
- If all sources disagree, confidence should be LOW and warnings should mention high uncertainty
- Consider functional groups in the substrate when choosing between conflicting recommendations
- User constraints MUST be respected in final recommendation
""",
    reaction_smiles="",
    ml_conditions="No ML recommendations available",
    rule_conditions="No rule-based recommendations available",
    protocol_conditions="No protocol recommendations available",
    constraints="None specified",
)

MULTI_SOURCE_SYNTHESIS_V2 = PromptTemplate(
    template="""Expert chemist: Synthesize recommendations for this reaction.

## REACTION: {reaction_smiles}

## SOURCES
ML (DRFP): {ml_conditions}
Rule (SCDB): {rule_conditions}
Protocol (Lit): {protocol_conditions}

## CONSTRAINTS: {constraints}

## TASK
Analyze sources 鈫?Synthesize ONE best recommendation with backups.

**Confidence Thresholds:**
- HIGH: All 3 sources agree on catalyst+solvent AND ML similarity >0.80
- MEDIUM: 2/3 sources agree OR ML similarity 0.65-0.80
- LOW: Sources disagree OR ML similarity <0.65

**Chemistry Guidelines:**
- Electron-poor aryl 鈫?electron-rich ligand (dppf, XPhos, SPhos)
- Heteroaryl halides 鈫?chelating ligands (watch Pd coordination)
- Nitro groups 鈫?avoid H2 (reduction risk); monitor proto-debromination
- Sterically hindered 鈫?bidentate/bulky ligands (dppf, SPhos)
- Aryl chlorides 鈫?strong base + high temp OR bulky ligand

**Decision Tree for Backups:**
- Backup 1: If main <30% conv after 6h (different ligand/solvent)
- Backup 2: If Backup 1 <20% after 12h (different catalyst family)

## OUTPUT (JSON only, no markdown)
{{
  "consensus_analysis": {{
    "catalyst": {{"agreement": "high|medium|low", "consensus_value": "...", "notes": "..."}},
    "solvent": {{"agreement": "...", "consensus_value": "...", "notes": "..."}},
    "temperature": {{"agreement": "...", "consensus_value": "...", "notes": "..."}},
    "base": {{"agreement": "...", "consensus_value": "...", "notes": "..."}}
  }},
  "confidence_level": "high|medium|low",
  "confidence_reasoning": "Source agreement + ML similarity assessment",
  "recommended_condition": {{
    "catalyst": "...",
    "ligand": "...",
    "solvent": "...",
    "temperature": "...",
    "base": "...",
    "additive": "...",
    "rationale": "Cite sources + substrate chemistry + constraints"
  }},
  "backup_conditions": [
    {{"catalyst": "...", "solvent": "...", "temperature": "...", "base": "...", 
      "when_to_use": "<30% conv: try this ligand/solvent"}}
  ],
  "warnings": ["Functional group risks", "Practical concerns"],
  "source_comparison": {{
    "ml_contribution": "...",
    "rule_contribution": "...",
    "protocol_contribution": "..."
  }}
}}

RULES: Respect constraints. Cite sources. Use chemistry guidelines for substrate-specific warnings.
""",
    reaction_smiles="",
    ml_conditions="No ML recommendations available",
    rule_conditions="No rule-based recommendations available",
    protocol_conditions="No protocol recommendations available",
    constraints="None specified",
)


# =============================================================================
# RULE BUILDER EXTRACTION
# =============================================================================

RULE_BUILDER_EXTRACTION = PromptTemplate(
    template="""You are a synthetic chemistry knowledge engineer tasked with
generating structured rule-database content for the "{family}" family.

Reference reactions:
{reference_block}

Protocol trends and notes:
{protocol_text}

Focus / constraints:
{focus}

Produce STRICT JSON with these keys:
{{
  "notes": "Concise rationale for this rule database",
  "scope": {{"scope_type": "...", "compatible_functional_groups": [], "incompatible_functional_groups": []}},
  "applies_if": {{"all": [], "any": []}},
  "default_rule": {{
    "id": "identifier",
    "description": "short text",
    "conditions": {{"pd_source": "...", "ligand": "...", "base": "...", "temperature_C": "...", "time_h": "..."}}
  }},
  "base_rules": [
    {{
      "id": "lower_snake_case_id",
      "name": "Readable name",
      "description": "1-2 sentence guidance",
      "reactant_features": {{"all": [], "any": []}},
      "conditions": {{"pd_source": "...", "ligand": "...", "base": "...", "temperature_C": "...", "time_h": "...", "solvent": "...", "additive": "... optional"}},
      "priority": 0
    }}
  ],
  "modifiers": [
    {{
      "id": "modifier_id",
      "when": ["feature_or_symptom", "symptom:hydrodehalogenation_observed"],
      "suggest": "Actionable recommendation",
      "rationale": "Optional explanation"
    }}
  ]
}}

Rules:
- Create at least 1 base rule and at most {max_base_rules}.
- Each conditions dict must include catalyst/precatalyst, ligand (if applicable), base, solvent, temperature, and time.
- Use lower_snake_case for IDs, camel case is not allowed.
- Use short lists for applies_if/reactant_features; do NOT embed prose.
- Modifiers should use "symptom:" prefixes when referencing lab observations.
- Return JSON ONLY. No markdown fences, commentary, or trailing text.
""",
    reference_block="- None provided",
    protocol_text="(no protocol text)",
    focus="None specified",
    max_base_rules=4,
)
