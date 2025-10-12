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
        "safety": SAFETY_ASSESSMENT,
        "troubleshooting": TROUBLESHOOTING,
        "spectroscopy": SPECTROSCOPY_INTERPRETATION,
        "protocol": PROTOCOL_GENERATION,
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
    }
