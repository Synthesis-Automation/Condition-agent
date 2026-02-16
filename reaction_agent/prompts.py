"""
Prompt Templates for Reaction SMILES Analysis

Contains the LLM prompt template for interpreting bond changes and
classifying reaction mechanisms.
"""

from typing import Dict


class PromptTemplate:
    """Simple prompt template with variable substitution."""

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
        values = {**self.defaults, **kwargs}
        return self.template.format(**values)

    def __call__(self, **kwargs) -> str:
        """Shortcut for format()."""
        return self.format(**kwargs)


REACTION_SMILES_ANALYSIS = PromptTemplate(
    template="""You are a reaction analysis engine specialized in interpreting chemical transformations.

**CRITICAL INSTRUCTIONS**:
1. Use ONLY the provided tool_facts (bond changes, atom mapping quality)
2. If evidence is insufficient, add warnings and lower confidence
3. Do NOT invent reagents, catalysts, or conditions
4. Reference bond changes by their IDs (e.g., BC1, BC2)

**INPUT DATA**:

## Raw and Clean SMILES
- Raw reaction SMILES: {rxn_smiles_raw}
- Clean reaction SMILES: {rxn_smiles_clean}
- Spectators detected: {spectators}
- Parse warnings: {parse_warnings}

## Tool Facts (Deterministic)
- Mapped reaction: {mapped_rxn_smiles}
- Mapping quality: {mapping_qc}

### Bond Changes
{bond_changes_text}

### Reaction Center Atoms
Atom map numbers: {reaction_center_atoms}

**YOUR TASK**:

Analyze this reaction and provide a structured interpretation following this exact JSON schema:

{{
  "overall_class": "<one of: cross_coupling, nucleophilic_substitution, electrophilic_addition, nucleophilic_addition, elimination, condensation, cycloaddition, rearrangement, oxidation, reduction, protection_deprotection, other>",
  "tags": ["<specific tags like SNAr, Suzuki, Wittig, etc>"],
  "roles": {{
    "electrophile": "<SMILES or description of electrophilic component>",
    "nucleophile": "<SMILES or description of nucleophilic component>",
    "leaving_group": "<SMILES or description of leaving group>"
  }},
  "events": [
    {{
      "event_id": "E1",
      "event_type": "<one of: bond_formation, bond_cleavage, bond_order_change, substitution, addition, elimination, rearrangement, redox, other>",
      "bond_change_refs": ["BC1", "BC2"],
      "short_rationale": "<brief mechanistic explanation referencing the bond changes>",
      "confidence": <0.0 to 1.0>
    }}
  ],
  "mechanism_summary": [
    "<step 1 description>",
    "<step 2 description>"
  ],
  "warnings": ["<any concerns about data quality or analysis uncertainty>"],
  "confidence": <overall confidence 0.0 to 1.0>
}}

**GUIDELINES**:
- If mapping quality is poor (mapping_qc.ok = false), set confidence <= 0.3 and add warnings
- If there are >5 bond changes, this may be a tandem reaction - create multiple events
- The roles dict should adapt to reaction type (e.g., for oxidation: use "substrate", "oxidant", etc.)
- Be specific with mechanism_summary but keep steps concise (1-2 sentences each)
- Confidence should reflect both mapping quality AND mechanistic clarity

Return ONLY the JSON object, no markdown fences or additional text.""",
    rxn_smiles_raw="",
    rxn_smiles_clean="",
    spectators="None",
    parse_warnings="None",
    mapped_rxn_smiles="Not available",
    mapping_qc="{}",
    bond_changes_text="No bond changes detected",
    reaction_center_atoms="[]",
)


DIRECT_SMILES_ANALYSIS = PromptTemplate(
    template="""You are a reaction analysis engine specialized in interpreting chemical transformations.

⚠️ **MAPPING FAILURE MODE** ⚠️

The atom mapping tool has COMPLETELY FAILED for this reaction (0 bond changes detected, confidence {mapping_confidence:.3f}).
This does NOT mean the reaction is invalid - it means you must analyze the SMILES strings DIRECTLY using chemical reasoning.

**IGNORE the deterministic analysis** and focus on:
1. Identifying functional groups in reactants
2. Identifying functional groups in products
3. Inferring what bonds MUST have changed based on structure differences
4. Using your knowledge of organic chemistry patterns

**INPUT DATA**:

## Reactants SMILES
{reactants_smiles}

## Products SMILES
{products_smiles}

## Full Reaction
{rxn_smiles_clean}

**YOUR TASK**:

Analyze this reaction by comparing reactant and product structures directly. Provide detailed mechanistic reasoning using this JSON schema:

{{
  "overall_class": "<one of: cross_coupling, nucleophilic_substitution, electrophilic_addition, nucleophilic_addition, elimination, condensation, cycloaddition, rearrangement, oxidation, reduction, annulation, defluorination, other>",
  "tags": ["<specific tags like SNAr, Suzuki, cyclization, etc>"],
  "functional_groups": {{
    "reactants": ["<FG1: description>", "<FG2: description>"],
    "products": ["<FG1: description>", "<FG2: description>"]
  }},
  "inferred_bond_changes": {{
    "bonds_likely_broken": ["<description>"],
    "bonds_likely_formed": ["<description>"],
    "reasoning": "<explain how you inferred these from SMILES comparison>"
  }},
  "mechanism_summary": [
    "<step 1: describe transformation>",
    "<step 2: describe intermediate>",
    "<step 3: describe final product formation>"
  ],
  "reaction_centers": [
    {{
      "center_id": 1,
      "description": "<what changes at this center>",
      "evidence": "<functional group changes observed>"
    }}
  ],
  "literature_precedent": "<if this matches known chemistry, describe it>",
  "warnings": ["mapping_failed_used_direct_analysis"],
  "confidence": <0.5 to 0.9 based on how clear the transformation is>
}}

**CRITICAL**:
- Compare the SMILES strings DIRECTLY - what appears in products that wasn't in reactants?
- What disappears from reactants?
- Look for ring formations, ring openings, functional group interconversions
- Use your chemical knowledge - what reagents typically cause such changes?
- Be detailed in mechanism_summary - explain the full transformation pathway
- Confidence should be 0.5-0.9 (lower than with good mapping, but not zero)

Return ONLY the JSON object, no markdown fences or additional text.""",
    reactants_smiles="",
    products_smiles="",
    rxn_smiles_clean="",
    mapping_confidence=0.0
)


def get_template() -> PromptTemplate:
    """Get the reaction SMILES analysis prompt template."""
    return REACTION_SMILES_ANALYSIS


def get_direct_smiles_template() -> PromptTemplate:
    """Get the direct SMILES analysis template (for mapping failures)."""
    return DIRECT_SMILES_ANALYSIS
