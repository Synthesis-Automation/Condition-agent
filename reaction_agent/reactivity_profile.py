"""
Rich reactivity profile — comprehensive output of the ReactionReasoningAgent.

Stores everything the agent discovered about a reaction: electrophile/nucleophile
characterization, mechanism analysis, selectivity risks, condition implications,
and taxonomy mapping. Some fields are used immediately for recommendation; others
are preserved for future explanation steps.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class ElectrophileProfile:
    """Characterization of the electrophilic center."""
    center_atom: str = ""                     # e.g. "C6 (ipso carbon)" or description
    hybridization: str = ""                   # "sp2_aryl", "sp2_vinyl", "sp3", "carbonyl"
    leaving_group: str = ""                   # "Br", "Cl", "OTf", "OH", "none"
    leaving_group_quality: str = ""           # "excellent", "good", "moderate", "poor"
    steric_score: float = 0.0                 # 0-10 from steric analysis
    steric_class: str = ""                    # "unhindered", "moderate", "hindered"
    electronic_score: float = 5.0             # 0-10 from electronic analysis
    electronic_class: str = ""                # "electron_poor", "neutral", "electron_rich"
    activation_tags: List[str] = field(default_factory=list)   # ["EWG_para_NO2", "heteroaryl"]
    raw_tool_data: Dict[str, Any] = field(default_factory=dict)


@dataclass
class NucleophileProfile:
    """Characterization of the nucleophilic species."""
    identity: str = ""                        # "primary_amine", "boronic_acid", "alkoxide"
    attacking_atom: str = ""                  # "N", "C", "O", "S"
    hardsoft: str = ""                        # "soft", "hard", "borderline"
    is_also_base: bool = False                # Could trigger elimination
    steric_bulk: str = ""                     # "small", "moderate", "bulky"
    functional_groups: List[str] = field(default_factory=list)
    raw_tool_data: Dict[str, Any] = field(default_factory=dict)


@dataclass
class MechanismAnalysis:
    """Mechanism hypothesis with evidence."""
    primary_class: str = ""                   # "oa_based_coupling", "snar", "sn2", "acyl_sub"
    evidence: List[str] = field(default_factory=list)          # Tool-grounded evidence
    confidence: float = 0.0
    alternative_mechanisms: List[str] = field(default_factory=list)
    snar_feasibility: Optional[Dict[str, Any]] = None
    requires_catalyst: bool = False
    likely_catalyst_metals: List[str] = field(default_factory=list)  # ["Pd", "Ni", "Cu"]
    key_intermediates: List[str] = field(default_factory=list)       # named intermediates
    stepwise: List[str] = field(default_factory=list)                # step-by-step pathway


@dataclass
class TransformationRecord:
    """What bonds changed and what functional groups transformed."""
    bonds_broken: List[str] = field(default_factory=list)   # ["C-Br (sp2, aryl, ipso)"]
    bonds_formed: List[str] = field(default_factory=list)   # ["C-C (sp2-sp2, biaryl)"]
    key_bond_type: str = ""                  # "C-C", "C-N", "C-O", "C-S"
    fg_removed: List[str] = field(default_factory=list)      # ["aryl_bromide"]
    fg_formed: List[str] = field(default_factory=list)       # ["biaryl"]
    redox_change: str = ""                   # "neutral", "oxidation", "reduction"
    raw_bond_analysis: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ReactivityProfile:
    """
    Comprehensive reactivity profile for a reaction.

    This is the main output of ReactionReasoningAgent.analyze().
    All fields are populated by the agent's tool calls and reasoning.
    """

    # --- Input ---
    reaction_smiles: str = ""
    reactant_smiles: List[str] = field(default_factory=list)
    product_smiles: List[str] = field(default_factory=list)

    # --- Core analysis (always populated) ---
    transformation: TransformationRecord = field(default_factory=TransformationRecord)
    electrophile: ElectrophileProfile = field(default_factory=ElectrophileProfile)
    nucleophile: NucleophileProfile = field(default_factory=NucleophileProfile)
    mechanism: MechanismAnalysis = field(default_factory=MechanismAnalysis)

    # --- Extended analysis (for future explanation steps) ---
    selectivity_risks: List[str] = field(default_factory=list)
    competing_pathways: List[str] = field(default_factory=list)
    driving_forces: List[str] = field(default_factory=list)
    condition_implications: Dict[str, Any] = field(default_factory=dict)
    functional_group_inventory: Dict[str, Any] = field(default_factory=dict)
    molecular_descriptors: Dict[str, Any] = field(default_factory=dict)

    # --- Taxonomy mapping (for recommendation) ---
    reaction_type: Optional[str] = None
    reaction_type_confidence: float = 0.0
    reaction_pattern_type: str = ""              # broad pattern: coupling_substitution, condensation, etc.
    named_reaction: str = ""                     # human-readable named reaction or description
    reacted_motifs: List[str] = field(default_factory=list)
    formed_motifs: List[str] = field(default_factory=list)
    taxonomy_reasoning: str = ""

    # --- All-component role assignment ---
    all_roles: Dict[str, str] = field(default_factory=dict)  # {smiles_fragment: role}

    # --- Product verification ---
    product_verification: Dict[str, Any] = field(default_factory=dict)  # expected/confirmed motifs + score

    # --- Missing conditions ---
    missing_conditions: List[str] = field(default_factory=list)  # conditions implied but absent from SMILES

    # --- Tandem / multi-step ---
    is_tandem: bool = False
    reactive_streams: List[Dict[str, Any]] = field(default_factory=list)  # bifurcating pathway breakdown
    tandem_steps: List[Dict[str, Any]] = field(default_factory=list)

    # --- Agent metadata ---
    reasoning_chain: List[Dict[str, Any]] = field(default_factory=list)
    tools_called: List[str] = field(default_factory=list)
    total_tool_calls: int = 0
    llm_model: str = ""
    confidence: float = 0.0
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dictionary."""
        from dataclasses import asdict
        return asdict(self)

    def summary(self) -> str:
        """One-line summary of the profile."""
        parts = []
        if self.reaction_pattern_type:
            parts.append(f"pattern={self.reaction_pattern_type}")
        if self.mechanism.primary_class:
            parts.append(f"mechanism={self.mechanism.primary_class}")
        if self.reaction_type:
            parts.append(f"type={self.reaction_type}")
        if self.reacted_motifs:
            parts.append(f"motifs={self.reacted_motifs}")
        if self.confidence:
            parts.append(f"conf={self.confidence:.2f}")
        return " | ".join(parts) if parts else "(empty profile)"


__all__ = [
    "ReactivityProfile",
    "ElectrophileProfile",
    "NucleophileProfile",
    "MechanismAnalysis",
    "TransformationRecord",
]
