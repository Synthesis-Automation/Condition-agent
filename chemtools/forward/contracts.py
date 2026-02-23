"""
Data contracts for the forward synthesis engine.

ForwardTemplateMatch  — a forward SMARTS template applicable to a reactant pair
ReactantProfile       — functional-group inventory of a single reactant
ProductPrediction     — a predicted product with scoring metadata
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class ReactantProfile:
    """Functional-group inventory and reactivity classification for one molecule."""

    smiles: str
    canonical_smiles: str
    molecular_weight: float
    heavy_atoms: int
    functional_groups: List[str]          # e.g. ["aryl_bromide", "ester"]
    fg_categories: Dict[str, List[str]]   # category → [fg names]
    electrophilic_sites: List[str]        # FG or SMARTS labels
    nucleophilic_sites: List[str]
    leaving_group: Optional[str] = None   # "Br", "F", "Cl", "OTf", etc.


@dataclass
class ForwardTemplateMatch:
    """A forward SMARTS template that matches the reactant functional groups."""

    template_name: str
    taxonomy_id: str
    category: str
    forward_smarts: str
    n_reactants: int                      # 1 or 2
    difficulty: float                     # 0.0–1.0
    description: str
    hte_families: List[str]
    notes: str = ""


@dataclass
class ProductPrediction:
    """A predicted product from applying a forward template to reactant(s)."""

    product_smiles: str                   # canonical SMILES
    reactant_a: str
    reactant_b: str                       # empty string for unimolecular reactions
    template_name: str
    taxonomy_id: str
    reaction_smiles: str                  # reactants>>product

    # Scoring components
    hte_yield_proxy: float = 0.0          # avg yield proxy from HTE database (0–100)
    chemoselectivity_penalty: float = 0.0 # deduction for competing templates
    structural_complexity: float = 0.0    # BertzCT of product
    overall_score: float = 0.0            # final ranking score

    # Metadata
    difficulty: float = 0.5
    description: str = ""
    notes: str = ""
    hte_families: List[str] = field(default_factory=list)
    new_stereocenters: int = 0            # flagged for LLM reasoning
    competing_templates: List[str] = field(default_factory=list)
    confidence_label: str = "medium"      # "high" / "medium" / "low"
    all_product_smiles: List[str] = field(default_factory=list)  # all regioisomers
