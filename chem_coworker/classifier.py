"""
TaskClassifier — lightweight rule-based task type detection.

No LLM call needed. Uses regex patterns and keyword heuristics to:
  - Classify the user's intent (ANALYZE, PREDICT, EXPLAIN, LOOKUP, COMPARE, TROUBLESHOOT)
  - Extract any SMILES strings from the query
  - Detect whether the query involves a reaction (has '>>')

The classified task type is passed to the LLM in REASON_PROMPT so it can
calibrate response depth and tool selection without re-deriving it.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional


class TaskType(Enum):
    ANALYZE = "analyze"
    PREDICT = "predict"
    FORWARD_SYNTHESIS = "forward_synthesis"
    EXPLAIN = "explain"
    LOOKUP = "lookup"
    COMPARE = "compare"
    TROUBLESHOOT = "troubleshoot"
    RETROSYNTHESIS = "retrosynthesis"
    GENERAL = "general"


# Regex to extract SMILES-like tokens (simplified: allows organic chemistry chars)
_SMILES_PATTERN = re.compile(
    r"(?<!\w)"               # not preceded by word char
    r"([A-Za-z\[\]()=#@+\-\\/\.:%0-9]{6,})"  # SMILES chars, min length 6
    r"(?!\w)",               # not followed by word char
)

# Reaction SMILES: must contain >>
_REACTION_SMILES_PATTERN = re.compile(r"[^\s]+>>[^\s]*")

# Keywords for each task type (order matters — checked top to bottom)
_TASK_KEYWORDS: List[tuple] = [
    # Retrosynthesis detected first — "synthesize", "make", "route" are strong signals.
    # We distinguish from forward: user provides a target (one molecule), not a reaction SMILES.
    (TaskType.RETROSYNTHESIS, [
        "retrosynthesis", "retrosynthetic", "retro synthesis", "retro-synthesis",
        "how to make", "how do i make", "how do i synthesize", "how to synthesize",
        "how can i make", "how can i prepare", "how can i synthesize",
        "synthesis of", "synthesize this", "synthesize the",
        "prepare this compound", "prepare this molecule",
        "total synthesis", "route to", "synthetic route",
        "starting material for", "from what starting material",
        "what precursors", "what starting materials",
        "build this molecule", "make this compound",
        "suggest a synthesis", "plan a synthesis", "plan a route",
        "disconnection", "disconnect", "synthon",
    ]),
    (TaskType.TROUBLESHOOT, [
        "low yield", "poor yield", "% yield", "side product", "side reaction",
        "didn't work", "did not work", "failed", "problem", "issue",
        "what went wrong", "improve yield", "optimize",
        "what could cause", "what might cause", "caused by",
        "why did", "why didn't", "why does it", "not working", "doesn't work",
    ]),
    (TaskType.COMPARE, [
        " vs ", " versus ", "compare", "difference between",
        "which is better", "prefer", "choose between",
    ]),
    (TaskType.LOOKUP, [
        "find", "list", "search", "what is the pka", "pka of",
        "what is the boiling", "show me", "give me a list",
        "which reagents", "available", "database", "catalog",
        "work for", "work with", "used for", "commonly used",
        "what bases", "what ligands", "what catalysts", "what solvents",
    ]),
    (TaskType.PREDICT, [
        "what condition", "recommend condition", "suggest condition",
        "what catalyst", "what solvent", "what base", "what ligand",
        "what reagent", "what will happen", "predict", "forecast",
        "which conditions", "best conditions",
    ]),
    (TaskType.EXPLAIN, [
        "explain", "why does", "why is", "how does", "mechanism of",
        "what is the mechanism", "describe the mechanism", "role of",
        "purpose of", "what is the role", "why use",
        "how does pd", "why pd", "why palladium",
    ]),
    (TaskType.ANALYZE, [
        "analyze", "analyse", "what reaction", "what is this reaction",
        "identify", "characterize", "what functional group",
        "drug-like", "properties of", "features of",
    ]),
]

_FORWARD_EXPLICIT_CUES = [
    "forward synthesis",
    "forward reaction prediction",
    "predict the product",
    "predict product",
    "major product",
    "expected product",
    "what is the product",
    "what would be the product",
    "product prediction",
]

_FORWARD_REACTANT_CUES = [
    "reactant a",
    "reactant b",
    "two reactants",
    "these reactants",
    "from these reagents",
    "couple these",
    "mixing",
]


@dataclass
class ClassificationResult:
    task_type: TaskType
    task_type_str: str
    has_reaction: bool
    has_molecule: bool
    reaction_smiles: List[str]
    molecule_smiles: List[str]
    all_smiles: List[str]
    primary_smiles: Optional[str]  # Best candidate SMILES for the main subject


class TaskClassifier:
    """
    Rule-based task type classifier and SMILES extractor.

    Usage:
        clf = TaskClassifier()
        result = clf.classify("Recommend conditions for Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
        print(result.task_type)         # TaskType.PREDICT
        print(result.reaction_smiles)   # ['Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1']
    """

    def classify(self, query: str) -> ClassificationResult:
        q_lower = query.lower()

        # Extract SMILES
        reaction_smiles = self._extract_reaction_smiles(query)
        molecule_smiles = self._extract_molecule_smiles(query, reaction_smiles)
        all_smiles = reaction_smiles + molecule_smiles

        has_reaction = bool(reaction_smiles)
        has_molecule = bool(all_smiles)

        # Classify task type
        task_type = self._classify_type(
            q_lower=q_lower,
            has_reaction=has_reaction,
            reaction_smiles=reaction_smiles,
            molecule_smiles=molecule_smiles,
        )

        # Primary SMILES: prefer reaction over molecule
        primary = reaction_smiles[0] if reaction_smiles else (molecule_smiles[0] if molecule_smiles else None)

        return ClassificationResult(
            task_type=task_type,
            task_type_str=task_type.value,
            has_reaction=has_reaction,
            has_molecule=has_molecule,
            reaction_smiles=reaction_smiles,
            molecule_smiles=molecule_smiles,
            all_smiles=all_smiles,
            primary_smiles=primary,
        )

    def _classify_type(
        self,
        q_lower: str,
        has_reaction: bool,
        reaction_smiles: List[str],
        molecule_smiles: List[str],
    ) -> TaskType:
        """Match keywords in priority order."""
        forward_like = self._is_forward_synthesis_query(
            q_lower=q_lower,
            has_reaction=has_reaction,
            reaction_smiles=reaction_smiles,
            molecule_smiles=molecule_smiles,
        )
        for task_type, keywords in _TASK_KEYWORDS:
            matched = any(kw in q_lower for kw in keywords)
            if matched:
                # Retrosynthesis requires a TARGET molecule, not a full reaction SMILES.
                # If the user already provided a reaction (>>), default to ANALYZE instead.
                if task_type == TaskType.RETROSYNTHESIS and has_reaction:
                    return TaskType.ANALYZE
                # Forward synthesis product prediction is chemistry-distinct from
                # generic PREDICT/ANALYZE and should route to its specialized workflow.
                if task_type in (TaskType.PREDICT, TaskType.ANALYZE) and forward_like:
                    return TaskType.FORWARD_SYNTHESIS
                return task_type

        # Default heuristics based on SMILES presence
        if forward_like:
            return TaskType.FORWARD_SYNTHESIS
        if has_reaction:
            return TaskType.ANALYZE
        if "?" in q_lower:
            return TaskType.GENERAL
        return TaskType.GENERAL

    def _is_forward_synthesis_query(
        self,
        q_lower: str,
        has_reaction: bool,
        reaction_smiles: List[str],
        molecule_smiles: List[str],
    ) -> bool:
        """
        Detect forward-synthesis/product-prediction queries using chemistry signals.

        Baseline policy:
        - never route full reaction SMILES (`>>`) to forward_synthesis
        - prefer forward_synthesis when user asks for product prediction and provides
          multiple reactants (SMILES or explicit reactant/reagent cues)
        """
        if has_reaction or reaction_smiles:
            return False

        explicit_forward = any(cue in q_lower for cue in _FORWARD_EXPLICIT_CUES)
        product_prediction_cue = any(
            cue in q_lower for cue in ("what will happen", "what happens", "forms when", "gives what")
        )
        multi_reactant_signal = (
            len(molecule_smiles) >= 2
            or any(cue in q_lower for cue in _FORWARD_REACTANT_CUES)
            or " + " in q_lower
        )

        if explicit_forward and (multi_reactant_signal or len(molecule_smiles) >= 1):
            return True
        if product_prediction_cue and len(molecule_smiles) >= 2:
            return True
        return False

    def _extract_reaction_smiles(self, query: str) -> List[str]:
        """Extract reaction SMILES (containing >>) from query text."""
        matches = _REACTION_SMILES_PATTERN.findall(query)
        return [m for m in matches if ">>" in m]

    def _extract_molecule_smiles(self, query: str, known_reactions: List[str]) -> List[str]:
        """Extract molecule SMILES (not part of a known reaction SMILES)."""
        # Remove known reaction SMILES from query before extracting molecules
        q_clean = query
        for rsmi in known_reactions:
            q_clean = q_clean.replace(rsmi, " ")

        candidates = _SMILES_PATTERN.findall(q_clean)
        results = []
        for candidate in candidates:
            # Must contain at least one organic atom
            if not any(c in candidate for c in ("C", "N", "O", "S", "P", "F", "B", "c", "n", "o", "s")):
                continue
            # Skip obvious non-SMILES tokens (URLs, file paths, etc.)
            if any(x in candidate for x in ("http", "//", "\\", ".py", ".json")):
                continue
            # Must look like SMILES: contain special chars OR ring closures OR atomic brackets
            # Plain English words (Recommend, conditions, Suzuki) are excluded this way
            has_smiles_chars = any(x in candidate for x in ("(", ")", "[", "]", "=", "#", "@", "+"))
            has_ring_closure = bool(re.search(r"[A-Za-z]\d", candidate))  # e.g. c1, C1, N2
            if has_smiles_chars or has_ring_closure:
                results.append(candidate)
        return results
