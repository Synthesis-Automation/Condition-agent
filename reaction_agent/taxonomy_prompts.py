"""
Taxonomy-aware LLM prompts for the reaction pipeline fallback stage.

Loads the controlled vocabulary (reaction types + motif labels) from the
chemtools taxonomy files and builds prompts that constrain the LLM to only
produce identifiers that exist in the taxonomy.

Stage 4 of ReactionPipeline — only invoked when the quality gate fails.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Dict, List, Optional, Set, Tuple

import logging

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# TaxonomyContext — loads and formats controlled vocabulary
# ---------------------------------------------------------------------------

class TaxonomyContext:
    """
    Loads the chemtools taxonomy vocabulary and formats it for LLM prompts.

    Uses lru_cache internally so the JSON files are read only once per process.

    Usage:
        ctx = TaxonomyContext()
        block = ctx.as_prompt_block()    # inject into LLM prompt
        valid = ctx.validate_reaction_type("Suzuki_miyaura")  # True
        valid = ctx.validate_motif("Ar-Br")                   # True
    """

    def __init__(self) -> None:
        self._reaction_type_ids: Optional[Set[str]] = None
        self._motif_labels: Optional[Set[str]] = None

    # ------------------------------------------------------------------
    # Public helpers
    # ------------------------------------------------------------------

    @property
    def reaction_type_ids(self) -> Set[str]:
        if self._reaction_type_ids is None:
            self._reaction_type_ids = _load_reaction_type_ids()
        return self._reaction_type_ids

    @property
    def motif_labels(self) -> Set[str]:
        if self._motif_labels is None:
            self._motif_labels = _load_motif_labels()
        return self._motif_labels

    def validate_reaction_type(self, reaction_type: Optional[str]) -> bool:
        """Return True if reaction_type is in the taxonomy vocabulary."""
        return bool(reaction_type and reaction_type in self.reaction_type_ids)

    def validate_motif(self, motif: Optional[str]) -> bool:
        """Return True if motif label is in the taxonomy vocabulary."""
        return bool(motif and motif in self.motif_labels)

    def filter_motifs(self, motifs: List[str]) -> Tuple[List[str], List[str]]:
        """
        Split motif list into (valid, invalid) based on taxonomy vocabulary.

        Args:
            motifs: List of motif label strings from LLM output

        Returns:
            (valid_motifs, invalid_motifs) — only valid ones should be used
        """
        valid, invalid = [], []
        for m in motifs:
            (valid if self.validate_motif(m) else invalid).append(m)
        return valid, invalid

    def as_prompt_block(self) -> str:
        """
        Format the taxonomy vocabulary as a compact LLM-injectable block.

        Returns a string ready to be embedded in an LLM prompt, listing all
        valid reaction type IDs and motif labels in a concise format.
        """
        rxn_ids = sorted(self.reaction_type_ids)
        motif_lbls = sorted(self.motif_labels)

        lines = [
            "=" * 70,
            "TAXONOMY CONTROLLED VOCABULARY — USE ONLY IDENTIFIERS FROM THESE LISTS",
            "=" * 70,
            "",
            "VALID REACTION TYPE IDs (pick exactly one):",
            "  " + ", ".join(rxn_ids),
            "",
            "VALID MOTIF LABELS (scaffold-substituent pairs, e.g. Ar-Br, Alkyl-NH2):",
            "  " + ", ".join(motif_lbls),
            "",
            "IMPORTANT: Any identifier NOT in these lists will be rejected.",
            "=" * 70,
        ]
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Cached loaders
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _load_reaction_type_ids() -> Set[str]:
    """Load all reaction type IDs from chemtools taxonomy."""
    try:
        from chemtools.taxonomy.loader import load_reaction_types_list
        rxn_list = load_reaction_types_list()
        return {entry["id"] for entry in rxn_list if isinstance(entry, dict) and "id" in entry}
    except Exception as e:
        logger.warning(f"TaxonomyContext: could not load reaction types: {e}")
        return set()


@lru_cache(maxsize=1)
def _load_motif_labels() -> Set[str]:
    """
    Load all motif labels (A-B format) from organic_compounds taxonomy.

    The organic_compounds JSON has entries like {"A": "Ar", "B": "-Br"}.
    The motif label is A + B, e.g. "Ar" + "-Br" = "Ar-Br".
    """
    try:
        from chemtools.taxonomy.loader import load_organic_compounds
        data = load_organic_compounds()
        compounds = data.get("compounds", [])
        labels: Set[str] = set()
        for entry in compounds:
            if not isinstance(entry, dict):
                continue
            a = entry.get("A", "")
            b = entry.get("B", "")
            if a and b:
                labels.add(f"{a}{b}")
        return labels
    except Exception as e:
        logger.warning(f"TaxonomyContext: could not load organic compounds: {e}")
        return set()


# ---------------------------------------------------------------------------
# LLM fallback prompt builder
# ---------------------------------------------------------------------------

def build_llm_fallback_prompt(
    normalized_smiles: str,
    taxonomy_block: str,
    feat_reaction_type: Optional[str],
    feat_confidence: float,
    feat_reacted_motifs: Tuple[str, ...],
    feat_formed_motifs: Tuple[str, ...],
    quality_issues: List[str],
) -> str:
    """
    Build the taxonomy-aware LLM prompt for Stage 4 fallback.

    Args:
        normalized_smiles: Canonical reaction SMILES (A.B>>P)
        taxonomy_block: Output of TaxonomyContext.as_prompt_block()
        feat_reaction_type: Reaction type from featurize_reaction() (may be None/uncertain)
        feat_confidence: Confidence from featurize_reaction()
        feat_reacted_motifs: Reacted motifs from featurize_reaction() (may be empty/wrong)
        feat_formed_motifs: Formed motifs from featurize_reaction() (may be empty/wrong)
        quality_issues: List of quality gate failure reasons

    Returns:
        Complete prompt string to send to LLM
    """
    issues_str = "\n".join(f"  - {issue}" for issue in quality_issues) if quality_issues else "  (none)"
    feat_rt = feat_reaction_type or "Unknown"
    feat_rm = list(feat_reacted_motifs) or []
    feat_fm = list(feat_formed_motifs) or []

    prompt = f"""You are a computational chemistry expert assigning taxonomy identifiers to reactions.

{taxonomy_block}

REACTION SMILES:
  {normalized_smiles}

DETERMINISTIC ANALYSIS (may be incomplete or low-confidence):
  reaction_type     : {feat_rt}  (confidence: {feat_confidence:.2f})
  reacted_motifs    : {feat_rm}
  formed_motifs     : {feat_fm}

QUALITY ISSUES DETECTED (why this fallback was triggered):
{issues_str}

YOUR TASK:
Analyze the reaction SMILES and assign the correct taxonomy identifiers.

Rules:
1. reaction_type MUST be exactly one ID from the VALID REACTION TYPE IDs list above.
2. reacted_motifs MUST be exactly 2 items, each from the VALID MOTIF LABELS list.
   - One motif per reactant; order does not matter.
   - These are the functional groups CONSUMED by the reaction.
3. formed_motifs MUST contain 1–3 items from the VALID MOTIF LABELS list.
   - These are the functional groups CREATED in the product.
4. If the deterministic analysis looks mostly correct, use it as a starting point.
5. confidence should reflect how certain you are (0.0–1.0).

Return ONLY a JSON object with no markdown fences:
{{
  "reaction_type": "<taxonomy ID>",
  "reacted_motifs": ["<motif1>", "<motif2>"],
  "formed_motifs": ["<motif1>"],
  "confidence": 0.0,
  "reasoning": "<one or two sentences explaining your assignment>"
}}"""

    return prompt


# ---------------------------------------------------------------------------
# LLM response parser
# ---------------------------------------------------------------------------

def parse_llm_fallback_response(
    raw_response: str,
    taxonomy: TaxonomyContext,
) -> Dict[str, Any]:
    """
    Parse and validate the JSON response from the LLM fallback.

    Validates all returned identifiers against the taxonomy vocabulary.
    Removes any invalid identifiers and records them as warnings.

    Args:
        raw_response: Raw LLM response string (expected to be JSON)
        taxonomy: TaxonomyContext for validation

    Returns:
        Dict with keys:
            reaction_type     : Optional[str]  — validated ID or None
            reacted_motifs    : List[str]       — validated motif labels
            formed_motifs     : List[str]       — validated motif labels
            confidence        : float
            reasoning         : str
            invalid_ids_found : List[str]       — rejected identifiers
            warnings          : List[str]
    """
    warnings: List[str] = []
    invalid_ids: List[str] = []

    # Strip markdown fences if present
    text = raw_response.strip()
    text = re.sub(r'^```(?:json)?\s*\n', '', text, flags=re.MULTILINE)
    text = re.sub(r'\n```\s*$', '', text, flags=re.MULTILINE)
    text = text.strip()

    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        warnings.append(f"LLM response is not valid JSON: {e}")
        return {
            "reaction_type": None,
            "reacted_motifs": [],
            "formed_motifs": [],
            "confidence": 0.0,
            "reasoning": "",
            "invalid_ids_found": [],
            "warnings": warnings,
        }

    # Validate reaction_type
    reaction_type = data.get("reaction_type")
    if reaction_type and not taxonomy.validate_reaction_type(reaction_type):
        invalid_ids.append(reaction_type)
        warnings.append(f"LLM reaction_type '{reaction_type}' not in taxonomy — discarded")
        reaction_type = None

    # Validate reacted_motifs
    raw_reacted = data.get("reacted_motifs", [])
    if not isinstance(raw_reacted, list):
        raw_reacted = []
    valid_reacted, inv_reacted = taxonomy.filter_motifs(raw_reacted)
    if inv_reacted:
        invalid_ids.extend(inv_reacted)
        warnings.append(f"Invalid reacted_motifs discarded: {inv_reacted}")

    # Validate formed_motifs
    raw_formed = data.get("formed_motifs", [])
    if not isinstance(raw_formed, list):
        raw_formed = []
    valid_formed, inv_formed = taxonomy.filter_motifs(raw_formed)
    if inv_formed:
        invalid_ids.extend(inv_formed)
        warnings.append(f"Invalid formed_motifs discarded: {inv_formed}")

    confidence = float(data.get("confidence", 0.0))
    confidence = max(0.0, min(1.0, confidence))
    reasoning = str(data.get("reasoning", ""))

    return {
        "reaction_type": reaction_type,
        "reacted_motifs": valid_reacted,
        "formed_motifs": valid_formed,
        "confidence": confidence,
        "reasoning": reasoning,
        "invalid_ids_found": invalid_ids,
        "warnings": warnings,
    }


__all__ = [
    "TaxonomyContext",
    "build_llm_fallback_prompt",
    "parse_llm_fallback_response",
]
