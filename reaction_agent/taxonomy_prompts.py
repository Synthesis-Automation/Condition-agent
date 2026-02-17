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
from pathlib import Path
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

    # ------------------------------------------------------------------
    # Search methods (used by reasoning agent tools)
    # ------------------------------------------------------------------

    def search_reaction_types(
        self,
        keywords: Optional[List[str]] = None,
        bond_type: Optional[str] = None,
        catalyst: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Search reaction types by keywords, bond type, or catalyst.

        Returns matching entries sorted by relevance score.
        """
        rxn_list = _load_reaction_types_list_full()
        if not rxn_list:
            return []

        kw_lower = [kw.lower() for kw in (keywords or [])]
        bond_lower = bond_type.lower().strip() if bond_type else ""
        cat_lower = catalyst.lower().strip() if catalyst else ""

        matches = []
        for entry in rxn_list:
            if not isinstance(entry, dict):
                continue
            score = 0

            # Build searchable text
            searchable_parts = [
                entry.get("id", ""),
                entry.get("description", ""),
            ]
            searchable_parts.extend(entry.get("aliases", []))
            searchable = " ".join(searchable_parts).lower()

            # Keyword matching
            for kw in kw_lower:
                if kw in searchable:
                    score += 1

            # Bond type matching (constraints.include_bond_formed)
            if bond_lower:
                constraints = entry.get("constraints", {})
                bond_formed = constraints.get("include_bond_formed", [])
                if isinstance(bond_formed, list):
                    for bf in bond_formed:
                        if bond_lower in str(bf).lower():
                            score += 2
                            break

            # Catalyst matching
            if cat_lower:
                catalysts = entry.get("catalysts", [])
                if isinstance(catalysts, list):
                    for cat in catalysts:
                        if cat_lower in str(cat).lower():
                            score += 2
                            break

            if score > 0:
                matches.append({**entry, "_match_score": score})

        matches.sort(key=lambda x: x.get("_match_score", 0), reverse=True)
        return matches

    def search_motifs(
        self,
        scaffold: str = "",
        substituent: str = "",
    ) -> List[str]:
        """Search motif labels by scaffold and/or substituent part.

        Returns matching motif label strings.
        """
        scaffold_lower = scaffold.lower().strip()
        sub_lower = substituent.lower().strip().lstrip("-")

        matches = []
        for label in sorted(self.motif_labels):
            # Label format: "Scaffold-Substituent" e.g. "Ar-Br"
            if scaffold_lower and scaffold_lower not in label.lower().split("-")[0]:
                continue
            if sub_lower:
                # Check if substituent appears in the part after first "-"
                parts = label.split("-", 1)
                if len(parts) < 2 or sub_lower not in parts[1].lower():
                    continue
            matches.append(label)
        return matches

    def get_motif_hierarchy(self, broad_motif: str) -> Dict[str, Any]:
        """Get specific sub-motifs for a broad motif from the scope index.

        Returns dict with 'specific_motifs' list and 'scaffold_parents' map.
        """
        scope_data = _load_motif_scope_index()
        scope_map = scope_data.get("scope_map", {})
        parent_map = scope_data.get("scaffold_parent_map", {})

        result: Dict[str, Any] = {"broad_motif": broad_motif}

        # Direct lookup
        if broad_motif in scope_map:
            result["specific_motifs"] = scope_map[broad_motif]
        else:
            # Case-insensitive fallback
            for key, val in scope_map.items():
                if key.lower() == broad_motif.lower():
                    result["specific_motifs"] = val
                    result["canonical_key"] = key
                    break
            else:
                result["specific_motifs"] = []

        # Relevant parent relationships
        relevant_parents = {}
        base_scaffold = broad_motif.split("-")[0]
        for child, parents in parent_map.items():
            if base_scaffold in parents or child == base_scaffold:
                relevant_parents[child] = parents
        if relevant_parents:
            result["scaffold_parents"] = relevant_parents

        return result


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


@lru_cache(maxsize=1)
def _load_reaction_types_list_full() -> List[Dict[str, Any]]:
    """Load full reaction types entries (not just IDs) for search."""
    try:
        from chemtools.taxonomy.loader import load_reaction_types_list
        return load_reaction_types_list()
    except Exception as e:
        logger.warning(f"Could not load reaction types list: {e}")
        return []


@lru_cache(maxsize=1)
def _load_motif_scope_index() -> Dict[str, Any]:
    """Load the motif scope index (broad -> specific motif mapping)."""
    try:
        # Try multiple paths to find the file
        candidates = [
            Path(__file__).resolve().parent.parent / "chemtools" / "taxonomy" / "data" / "motif_scope_index.v1.json",
        ]
        # Also try via chemtools package location
        try:
            import chemtools.taxonomy
            pkg_dir = Path(chemtools.taxonomy.__file__).parent
            candidates.insert(0, pkg_dir / "data" / "motif_scope_index.v1.json")
        except Exception:
            pass

        for path in candidates:
            if path.exists():
                with open(path, encoding="utf-8") as f:
                    return json.load(f)

        logger.warning("motif_scope_index.v1.json not found")
        return {}
    except Exception as e:
        logger.warning(f"Could not load motif scope index: {e}")
        return {}


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
    # Cached loaders (used by reasoning_tools.py)
    "_load_reaction_types_list_full",
    "_load_motif_scope_index",
]
