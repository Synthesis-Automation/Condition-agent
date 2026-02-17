"""
Modular 5-stage reaction SMILES pipeline.

Orchestrates SMILES normalization, deterministic featurization, quality
evaluation, LLM-taxonomy fallback, and result merging into a single coherent
workflow.  Each stage is a public method so it can be tested and debugged
independently.

Stages:
    1. normalize()   — SMILES pre-check and canonicalization (chemtools)
    2. featurize()   — Deterministic taxonomy featurization (featurize_reaction)
    3. evaluate()    — Quality gate (4 hard criteria)
    4. llm_fallback()— LLM-assisted recovery constrained to taxonomy vocabulary
    5. merge()       — Patch featurize_reaction() fields with LLM where needed

Usage:
    from llmtools.clients import LLMClient
    from reaction_agent.smiles_pipeline import ReactionPipeline
    from reaction_agent.pipeline_eval import QualityConfig

    client = LLMClient(provider="openai", model="gpt-4o")
    pipeline = ReactionPipeline(llm_client=client)

    result = pipeline.run("Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
    print(result.reaction_type)          # "Suzuki_miyaura"
    print(result.reacted_motifs)         # ("Ar-Br", "Ar-B(OH)2")
    print(result.used_llm_fallback)      # False  (quality gate passed)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from llmtools.clients import LLMClient

from reaction_agent.pipeline_eval import QualityConfig, QualityEvaluator, QualityReport
from reaction_agent.taxonomy_prompts import (
    TaxonomyContext,
    build_llm_fallback_prompt,
    parse_llm_fallback_response,
)

logger = logging.getLogger(__name__)

_UNCLASSIFIED_MARKER = "Unclassified-Reactant"


# ===========================================================================
# Stage dataclasses — each stage's full output is preserved in PipelineResult
# ===========================================================================

@dataclass
class NormalizationResult:
    """Output of Stage 1 (SMILES normalization and pre-check)."""

    success: bool
    normalized_smiles: str           # "A.B>>P" in canonical form
    reactants: List[str]             # Individual reactant SMILES (canonical)
    product: Optional[str]           # First product SMILES (canonical), or None
    warnings: List[str] = field(default_factory=list)
    error: Optional[str] = None      # Human-readable error if success=False


@dataclass
class FeaturizationResult:
    """Output of Stage 2 (deterministic featurize_reaction)."""

    success: bool
    reaction_key: Optional[str]                    # CRK-v1 key, e.g. "bond_formed: [C-C] | ..."
    reaction_type: Optional[str]                   # Taxonomy reaction type ID
    reaction_type_confidence: float                # 0.0–1.0
    reacted_motifs: Tuple[str, ...]                # Motifs consumed (one per reactant)
    formed_motifs: Tuple[str, ...]                 # Motifs created in product
    spectator_motifs: Tuple[str, ...]              # Unchanged motifs
    has_unclassified_reactant: bool                # True if coverage guard fired
    reactant_motif_count: int                      # len(reacted_motifs) excl. unclassified
    raw: Dict[str, Any] = field(default_factory=dict)   # Full featurize_reaction() output
    warnings: List[str] = field(default_factory=list)


@dataclass
class LLMFallbackResult:
    """Output of Stage 4 (LLM taxonomy fallback)."""

    success: bool
    reaction_type: Optional[str]        # Taxonomy-validated reaction type ID
    reacted_motifs: Tuple[str, ...]     # Taxonomy-validated motif labels
    formed_motifs: Tuple[str, ...]      # Taxonomy-validated motif labels
    confidence: float                   # LLM self-reported confidence (0.0–1.0)
    reasoning: str                      # LLM explanation
    invalid_ids_found: List[str] = field(default_factory=list)   # Rejected IDs
    raw_response: str = ""              # Full LLM JSON string for debugging
    warnings: List[str] = field(default_factory=list)


@dataclass
class PipelineResult:
    """
    Final output of the full pipeline.

    All intermediate stage results are preserved for debugging.
    Use the top-level merged fields (reaction_type, reacted_motifs, etc.)
    for downstream consumption.
    """

    # Input
    raw_smiles: str

    # Stage outputs (always present; check success flag before using)
    normalization: NormalizationResult
    featurization: Optional[FeaturizationResult] = None
    quality: Optional[QualityReport] = None
    llm_fallback: Optional[LLMFallbackResult] = None

    # Final merged / patched fields — use these downstream
    reaction_type: Optional[str] = None
    reaction_type_confidence: float = 0.0
    reacted_motifs: Tuple[str, ...] = ()
    formed_motifs: Tuple[str, ...] = ()
    spectator_motifs: Tuple[str, ...] = ()
    reaction_key: Optional[str] = None
    used_llm_fallback: bool = False

    # Aggregated warnings from all stages
    pipeline_warnings: List[str] = field(default_factory=list)

    @property
    def success(self) -> bool:
        """True if normalization succeeded (minimum requirement for any result)."""
        return self.normalization.success

    @classmethod
    def from_failed_normalization(cls, norm: NormalizationResult) -> "PipelineResult":
        """Build a failed result when normalization itself fails."""
        return cls(
            raw_smiles=norm.normalized_smiles or "",
            normalization=norm,
            pipeline_warnings=norm.warnings + ([norm.error] if norm.error else []),
        )


# ===========================================================================
# Pipeline orchestrator
# ===========================================================================

class ReactionPipeline:
    """
    5-stage modular pipeline: normalize → featurize → evaluate → fallback → merge.

    Each stage is a public method and returns a typed dataclass so any stage
    can be called in isolation for unit testing or debugging.

    Args:
        llm_client: LLM client for Stage 4 fallback (optional — if None,
            the pipeline skips Stage 4 and returns featurize_reaction() output
            even if the quality gate fails).
        quality_config: Configures the four quality criteria for Stage 3.
        llm_model_override: Override the LLM model for Stage 4 (e.g. "gpt-4o").
            If None, the client's current model is used.

    Example:
        >>> from llmtools.clients import LLMClient
        >>> from reaction_agent.smiles_pipeline import ReactionPipeline
        >>> client = LLMClient(provider="openai", model="gpt-4o")
        >>> pipeline = ReactionPipeline(llm_client=client)
        >>> result = pipeline.run("Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
        >>> print(result.reaction_type)
    """

    def __init__(
        self,
        llm_client: Optional["LLMClient"] = None,
        quality_config: Optional[QualityConfig] = None,
        llm_model_override: Optional[str] = None,
    ) -> None:
        self.llm_client = llm_client
        self._evaluator = QualityEvaluator(quality_config)
        self._taxonomy = TaxonomyContext()
        self._llm_model_override = llm_model_override

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, smiles: str) -> PipelineResult:
        """
        Execute the full pipeline for a single reaction SMILES.

        Args:
            smiles: Raw reaction SMILES string (reactants>>product)

        Returns:
            PipelineResult with all stage outputs and merged final fields
        """
        # Stage 1 — normalize
        norm = self.normalize(smiles)
        if not norm.success:
            logger.warning(f"Pipeline Stage 1 failed: {norm.error}")
            return PipelineResult.from_failed_normalization(norm)

        # Stage 2 — featurize
        feat = self.featurize(norm.normalized_smiles)

        # Stage 3 — evaluate quality
        quality = self.evaluate(feat)

        # Stage 4 — LLM fallback (conditional)
        llm_result: Optional[LLMFallbackResult] = None
        if quality.needs_llm_fallback:
            if self.llm_client is not None:
                logger.info(
                    f"Pipeline Stage 4: LLM fallback triggered. Issues: {quality.issues}"
                )
                llm_result = self.llm_fallback(norm.normalized_smiles, feat, quality)
            else:
                logger.warning(
                    "Pipeline Stage 3 failed but no LLM client provided — "
                    "skipping Stage 4. Quality issues: %s",
                    quality.issues,
                )

        # Stage 5 — merge
        return self.merge(norm, feat, quality, llm_result)

    # ------------------------------------------------------------------
    # Stage 1 — SMILES normalization and pre-check
    # ------------------------------------------------------------------

    def normalize(self, smiles: str) -> NormalizationResult:
        """
        Normalize and validate a reaction SMILES string.

        Uses chemtools.smiles.normalize_reaction() for RDKit canonicalization.
        Runs pre-checks: non-empty, exactly one '>>', at least one reactant,
        product present, no individual SMILES tokens with parse errors.

        Args:
            smiles: Raw reaction SMILES (any format)

        Returns:
            NormalizationResult — check .success before proceeding
        """
        warnings: List[str] = []

        # Basic string checks
        smiles = (smiles or "").strip()
        if not smiles:
            return NormalizationResult(
                success=False,
                normalized_smiles="",
                reactants=[],
                product=None,
                error="Empty SMILES input",
            )

        if smiles.count(">>") != 1:
            return NormalizationResult(
                success=False,
                normalized_smiles=smiles,
                reactants=[],
                product=None,
                error=f"Reaction SMILES must contain exactly one '>>'. Got: '{smiles}'",
            )

        try:
            from chemtools.smiles import normalize_reaction
        except ImportError:
            return NormalizationResult(
                success=False,
                normalized_smiles=smiles,
                reactants=[],
                product=None,
                error="chemtools.smiles not available",
            )

        try:
            record = normalize_reaction(smiles)
        except Exception as e:
            return NormalizationResult(
                success=False,
                normalized_smiles=smiles,
                reactants=[],
                product=None,
                error=f"normalize_reaction() raised: {e}",
            )

        # Check for parse errors in any component
        errors = record.get("errors", [])
        if errors:
            warnings.append(f"Some SMILES tokens failed to parse: {errors}")

        # Extract reactant and product SMILES
        reactant_payloads = record.get("reactants", [])
        product_payloads = record.get("products", [])

        reactants = [
            p.get("smiles_norm") or p.get("largest_smiles") or p.get("input", "")
            for p in reactant_payloads
            if isinstance(p, dict)
        ]
        reactants = [r for r in reactants if r]

        products_raw = [
            p.get("smiles_norm") or p.get("largest_smiles") or p.get("input", "")
            for p in product_payloads
            if isinstance(p, dict)
        ]
        products_raw = [p for p in products_raw if p]
        product = products_raw[0] if products_raw else None

        # Validate: at least one reactant
        if not reactants:
            return NormalizationResult(
                success=False,
                normalized_smiles=smiles,
                reactants=[],
                product=product,
                error="No valid reactant SMILES found after normalization",
                warnings=warnings,
            )

        # Validate: product required for downstream taxonomy matching
        if not product:
            return NormalizationResult(
                success=False,
                normalized_smiles=smiles,
                reactants=reactants,
                product=None,
                error=(
                    "Product SMILES is required for featurize_reaction() and HTE matching. "
                    "Provide reaction SMILES in format: reactants>>product"
                ),
                warnings=warnings,
            )

        # Reconstruct canonical reaction SMILES: "A.B>>P"
        normalized_smiles = ".".join(reactants) + ">>" + product

        return NormalizationResult(
            success=True,
            normalized_smiles=normalized_smiles,
            reactants=reactants,
            product=product,
            warnings=warnings,
        )

    # ------------------------------------------------------------------
    # Stage 2 — Deterministic featurization
    # ------------------------------------------------------------------

    def featurize(self, normalized_smiles: str) -> FeaturizationResult:
        """
        Run featurize_reaction() with CRK options on the normalized SMILES.

        Extracts reaction_key, reaction_type, reacted_motifs, formed_motifs,
        and spectator_motifs from the returned payload.  Also detects whether
        the coverage guard inserted 'Unclassified-Reactant'.

        Args:
            normalized_smiles: Canonical reaction SMILES from Stage 1

        Returns:
            FeaturizationResult — check .success before proceeding
        """
        warnings: List[str] = []

        try:
            from chemtools.featurizers.unified import featurize_reaction
            from chemtools.featurizers.formatters.reaction import get_crk_options
        except ImportError as e:
            return FeaturizationResult(
                success=False,
                reaction_key=None,
                reaction_type=None,
                reaction_type_confidence=0.0,
                reacted_motifs=(),
                formed_motifs=(),
                spectator_motifs=(),
                has_unclassified_reactant=False,
                reactant_motif_count=0,
                warnings=[f"chemtools import failed: {e}"],
            )

        try:
            options = get_crk_options()
            # Also request reaction type detection for the quality gate
            options["include_reaction_type"] = True
            payload = featurize_reaction(normalized_smiles, options=options)
        except Exception as e:
            return FeaturizationResult(
                success=False,
                reaction_key=None,
                reaction_type=None,
                reaction_type_confidence=0.0,
                reacted_motifs=(),
                formed_motifs=(),
                spectator_motifs=(),
                has_unclassified_reactant=False,
                reactant_motif_count=0,
                warnings=[f"featurize_reaction() raised: {e}"],
                raw={},
            )

        if not isinstance(payload, dict):
            return FeaturizationResult(
                success=False,
                reaction_key=None,
                reaction_type=None,
                reaction_type_confidence=0.0,
                reacted_motifs=(),
                formed_motifs=(),
                spectator_motifs=(),
                has_unclassified_reactant=False,
                reactant_motif_count=0,
                warnings=["featurize_reaction() returned non-dict"],
                raw={},
            )

        # Unwrap the "reaction" wrapper if present (legacy compat shim in unified.py)
        reaction = payload.get("reaction") if isinstance(payload.get("reaction"), dict) else payload

        aggregates = reaction.get("aggregates") or {}
        reaction_key = str(reaction.get("reaction_key") or "").strip()

        reacted_motifs = tuple(aggregates.get("reacted_motifs") or [])
        formed_motifs = tuple(aggregates.get("formed_motifs") or [])
        spectator_motifs = tuple(aggregates.get("spectator_motifs") or [])

        # Reaction type detection
        # featurize_reaction() returns reaction_type as either:
        #   str  — plain taxonomy ID, e.g. "Suzuki_miyaura" (high confidence, treat as 1.0)
        #   dict — {"reaction_type": "...", "confidence": 0.92}  (explicit confidence)
        #   None — detection failed / not requested
        rxn_type_data = reaction.get("reaction_type")
        if isinstance(rxn_type_data, str) and rxn_type_data:
            reaction_type = rxn_type_data
            reaction_type_confidence = 1.0   # string form = deterministically assigned
        elif isinstance(rxn_type_data, dict):
            reaction_type = rxn_type_data.get("reaction_type") or rxn_type_data.get("type")
            reaction_type_confidence = float(rxn_type_data.get("confidence", 0.0))
        else:
            reaction_type = None
            reaction_type_confidence = 0.0

        # Detect coverage guard marker and count real reactant motifs
        has_unclassified = _UNCLASSIFIED_MARKER in reacted_motifs
        real_reactant_motifs = [m for m in reacted_motifs if m != _UNCLASSIFIED_MARKER]
        reactant_motif_count = len(real_reactant_motifs)

        return FeaturizationResult(
            success=True,
            reaction_key=reaction_key or None,
            reaction_type=reaction_type,
            reaction_type_confidence=reaction_type_confidence,
            reacted_motifs=reacted_motifs,
            formed_motifs=formed_motifs,
            spectator_motifs=spectator_motifs,
            has_unclassified_reactant=has_unclassified,
            reactant_motif_count=reactant_motif_count,
            raw=payload,
            warnings=warnings,
        )

    # ------------------------------------------------------------------
    # Stage 3 — Quality gate
    # ------------------------------------------------------------------

    def evaluate(self, feat: FeaturizationResult) -> QualityReport:
        """
        Evaluate featurization quality against the four hard criteria.

        Delegates to QualityEvaluator — see pipeline_eval.py for criteria.

        Args:
            feat: Output of Stage 2

        Returns:
            QualityReport with pass/fail and per-criterion detail
        """
        return self._evaluator.evaluate(feat)

    # ------------------------------------------------------------------
    # Stage 4 — LLM taxonomy fallback
    # ------------------------------------------------------------------

    def llm_fallback(
        self,
        normalized_smiles: str,
        feat: FeaturizationResult,
        quality: QualityReport,
    ) -> LLMFallbackResult:
        """
        Call LLM with taxonomy-controlled vocabulary to recover missing/wrong fields.

        The LLM is shown the full taxonomy vocabulary and constrained to only
        produce identifiers that exist in the taxonomy.  Any out-of-vocabulary
        identifiers in the response are removed with a warning.

        Args:
            normalized_smiles: Canonical reaction SMILES from Stage 1
            feat: Stage 2 output (context for LLM prompt)
            quality: Stage 3 report (used to explain why fallback is needed)

        Returns:
            LLMFallbackResult — check .success before merging
        """
        assert self.llm_client is not None, "llm_client required for Stage 4"

        taxonomy_block = self._taxonomy.as_prompt_block()
        prompt = build_llm_fallback_prompt(
            normalized_smiles=normalized_smiles,
            taxonomy_block=taxonomy_block,
            feat_reaction_type=feat.reaction_type,
            feat_confidence=feat.reaction_type_confidence,
            feat_reacted_motifs=feat.reacted_motifs,
            feat_formed_motifs=feat.formed_motifs,
            quality_issues=quality.issues,
        )

        # Use model override if specified
        original_model = self.llm_client.model
        if self._llm_model_override:
            self.llm_client.model = self._llm_model_override

        try:
            response = self.llm_client.chat(
                prompt=prompt,
                temperature=0.0,
                max_tokens=1000,
            )
            raw_response = response.content
        except Exception as e:
            if self._llm_model_override:
                self.llm_client.model = original_model
            return LLMFallbackResult(
                success=False,
                reaction_type=None,
                reacted_motifs=(),
                formed_motifs=(),
                confidence=0.0,
                reasoning="",
                warnings=[f"LLM call failed: {e}"],
                raw_response="",
            )
        finally:
            if self._llm_model_override:
                self.llm_client.model = original_model

        parsed = parse_llm_fallback_response(raw_response, self._taxonomy)

        success = bool(
            parsed.get("reaction_type")
            and len(parsed.get("reacted_motifs", [])) == 2
        )

        if not success:
            parsed["warnings"].append(
                "LLM fallback did not return valid reaction_type + 2 reacted_motifs"
            )

        return LLMFallbackResult(
            success=success,
            reaction_type=parsed.get("reaction_type"),
            reacted_motifs=tuple(parsed.get("reacted_motifs", [])),
            formed_motifs=tuple(parsed.get("formed_motifs", [])),
            confidence=parsed.get("confidence", 0.0),
            reasoning=parsed.get("reasoning", ""),
            invalid_ids_found=parsed.get("invalid_ids_found", []),
            raw_response=raw_response,
            warnings=parsed.get("warnings", []),
        )

    # ------------------------------------------------------------------
    # Stage 5 — Merge / patch
    # ------------------------------------------------------------------

    def merge(
        self,
        norm: NormalizationResult,
        feat: FeaturizationResult,
        quality: QualityReport,
        llm_result: Optional[LLMFallbackResult],
    ) -> PipelineResult:
        """
        Merge Stage 2 (featurize) and Stage 4 (LLM) outputs.

        Strategy:
        - Keep all featurize_reaction() fields that are structural/reliable:
          reaction_key, spectator_motifs, full raw output.
        - Override reaction_type, reacted_motifs, formed_motifs from LLM
          when LLM fallback succeeded AND LLM confidence > featurize confidence.

        Args:
            norm: Stage 1 result
            feat: Stage 2 result
            quality: Stage 3 result
            llm_result: Stage 4 result (None if fallback was skipped)

        Returns:
            PipelineResult with merged fields populated
        """
        all_warnings: List[str] = list(norm.warnings)
        if feat:
            all_warnings.extend(feat.warnings)
        if quality:
            all_warnings.extend(quality.issues)
        if llm_result:
            all_warnings.extend(llm_result.warnings)

        # Start with featurize_reaction() values as the base
        reaction_type = feat.reaction_type if feat.success else None
        reaction_type_confidence = feat.reaction_type_confidence if feat.success else 0.0
        reacted_motifs = feat.reacted_motifs if feat.success else ()
        formed_motifs = feat.formed_motifs if feat.success else ()
        spectator_motifs = feat.spectator_motifs if feat.success else ()
        reaction_key = feat.reaction_key if feat.success else None
        used_llm_fallback = False

        # Patch with LLM output if fallback succeeded
        if llm_result and llm_result.success:
            llm_conf = llm_result.confidence
            feat_conf = feat.reaction_type_confidence if feat.success else 0.0

            # Override reaction_type if LLM is more confident
            if llm_result.reaction_type and llm_conf >= feat_conf:
                reaction_type = llm_result.reaction_type
                reaction_type_confidence = llm_conf
                used_llm_fallback = True

            # Always override motifs when LLM fallback ran (this was the point)
            if llm_result.reacted_motifs:
                reacted_motifs = llm_result.reacted_motifs
                used_llm_fallback = True
            if llm_result.formed_motifs:
                formed_motifs = llm_result.formed_motifs
                used_llm_fallback = True

        elif llm_result and not llm_result.success:
            all_warnings.append(
                "LLM fallback was triggered but did not produce a valid result — "
                "falling back to featurize_reaction() output"
            )

        return PipelineResult(
            raw_smiles=norm.normalized_smiles or "",
            normalization=norm,
            featurization=feat,
            quality=quality,
            llm_fallback=llm_result,
            reaction_type=reaction_type,
            reaction_type_confidence=reaction_type_confidence,
            reacted_motifs=reacted_motifs,
            formed_motifs=formed_motifs,
            spectator_motifs=spectator_motifs,
            reaction_key=reaction_key,
            used_llm_fallback=used_llm_fallback,
            pipeline_warnings=all_warnings,
        )


__all__ = [
    "ReactionPipeline",
    "PipelineResult",
    "NormalizationResult",
    "FeaturizationResult",
    "LLMFallbackResult",
]
