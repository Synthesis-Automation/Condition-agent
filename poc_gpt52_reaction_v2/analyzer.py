"""General-purpose reaction analyzer PoC (taxonomy-first + optional LLM reranking)."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from chemtools.detection import detect_reaction_type
from chemtools.taxonomy.loader import load_reaction_types_dict


class GeneralReactionAnalysisError(ValueError):
    """Raised when analysis cannot proceed."""


@dataclass
class ReactionCandidate:
    """One taxonomy candidate from deterministic detection."""

    reaction_type: str
    deterministic_score: float
    detector_confidence: float
    taxonomy_name: str
    taxonomy_category: str
    evidence: Dict[str, Any] = field(default_factory=dict)


@dataclass
class FinalDecision:
    """Final selected reaction type."""

    reaction_type: str
    confidence: float
    source: str
    rationale: str


@dataclass
class ValidationGate:
    """Validator output for chemistry and consistency checks."""

    passed: bool
    issues: List[str] = field(default_factory=list)


@dataclass
class GeneralReactionAnalysisResult:
    """Structured output for the general-purpose v2 PoC."""

    reaction_smiles: str
    reactant_smiles_list: List[str]
    product_smiles_list: List[str]
    canonical_reactants: List[str]
    canonical_products: List[str]
    principal_pair: Dict[str, Any]
    mcs_smarts: str
    mcs_atoms: int
    mcs_ratio: float
    core_reactant_formula: str
    core_product_formula: str
    core_formula_delta: Dict[str, int]
    side_reactant_formula: str
    side_product_formula: str
    side_formula_delta: Dict[str, int]
    detection_error: Optional[str]
    reacted_motifs: List[str]
    formed_motifs: List[str]
    reaction_key: str
    taxonomy_candidates: List[ReactionCandidate]
    decision: FinalDecision
    validation: ValidationGate
    llm_rerank: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-ready dict."""
        return asdict(self)


def _split_side(side: str) -> List[str]:
    parts = [token.strip() for token in side.strip().split(".")]
    return [token for token in parts if token]


def _require_mol(smiles: str, label: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise GeneralReactionAnalysisError(f"Invalid {label} SMILES: {smiles}")
    return mol


def _formula_to_counts(formula: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[symbol] = counts.get(symbol, 0) + (int(count) if count else 1)
    return counts


def _count_delta(reactant_formula: str, product_formula: str) -> Dict[str, int]:
    reactant = _formula_to_counts(reactant_formula)
    product = _formula_to_counts(product_formula)
    keys = sorted(set(reactant) | set(product))
    return {key: product.get(key, 0) - reactant.get(key, 0) for key in keys}


def _counts_to_formula(counts: Dict[str, int]) -> str:
    if not counts:
        return ""
    pieces: List[str] = []
    if "C" in counts:
        pieces.append(f"C{counts['C']}" if counts["C"] != 1 else "C")
    if "H" in counts:
        pieces.append(f"H{counts['H']}" if counts["H"] != 1 else "H")
    for symbol in sorted(k for k in counts if k not in {"C", "H"}):
        number = counts[symbol]
        pieces.append(f"{symbol}{number}" if number != 1 else symbol)
    return "".join(pieces)


def _merge_formula_counts(formulas: Sequence[str]) -> Dict[str, int]:
    merged: Dict[str, int] = {}
    for formula in formulas:
        for symbol, count in _formula_to_counts(formula).items():
            merged[symbol] = merged.get(symbol, 0) + count
    return merged


def _select_principal_pair(
    reactant_mols: Sequence[Chem.Mol],
    product_mols: Sequence[Chem.Mol],
) -> Tuple[int, int, str, int, float]:
    best_i = 0
    best_j = 0
    best_smarts = ""
    best_atoms = -1
    best_ratio = -1.0

    for i, reactant_mol in enumerate(reactant_mols):
        for j, product_mol in enumerate(product_mols):
            mcs = rdFMCS.FindMCS([reactant_mol, product_mol], timeout=4, ringMatchesRingOnly=True)
            atoms = int(getattr(mcs, "numAtoms", 0))
            denominator = max(reactant_mol.GetNumAtoms(), product_mol.GetNumAtoms())
            ratio = (atoms / denominator) if denominator else 0.0
            if atoms > best_atoms or (atoms == best_atoms and ratio > best_ratio):
                best_i = i
                best_j = j
                best_smarts = getattr(mcs, "smartsString", "")
                best_atoms = atoms
                best_ratio = ratio

    return best_i, best_j, best_smarts, best_atoms, float(best_ratio)


def _strip_code_fence(text: str) -> str:
    content = text.strip()
    if content.startswith("```") and content.endswith("```"):
        lines = content.splitlines()
        if len(lines) >= 3:
            return "\n".join(lines[1:-1]).strip()
    return content


def _parse_json_object(text: str) -> Optional[Dict[str, Any]]:
    cleaned = _strip_code_fence(text)
    try:
        payload = json.loads(cleaned)
        return payload if isinstance(payload, dict) else None
    except Exception:
        pass

    left = cleaned.find("{")
    right = cleaned.rfind("}")
    if left >= 0 and right > left:
        snippet = cleaned[left : right + 1]
        try:
            payload = json.loads(snippet)
            return payload if isinstance(payload, dict) else None
        except Exception:
            return None
    return None


class GeneralReactionAnalyzerV2:
    """General-purpose reaction analyzer with taxonomy-grounded reasoning."""

    def __init__(
        self,
        *,
        use_llm: bool = False,
        provider: str = "openai",
        model: str = "gpt-5.2",
        max_candidates: int = 8,
        min_confidence: float = 0.5,
    ):
        self.use_llm = use_llm
        self.provider = provider
        self.model = model
        self.max_candidates = max_candidates
        self.min_confidence = min_confidence

    def analyze(self, reaction_smiles: str) -> GeneralReactionAnalysisResult:
        """Analyze a reaction in a general-purpose, taxonomy-first way."""
        if ">>" not in reaction_smiles:
            raise GeneralReactionAnalysisError("Reaction SMILES must contain '>>'.")

        reactants_side, products_side = reaction_smiles.split(">>", 1)
        reactant_smiles_list = _split_side(reactants_side)
        product_smiles_list = _split_side(products_side)
        if not reactant_smiles_list or not product_smiles_list:
            raise GeneralReactionAnalysisError("Reaction SMILES must contain at least one reactant and one product.")

        reactant_mols = [_require_mol(smiles, "reactant") for smiles in reactant_smiles_list]
        product_mols = [_require_mol(smiles, "product") for smiles in product_smiles_list]

        canonical_reactants = [Chem.MolToSmiles(mol) for mol in reactant_mols]
        canonical_products = [Chem.MolToSmiles(mol) for mol in product_mols]

        pair_i, pair_j, mcs_smarts, mcs_atoms, mcs_ratio = _select_principal_pair(reactant_mols, product_mols)
        core_reactant_mol = reactant_mols[pair_i]
        core_product_mol = product_mols[pair_j]

        core_reactant_formula = CalcMolFormula(core_reactant_mol)
        core_product_formula = CalcMolFormula(core_product_mol)
        core_formula_delta = _count_delta(core_reactant_formula, core_product_formula)

        side_reactant_formulas = [CalcMolFormula(mol) for mol in reactant_mols]
        side_product_formulas = [CalcMolFormula(mol) for mol in product_mols]
        side_reactant_formula = _counts_to_formula(_merge_formula_counts(side_reactant_formulas))
        side_product_formula = _counts_to_formula(_merge_formula_counts(side_product_formulas))
        side_formula_delta = _count_delta(side_reactant_formula, side_product_formula)

        principal_pair = {
            "reactant_index": pair_i,
            "product_index": pair_j,
            "reactant_smiles": reactant_smiles_list[pair_i],
            "product_smiles": product_smiles_list[pair_j],
            "canonical_reactant": canonical_reactants[pair_i],
            "canonical_product": canonical_products[pair_j],
        }

        detection_result = detect_reaction_type(reaction_smiles)
        reaction_types = load_reaction_types_dict()
        reaction_types_lc = {key.lower(): value for key, value in reaction_types.items()}

        candidates = self._build_candidates(
            matches=detection_result.matches,
            reaction_types=reaction_types,
            reaction_types_lc=reaction_types_lc,
            mcs_ratio=mcs_ratio,
        )
        deterministic_decision = self._deterministic_decision(candidates, detection_result.error)
        llm_rerank = None
        decision = deterministic_decision

        if self.use_llm and candidates:
            llm_rerank = self._llm_rerank(
                reaction_smiles=reaction_smiles,
                candidates=candidates,
                deterministic_decision=deterministic_decision,
                principal_pair=principal_pair,
                core_formula_delta=core_formula_delta,
                side_formula_delta=side_formula_delta,
                reacted_motifs=detection_result.reacted_motifs,
                formed_motifs=detection_result.formed_motifs,
            )
            decision = self._merge_llm_decision(
                deterministic=deterministic_decision,
                llm_rerank=llm_rerank,
                allowed=[candidate.reaction_type for candidate in candidates],
            )

        validation = self._validate_decision(
            decision=decision,
            candidates=candidates,
            detection_error=detection_result.error,
        )
        if not validation.passed and decision.reaction_type != "unknown":
            decision = FinalDecision(
                reaction_type="unknown",
                confidence=0.0,
                source="validator_fallback",
                rationale="Validator rejected non-unknown classification: " + "; ".join(validation.issues),
            )

        return GeneralReactionAnalysisResult(
            reaction_smiles=reaction_smiles,
            reactant_smiles_list=reactant_smiles_list,
            product_smiles_list=product_smiles_list,
            canonical_reactants=canonical_reactants,
            canonical_products=canonical_products,
            principal_pair=principal_pair,
            mcs_smarts=mcs_smarts,
            mcs_atoms=mcs_atoms,
            mcs_ratio=round(mcs_ratio, 3),
            core_reactant_formula=core_reactant_formula,
            core_product_formula=core_product_formula,
            core_formula_delta=core_formula_delta,
            side_reactant_formula=side_reactant_formula,
            side_product_formula=side_product_formula,
            side_formula_delta=side_formula_delta,
            detection_error=detection_result.error,
            reacted_motifs=detection_result.reacted_motifs,
            formed_motifs=detection_result.formed_motifs,
            reaction_key=detection_result.reaction_key,
            taxonomy_candidates=candidates,
            decision=decision,
            validation=validation,
            llm_rerank=llm_rerank,
        )

    def _build_candidates(
        self,
        *,
        matches: Sequence[Any],
        reaction_types: Dict[str, Dict[str, Any]],
        reaction_types_lc: Dict[str, Dict[str, Any]],
        mcs_ratio: float,
    ) -> List[ReactionCandidate]:
        candidates: List[ReactionCandidate] = []
        for match in list(matches)[: self.max_candidates]:
            reaction_id = str(getattr(match, "reaction_type", "")).strip()
            if not reaction_id:
                continue
            reaction_def = reaction_types.get(reaction_id) or reaction_types_lc.get(reaction_id.lower()) or {}
            detector_confidence = float(getattr(match, "confidence", 0.0))
            cross_slot_bonus = 0.03 if getattr(match, "electrophile", []) and getattr(match, "nucleophile", []) else 0.0
            scaffold_bonus = 0.02 if mcs_ratio >= 0.5 else 0.0
            deterministic_score = min(1.0, detector_confidence + cross_slot_bonus + scaffold_bonus)

            candidates.append(
                ReactionCandidate(
                    reaction_type=reaction_id,
                    deterministic_score=round(deterministic_score, 2),
                    detector_confidence=round(detector_confidence, 2),
                    taxonomy_name=str(reaction_def.get("name") or reaction_id),
                    taxonomy_category=str(reaction_def.get("category") or "unknown"),
                    evidence={
                        "electrophile": list(getattr(match, "electrophile", []) or []),
                        "nucleophile": list(getattr(match, "nucleophile", []) or []),
                        "product": list(getattr(match, "product", []) or []),
                        "slot_sources": dict(getattr(match, "slot_sources", {}) or {}),
                    },
                )
            )

        candidates.sort(
            key=lambda item: (
                -item.deterministic_score,
                -item.detector_confidence,
                item.reaction_type,
            )
        )
        return candidates

    def _deterministic_decision(
        self,
        candidates: Sequence[ReactionCandidate],
        detection_error: Optional[str],
    ) -> FinalDecision:
        if not candidates:
            reason = "No taxonomy candidates from deterministic detector."
            if detection_error:
                reason += f" detector_error={detection_error}"
            return FinalDecision(
                reaction_type="unknown",
                confidence=0.0,
                source="deterministic",
                rationale=reason,
            )

        top = candidates[0]
        if top.deterministic_score < self.min_confidence:
            return FinalDecision(
                reaction_type="unknown",
                confidence=top.deterministic_score,
                source="deterministic",
                rationale=(
                    f"Top candidate {top.reaction_type} is below min confidence "
                    f"({top.deterministic_score:.2f} < {self.min_confidence:.2f})."
                ),
            )

        return FinalDecision(
            reaction_type=top.reaction_type,
            confidence=top.deterministic_score,
            source="deterministic",
            rationale=f"Selected top taxonomy candidate {top.reaction_type} from deterministic evidence.",
        )

    def _llm_rerank(
        self,
        *,
        reaction_smiles: str,
        candidates: Sequence[ReactionCandidate],
        deterministic_decision: FinalDecision,
        principal_pair: Dict[str, Any],
        core_formula_delta: Dict[str, int],
        side_formula_delta: Dict[str, int],
        reacted_motifs: Sequence[str],
        formed_motifs: Sequence[str],
    ) -> Dict[str, Any]:
        try:
            from llmtools.clients import LLMClient
        except Exception as exc:  # pragma: no cover
            return {"enabled": False, "error": f"LLM client import failed: {exc}"}

        allowed = [candidate.reaction_type for candidate in candidates]
        prompt = (
            "You are a chemistry classifier reranker.\n"
            "Select only from the allowed taxonomy reaction_type list.\n"
            "If evidence is insufficient, return unknown.\n"
            "Return JSON only with keys: final_reaction_type, confidence, rationale, alternatives.\n\n"
            f"Allowed reaction_type values: {json.dumps(allowed, ensure_ascii=True)}\n"
            f"Deterministic top decision: {json.dumps(asdict(deterministic_decision), ensure_ascii=True)}\n"
            f"Reaction: {reaction_smiles}\n"
            f"Principal pair: {json.dumps(principal_pair, ensure_ascii=True)}\n"
            f"Core formula delta: {json.dumps(core_formula_delta, ensure_ascii=True)}\n"
            f"Side formula delta: {json.dumps(side_formula_delta, ensure_ascii=True)}\n"
            f"Reacted motifs: {json.dumps(list(reacted_motifs), ensure_ascii=True)}\n"
            f"Formed motifs: {json.dumps(list(formed_motifs), ensure_ascii=True)}\n"
            f"Candidates: {json.dumps([asdict(candidate) for candidate in candidates], ensure_ascii=True)}\n"
        )

        try:
            client = LLMClient(provider=self.provider, model=self.model, max_tokens=1200)
            response = client.chat(prompt=prompt, temperature=0.2)
            payload = _parse_json_object(response.content)
            return {
                "enabled": True,
                "provider": response.provider,
                "model": response.model,
                "tokens": response.total_tokens,
                "raw_content": response.content,
                "parsed": payload,
            }
        except Exception as exc:  # pragma: no cover
            return {"enabled": True, "error": str(exc)}

    def _merge_llm_decision(
        self,
        *,
        deterministic: FinalDecision,
        llm_rerank: Dict[str, Any],
        allowed: Sequence[str],
    ) -> FinalDecision:
        parsed = llm_rerank.get("parsed")
        if not isinstance(parsed, dict):
            return deterministic

        llm_type = str(parsed.get("final_reaction_type") or "").strip()
        if not llm_type:
            return deterministic

        try:
            llm_confidence = max(0.0, min(1.0, float(parsed.get("confidence", 0.0))))
        except Exception:
            llm_confidence = 0.0
        rationale = str(parsed.get("rationale") or "LLM reranked candidate list.")

        if llm_type.lower() == "unknown":
            if llm_confidence >= 0.8:
                return FinalDecision(
                    reaction_type="unknown",
                    confidence=llm_confidence,
                    source="llm_rerank",
                    rationale=rationale,
                )
            return deterministic

        if llm_type not in allowed:
            return deterministic

        if llm_confidence >= max(self.min_confidence, deterministic.confidence + 0.05):
            return FinalDecision(
                reaction_type=llm_type,
                confidence=llm_confidence,
                source="llm_rerank",
                rationale=rationale,
            )
        return deterministic

    def _validate_decision(
        self,
        *,
        decision: FinalDecision,
        candidates: Sequence[ReactionCandidate],
        detection_error: Optional[str],
    ) -> ValidationGate:
        issues: List[str] = []
        candidate_ids = {candidate.reaction_type for candidate in candidates}

        if decision.reaction_type != "unknown" and not candidates:
            issues.append("no_candidates_for_non_unknown_decision")
        if decision.reaction_type != "unknown" and decision.reaction_type not in candidate_ids:
            issues.append("decision_not_in_candidate_list")
        if decision.reaction_type != "unknown" and decision.confidence < self.min_confidence:
            issues.append("decision_confidence_below_threshold")
        if detection_error == "no_motif_changes" and decision.reaction_type != "unknown":
            issues.append("detector_reported_no_motif_changes")

        return ValidationGate(passed=not issues, issues=issues)


def analyze_reaction_general(
    reaction_smiles: str,
    *,
    use_llm: bool = False,
    provider: str = "openai",
    model: str = "gpt-5.2",
    max_candidates: int = 8,
    min_confidence: float = 0.5,
) -> GeneralReactionAnalysisResult:
    """Convenience function for one-shot analysis."""
    analyzer = GeneralReactionAnalyzerV2(
        use_llm=use_llm,
        provider=provider,
        model=model,
        max_candidates=max_candidates,
        min_confidence=min_confidence,
    )
    return analyzer.analyze(reaction_smiles)
