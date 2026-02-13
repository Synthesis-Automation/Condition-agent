"""Proof-of-concept reaction analyzer with deterministic reasoning and optional LLM refinement."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from chemtools.util.smarts_cache import compile_smarts


_SMARTS_PATTERNS = {
    "aryl_halide_like": "[$([c,n][Cl,Br,I,F])]",
    "nn_single": "[N]-[N]",
    "nn_double": "[N]=[N]",
}


class ReactionAnalysisError(ValueError):
    """Raised when reaction analysis cannot proceed."""


@dataclass
class AnalysisStep:
    """One visible reasoning step in the PoC pipeline."""

    name: str
    evidence: Dict[str, Any]
    conclusion: str


@dataclass
class ReactionAnalysisResult:
    """Structured output for deterministic + optional LLM analysis."""

    reaction_smiles: str
    reactant_smiles_list: List[str]
    product_smiles_list: List[str]
    canonical_reactants: List[str]
    canonical_products: List[str]
    reactant_smiles: str
    product_smiles: str
    canonical_reactant: str
    canonical_product: str
    reactant_formula: str
    product_formula: str
    formula_delta: Dict[str, int]
    side_reactant_formula: str
    side_product_formula: str
    side_formula_delta: Dict[str, int]
    selected_pair: Dict[str, Any]
    mcs_smarts: str
    mcs_atoms: int
    features: Dict[str, Any]
    hypotheses: List[Dict[str, Any]]
    best_hypothesis: Dict[str, Any]
    confidence: float
    summary: str
    reasoning_steps: List[AnalysisStep] = field(default_factory=list)
    llm_refinement: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-ready dict."""
        return asdict(self)


def _split_side(side: str) -> List[str]:
    parts = [item.strip() for item in side.strip().split(".")]
    return [item for item in parts if item]


def _formula_to_counts(formula: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[symbol] = counts.get(symbol, 0) + (int(count) if count else 1)
    return counts


def _counts_to_formula(counts: Dict[str, int]) -> str:
    if not counts:
        return ""
    pieces: List[str] = []
    if "C" in counts:
        pieces.append(f"C{counts['C']}" if counts["C"] != 1 else "C")
    if "H" in counts:
        pieces.append(f"H{counts['H']}" if counts["H"] != 1 else "H")
    for symbol in sorted(k for k in counts if k not in {"C", "H"}):
        n = counts[symbol]
        pieces.append(f"{symbol}{n}" if n != 1 else symbol)
    return "".join(pieces)


def _count_delta(reactant_formula: str, product_formula: str) -> Dict[str, int]:
    reactant = _formula_to_counts(reactant_formula)
    product = _formula_to_counts(product_formula)
    all_symbols = sorted(set(reactant) | set(product))
    return {sym: product.get(sym, 0) - reactant.get(sym, 0) for sym in all_symbols}


def _merge_counts(formulas: Sequence[str]) -> Dict[str, int]:
    merged: Dict[str, int] = {}
    for formula in formulas:
        for symbol, count in _formula_to_counts(formula).items():
            merged[symbol] = merged.get(symbol, 0) + count
    return merged


def _require_mol(smiles: str, label: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ReactionAnalysisError(f"Invalid {label} SMILES: {smiles}")
    return mol


def _count_halogens(mol: Chem.Mol) -> int:
    return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() in {"F", "Cl", "Br", "I"})


def _count_aromatic_ring_n(mol: Chem.Mol) -> int:
    return sum(
        1
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "N" and atom.GetIsAromatic() and atom.IsInRing()
    )


def _has_substructure(mol: Chem.Mol, key: str) -> bool:
    pattern = compile_smarts(_SMARTS_PATTERNS[key], validate=False)
    return bool(pattern and mol.HasSubstructMatch(pattern))


def _has_substructure_any(mols: Sequence[Chem.Mol], key: str) -> bool:
    pattern = compile_smarts(_SMARTS_PATTERNS[key], validate=False)
    if not pattern:
        return False
    return any(mol.HasSubstructMatch(pattern) for mol in mols)


def _select_principal_pair(
    reactant_mols: Sequence[Chem.Mol],
    product_mols: Sequence[Chem.Mol],
) -> Tuple[int, int, str, int]:
    best_i = 0
    best_j = 0
    best_smarts = ""
    best_atoms = -1
    best_ratio = -1.0
    for i, reactant_mol in enumerate(reactant_mols):
        for j, product_mol in enumerate(product_mols):
            mcs = rdFMCS.FindMCS([reactant_mol, product_mol], timeout=4, ringMatchesRingOnly=True)
            atoms = int(getattr(mcs, "numAtoms", 0))
            largest = max(reactant_mol.GetNumAtoms(), product_mol.GetNumAtoms())
            ratio = (atoms / largest) if largest else 0.0
            if atoms > best_atoms or (atoms == best_atoms and ratio > best_ratio):
                best_i = i
                best_j = j
                best_smarts = getattr(mcs, "smartsString", "")
                best_atoms = atoms
                best_ratio = ratio
    return best_i, best_j, best_smarts, best_atoms


class ReactionAnalyzerPOC:
    """Standalone PoC analyzer with 'thinking-like' visible reasoning stages."""

    def __init__(self, *, use_llm: bool = False, provider: str = "openai", model: str = "gpt-5.2"):
        self.use_llm = use_llm
        self.provider = provider
        self.model = model

    def analyze(self, reaction_smiles: str) -> ReactionAnalysisResult:
        """Analyze single- or multi-component reaction SMILES."""
        if ">>" not in reaction_smiles:
            raise ReactionAnalysisError("Reaction SMILES must contain '>>'.")

        reactants_side, products_side = reaction_smiles.split(">>", 1)
        reactant_smiles_list = _split_side(reactants_side)
        product_smiles_list = _split_side(products_side)
        if not reactant_smiles_list or not product_smiles_list:
            raise ReactionAnalysisError("Reaction SMILES must contain at least one reactant and one product.")

        reactant_mols = [_require_mol(smiles, "reactant") for smiles in reactant_smiles_list]
        product_mols = [_require_mol(smiles, "product") for smiles in product_smiles_list]
        canonical_reactants = [Chem.MolToSmiles(mol) for mol in reactant_mols]
        canonical_products = [Chem.MolToSmiles(mol) for mol in product_mols]

        pair_i, pair_j, mcs_smarts, mcs_atoms = _select_principal_pair(reactant_mols, product_mols)
        core_reactant_mol = reactant_mols[pair_i]
        core_product_mol = product_mols[pair_j]
        core_reactant_smiles = reactant_smiles_list[pair_i]
        core_product_smiles = product_smiles_list[pair_j]
        core_canonical_reactant = canonical_reactants[pair_i]
        core_canonical_product = canonical_products[pair_j]

        reactant_formula = CalcMolFormula(core_reactant_mol)
        product_formula = CalcMolFormula(core_product_mol)
        formula_delta = _count_delta(reactant_formula, product_formula)

        side_reactant_formulas = [CalcMolFormula(mol) for mol in reactant_mols]
        side_product_formulas = [CalcMolFormula(mol) for mol in product_mols]
        side_reactant_formula = _counts_to_formula(_merge_counts(side_reactant_formulas))
        side_product_formula = _counts_to_formula(_merge_counts(side_product_formulas))
        side_formula_delta = _count_delta(side_reactant_formula, side_product_formula)

        features = {
            "reactant_count": len(reactant_mols),
            "product_count": len(product_mols),
            "core_reactant_halogen_count": _count_halogens(core_reactant_mol),
            "core_product_halogen_count": _count_halogens(core_product_mol),
            "side_reactant_halogen_count": sum(_count_halogens(mol) for mol in reactant_mols),
            "side_product_halogen_count": sum(_count_halogens(mol) for mol in product_mols),
            "core_reactant_aromatic_ring_n_count": _count_aromatic_ring_n(core_reactant_mol),
            "core_product_has_aryl_halide_like": _has_substructure(core_product_mol, "aryl_halide_like"),
            "core_reactant_has_aryl_halide_like": _has_substructure(core_reactant_mol, "aryl_halide_like"),
            "core_product_has_nn_single": _has_substructure(core_product_mol, "nn_single"),
            "core_product_has_nn_double": _has_substructure(core_product_mol, "nn_double"),
            "side_reactants_have_nn": _has_substructure_any(reactant_mols, "nn_single")
            or _has_substructure_any(reactant_mols, "nn_double"),
            "core_product_has_nh_aromatic": any(
                atom.GetSymbol() == "N"
                and atom.GetIsAromatic()
                and atom.GetTotalNumHs(includeNeighbors=True) > 0
                for atom in core_product_mol.GetAtoms()
            ),
        }

        selected_pair = {
            "reactant_index": pair_i,
            "product_index": pair_j,
            "reactant_smiles": core_reactant_smiles,
            "product_smiles": core_product_smiles,
        }

        hypotheses, best_hypothesis, confidence, steps = self._build_hypotheses(
            formula_delta=formula_delta,
            side_formula_delta=side_formula_delta,
            features=features,
            mcs_atoms=mcs_atoms,
            selected_pair=selected_pair,
        )

        summary = self._summarize(
            best_hypothesis=best_hypothesis,
            formula_delta=formula_delta,
            side_formula_delta=side_formula_delta,
            features=features,
        )
        result = ReactionAnalysisResult(
            reaction_smiles=reaction_smiles,
            reactant_smiles_list=reactant_smiles_list,
            product_smiles_list=product_smiles_list,
            canonical_reactants=canonical_reactants,
            canonical_products=canonical_products,
            reactant_smiles=core_reactant_smiles,
            product_smiles=core_product_smiles,
            canonical_reactant=core_canonical_reactant,
            canonical_product=core_canonical_product,
            reactant_formula=reactant_formula,
            product_formula=product_formula,
            formula_delta=formula_delta,
            side_reactant_formula=side_reactant_formula,
            side_product_formula=side_product_formula,
            side_formula_delta=side_formula_delta,
            selected_pair=selected_pair,
            mcs_smarts=mcs_smarts,
            mcs_atoms=mcs_atoms,
            features=features,
            hypotheses=hypotheses,
            best_hypothesis=best_hypothesis,
            confidence=confidence,
            summary=summary,
            reasoning_steps=steps,
        )

        if self.use_llm:
            result.llm_refinement = self._llm_refine(result)

        return result

    def _build_hypotheses(
        self,
        *,
        formula_delta: Dict[str, int],
        side_formula_delta: Dict[str, int],
        features: Dict[str, Any],
        mcs_atoms: int,
        selected_pair: Dict[str, Any],
    ) -> tuple[List[Dict[str, Any]], Dict[str, Any], float, List[AnalysisStep]]:
        nitrogen_gain = formula_delta.get("N", 0)
        halogen_loss = features["core_reactant_halogen_count"] > features["core_product_halogen_count"]
        has_nn = features["core_product_has_nn_single"] or features["core_product_has_nn_double"]

        snar_score = 0.0
        if halogen_loss:
            snar_score += 0.25
        if nitrogen_gain >= 2:
            snar_score += 0.20
        elif features["side_reactants_have_nn"] and has_nn:
            snar_score += 0.10
        if features["core_reactant_aromatic_ring_n_count"] >= 2:
            snar_score += 0.20
        if has_nn:
            snar_score += 0.20
        if mcs_atoms >= 10:
            snar_score += 0.10
        if features["side_reactants_have_nn"]:
            snar_score += 0.05
        snar_score = min(1.0, snar_score)

        hypotheses: List[Dict[str, Any]] = [
            {
                "reaction_class": "SNAr_hydrazinolysis",
                "score": round(snar_score, 2),
                "why": (
                    "Principal scaffold loses halogen on an electron-poor aza ring and installs an N-N fragment."
                ),
            },
            {
                "reaction_class": "other_amination_or_annotation_error",
                "score": round(max(0.0, 1.0 - snar_score), 2),
                "why": "Fallback when evidence for SNAr hydrazinolysis is incomplete.",
            },
        ]
        hypotheses.sort(key=lambda item: item["score"], reverse=True)
        best = hypotheses[0]

        steps = [
            AnalysisStep(
                name="Parse and normalize",
                evidence={
                    "valid_reaction": True,
                    "reactant_count": features["reactant_count"],
                    "product_count": features["product_count"],
                    "selected_pair": selected_pair,
                    "mcs_atoms": mcs_atoms,
                },
                conclusion="Reaction parsed; selected principal transformed pair via maximum shared scaffold.",
            ),
            AnalysisStep(
                name="Stoichiometric delta",
                evidence={"core_formula_delta": formula_delta, "side_formula_delta": side_formula_delta},
                conclusion="Compared both core-pair and whole-side stoichiometry for robustness.",
            ),
            AnalysisStep(
                name="Functional signal extraction",
                evidence={
                    "core_reactant_has_aryl_halide_like": features["core_reactant_has_aryl_halide_like"],
                    "core_product_has_nn_single_or_double": has_nn,
                    "core_reactant_aromatic_ring_n_count": features["core_reactant_aromatic_ring_n_count"],
                    "side_reactants_have_nn": features["side_reactants_have_nn"],
                    "core_product_has_nh_aromatic": features["core_product_has_nh_aromatic"],
                },
                conclusion="Signals are consistent with hydrazinyl substitution and tautomeric depiction.",
            ),
            AnalysisStep(
                name="Mechanistic ranking",
                evidence={"top_hypothesis": best["reaction_class"], "score": best["score"]},
                conclusion="SNAr hydrazinolysis is the most plausible mechanism class.",
            ),
        ]
        return hypotheses, best, float(best["score"]), steps

    def _summarize(
        self,
        *,
        best_hypothesis: Dict[str, Any],
        formula_delta: Dict[str, int],
        side_formula_delta: Dict[str, int],
        features: Dict[str, Any],
    ) -> str:
        if best_hypothesis["reaction_class"] != "SNAr_hydrazinolysis":
            return "Could not confidently assign SNAr hydrazinolysis from deterministic evidence."

        tautomer_note = (
            " Product aromatic [nH] indicates the reported structure is likely an imino/amino tautomer form."
            if features["core_product_has_nh_aromatic"]
            else ""
        )
        side_note = (
            " Side-level stoichiometry differs from core-pair delta, likely due to explicit reagents/byproducts."
            if side_formula_delta != formula_delta
            else ""
        )
        return (
            "Most likely transformation: aryl-halide substitution on an aza-aromatic ring by a hydrazine fragment "
            f"(core delta N={formula_delta.get('N', 0)}, core delta halogens="
            f"{features['core_product_halogen_count'] - features['core_reactant_halogen_count']})."
            + tautomer_note
            + side_note
        )

    def _llm_refine(self, result: ReactionAnalysisResult) -> Dict[str, Any]:
        try:
            from llmtools.clients import LLMClient
        except Exception as exc:  # pragma: no cover
            return {"enabled": False, "error": f"LLM client import failed: {exc}"}

        try:
            client = LLMClient(provider=self.provider, model=self.model, max_tokens=1200)
            prompt = (
                "You are refining a reaction-mechanism analysis.\n"
                "Use the deterministic evidence exactly as constraints.\n"
                "Return concise JSON with keys: final_call, confidence, alternatives, rationale.\n\n"
                f"Reaction: {result.reaction_smiles}\n"
                f"Selected pair: {result.selected_pair}\n"
                f"Core formula delta: {result.formula_delta}\n"
                f"Side formula delta: {result.side_formula_delta}\n"
                f"Features: {result.features}\n"
                f"Current best: {result.best_hypothesis}\n"
            )
            response = client.chat(prompt=prompt, temperature=0.2)
            return {
                "enabled": True,
                "provider": response.provider,
                "model": response.model,
                "tokens": response.total_tokens,
                "content": response.content,
            }
        except Exception as exc:  # pragma: no cover
            return {"enabled": True, "error": str(exc)}


def analyze_reaction(
    reaction_smiles: str,
    *,
    use_llm: bool = False,
    provider: str = "openai",
    model: str = "gpt-5.2",
) -> ReactionAnalysisResult:
    """Convenience function for one-shot analysis."""
    analyzer = ReactionAnalyzerPOC(use_llm=use_llm, provider=provider, model=model)
    return analyzer.analyze(reaction_smiles)
