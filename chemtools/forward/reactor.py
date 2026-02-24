"""
Forward synthesis reaction engine.

ReactantAnalyzer  — analyze a reactant pair for compatible forward templates
ForwardReactor    — apply forward SMARTS templates via AllChem.RunReactants()

Design mirrors chemtools.retro.disconnector:
  disconnector: target_mol → RunReactants(retro_smarts) → precursor SMILES
  reactor:      reactant_mol(s) → RunReactants(fwd_smarts) → product SMILES
"""
from __future__ import annotations

import logging
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple

from .contracts import ForwardTemplateMatch, ProductPrediction, ReactantProfile

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Module-level SMARTS cache — avoids recompiling the same pattern repeatedly
# ---------------------------------------------------------------------------

@lru_cache(maxsize=256)
def _compile_smarts(smarts: str) -> Optional[Any]:
    """Compile and cache a SMARTS pattern (returns RDKit Mol or None)."""
    try:
        from rdkit import Chem, rdBase
        with rdBase.BlockLogs():
            return Chem.MolFromSmarts(smarts)
    except Exception:
        return None


@lru_cache(maxsize=256)
def _compile_reaction(smarts: str) -> Optional[Any]:
    """Compile and cache a forward reaction SMARTS."""
    try:
        from rdkit.Chem import AllChem, rdBase
        from rdkit import rdBase as _rdBase
        with _rdBase.BlockLogs():
            rxn = AllChem.ReactionFromSmarts(smarts)
        return rxn
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Functional group SMARTS for reactant analysis
# ---------------------------------------------------------------------------

#: Primary FG patterns keyed by label — used for detect_reactive_groups()
_FG_SMARTS: Dict[str, str] = {
    # Electrophilic leaving groups / electrophiles
    "aryl_bromide":      "[c][Br]",
    "aryl_chloride":     "[c][Cl]",
    "aryl_iodide":       "[c][I]",
    "aryl_fluoride":     "[c][F]",
    "aryl_triflate":     "[c]OC(F)(F)F",
    "alkyl_bromide":     "[CX4][Br]",
    "alkyl_chloride":    "[CX4][Cl]",
    "acyl_chloride":     "[C](=O)[Cl]",
    "sulfonyl_chloride": "[S](=O)(=O)[Cl]",
    "chloroformate":     "[C](=O)O[Cl]",
    "arylboronic_acid":  "[c][B](O)O",
    "boronic_ester":     "[c][B]1OC(C)(C)C(C)(C)O1",
    "arylstannane":      "[c][Sn]",
    "arylsilane":        "[c][Si]",
    "vinyl_halide":      "[C]=[C][Br,Cl,I]",
    "isocyanate":        "[N]=[C]=O",
    "isothiocyanate":    "[N]=[C]=S",
    "aldehyde":          "[CH1](=O)",
    "ketone":            "[CX3;!$(C=O)O](=O)[#6]",
    "carboxylic_acid":   "[C](=O)[OH]",
    "ester":             "[C](=O)O[CX4]",
    "michael_acceptor":  "[C]=[C][C](=O)",
    "thioester":         "[C](=O)[S]",
    "aryl_azide":        "[c][N3]",
    "alkyl_azide":       "[CX4][N3]",
    "terminal_alkyne":   "[C]#[CH]",
    "terminal_alkene":   "[C]=[CH2]",
    "alkyl_halide_gen":  "[CX4][Br,Cl,I]",
    # Nucleophiles
    "primary_amine":     "[NH2][#6;!$(C=O)]",
    "secondary_amine":   "[NH1]([#6])[#6;!$(C=O)]",
    "aniline":           "[c][NH2]",
    "aryl_amine":        "[c][NH]",
    "alcohol":           "[CX4][OH]",
    "phenol":            "[c][OH]",
    "thiol":             "[#16H]",
    "boronic_acid":      "[B](O)O",
    "alkyne":            "[C]#[C]",
    "active_methylene":  "[CH2]([C]=O)[C]=O",
    "organozinc":        "[CX4][Zn]",
    "grignard":          "[CX4][Mg]",
}

#: Map FG label → role (electrophile / nucleophile / ambiphile)
_FG_ROLE: Dict[str, str] = {
    "aryl_bromide": "electrophile",
    "aryl_chloride": "electrophile",
    "aryl_iodide": "electrophile",
    "aryl_fluoride": "electrophile",
    "aryl_triflate": "electrophile",
    "alkyl_bromide": "electrophile",
    "alkyl_halide_gen": "electrophile",
    "alkyl_chloride": "electrophile",
    "acyl_chloride": "electrophile",
    "sulfonyl_chloride": "electrophile",
    "chloroformate": "electrophile",
    "arylboronic_acid": "ambiphile",
    "boronic_ester": "ambiphile",
    "arylstannane": "ambiphile",
    "arylsilane": "ambiphile",
    "vinyl_halide": "electrophile",
    "isocyanate": "electrophile",
    "isothiocyanate": "electrophile",
    "aldehyde": "electrophile",
    "ketone": "electrophile",
    "carboxylic_acid": "electrophile",
    "ester": "electrophile",
    "michael_acceptor": "electrophile",
    "thioester": "ambiphile",
    "aryl_azide": "ambiphile",
    "alkyl_azide": "ambiphile",
    "terminal_alkyne": "ambiphile",
    "terminal_alkene": "ambiphile",
    "primary_amine": "nucleophile",
    "secondary_amine": "nucleophile",
    "aniline": "nucleophile",
    "aryl_amine": "nucleophile",
    "alcohol": "nucleophile",
    "phenol": "nucleophile",
    "thiol": "nucleophile",
    "boronic_acid": "ambiphile",
    "alkyne": "ambiphile",
    "active_methylene": "ambiphile",
    "organozinc": "nucleophile",
    "grignard": "nucleophile",
}


def _canonical(smiles: str) -> str:
    try:
        from rdkit import Chem, rdBase
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol)
    except Exception:
        return smiles


def _parse_mol(smiles: str) -> Optional[Any]:
    try:
        from rdkit import Chem, rdBase
        with rdBase.BlockLogs():
            return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# ReactantAnalyzer
# ---------------------------------------------------------------------------

class ReactantAnalyzer:
    """Analyze one or two reactant molecules for forward synthesis compatibility.

    Usage::

        analyzer = ReactantAnalyzer()
        profile = analyzer.profile("Brc1ccccc1")
        compatible = analyzer.find_compatible_templates("Brc1ccccc1", "Nc1ccccc1")
    """

    def profile(self, smiles: str) -> ReactantProfile:
        """Build a ReactantProfile for one molecule."""
        from rdkit import Chem, rdBase
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol_smiles = smiles.split(">>")[0].split(">")[0].strip() if ">" in smiles else smiles
        mol = _parse_mol(mol_smiles)
        if mol is None:
            raise ValueError(f"Cannot parse SMILES: {mol_smiles!r}")

        # Canonical SMILES
        with rdBase.BlockLogs():
            canonical = Chem.MolToSmiles(mol)

        # Properties
        mw = round(Descriptors.MolWt(mol), 1)
        ha = mol.GetNumHeavyAtoms()

        # FG detection
        fgs: List[str] = []
        for label, smarts in _FG_SMARTS.items():
            pat = _compile_smarts(smarts)
            if pat is not None and mol.HasSubstructMatch(pat):
                fgs.append(label)

        # Group by role
        elec, nuc = [], []
        for fg in fgs:
            role = _FG_ROLE.get(fg, "ambiphile")
            if role == "electrophile":
                elec.append(fg)
            elif role == "nucleophile":
                nuc.append(fg)

        # Category map (simplified)
        fg_cats: Dict[str, List[str]] = {}
        for fg in fgs:
            role = _FG_ROLE.get(fg, "ambiphile")
            fg_cats.setdefault(role, []).append(fg)

        # Leaving group
        lg = None
        if "aryl_bromide" in fgs:
            lg = "Br"
        elif "aryl_chloride" in fgs:
            lg = "Cl"
        elif "aryl_iodide" in fgs:
            lg = "I"
        elif "aryl_fluoride" in fgs:
            lg = "F"
        elif "alkyl_bromide" in fgs or "alkyl_halide_gen" in fgs:
            lg = "Br/Cl"

        return ReactantProfile(
            smiles=mol_smiles,
            canonical_smiles=canonical,
            molecular_weight=mw,
            heavy_atoms=ha,
            functional_groups=fgs,
            fg_categories=fg_cats,
            electrophilic_sites=elec,
            nucleophilic_sites=nuc,
            leaving_group=lg,
        )

    def find_compatible_templates(
        self,
        smiles_a: str,
        smiles_b: str = "",
    ) -> List[ForwardTemplateMatch]:
        """Return all HTE forward templates compatible with the reactant(s).

        Compatibility is checked by verifying that the forward SMARTS LHS
        patterns match the supplied reactant SMILES using RDKit substructure.
        """
        from rdkit import Chem, rdBase

        mol_a = _parse_mol(smiles_a)
        mol_b = _parse_mol(smiles_b) if smiles_b else None
        if mol_a is None:
            return []

        # Load templates
        try:
            from chemtools.retro.hte_templates import HTE_TEMPLATES, get_template_taxonomy_id
        except ImportError:
            return []

        results: List[ForwardTemplateMatch] = []
        for tmpl in HTE_TEMPLATES:
            fwd_smarts = tmpl.get("forward_smarts")
            if not fwd_smarts:
                continue
            n_r = tmpl.get("n_reactants", tmpl.get("n_precursors", 2))

            # When a second reactant is provided, only test bimolecular templates.
            # Unimolecular templates (n_r==1) should only fire when smiles_b is
            # empty — otherwise the user supplied two reactants and expects a
            # coupling reaction, not a single-reactant transformation.
            if mol_b is not None and n_r == 1:
                continue

            # Try to compile the reaction and test with RunReactants
            rxn = _compile_reaction(fwd_smarts)
            if rxn is None:
                continue

            try:
                if n_r == 1:
                    with rdBase.BlockLogs():
                        out = rxn.RunReactants((mol_a,))
                    match = bool(out)
                else:
                    # Try both orderings of A and B
                    match = False
                    if mol_b is not None:
                        with rdBase.BlockLogs():
                            out = rxn.RunReactants((mol_a, mol_b))
                        if out:
                            match = True
                        if not match:
                            with rdBase.BlockLogs():
                                out = rxn.RunReactants((mol_b, mol_a))
                            if out:
                                match = True
            except Exception:
                continue

            if match:
                taxonomy_id = get_template_taxonomy_id(tmpl)
                results.append(ForwardTemplateMatch(
                    template_name=tmpl["name"],
                    taxonomy_id=taxonomy_id,
                    category=tmpl.get("category", ""),
                    forward_smarts=fwd_smarts,
                    n_reactants=n_r,
                    difficulty=tmpl.get("difficulty", 0.5),
                    description=tmpl.get("description", ""),
                    hte_families=tmpl.get("hte_families", []),
                    notes=tmpl.get("notes", ""),
                ))

        # Sort by difficulty (easiest first)
        results.sort(key=lambda t: t.difficulty)
        return results


# ---------------------------------------------------------------------------
# ForwardReactor
# ---------------------------------------------------------------------------

class ForwardReactor:
    """Apply forward SMARTS templates to reactant(s) and collect product SMILES.

    Usage::

        reactor = ForwardReactor()
        predictions = reactor.generate("Brc1ccccc1", "Nc1ccccc1", top_k=5)
    """

    _analyzer = ReactantAnalyzer()

    def generate(
        self,
        smiles_a: str,
        smiles_b: str = "",
        reaction_name: str = "",
        top_k: int = 5,
    ) -> List[ProductPrediction]:
        """Generate product predictions for the given reactant(s).

        Args:
            smiles_a: First reactant SMILES.
            smiles_b: Second reactant SMILES (empty for unimolecular reactions).
            reaction_name: Optional filter — only apply templates matching
                           this name or hte_family (e.g. "suzuki_miyaura").
            top_k: Maximum number of predictions to return.

        Returns:
            List[ProductPrediction] sorted by overall_score descending.
        """
        from rdkit import Chem, rdBase

        mol_a = _parse_mol(smiles_a)
        mol_b = _parse_mol(smiles_b) if smiles_b else None
        if mol_a is None:
            return []

        # Load templates
        try:
            from chemtools.retro.hte_templates import HTE_TEMPLATES, get_template_taxonomy_id
        except ImportError:
            return []

        templates = HTE_TEMPLATES
        if reaction_name:
            rn_lower = reaction_name.lower()
            templates = [
                t for t in templates
                if t["name"].lower() == rn_lower
                or any(rn_lower in f.lower() for f in t.get("hte_families", []))
            ]

        # Find templates that match reactant FGs (via compatibility check)
        compatible_names = {
            t.template_name
            for t in self._analyzer.find_compatible_templates(smiles_a, smiles_b)
        }

        predictions: List[ProductPrediction] = []

        for tmpl in templates:
            # Skip templates not compatible unless explicit reaction_name given
            if not reaction_name and tmpl["name"] not in compatible_names:
                continue

            fwd_smarts = tmpl.get("forward_smarts")
            if not fwd_smarts:
                continue
            n_r = tmpl.get("n_reactants", tmpl.get("n_precursors", 2))

            rxn = _compile_reaction(fwd_smarts)
            if rxn is None:
                continue

            products = self._run_forward_rxn(rxn, mol_a, mol_b, n_r)
            if not products:
                continue

            # Determine canonical reactants string
            can_a = _canonical(smiles_a)
            can_b = _canonical(smiles_b) if smiles_b else ""
            reactants_str = f"{can_a}.{can_b}" if can_b else can_a

            # Build ProductPrediction for the best (first) product;
            # store all regioisomers in all_product_smiles
            primary = products[0]
            rxn_smi = f"{reactants_str}>>{primary}"

            pred = ProductPrediction(
                product_smiles=primary,
                reactant_a=can_a,
                reactant_b=can_b,
                template_name=tmpl["name"],
                taxonomy_id=get_template_taxonomy_id(tmpl),
                reaction_smiles=rxn_smi,
                difficulty=tmpl.get("difficulty", 0.5),
                description=tmpl.get("description", ""),
                notes=tmpl.get("notes", ""),
                hte_families=tmpl.get("hte_families", []),
                all_product_smiles=products,
            )

            # Count new stereocenters in product vs reactants
            pred.new_stereocenters = self._count_new_stereocenters(
                primary, smiles_a, smiles_b
            )

            predictions.append(pred)
            if len(predictions) >= top_k * 3:   # over-generate, then trim after scoring
                break

        return predictions

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _run_forward_rxn(
        self,
        rxn: Any,
        mol_a: Any,
        mol_b: Optional[Any],
        n_reactants: int,
    ) -> List[str]:
        """Run RunReactants and return deduplicated product SMILES."""
        from rdkit import Chem, rdBase

        seen: set = set()
        results: List[str] = []

        reactant_sets: List[Tuple] = []
        if n_reactants == 1:
            reactant_sets = [(mol_a,)]
        else:
            if mol_b is not None:
                reactant_sets = [(mol_a, mol_b), (mol_b, mol_a)]
            else:
                reactant_sets = [(mol_a,)]

        for reactants in reactant_sets:
            try:
                with rdBase.BlockLogs():
                    out_tuples = rxn.RunReactants(reactants)
            except Exception:
                continue
            for tup in out_tuples:
                if not tup:
                    continue
                # Take the first component (primary product)
                prod_mol = tup[0]
                if prod_mol is None:
                    continue
                try:
                    with rdBase.BlockLogs():
                        smi = Chem.MolToSmiles(prod_mol)
                    if not smi:
                        continue
                    # Validate round-trip — try full sanitize first, then
                    # fall back to sanitize=False to accept charged aromatic
                    # systems (e.g. CuAAC triazole [n+][n-]) that MolToSmiles
                    # can generate but MolFromSmiles may reject on kekulization.
                    with rdBase.BlockLogs():
                        check_mol = Chem.MolFromSmiles(smi)
                    if check_mol is None:
                        with rdBase.BlockLogs():
                            check_mol = Chem.MolFromSmiles(smi, sanitize=False)
                        if check_mol is None:
                            continue
                except Exception:
                    continue
                if smi not in seen:
                    seen.add(smi)
                    results.append(smi)
                if len(results) >= 8:
                    return results

        return results

    def _count_new_stereocenters(
        self,
        product_smiles: str,
        smiles_a: str,
        smiles_b: str,
    ) -> int:
        """Count stereocenters in product not present in either reactant."""
        try:
            from rdkit.Chem import rdMolDescriptors
            prod_sc = rdMolDescriptors.CalcNumAtomStereoCenters(_parse_mol(product_smiles) or _parse_mol("C"))
            a_sc = rdMolDescriptors.CalcNumAtomStereoCenters(_parse_mol(smiles_a) or _parse_mol("C"))
            b_sc = rdMolDescriptors.CalcNumAtomStereoCenters(_parse_mol(smiles_b) or _parse_mol("C")) if smiles_b else 0
            return max(0, prod_sc - a_sc - b_sc)
        except Exception:
            return 0
