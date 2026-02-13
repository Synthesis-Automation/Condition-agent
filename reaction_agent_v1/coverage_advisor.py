"""Coverage expansion advisor for taxonomy and tool gaps."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any, Dict, List

from rdkit import Chem

from chemtools.taxonomy.loader import load_reaction_types_dict

from .contracts import CoverageSuggestionCard


def _count_aromatic_ring_n(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return sum(
        1
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "N" and atom.GetIsAromatic() and atom.IsInRing()
    )


class CoverageAdvisor:
    """Suggest high-impact taxonomy and tool improvements from unknown decisions."""

    def suggest(self, reaction_smiles: str, analysis: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Generate ranked suggestion cards for coverage expansion."""
        decision = analysis.get("decision") or {}
        if str(decision.get("reaction_type")) != "unknown":
            return []

        suggestions: List[CoverageSuggestionCard] = []
        core_delta = analysis.get("core_formula_delta") or {}
        reacted_motifs = analysis.get("reacted_motifs") or []
        formed_motifs = analysis.get("formed_motifs") or []
        candidates = analysis.get("taxonomy_candidates") or []
        principal_pair = analysis.get("principal_pair") or {}
        reactant_smiles = str(principal_pair.get("reactant_smiles") or "")

        halogen_loss = any(int(core_delta.get(symbol, 0)) < 0 for symbol in ("F", "Cl", "Br", "I"))
        nitrogen_gain = int(core_delta.get("N", 0)) > 0
        ring_n = _count_aromatic_ring_n(reactant_smiles) if reactant_smiles else 0

        taxonomy = load_reaction_types_dict()
        taxonomy_ids = {key.lower() for key in taxonomy}

        if halogen_loss and nitrogen_gain and ring_n >= 2 and "snar_cn" in taxonomy_ids:
            suggestions.append(
                CoverageSuggestionCard(
                    suggestion_id="taxonomy-snar-cn-gap",
                    kind="taxonomy",
                    priority=100,
                    title="Improve SNAr_CN detector coverage for aza-aryl hydrazinolysis-like reactions",
                    rationale=(
                        "Core delta shows halogen loss and nitrogen gain on an electron-poor aza aromatic substrate, "
                        "but taxonomy detector returned no candidate."
                    ),
                    proposed_changes=[
                        "Add/expand motif constraints in SNAr_CN for diazine-like activated heteroaryl halides.",
                        "Add regression reaction examples from unknown queue into SNAr_CN validation set.",
                    ],
                    evidence={
                        "reaction_smiles": reaction_smiles,
                        "core_formula_delta": core_delta,
                        "reactant_aromatic_ring_n_count": ring_n,
                    },
                )
            )

        if not reacted_motifs and not formed_motifs:
            suggestions.append(
                CoverageSuggestionCard(
                    suggestion_id="tool-motif-extraction-gap",
                    kind="tool",
                    priority=90,
                    title="Add bond-change fallback when motif extraction returns empty",
                    rationale="Unknown decisions with empty motif deltas indicate upstream detection blind spots.",
                    proposed_changes=[
                        "Add secondary bond-change signature extraction tool for candidate retrieval fallback.",
                        "Route empty-motif cases to fallback detector before final unknown decision.",
                    ],
                    evidence={
                        "reaction_smiles": reaction_smiles,
                        "reacted_motifs": reacted_motifs,
                        "formed_motifs": formed_motifs,
                    },
                )
            )

        if not candidates:
            suggestions.append(
                CoverageSuggestionCard(
                    suggestion_id="tool-candidate-retrieval-gap",
                    kind="tool",
                    priority=85,
                    title="Add secondary candidate retrieval path for zero-hit taxonomy detections",
                    rationale=(
                        "Detector returned zero taxonomy candidates; add fallback retrieval from "
                        "bond-change signatures and principal-pair formula deltas."
                    ),
                    proposed_changes=[
                        "Implement fallback candidate retrieval tool using bond-change and delta signatures.",
                        "Feed fallback candidates into validator-gated ranking before final unknown.",
                    ],
                    evidence={
                        "reaction_smiles": reaction_smiles,
                        "candidate_count": len(candidates),
                    },
                )
            )

        if not suggestions:
            suggestions.append(
                CoverageSuggestionCard(
                    suggestion_id="taxonomy-generic-unknown-cluster",
                    kind="taxonomy",
                    priority=50,
                    title="Cluster unknown reactions for taxonomy expansion",
                    rationale="No strong single-family signature; prioritize unknown-cluster triage workflow.",
                    proposed_changes=[
                        "Cluster by principal-pair delta + reagent role patterns + reaction key fragments.",
                        "Promote top cluster to candidate taxonomy rule/tool update.",
                    ],
                    evidence={"reaction_smiles": reaction_smiles},
                )
            )

        suggestions.sort(key=lambda row: (-row.priority, row.suggestion_id))
        return [asdict(row) for row in suggestions]
