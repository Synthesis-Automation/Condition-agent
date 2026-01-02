"""
Unified feature bundles for molecules and reactions.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from chemtools.smiles import normalize_reaction
from chemtools.taxonomy.reagent_v2 import ReagentTaxonomyV2
from chemtools.util import rdkit_helpers
from chemtools.util.functional_groups import detect_all, get_functional_groups

from .analysis.reaction_context import classify_reactants_with_context, get_reactant_summary
from .molecule import featurize_molecule as _featurize_molecule
from .reaction_detection import detect_reaction_types

_AGENT_ROLE_FLAGS = (
    "metal_catalyst",
    "organo_catalyst",
    "enzyme",
    "ligand",
    "base",
    "acid",
    "solvent",
    "additive",
    "oxidant",
    "reductant",
    "condensation_agent",
    "other_reagent",
)


@lru_cache(maxsize=1)
def _load_reagent_taxonomy() -> Optional[ReagentTaxonomyV2]:
    try:
        return ReagentTaxonomyV2.from_path()
    except Exception:
        return None


def _rdkit_props(smiles: str) -> Dict[str, Any]:
    mol = rdkit_helpers.parse_smiles(smiles)
    if mol is None:
        return {}
    try:
        from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    except Exception:
        return {}
    return {
        "molecular_weight": float(Descriptors.MolWt(mol)),
        "logP": float(Descriptors.MolLogP(mol)),
        "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
        "HBA": float(Lipinski.NumHAcceptors(mol)),
        "HBD": float(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": float(Lipinski.NumRotatableBonds(mol)),
        "fraction_csp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "ring_count": int(rdMolDescriptors.CalcNumRings(mol)),
        "aromatic_ring_count": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "heavy_atom_count": int(Descriptors.HeavyAtomCount(mol)),
    }


def _to_bool(value: Any, *, default: bool) -> bool:
    if value is None:
        return default
    return bool(value)


def _molecule_bundle(
    smiles: str,
    *,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    options = options or {}
    include_rdkit_props = _to_bool(options.get("include_rdkit_props"), default=True)
    include_functional_groups = _to_bool(options.get("include_functional_groups"), default=True)

    analysis = _featurize_molecule(smiles, registry_paths=registry_paths, options=options)

    fg_map = detect_all(smiles) if include_functional_groups else {}
    fg_list = get_functional_groups(smiles) if include_functional_groups else []
    rdkit_props = _rdkit_props(smiles) if include_rdkit_props else {}

    molecule = dict(analysis)
    molecule["functional_group_map"] = fg_map
    molecule["functional_groups"] = sorted(fg_list)
    molecule["rdkit_props"] = rdkit_props
    return molecule


def _agent_smiles(agent: Dict[str, Any]) -> str:
    return (
        agent.get("smiles_norm")
        or agent.get("largest_smiles")
        or agent.get("input")
        or ""
    )


def _classify_agents(agents: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
    taxonomy = _load_reagent_taxonomy()
    role_counts: Dict[str, int] = {}
    family_counts: Dict[str, int] = {}
    role_flags = {role_id: False for role_id in _AGENT_ROLE_FLAGS}
    entries: List[Dict[str, Any]] = []

    for agent in agents:
        smiles = _agent_smiles(agent)
        if not smiles:
            continue
        entry: Dict[str, Any] = {
            "smiles": smiles,
            "input": agent.get("input") or smiles,
        }
        match = None
        if taxonomy is not None:
            match = taxonomy.classify({"smiles": smiles, "name": None, "cas": None})
        if match is not None:
            role_flags[match.role_id] = True
            role_counts[match.role_id] = role_counts.get(match.role_id, 0) + 1
            family_counts[match.family_id] = family_counts.get(match.family_id, 0) + 1
            role = taxonomy.get_role(match.role_id) if taxonomy else None
            family = taxonomy.get_family(match.family_id) if taxonomy else None
            entry.update(
                {
                    "role_id": match.role_id,
                    "role_name": role.name if role else match.role_id,
                    "family_id": match.family_id,
                    "family_name": family.name if family else match.family_id,
                    "match_kind": match.match_kind,
                }
            )
        entries.append(entry)

    flags = {
        "has_catalyst": any(
            role_flags.get(role_id)
            for role_id in ("metal_catalyst", "organo_catalyst", "enzyme")
        ),
        "has_metal_catalyst": role_flags.get("metal_catalyst", False),
        "has_organo_catalyst": role_flags.get("organo_catalyst", False),
        "has_enzyme": role_flags.get("enzyme", False),
        "has_ligand": role_flags.get("ligand", False),
        "has_base": role_flags.get("base", False),
        "has_acid": role_flags.get("acid", False),
        "has_solvent": role_flags.get("solvent", False),
        "has_additive": role_flags.get("additive", False),
        "has_oxidant": role_flags.get("oxidant", False),
        "has_reductant": role_flags.get("reductant", False),
        "has_condensation_agent": role_flags.get("condensation_agent", False),
    }

    return {
        "agents": entries,
        "agent_count": len(entries),
        "role_counts": role_counts,
        "family_counts": family_counts,
        "role_flags": role_flags,
        "flags": flags,
    }


def featurize_molecule(
    smiles: str,
    *,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Return a unified molecule feature bundle.
    """
    molecule = _molecule_bundle(
        smiles,
        registry_paths=registry_paths,
        options=options,
    )
    meta = {
        "rdkit_available": rdkit_helpers.rdkit_available(),
        "errors": [molecule.get("meta", {}).get("error")] if molecule.get("meta", {}).get("error") else [],
    }
    return {
        "schema_version": "v1",
        "kind": "molecule",
        "molecule": molecule,
        "meta": meta,
    }


def _reaction_type_summary(detection: Any) -> Dict[str, Any]:
    matches = detection.matches if detection else []
    if not matches:
        return {"reaction_type": "Unknown", "confidence": 0.0, "slot_evidence": {}}
    best = matches[0]
    required = int(best.required_slots or 0)
    matched = int(best.matched_slots or 0)
    confidence = float(matched) / float(required) if required else 0.0
    return {
        "reaction_type": best.reaction_type,
        "name": best.name,
        "category": best.category,
        "confidence": confidence,
        "slot_evidence": {slot: list(values) for slot, values in best.slot_evidence.items()},
    }


def _extract_scores(result: Any) -> List[float]:
    scores: List[float] = []
    if isinstance(result, dict):
        score = result.get("score_0_10")
        if isinstance(score, (int, float)):
            scores.append(float(score))
    elif isinstance(result, list):
        for entry in result:
            if isinstance(entry, dict):
                score = entry.get("score_0_10")
                if isinstance(score, (int, float)):
                    scores.append(float(score))
    return scores


def _aggregate_reaction_features(reactants: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
    reactant_list = list(reactants)
    aryl_scores: List[float] = []
    alkyl_scores: List[float] = []
    electronic_scores: List[float] = []
    motifs: set[str] = set()
    functional_groups: set[str] = set()

    for reactant in reactant_list:
        for entry in reactant.get("steric", {}).get("aryl", []):
            aryl_scores.extend(_extract_scores(entry.get("result")))
        for entry in reactant.get("steric", {}).get("alkyl", []):
            alkyl_scores.extend(_extract_scores(entry.get("result")))
        for entry in reactant.get("electronics", {}).get("aryl", []):
            electronic_scores.extend(_extract_scores(entry.get("result")))
        for motif in reactant.get("motifs", []):
            compound_id = motif.get("compound_id")
            if compound_id:
                motifs.add(str(compound_id))
        fg_map = reactant.get("functional_group_map") or {}
        for token, present in fg_map.items():
            if present:
                functional_groups.add(token)

    avg_electronic = None
    if electronic_scores:
        avg_electronic = round(sum(electronic_scores) / len(electronic_scores), 2)

    return {
        "reactant_count": len(reactant_list),
        "motif_ids": sorted(motifs),
        "functional_group_ids": sorted(functional_groups),
        "max_aryl_steric": max(aryl_scores) if aryl_scores else 0.0,
        "max_alkyl_steric": max(alkyl_scores) if alkyl_scores else 0.0,
        "avg_aryl_electronic": avg_electronic if avg_electronic is not None else 5.0,
        "electron_poor_aryl": any(score > 6.5 for score in electronic_scores),
    }


def featurize_reaction(
    reaction_smiles: str,
    *,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Return a unified reaction feature bundle (reaction-level + per-reactant bundles).
    """
    options = options or {}
    include_roles = _to_bool(options.get("include_roles"), default=True)
    include_agent_roles = _to_bool(options.get("include_agent_roles"), default=True)

    detection = detect_reaction_types(
        reaction_smiles,
        max_hits_per_compound=(options or {}).get("max_hits_per_compound"),
    )
    detection_payload = detection.to_dict()
    reaction_type = _reaction_type_summary(detection)

    normalized = normalize_reaction(reaction_smiles)
    reactant_smiles = [
        item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
        for item in (normalized.get("reactants") or [])
    ]
    reactant_smiles = [s for s in reactant_smiles if s]

    agent_roles = None
    if include_agent_roles:
        agent_roles = _classify_agents(normalized.get("agents") or [])

    reactant_bundles = [
        _molecule_bundle(smiles, registry_paths=registry_paths, options=options)
        for smiles in reactant_smiles
    ]

    aggregates = _aggregate_reaction_features(reactant_bundles)

    roles_summary = None
    if include_roles:
        try:
            roles = classify_reactants_with_context(reaction_smiles)
            roles_summary = get_reactant_summary(roles)
        except Exception:
            roles_summary = None

    reaction = {
        "reaction_smiles": reaction_smiles,
        "normalized": normalized,
        "detection": detection_payload,
        "reaction_type": reaction_type,
        "reactants": reactant_bundles,
        "aggregates": aggregates,
        "roles": roles_summary,
        "agent_roles": agent_roles,
    }

    meta = {
        "rdkit_available": rdkit_helpers.rdkit_available(),
        "errors": [detection_payload.get("error")] if detection_payload.get("error") else [],
    }
    return {
        "schema_version": "v1",
        "kind": "reaction",
        "reaction": reaction,
        "meta": meta,
    }


__all__ = ["featurize_molecule", "featurize_reaction"]
