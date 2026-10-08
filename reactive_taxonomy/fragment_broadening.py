"""Explicit target-validated alternatives to a chosen fragment query.

Generated choices describe specific relaxations, not a universal SMARTS subset
proof. Atom selections refer to the input query order; target alignments remain
structural hypotheses, including symmetry alternatives.
"""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any

from rdkit import Chem

from .chemistry.smarts_cache import compile_smarts
from .fragment_search import compile_fragment_query, fragment_embeddings, validate_fragment_target
from .search_fragments import suggest_search_fragments


def broadening_policy() -> dict[str, Any]:
    """Load the versioned bounded query-edit vocabulary."""
    value = json.loads((Path(__file__).parent / "definitions/fragment_broadening.v1.json").read_text("utf-8"))
    if (value.get("schema_version") != "fragment_broadening_policy.v1"
            or value.get("definition_version") != "fragment_broadening.v1@1.0"
            or value.get("relaxations") != ["peripheral_context", "ring_boundary", "aromatic_carbon_nitrogen"]
            or value.get("aromatic_alternative") != "[c,n;+0]"):
        raise ValueError("Unsupported fragment broadening policy")
    for key, maximum in (("max_variants", 7), ("max_target_alignments", 64)):
        if type(value.get(key)) is not int or not 1 <= value[key] <= maximum:
            raise ValueError(f"Invalid fragment broadening limit: {key}")
    compile_smarts(value["aromatic_alternative"], validate=True)
    return value


def propose_fragment_queries(
    target_smiles: str, query: str, query_format: str = "smiles",
    topology: str = "preserve_rings", aromatic_atom_ids: list[int] | None = None,
) -> dict[str, Any]:
    """Preview bounded query alternatives without searching or inferring reactions.

Automatic peripheral edits require SMILES. C/N alternatives require explicit
neutral aromatic c/n positions without isotope, explicit H or stereo constraints.
Custom SMARTS remains supported, but is never rewritten automatically.
"""
    policy = broadening_policy()
    parent = compile_fragment_query(query, query_format, topology)
    validation = validate_fragment_target(parent, target_smiles)
    if not validation.matches_target:
        raise ValueError("Parent query does not match target_smiles; correct it before broadening")
    target = Chem.MolFromSmiles(validation.target_smiles)
    choices: list[dict[str, Any]] = []
    seen = {parent.query_id}

    def add(expression: str, fmt: str, mode: str, edits: list[str], reason: str) -> None:
        compiled = compile_fragment_query(expression, fmt, mode)
        check = validate_fragment_target(compiled, validation.target_smiles)
        if not check.matches_target or compiled.query_id in seen:
            return
        alignments, truncated = fragment_embeddings(compiled, target, maximum=policy["max_target_alignments"])
        seen.add(compiled.query_id)
        identity = [parent.query_id, compiled.query_id, policy["definition_version"], validation.target_smiles]
        choices.append({
            "variant_id": "FBQ1:" + hashlib.sha256(json.dumps(identity).encode()).hexdigest(),
            "parent_query_id": parent.query_id, "query": expression, "query_format": fmt,
            "topology": mode, "relaxations": edits, "reason": reason,
            "target_validation": asdict(check), "target_alignments": [list(a) for a in alignments],
            "target_atom_ids": sorted(alignments[0]), "target_alignments_truncated": truncated,
            "alignment_ambiguous": len(alignments) > 1 or truncated,
        })

    if topology == "preserve_rings":
        add(query, query_format, "subgraph", ["ring_boundary"],
            "Allow additional ring fusion and relax ring membership at query boundaries; retain atom and bond constraints.")
    eligible: list[int] = []
    atoms: list[dict[str, Any]] = []
    if query_format == "smiles":
        original = Chem.MolFromSmiles(query)
        for atom in original.GetAtoms():
            allowed = (atom.GetIsAromatic() and atom.GetAtomicNum() in (6, 7)
                       and not atom.GetFormalCharge() and not atom.GetIsotope()
                       and not atom.GetNumExplicitHs()
                       and atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED)
            if allowed:
                eligible.append(atom.GetIdx())
            atoms.append({"query_atom_id": atom.GetIdx(), "element": atom.GetSymbol(),
                          "aromatic": atom.GetIsAromatic(), "allows_carbon_nitrogen": allowed})
        for suggestion in suggest_search_fragments(query).candidates:
            if len(suggestion.target_atom_ids) < original.GetNumAtoms():
                add(suggestion.query, "smiles", topology, ["peripheral_context"],
                    "Remove peripheral context while retaining complete selected rings and required valence/stereo context; omitted groups become unconstrained.")
    if aromatic_atom_ids is not None:
        if (not isinstance(aromatic_atom_ids, list) or not aromatic_atom_ids
                or any(type(i) is not int or i not in eligible for i in aromatic_atom_ids)
                or len(set(aromatic_atom_ids)) != len(aromatic_atom_ids)):
            raise ValueError("Select unique eligible aromatic c/n atom IDs from the SMILES query")
        mol = Chem.RWMol(compile_fragment_query(query, "smiles", "subgraph").molecule)
        for index in sorted(aromatic_atom_ids):
            pattern = compile_smarts(policy["aromatic_alternative"], validate=True)
            mol.ReplaceAtom(index, pattern.GetAtomWithIdx(0))
        edits = (["ring_boundary"] if topology == "preserve_rings" else []) + ["aromatic_carbon_nitrogen"]
        add(Chem.MolToSmarts(mol), "smarts", "subgraph", edits,
            f"Allow neutral aromatic carbon or nitrogen at query atoms {sorted(aromatic_atom_ids)}; retain other atoms and bonds. Ring-boundary restrictions are relaxed.")
    return {
        "schema_version": "fragment_query_alternatives.v1", "definition_version": policy["definition_version"],
        "target_smiles": validation.target_smiles, "parent_query": parent.describe(),
        "target_validation": asdict(validation), "query_atoms": atoms,
        "aromatic_atom_ids": sorted(aromatic_atom_ids or []),
        "variants": choices[:policy["max_variants"]], "output_truncated": len(choices) > policy["max_variants"],
        "limitations": ["Choices are explicit query edits, not automatic route recommendations.",
                        "The first target alignment is highlighted; alternative alignments remain hypotheses.",
                        "Custom SMARTS is not automatically simplified; edit its constraints explicitly."],
    }
