"""Explicit fragment queries and evidence projection onto product embeddings.

No correspondence is inferred here. Relationship evidence uses validated supplied
maps only; other observations remain discoverable with unresolved relationships.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from rdkit import Chem

from .chemistry.smarts_cache import compile_smarts
from .reaction_models import ReactionAtomReference


QUERY_COMPILER_VERSION = "fragment_query_compiler.v3"


def fragment_search_policy() -> dict[str, Any]:
    """Load and validate the versioned search limits and evidence vocabulary."""
    value = json.loads((Path(__file__).parent / "definitions/fragment_search.v1.json").read_text("utf-8"))
    if value.get("schema_version") != "fragment_search_policy.v1":
        raise ValueError("Unsupported fragment search policy")
    for key in ("max_query_characters", "max_query_atoms", "max_embeddings",
                "max_matched_products", "max_observations", "search_batch_size",
                "max_result_bytes", "max_procedure_characters"):
        if type(value.get(key)) is not int or value[key] < 1:
            raise ValueError(f"Invalid fragment search policy: {key}")
    if value["relationship_order"] != ["constructed", "modified", "boundary_changed",
                                        "carried_through", "reported_use", "unresolved"]:
        raise ValueError("Invalid fragment relationship vocabulary")
    return value


@dataclass(frozen=True)
class FragmentQuery:
    """Validated query expression; atom indices refer to this compiled order."""

    expression: str
    query_format: str
    topology: str
    query_id: str
    definition_version: str
    molecule: Any
    compiler_version: str = QUERY_COMPILER_VERSION

    def describe(self) -> dict[str, Any]:
        """Return portable query semantics without serializing the RDKit object."""
        return {"expression": self.expression, "query_format": self.query_format,
                "topology": self.topology, "query_id": self.query_id,
                "definition_version": self.definition_version,
                "compiler_version": self.compiler_version,
                "compiled_smarts": Chem.MolToSmarts(self.molecule),
                "atom_count": self.molecule.GetNumAtoms(), "specified_stereo_required": True,
                "automatic_relaxations": []}


def compile_fragment_query(
    query: str, query_format: str = "smiles", topology: str = "preserve_rings",
) -> FragmentQuery:
    """Compile a connected fragment, rejecting unsupported or contradictory queries."""
    policy = fragment_search_policy()
    if not isinstance(query, str) or not query.strip() or len(query) > policy["max_query_characters"]:
        raise ValueError("query must be nonempty text of at most 2000 characters")
    if query_format not in {"smiles", "smarts"} or topology not in {"preserve_rings", "subgraph"}:
        raise ValueError("Use smiles/smarts and preserve_rings/subgraph")
    query = query.strip()
    if query_format == "smarts" and any(token in query for token in ("$", "|", ">")):
        raise ValueError("Recursive, extended, and reaction SMARTS are unsupported; use a connected core")
    mol = Chem.MolFromSmiles(query) if query_format == "smiles" else compile_smarts(query, validate=True)
    if mol is None or not 0 < mol.GetNumAtoms() <= policy["max_query_atoms"]:
        raise ValueError("Invalid fragment or fragment exceeds 100 atoms")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Select one connected fragment; disconnected queries are unsupported")
    if any(atom.GetAtomMapNum() for atom in mol.GetAtoms()):
        raise ValueError("Remove atom-map labels from the query; returned query indices identify atoms")
    mol = Chem.Mol(mol)
    Chem.GetSymmSSSR(mol)
    if query_format == "smarts" and topology == "preserve_rings":
        for atom in mol.GetAtoms():
            description = atom.DescribeQuery()
            if (not atom.IsInRing() and (atom.GetIsAromatic() or "Ring" in description)) or "AtomOr" in description:
                raise ValueError("Partial/alternative ring SMARTS require topology='subgraph'")
        if any(any(line.strip() == "BondOr" for line in bond.DescribeQuery().splitlines()) for bond in mol.GetBonds()):
            raise ValueError("Alternative bond SMARTS require topology='subgraph'")
    if query_format == "smiles":
        if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
            raise ValueError("Radical SMILES queries are unsupported; select a closed-shell core")
        # Mol queries do not constrain neutral charge or aliphaticity by default.
        # Convert atoms to explicit query atoms while retaining query-order stereo.
        for atom in list(mol.GetAtoms()):
            # In a decorated SMARTS atom, H can mean hydrogen count rather than
            # element 1. Use an atomic number for explicit H / D / T atoms.
            symbol = ("#1" if atom.GetAtomicNum() == 1 else
                      atom.GetSymbol().lower() if atom.GetIsAromatic() else atom.GetSymbol())
            isotope = str(atom.GetIsotope()) if atom.GetIsotope() else ""
            hydrogen = f";H{atom.GetNumExplicitHs()}" if atom.GetNoImplicit() or atom.GetNumExplicitHs() else ""
            pattern = Chem.Mol(compile_smarts(f"[{isotope}{symbol}{hydrogen};{atom.GetFormalCharge():+d}]", validate=True))
            replacement = pattern.GetAtomWithIdx(0)
            replacement.SetChiralTag(atom.GetChiralTag())
            editable = Chem.RWMol(mol)
            editable.ReplaceAtom(atom.GetIdx(), replacement)
            mol = editable.GetMol()
        Chem.GetSymmSSSR(mol)
    identity = json.dumps([query, query_format, topology, policy["definition_version"], QUERY_COMPILER_VERSION])
    return FragmentQuery(query, query_format, topology,
                         "FQ1:" + hashlib.sha256(identity.encode()).hexdigest(),
                         policy["definition_version"], mol)


def fragment_embeddings(
    query: FragmentQuery, molecule: Any, *, maximum: int = 64,
) -> tuple[tuple[tuple[int, ...], ...], bool]:
    """Return query-ordered matches and an explicit enumeration truncation flag."""
    parameters = Chem.SubstructMatchParameters()
    parameters.useChirality = True
    parameters.uniquify = False  # Query edge assignments can differ on the same atom set.
    parameters.maxMatches = maximum + 1
    if query.topology == "preserve_rings":
        def ring_match(product: Any, match: Any) -> bool:
            for atom in query.molecule.GetAtoms():
                other = product.GetAtomWithIdx(match[atom.GetIdx()])
                if atom.IsInRing() != other.IsInRing():
                    return False
                query_ring_degree = sum(b.IsInRing() for b in atom.GetBonds())
                if query_ring_degree != sum(b.IsInRing() for b in other.GetBonds()):
                    return False
            return all(b.IsInRing() == product.GetBondBetweenAtoms(
                match[b.GetBeginAtomIdx()], match[b.GetEndAtomIdx()]).IsInRing()
                       for b in query.molecule.GetBonds())
        parameters.setExtraFinalCheck(ring_match)
    matches = molecule.GetSubstructMatches(query.molecule, parameters)
    return tuple(matches[:maximum]), len(matches) > maximum


@dataclass(frozen=True)
class FragmentTargetValidation:
    """Target membership under the exact corpus query semantics, not feasibility."""

    target_smiles: str
    query_id: str
    matches_target: bool
    definition_version: str
    compiler_version: str
    schema_version: str = "fragment_target_validation.v1"


def validate_fragment_target(query: FragmentQuery, target_smiles: str) -> FragmentTargetValidation:
    """Check the target with the same stereo and topology rules as corpus hits.

    A mismatch is evidence about the query, not absence from the corpus. No query
    relaxation or atom correspondence is inferred.
    """
    if not isinstance(target_smiles, str) or not target_smiles.strip():
        raise ValueError("Provide one connected target SMILES")
    canonical, _ = indexed_product(target_smiles)
    molecule = Chem.MolFromSmiles(canonical)
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("Provide one connected target SMILES")
    matches, _ = fragment_embeddings(query, molecule, maximum=1)
    return FragmentTargetValidation(canonical, query.query_id, bool(matches),
                                    query.definition_version, query.compiler_version)


def indexed_product(smiles: str) -> tuple[str, tuple[int, ...]]:
    """Normalize maps away, preserving canonical-index -> original-index identity."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Invalid product component")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    order = tuple(json.loads(mol.GetProp("_smilesAtomOutputOrder")))
    if len(order) != mol.GetNumAtoms():
        raise ValueError("Product atom-order projection is incomplete")
    return canonical, order


def _reference(atom: Any, side: str, component: int) -> dict[str, Any]:
    return asdict(ReactionAtomReference(
        side, component, atom.GetIdx(), atom.GetAtomMapNum() or None,
        atom.GetSymbol(), atom.GetFormalCharge(), atom.GetIsAromatic(),
        str(atom.GetHybridization()), "not_computed",
        str(atom.GetChiralTag()), atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None,
    ))


def project_fragment_evidence(reaction_smiles: str, observation: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Precompute graph changes using only already validated supplied correspondence.

Unmapped/reconstructed observations retain unresolved evidence in this initial
implementation. Independent graph checks prevent stale supplied-map claims.
"""
    parts = reaction_smiles.split(">")
    if len(parts) != 3:
        return []
    sides = [[Chem.MolFromSmiles(s) for s in parts[index].split(".") if s] for index in (0, 2)]
    reactants, products = sides
    projections = [{"component_index": i, "origins": {}, "changes": [],
                    "uncertain_atoms": list(range(m.GetNumAtoms())) if m else [],
                    "evidence_status": "unresolved", "warnings": list(observation.get("warnings") or ())}
                   for i, m in enumerate(products)]
    if not observation.get("valid") or observation.get("evidence_quality") != "validated_atom_mapping":
        return projections
    if observation.get("input_reaction_smiles") != reaction_smiles:
        return projections
    if observation.get("edit_hypotheses") or any(
        c.get("status") in {"ambiguous", "invalid"} for c in (observation.get("evidence_candidates") or ())
    ) or any("CONFLICT" in str(w).upper() for w in (observation.get("warnings") or ())):
        return projections
    maps = []
    for molecules in sides:
        lookup = {}
        for ci, mol in enumerate(molecules):
            if mol is None:
                return projections
            for atom in mol.GetAtoms():
                number = atom.GetAtomMapNum()
                if number:
                    if number in lookup:
                        return projections
                    lookup[number] = (ci, atom)
        maps.append(lookup)
    before, after = maps
    for number in before.keys() & after.keys():
        a, b = before[number][1], after[number][1]
        if (a.GetAtomicNum(), a.GetIsotope()) != (b.GetAtomicNum(), b.GetIsotope()):
            return projections
    for ci, mol in enumerate(products):
        entry = projections[ci]
        entry["evidence_status"] = "validated_supplied_mapping"
        uncertain = set()
        for atom in mol.GetAtoms():
            ai, number = atom.GetIdx(), atom.GetAtomMapNum()
            if not number or number not in before:
                uncertain.add(ai)
                continue
            rci, origin = before[number]
            entry["origins"][str(ai)] = {
                "reactant": _reference(origin, "reactant", rci),
                "product": _reference(atom, "product", ci),
            }
            state = lambda a: (a.GetFormalCharge(), a.GetTotalNumHs(), a.GetIsAromatic(),
                               a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else None)
            if state(origin) != state(atom):
                entry["changes"].append({"kind": "atom_state", "product_atoms": [ai],
                                         "before": list(state(origin)), "after": list(state(atom))})
            if any(not n.GetAtomMapNum() for n in origin.GetNeighbors()) or any(
                not n.GetAtomMapNum() or n.GetAtomMapNum() not in before for n in atom.GetNeighbors()
            ):
                uncertain.add(ai)
        # Compare mapped incident edges; absent endpoints are boundary losses, not no_bond formation.
        keys = set()
        for molecules in sides:
            for component in molecules:
                for bond in component.GetBonds():
                    nums = (bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum())
                    if all(nums):
                        keys.add(tuple(sorted(nums)))
        for a, b in sorted(keys):
            endpoints = [after[n][1].GetIdx() if n in after and after[n][0] == ci else None for n in (a, b)]
            if all(i is None for i in endpoints):
                continue
            if a not in before or b not in before:
                uncertain.update(i for i in endpoints if i is not None)
                continue
            def bond_state(lookup: dict, molecules: list) -> str | None:
                if a not in lookup or b not in lookup or lookup[a][0] != lookup[b][0]:
                    return None
                bond = molecules[lookup[a][0]].GetBondBetweenAtoms(lookup[a][1].GetIdx(), lookup[b][1].GetIdx())
                return str(bond.GetBondType()) if bond else None
            old, new = bond_state(before, reactants), bond_state(after, products)
            if old != new:
                entry["changes"].append({"kind": "formed" if old is None else "broken" if new is None else "order_changed",
                                         "product_atoms": endpoints, "atom_maps": [a, b],
                                         "before": {"state_kind": "bond" if old else "no_bond", "order": old},
                                         "after": {"state_kind": "bond" if new else "endpoint_absent" if None in endpoints else "no_bond", "order": new}})
            elif old is not None:
                old_bond = reactants[before[a][0]].GetBondBetweenAtoms(before[a][1].GetIdx(), before[b][1].GetIdx())
                new_bond = products[after[a][0]].GetBondBetweenAtoms(after[a][1].GetIdx(), after[b][1].GetIdx())
                if old_bond.GetStereo() != new_bond.GetStereo():
                    entry["changes"].append({"kind": "bond_stereo", "product_atoms": endpoints,
                                             "before": str(old_bond.GetStereo()), "after": str(new_bond.GetStereo())})
        entry["uncertain_atoms"] = sorted(uncertain)
    return projections


def classify_fragment_embedding(
    query: FragmentQuery, embedding: tuple[int, ...], original_order: tuple[int, ...],
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    """Explain local changes for one embedding without asserting route feasibility."""
    matched = tuple(original_order[i] for i in embedding)
    selected = set(matched)
    query_edges = {frozenset((matched[b.GetBeginAtomIdx()], matched[b.GetEndAtomIdx()]))
                   for b in query.molecule.GetBonds()}
    relationships, witnesses = set(), []
    for change in evidence.get("changes", ()):
        endpoints = change["product_atoms"]
        inside = [i in selected for i in endpoints]
        if not any(inside):
            continue
        relation = "modified"
        if len(endpoints) == 2:
            if not all(inside):
                relation = "boundary_changed"
            elif change["kind"] == "formed" and frozenset(endpoints) in query_edges:
                relation = "constructed"
        relationships.add(relation)
        witnesses.append({**change, "relationship": relation})
    unresolved = (evidence.get("evidence_status") != "validated_supplied_mapping"
                  or bool(selected.intersection(evidence.get("uncertain_atoms", ()))))
    if unresolved:
        relationships.add("unresolved")
    if not relationships:
        relationships.add("carried_through")
    return {"query_to_original_product_atoms": list(matched),
            "relationships": sorted(relationships), "witnesses": witnesses,
            "atom_origins": [evidence.get("origins", {}).get(str(i)) for i in matched],
            "evidence_status": evidence.get("evidence_status", "unresolved")}
