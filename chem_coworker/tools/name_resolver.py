"""
Chemical name → SMILES resolver tools for ChemCoworker.

Two ToolPlugin entries:

  resolve_name_to_smiles  — convert a single chemical name to SMILES
  resolve_names_to_smiles — convert a comma-separated list of names (batch)

Both use a three-backend cascade: PubChem → OPSIN → CACTUS, with in-memory
caching and RDKit validation. These tools are especially useful for retrosynthesis
queries where the user types a molecule name rather than a SMILES string.

Typical usage in retrosynthesis workflow:
  User: "How do I make tamoxifen?"
  → Classifier finds no SMILES in query
  → LLM plan includes resolve_name_to_smiles as the first tool group
  → Resolved SMILES is passed to inspect_target / generate_disconnections
"""
from __future__ import annotations

from typing import Any, Dict, List

from ._helpers import _error, _success
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Tool A: resolve_name_to_smiles
# ---------------------------------------------------------------------------

def _resolve_name_to_smiles(name: str) -> Dict[str, Any]:
    """Convert a chemical name to its canonical SMILES using PubChem/OPSIN/CACTUS.

    Tries three external services in cascade:
      1. PubChem  — handles IUPAC names, trivial names, trade names, CAS numbers
      2. OPSIN    — specialized in systematic IUPAC nomenclature
      3. CACTUS   — NCI broad-coverage fallback

    The returned SMILES is validated by RDKit and canonicalized. Also returns
    the CID (PubChem Compound ID) and molecular formula when available.

    Args:
        name: Chemical name to look up (IUPAC, trivial, CAS, or trade name).
              Examples: "tamoxifen", "4-methylbiphenyl", "50-78-2" (aspirin CAS)

    Returns:
        dict with smiles, canonical_smiles, source, name_used, and optionally
        pubchem_cid and molecular_formula.
    """
    if not name or not name.strip():
        return _error("Name cannot be empty")

    name = name.strip()

    try:
        from chemtools.util.name_to_smiles import resolve_name
        smiles, source = resolve_name(name)

        if not smiles:
            return _error(
                f"Could not resolve '{name}' to a SMILES string. "
                f"Tried PubChem, OPSIN, and CACTUS. "
                f"Check spelling, try a synonym, or provide the SMILES directly.",
                {"name_tried": name, "source": "not_found"},
            )

        # Optionally enrich with PubChem CID + formula
        pubchem_cid = None
        mol_formula = None
        if source == "pubchem":
            try:
                pubchem_cid, mol_formula = _pubchem_enrich(name)
            except Exception:
                pass

        result: Dict[str, Any] = {
            "name": name,
            "smiles": smiles,
            "source": source,
            "validated": True,
        }
        if pubchem_cid:
            result["pubchem_cid"] = pubchem_cid
        if mol_formula:
            result["molecular_formula"] = mol_formula

        return _success(result)

    except Exception as exc:
        return _error(f"Name resolution failed: {exc}")


def _pubchem_enrich(name: str):
    """Fetch PubChem CID and molecular formula as enrichment data."""
    import json
    from urllib.parse import quote_plus

    try:
        import requests
        encoded = quote_plus(name)
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{encoded}/property/MolecularFormula,CID/JSON"
        )
        resp = requests.get(url, timeout=8, headers={"User-Agent": "ChemCoworker/1.0"})
        if resp.status_code == 200:
            data = json.loads(resp.text)
            props = data.get("PropertyTable", {}).get("Properties", [])
            if props:
                return props[0].get("CID"), props[0].get("MolecularFormula")
    except Exception:
        pass
    return None, None


resolve_name_to_smiles_tool = ToolPlugin(
    name="resolve_name_to_smiles",
    category="name_resolver",
    description=(
        "Convert a chemical name (IUPAC, trivial, trade name, or CAS number) to its "
        "SMILES string using PubChem → OPSIN → CACTUS cascade with RDKit validation. "
        "Use this FIRST when a retrosynthesis or analysis query contains a molecule name "
        "rather than a SMILES string. Examples: 'tamoxifen', '4-methylbiphenyl', 'aspirin'."
    ),
    prerequisites=[],
    fn=_resolve_name_to_smiles,
)


# ---------------------------------------------------------------------------
# Tool B: resolve_names_to_smiles (batch)
# ---------------------------------------------------------------------------

def _resolve_names_to_smiles(names: str) -> Dict[str, Any]:
    """Convert multiple chemical names to SMILES in a single call.

    Accepts a comma-separated string of chemical names. Resolves each
    sequentially using the PubChem → OPSIN → CACTUS cascade. Useful when
    a query mentions several reagents or substrates by name.

    Args:
        names: Comma-separated chemical names.
               Example: "phenylboronic acid, 4-bromotoluene, palladium acetate"

    Returns:
        dict with resolved list (name, smiles, source, resolved: bool),
        total_count, resolved_count, failed (list of names not found).
    """
    if not names or not names.strip():
        return _error("Names cannot be empty")

    # Parse comma-separated input
    name_list = [n.strip() for n in names.split(",") if n.strip()]
    if not name_list:
        return _error("No valid names found after splitting on commas")

    try:
        from chemtools.util.name_to_smiles import resolve_names
        results_map = resolve_names(name_list)

        resolved: List[Dict[str, Any]] = []
        failed: List[str] = []

        for name in name_list:
            smiles, source = results_map.get(name, (None, "not_found"))
            if smiles:
                resolved.append({
                    "name": name,
                    "smiles": smiles,
                    "source": source,
                    "resolved": True,
                })
            else:
                resolved.append({
                    "name": name,
                    "smiles": None,
                    "source": "not_found",
                    "resolved": False,
                })
                failed.append(name)

        return _success({
            "total_count": len(name_list),
            "resolved_count": len(name_list) - len(failed),
            "failed_count": len(failed),
            "results": resolved,
            "failed": failed,
            "message": (
                f"Resolved {len(name_list) - len(failed)}/{len(name_list)} names. "
                + (f"Could not resolve: {', '.join(failed)}" if failed else "All names resolved.")
            ),
        })

    except Exception as exc:
        return _error(f"Batch name resolution failed: {exc}")


resolve_names_to_smiles_tool = ToolPlugin(
    name="resolve_names_to_smiles",
    category="name_resolver",
    description=(
        "Convert multiple chemical names (comma-separated) to SMILES strings in batch. "
        "Uses PubChem → OPSIN → CACTUS cascade for each name. Useful when a query "
        "mentions several reagents, substrates, or intermediates by name. "
        "Example input: 'phenylboronic acid, 4-bromotoluene, Pd(OAc)2'"
    ),
    prerequisites=[],
    fn=_resolve_names_to_smiles,
)


# ---------------------------------------------------------------------------
# Tool C: resolve_smiles_to_name (single SMILES → names)
# ---------------------------------------------------------------------------

def _resolve_smiles_to_name(smiles: str) -> Dict[str, Any]:
    """Convert a SMILES string to its chemical names using PubChem → CACTUS.

    Canonicalizes the input SMILES with RDKit before lookup, then queries:
      1. PubChem — returns IUPAC name, preferred/common name, synonyms,
                   molecular formula, and PubChem CID
      2. CACTUS  — fallback for names and IUPAC when PubChem misses

    Useful for labelling precursor SMILES from retrosynthetic disconnections,
    or identifying unknown molecules by structure.

    Args:
        smiles: Molecule SMILES. Reaction SMILES (with >>) are stripped to
                the first component automatically.
                Example: "CC(=O)Oc1ccccc1C(=O)O"

    Returns:
        dict with iupac_name, preferred_name, synonyms (up to 8),
        molecular_formula, pubchem_cid, source, and canonical_smiles.
    """
    if not smiles or not smiles.strip():
        return _error("SMILES cannot be empty")

    smiles = smiles.strip()

    try:
        from chemtools.util.smiles_to_name import resolve_smiles
        result = resolve_smiles(smiles)

        if result.get("source") == "not_found":
            return _error(
                f"Could not identify names for SMILES '{smiles}'. "
                "The compound may not be in PubChem or CACTUS databases. "
                "Try canonicalizing the SMILES or check that it is valid.",
                {"smiles": smiles, "canonical_smiles": result.get("canonical_smiles", smiles)},
            )

        return _success(result)

    except Exception as exc:
        return _error(f"SMILES-to-name resolution failed: {exc}")


resolve_smiles_to_name_tool = ToolPlugin(
    name="resolve_smiles_to_name",
    category="name_resolver",
    description=(
        "Convert a SMILES string to its chemical names (IUPAC name, preferred/common name, "
        "synonyms, molecular formula, PubChem CID) using PubChem → CACTUS cascade. "
        "Useful for labelling retrosynthetic precursors, identifying unknown structures, "
        "or annotating reaction SMILES with human-readable names."
    ),
    prerequisites=[],
    fn=_resolve_smiles_to_name,
)


# ---------------------------------------------------------------------------
# Tool D: resolve_smiles_to_names_batch (multiple SMILES → names)
# ---------------------------------------------------------------------------

def _resolve_smiles_to_names_batch(smiles_list: str) -> Dict[str, Any]:
    """Convert multiple SMILES strings to chemical names in a single call.

    Accepts a comma-separated string of SMILES. Resolves each via the
    PubChem → CACTUS cascade. Useful for annotating all precursors returned
    by generate_disconnections in one step.

    Args:
        smiles_list: Comma-separated SMILES strings.
                     Example: "Brc1ccccc1, OB(O)c1ccccc1, c1ccc(-c2ccccc2)cc1"

    Returns:
        dict with results list (smiles, iupac_name, preferred_name, source,
        resolved: bool), total_count, resolved_count, failed.
    """
    if not smiles_list or not smiles_list.strip():
        return _error("SMILES list cannot be empty")

    items = [s.strip() for s in smiles_list.split(",") if s.strip()]
    if not items:
        return _error("No valid SMILES found after splitting on commas")

    try:
        from chemtools.util.smiles_to_name import resolve_smiles
        results = []
        failed = []

        for smi in items:
            result = resolve_smiles(smi)
            resolved = result.get("source") != "not_found"
            entry: Dict[str, Any] = {
                "smiles": smi,
                "canonical_smiles": result.get("canonical_smiles", smi),
                "iupac_name": result.get("iupac_name", ""),
                "preferred_name": result.get("preferred_name", ""),
                "molecular_formula": result.get("molecular_formula", ""),
                "pubchem_cid": result.get("pubchem_cid"),
                "source": result.get("source", "not_found"),
                "resolved": resolved,
            }
            results.append(entry)
            if not resolved:
                failed.append(smi)

        return _success({
            "total_count": len(items),
            "resolved_count": len(items) - len(failed),
            "failed_count": len(failed),
            "results": results,
            "failed": failed,
            "message": (
                f"Named {len(items) - len(failed)}/{len(items)} SMILES. "
                + (f"Could not identify: {', '.join(failed)}" if failed else "All SMILES identified.")
            ),
        })

    except Exception as exc:
        return _error(f"Batch SMILES-to-name resolution failed: {exc}")


resolve_smiles_to_names_batch_tool = ToolPlugin(
    name="resolve_smiles_to_names_batch",
    category="name_resolver",
    description=(
        "Convert multiple SMILES strings (comma-separated) to chemical names in batch. "
        "Returns IUPAC name, preferred name, molecular formula, and PubChem CID for each. "
        "Use after generate_disconnections to annotate precursor SMILES with human-readable names."
    ),
    prerequisites=[],
    fn=_resolve_smiles_to_names_batch,
)


# ---------------------------------------------------------------------------
# Module export
# ---------------------------------------------------------------------------

NAME_RESOLVER_TOOLS = [
    resolve_name_to_smiles_tool,
    resolve_names_to_smiles_tool,
    resolve_smiles_to_name_tool,
    resolve_smiles_to_names_batch_tool,
]
