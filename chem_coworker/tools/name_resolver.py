"""
Chemical identifier resolver tools for ChemCoworker.

Two consolidated ToolPlugin entries:

  resolve_to_smiles  — convert any identifier (name, CAS, InChI, trade name) → SMILES
  smiles_to_info     — convert SMILES → name + CAS + formula + PubChem CID

Both use a three-backend cascade: PubChem → OPSIN → CACTUS, with in-memory
caching and RDKit validation. These tools are especially useful for retrosynthesis
queries where the user types a molecule name or CAS number rather than a SMILES string.

Typical usage in retrosynthesis workflow:
  User: "How do I make tamoxifen?"
  → Classifier finds no SMILES in query
  → LLM plan includes resolve_to_smiles as the first tool group
  → Resolved SMILES is passed to inspect_target / generate_disconnections
"""
from __future__ import annotations

import re
from typing import Any, Dict, List

from ._helpers import _error, _success
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Internal helpers (kept for reuse; not exported as ToolPlugins)
# ---------------------------------------------------------------------------

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


# CAS pattern: 1-7 digits, dash, 2 digits, dash, 1 digit checksum
_CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")


# ---------------------------------------------------------------------------
# Tool 1: resolve_to_smiles
# Replaces: resolve_name_to_smiles, resolve_names_to_smiles,
#           resolve_cas (single identifier path)
# ---------------------------------------------------------------------------

def _resolve_to_smiles(identifier: str) -> Dict[str, Any]:
    """Convert any chemical identifier to its canonical SMILES string.

    Accepts any of:
      • Chemical names — IUPAC, trivial, trade names (e.g. "tamoxifen", "aspirin")
      • CAS Registry Numbers (e.g. "50-78-2", "58-08-2")
      • InChI strings (e.g. "InChI=1S/C9H8O4/...")
      • Comma-separated list of the above

    Uses a PubChem → OPSIN → CACTUS cascade with RDKit validation.
    This is the FIRST tool to call when the user provides a molecule name or
    CAS number instead of a SMILES string.

    Args:
        identifier: Chemical name, CAS number, InChI, or comma-separated list.
                    Examples: "tamoxifen"
                              "50-78-2"
                              "phenylboronic acid, 4-bromotoluene"

    Returns:
        dict with smiles (primary resolved SMILES), source, name_used,
        and optionally pubchem_cid, molecular_formula. For multi-input:
        also includes 'results' list with per-identifier outcomes.
    """
    if not identifier or not identifier.strip():
        return _error("Identifier cannot be empty")

    identifier = identifier.strip()

    # Multi-identifier: comma-separated list
    parts = [p.strip() for p in identifier.split(",") if p.strip()]
    if len(parts) > 1:
        return _resolve_batch(parts)

    # Single identifier — auto-detect type
    return _resolve_single(identifier)


def _resolve_single(identifier: str) -> Dict[str, Any]:
    """Resolve a single identifier (name, CAS, or InChI) to SMILES."""
    try:
        # CAS number path
        if _CAS_RE.match(identifier):
            try:
                from chemtools.util.cas_resolver import resolve_cas as _cas_fn
                result = _cas_fn(identifier)
                smiles = result.get("smiles", "")
                if smiles:
                    return _success({
                        "name": identifier,
                        "smiles": smiles,
                        "source": result.get("source", "pubchem_cas"),
                        "validated": True,
                        "name_used": result.get("name", identifier),
                        **{k: result[k] for k in ("molecular_formula", "pubchem_cid", "iupac_name")
                           if k in result and result[k]},
                    })
            except Exception:
                pass  # fall through to name resolver

        # InChI path — PubChem handles InChI as name query
        # Name path (also handles InChI via PubChem)
        from chemtools.util.name_to_smiles import resolve_name
        smiles, source = resolve_name(identifier)

        if not smiles:
            return _error(
                f"Could not resolve '{identifier}' to a SMILES string. "
                f"Tried PubChem, OPSIN, and CACTUS. "
                f"Check spelling, try a synonym, or provide the SMILES directly.",
                {"name_tried": identifier, "source": "not_found"},
            )

        # Optionally enrich with PubChem CID + formula
        pubchem_cid = None
        mol_formula = None
        if source == "pubchem":
            try:
                pubchem_cid, mol_formula = _pubchem_enrich(identifier)
            except Exception:
                pass

        result: Dict[str, Any] = {
            "name": identifier,
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
        return _error(f"Identifier resolution failed: {exc}")


def _resolve_batch(parts: List[str]) -> Dict[str, Any]:
    """Resolve a list of identifiers; return summary + per-identifier results."""
    resolved_list: List[Dict[str, Any]] = []
    failed: List[str] = []
    primary_smiles = None

    for part in parts:
        r = _resolve_single(part)
        if r.get("success") and r.get("smiles"):
            entry = {"identifier": part, "smiles": r["smiles"],
                     "source": r.get("source", ""), "resolved": True}
            if primary_smiles is None:
                primary_smiles = r["smiles"]
        else:
            entry = {"identifier": part, "smiles": None, "source": "not_found", "resolved": False}
            failed.append(part)
        resolved_list.append(entry)

    return _success({
        "smiles": primary_smiles or "",
        "total_count": len(parts),
        "resolved_count": len(parts) - len(failed),
        "results": resolved_list,
        "failed": failed,
        "message": (
            f"Resolved {len(parts) - len(failed)}/{len(parts)} identifiers. "
            + (f"Could not resolve: {', '.join(failed)}" if failed else "All resolved.")
        ),
    })


resolve_to_smiles_tool = ToolPlugin(
    name="resolve_to_smiles",
    category="name_resolver",
    description=(
        "Convert any chemical identifier to SMILES: accepts chemical names (IUPAC, trivial, "
        "trade names), CAS Registry Numbers (e.g. '50-78-2'), InChI strings, or a "
        "comma-separated list of any of the above. Uses PubChem → OPSIN → CACTUS cascade "
        "with RDKit validation. Call this FIRST when the query contains a molecule name or "
        "CAS number rather than a SMILES. Args: identifier (str)."
    ),
    prerequisites=[],
    fn=_resolve_to_smiles,
    provides=["resolved_smiles", "smiles"],
)


# ---------------------------------------------------------------------------
# Tool 2: smiles_to_info
# Replaces: resolve_smiles_to_name, resolve_smiles_to_names_batch, smiles_to_cas
# ---------------------------------------------------------------------------

def _smiles_to_info(smiles: str) -> Dict[str, Any]:
    """Look up chemical identifiers and names for a SMILES string.

    Returns IUPAC name, preferred/common name, synonyms, CAS number(s),
    molecular formula, and PubChem CID. Accepts a single SMILES or a
    comma-separated list of SMILES strings.

    Useful for:
      • Labelling retrosynthetic precursors with human-readable names
      • Finding CAS numbers for ordering / safety data lookup
      • Identifying unknown structures by SMILES

    Args:
        smiles: Molecule SMILES, or comma-separated SMILES strings.
                Reaction SMILES (with >>) are stripped to the first component.
                Example: "CC(=O)Oc1ccccc1C(=O)O"
                Example: "Brc1ccccc1, OB(O)c1ccccc1"

    Returns:
        dict with iupac_name, preferred_name, synonyms (up to 8),
        molecular_formula, cas_numbers, primary_cas, pubchem_cid,
        canonical_smiles, and source.
        For multi-SMILES input: also includes 'results' list.
    """
    if not smiles or not smiles.strip():
        return _error("SMILES cannot be empty")

    smiles = smiles.strip()

    # Multi-SMILES: comma-separated list
    parts = [p.strip() for p in smiles.split(",") if p.strip()]
    if len(parts) > 1:
        return _smiles_to_info_batch(parts)

    return _smiles_to_info_single(smiles)


def _smiles_to_info_single(smiles: str) -> Dict[str, Any]:
    """Look up name + CAS for a single SMILES."""
    try:
        from chemtools.util.smiles_to_name import resolve_smiles

        result = resolve_smiles(smiles)
        if result.get("source") == "not_found":
            return _error(
                f"Could not identify names for SMILES '{smiles}'. "
                "The compound may not be in PubChem or CACTUS databases.",
                {"smiles": smiles, "canonical_smiles": result.get("canonical_smiles", smiles)},
            )

        # Enrich with CAS numbers
        cas_numbers: List[str] = []
        primary_cas = ""
        try:
            from chemtools.util.cas_resolver import smiles_to_cas as _stoc
            cas_numbers = _stoc(smiles) or []
            primary_cas = cas_numbers[0] if cas_numbers else ""
        except Exception:
            pass

        enriched = dict(result)
        enriched["cas_numbers"] = cas_numbers
        enriched["primary_cas"] = primary_cas
        return _success(enriched)

    except Exception as exc:
        return _error(f"SMILES info lookup failed: {exc}")


def _smiles_to_info_batch(smiles_list: List[str]) -> Dict[str, Any]:
    """Look up name/CAS for each SMILES in the list."""
    results = []
    failed = []

    for smi in smiles_list:
        r = _smiles_to_info_single(smi)
        resolved = r.get("success", False)
        entry: Dict[str, Any] = {
            "smiles": smi,
            "canonical_smiles": r.get("canonical_smiles", smi),
            "iupac_name": r.get("iupac_name", ""),
            "preferred_name": r.get("preferred_name", ""),
            "molecular_formula": r.get("molecular_formula", ""),
            "cas_numbers": r.get("cas_numbers", []),
            "primary_cas": r.get("primary_cas", ""),
            "pubchem_cid": r.get("pubchem_cid"),
            "source": r.get("source", "not_found"),
            "resolved": resolved,
        }
        results.append(entry)
        if not resolved:
            failed.append(smi)

    return _success({
        "total_count": len(smiles_list),
        "resolved_count": len(smiles_list) - len(failed),
        "results": results,
        "failed": failed,
        "message": (
            f"Named {len(smiles_list) - len(failed)}/{len(smiles_list)} SMILES. "
            + (f"Could not identify: {', '.join(failed)}" if failed else "All SMILES identified.")
        ),
    })


smiles_to_info_tool = ToolPlugin(
    name="smiles_to_info",
    category="name_resolver",
    description=(
        "Look up chemical names, CAS number, molecular formula, and PubChem CID for a SMILES "
        "string. Accepts a single SMILES or comma-separated list. Uses PubChem → CACTUS. "
        "Useful for labelling retrosynthetic precursors, annotating reaction SMILES with "
        "human-readable names, or finding CAS numbers for purchasing. Args: smiles (str)."
    ),
    prerequisites=[],
    fn=_smiles_to_info,
)


# ---------------------------------------------------------------------------
# Module export — 2 consolidated tools replace the original 6
# ---------------------------------------------------------------------------

NAME_RESOLVER_TOOLS = [
    resolve_to_smiles_tool,
    smiles_to_info_tool,
]
