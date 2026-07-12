"""
High-level helpers for adding reagents to the taxonomy registry.

This module centralizes the logic that was previously scattered across the
CLI and UI so automated clients (LangChain agents, GUIs, tests) can add
reagents through a single, documented API.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from .registry_store import (
    RegistryStore,
    build_registry_entry,
    infer_abbreviations,
)
from .taxonomy_store import RoleHeuristics, TaxonomyStore
from .taxonomy_utils import (
    DEFAULT_RESOLVER_TIMEOUT,
    dedupe_synonyms,
    normalize_cas,
    resolve_identity_from_cas,
    tokenize_all,
)
import re

DEFAULT_FLAT_REGISTRY_DIR = Path("condition_registry/definitions")


class ReagentAdditionError(RuntimeError):
    """Raised when automatic reagent addition fails."""


def _coerce_synonyms(values: Optional[Sequence[str]]) -> List[str]:
    if not values:
        return []
    cleaned: List[str] = []
    for value in values:
        if not value:
            continue
        text = str(value).strip()
        if text and text not in cleaned:
            cleaned.append(text)
    return cleaned


def _coerce_numeric(value: Optional[str | int | float]) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value).strip()
    if not text:
        return None
    match = re.search(r"[-+]?\d+(?:\.\d+)?", text)
    if not match:
        return None
    return match.group(0)


def _coerce_float(value: Optional[str | int | float]) -> Optional[float]:
    text = _coerce_numeric(value)
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _coerce_temperature_c(value: Optional[str | int | float]) -> Optional[str]:
    """Coerce a temperature value to Celsius as a numeric string."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value).strip()
    if not text:
        return None
    number = _coerce_float(text)
    if number is None:
        return None
    lower = text.lower()
    if "fahrenheit" in lower or "°f" in lower or "deg f" in lower or "degree f" in lower:
        number = (number - 32.0) * 5.0 / 9.0
    return f"{number:.1f}"


def _coerce_temperature_float(value: Optional[str | int | float]) -> Optional[float]:
    text = _coerce_temperature_c(value)
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _guess_reagent_type(mp: Optional[float], bp: Optional[float]) -> str:
    """Best-effort phase guess at room temperature (25C)."""
    if bp is not None and bp <= 25.0:
        return "gas"
    if mp is not None and mp <= 25.0:
        return "liquid"
    if bp is not None and bp <= 150.0 and (mp is None or mp <= 25.0):
        return "liquid"
    return "powder"


def _resolve_identity(cas: str, *, enabled: bool, timeout: float) -> Dict[str, Any]:
    if not enabled:
        return {}
    try:
        resolved = resolve_identity_from_cas(cas, timeout=timeout)
    except Exception:
        return {}
    return resolved or {}


def _load_taxonomy(taxonomy_dir: Optional[Path]) -> TaxonomyStore:
    """Load reagent taxonomy v2 for family/role inference."""
    target = taxonomy_dir.expanduser() if taxonomy_dir is not None else None
    try:
        return TaxonomyStore(target)
    except FileNotFoundError as exc:
        raise ReagentAdditionError(f"Unable to load reagent taxonomy v2: {exc}") from exc


def _load_registry(registry_dir: Optional[Path]) -> RegistryStore:
    """Load the flattened reagent registry used for storage."""
    if registry_dir is not None:
        return RegistryStore(registry_dir.expanduser())
    return RegistryStore(DEFAULT_FLAT_REGISTRY_DIR)


def add_reagent_entry(
    *,
    cas: str,
    taxonomy_dir: Optional[Path | str] = None,
    registry_dir: Optional[Path | str] = None,
    name: Optional[str] = None,
    formula: Optional[str] = None,
    molecular_weight: Optional[str] = None,
    density: Optional[str] = None,
    boiling_point: Optional[str] = None,
    melting_point: Optional[str] = None,
    reagent_type: Optional[str] = None,
    synonyms: Optional[Sequence[str]] = None,
    abbreviation: Optional[str] = None,
    tag: Optional[str] = None,
    role: Optional[str] = None,
    family_id: Optional[str] = None,
    smiles: Optional[str] = None,
    allow_default_family: bool = False,
    dry_run: bool = False,
    auto_resolve: bool = True,
    resolver_timeout: float = DEFAULT_RESOLVER_TIMEOUT,
) -> Dict[str, Any]:
    """
    Add a reagent entry to the flattened registry using v2 taxonomy for classification.

    Args:
        cas: CAS registry identifier (required).
        taxonomy_dir: Optional path to the reagent taxonomy v2 directory. When not supplied,
            the packaged v2 taxonomy is used.
        registry_dir: Optional path to the flattened reagent registry directory. When not supplied,
            the standalone condition_registry definitions directory is used.
        name: Preferred reagent name. If omitted, an online resolver is used when permitted.
        formula: Molecular formula override (e.g., C6H6).
        molecular_weight: Molecular weight override (e.g., 78.11).
        density: Density override (e.g., 0.867 g/mL).
        boiling_point: Boiling point override (e.g., 110.6 C).
        melting_point: Melting point override (e.g., -95 C).
        reagent_type: Physical form hint (e.g., powder, liquid, gas).
        synonyms: Optional list of synonyms/aliases.
        abbreviation: Explicit abbreviation override. Defaults to the reagent name.
        tag: Optional short tag for additional notes.
        role: Reagent role (ligand, base, metal_catalyst, ...). When omitted, the role
            heuristics attempt to infer it.
        family_id: Target family identifier. When omitted, heuristics attempt to infer it
            (optionally falling back to default families when ``allow_default_family`` is True).
        smiles: Optional SMILES annotation stored with the entry.
        allow_default_family: Allow default family fallback when inference fails.
        dry_run: When True, do not update files; return a preview instead.
        auto_resolve: Whether to resolve identity from CAS when data is missing.
        resolver_timeout: Timeout for the resolver (seconds).

    Returns:
        JSON-serialisable dictionary describing the action performed.

    Raises:
        ReagentAdditionError: When mandatory information is missing or inconsistent.
    """
    if not cas:
        raise ReagentAdditionError("CAS number is required.")

    normalized_cas = normalize_cas(cas)
    taxonomy_path = Path(taxonomy_dir).expanduser() if taxonomy_dir else None
    registry_path = Path(registry_dir).expanduser() if registry_dir else None

    taxonomy = _load_taxonomy(taxonomy_path)
    registry = _load_registry(registry_path)
    heuristics = RoleHeuristics(taxonomy)

    resolved_identity = _resolve_identity(normalized_cas, enabled=auto_resolve, timeout=resolver_timeout)
    resolved_name = resolved_identity.get("name")
    resolved_synonyms = _coerce_synonyms(resolved_identity.get("synonyms"))
    resolved_smiles = resolved_identity.get("smiles")
    resolved_formula = resolved_identity.get("molecular_formula")
    resolved_mw = resolved_identity.get("molecular_weight")
    resolved_density = resolved_identity.get("density")
    resolved_bp = resolved_identity.get("boiling_point")
    resolved_mp = resolved_identity.get("melting_point")
    auto_resolve_source = resolved_identity.get("source")

    final_name = (name or resolved_name or "").strip()
    if not final_name:
        raise ReagentAdditionError(
            "Unable to determine reagent name. Provide a name or enable CAS identity resolution."
        )

    provided_syns = _coerce_synonyms(synonyms)
    abbr = (abbreviation or final_name).strip()
    combined_synonyms = dedupe_synonyms([final_name, abbr, *provided_syns, *resolved_synonyms])
    token_list = tokenize_all(combined_synonyms)
    token_set = {tok for tok in token_list if tok}

    role_reason: Optional[str] = None
    family_tokens: Optional[List[str]] = None
    default_rejection: Optional[str] = None
    family_candidates: List[Dict[str, Any]] = []
    used_default = False

    if family_id:
        family_role = taxonomy.role_for_family(family_id)
        if not family_role:
            raise ReagentAdditionError(f"Unknown family '{family_id}'.")
        if role and role != family_role:
            raise ReagentAdditionError(
                f"Provided role '{role}' conflicts with family '{family_id}' (expected role '{family_role}')."
            )
        role = family_role

    if not role:
        inferred_role = heuristics.infer_role(final_name, combined_synonyms)
        if inferred_role:
            role, pattern = inferred_role
            role_reason = pattern

    if not role:
        raise ReagentAdditionError(
            "Unable to determine reagent role. Provide it explicitly or supply additional synonyms."
        )
    valid_roles = set(taxonomy.role_data.keys())
    if role not in valid_roles:
        raise ReagentAdditionError(f"Unsupported role '{role}'.")

    if not family_id:
        family_inference = heuristics.infer_family(final_name, combined_synonyms)
        if family_inference:
            inferred_role, inferred_family, matches = family_inference
            if inferred_role == role:
                family_id = inferred_family
                family_tokens = matches
            else:
                family_tokens = matches  # provide hint even if mismatch

    if not family_id:
        default_family = heuristics.default_family_for_role(role) if allow_default_family else None
        if default_family:
            use_default = False
            if token_set and taxonomy.family_token_overlap(role, default_family, token_set):
                use_default = True
            if use_default:
                family_id = default_family
                used_default = True
            else:
                default_rejection = (
                    f"default family '{default_family}' rejected: no token overlap with "
                    f"input tokens ({', '.join(sorted(token_set)[:6]) or 'none'})."
                )

    if not family_id:
        raise ReagentAdditionError(
            "Unable to determine family. Provide it explicitly or enable default fallback."
        )

    suggestions = taxonomy.suggest_families(role, token_set, limit=5)
    if suggestions:
        family_candidates = suggestions

    family_entry = registry.family_entry(role, family_id)
    if not family_entry:
        family_entry = taxonomy.family_data(role, family_id)
    if not family_entry:
        raise ReagentAdditionError(f"Family '{family_id}' not found for role '{role}'.")

    existing = registry.find_by_cas(normalized_cas, role=role)
    if existing:
        existing_role, existing_family, data = existing
        result = {
            "status": "exists",
            "cas": normalized_cas,
            "name": data.get("name"),
            "role": existing_role,
            "family_id": existing_family,
            "registry_type": "flat",
            "registry_file": str(registry.file_for_role(existing_role)),
        }
        if family_tokens:
            result["family_tokens"] = family_tokens
        if role_reason:
            result["role_reason"] = role_reason
        return result

    final_smiles = smiles or resolved_smiles
    abbreviations = infer_abbreviations(final_name, combined_synonyms)
    abbr_value = (abbreviation or (abbreviations[0] if abbreviations else "")).strip()
    entry = build_registry_entry(
        name=final_name,
        abbreviation=abbr_value,
        cas=normalized_cas,
        smile=final_smiles,
        role=role,
        family_id=family_id,
        tag=tag,
        family_entry=family_entry,
        synonyms=list(combined_synonyms),
    )
    properties: Dict[str, Any] = {}
    if formula or resolved_formula:
        properties["formula"] = str(formula or resolved_formula).strip()
    if molecular_weight is not None or resolved_mw is not None:
        properties["mw"] = str(molecular_weight or resolved_mw).strip()
    if density or resolved_density:
        coerced = _coerce_numeric(density or resolved_density)
        if coerced is not None:
            properties["density"] = coerced
    if boiling_point or resolved_bp:
        coerced = _coerce_temperature_c(boiling_point or resolved_bp)
        if coerced is not None:
            properties["bp"] = coerced
    if melting_point or resolved_mp:
        coerced = _coerce_temperature_c(melting_point or resolved_mp)
        if coerced is not None:
            properties["mp"] = coerced
    if reagent_type:
        properties["type"] = str(reagent_type).strip()
    else:
        mp_val = _coerce_temperature_float(melting_point or resolved_mp)
        bp_val = _coerce_temperature_float(boiling_point or resolved_bp)
        properties["type"] = _guess_reagent_type(mp_val, bp_val)
    properties.setdefault("volatile", "0")
    properties.setdefault("viscose", "0")
    if properties:
        entry["properties"] = properties

    result: Dict[str, Any] = {
        "cas": normalized_cas,
        "name": final_name,
        "role": role,
        "family_id": family_id,
        "registry_type": "flat",
        "registry_file": str(registry.file_for_role(role)),
        "dry_run": dry_run,
        "used_default_family": used_default,
        "entry_preview": entry,
    }
    if tag:
        result["tag"] = tag
    if family_candidates:
        result["family_candidates"] = family_candidates
    if auto_resolve_source:
        result["auto_resolve_source"] = auto_resolve_source
    if resolved_smiles and not smiles:
        result["smiles_source"] = auto_resolve_source or "resolver"
    if family_tokens:
        result["family_tokens"] = family_tokens
    if role_reason:
        result["role_reason"] = role_reason
    if default_rejection:
        result["default_rejection"] = default_rejection

    if dry_run:
        result["status"] = "dry_run"
        return result

    registry.add_entry(role, entry)
    path = registry.save_role(role)
    result["status"] = "written"
    result["written_to"] = str(path)
    return result


def list_available_families(taxonomy_dir: Optional[Path | str] = None) -> List[Dict[str, Any]]:
    """
    Return metadata for all known families.

    Args:
        taxonomy_dir: Optional path to reagent taxonomy v2 directory.

    Returns:
        List of dictionaries describing families (role, id, label).
    """
    taxonomy_path = Path(taxonomy_dir).expanduser() if taxonomy_dir else None
    store = _load_taxonomy(taxonomy_path)
    try:
        families = store.list_families()
        if families and isinstance(families[0], tuple):
            families = [
                {"role": role, "family_id": family, "label": label}
                for role, family, label in families
            ]
        return families
    except Exception as exc:  # pragma: no cover - defensive fallback
        raise ReagentAdditionError(f"Unable to list families: {exc}") from exc
