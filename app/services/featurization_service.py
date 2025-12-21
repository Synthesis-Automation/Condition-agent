"""
Featurization Service - Molecular and role-aware featurization.

This service handles:
- Molecular featurization for C-N coupling substrates
- Role-aware molecule featurization
- Role-aware reaction featurization
- Field registry information for role-aware features
"""

from typing import Dict, Any, List, Optional
from chemtools import featurizers
from chemtools.contracts import FeaturizeUllmannRequest, RoleAwareMolRequest, RoleAwareReactionRequest
from chemtools.exceptions import ValidationError, ResourceNotAvailableError

try:
    from chemtools.features.role import (
        featurize_mol as role_featurize_mol,
        featurize_reaction as role_featurize_reaction,
    )
    from chemtools.features.role.registry import REGISTRY as ROLE_REGISTRY
    _HAS_ROLE_FEATS = True
except Exception:
    _HAS_ROLE_FEATS = False


def featurize_molecular(req: FeaturizeUllmannRequest) -> Dict[str, Any]:
    """
    Featurize electrophile and nucleophile for C-N coupling predictions.
    
    Args:
        req: FeaturizeUllmannRequest with electrophile and nucleophile SMILES
        
    Returns:
        Dictionary with molecular features
        
    Raises:
        ValidationError: If SMILES are invalid
    """
    if not req.electrophile or not req.electrophile.strip():
        raise ValidationError("Electrophile SMILES cannot be empty")
    if not req.nucleophile or not req.nucleophile.strip():
        raise ValidationError("Nucleophile SMILES cannot be empty")
    
    return featurizers.reaction_pair.featurize_pair(req.electrophile, req.nucleophile)


def featurize_role_aware_molecule(req: RoleAwareMolRequest) -> Dict[str, Any]:
    """
    Featurize a molecule with role-aware features.
    
    Args:
        req: RoleAwareMolRequest with SMILES and optional roles
        
    Returns:
        Dictionary with role-aware features including vector
        
    Raises:
        ResourceNotAvailableError: If role-aware featurization is not available
        ValidationError: If SMILES is invalid
    """
    if not _HAS_ROLE_FEATS:
        raise ResourceNotAvailableError("role-aware featurization unavailable")
    
    if not req.smiles or not req.smiles.strip():
        raise ValidationError("SMILES cannot be empty")
    
    out = role_featurize_mol(req.smiles, roles=req.roles or None)
    
    # Convert numpy array to JSON-serializable list
    vec = out.get("vector")
    try:
        out["vector"] = vec.tolist()  # type: ignore
    except Exception:
        pass
    
    return out


def featurize_role_aware_reaction(req: RoleAwareReactionRequest) -> Dict[str, Any]:
    """
    Featurize a reaction with role-aware features for all reactants.
    
    Args:
        req: RoleAwareReactionRequest with reaction SMILES
        
    Returns:
        Dictionary with role-aware features for each reactant
        
    Raises:
        ResourceNotAvailableError: If role-aware featurization is not available
        ValidationError: If reaction SMILES is invalid
    """
    if not _HAS_ROLE_FEATS:
        raise ResourceNotAvailableError("role-aware featurization unavailable")
    
    if not req.reaction or not req.reaction.strip():
        raise ValidationError("Reaction SMILES cannot be empty")
    
    out = role_featurize_reaction(req.reaction)
    
    # Ensure vectors are JSON-serializable lists
    try:
        for item in out.get("reactants") or []:  # type: ignore[union-attr]
            vec = item.get("vector")
            try:
                item["vector"] = vec.tolist()  # type: ignore
            except Exception:
                pass
    except Exception:
        pass
    
    return out


def get_role_aware_fields(roles: Optional[str] = None) -> Dict[str, Any]:
    """
    Get role-aware field order and registry information.
    
    Args:
        roles: Comma-separated list of roles (default: "amine,alcohol,aryl_halide")
        
    Returns:
        Dictionary containing:
        - roles: List of roles used
        - fields: List of field names in order
        - counts: Field counts by category
        - registry: Complete role registry
        
    Raises:
        ResourceNotAvailableError: If role-aware featurization is not available
    """
    if not _HAS_ROLE_FEATS:
        raise ResourceNotAvailableError("role-aware featurization unavailable")
    
    # Parse and normalize roles
    default_roles = ["amine", "alcohol", "aryl_halide"]
    if roles is None or not str(roles).strip():
        use_roles = default_roles
    else:
        use_roles = [r.strip() for r in str(roles).split(",") if r.strip()]
        # Keep only known roles, preserve order
        known = {"amine", "alcohol", "aryl_halide"}
        use_roles = [r for r in use_roles if r in known]
        if not use_roles:
            use_roles = default_roles
    
    # Assemble fields: global -> role fields (in order) -> fingerprints per role
    fields: List[str] = []
    fields.extend([f.get("name", "") for f in ROLE_REGISTRY.get("global", [])])
    for r in use_roles:
        fields.extend([f.get("name", "") for f in ROLE_REGISTRY.get(r, [])])
    for r in use_roles:
        bits = int(ROLE_REGISTRY.get("fingerprints", {}).get(r, {}).get("bits", 512))
        fields.extend([f"fp.{r}.{i}" for i in range(bits)])
    
    counts = {
        "global": len(ROLE_REGISTRY.get("global", [])),
        "by_role": {r: len(ROLE_REGISTRY.get(r, [])) for r in use_roles},
        "fingerprints": {r: int(ROLE_REGISTRY.get("fingerprints", {}).get(r, {}).get("bits", 512)) for r in use_roles},
    }
    
    return {
        "roles": use_roles,
        "fields": fields,
        "counts": {**counts, "total": len(fields)},
        "registry": ROLE_REGISTRY,
    }


def is_role_aware_available() -> bool:
    """Check if role-aware featurization is available."""
    return _HAS_ROLE_FEATS
