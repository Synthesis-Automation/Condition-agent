"""MolPipeline integration for advanced featurization.

Provides optional MolPipeline feature extraction for precedent results.
"""
from typing import Dict, Any, Mapping, Sequence

try:
    from ..integrations import molpipeline as molpipe_integration
    _HAS_MOLPIPELINE = molpipe_integration.is_available()
except Exception:
    molpipe_integration = None  # type: ignore
    _HAS_MOLPIPELINE = False


def _attach_molpipeline_features(pack: Dict[str, Any], cfg: Any) -> Dict[str, Any]:
    """Attach MolPipeline features to precedent pack.
    
    Args:
        pack: Precedent result dict with 'precedents' list
        cfg: MolPipeline configuration dict
        
    Returns:
        Updated pack with molpipeline features added
    """
    if not isinstance(pack, dict):
        return pack
    if not isinstance(cfg, Mapping):
        return pack
    if not _HAS_MOLPIPELINE or molpipe_integration is None:
        return pack
    suppress = bool(cfg.get('suppress_errors', True))
    try:
        aggregator = cfg.get('aggregator') if isinstance(cfg, Mapping) else None
        if aggregator is None:
            roles = cfg.get('roles') if isinstance(cfg.get('roles'), Sequence) else None
            aggregator = molpipe_integration.build_default_role_aggregator(
                roles=roles,
                aggregate=str(cfg.get('aggregate', 'mean')),
                missing_strategy=str(cfg.get('missing_strategy', 'zeros')),
                n_jobs=int(cfg.get('n_jobs', 1)),
                ligand_n_bits=int(cfg.get('ligand_n_bits', 512)),
                ligand_radius=int(cfg.get('ligand_radius', 2)),
            )
    except Exception as exc:
        if suppress:
            pack.setdefault('molpipeline_warnings', []).append(str(exc))
            return pack
        raise
    include_role = bool(cfg.get('include_role_features', True))
    include_concat = bool(cfg.get('include_concat', True))
    role_key = str(cfg.get('role_features_key', 'molpipeline_role_features'))
    concat_key = str(cfg.get('concat_key', 'molpipeline_feature_vector'))
    precedents = pack.get('precedents')
    if isinstance(precedents, list):
        for row in precedents:
            if not isinstance(row, dict):
                continue
            try:
                role_feats = aggregator.featurize_roles(reaction=row)
                if include_role:
                    row[role_key] = {r: (vec.tolist() if hasattr(vec, 'tolist') else vec) if vec is not None else None for r, vec in role_feats.items()}
                if include_concat:
                    row[concat_key] = aggregator.concatenate(reaction=row).tolist()
            except Exception as exc:
                if suppress:
                    row.setdefault('molpipeline_error', str(exc))
                else:
                    raise
    try:
        query_role_smiles = cfg.get('query_role_smiles') if isinstance(cfg, Mapping) else None
        if query_role_smiles is None and isinstance(cfg.get('query_reaction'), Mapping):
            query_role_smiles = molpipe_integration.collect_role_smiles(cfg['query_reaction'])
        if query_role_smiles:
            pack['molpipeline_query_vector'] = aggregator.concatenate(role_smiles=query_role_smiles).tolist()
    except Exception as exc:
        if suppress:
            pack.setdefault('molpipeline_warnings', []).append(str(exc))
        else:
            raise
    return pack
