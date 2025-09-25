from pathlib import Path
text = Path('app/ui_gradio.py').read_text(encoding='utf-8')
start = text.index('def ui_precedent_search(')
end = text.index('\n    return pack, html\n', start) + len('\n    return pack, html\n')
new_block = '''def ui_precedent_search(
    reaction: str,
    k: int,
    use_drfp: bool,
    drfp_weight: float,
    drfp_bits: int,
    drfp_radius: int,
    precompute_scope: str,
    use_molpipeline: bool,
    molpipeline_cfg_json: str,
    molpipeline_query_json: str,
) -> Tuple[Dict[str, Any], str, Dict[str, Any]]:
    norm = smiles.normalize_reaction(reaction or "")
    elec, nuc, reactants = _pick_elec_nuc_from_reaction(reaction or "")
    fam = router.detect_family(reactants).get("family") or "Unknown"
    feat = featurizers.molecular.featurize(elec, nuc)
    if isinstance(feat, dict) and "role_aware" in feat:
        try:
            feat = {k: v for k, v in feat.items() if k != "role_aware"}
        except Exception:
            feat.pop("role_aware", None)  # type: ignore

    relax: Dict[str, Any] = {
        "reaction_smiles": reaction or "",
        "use_drfp": bool(use_drfp),
        "drfp_weight": float(drfp_weight),
        "drfp_n_bits": int(drfp_bits),
        "drfp_radius": int(drfp_radius),
    }
    if precompute_scope in {"candidates", "all"}:
        relax["precompute_drfp"] = True
        relax["precompute_scope"] = precompute_scope
    else:
        relax["precompute_drfp"] = False

    molpipeline_summary: Dict[str, Any] = {
        "available": bool(_HAS_MOLPIPELINE),
    }
    if _HAS_MOLPIPELINE and _MOLPIPELINE_ENV:
        molpipeline_summary.update(
            {
                "version": _MOLPIPELINE_ENV.version,
                "rdkit": _MOLPIPELINE_ENV.rdkit_version,
            }
        )

    if use_molpipeline:
        if not _HAS_MOLPIPELINE:
            molpipeline_summary["enabled"] = False
            molpipeline_summary["error"] = "MolPipeline integration not available."
        else:
            cfg = _safe_json_loads(molpipeline_cfg_json)
            if not isinstance(cfg, dict):
                cfg = {}
            roles = cfg.get("roles")
            if isinstance(roles, list):
                cfg["roles"] = [str(r).strip().upper() for r in roles if isinstance(r, str) and r.strip()]
            else:
                cfg.pop("roles", None)
            cfg["include_role_features"] = bool(cfg.get("include_role_features", True))
            cfg["include_concat"] = bool(cfg.get("include_concat", True))
            cfg["suppress_errors"] = bool(cfg.get("suppress_errors", True))
            query_map = _safe_json_loads(molpipeline_query_json)
            if isinstance(query_map, dict) and query_map:
                cfg["query_role_smiles"] = query_map
            molpipeline_summary.update(
                {
                    "enabled": True,
                    "config": {
                        k: cfg[k]
                        for k in ("roles", "include_role_features", "include_concat")
                        if k in cfg
                    },
                }
            )
            relax["molpipeline"] = cfg
    else:
        molpipeline_summary["enabled"] = False

    pack = precedent.knn(fam, feat, k=int(k or 25), relax=relax)

    if use_molpipeline and _HAS_MOLPIPELINE:
        if pack.get("molpipeline_warnings"):
            molpipeline_summary["warnings"] = pack.get("molpipeline_warnings")
        if pack.get("molpipeline_query_vector"):
            vector = pack["molpipeline_query_vector"]
            molpipeline_summary["query_vector_length"] = len(vector) if isinstance(vector, list) else None
        first_vec = None
        for row in pack.get("precedents") or []:
            vec = row.get("molpipeline_feature_vector")
            if isinstance(vec, list):
                first_vec = vec
                break
        if first_vec is not None:
            molpipeline_summary["feature_length"] = len(first_vec)

    precs = list(pack.get("precedents") or [])
    html_rows: List[str] = []
    html_rows.append('<table style="border-collapse:collapse; width:100%; table-layout:fixed">')
    html_rows.append(
        "<colgroup>"
        "<col style='width:340px'/>"
        "<col style='width:110px'/>"
        "<col/>"
        "<col style='width:70px'/>"
        "<col style='width:160px'/>"
        "<col style='width:120px'/>"
        "<col style='width:120px'/>"
        "<col style='width:70px'/>"
        "<col style='width:70px'/>"
        "</colgroup>"
    )
    html_rows.append(
        "<tr>"
        "<th style='text-align:left'>image</th>"
        "<th style='text-align:left'>reaction_id</th>"
        "<th style='text-align:left'>reaction_smiles</th>"
        "<th style='text-align:left'>yield</th>"
        "<th style='text-align:left'>core</th>"
        "<th style='text-align:left'>base</th>"
        "<th style='text-align:left'>solvent</th>"
        "<th style='text-align:left'>T_C</th>"
        "<th style='text-align:left'>time_h</th>"
        "</tr>"
    )
    def _img_data_uri(img_obj: Any) -> str | None:
        try:
            import base64, io
            buf = io.BytesIO()
            img_obj.save(buf, format="PNG")  # type: ignore[attr-defined]
            b64 = base64.b64encode(buf.getvalue()).decode("ascii")
            return f"data:image/png;base64,{b64}"
        except Exception:
            return None

    def _render_img(reactants_smi: str, products_smi: str) -> str:
        try:
            from chemtools.util.rdkit_helpers import rdkit_available as _rd_avail
            if not _rd_avail():
                return ""
            from rdkit import Chem  # type: ignore
            from rdkit.Chem import Draw  # type: ignore

            def _grid(smi: str):
                ms = []
                for s in [x for x in (smi or '').split('.') if x]:
                    try:
                        m = Chem.MolFromSmiles(s)
                        if m is not None:
                            ms.append(m)
                    except Exception:
                        continue
                if not ms:
                    return None
                return Draw.MolsToGridImage(ms, molsPerRow=min(3, len(ms)), subImgSize=(220,220))

            l = _grid(reactants_smi)
            r = _grid(products_smi)
            if l is None and r is None:
                return ""
            if l is None:
                uri = _img_data_uri(r)
                return f"<img src='{uri}' width='320'/>" if uri else ""
            if r is None:
                uri = _img_data_uri(l)
                return f"<img src='{uri}' width='320'/>" if uri else ""
            try:
                from PIL import Image as _Image, ImageDraw as _ImageDraw  # type: ignore
                w = l.width + r.width + 60
                h = max(l.height, r.height)
                canvas = _Image.new("RGB", (w, h), (255,255,255))
                canvas.paste(l, (0, (h - l.height)//2))
                canvas.paste(r, (l.width + 60, (h - r.height)//2))
                dr = _ImageDraw.Draw(canvas)
                y = h//2
                dr.line((l.width + 10, y, l.width + 50, y), fill=(0,0,0), width=3)
                dr.polygon([(l.width + 50, y), (l.width + 38, y - 6), (l.width + 38, y + 6)], fill=(0,0,0))
                uri = _img_data_uri(canvas)
                return f"<img src='{uri}' width='320'/>" if uri else ""
            except Exception:
                uri = _img_data_uri(l)
                return f"<img src='{uri}' width='320'/>" if uri else ""
        except Exception:
            return ""

    def _resolve_name(uid: str) -> str:
        return _compound_display_label(str(uid or ""))

    for p in precs:
        img_html = _render_img(p.get("reactants_smiles", ""), p.get("products_smiles", ""))
        base_uid = str(p.get("base_uid") or "")
        solv_uid = str(p.get("solvent_uid") or "")
        base_name = _resolve_name(base_uid)
        solv_name = _resolve_name(solv_uid)
        html_rows.append(
            "<tr>"
            f"<td style='vertical-align:top'>{img_html}</td>"
            f"<td style='vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{p.get('reaction_id','')}</td>"
            f"<td style='vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{p.get('reaction_smiles','')}</td>"
            f"<td style='vertical-align:top; white-space:nowrap'>{p.get('yield','')}</td>"
            f"<td style='vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{p.get('core','')}</td>"
            f"<td title='{base_uid}' style='vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{base_name}</td>"
            f"<td title='{solv_uid}' style='vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{solv_name}</td>"
            f"<td style='vertical-align:top; white-space:nowrap'>{p.get('T_C','')}</td>"
            f"<td style='vertical-align:top; white-space:nowrap'>{p.get('time_h','')}</td>"
            "</tr>"
        )
    html_rows.append("</table>")
    html = "\n".join(html_rows)
    return pack, html, molpipeline_summary
'''
text = text[:start] + new_block + text[end:]
Path('app/ui_gradio.py').write_text(text, encoding='utf-8')
