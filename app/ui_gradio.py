# -*- coding: utf-8 -*-
"""
Gradio UI to interactively test chemtools functions.

Run:
  python app/ui_gradio.py

Then open the URL printed by Gradio (default http://127.0.0.1:7860).
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple
from functools import lru_cache
from types import SimpleNamespace
import json
import os, sys

# Ensure project root is importable when running as a script
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import gradio as gr

# Theme: professional, compact spacing, subtle radius, readable fonts
THEME = gr.themes.Soft(
    primary_hue="slate",
    secondary_hue="zinc",
    neutral_hue="slate",
    font=[gr.themes.GoogleFont("Inter"), "ui-sans-serif", "system-ui", "-apple-system"],
    font_mono=[gr.themes.GoogleFont("JetBrains Mono"), "ui-monospace", "SFMono-Regular"],
    radius_size="sm",
    spacing_size="sm",
)

# Constrain overall width for a cleaner layout
CSS = """
.gradio-container{max-width:1600px !important;margin:0 auto !important;}
.gradio-container [role='tablist']{flex-wrap:wrap;overflow-x:auto;gap:0.25rem;}
.gradio-container [role='tab']{white-space:nowrap;}
"""

# Enable RDKit by default unless explicitly disabled by the environment
os.environ.setdefault("CHEMTOOLS_DISABLE_RDKIT", "0")
# Prefer rxn-insight enabled by default (can be disabled via env or UI toggle)
os.environ.setdefault("CHEMTOOLS_DISABLE_RXN_INSIGHT", "0")

from chemtools import smiles, router, properties, featurizers, recommend, precedent, reaction_similarity as rs
from chemtools import registry as creg
from chemtools.util.rdkit_helpers import rdkit_available

try:
    from chemtools.integrations import molpipeline as molpipe

    _MOLPIPELINE_ENV = molpipe.environment_snapshot()
    _HAS_MOLPIPELINE = _MOLPIPELINE_ENV.available
except Exception:
    molpipe = None  # type: ignore
    _MOLPIPELINE_ENV = SimpleNamespace(available=False, version=None, rdkit_version=None, sklearn_version=None, shap_version=None, import_error=None)  # type: ignore
    _HAS_MOLPIPELINE = False

# Optionally silence RDKit console logs (parse errors, warnings) to keep UI clean.
# Enabled by default; set CHEMTOOLS_RDKIT_QUIET=0 to show RDKit logs.
if os.environ.get("CHEMTOOLS_RDKIT_QUIET", "1").strip().lower() not in {"0", "false", "no", "off"}:
    try:
        from rdkit import RDLogger  # type: ignore
        # Prefer strict level; also disable common channels for wider compatibility
        try:
            RDLogger.logger().setLevel(RDLogger.CRITICAL)
        except Exception:
            pass
        for channel in ("rdApp.error", "rdApp.warning", "rdApp.info"):
            try:
                RDLogger.DisableLog(channel)
            except Exception:
                pass
    except Exception:
        pass
# Optional role-aware single-molecule featurizer
try:
    from chem_feats import featurize_mol as role_feat_mol  # type: ignore
except Exception:
    role_feat_mol = None  # type: ignore


@lru_cache(maxsize=512)
def _registry_resolve_cached(query: str) -> Dict[str, Any]:
    q = (query or "").strip()
    if not q:
        return {"error": "EMPTY_QUERY"}
    try:
        res = creg.resolve(q)
        if isinstance(res, dict):
            return res
    except Exception as exc:  # pragma: no cover - UI helper
        return {"error": str(exc)}
    return {"error": "NOT_FOUND"}


@lru_cache(maxsize=512)
def _compound_display_label(value: str) -> str:
    candidate = (value or "").strip()
    if not candidate:
        return ""
    try:
        prop_res = properties.lookup(candidate)
        if isinstance(prop_res, dict) and prop_res.get("found"):
            rec = prop_res.get("record")
            if isinstance(rec, dict):
                for key in ("token", "name", "uid"):
                    val = rec.get(key)
                    if isinstance(val, str) and val.strip():
                        return val.strip()
    except Exception:  # pragma: no cover - defensive UI helper
        pass
    detail = _registry_resolve_cached(candidate)
    if isinstance(detail, dict) and not detail.get("error"):
        props = detail.get("props") if isinstance(detail.get("props"), dict) else {}
        if isinstance(props, dict):
            for key in ("token", "abbreviation", "name", "generic_core"):
                val = props.get(key)
                if isinstance(val, str) and val.strip():
                    return val.strip()
        for key in ("name", "uid"):
            val = detail.get(key)
            if isinstance(val, str) and val.strip():
                return val.strip()
        aliases = detail.get("aliases")
        if isinstance(aliases, list):
            for alias in aliases:
                if isinstance(alias, str):
                    alias = alias.strip()
                    if alias:
                        return alias
    return candidate


def _safe_json_loads(s: str) -> Dict[str, Any]:
    s = (s or "").strip()
    if not s:
        return {}
    try:
        obj = json.loads(s)
        return obj if isinstance(obj, dict) else {}
    except Exception:
        return {}


# --- Handlers for tabs ---

def ui_normalize_smiles(smi: str) -> Dict[str, Any]:
    # If input looks like a reaction (contains '>' or '>>'), route to reaction normalizer.
    # Do not fall back to single-molecule parsing on failure, to avoid RDKit
    # attempting to parse reaction strings like 'A>>B' as molecule SMILES.
    t = (smi or "").strip()
    if ">" in t:
        return smiles.normalize_reaction(t)
    return smiles.normalize(t)


def _detect_family_md(use_rxi_requested: bool, out: Dict[str, Any]) -> str:
    lines: list[str] = []
    # Comparison of rule-based vs rxn-insight
    rule = out.get("rule") if isinstance(out, dict) else None
    rxn = out.get("rxn") if isinstance(out, dict) else None
    rf = (rule or {}).get("family") if isinstance(rule, dict) else None
    rxf = (rxn or {}).get("family") if isinstance(rxn, dict) else None
    agree = bool(rf and rxf and (str(rf) == str(rxf)))
    if agree and rf:
        lines.append(f"Detected family: {rf}")
    elif rf and rxf:
        lines.append(f"Conflict: rule-based={rf} vs rxn-insight={rxf}")
    fam = out.get("family")
    if fam and not agree:
        lines.append(f"final choice: {fam}")
    # Details
    if isinstance(rule, dict) and rf:
        rc = rule.get("confidence")
        try:
            lines.append(f"rule-based: {rf} (conf {float(rc):.3f})")
        except Exception:
            lines.append(f"rule-based: {rf} (conf {rc})")
    if isinstance(rxn, dict) and rxf:
        rxc = rxn.get("confidence")
        try:
            lines.append(f"rxn-insight mapped: {rxf} (conf {float(rxc):.3f} if present)")
        except Exception:
            lines.append(f"rxn-insight mapped: {rxf}")
    if rf and rxf:
        lines.append(f"agreement: {'yes' if agree else 'no'}")
    auto = out.get("auto") if isinstance(out, dict) else None
    if isinstance(auto, dict):
        avail = bool(auto.get("available"))
        status = "success" if auto.get("success") else "fallback"
        name = auto.get("rxn_name")
        klass = auto.get("rxn_class")
        conf = auto.get("confidence")
        lines.append(f"rxn-insight available: {avail}")
        lines.append(f"auto-detection status: {status}")
        if name:
            lines.append(f"rxn-insight name: {name}")
        if klass:
            lines.append(f"rxn-insight class: {klass}")
        if conf is not None:
            try:
                lines.append(f"rxn-insight confidence: {float(conf):.3f}")
            except Exception:
                lines.append(f"rxn-insight confidence: {conf}")
    else:
        if use_rxi_requested:
            lines.append("rxn-insight: not available; using rule-based fallback")
        else:
            lines.append("rxn-insight: disabled; using rule-based detection")
    return "\n".join(f"- {x}" for x in lines) if lines else "- No detection details"


def ui_detect_family(reactants_text: str, use_rxn_insight: bool) -> Tuple[Dict[str, Any], str]:
    # Accept either a dot-separated list of reactant SMILES or a full reaction SMILES.
    # If a reaction arrow is present, prefer the rxn-insight path when enabled.
    t = (reactants_text or "").strip()
    if not t:
        return {"family": "Unknown", "confidence": 0.0, "hits": {}}
    if ">" in t:
        try:
            out = router.detect_family_from_reaction(t, use_rxn_insight=bool(use_rxn_insight))
            # If both methods are present and disagree, set family to CONFLICT
            rule = out.get("rule") if isinstance(out, dict) else None
            rxn = out.get("rxn") if isinstance(out, dict) else None
            rf = (rule or {}).get("family") if isinstance(rule, dict) else None
            rxf = (rxn or {}).get("family") if isinstance(rxn, dict) else None
            if rf and rxf and str(rf) != str(rxf):
                out = dict(out)
                out["family"] = "CONFLICT"
                out["conflict_between"] = {"rule": rf, "rxn": rxf}
                out["status"] = "conflict"
            return out, _detect_family_md(bool(use_rxn_insight), out)
        except Exception:
            pass
    # Split reactants by '.' or newlines
    if "\n" in t:
        reactants = [s.strip() for s in t.splitlines() if s.strip()]
    else:
        reactants = [s.strip() for s in t.split(".") if s.strip()]
    if use_rxn_insight:
        # Build a reaction string with empty products to enable rxn-insight route
        rxn_text = ".".join(reactants) + ">>"
        try:
            out = router.detect_family_from_reaction(rxn_text, use_rxn_insight=True)
            rule = out.get("rule") if isinstance(out, dict) else None
            rxn = out.get("rxn") if isinstance(out, dict) else None
            rf = (rule or {}).get("family") if isinstance(rule, dict) else None
            rxf = (rxn or {}).get("family") if isinstance(rxn, dict) else None
            if rf and rxf and str(rf) != str(rxf):
                out = dict(out)
                out["family"] = "CONFLICT"
                out["conflict_between"] = {"rule": rf, "rxn": rxf}
                out["status"] = "conflict"
            return out, _detect_family_md(True, out)
        except Exception:
            pass
    out = router.detect_family(reactants)
    return out, _detect_family_md(False, out)


def ui_featurize_molecular(electrophile: str, nucleophile: str) -> Dict[str, Any]:
    return featurizers.molecular.featurize(electrophile or "", nucleophile or "")


def ui_properties_lookup(query: str) -> Dict[str, Any]:
    q = (query or "").strip()
    result: Dict[str, Any] = {"query": q}
    props_out = properties.lookup(q)
    result["properties"] = props_out
    registry_out = None
    lookup_keys: List[str] = []
    if isinstance(props_out, dict):
        rec = props_out.get("record")
        if isinstance(rec, dict) and rec.get("uid"):
            lookup_keys.append(str(rec["uid"]))
    if q:
        lookup_keys.append(q)
    for key in lookup_keys:
        detail = _registry_resolve_cached(key)
        if isinstance(detail, dict) and not detail.get("error"):
            registry_out = detail
            break
    if registry_out is not None:
        result["registry"] = registry_out
    elif q:
        result["registry"] = {"error": "NOT_FOUND"}
    else:
        result["registry"] = {"error": "EMPTY_QUERY"}
    return result


def ui_recommend(
    reaction: str,
    k: int,
    use_rxn_insight: bool,
    relax_json: str,
    constraints_json: str,
) -> Dict[str, Any]:
    relax = _safe_json_loads(relax_json)
    relax["use_rxn_insight"] = bool(use_rxn_insight)
    constraints = _safe_json_loads(constraints_json)
    return recommend.recommend_from_reaction(
        reaction=reaction or "",
        k=int(k or 25),
        relax=relax,
        constraint_rules=constraints,
    )


def ui_recommend_both(
    reaction: str,
    k: int,
    use_rxn_insight: bool,
    relax_json: str,
    constraints_json: str,
) -> Tuple[Dict[str, Any], List[List[Any]], str, str]:
    out = ui_recommend(reaction, k, use_rxn_insight, relax_json, constraints_json)
    # Build a compact 1-row table for the main recommendation
    hdr = ["core", "base_uid", "solvent_uid", "T_C", "time_h", "confidence"]
    rec = out.get("recommendation") or {}
    row = [
        rec.get("core", ""),
        rec.get("base_uid", ""),
        rec.get("solvent_uid", ""),
        rec.get("T_C", ""),
        rec.get("time_h", ""),
        rec.get("confidence", ""),
    ]
    # Human-readable summary using names/tokens looked up from registry/properties
    def _resolve_name(uid: str) -> str:
        return _compound_display_label(str(uid or ""))
    core_txt = str(rec.get("core") or "")
    # strip trailing /none if present
    if core_txt.endswith("/none"):
        core_txt = core_txt[:-5]
    base_txt = _resolve_name(str(rec.get("base_uid") or ""))
    solv_txt = _resolve_name(str(rec.get("solvent_uid") or ""))
    human = "/".join([s for s in [core_txt, base_txt, solv_txt] if s]) or ""

    # Detection markdown (rxn-insight and fallback details)
    det = (out.get("formatted") or {}).get("detection") or {}
    auto = det.get("auto") if isinstance(det, dict) else None
    det_lines: list[str] = []
    if isinstance(det, dict):
        rt = det.get("reaction_type")
        if rt:
            det_lines.append(f"Detected family: {rt}")
    if isinstance(auto, dict):
        avail = bool(auto.get("rxn_insight_available"))
        name = auto.get("rxn_insight_name")
        klass = auto.get("rxn_insight_class")
        conf = auto.get("rxn_insight_confidence")
        status = auto.get("status")
        det_lines.append(f"rxn-insight available: {avail}")
        if status:
            det_lines.append(f"auto-detection status: {status}")
        if name:
            det_lines.append(f"rxn-insight name: {name}")
        if klass:
            det_lines.append(f"rxn-insight class: {klass}")
        if conf is not None:
            try:
                det_lines.append(f"rxn-insight confidence: {float(conf):.3f}")
            except Exception:
                det_lines.append(f"rxn-insight confidence: {conf}")
    elif det and not auto:
        det_lines.append("rxn-insight: not available; using rule-based fallback")
    detection_md = "\n".join(f"- {x}" for x in det_lines) if det_lines else "- No detection details"

    return out, [hdr, row], human, detection_md


def ui_design_plate(
    reaction: str,
    plate_size: int,
    use_rxn_insight: bool,
    relax_json: str,
    constraints_json: str,
) -> Tuple[str, List[List[Any]], Dict[str, Any]]:
    relax = _safe_json_loads(relax_json)
    relax["use_rxn_insight"] = bool(use_rxn_insight)
    constraints = _safe_json_loads(constraints_json)
    out = recommend.design_plate_from_reaction(
        reaction=reaction or "",
        plate_size=int(plate_size or 24),
        relax=relax,
        constraint_rules=constraints,
    )
    # Prepare a simple 2D list for Dataframe display (header + rows)
    header = ["well_id", "core", "base_uid", "solvent_uid", "additive_uids", "T_C", "time_h"]
    rows = out.get("rows") or []
    table = [header]
    for r in rows:
        table.append([
            r.get("well_id", ""),
            r.get("core", ""),
            r.get("base_uid", ""),
            r.get("solvent_uid", ""),
            r.get("additive_uids", ""),
            r.get("T_C", ""),
            r.get("time_h", ""),
        ])
    meta = out.get("meta") or {}
    return out.get("csv", ""), table, meta


def _pick_elec_nuc_from_reaction(rsmi: str) -> Tuple[str, str, List[str]]:
    from chemtools.smiles import normalize_reaction
    norm = normalize_reaction(rsmi or "")
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]
    def is_electrophile(s: str) -> bool:
        t = (s or "").lower()
        return ("br" in t) or ("cl" in t) or (" i" in t) or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
    elec, nuc = "", ""
    if reactants:
        if len(reactants) == 1:
            elec, nuc = reactants[0], ""
        else:
            r0, r1 = reactants[0], reactants[1]
            elec, nuc = (r0, r1) if is_electrophile(r0) else ((r1, r0) if is_electrophile(r1) else (r0, r1))
    return elec, nuc, reactants


def ui_precedent_search(
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
    # Detect family and featurize from reaction
    elec, nuc, reactants = _pick_elec_nuc_from_reaction(reaction or "")
    fam = router.detect_family(reactants).get("family") or "Unknown"
    feat = featurizers.molecular.featurize(elec, nuc)
    # Drop nested role_aware to keep features hashable for caching
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
    if _MOLPIPELINE_ENV is not None:
        molpipeline_summary.update(
            {
                "version": getattr(_MOLPIPELINE_ENV, "version", None),
                "rdkit": getattr(_MOLPIPELINE_ENV, "rdkit_version", None),
                "sklearn": getattr(_MOLPIPELINE_ENV, "sklearn_version", None),
                "shap": getattr(_MOLPIPELINE_ENV, "shap_version", None),
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
            sanitized_cfg: Dict[str, Any] = {}
            roles = cfg.get("roles")
            if isinstance(roles, list):
                clean_roles = []
                for role in roles:
                    if isinstance(role, str) and role.strip():
                        clean_roles.append(role.strip().upper())
                if clean_roles:
                    sanitized_cfg["roles"] = clean_roles
            include_role = cfg.get("include_role_features")
            sanitized_cfg["include_role_features"] = bool(True if include_role is None else include_role)
            include_concat = cfg.get("include_concat")
            sanitized_cfg["include_concat"] = bool(True if include_concat is None else include_concat)
            sanitized_cfg["suppress_errors"] = bool(cfg.get("suppress_errors", True))
            if isinstance(cfg.get("aggregate"), str) and cfg["aggregate"].strip():
                sanitized_cfg["aggregate"] = cfg["aggregate"].strip()
            if isinstance(cfg.get("missing_strategy"), str) and cfg["missing_strategy"].strip():
                sanitized_cfg["missing_strategy"] = cfg["missing_strategy"].strip()
            for key in ("n_jobs", "ligand_n_bits", "ligand_radius"):
                if key in cfg:
                    try:
                        sanitized_cfg[key] = int(cfg[key])
                    except Exception:
                        continue
            query_raw = _safe_json_loads(molpipeline_query_json)
            if isinstance(query_raw, dict) and query_raw:
                query_map: Dict[str, List[str]] = {}
                for role, value in query_raw.items():
                    if not isinstance(role, str):
                        continue
                    role_norm = role.strip().upper()
                    if not role_norm:
                        continue
                    smiles_list: List[str] = []
                    if isinstance(value, str):
                        trimmed = value.strip()
                        if trimmed:
                            smiles_list = [trimmed]
                    elif isinstance(value, (list, tuple, set)):
                        for item in value:
                            if isinstance(item, str) and item.strip():
                                smiles_list.append(item.strip())
                    if smiles_list:
                        query_map[role_norm] = smiles_list
                if query_map:
                    sanitized_cfg["query_role_smiles"] = query_map
            molpipeline_summary.update(
                {
                    "enabled": True,
                    "config": {
                        key: sanitized_cfg[key]
                        for key in (
                            "roles",
                            "include_role_features",
                            "include_concat",
                            "aggregate",
                            "missing_strategy",
                            "n_jobs",
                            "ligand_n_bits",
                            "ligand_radius",
                        )
                        if key in sanitized_cfg
                    },
                }
            )
            if "query_role_smiles" in sanitized_cfg:
                molpipeline_summary.setdefault("config", {})["query_roles"] = sorted(sanitized_cfg["query_role_smiles"].keys())
            relax["molpipeline"] = sanitized_cfg
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

    # Build HTML table with embedded reaction images; omit reactant/product SMILES columns
    precs = list(pack.get("precedents") or [])
    html_rows: List[str] = []
    html_rows.append('<table style="border-collapse:collapse; width:100%; table-layout:fixed">')
    html_rows.append(
        "<colgroup>"
        "<col style='width:340px'/>"   # image
        "<col style='width:110px'/>"   # reaction_id
        "<col/>"                        # reaction_smiles (flex)
        "<col style='width:70px'/>"    # yield
        "<col style='width:160px'/>"   # core
        "<col style='width:120px'/>"   # base_uid
        "<col style='width:120px'/>"   # solvent_uid
        "<col style='width:70px'/>"    # T_C
        "<col style='width:70px'/>"    # time_h
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
    # Image helpers
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


def ui_similarity_tanimoto(q: str, r: str, n_bits: int, radius: int) -> Dict[str, Any]:
    if not q or not r:
        return {"ok": False, "error": "provide two reaction SMILES"}
    if not rs.drfp_available():
        return {"ok": False, "error": "drfp/numpy not available in this environment"}
    a = rs.encode_drfp_cached(q, n_bits=int(n_bits or 4096), radius=int(radius or 3))
    b = rs.encode_drfp_cached(r, n_bits=int(n_bits or 4096), radius=int(radius or 3))
    sim = rs.tanimoto(a, b)
    try:
        sim = float(sim)
    except Exception:
        sim = 0.0
    return {"ok": True, "tanimoto": round(sim, 4)}


def ui_core_search(core: str, family: str, fuzzy: bool, limit: int):
    fam = (family or "").strip() or None
    rows = precedent.find_reactions_by_core(core or "", family=fam, fuzzy=bool(fuzzy), limit=int(limit or 25))
    payload = {
        "query": {"core": core, "family": fam, "fuzzy": bool(fuzzy), "limit": int(limit or 25)},
        "count": len(rows),
    }
    # Build table (header + rows)
    header = [
        "reaction_id",
        "rxn_type",
        "condition_core",
        "yield_value",
        "base_uid",
        "solvent_uid",
        "T_C",
        "time_h",
        "reaction_smiles",
    ]
    table: List[List[Any]] = [header]
    for r in rows:
        table.append([
            r.get("reaction_id", ""),
            r.get("rxn_type", ""),
            r.get("condition_core", ""),
            r.get("yield_value", ""),
            r.get("base_uid", ""),
            r.get("solvent_uid", ""),
            r.get("T_C", ""),
            r.get("time_h", ""),
            r.get("reaction_smiles", ""),
        ])
    return payload, table


def build_demo() -> gr.Blocks:
    default_roles = ["ADDITIVE", "BASE", "CATALYST", "LIGAND", "SOLVENT"]
    roles = list(default_roles)
    compound_types: List[str] = []
    try:
        cats = creg.categories()
        if isinstance(cats, dict):
            role_list = cats.get("roles")
            if isinstance(role_list, list) and role_list:
                roles = [str(r) for r in role_list if isinstance(r, str)] or roles
            ct_list = cats.get("compound_types")
            if isinstance(ct_list, list):
                compound_types = [str(c) for c in ct_list if isinstance(c, str)]
    except Exception:  # pragma: no cover - UI helper
        pass
    # Preserve order while removing duplicates
    roles = list(dict.fromkeys(roles))
    compound_types = list(dict.fromkeys(compound_types))
    role_choices = [""] + roles
    compound_type_choices = [""] + compound_types
    with gr.Blocks(title="ChemTools UI", theme=THEME, css=CSS) as demo:
        gr.Markdown("""
        # ChemTools - Interactive UI
        Try common chemistry tools without writing code. Use the tabs below.
        """)

        with gr.Tab("SMILES Normalize"):
            smi_in = gr.Textbox(label="SMILES", value="c1ccccc1O")
            smi_btn = gr.Button("Normalize", variant="primary")
            smi_out = gr.JSON(label="Result")
            smi_btn.click(ui_normalize_smiles, inputs=[smi_in], outputs=[smi_out])

        with gr.Tab("Detect Family"):
            gr.Markdown("Enter reactants separated by '.' or new lines, or a full reaction SMILES.")
            react_in = gr.Textbox(label="Reactants or Reaction", value="Clc1ccccc1.OB(O)c1ccccc1>>")
            react_use_rxi = gr.Checkbox(value=True, label="Use rxn-insight auto-detect (if available)")
            react_btn = gr.Button("Detect", variant="primary")
            react_out = gr.JSON(label="Detected Family")
            react_md = gr.Markdown(label="Auto-Detection")
            react_btn.click(ui_detect_family, inputs=[react_in, react_use_rxi], outputs=[react_out, react_md])

        # Removed the legacy multi-molecule featurizer tab to keep only
        # the single-molecule (basic) featurizer in the UI.

        with gr.Tab("Single Molecule (basic)"):
            smi_in = gr.Textbox(label="SMILES", value="Clc1ccccc1")
            roles_in = gr.CheckboxGroup(
                choices=["amine", "alcohol", "aryl_halide"],
                label="Roles (optional; leave empty for globals-only)",
                value=[],
            )
            show_full = gr.Checkbox(value=False, label="Show full fields list (unchecked shows preview)")
            with gr.Row():
                smi_btn = gr.Button("Featurize molecule", variant="primary")
                dl_btn = gr.DownloadButton(label="Download CSV")
            smi_out = gr.JSON(label="Result (vector length, masks, fields)")
            rdkit_status = gr.Markdown(label="Status")
            csv_code = gr.Code(label="CSV preview")

            def _ui_single_molecule(smi: str, roles: list[str], show_full_fields: bool):
                if role_feat_mol is None:
                    return {"error": "role-aware featurizer unavailable"}, ""
                roles = roles or []
                out = role_feat_mol(smi or "", roles=roles)
                vec = out.get("vector")
                try:
                    out["vector"] = vec.tolist()  # type: ignore
                except Exception:
                    pass
                out["length"] = len(out.get("vector") or [])
                # fields handling
                flds = out.get("fields") or []
                if not show_full_fields:
                    out["fields_preview"] = flds[:12]
                    try:
                        del out["fields"]
                    except Exception:
                        pass
                # RDKit status hint
                status_msg = ""
                if not rdkit_available():
                    status_msg = (
                        "?? RDKit is not available or disabled. Global descriptors, role-specific details, "
                        "and fingerprints default to 0. Set `CHEMTOOLS_DISABLE_RDKIT=0` (or unset it) and install RDKit "
                        "to enable numeric features."
                    )
                else:
                    status_msg = "RDKit detected ?"
                return out, status_msg

            def _to_csv(smi: str, roles: list[str]):
                if role_feat_mol is None:
                    return b"error,role-aware featurizer unavailable\n"
                roles = roles or []
                out = role_feat_mol(smi or "", roles=roles)
                vec = out.get("vector")
                try:
                    vec_list = vec.tolist()  # type: ignore
                except Exception:
                    vec_list = list(vec or [])  # type: ignore
                fields = out.get("fields") or []
                # Build CSV: field,value
                lines = ["field,value"]
                n = min(len(fields), len(vec_list))
                for i in range(n):
                    name = str(fields[i]).replace('"', '""')
                    val = vec_list[i]
                    lines.append(f'"{name}",{val}')
                csv_text = "\n".join(lines) + "\n"
                return csv_text.encode("utf-8")

            smi_btn.click(_ui_single_molecule, inputs=[smi_in, roles_in, show_full], outputs=[smi_out, rdkit_status])
            # Update CSV preview and download
            def _update_csv_preview(smi: str, roles: list[str]) -> str:
                data = _to_csv(smi, roles)
                try:
                    return data.decode("utf-8") if isinstance(data, (bytes, bytearray)) else str(data)
                except Exception:
                    return str(data)
            smi_btn.click(_update_csv_preview, inputs=[smi_in, roles_in], outputs=[csv_code])
            dl_btn.click(_to_csv, inputs=[smi_in, roles_in], outputs=[dl_btn])

        with gr.Tab("Properties Lookup"):
            prop_in = gr.Textbox(label="Query (name, CAS, token)", value="Water")
            prop_btn = gr.Button("Lookup", variant="primary")
            prop_out = gr.JSON(label="Record")
            prop_btn.click(ui_properties_lookup, inputs=[prop_in], outputs=[prop_out])

            gr.Markdown("""
            ### Registry Categories and Listing
            Use the controls below to list items by role or compound type (e.g., solvents).
            """)
            # Handlers for categories and search
            def ui_registry_categories() -> Dict[str, Any]:
                try:
                    cats = creg.categories()
                except Exception as e:
                    return {"error": str(e)}
                if not isinstance(cats, dict):
                    cats = {}
                result: Dict[str, Any] = dict(cats)
                try:
                    if "taxonomy_dir" not in result:
                        env_dir = os.environ.get("CHEMTOOLS_TAXONOMY_DIR")
                        if env_dir and os.path.isdir(env_dir):
                            result["taxonomy_dir"] = env_dir
                        else:
                            fallback = os.path.join(ROOT, "data", "compound_taxonomy")
                            if os.path.isdir(fallback):
                                result["taxonomy_dir"] = fallback
                except Exception:
                    pass
                return result

            def ui_registry_search(q: str, role: str, compound_type: str, limit: int) -> List[List[Any]]:
                rows = creg.search(q=(q or None), role=(role or None), compound_type=(compound_type or None), limit=int(limit or 50))
                header = ["uid", "role", "name", "compound_type", "token", "generic_core", "aliases (sample)"]
                table: List[List[Any]] = [header]
                for r in rows:
                    uid = str(r.get("uid") or "")
                    token = ""
                    generic = ""
                    alias_text = ""
                    if uid:
                        detail = _registry_resolve_cached(uid)
                        if isinstance(detail, dict) and not detail.get("error"):
                            props = detail.get("props") if isinstance(detail.get("props"), dict) else {}
                            if isinstance(props, dict):
                                token = str(props.get("token") or props.get("abbreviation") or props.get("name") or "").strip()
                                generic = str(props.get("generic_core") or "").strip()
                            aliases = detail.get("aliases") if isinstance(detail.get("aliases"), list) else []
                            if aliases:
                                alias_parts: List[str] = []
                                for alias in aliases:
                                    if isinstance(alias, str):
                                        alias = alias.strip()
                                        if alias:
                                            alias_parts.append(alias)
                                    if len(alias_parts) >= 5:
                                        break
                                alias_text = ", ".join(alias_parts)
                    table.append([uid, r.get("role", ""), r.get("name", ""), r.get("compound_type", ""), token, generic, alias_text])
                return table

            with gr.Row():
                cat_btn = gr.Button("List categories", variant="secondary")
                cat_json = gr.JSON(label="Categories (roles, compound_types)")
                cat_btn.click(ui_registry_categories, outputs=[cat_json])
            with gr.Row():
                rs_q = gr.Textbox(label="Filter (optional substring)", value="")
                rs_role = gr.Dropdown(label="Role", choices=role_choices, value="", allow_custom_value=True)
                rs_ct = gr.Dropdown(label="Compound type (optional)", choices=compound_type_choices, value="", allow_custom_value=True)
                rs_limit = gr.Slider(label="Limit", minimum=1, maximum=500, value=50, step=1)
            rs_btn = gr.Button("List by filter", variant="secondary")
            rs_tbl = gr.Dataframe(label="Registry items", interactive=False)
            rs_btn.click(ui_registry_search, inputs=[rs_q, rs_role, rs_ct, rs_limit], outputs=[rs_tbl])

        with gr.Tab("Recommend Conditions"):
            rec_in = gr.Textbox(label="Reaction SMILES", value="Brc1ccccc1.Nc1ccccc1>>")
            rec_k = gr.Slider(label="k (neighbors)", minimum=5, maximum=100, value=25, step=1)
            rec_use_rxi = gr.Checkbox(label="Use rxn-insight auto-detect (if available)", value=True)
            rec_relax = gr.Textbox(label="Relax (JSON)", value="")
            rec_constraints = gr.Textbox(label="Constraints (JSON)", value="")
            rec_btn = gr.Button("Recommend", variant="primary")
            rec_out = gr.JSON(label="Recommendation Pack")
            rec_tbl = gr.Dataframe(label="Recommendation (table)", interactive=False)
            rec_human = gr.Textbox(label="Recommended (human?readable core/base/solvent)", interactive=False)
            rec_detect = gr.Markdown(label="Auto-Detection")
            rec_btn.click(ui_recommend_both, inputs=[rec_in, rec_k, rec_use_rxi, rec_relax, rec_constraints], outputs=[rec_out, rec_tbl, rec_human, rec_detect])

        with gr.Tab("Design Plate"):
            plate_in = gr.Textbox(label="Reaction SMILES", value="Brc1ccccc1.Nc1ccccc1>>")
            plate_n = gr.Slider(label="Plate size", minimum=6, maximum=96, value=24, step=1)
            plate_use_rxi = gr.Checkbox(label="Use rxn-insight auto-detect (if available)", value=True)
            plate_relax = gr.Textbox(label="Relax (JSON)", value="")
            plate_constraints = gr.Textbox(label="Constraints (JSON)", value="")
            plate_btn = gr.Button("Design", variant="primary")
            plate_csv = gr.Code(label="CSV")
            plate_tbl = gr.Dataframe(label="Preview", interactive=False)
            plate_meta = gr.JSON(label="Meta")
            plate_btn.click(
                ui_design_plate,
                inputs=[plate_in, plate_n, plate_use_rxi, plate_relax, plate_constraints],
                outputs=[plate_csv, plate_tbl, plate_meta],
            )

        with gr.Tab("Precedent Search"):
            ps_in = gr.Textbox(label="Reaction SMILES", value="Brc1ccccc1.Nc1ccccc1>>")
            ps_k = gr.Slider(label="k (neighbors)", minimum=5, maximum=200, value=50, step=1)
            with gr.Row():
                ps_use_drfp = gr.Checkbox(label="Use DRFP re-ranking", value=True)
                ps_drfp_w = gr.Slider(label="DRFP weight", minimum=0.0, maximum=1.0, value=0.4, step=0.05)
            with gr.Row():
                ps_bits = gr.Number(label="DRFP bits", value=4096, precision=0)
                ps_radius = gr.Number(label="DRFP radius", value=3, precision=0)
                ps_prec_scope = gr.Radio(label="Precompute scope", choices=["none", "candidates", "all"], value="candidates")
            if _HAS_MOLPIPELINE:
                env_lines = ["MolPipeline detected."]
                if getattr(_MOLPIPELINE_ENV, "version", None):
                    env_lines.append(f"version: {_MOLPIPELINE_ENV.version}")
                extras = []
                for attr, label in (("rdkit_version", "rdkit"), ("sklearn_version", "sklearn"), ("shap_version", "shap")):
                    val = getattr(_MOLPIPELINE_ENV, attr, None)
                    if val:
                        extras.append(f"{label}: {val}")
                if extras:
                    env_lines.extend(extras)
                molpipe_status_md = "\n".join(env_lines)
            else:
                molpipe_status_md = ("MolPipeline extras not installed. Run `pip install chemtools-project[molpipeline]` to enable role features.")
            with gr.Accordion("MolPipeline (optional)", open=False):
                gr.Markdown(molpipe_status_md)
                ps_use_molpipe = gr.Checkbox(label="Attach MolPipeline role features", value=False)
                ps_mol_cfg = gr.Textbox(
                    label="MolPipeline config (JSON)",
                    placeholder="{\"roles\": [\"LIGAND\", \"BASE\", \"SOLVENT\"], \"aggregate\": \"mean\"}",
                    lines=4,
                )
                ps_mol_query = gr.Textbox(
                    label="Query role SMILES override (JSON)",
                    placeholder="{\"LIGAND\": [\"PPh3\"], \"BASE\": [\"K3PO4\"]}",
                    lines=3,
                )
            ps_btn = gr.Button("Search", variant="primary")
            ps_pack = gr.JSON(label="Pack (prototype, support, precedents)")
            ps_tbl = gr.HTML(label="Top precedents")
            ps_mol_summary = gr.JSON(label="MolPipeline summary")
            ps_btn.click(
                ui_precedent_search,
                inputs=[
                    ps_in,
                    ps_k,
                    ps_use_drfp,
                    ps_drfp_w,
                    ps_bits,
                    ps_radius,
                    ps_prec_scope,
                    ps_use_molpipe,
                    ps_mol_cfg,
                    ps_mol_query,
                ],
                outputs=[ps_pack, ps_tbl, ps_mol_summary],
            )

        with gr.Tab("DRFP Similarity"):
            s_q = gr.Textbox(label="Query reaction SMILES", value="Brc1ccccc1.Nc1ccccc1>>")
            s_r = gr.Textbox(label="Reference reaction SMILES", value="Clc1ccccc1.Nc1ccccc1>>")
            s_bits = gr.Number(label="DRFP bits", value=4096, precision=0)
            s_radius = gr.Number(label="DRFP radius", value=3, precision=0)
            s_btn = gr.Button("Compute Tanimoto", variant="primary")
            s_out = gr.JSON(label="Result")
            s_btn.click(ui_similarity_tanimoto, inputs=[s_q, s_r, s_bits, s_radius], outputs=[s_out])

        with gr.Tab("Core Search"):
            cs_core = gr.Textbox(label="Condition core (e.g., Pd/XPhos or XPhos)", value="Pd/XPhos")
            cs_family = gr.Dropdown(
                label="Family (optional)",
                choices=["", "Ullmann C-N", "Buchwald C-N", "Suzuki_CC", "Amide_Coupling"],
                value="",
            )
            with gr.Row():
                cs_fuzzy = gr.Checkbox(label="Fuzzy ligand matching", value=True)
                cs_limit = gr.Slider(label="Limit", minimum=1, maximum=200, value=25, step=1)
            cs_btn = gr.Button("Search", variant="primary")
            cs_json = gr.JSON(label="Summary")
            cs_tbl = gr.Dataframe(label="Matches", interactive=False)
            cs_btn.click(ui_core_search, inputs=[cs_core, cs_family, cs_fuzzy, cs_limit], outputs=[cs_json, cs_tbl])

            gr.Markdown("""
            ### List Available Cores
            List unique condition cores present in the loaded reaction dataset (optionally by family).
            """)
            def ui_list_cores(family: str, limit: int) -> List[List[Any]]:
                fam = (family or "").strip() or None
                rows = precedent.list_cores(family=fam, top_n=int(limit or 200), include_counts=True)  # type: ignore[arg-type]
                header = ["core", "count"]
                table: List[List[Any]] = [header]
                for item in rows:  # type: ignore[assignment]
                    table.append([item.get("core", ""), item.get("count", 0)])  # type: ignore[union-attr]
                return table
            with gr.Row():
                lc_family = gr.Dropdown(
                    label="Family (optional)",
                choices=["", "Ullmann C-N", "Buchwald C-N", "Suzuki_CC", "Amide_Coupling"],
                    value="",
                )
                lc_limit = gr.Slider(label="Top N", minimum=5, maximum=500, value=200, step=5)
                lc_btn = gr.Button("List cores", variant="secondary")
            lc_tbl = gr.Dataframe(label="Cores", interactive=False)
            lc_btn.click(ui_list_cores, inputs=[lc_family, lc_limit], outputs=[lc_tbl])

    return demo


if __name__ == "__main__":
    demo = build_demo()
    # Bind to localhost; use share=True if you need a public link
    host = os.environ.get("GRADIO_SERVER_NAME", "127.0.0.1")
    try:
        port = int(os.environ.get("GRADIO_SERVER_PORT", "7860"))
    except Exception:
        port = 7860
    demo.launch(server_name=host, server_port=port)

