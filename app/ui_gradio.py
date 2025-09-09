# -*- coding: utf-8 -*-
"""
Gradio UI to interactively test chemtools functions.

Run:
  python app/ui_gradio.py

Then open the URL printed by Gradio (default http://127.0.0.1:7860).
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple
import json
import os, sys

# Ensure project root is importable when running as a script
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import gradio as gr

from chemtools import smiles, router, properties, featurizers, recommend, precedent, reaction_similarity as rs
# Optional role-aware single-molecule featurizer
try:
    from chem_feats import featurize_mol as role_feat_mol  # type: ignore
except Exception:
    role_feat_mol = None  # type: ignore


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
    return smiles.normalize(smi or "")


def ui_detect_family(reactants_text: str) -> Dict[str, Any]:
    # Split reactants by '.' (common SMILES fragment separator) or newlines
    t = (reactants_text or "").strip()
    if not t:
        return {"family": "Unknown", "confidence": 0.0, "hits": {}}
    if "\n" in t:
        reactants = [s.strip() for s in t.splitlines() if s.strip()]
    else:
        reactants = [s.strip() for s in t.split(".") if s.strip()]
    return router.detect_family(reactants)


def ui_featurize_molecular(electrophile: str, nucleophile: str) -> Dict[str, Any]:
    return featurizers.molecular.featurize(electrophile or "", nucleophile or "")


def ui_properties_lookup(query: str) -> Dict[str, Any]:
    return properties.lookup(query or "")


def ui_recommend(
    reaction: str,
    k: int,
    relax_json: str,
    constraints_json: str,
) -> Dict[str, Any]:
    relax = _safe_json_loads(relax_json)
    constraints = _safe_json_loads(constraints_json)
    return recommend.recommend_from_reaction(
        reaction=reaction or "",
        k=int(k or 25),
        relax=relax,
        constraint_rules=constraints,
    )


def ui_design_plate(
    reaction: str,
    plate_size: int,
    relax_json: str,
    constraints_json: str,
) -> Tuple[str, List[List[Any]], Dict[str, Any]]:
    relax = _safe_json_loads(relax_json)
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
) -> Tuple[Dict[str, Any], List[List[Any]]]:
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

    pack = precedent.knn(fam, feat, k=int(k or 25), relax=relax)

    # Build neighbors table (header + rows)
    header = [
        "reaction_id",
        "reaction_smiles",
        "reactants_smiles",
        "products_smiles",
        "yield",
        "core",
        "base_uid",
        "solvent_uid",
        "T_C",
        "time_h",
    ]
    precs = list(pack.get("precedents") or [])
    table: List[List[Any]] = [header]
    for p in precs:
        table.append([
            p.get("reaction_id", ""),
            p.get("reaction_smiles", ""),
            p.get("reactants_smiles", ""),
            p.get("products_smiles", ""),
            p.get("yield", ""),
            p.get("core", ""),
            p.get("base_uid", ""),
            p.get("solvent_uid", ""),
            p.get("T_C", ""),
            p.get("time_h", ""),
        ])
    return pack, table


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


def build_demo() -> gr.Blocks:
    with gr.Blocks(title="ChemTools UI") as demo:
        gr.Markdown("""
        # ChemTools - Interactive UI
        Try common chemistry tools without writing code. Use the tabs below.
        """)

        with gr.Tab("SMILES Normalize"):
            smi_in = gr.Textbox(label="SMILES", value="c1ccccc1O")
            smi_btn = gr.Button("Normalize")
            smi_out = gr.JSON(label="Result")
            smi_btn.click(ui_normalize_smiles, inputs=[smi_in], outputs=[smi_out])

        with gr.Tab("Detect Family"):
            gr.Markdown("Enter reactants separated by '.' or new lines.")
            react_in = gr.Textbox(label="Reactants", value="Clc1ccccc1.Nc1ccccc1")
            react_btn = gr.Button("Detect")
            react_out = gr.JSON(label="Family Result")
            react_btn.click(ui_detect_family, inputs=[react_in], outputs=[react_out])

        with gr.Tab("Molecular Featurizer"):
            elec_in = gr.Textbox(label="Electrophile", value="Clc1ccccc1")
            nuc_in = gr.Textbox(label="Nucleophile", value="Nc1ccccc1")
            feat_btn = gr.Button("Featurize")
            feat_out = gr.JSON(label="Features")
            feat_btn.click(ui_featurize_molecular, inputs=[elec_in, nuc_in], outputs=[feat_out])

        with gr.Tab("Single Molecule (basic)"):
            smi_in = gr.Textbox(label="SMILES", value="Clc1ccccc1")
            roles_in = gr.CheckboxGroup(
                choices=["amine", "alcohol", "aryl_halide"],
                label="Roles (optional; leave empty for globals-only)",
                value=[],
            )
            show_full = gr.Checkbox(value=False, label="Show full fields list (unchecked shows preview)")
            with gr.Row():
                smi_btn = gr.Button("Featurize molecule")
                dl_btn = gr.DownloadButton(label="Download CSV")
            smi_out = gr.JSON(label="Result (vector length, masks, fields)")
            csv_code = gr.Code(label="CSV preview")

            def _ui_single_molecule(smi: str, roles: list[str], show_full_fields: bool) -> dict:
                if role_feat_mol is None:
                    return {"error": "role-aware featurizer unavailable"}
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
                return out

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

            smi_btn.click(_ui_single_molecule, inputs=[smi_in, roles_in, show_full], outputs=[smi_out])
            # Update CSV preview and download
            def _update_csv_preview(smi: str, roles: list[str]) -> str:
                return _to_csv(smi, roles)
            smi_btn.click(_update_csv_preview, inputs=[smi_in, roles_in], outputs=[csv_code])
            dl_btn.click(_to_csv, inputs=[smi_in, roles_in], outputs=[dl_btn])

        with gr.Tab("Properties Lookup"):
            prop_in = gr.Textbox(label="Query (name, CAS, token)", value="Water")
            prop_btn = gr.Button("Lookup")
            prop_out = gr.JSON(label="Record")
            prop_btn.click(ui_properties_lookup, inputs=[prop_in], outputs=[prop_out])

        with gr.Tab("Recommend Conditions"):
            rec_in = gr.Textbox(label="Reaction SMILES", value="Brc1ccccc1.Nc1ccccc1>>")
            rec_k = gr.Slider(label="k (neighbors)", minimum=5, maximum=100, value=25, step=1)
            rec_relax = gr.Textbox(label="Relax (JSON)", value="")
            rec_constraints = gr.Textbox(label="Constraints (JSON)", value="")
            rec_btn = gr.Button("Recommend")
            rec_out = gr.JSON(label="Recommendation Pack")
            rec_btn.click(ui_recommend, inputs=[rec_in, rec_k, rec_relax, rec_constraints], outputs=[rec_out])

        with gr.Tab("Design Plate"):
            plate_in = gr.Textbox(label="Reaction SMILES", value="Brc1ccccc1.Nc1ccccc1>>")
            plate_n = gr.Slider(label="Plate size", minimum=6, maximum=96, value=24, step=1)
            plate_relax = gr.Textbox(label="Relax (JSON)", value="")
            plate_constraints = gr.Textbox(label="Constraints (JSON)", value="")
            plate_btn = gr.Button("Design")
            plate_csv = gr.Code(label="CSV")
            plate_tbl = gr.Dataframe(label="Preview", interactive=False)
            plate_meta = gr.JSON(label="Meta")
            plate_btn.click(
                ui_design_plate,
                inputs=[plate_in, plate_n, plate_relax, plate_constraints],
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
            ps_btn = gr.Button("Search")
            ps_pack = gr.JSON(label="Pack (prototype, support, precedents)")
            ps_tbl = gr.Dataframe(label="Top precedents", interactive=False)
            ps_btn.click(
                ui_precedent_search,
                inputs=[ps_in, ps_k, ps_use_drfp, ps_drfp_w, ps_bits, ps_radius, ps_prec_scope],
                outputs=[ps_pack, ps_tbl],
            )

        with gr.Tab("DRFP Similarity"):
            s_q = gr.Textbox(label="Query reaction SMILES", value="Brc1ccccc1.Nc1ccccc1>>")
            s_r = gr.Textbox(label="Reference reaction SMILES", value="Clc1ccccc1.Nc1ccccc1>>")
            s_bits = gr.Number(label="DRFP bits", value=4096, precision=0)
            s_radius = gr.Number(label="DRFP radius", value=3, precision=0)
            s_btn = gr.Button("Compute Tanimoto")
            s_out = gr.JSON(label="Result")
            s_btn.click(ui_similarity_tanimoto, inputs=[s_q, s_r, s_bits, s_radius], outputs=[s_out])

    return demo


if __name__ == "__main__":
    demo = build_demo()
    # Bind to localhost; use share=True if you need a public link
    demo.launch(server_name="127.0.0.1", server_port=7860)
