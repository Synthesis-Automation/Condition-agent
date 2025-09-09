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

from chemtools import smiles, router, properties, featurizers, recommend


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


def ui_featurize_ullmann(electrophile: str, nucleophile: str) -> Dict[str, Any]:
    return featurizers.ullmann.featurize(electrophile or "", nucleophile or "")


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


def build_demo() -> gr.Blocks:
    with gr.Blocks(title="ChemTools UI") as demo:
        gr.Markdown("""
        # ChemTools – Interactive UI
        Try common chemistry tools without writing code. Use the tabs below.
        """)

        with gr.Tab("SMILES Normalize"):
            smi_in = gr.Textbox(label="SMILES", value="c1ccccc1O")
            smi_btn = gr.Button("Normalize")
            smi_out = gr.JSON(label="Result")
            smi_btn.click(ui_normalize_smiles, inputs=[smi_in], outputs=[smi_out])

        with gr.Tab("Detect Family"):
            gr.Markdown("""Enter reactants separated by '.' or new lines.""")
            react_in = gr.Textbox(label="Reactants", value="Clc1ccccc1.Nc1ccccc1")
            react_btn = gr.Button("Detect")
            react_out = gr.JSON(label="Family Result")
            react_btn.click(ui_detect_family, inputs=[react_in], outputs=[react_out])

        with gr.Tab("Ullmann Featurizer"):
            elec_in = gr.Textbox(label="Electrophile", value="Clc1ccccc1")
            nuc_in = gr.Textbox(label="Nucleophile", value="Nc1ccccc1")
            feat_btn = gr.Button("Featurize")
            feat_out = gr.JSON(label="Features")
            feat_btn.click(ui_featurize_ullmann, inputs=[elec_in, nuc_in], outputs=[feat_out])

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

    return demo


if __name__ == "__main__":
    demo = build_demo()
    # Bind to localhost; use share=True if you need a public link
    demo.launch(server_name="127.0.0.1", server_port=7860)
