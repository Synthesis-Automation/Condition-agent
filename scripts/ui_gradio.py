from __future__ import annotations

import os
import sys
from typing import Any, Dict, List
import tempfile
from pathlib import Path
import json

import gradio as gr

# Theme: emphasize action buttons with a distinct primary hue
THEME = gr.themes.Soft(primary_hue="indigo", secondary_hue="amber", neutral_hue="slate")

# Ensure project root on sys.path for local execution
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Ensure RDKit enabled by default for UI runs (respect explicit disable via env)
os.environ.setdefault("CHEMTOOLS_DISABLE_RDKIT", os.environ.get("CHEMTOOLS_DISABLE_RDKIT", "0") or "0")

# Ensure dataset auto-load for UI runs (respects explicit disable via env)
os.environ.setdefault("CHEMTOOLS_LOAD_DATASET", os.environ.get("CHEMTOOLS_LOAD_DATASET", "1"))

from chemtools.recommend import recommend_from_reaction, design_plate_from_reaction


def _build_relax(use_drfp: bool, drfp_weight: float) -> Dict[str, Any]:
    return {
        "use_drfp": bool(use_drfp),
        "drfp_weight": float(drfp_weight),
        "precompute_drfp": True,
        "precompute_scope": "candidates",
    }


def _build_constraints(no_chloro: bool, no_hmpa: bool, aqueous_only: bool, min_bp: float | None) -> Dict[str, Any]:
    rules: Dict[str, Any] = {
        "no_chlorinated": bool(no_chloro),
        "no_HMPA": bool(no_hmpa),
        "aqueous_only": bool(aqueous_only),
    }
    try:
        if min_bp is not None and float(min_bp) > 0:
            rules["min_bp_C"] = float(min_bp)
    except Exception:
        pass
    return rules


def ui_recommend(reaction: str, k: int, use_drfp: bool, drfp_weight: float,
                 no_chloro: bool, no_hmpa: bool, aqueous_only: bool, min_bp: float | None):
    relax = _build_relax(use_drfp, drfp_weight)
    rules = _build_constraints(no_chloro, no_hmpa, aqueous_only, min_bp)
    out = recommend_from_reaction(reaction, k=int(k), relax=relax, constraint_rules=rules)
    rec = out.get("recommendation") or {}
    alt = out.get("alternatives") or {}
    reasons = out.get("reasons") or []
    # Format a compact summary
    summary = {
        "core": rec.get("core"),
        "base_uid": rec.get("base_uid"),
        "solvent_uid": rec.get("solvent_uid"),
        "T_C": rec.get("T_C"),
        "time_h": rec.get("time_h"),
        "confidence": rec.get("confidence"),
        "bin": out.get("bin"),
        "family": out.get("family"),
        "precedent_support": (out.get("precedent_pack") or {}).get("support"),
    }
    # Render top alternatives as markdown
    md = []
    if reasons:
        md.append("\n".join(f"- {r}" for r in reasons))
    if alt:
        cores = ", ".join(f"{c} ({n})" for c, n in (alt.get("cores") or []) if c)
        bases = ", ".join(f"{b} ({n})" for b, n in (alt.get("bases") or []) if b)
        solv = ", ".join(f"{s} ({n})" for s, n in (alt.get("solvents") or []) if s)
        md.append(f"Alternatives — cores: {cores or 'n/a'}; bases: {bases or 'n/a'}; solvents: {solv or 'n/a'}")
    json_text = json.dumps(out, ensure_ascii=False, indent=2)
    return summary, "\n\n".join(md), json_text


def ui_plate(reaction: str, use_drfp: bool, drfp_weight: float,
             no_chloro: bool, no_hmpa: bool, aqueous_only: bool, min_bp: float | None):
    relax = _build_relax(use_drfp, drfp_weight)
    rules = _build_constraints(no_chloro, no_hmpa, aqueous_only, min_bp)
    out = design_plate_from_reaction(reaction, plate_size=24, relax=relax, constraint_rules=rules)
    rows = out.get("rows") or []
    csv_text = out.get("csv") or ""
    return rows, csv_text, csv_text


def ui_download_csv(csv_text: str | None):
    txt = csv_text or ""
    if not txt.strip():
        return None
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".csv", mode="w", encoding="utf-8") as f:
            f.write(txt)
            path = f.name
        return Path(path)
    except Exception:
        return None


def ui_download_json(json_text: str | None):
    txt = json_text or ""
    if not txt.strip():
        return None
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".json", mode="w", encoding="utf-8") as f:
            f.write(txt)
            path = f.name
        return Path(path)
    except Exception:
        return None


def ui_precedents(reaction: str, k: int, use_drfp: bool, drfp_weight: float):
    relax = _build_relax(use_drfp, drfp_weight)
    out = recommend_from_reaction(reaction, k=int(k), relax=relax, constraint_rules={})
    pack = out.get("precedent_pack") or {}
    precs = list(pack.get("precedents") or [])
    # Build rows for the Dataframe and CSV
    header = ["reaction_id","yield","core","base_uid","solvent_uid","T_C","time_h"]
    rows = []
    for p in precs:
        rows.append([
            p.get("reaction_id"),
            p.get("yield"),
            p.get("core"),
            p.get("base_uid"),
            p.get("solvent_uid"),
            p.get("T_C"),
            p.get("time_h"),
        ])
    def esc(v: Any) -> str:
        s = "" if v is None else str(v)
        return ('"' + s.replace('"','""') + '"') if ("," in s or '"' in s) else s
    lines = [",".join(header)]
    for r in rows:
        lines.append(
            ",".join(esc(x) for x in r)
        )
    csv_text = "\n".join(lines) + ("\n" if lines else "")
    json_text = json.dumps(precs, ensure_ascii=False, indent=2)
    # Return a 2D list for the Dataframe component instead of list-of-dicts
    return rows, csv_text, json_text


with gr.Blocks(title="Condition Recommender", theme=THEME) as demo:
    gr.Markdown("""
    # Condition Recommendation (Ullmann-first)
    Enter a reaction SMILES (reactants>>products). The app recommends Core/Base/Solvent/T/time and can design a 24‑well plate diversified across cores.
    """)

    with gr.Tab("Recommend"):
        rxn = gr.Textbox(label="Reaction SMILES", placeholder="Brc1ccc(F)cc1.Nc1ccccc1>>Nc1ccc(F)cc1")
        k = gr.Slider(5, 200, value=25, step=1, label="Neighbors (k)")
        with gr.Row():
            use_drfp = gr.Checkbox(True, label="Use DRFP re‑ranking")
            drfp_w = gr.Slider(0.0, 1.0, value=0.4, step=0.05, label="DRFP weight")
        gr.Markdown("Constraints")
        with gr.Row():
            no_chloro = gr.Checkbox(False, label="No chlorinated solvents")
            no_hmpa = gr.Checkbox(True, label="No HMPA")
            aq_only = gr.Checkbox(False, label="Aqueous only")
            min_bp = gr.Slider(0, 200, value=0, step=5, label="Min solvent bp (C)")
        btn = gr.Button("Recommend", variant="primary")
        summary = gr.JSON(label="Recommendation")
        reasons = gr.Markdown()
        json_state = gr.State("")
        btn.click(ui_recommend, [rxn, k, use_drfp, drfp_w, no_chloro, no_hmpa, aq_only, min_bp], [summary, reasons, json_state])
        with gr.Row():
            btn_json = gr.Button("Download JSON", variant="secondary")
            file_json = gr.File(label="recommendation.json")
        btn_json.click(ui_download_json, [json_state], [file_json])

    with gr.Tab("Design 24‑well Plate"):
        rxn2 = gr.Textbox(label="Reaction SMILES", placeholder="Brc1ccc(F)cc1.Nc1ccccc1>>Nc1ccc(F)cc1")
        with gr.Row():
            use_drfp2 = gr.Checkbox(True, label="Use DRFP re‑ranking")
            drfp_w2 = gr.Slider(0.0, 1.0, value=0.4, step=0.05, label="DRFP weight")
        gr.Markdown("Constraints")
        with gr.Row():
            no_chloro2 = gr.Checkbox(False, label="No chlorinated solvents")
            no_hmpa2 = gr.Checkbox(True, label="No HMPA")
            aq_only2 = gr.Checkbox(False, label="Aqueous only")
            min_bp2 = gr.Slider(0, 200, value=0, step=5, label="Min solvent bp (C)")
        btn2 = gr.Button("Design Plate", variant="primary")
        grid = gr.Dataframe(headers=["well_id","core","base_uid","solvent_uid","additive_uids","T_C","time_h"],
                            datatype=["str","str","str","str","str","str","str"],
                            interactive=False, wrap=True)
        csv_box = gr.Textbox(label="CSV", lines=10)
        csv_state = gr.State("")
        btn2.click(ui_plate, [rxn2, use_drfp2, drfp_w2, no_chloro2, no_hmpa2, aq_only2, min_bp2], [grid, csv_box, csv_state])
        with gr.Row():
            btn_dl = gr.Button("Download CSV", variant="secondary")
            file_out = gr.File(label="plate.csv")
        btn_dl.click(ui_download_csv, [csv_state], [file_out])

    with gr.Tab("Precedents"):
        rxn3 = gr.Textbox(label="Reaction SMILES", placeholder="Brc1ccc(F)cc1.Nc1ccccc1>>Nc1ccc(F)cc1")
        k3 = gr.Slider(5, 200, value=25, step=1, label="Neighbors (k)")
        with gr.Row():
            use_drfp3 = gr.Checkbox(True, label="Use DRFP re‑ranking")
            drfp_w3 = gr.Slider(0.0, 1.0, value=0.4, step=0.05, label="DRFP weight")
        btn3 = gr.Button("Fetch Precedents", variant="primary")
        tbl = gr.Dataframe(headers=["reaction_id","yield","core","base_uid","solvent_uid","T_C","time_h"],
                           datatype=["str","number","str","str","str","number","number"],
                           interactive=False, wrap=True)
        prec_csv_state = gr.State("")
        prec_json_state = gr.State("")
        btn3.click(ui_precedents, [rxn3, k3, use_drfp3, drfp_w3], [tbl, prec_csv_state, prec_json_state])
        with gr.Row():
            btn_prec_csv = gr.Button("Download CSV", variant="secondary")
            file_prec_csv = gr.File(label="precedents.csv")
            btn_prec_json = gr.Button("Download JSON", variant="secondary")
            file_prec_json = gr.File(label="precedents.json")
        btn_prec_csv.click(ui_download_csv, [prec_csv_state], [file_prec_csv])
        btn_prec_json.click(ui_download_json, [prec_json_state], [file_prec_json])


if __name__ == "__main__":
    # Launch on http://127.0.0.1:7860
    demo.launch()
