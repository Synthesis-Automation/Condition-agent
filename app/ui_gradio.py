# -*- coding: utf-8 -*-
"""
Gradio UI to interactively test chemtools functions.

Run:
  python app/ui_gradio.py

Then open the URL printed by Gradio (default http://127.0.0.1:7860).
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple, Optional
from functools import lru_cache
from types import SimpleNamespace
import json
import os, sys
from pathlib import Path

# Ensure project root is importable when running as a script
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import gradio as gr

from condition_mcp import DEFAULT_SCHEMA_VERSION
from condition_mcp.tools import (
    detect_family as mcp_detect_family,
    featurize_substrates as mcp_featurize_substrates,
    normalize_reaction as mcp_normalize_reaction,
)

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

# Unified mapping between UI labels, structured datasets, and rule families.
# Each entry ensures the Structured (ML) and Rule-based recommenders operate on
# the same reaction family definition.
RECOMMEND_REACTION_TYPE_CONFIGS: List[Dict[str, str]] = [
    {
        "label": "Auto-detect (default)",
        "structured_family": "",
        "rule_family": "",
        "dataset_label": "Auto-detect via router",
    },
    {
        "label": "Suzuki (Suzuki_CC)",
        "structured_family": "Suzuki_CC",
        "rule_family": "Suzuki_Coupling",
        "dataset_label": "Suzuki_CC",
    },
    {
        "label": "Ullmann C–N",
        "structured_family": "Ullmann C?N",
        "rule_family": "Ullmann_CN",
        "dataset_label": "Ullmann C?N",
    },
    {
        "label": "Buchwald C–N",
        "structured_family": "Buchwald C?N",
        "rule_family": "Buchwald_CN",
        "dataset_label": "Buchwald C?N",
    },
    {
        "label": "Amide Coupling",
        "structured_family": "Amide_Coupling",
        "rule_family": "Amide_Formation",
        "dataset_label": "Amide_Coupling",
    },
]

RECOMMEND_FAMILY_OPTIONS = [cfg["label"] for cfg in RECOMMEND_REACTION_TYPE_CONFIGS]
RECOMMEND_FAMILY_VALUE_MAP = {cfg["label"]: cfg["structured_family"] for cfg in RECOMMEND_REACTION_TYPE_CONFIGS}
RECOMMEND_RULE_FAMILY_MAP = {cfg["label"]: cfg["rule_family"] for cfg in RECOMMEND_REACTION_TYPE_CONFIGS}
RECOMMEND_DATASET_LABEL_MAP = {cfg["label"]: cfg.get("dataset_label", cfg["structured_family"]) for cfg in RECOMMEND_REACTION_TYPE_CONFIGS}
RECOMMEND_REACTION_TYPE_CONFIG_MAP = {cfg["label"]: cfg for cfg in RECOMMEND_REACTION_TYPE_CONFIGS}

# Catalyst class options (optional filter)
CATALYST_CLASS_DISPLAY = [
    ("Any (no filter)", ""),
    ("Pd (palladium)", "Pd"),
    ("Ni (nickel)", "Ni"),
    ("Co (cobalt)", "Co"),
    ("Cu (copper)", "Cu"),
    ("Ru (ruthenium)", "Ru"),
    ("Rh (rhodium)", "Rh"),
    ("Ir (iridium)", "Ir"),
    ("Fe (iron)", "Fe"),
    ("Ag (silver)", "Ag"),
    ("Au (gold)", "Au"),
    ("Zn (zinc)", "Zn"),
    ("Organo catalyst", "organo_catalyst"),
    ("Enzyme", "enzyme"),
    ("Other", "other"),
]
CATALYST_CLASS_OPTIONS = [label for label, _ in CATALYST_CLASS_DISPLAY]
CATALYST_CLASS_VALUE_MAP = {label: value for label, value in CATALYST_CLASS_DISPLAY}

CONDITION_RULE_FAMILIES: List[str] = []

_RULE_SAMPLE_FEATURES: Dict[str, Dict[str, Any]] = {
    "Buchwald_CN": {
        "electrophile": {
            "class": "aryl chloride",
            "electronics": "electron-poor",
            "ortho_sub_count": 1,
        },
        "nucleophile": {
            "class": "primary aniline",
        },
        "solvent": "DMF",
    },
    "Amide_Formation": {
        "reaction_family": "amide_formation",
        "category": "direct_dehydrative_boron_catalysis",
        "substrate_class": "Aliphatic/aromatic carboxylic acids with primary/secondary amines",
        "electrophile": {
            "class": "carboxylic acid",
            "state": "liquid",
        },
        "nucleophile": {
            "class": "primary amine",
        },
        "solvent": "toluene",
        "water_management": "Dean-Stark",
    },
    "Suzuki_Coupling": {
        "reaction_family": "suzuki_coupling",
        "category": "aryl_halides",
        "substrate_class": "Unhindered aryl iodides and bromides",
        "electrophile": {
            "class": "aryl halide",
            "halide": "Br",
            "electronics": "neutral",
        },
        "boron_partner": {
            "class": "aryl boron",
            "identity": "aryl-Bpin",
        },
        "base": "K2CO3",
        "solvent": "1,4-dioxane/H2O",
    },
    "Ullmann_CN": {
        "reaction_family": "ullmann_cn",
        "electrophile": {
            "class": "aryl bromide or iodide",
            "halide": "Br",
        },
        "nucleophile": {
            "class": "primary aniline",
        },
        "base": "K3PO4",
        "ligand": "N,N-dimethylglycine",
        "cu_source": "CuI",
        "solvent": "DMSO",
        "temperature_C": 90,
    },
}

_RULE_SAMPLE_JOB_CTX: Dict[str, Dict[str, Any]] = {
    "Buchwald_CN": {"base": "DBU", "mode": "batch"},
    "Amide_Formation": {
        "mode": "Dean-Stark",
        "water_removal": "molecular_sieves",
    },
    "Suzuki_Coupling": {
        "mode": "batch",
        "base": "K2CO3",
    },
    "Ullmann_CN": {
        "mode": "batch",
        "temperature_C": 100,
        "base": "K3PO4",
        "cu_source": "CuI",
    },
}

_RULE_SAMPLE_FEATURES_TEXT: Dict[str, str] = {
    key: json.dumps(value, indent=2, sort_keys=True) for key, value in _RULE_SAMPLE_FEATURES.items()
}
_RULE_SAMPLE_JOB_CTX_TEXT: Dict[str, str] = {
    key: json.dumps(value, indent=2, sort_keys=True) for key, value in _RULE_SAMPLE_JOB_CTX.items()
}

# Constrain overall width for a cleaner layout
CSS = """
.gradio-container{max-width:1600px !important;margin:0 auto !important;}
.gradio-container [role='tablist']{flex-wrap:wrap;overflow-x:auto;gap:0.25rem;}
.gradio-container [role='tab']{white-space:nowrap;}
.gradio-container [role='tab'][aria-selected='true']{
    background-color:#ea580c !important;
    color:#fff !important;
    border-color:#ea580c !important;
    box-shadow:0 0 0 1px #ea580c inset;
}
.gradio-container [role='tab'][aria-selected='true'] svg{
    color:#fff !important;
    fill:#fff !important;
}
.gradio-container button:not([role='tab']){
    background-color:#2563eb !important;
    border-color:#1d4ed8 !important;
    color:#fff !important;
    box-shadow:0 1px 2px rgba(15, 23, 42, 0.2);
    transition:background-color 0.15s ease, box-shadow 0.15s ease;
}
.gradio-container button:not([role='tab']):hover{
    background-color:#1d4ed8 !important;
    box-shadow:0 2px 6px rgba(15, 23, 42, 0.25);
}
.gradio-container button:not([role='tab']):focus-visible{
    outline:2px solid #1e3a8a !important;
    outline-offset:2px;
}
"""

CONDITION_SET_SCHEMA_PATH = Path(ROOT) / "condition_mcp" / "resources" / "schemas" / "condition_set.json"

# Enable RDKit by default unless explicitly disabled by the environment
os.environ.setdefault("CHEMTOOLS_DISABLE_RDKIT", "0")
# Prefer rxn-insight enabled by default (can be disabled via env or UI toggle)
os.environ.setdefault("CHEMTOOLS_DISABLE_RXN_INSIGHT", "0")

from chemtools import smiles, router, properties, featurizers, recommend, precedent, reaction_similarity as rs, condition_rules
from chemtools import registry as creg
from chemtools.util.rdkit_helpers import rdkit_available

_CRL_DEFAULT: Dict[str, Any] | None = None

try:
    from condition_rule_library import api as crl_api
except Exception:  # pragma: no cover - optional dependency guard
    crl_api = None  # type: ignore
else:
    try:
        _CRL_DEFAULT = crl_api.load_default_crl()
        CONDITION_RULE_FAMILIES = sorted(list((_CRL_DEFAULT.get("families") or {}).keys()))
    except Exception:  # pragma: no cover - defensive fallback
        _CRL_DEFAULT = None
        CONDITION_RULE_FAMILIES = []

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
# Optional role-aware featurizers
try:
    from chem_feats import featurize_mol as role_feat_mol  # type: ignore
except Exception:
    role_feat_mol = None  # type: ignore

try:
    from chem_feats import featurize_reaction as role_feat_rxn  # type: ignore
except Exception:
    role_feat_rxn = None  # type: ignore


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


def _rule_default_payloads(family: str) -> Tuple[str, str]:
    feat = _RULE_SAMPLE_FEATURES_TEXT.get(family)
    if feat is None:
        feat = json.dumps({"TODO": "provide features"}, indent=2, sort_keys=True)
    job = _RULE_SAMPLE_JOB_CTX_TEXT.get(family)
    if job is None:
        job = json.dumps({}, indent=2)
    return feat, job


def _format_detection_details(det: Dict[str, Any] | None) -> str:
    if not isinstance(det, dict):
        return '- No detection details'
    lines: List[str] = []
    rt = det.get('reaction_type') if isinstance(det, dict) else None
    if rt:
        lines.append(f'Detected family: {rt}')
    source = det.get('source') if isinstance(det, dict) else None
    if source:
        lines.append(f'Source: {source}')
    provided = det.get('provided_reaction_type') if isinstance(det, dict) else None
    if provided:
        lines.append(f'Provided reaction type: {provided}')
    auto = det.get('auto') if isinstance(det, dict) else None
    if isinstance(auto, dict):
        lines.append(f"rxn-insight available: {bool(auto.get('rxn_insight_available'))}")
        status = auto.get('status')
        if status:
            lines.append(f"rxn-insight status: {status}")
        auto_rt = auto.get('reaction_type')
        if auto_rt:
            lines.append(f"rxn-insight family: {auto_rt}")
        auto_name = auto.get('rxn_insight_name')
        if auto_name:
            lines.append(f"rxn-insight name: {auto_name}")
        auto_class = auto.get('rxn_insight_class')
        if auto_class:
            lines.append(f"rxn-insight class: {auto_class}")
        auto_conf = auto.get('rxn_insight_confidence')
        if auto_conf is not None:
            try:
                lines.append(f"rxn-insight confidence: {float(auto_conf):.3f}")
            except Exception:
                lines.append(f"rxn-insight confidence: {auto_conf}")
    rule = det.get('rule_based') if isinstance(det, dict) else None
    if isinstance(rule, dict):
        rule_rt = rule.get('reaction_type')
        if rule_rt:
            lines.append(f"Rule-based family: {rule_rt}")
    if not lines:
        return '- No detection details'
    return '\n'.join(f'- {line}' for line in lines)


def _format_starting_materials_summary(starting: Dict[str, Any] | None) -> str:
    if not isinstance(starting, dict):
        return '- Starting-material featurization unavailable.'

    lines: List[str] = []
    role_pack = starting.get('role_featurization')
    if isinstance(role_pack, dict):
        assignments = role_pack.get('assignments') or {}
        reactants = role_pack.get('reactants') or []
        for label, info in assignments.items():
            idx = info.get('index')
            smi = None
            if isinstance(idx, int) and 0 <= idx < len(reactants):
                smi = reactants[idx].get('smiles')
            role_token = info.get('role') or 'unspecified'
            reason = 'mask' if info.get('matched_via_mask') else 'fallback'
            summary = smi or 'n/a'
            lines.append(f"- {label}: {summary} ({role_token}; {reason})")
        if not assignments:
            lines.append('- Role-aware vectors present but no assignments matched.')
    else:
        lines.append('- Role-aware vectors unavailable (install `chem-feats` to enable).')

    rule_feats = starting.get('rule_features')
    if isinstance(rule_feats, dict):
        elec = rule_feats.get('electrophile') or {}
        nuc = rule_feats.get('nucleophile') or {}
        if rule_feats.get('reaction_family'):
            lines.append(f"- Reaction family: {rule_feats.get('reaction_family')}")
        if rule_feats.get('category'):
            lines.append(f"- Rule category: {rule_feats.get('category')}")
        if rule_feats.get('substrate_class'):
            lines.append(f"- Substrate class: {rule_feats.get('substrate_class')}")
        lines.append(f"- Electrophile class: {elec.get('class') or 'unknown'}")
        if elec.get('description'):
            lines.append(f"  • detail: {elec.get('description')}")
        if elec.get('alpha_branching'):
            lines.append(f"  • alpha branching: {elec.get('alpha_branching')}")
        if elec.get('state'):
            lines.append(f"  • state: {elec.get('state')}")
        if elec.get('smiles'):
            lines.append(f"  • electrophile smiles: {elec.get('smiles')}")
        lines.append(f"- Nucleophile class: {nuc.get('class') or 'unknown'}")
        if nuc.get('basicity'):
            lines.append(f"  • basicity: {nuc.get('basicity')}")
        if nuc.get('steric_alpha'):
            lines.append(f"  • steric alpha: {nuc.get('steric_alpha')}")
        if nuc.get('smiles'):
            lines.append(f"  • nucleophile smiles: {nuc.get('smiles')}")
        if rule_feats.get('ortho_count') is not None:
            lines.append(f"- ortho substitution: {rule_feats.get('ortho_count')}")
        if rule_feats.get('para_EWG') is not None:
            lines.append(f"- para EWG: {rule_feats.get('para_EWG')}")
        if rule_feats.get('n_basicity') is not None:
            lines.append(f"- nucleophile basicity: {rule_feats.get('n_basicity')}")
        if rule_feats.get('water_management'):
            lines.append(f"- Water management: {rule_feats.get('water_management')}")
    else:
        lines.append('- Rule-based substrate features missing from response.')

    return '\n'.join(lines)


_STARTING_TABLE_HEADERS: List[str] = [
    "reactant",
    "smiles",
    "assigned_roles",
    "amine.present",
    "amine.class_ps3",
    "amine.aniline_flag",
    "amine.alpha_branch",
    "amine.h_count_on_N",
    "aryl_halide.present",
    "aryl_halide.halide",
    "aryl_halide.ortho_block",
    "aryl_halide.ipso_degree",
    "aryl_halide.para_EWG",
    "aryl_halide.heteroaryl",
]


def _starting_materials_table(starting: Dict[str, Any] | None) -> List[List[Any]]:
    rows: List[List[Any]] = []
    if not isinstance(starting, dict):
        return rows

    role_pack = starting.get('role_featurization') if isinstance(starting.get('role_featurization'), dict) else None
    if not isinstance(role_pack, dict):
        return rows

    reactants = role_pack.get('reactants')
    if not isinstance(reactants, list):
        return rows

    assignments = role_pack.get('assignments') if isinstance(role_pack.get('assignments'), dict) else {}
    idx_roles: Dict[int, List[str]] = {}
    if isinstance(assignments, dict):
        for info in assignments.values():
            if not isinstance(info, dict):
                continue
            idx = info.get('index')
            role = info.get('role') or info.get('label')
            if isinstance(idx, int):
                idx_roles.setdefault(idx, []).append(str(role or ''))

    feature_keys = [
        'amine.present',
        'amine.class_ps3',
        'amine.aniline_flag',
        'amine.alpha_branch',
        'amine.h_count_on_N',
        'aryl_halide.present',
        'aryl_halide.halide',
        'aryl_halide.ortho_block',
        'aryl_halide.ipso_degree',
        'aryl_halide.para_EWG',
        'aryl_halide.heteroaryl',
    ]

    for idx, react in enumerate(reactants):
        if not isinstance(react, dict):
            continue
        smi = react.get('smiles') or react.get('input') or ''
        fields = react.get('fields') or []
        vector = react.get('vector')
        if hasattr(vector, 'tolist'):
            try:
                vector = vector.tolist()
            except Exception:
                vector = list(vector)
        if not isinstance(vector, list):
            try:
                vector = list(vector or [])
            except Exception:
                vector = []
        field_map: Dict[str, Any] = {}
        for name, value in zip(fields, vector):
            if isinstance(name, str):
                field_map[name] = value

        amine_present = react.get('masks', {}).get('amine') if isinstance(react.get('masks'), dict) else None
        aryl_present = react.get('masks', {}).get('aryl_halide') if isinstance(react.get('masks'), dict) else None

        row: List[Any] = [
            idx,
            smi,
            ', '.join(idx_roles.get(idx, [])) or '',
        ]

        for key in feature_keys:
            if key == 'amine.present' and amine_present is not None:
                row.append(amine_present)
                continue
            if key == 'aryl_halide.present' and aryl_present is not None:
                row.append(aryl_present)
                continue
            row.append(field_map.get(key))

        rows.append(row)

    if not rows:
        rows.append([None] + [''] * (len(_STARTING_TABLE_HEADERS) - 1))

    return rows


def _rule_feature_component_updates(features_obj: Dict[str, Any] | None, fallback_family: str | None = None) -> Tuple[Any, Any, Any]:
    if not isinstance(features_obj, dict) or not features_obj:
        return gr.update(), gr.update(), gr.update()
    try:
        features_pretty = json.dumps(features_obj, indent=2, sort_keys=True)
    except Exception:
        features_pretty = json.dumps(features_obj, indent=2)
    try:
        features_compact = json.dumps(features_obj, sort_keys=True)
    except Exception:
        features_compact = json.dumps(features_obj)
    family_val = str(features_obj.get("family") or fallback_family or "").strip()
    family_update = gr.update(value=family_val) if family_val else gr.update()
    return family_update, gr.update(value=features_pretty), gr.update(value=features_compact)


def _rule_payload_from_starting(starting: Dict[str, Any] | None) -> Tuple[Any, Any, Any]:
    if not isinstance(starting, dict):
        return gr.update(), gr.update(), gr.update()
    features_obj = starting.get("rule_features")
    fallback_family = None
    if isinstance(features_obj, dict):
        if isinstance(features_obj.get("family"), str):
            fallback_family = features_obj.get("family")
    return _rule_feature_component_updates(features_obj, fallback_family)


CONDITION_RULE_FAMILY_SET = set(CONDITION_RULE_FAMILIES)


def _get_reaction_type_config(label: str) -> Dict[str, str]:
    return RECOMMEND_REACTION_TYPE_CONFIG_MAP.get(label or "", {})


def ui_recommend_reaction_type_meta(label: str) -> Tuple[str, str, str, str]:
    cfg = _get_reaction_type_config(label)
    structured_family = cfg.get("structured_family", "")
    dataset_label = cfg.get("dataset_label") or structured_family or "Auto-detect via router"
    rule_family = cfg.get("rule_family", "")

    summary_lines: List[str] = []
    if structured_family:
        if dataset_label and dataset_label != structured_family:
            summary_lines.append(f"- **Structured dataset**: `{structured_family}` ({dataset_label})")
        else:
            summary_lines.append(f"- **Structured dataset**: `{structured_family}`")
    else:
        summary_lines.append(f"- **Structured dataset**: {dataset_label}")

    if rule_family:
        availability = "available" if rule_family in CONDITION_RULE_FAMILY_SET else "missing"
        summary_lines.append(f"- **Rule family**: `{rule_family}` ({availability})")
    else:
        summary_lines.append("- **Rule family**: none (rule recommender disabled)")

    if structured_family and rule_family:
        summary_lines.append(f"- Mapping: `{structured_family}` ↔ `{rule_family}`")

    features_text, job_ctx_text = "{}", "{}"
    if rule_family:
        features_text, job_ctx_text = _rule_default_payloads(rule_family)

    return "\n".join(summary_lines), rule_family, features_text, job_ctx_text


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


def _split_reactant_block(text: str) -> List[str]:
    t = (text or "").replace("\r\n", "\n").strip()
    if not t:
        return []
    if "\n" in t:
        return [s.strip() for s in t.splitlines() if s.strip()]
    return [s.strip() for s in t.split(".") if s.strip()]


def ui_mcp_normalize(smiles_rxn: str, include_agents: bool) -> Tuple[Dict[str, Any], str]:
    rsmi = (smiles_rxn or "").strip()
    if not rsmi:
        return {"error": "Provide a reaction SMILES string."}, ""
    try:
        out = mcp_normalize_reaction({"smiles_rxn": rsmi, "include_agents": bool(include_agents)})
    except Exception as exc:  # pragma: no cover - UI surface
        return {"error": str(exc)}, ""
    reactants = out.get("reactants") if isinstance(out, dict) else []
    reactant_text = "\n".join(str(r) for r in (reactants or []) if r)
    return out, reactant_text


def ui_mcp_detect(reactants_text: str) -> Dict[str, Any]:
    reactants = _split_reactant_block(reactants_text)
    if not reactants:
        return {"error": "Provide one or more reactant SMILES (separate with new lines or '.')."}
    try:
        return mcp_detect_family({"reactants": reactants})
    except Exception as exc:  # pragma: no cover
        return {"error": str(exc)}


def ui_mcp_featurize(reactants_text: str) -> Dict[str, Any]:
    reactants = _split_reactant_block(reactants_text)
    if not reactants:
        return {"error": "Provide reactant SMILES first (one per line or separated by '.')."}
    try:
        return mcp_featurize_substrates({"reactants": reactants})
    except Exception as exc:  # pragma: no cover
        return {"error": str(exc)}


def ui_mcp_condition_set_schema() -> Dict[str, Any]:
    try:
        text = CONDITION_SET_SCHEMA_PATH.read_text(encoding="utf-8")
    except FileNotFoundError:
        return {
            "error": "ConditionSet schema not found",
            "path": str(CONDITION_SET_SCHEMA_PATH),
        }
    except Exception as exc:  # pragma: no cover
        return {"error": f"Failed to read schema: {exc}"}
    try:
        return json.loads(text)
    except Exception as exc:  # pragma: no cover
        return {"error": f"Failed to parse schema JSON: {exc}"}


def ui_mcp_rule_resource(family: str) -> Dict[str, Any]:
    if not crl_api or not _CRL_DEFAULT:
        return {"error": "Condition Rule Library not available in this environment."}
    fam = (family or "").strip()
    if not fam:
        return {"error": "Choose a rule family."}
    families = _CRL_DEFAULT.get("families") if isinstance(_CRL_DEFAULT, dict) else None
    entry = (families or {}).get(fam)
    if entry is None:
        return {"error": f"No rule entry found for family '{fam}'."}
    return {
        "family": fam,
        "metadata": {k: v for k, v in entry.items() if k != "rules"},
        "rules": entry.get("rules"),
    }


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
    det = (out.get("formatted") or {}).get("detection")
    detection_md = _format_detection_details(det)

    return out, [hdr, row], human, detection_md



def ui_recommend_structured(
    reaction: str,
    reaction_type: Optional[str] = None,
    k: int = 50,
    limit: int = 5,
    catalyst_class: Optional[str] = None,
    use_rxn_insight: bool = True,
    relax_json: str = "",
    constraints_json: str = "",
) -> Tuple[Dict[str, Any], List[List[Any]], str, str]:
    # Back-compat for legacy UI signature without reaction_type/limit.
    if isinstance(reaction_type, (int, float)) and isinstance(k, bool):
        legacy_k = int(reaction_type)
        legacy_use = bool(k)
        legacy_relax = limit if isinstance(limit, str) else ""
        legacy_constraints = use_rxn_insight if isinstance(use_rxn_insight, str) else ""
        reaction_type = None
        k = legacy_k
        limit = 5
        use_rxn_insight = legacy_use
        relax_json = legacy_relax
        constraints_json = legacy_constraints

    family_override = reaction_type
    if isinstance(family_override, str):
        family_override = RECOMMEND_FAMILY_VALUE_MAP.get(family_override, family_override)
        if isinstance(family_override, str) and not family_override.strip():
            family_override = None

    relax = _safe_json_loads(relax_json)
    if not isinstance(relax, dict):
        relax = {}
    # Thread catalyst_class filter (optional)
    if isinstance(catalyst_class, str) and catalyst_class.strip():
        cc_val = CATALYST_CLASS_VALUE_MAP.get(catalyst_class, catalyst_class).strip()
        if cc_val:
            relax["catalyst_class"] = cc_val
    if not bool(use_rxn_insight):
        relax["use_rxn_insight"] = False
    constraints = _safe_json_loads(constraints_json)
    if not isinstance(constraints, dict):
        constraints = {}

    data = recommend.recommend_conditions_structured(
        reaction=reaction or "",
        reaction_type=family_override,
        k=int(k or 50),
        limit=int(limit or 5),
        relax=relax,
        constraints=constraints or None,
    )

    recs = data.get("recommendations") or []
    rows: List[List[Any]] = []
    human_parts: List[str] = []
    for rec in recs:
        summary = rec.get("summary") or {}
        base = summary.get("base") or {}
        solvent = summary.get("solvent") or {}
        support = summary.get("support") or {}
        rows.append([
            rec.get("rank"),
            summary.get("core"),
            base.get("name") or base.get("cas"),
            solvent.get("name") or solvent.get("cas"),
            summary.get("confidence"),
            support.get("count"),
        ])
        label_parts = [
            str(summary.get("core") or ""),
            str(base.get("name") or base.get("cas") or ""),
            str(solvent.get("name") or solvent.get("cas") or ""),
        ]
        label = "/".join(part for part in label_parts if part)
        if label:
            human_parts.append(label)
    human = "; ".join(human_parts) if human_parts else ""
    detection_md = _format_detection_details(data.get("detection"))
    starting_payload = data.get("starting_materials")
    starting_md = _format_starting_materials_summary(starting_payload)
    starting_table = _starting_materials_table(starting_payload)
    return (
        data,
        rows,
        human or "n/a",
        detection_md,
        starting_md,
        starting_table,
        starting_payload or {},
    )


def ui_featurize_reactants_only(
    reaction: str,
    reaction_type: Optional[str] = None,
    use_rxn_insight: bool = True,
    relax_json: str = "",
) -> Tuple[str, str, List[List[Any]], Dict[str, Any]]:
    reaction_text = (reaction or "").strip()
    if not reaction_text:
        placeholder = _starting_materials_table(None)
        return (
            "- Provide a reaction SMILES first.",
            "- Starting-material featurization unavailable.",
            placeholder,
            {},
        )

    fam_override = reaction_type
    if isinstance(fam_override, str):
        fam_override = RECOMMEND_FAMILY_VALUE_MAP.get(fam_override, fam_override)
        if isinstance(fam_override, str) and not fam_override.strip():
            fam_override = None
    relax = _safe_json_loads(relax_json)
    if not isinstance(relax, dict):
        relax = {}
    if "use_rxn_insight" in relax:
        use_rxi = bool(relax.get("use_rxn_insight"))
    else:
        use_rxi = bool(use_rxn_insight)

    try:
        norm = recommend.normalize_reaction(reaction_text)
    except Exception as exc:
        placeholder = _starting_materials_table(None)
        message = f"- Normalization failed: {exc}"
        return (message, "- Starting-material featurization unavailable.", placeholder, {"error": str(exc)})

    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]

    fam_override_clean = (fam_override.strip() if isinstance(fam_override, str) else None)
    rxn_smiles_norm = norm.get("normalized") or reaction_text
    rxn_auto: Dict[str, Any] | None = None
    if use_rxi and getattr(recommend, "_HAS_RXN_INSIGHT", False):
        try:
            rxn_auto = recommend._rxn_detect(rxn_smiles_norm)  # type: ignore[attr-defined]
        except Exception:
            rxn_auto = None
    auto_family = None
    if rxn_auto and rxn_auto.get("success") and rxn_auto.get("mapped_family"):
        auto_family = str(rxn_auto.get("mapped_family") or "Unknown")

    rule_info = recommend.detect_family(reactants)
    rule_family = rule_info.get("family") if isinstance(rule_info, dict) else None

    if fam_override_clean:
        fam = fam_override_clean or auto_family or rule_family or "Unknown"
    else:
        fam = auto_family or rule_family or "Unknown"

    detection = {
        "reaction_type": fam,
        "source": "user_supplied" if fam_override_clean else "auto",
        "provided_reaction_type": fam_override_clean,
        "auto": None,
        "rule_based": {
            "reaction_type": rule_family,
            "raw": rule_info,
        },
    }
    if rxn_auto is not None:
        detection["auto"] = {
            "rxn_insight_available": bool(rxn_auto.get("available")),
            "status": "success" if rxn_auto.get("success") else "fallback",
            "rxn_insight_name": rxn_auto.get("rxn_name"),
            "rxn_insight_class": rxn_auto.get("rxn_class"),
            "rxn_insight_confidence": rxn_auto.get("confidence"),
            "reaction_type": rxn_auto.get("mapped_family"),
        }

    try:
        pick = getattr(recommend, "_pick_electrophile_nucleophile", None)
        if callable(pick):
            elec, nuc = pick(reactants)
        else:
            elec = reactants[0] if reactants else ""
            nuc = reactants[1] if len(reactants) > 1 else ""
        features = recommend.feat_molecular.featurize(elec, nuc)
        role_pack = recommend._role_featurization_for_reactants(fam, reactants)
        rule_features = recommend._compose_rule_features(
            fam,
            features,
            role_pack,
            reactants=reactants,
            detection=detection,
        )
        starting_payload = {
            "role_featurization": role_pack,
            "rule_features": rule_features,
        }
        starting_table = _starting_materials_table(starting_payload)
        starting_md = _format_starting_materials_summary(starting_payload)
    except Exception as exc:
        starting_payload = {"error": str(exc)}
        starting_table = _starting_materials_table(None)
        starting_md = f"- Reactant featurization failed: {exc}"

    detection_md = _format_detection_details(detection)
    return detection_md, starting_md, starting_table, starting_payload
def ui_condition_rule_recommend(
    family: str,
    features_json: str,
    job_ctx_json: str,
    top_n: int,
    *_,
    rule_features_override: Optional[str] = None,
) -> Tuple[Dict[str, Any], List[List[Any]], str]:
    fam = (family or "").strip()
    header = [
        "rank",
        "playbook",
        "name",
        "score",
        "ok",
        "ligand",
        "pd_source",
        "base",
        "solvent",
        "temperature_C",
    ]
    if not fam:
        return {"error": "Select a reaction family."}, [header], "- No family selected -"

    features_payload = rule_features_override if rule_features_override else features_json
    features = _safe_json_loads(features_payload)
    if not isinstance(features, dict) or not features:
        return {
            "error": "Features payload must be a non-empty JSON object.",
            "family": fam,
        }, [header], "- Provide features to evaluate rules -"

    job_ctx = _safe_json_loads(job_ctx_json)

    try:
        recs = condition_rules.recommend_rule_based(
            fam,
            features,
            job_ctx=job_ctx or None,
            top_n=max(1, int(top_n or 3)),
        )
    except Exception as exc:  # pragma: no cover - UI error path
        return {
            "error": str(exc),
            "family": fam,
        }, [header], f"⚠️ {exc}"

    result_payload: Dict[str, Any] = {
        "family": fam,
        "input_features": features,
        "job_ctx": job_ctx or {},
        "recommendations": recs,
    }

    table: List[List[Any]] = [header]
    md_lines: List[str] = []

    defaults: Dict[str, Any] = {}
    try:
        defaults = condition_rules.get_family_defaults(fam)
    except Exception:
        defaults = {}

    if defaults:
        result_payload["family_defaults"] = defaults
        defaults_preview = {
            key: defaults[key]
            for key in list(defaults.keys())[:6]
        }
        md_lines.append("**Family defaults (truncated preview)**")
        md_lines.append("```json")
        md_lines.append(json.dumps(defaults_preview, indent=2, sort_keys=True))
        md_lines.append("```")

    if not recs:
        md_lines.append("No rule-based playbooks matched the provided features.")
        return result_payload, table, "\n".join(md_lines) if md_lines else "No results"

    for idx, rec in enumerate(recs, start=1):
        suggestion = rec.get("suggested_recipe") or {}
        lig = suggestion.get("ligands") or suggestion.get("ligand") or ""
        pd_src = suggestion.get("pd_source") or suggestion.get("pd_source_uid") or ""
        base = suggestion.get("base") or suggestion.get("base_uid") or ""
        solv = suggestion.get("solvent") or suggestion.get("solvent_uid") or ""
        temp = suggestion.get("temperature_C") or suggestion.get("temperature") or ""

        table.append([
            idx,
            rec.get("playbook_id", ""),
            rec.get("name", ""),
            rec.get("score", ""),
            "✅" if rec.get("ok", True) else "⚠️",
            lig,
            pd_src,
            base,
            solv,
            temp,
        ])

        breakdown = rec.get("score_breakdown") or {}
        guard_msgs = rec.get("guard_messages") or []
        md_lines.append(
            f"**{rec.get('playbook_id', rec.get('name', 'Playbook'))}** — score {rec.get('score', '?')} (guards {breakdown.get('guard_adjust', 0)}, priors {breakdown.get('prior', 0)})"
        )
        if guard_msgs:
            for gm in guard_msgs:
                action = str(gm.get("action", "")).upper()
                rationale = gm.get("rationale") or ""
                md_lines.append(f"- `{action}` {rationale}")
        else:
            md_lines.append("- No guard messages.")

    return result_payload, table, "\n".join(md_lines)


def ui_design_plate(
    reaction: str,
    plate_size: int,
    catalyst_class: Optional[str],
    use_rxn_insight: bool,
    relax_json: str,
    constraints_json: str,
) -> Tuple[str, List[List[Any]], Dict[str, Any]]:
    relax = _safe_json_loads(relax_json)
    relax["use_rxn_insight"] = bool(use_rxn_insight)
    if isinstance(catalyst_class, str) and catalyst_class.strip():
        cc_val = CATALYST_CLASS_VALUE_MAP.get(catalyst_class, catalyst_class).strip()
        if cc_val:
            relax["catalyst_class"] = cc_val
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
    catalyst_class: Optional[str],
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
    if isinstance(catalyst_class, str) and catalyst_class.strip():
        cc_val = CATALYST_CLASS_VALUE_MAP.get(catalyst_class, catalyst_class).strip()
        if cc_val:
            relax["catalyst_class"] = cc_val
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

        with gr.Tab("Condition MCP"):
            gr.Markdown(
                f"""
                ### Condition MCP playground
                Exercise the thin MCP tool wrappers and inspect read-only resources.
                Current schema baseline: `{DEFAULT_SCHEMA_VERSION}`.
                """
            )
            with gr.Row():
                with gr.Column():
                    gr.Markdown("#### Tool wrappers")
                    mcp_rxn = gr.Textbox(
                        label="Reaction SMILES",
                        value="Brc1ccccc1.Nc1ccccc1>>",
                        lines=3,
                    )
                    mcp_agents = gr.Checkbox(label="Include agents in normalize output", value=True)
                    mcp_norm_btn = gr.Button("Normalize reaction", variant="primary")
                    mcp_norm_json = gr.JSON(label="normalize_reaction output")
                    mcp_reactants = gr.Textbox(
                        label="Reactant SMILES (auto-filled)",
                        value="",
                        lines=4,
                    )
                    mcp_norm_btn.click(
                        ui_mcp_normalize,
                        inputs=[mcp_rxn, mcp_agents],
                        outputs=[mcp_norm_json, mcp_reactants],
                    )

                    with gr.Row():
                        mcp_detect_btn = gr.Button("Detect family", variant="secondary")
                        mcp_feat_btn = gr.Button("Featurize substrates", variant="secondary")
                    mcp_detect_json = gr.JSON(label="detect_family output")
                    mcp_feat_json = gr.JSON(label="featurize_substrates output")
                    mcp_detect_btn.click(ui_mcp_detect, inputs=[mcp_reactants], outputs=[mcp_detect_json])
                    mcp_feat_btn.click(ui_mcp_featurize, inputs=[mcp_reactants], outputs=[mcp_feat_json])

                with gr.Column():
                    gr.Markdown("#### Resources")
                    mcp_schema_btn = gr.Button("Show ConditionSet schema", variant="secondary")
                    mcp_schema_json = gr.JSON(label="doc://schemas/condition_set.json")
                    mcp_schema_btn.click(ui_mcp_condition_set_schema, outputs=[mcp_schema_json])

                    if CONDITION_RULE_FAMILIES:
                        mcp_rule_family = gr.Dropdown(
                            label="Condition Rule Library family",
                            choices=CONDITION_RULE_FAMILIES,
                            value=CONDITION_RULE_FAMILIES[0],
                            interactive=True,
                        )
                        mcp_rule_btn = gr.Button("Show rule resource", variant="secondary")
                        mcp_rule_json = gr.JSON(label="cond://rules/<family>.json")
                        mcp_rule_btn.click(ui_mcp_rule_resource, inputs=[mcp_rule_family], outputs=[mcp_rule_json])
                        mcp_rule_family.change(ui_mcp_rule_resource, inputs=[mcp_rule_family], outputs=[mcp_rule_json])
                    else:
                        gr.Markdown("⚠️ Condition Rule Library not available; install `condition-rule-library` data bundle.")

        with gr.Tab("Featurize"):
            gr.Markdown(
                """
                ### Role-aware single input featurization
                Provide a molecule SMILES, a reaction SMILES (we'll featurize the reactants), or a dot/newline-separated list of starting materials.
                """
            )
            smi_in = gr.Textbox(label="SMILES or Reaction SMILES", value="Clc1ccccc1", lines=3)
            roles_in = gr.CheckboxGroup(
                choices=["amine", "alcohol", "aryl_halide"],
                label="Roles (optional; applied to each reactant)",
                value=[],
            )
            show_full = gr.Checkbox(value=False, label="Show full fields list (unchecked shows preview)")
            with gr.Row():
                smi_btn = gr.Button("Featurize", variant="primary")
                dl_btn = gr.DownloadButton(label="Download CSV")
            smi_out = gr.JSON(label="Result (vector lengths, masks, fields)")
            rdkit_status = gr.Markdown(label="Status")
            csv_code = gr.Code(label="CSV preview")

            def _rdkit_status_message() -> str:
                if not rdkit_available():
                    return (
                        "⚠️ RDKit isn't available or is disabled. Descriptor blocks default to zeros. "
                        "Set `CHEMTOOLS_DISABLE_RDKIT=0` (or unset it) and install RDKit to enable numeric features."
                    )
                return "✅ RDKit detected; numeric descriptors enabled."

            def _format_feature_result(raw: Dict[str, Any] | None, show_full_fields: bool) -> Dict[str, Any]:
                if not isinstance(raw, dict):
                    return {"error": "Featurizer returned an unexpected payload."}
                result = dict(raw)
                vec = result.get("vector")
                if isinstance(vec, list):
                    vec_list = vec
                elif hasattr(vec, "tolist"):
                    try:
                        vec_list = list(vec.tolist())  # type: ignore[arg-type]
                    except Exception:
                        vec_list = list(vec) if isinstance(vec, (tuple, list)) else ([] if vec is None else [vec])
                else:
                    vec_list = list(vec) if isinstance(vec, (tuple, list)) else ([] if vec is None else [vec])
                result["vector"] = vec_list
                result["length"] = len(vec_list)
                fields = result.get("fields")
                if isinstance(fields, list):
                    field_list = list(fields)
                elif isinstance(fields, tuple):
                    field_list = list(fields)
                else:
                    field_list = []
                if not show_full_fields:
                    result["fields_preview"] = field_list[:12]
                    result.pop("fields", None)
                else:
                    result["fields"] = field_list
                return result

            def _extract_reactants(text: str, treat_collection: bool) -> List[str]:
                if treat_collection:
                    return _split_reactant_block(text)
                head = text.split(">", 1)[0]
                return _split_reactant_block(head)

            def _ui_featurize(input_text: str, roles: List[str], show_full_fields: bool):
                if role_feat_mol is None:
                    return {"error": "role-aware featurizer unavailable"}, ""
                cleaned = (input_text or "").strip()
                if not cleaned:
                    return {"error": "Provide a SMILES string or reaction SMILES."}, ""
                resolved_roles: Optional[List[str]] = roles or None
                is_reaction = ">" in cleaned
                treat_collection = not is_reaction and len(_split_reactant_block(cleaned)) > 1

                try:
                    if is_reaction or treat_collection:
                        reactants = _extract_reactants(cleaned, treat_collection)
                        if not reactants:
                            return {"error": "No reactant SMILES detected in the input."}, ""
                        payloads: List[Dict[str, Any]] = []
                        for idx, smi in enumerate(reactants):
                            raw = role_feat_mol(smi, roles=resolved_roles)
                            formatted = _format_feature_result(raw, show_full_fields)
                            formatted["smiles"] = smi
                            formatted["reactant_index"] = idx
                            payloads.append(formatted)
                        output: Dict[str, Any] = {
                            "input": cleaned,
                            "input_type": "reaction" if is_reaction else "collection",
                            "reactant_count": len(payloads),
                            "reactants": payloads,
                        }
                        if is_reaction and role_feat_rxn and resolved_roles is None:
                            try:
                                raw_rxn = role_feat_rxn(cleaned)
                                if isinstance(raw_rxn, dict):
                                    output["raw_reaction_payload"] = raw_rxn
                            except Exception:
                                pass
                    else:
                        raw = role_feat_mol(cleaned, roles=resolved_roles)
                        formatted = _format_feature_result(raw, show_full_fields)
                        formatted["smiles"] = cleaned
                        output = {
                            "input": cleaned,
                            "input_type": "molecule",
                            **formatted,
                        }
                except Exception as exc:
                    return {"error": str(exc)}, ""

                return output, _rdkit_status_message()

            def _to_csv(input_text: str, roles: List[str]):
                if role_feat_mol is None:
                    return b"error,role-aware featurizer unavailable\n"
                cleaned = (input_text or "").strip()
                if not cleaned:
                    return b"error,Provide a SMILES string first\n"
                resolved_roles: Optional[List[str]] = roles or None
                is_reaction = ">" in cleaned
                treat_collection = not is_reaction and len(_split_reactant_block(cleaned)) > 1
                reactants = _extract_reactants(cleaned, treat_collection) if (is_reaction or treat_collection) else [cleaned]
                if not reactants:
                    return b"error,No reactants detected\n"
                lines = ["target,smiles,field,value"]
                for idx, smi in enumerate(reactants):
                    target = f"reactant_{idx + 1}" if (is_reaction or treat_collection) else "input"
                    try:
                        raw = role_feat_mol(smi, roles=resolved_roles)
                    except Exception as exc:
                        err = f"error,failed to featurize {smi},{exc}\n"
                        return err.encode("utf-8")
                    vec = raw.get("vector")
                    try:
                        vec_list = vec.tolist()  # type: ignore[attr-defined]
                    except Exception:
                        vec_list = list(vec or [])  # type: ignore[arg-type]
                    fields = raw.get("fields") or []
                    n = min(len(fields), len(vec_list))
                    escaped_smi = smi.replace('"', '""')
                    for i in range(n):
                        name = str(fields[i]).replace('"', '""')
                        val = json.dumps(vec_list[i])
                        lines.append(f'"{target}","{escaped_smi}","{name}",{val}')
                csv_text = "\n".join(lines) + "\n"
                return csv_text.encode("utf-8")

            smi_btn.click(_ui_featurize, inputs=[smi_in, roles_in, show_full], outputs=[smi_out, rdkit_status])

            def _update_csv_preview(input_text: str, roles: List[str]) -> str:
                data = _to_csv(input_text, roles)
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
            gr.Markdown("### Recommendation modes\nSelect between data-driven suggestions and rule-based playbooks.")

            default_reaction_type = RECOMMEND_FAMILY_OPTIONS[0] if RECOMMEND_FAMILY_OPTIONS else ""
            default_summary, default_rule_family, default_rule_features, default_rule_job_ctx = ui_recommend_reaction_type_meta(default_reaction_type)

            with gr.Row():
                rec_in = gr.Textbox(
                    label="Reaction SMILES",
                    value="Brc1ccccc1.Nc1ccccc1>>",
                    lines=3,
                )
                rec_type = gr.Dropdown(
                    label="Reaction Type",
                    choices=RECOMMEND_FAMILY_OPTIONS,
                    value=default_reaction_type,
                    interactive=True,
                )

            rec_type_summary = gr.Markdown(default_summary or "- No reaction type selected -", label="Reaction family mapping")

            with gr.Row():
                rec_detect = gr.Markdown(label="Detection Details")
                rec_starting_md = gr.Markdown(label="Starting materials summary")

            with gr.Row():
                feat_btn = gr.Button("Featurize Reactants", variant="secondary")

            with gr.Accordion("Starting materials featurization", open=False):
                rec_starting_table = gr.Dataframe(
                    label="Reactant feature highlights",
                    headers=_STARTING_TABLE_HEADERS,
                    interactive=False,
                )
                rec_starting_json = gr.JSON(label="Raw payload")
            rec_rule_features_json = gr.Textbox(visible=False)

            rule_family_component = None
            rule_features_component = None
            rule_job_component = None
            feat_click = None

            with gr.Tabs():
                with gr.Tab("Structured (ML)"):
                    rec_k = gr.Slider(label="k (neighbors)", minimum=5, maximum=100, value=50, step=1)
                    rec_limit = gr.Slider(label="Number of recommendations", minimum=1, maximum=10, value=5, step=1)
                    rec_cat = gr.Dropdown(
                        label="Catalyst class (optional filter)",
                        choices=CATALYST_CLASS_OPTIONS,
                        value=CATALYST_CLASS_OPTIONS[0],
                        interactive=True,
                    )
                    rec_use_rxi = gr.Checkbox(label="Use rxn-insight auto-detect (if available)", value=True)
                    rec_relax = gr.Textbox(label="Relax (JSON)", value="")
                    rec_constraints = gr.Textbox(label="Constraints (JSON)", value="")
                    rec_btn = gr.Button("Recommend", variant="primary")
                    rec_out = gr.JSON(label="Structured Recommendations")
                    rec_tbl = gr.Dataframe(headers=["rank", "core", "base", "solvent", "confidence", "support"], label="Recommendation Summary", interactive=False)
                    rec_human = gr.Textbox(label="Top recommended (core/base/solvent)", interactive=False)
                    rec_click = rec_btn.click(
                        ui_recommend_structured,
                        inputs=[rec_in, rec_type, rec_k, rec_limit, rec_cat, rec_use_rxi, rec_relax, rec_constraints],
                        outputs=[rec_out, rec_tbl, rec_human, rec_detect, rec_starting_md, rec_starting_table, rec_starting_json],
                    )
                    feat_click = feat_btn.click(
                        ui_featurize_reactants_only,
                        inputs=[rec_in, rec_type, rec_use_rxi, rec_relax],
                        outputs=[rec_detect, rec_starting_md, rec_starting_table, rec_starting_json],
                    )

                with gr.Tab("Rule-based (CRL)"):
                    if CONDITION_RULE_FAMILIES:
                        gr.Markdown(
                            "Use the shared reaction inputs above, then inspect deterministic playbooks for the mapped rule family. "
                            "Edit the JSON payloads to explore guard and prior behavior."
                        )
                        rule_family_component = gr.Textbox(
                            label="Rule family (auto-mapped)",
                            value=default_rule_family,
                            interactive=False,
                        )
                        rule_features_component = gr.Textbox(
                            label="Features (JSON)",
                            value=default_rule_features,
                            lines=12,
                        )
                        rule_job_component = gr.Textbox(
                            label="Job context (JSON, optional)",
                            value=default_rule_job_ctx,
                            lines=6,
                        )
                        rule_topn = gr.Slider(
                            label="Top playbooks",
                            minimum=1,
                            maximum=10,
                            value=3,
                            step=1,
                        )
                        rule_btn = gr.Button("Run rule-based recommender", variant="secondary")
                        rule_json = gr.JSON(label="Rule recommendations")
                        rule_table = gr.Dataframe(label="Playbook summary", interactive=False)
                        rule_md = gr.Markdown(label="Details")
                        rule_btn.click(
                            ui_condition_rule_recommend,
                            inputs=[rule_family_component, rule_features_component, rule_job_component, rule_topn, rec_rule_features_json],
                            outputs=[rule_json, rule_table, rule_md],
                        )
                    else:
                        gr.Markdown("⚠️ Rule-based condition library not available in this build.")

            def _auto_rule_payload(structured: Dict[str, Any]) -> Tuple[Any, Any, Any]:
                if not isinstance(structured, dict):
                    return gr.update(), gr.update(), gr.update()
                features_obj = structured.get("rule_features")
                fallback_family = structured.get("family") if isinstance(structured.get("family"), str) else None
                return _rule_feature_component_updates(features_obj, fallback_family)

            if CONDITION_RULE_FAMILIES and rule_family_component and rule_features_component and rule_job_component:
                rec_type.change(
                    ui_recommend_reaction_type_meta,
                    inputs=[rec_type],
                    outputs=[rec_type_summary, rule_family_component, rule_features_component, rule_job_component],
                )
                if rec_click is not None:
                    rec_click.then(
                        _auto_rule_payload,
                        inputs=[rec_out],
                        outputs=[rule_family_component, rule_features_component, rec_rule_features_json],
                    )
                if feat_click is not None:
                    feat_click.then(
                        _rule_payload_from_starting,
                        inputs=[rec_starting_json],
                        outputs=[rule_family_component, rule_features_component, rec_rule_features_json],
                    )
            else:
                def _meta_summary_only(selection: str) -> str:
                    summary, _, _, _ = ui_recommend_reaction_type_meta(selection)
                    return summary or "- No reaction type selected -"

                rec_type.change(_meta_summary_only, inputs=[rec_type], outputs=[rec_type_summary])


        with gr.Tab("Design Plate"):
            plate_in = gr.Textbox(label="Reaction SMILES", value="Brc1ccccc1.Nc1ccccc1>>")
            plate_n = gr.Slider(label="Plate size", minimum=6, maximum=96, value=24, step=1)
            plate_cat = gr.Dropdown(
                label="Catalyst class (optional filter)",
                choices=CATALYST_CLASS_OPTIONS,
                value=CATALYST_CLASS_OPTIONS[0],
                interactive=True,
            )
            plate_use_rxi = gr.Checkbox(label="Use rxn-insight auto-detect (if available)", value=True)
            plate_relax = gr.Textbox(label="Relax (JSON)", value="")
            plate_constraints = gr.Textbox(label="Constraints (JSON)", value="")
            plate_btn = gr.Button("Design", variant="primary")
            plate_csv = gr.Code(label="CSV")
            plate_tbl = gr.Dataframe(label="Preview", interactive=False)
            plate_meta = gr.JSON(label="Meta")
            plate_btn.click(
                ui_design_plate,
                inputs=[plate_in, plate_n, plate_cat, plate_use_rxi, plate_relax, plate_constraints],
                outputs=[plate_csv, plate_tbl, plate_meta],
            )

        with gr.Tab("Precedent Search"):
            ps_in = gr.Textbox(label="Reaction SMILES", value="Brc1ccccc1.Nc1ccccc1>>")
            ps_k = gr.Slider(label="k (neighbors)", minimum=5, maximum=200, value=50, step=1)
            ps_cat = gr.Dropdown(
                label="Catalyst class (optional filter)",
                choices=CATALYST_CLASS_OPTIONS,
                value=CATALYST_CLASS_OPTIONS[0],
                interactive=True,
            )
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
                    ps_cat,
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

