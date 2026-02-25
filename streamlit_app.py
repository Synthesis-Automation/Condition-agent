"""
ChemCoworker — Streamlit GUI
=============================
Run with:
    streamlit run streamlit_app.py

Requires:
    pip install streamlit
"""
from __future__ import annotations

import json
import os
import re
import tempfile
import uuid
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import streamlit as st

# ---------------------------------------------------------------------------
# Page config — must be first Streamlit call
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="ChemCoworker",
    page_icon="⚗",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------------------------------------------------------------------------
# Local imports (lazy so import errors surface as readable Streamlit messages)
# ---------------------------------------------------------------------------
try:
    from chem_coworker.agent import ChemCoworker
    from chem_coworker.response import ChemResponse
except ImportError as e:
    st.error(f"Could not import chem_coworker: {e}")
    st.stop()

try:
    from chemtools.visualization.rendering import (
        render_molecule_image,
        render_reaction_image,
    )
    _RENDERING_AVAILABLE = True
except Exception:
    _RENDERING_AVAILABLE = False

# ---------------------------------------------------------------------------
# Model catalogue (mirrors chem_coworker/cli.py SELECTABLE_MODELS)
# ---------------------------------------------------------------------------
SELECTABLE_MODELS: List[Dict[str, str]] = [
    {"name": "o4-mini",       "provider": "openai"},
    {"name": "gpt-5.2",       "provider": "openai"},
    {"name": "glm-5",         "provider": "aliyun"},
    {"name": "glm-4.7",       "provider": "aliyun"},
    {"name": "MiniMax-M2.1",  "provider": "aliyun"},
    {"name": "deepseek-v3.2", "provider": "aliyun"},
]
_MODELS_BY_PROVIDER: Dict[str, List[str]] = {}
for _m in SELECTABLE_MODELS:
    _MODELS_BY_PROVIDER.setdefault(_m["provider"], []).append(_m["name"])

# ---------------------------------------------------------------------------
# SMILES regex (reactions contain >>, molecules don't)
# ---------------------------------------------------------------------------
_RXNSMILES_RE = re.compile(
    r"[A-Za-z0-9@+\-\[\]()\\/%=#$!:.]{4,}"   # reactants side
    r">>"
    r"[A-Za-z0-9@+\-\[\]()\\/%=#$!:.]{3,}"   # products side
)
_MOLSMILES_RE = re.compile(
    r"(?<!\w)([A-Za-z0-9@+\-\[\]()\\/%=#$!:.]{6,})(?!\w)"
)


# ===========================================================================
# Helper: SMILES extraction
# ===========================================================================

def _extract_smiles(
    query: str,
    tool_results: Dict[str, Any],
) -> Optional[Tuple[str, str]]:
    """
    Return (smiles, kind) for the primary input SMILES, or None.

    Priority:
      1. Normalised reaction SMILES from normalize_reaction / analyze_bond_changes
      2. Reaction SMILES from the user query
      3. Target molecule SMILES from inspect_target
    """
    # 1. Normalised reaction SMILES from tools
    for key in ("normalize_reaction", "analyze_bond_changes"):
        res = tool_results.get(key, {})
        for field in ("reaction_smiles", "canonical_smiles", "smiles"):
            val = res.get(field)
            if val and ">>" in str(val):
                return str(val), "reaction"

    for key in ("inspect_target",):
        res = tool_results.get(key, {})
        val = res.get("smiles") or res.get("canonical_smiles")
        if val and ">>" not in str(val):
            return str(val), "molecule"

    # 2. Reaction SMILES from user query
    m = _RXNSMILES_RE.search(query)
    if m:
        return m.group(0), "reaction"

    return None


def _extract_reaction_smiles_from_text(text: str) -> List[str]:
    """
    Extract all unique reaction SMILES from answer text.

    The system prompts instruct the model to always wrap generated SMILES in
    backticks, e.g. `A.B>>C`, so backtick-delimited patterns are checked first
    (higher confidence).  Plain reaction SMILES in the text are collected second.

    Returns an ordered list of unique SMILES strings.
    """
    found: List[str] = []
    seen: set = set()

    # 1. Backtick-wrapped reaction SMILES  e.g. `Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(...)cc1`
    for m in re.finditer(r"`([^`\s]{4,}>>[^`\s]{3,})`", text):
        s = m.group(1).strip()
        if s not in seen:
            seen.add(s)
            found.append(s)

    # 2. Plain reaction SMILES not already captured
    for m in _RXNSMILES_RE.finditer(text):
        s = m.group(0).strip()
        if s not in seen:
            seen.add(s)
            found.append(s)

    return found


# ===========================================================================
# Rendering helpers
# ===========================================================================

def _render_smiles_image(smiles: str, kind: str) -> None:
    """Render a SMILES string as an image and display it with st.image()."""
    if not _RENDERING_AVAILABLE:
        st.info(f"RDKit rendering unavailable. SMILES: `{smiles}`")
        return

    tmp_path = os.path.join(
        tempfile.gettempdir(), f"chemcoworker_{uuid.uuid4().hex}.png"
    )
    # Track temp files for cleanup on new conversation
    st.session_state.setdefault("_tmp_files", []).append(tmp_path)

    try:
        if kind == "reaction":
            render_reaction_image(smiles, tmp_path, size=(900, 280))
            st.image(tmp_path, use_container_width=True, caption=smiles)
        else:
            render_molecule_image(smiles, tmp_path, size=(400, 260))
            st.image(tmp_path, caption=smiles)
    except Exception as exc:
        st.warning(f"Could not render SMILES image: {exc}\n\nSMILES: `{smiles}`")


def _render_conditions(conditions: List[Dict[str, Any]]) -> None:
    """Render a conditions list as a styled dataframe."""
    if not conditions:
        return

    PREFERRED_COLS = [
        "rank", "catalyst", "ligand", "base", "solvent",
        "temperature", "yield_pct", "yield_percent", "source",
    ]

    df = pd.DataFrame(conditions)
    # Keep only columns that are present, in preferred order
    ordered = [c for c in PREFERRED_COLS if c in df.columns]
    remainder = [c for c in df.columns if c not in ordered]
    df = df[ordered + remainder]

    # Rename yield_percent → yield_pct for consistency
    df = df.rename(columns={"yield_percent": "yield_pct"}, errors="ignore")

    st.markdown("**Recommended conditions**")

    # Highlight top row
    top = conditions[0]
    cols = st.columns(min(4, len([k for k in ("catalyst", "ligand", "base", "solvent") if k in top])))
    display_keys = [k for k in ("catalyst", "ligand", "base", "solvent") if k in top]
    for col, key in zip(cols, display_keys):
        col.metric(key.capitalize(), top.get(key, "—"))

    st.dataframe(df, use_container_width=True, hide_index=True)


def _render_tool_details(
    tool_results: Dict[str, Any],
    tools_called: List[str],
) -> None:
    """Render tool call payloads inside a collapsed expander."""
    if not tool_results:
        return
    n = len(tools_called)
    with st.expander(f"🔧 {n} tool{'s' if n != 1 else ''} used", expanded=False):
        for tool_name in tools_called:
            result = tool_results.get(tool_name, {})
            st.markdown(f"**`{tool_name}`**")
            try:
                st.code(json.dumps(result, indent=2, default=str), language="json")
            except Exception:
                st.code(str(result))


def _render_metadata(response: "ChemResponse") -> None:
    """Render timing, confidence, and warnings in a compact footer."""
    parts = []
    if response.model:
        parts.append(f"`{response.model}`")
    if response.provider:
        parts.append(response.provider)
    if response.elapsed_s:
        parts.append(f"{response.elapsed_s:.1f} s")
    if response.llm_calls:
        parts.append(f"{response.llm_calls} LLM call{'s' if response.llm_calls != 1 else ''}")
    if response.confidence:
        parts.append(f"confidence {response.confidence:.0%}")

    st.caption("  ·  ".join(parts))

    for w in response.warnings:
        if w.startswith("[critic]"):
            continue  # critic findings already inline in answer text
        st.warning(w)


# ===========================================================================
# Main response renderer
# ===========================================================================

def render_response(response: "ChemResponse") -> None:
    """Render a full ChemResponse inside the current st.chat_message context."""
    # Hypothesis badge
    if response.hypothesis:
        task_label = response.task_type.replace("_", " ").title() if response.task_type else ""
        badge = f"**{task_label}** · {response.hypothesis}" if task_label else response.hypothesis
        if response.confidence:
            badge += f" · {response.confidence:.0%} confidence"
        st.caption(badge)

    # Primary input SMILES image (from tools / user query)
    primary_smiles_info = _extract_smiles(response.query, response.tool_results)
    if primary_smiles_info:
        smiles, kind = primary_smiles_info
        _render_smiles_image(smiles, kind)

    # Main answer text
    if response.answer:
        st.markdown(response.answer)

    # Model-generated reaction SMILES found in the answer text
    if response.answer:
        answer_smiles = _extract_reaction_smiles_from_text(response.answer)
        # Exclude the primary SMILES already shown above to avoid duplication
        primary_smiles = primary_smiles_info[0] if primary_smiles_info else None
        step_smiles = [s for s in answer_smiles if s != primary_smiles]
        if step_smiles:
            label = f"Reaction images ({len(step_smiles)} step{'s' if len(step_smiles) != 1 else ''})"
            with st.expander(label, expanded=True):
                for i, s in enumerate(step_smiles, 1):
                    st.caption(f"Step {i}")
                    _render_smiles_image(s, "reaction")

    # Conditions table
    _render_conditions(response.conditions)

    # Tool call details
    _render_tool_details(response.tool_results, response.tools_called)

    # Footer
    _render_metadata(response)


# ===========================================================================
# Session state initialisation
# ===========================================================================

def _init_session(provider: str, model: str) -> None:
    """Initialise or reinitialise session state."""
    if (
        "coworker" not in st.session_state
        or st.session_state.get("_model") != model
        or st.session_state.get("_provider") != provider
    ):
        st.session_state["coworker"] = ChemCoworker(provider=provider, model=model)
        st.session_state["_model"] = model
        st.session_state["_provider"] = provider

    st.session_state.setdefault("history", [])
    st.session_state.setdefault("messages", [])
    st.session_state.setdefault("_total_llm_calls", 0)
    st.session_state.setdefault("_total_elapsed", 0.0)
    st.session_state.setdefault("_tmp_files", [])


def _new_conversation() -> None:
    """Clear conversation and clean up temp image files."""
    for path in st.session_state.get("_tmp_files", []):
        try:
            os.unlink(path)
        except OSError:
            pass
    st.session_state["history"] = []
    st.session_state["messages"] = []
    st.session_state["_total_llm_calls"] = 0
    st.session_state["_total_elapsed"] = 0.0
    st.session_state["_tmp_files"] = []


# ===========================================================================
# Sidebar
# ===========================================================================

def _sidebar() -> Tuple[str, str]:
    """Render sidebar controls; return (provider, model)."""
    with st.sidebar:
        st.title("⚗ ChemCoworker")
        st.divider()

        provider = st.radio(
            "Provider",
            options=list(_MODELS_BY_PROVIDER.keys()),
            horizontal=True,
            key="sidebar_provider",
        )
        model_options = _MODELS_BY_PROVIDER[provider]
        model = st.selectbox(
            "Model",
            options=model_options,
            key="sidebar_model",
        )

        st.divider()

        if st.button("New conversation", use_container_width=True):
            _new_conversation()
            st.rerun()

        st.divider()
        st.markdown("**Session stats**")
        st.metric("LLM calls", st.session_state.get("_total_llm_calls", 0))
        st.metric("Total time", f"{st.session_state.get('_total_elapsed', 0.0):.1f} s")

        if not _RENDERING_AVAILABLE:
            st.warning(
                "RDKit not available — SMILES images disabled. "
                "Install rdkit to enable rendering."
            )

        st.divider()
        st.caption("Powered by LangChain + RDKit")

    return provider, model


# ===========================================================================
# Main app
# ===========================================================================

def main() -> None:
    provider, model = _sidebar()
    _init_session(provider, model)

    st.header("ChemCoworker", divider="gray")

    # Replay existing messages
    for msg in st.session_state["messages"]:
        with st.chat_message(msg["role"]):
            if msg["role"] == "user":
                st.markdown(msg["text"])
            else:
                render_response(msg["response"])

    # Chat input
    prompt = st.chat_input("Enter a reaction SMILES or chemistry question…")
    if not prompt:
        return

    # Show user message
    with st.chat_message("user"):
        st.markdown(prompt)
    st.session_state["messages"].append({"role": "user", "text": prompt})

    # Run agent
    with st.chat_message("assistant"):
        with st.spinner("Thinking…"):
            try:
                resp, new_history = st.session_state["coworker"].chat(
                    prompt,
                    st.session_state["history"],
                )
            except Exception as exc:
                st.error(f"Agent error: {exc}")
                return

        st.session_state["history"] = new_history
        st.session_state["_total_llm_calls"] += resp.llm_calls
        st.session_state["_total_elapsed"] += resp.elapsed_s

        render_response(resp)

    st.session_state["messages"].append(
        {"role": "assistant", "text": resp.answer, "response": resp}
    )


if __name__ == "__main__":
    main()
