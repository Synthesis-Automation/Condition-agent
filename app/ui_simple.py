# -*- coding: utf-8 -*-
"""
Simplified Gradio UI for Condition Agent - Recommendation Focus

A clean, focused interface for reaction condition recommendations using both
rule-based (SchemeConditionDB) and ML-based approaches.

Run:
  python app/ui_simple.py

Then open http://127.0.0.1:7860
"""

from __future__ import annotations

import json
import os
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Suppress rxnmapper pkg_resources deprecation warning
warnings.filterwarnings("ignore", message="pkg_resources is deprecated", category=UserWarning)
warnings.filterwarnings("ignore", module="rxnmapper")

# Ensure project root is importable
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import gradio as gr

# Import chemtools modules
from chemtools import recommend, router, smiles, featurizers
from chemtools import reagent_lookup  # Reagent database lookup
from chemtools import output_formatter  # Enhanced output formatting
from chemtools import precedent  # For dataset pre-loading

# Optional: Reaction type detector (rxn-insight)
try:
    from chemtools.reaction_type_detector import detect_reaction_type, is_available as detector_available
    DETECTOR_AVAILABLE = detector_available()
except Exception:
    detect_reaction_type = None
    DETECTOR_AVAILABLE = False

# Optional: SchemeConditionDB matcher
try:
    from chemtools.scdb_matcher import load_db as scdb_load_db, match as scdb_match
    SCDB_AVAILABLE = True
except Exception:
    scdb_load_db = None
    scdb_match = None
    SCDB_AVAILABLE = False

# ============================================================================
# Dataset Pre-loading & Caching
# ============================================================================

def preload_datasets():
    """Pre-load all reaction datasets into memory at startup.
    
    This function loads datasets eagerly to avoid the 157-second delay
    on first query. The precedent._load() function uses @lru_cache,
    so this data stays in memory for subsequent queries.
    """
    print("\n" + "="*70)
    print("PRE-LOADING REACTION DATASETS INTO MEMORY")
    print("="*70)
    
    start_time = time.time()
    
    # Trigger dataset loading by calling precedent._load()
    # This caches all datasets in memory via @lru_cache decorator
    try:
        print("Loading all reaction datasets...", end=" ", flush=True)
        
        # Call _load() which loads ALL datasets from data/reaction_dataset/
        # and caches them in memory
        all_precedents = precedent._load()
        total_count = len(all_precedents)
        
        load_time = time.time() - start_time
        print(f"✓ Done!")
        print(f"Loaded {total_count:,} total precedents in {load_time:.1f}s")
        
        # Show breakdown by reaction type
        type_counts = {}
        for prec in all_precedents:
            rxn_type = prec.get('rxn_type', 'Unknown')
            type_counts[rxn_type] = type_counts.get(rxn_type, 0) + 1
        
        print("\nBreakdown by reaction type:")
        for rxn_type, count in sorted(type_counts.items(), key=lambda x: -x[1]):
            print(f"  - {rxn_type}: {count:,} precedents")
        
        print("="*70)
        print(f"✓ Datasets ready! First query will be {load_time:.0f}s faster.\n")
        
        return True
        
    except Exception as e:
        load_time = time.time() - start_time
        print(f"✗ Failed after {load_time:.1f}s")
        print(f"Error: {e}")
        print("="*70)
        print("⚠ Will load datasets on-demand (first query will be slower)\n")
        return False

# ============================================================================
# Configuration
# ============================================================================

# Paths to rule databases (SchemeConditionDB)
SCDB_DIR = Path(ROOT) / "data" / "conditionDB"
RULE_DATABASES = {
    "C-N Coupling (Cu)": str(SCDB_DIR / "C_N_Coupling_Cu_db.json"),
    "C-N Coupling (Pd)": str(SCDB_DIR / "C_N_Coupling_Pd_db.json"),
    "C-N Coupling (Ni)": str(SCDB_DIR / "C_N_Coupling_Ni_db.json"),
    "Amide Formation": str(SCDB_DIR / "Amide_formation_db.json"),
    "Suzuki Coupling": str(SCDB_DIR / "Suzuki_db.json"),
}

# ML family name mappings
ML_FAMILY_MAP = {
    "Auto-detect": None,
    "C-N Coupling (Cu)": "C_N_Coupling_Cu",
    "C-N Coupling (Pd)": "C_N_Coupling_Pd",
    "C-N Coupling (Ni)": "C_N_Coupling_Ni",
    "Amide Formation": "Amide_Coupling",
    "Suzuki Coupling": "Suzuki_CC",
}

# Gradio theme
THEME = gr.themes.Soft(
    primary_hue="blue",
    secondary_hue="slate",
    neutral_hue="slate",
    font=[gr.themes.GoogleFont("Inter"), "ui-sans-serif", "system-ui"],
    font_mono=[gr.themes.GoogleFont("JetBrains Mono"), "ui-monospace"],
    radius_size="sm",
    spacing_size="sm",
)

# ============================================================================
# Helper Functions
# ============================================================================

def format_ml_recommendations(data: Dict[str, Any]) -> Tuple[str, List[List[Any]]]:
    """Format ML-based recommendation results."""
    
    # Extract detection info
    detection = data.get("detection", {})
    family = detection.get("family", "Unknown")
    confidence = detection.get("confidence", 0.0)
    
    # Build detection summary
    lines = [
        f"**Detected Reaction Type:** {family}",
        f"**Confidence:** {confidence:.2%}",
        ""
    ]
    
    # Extract recommendations
    recs = data.get("recommendations", [])
    
    if not recs:
        lines.append("*No recommendations found.*")
        return "\n".join(lines), []
    
    lines.append(f"**Found {len(recs)} recommendation(s):**")
    
    # Build table rows
    table_data = []
    for i, rec in enumerate(recs, 1):
        summary = rec.get("summary", {})
        
        # Extract components from summary
        core = summary.get("core", "")
        base = summary.get("base", "")
        solvent = summary.get("solvent", "")
        confidence = summary.get("confidence", "")
        support = summary.get("support", "")
        
        # Extract catalytic system from chemicals array
        chemicals = rec.get("chemicals", [])
        metal_precursors = [c for c in chemicals if c.get("role") == "metal_precursor"]
        ligands = [c for c in chemicals if c.get("role") == "ligand"]
        
        # Format core (remove /none suffix if present)
        if isinstance(core, str) and core.endswith("/none"):
            core = core[:-5]
        
        # Format component display with reagent lookup
        def fmt_component(comp: Any, comp_type: str = "") -> str:
            if isinstance(comp, dict):
                name = comp.get("name", "")
            elif isinstance(comp, str):
                name = comp
            else:
                name = str(comp) if comp else ""
            
            if not name or not comp_type:
                return name
            
            # Try to enrich with reagent info
            try:
                reagent_info = reagent_lookup.enrich_reagent_info(name, comp_type)
                if reagent_info.get('found'):
                    cas = reagent_info.get('cas')
                    abbr = reagent_info.get('abbreviation')
                    full_name = reagent_info.get('name', name)
                    
                    # Format: Name (abbr) [CAS: xxx]
                    parts = [full_name]
                    if abbr and abbr != full_name:
                        parts.append(f"({abbr})")
                    if cas:
                        parts.append(f"[CAS: {cas}]")
                    return " ".join(parts)
            except Exception:
                pass
            
            return name
        
        # Format catalytic system components
        catalyst_parts = []
        for metal in metal_precursors:
            metal_str = fmt_component(metal, "metal_precursor")
            if metal_str:
                catalyst_parts.append(metal_str)
        for lig in ligands:
            lig_str = fmt_component(lig, "ligand")
            if lig_str:
                catalyst_parts.append(lig_str)
        
        catalyst_str = " + ".join(catalyst_parts) if catalyst_parts else ""
        
        base_str = fmt_component(base, "base")
        solvent_str = fmt_component(solvent, "solvent")
        
        table_data.append([
            i,
            core,
            catalyst_str,
            base_str,
            solvent_str,
            f"{confidence:.2%}" if isinstance(confidence, (int, float)) else str(confidence),
            str(support) if support else ""
        ])
    
    return "\n".join(lines), table_data


def format_rule_recommendations(result: Any, db_name: str) -> Tuple[str, List[List[Any]]]:
    """Format rule-based (SchemeConditionDB) recommendation results."""
    
    if not SCDB_AVAILABLE:
        return "**Error:** SchemeConditionDB matcher not available.", []
    
    # Extract result data - handle both object attributes and dict
    try:
        # Try to get as dict first
        if hasattr(result, 'to_json_dict'):
            result_dict = result.to_json_dict()
        elif isinstance(result, dict):
            result_dict = result
        else:
            result_dict = {}
            
        # Get attributes directly from object if dict is empty
        entry_id = getattr(result, 'entry_id', result_dict.get('entry_id', ''))
        entry_name = getattr(result, 'entry_name', result_dict.get('entry_name', ''))
        match_type = getattr(result, 'match_type', result_dict.get('match_type', ''))
        
        # Get conditions - try multiple approaches
        conditions = None
        if hasattr(result, 'conditions'):
            conditions = result.conditions
        elif 'conditions' in result_dict:
            conditions = result_dict['conditions']
            
    except Exception as e:
        return f"**Error:** Failed to process result: {e}", []
    
    # Build summary
    lines = [
        f"**Database:** {db_name}",
        f"**Matched Rule:** {entry_id}",
        f"**Rule Name:** {entry_name}",
        f"**Match Type:** {match_type}",
        ""
    ]
    
    # Check if conditions exist
    if not conditions:
        lines.append("*No conditions returned.*")
        return "\n".join(lines), []
    
    lines.append("**Recommended Conditions:**")
    
    # Build conditions table
    # The conditions structure is: {"category": [option1, option2, ...], ...}
    # or nested: {"loadings": {"Pd_mol%": [0.5, 1.5], ...}, ...}
    table_data = []
    
    # Map condition keys to reagent types for lookup
    reagent_type_map = {
        'pd_source': 'metal_precursor',
        'cu_source': 'metal_precursor',
        'ni_source': 'metal_precursor',
        'metal_source': 'metal_precursor',
        'catalyst': 'metal_precursor',
        'ligand': 'ligand',
        'ligands': 'ligand',
        'base': 'base',
        'solvent': 'solvent',
        'solvents': 'solvent',
        'oxidant': 'oxidant',
        'reductant': 'reductant',
    }
    
    def enrich_option(option_str: str, category_key: str) -> str:
        """Enrich a single option with reagent database info."""
        reagent_type = reagent_type_map.get(category_key.lower())
        if not reagent_type:
            return option_str
        
        try:
            info = reagent_lookup.enrich_reagent_info(option_str, reagent_type)
            if info.get('found'):
                parts = [info.get('name', option_str)]
                
                abbr = info.get('abbreviation')
                if abbr and abbr != option_str:
                    parts.append(f"({abbr})")
                
                cas = info.get('cas')
                if cas:
                    parts.append(f"[CAS: {cas}]")
                
                return " ".join(parts)
        except Exception:
            pass
        
        return option_str
    
    if isinstance(conditions, dict):
        for key, value in sorted(conditions.items()):
            # Format the key nicely
            display_key = key.replace("_", " ").title()
            
            # Handle nested dict (e.g., loadings)
            if isinstance(value, dict):
                for sub_key, sub_value in sorted(value.items()):
                    sub_display_key = sub_key.replace("_", " ").title()
                    if isinstance(sub_value, list):
                        options = " or ".join(str(v) for v in sub_value)
                        table_data.append([sub_display_key, options, ""])
                    else:
                        table_data.append([sub_display_key, str(sub_value), ""])
            
            # Handle list of options - enrich each option
            elif isinstance(value, list):
                if len(value) == 0:
                    options = "(none specified)"
                elif len(value) == 1:
                    options = enrich_option(str(value[0]), key)
                else:
                    # Enrich each option
                    enriched_options = [enrich_option(str(v), key) for v in value]
                    
                    # For better readability, use newlines if options are long
                    if any(len(s) > 40 for s in enriched_options):
                        options = "\n".join(f"• {s}" for s in enriched_options)
                    else:
                        options = " or ".join(enriched_options)
                
                table_data.append([display_key, options, ""])
            
            # Handle single values - enrich if possible
            else:
                enriched_value = enrich_option(str(value), key)
                table_data.append([display_key, enriched_value, ""])
    
    elif isinstance(conditions, str):
        table_data.append(["Conditions", conditions, ""])
    
    else:
        lines.append(f"*Conditions format not recognized (type: {type(conditions)})*")
    
    if not table_data:
        lines.append("*No condition data to display.*")
    
    return "\n".join(lines), table_data


def auto_detect_reaction_type(reaction_smiles: str) -> str:
    """Auto-detect reaction type from SMILES (for display only)."""
    try:
        result = router.detect_family(reaction_smiles)
        family = result.get("family", "Unknown")
        confidence = result.get("confidence", 0.0)
        return f"**Auto-detected:** {family} (confidence: {confidence:.2%})"
    except Exception as e:
        return f"**Auto-detection failed:** {e}"


def detect_and_map_reaction_type(reaction_smiles: str, requested_type: str) -> Dict[str, Any]:
    """
    Detect reaction type and map to both ML and rule-based naming conventions.
    
    Returns dict with:
        - detected_family: str - The detected chemtools family (e.g., "Buchwald_CN", "Ullmann_CN")
        - ml_family: str - ML API family name (e.g., "C_N_Coupling_Pd", "C_N_Coupling_Cu")
        - rule_db_name: str - Rule database name (e.g., "C-N Coupling (Pd)", "C-N Coupling (Cu)")
        - confidence: float - Detection confidence
        - method: str - Detection method used
        - success: bool - Whether detection succeeded
        - message: str - Human-readable detection message
    """
    
    # If not auto-detect, return the requested type
    if requested_type != "Auto-detect":
        return {
            "detected_family": requested_type,
            "ml_family": ML_FAMILY_MAP.get(requested_type),
            "rule_db_name": requested_type,
            "confidence": 1.0,
            "method": "user_specified",
            "success": True,
            "message": f"Using user-specified type: {requested_type}",
        }
    
    # Try rxn-insight detector first (more accurate)
    if DETECTOR_AVAILABLE and detect_reaction_type:
        try:
            result = detect_reaction_type(reaction_smiles)
            if result.get("success"):
                mapped_family = result.get("mapped_family")  # e.g., "Buchwald_CN", "Ullmann_CN", "Suzuki_CC"
                rxn_class = result.get("rxn_class")  # e.g., "Heteroatom Alkylation and Arylation"
                rxn_name = result.get("rxn_name")  # e.g., "Buchwald-Hartwig C-N coupling"
                confidence = result.get("confidence", 0.8)
                catalysts = result.get("catalysts", [])
                
                # Map to ML and rule database names
                ml_family = None
                rule_db_name = None
                
                # Handle C-N coupling variations based on catalyst
                if mapped_family in ["Buchwald_CN", "C_N_Coupling_Pd"]:
                    ml_family = "C_N_Coupling_Pd"
                    rule_db_name = "C-N Coupling (Pd)"
                    detected_family = "Buchwald_CN (Pd-catalyzed)"
                elif mapped_family == "Ullmann_CN" or mapped_family == "C_N_Coupling_Cu":
                    ml_family = "C_N_Coupling_Cu"
                    rule_db_name = "C-N Coupling (Cu)"
                    detected_family = "Ullmann_CN (Cu-catalyzed)"
                elif mapped_family == "C_N_Coupling_Ni":
                    ml_family = "C_N_Coupling_Ni"
                    rule_db_name = "C-N Coupling (Ni)"
                    detected_family = "C_N_Coupling_Ni"
                elif mapped_family in ["Suzuki_CC", "Suzuki"]:
                    ml_family = "Suzuki_CC"
                    rule_db_name = "Suzuki Coupling"
                    detected_family = "Suzuki_CC"
                elif mapped_family == "Amide_Coupling":
                    ml_family = "Amide_Coupling"
                    rule_db_name = "Amide Formation"
                    detected_family = "Amide_Coupling"
                else:
                    detected_family = mapped_family or "Unknown"
                
                # Build message
                msg_parts = []
                if rxn_class:
                    msg_parts.append(f"Class: {rxn_class}")
                if rxn_name:
                    msg_parts.append(f"Name: {rxn_name}")
                if catalysts:
                    msg_parts.append(f"Catalysts: {', '.join(catalysts)}")
                
                message = f"**Auto-detected (rxn-insight):** {detected_family}\n"
                if msg_parts:
                    message += "  " + " | ".join(msg_parts)
                if confidence:
                    message += f"\n  Confidence: {confidence:.2%}"
                
                return {
                    "detected_family": detected_family,
                    "ml_family": ml_family,
                    "rule_db_name": rule_db_name,
                    "confidence": confidence or 0.8,
                    "method": "rxn_insight",
                    "success": True,
                    "message": message,
                    "raw_detection": result,
                }
        except Exception as e:
            print(f"⚠ rxn-insight detection failed: {e}")
    
    # Fallback to router-based detection
    try:
        result = router.detect_family(reaction_smiles)
        detected_family = result.get("family", "Unknown")
        confidence = result.get("confidence", 0.0)
        
        # Map router family names to ML and rule DB names
        family_to_ml = {
            "C_N_Coupling_Cu": "C_N_Coupling_Cu",
            "C_N_Coupling_Pd": "C_N_Coupling_Pd",
            "C_N_Coupling_Ni": "C_N_Coupling_Ni",
            "Ullmann_CN": "C_N_Coupling_Cu",
            "Buchwald_CN": "C_N_Coupling_Pd",
            "Amide_Coupling": "Amide_Coupling",
            "Suzuki_CC": "Suzuki_CC",
            "Suzuki": "Suzuki_CC",
        }
        
        family_to_db = {
            "C_N_Coupling_Cu": "C-N Coupling (Cu)",
            "C_N_Coupling_Pd": "C-N Coupling (Pd)",
            "C_N_Coupling_Ni": "C-N Coupling (Ni)",
            "Ullmann_CN": "C-N Coupling (Cu)",
            "Buchwald_CN": "C-N Coupling (Pd)",
            "Amide_Coupling": "Amide Formation",
            "Suzuki_CC": "Suzuki Coupling",
            "Suzuki": "Suzuki Coupling",
        }
        
        ml_family = family_to_ml.get(detected_family)
        rule_db_name = family_to_db.get(detected_family)
        
        message = f"**Auto-detected (rule-based):** {detected_family}"
        if confidence > 0:
            message += f"\n  Confidence: {confidence:.2%}"
        
        return {
            "detected_family": detected_family,
            "ml_family": ml_family,
            "rule_db_name": rule_db_name,
            "confidence": confidence,
            "method": "router",
            "success": bool(ml_family or rule_db_name),
            "message": message,
            "raw_detection": result,
        }
        
    except Exception as e:
        return {
            "detected_family": "Unknown",
            "ml_family": None,
            "rule_db_name": None,
            "confidence": 0.0,
            "method": "failed",
            "success": False,
            "message": f"**Auto-detection failed:** {e}",
        }


def warmup_ml_system():
    """Pre-load datasets and caches to avoid slow first recommendation."""
    print("\n" + "="*60)
    print("WARMING UP ML RECOMMENDATION SYSTEM")
    print("="*60)
    
    start_time = time.time()
    
    # Simple dummy reaction to trigger dataset loading
    dummy_reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    
    try:
        print("Loading datasets and initializing caches...")
        _ = recommend.recommend_conditions_structured(
            reaction=dummy_reaction,
            reaction_type="C_N_Coupling_Cu",  # Most common type
            k=5,  # Small k for faster warmup
            limit=1,
            relax={
                "use_drfp": False,
                "precompute_drfp": False,
                "use_rxn_insight": True,
            },
            constraints=None,
        )
        
        elapsed = time.time() - start_time
        print(f"✓ Warmup complete in {elapsed:.1f}s")
        print(f"  Subsequent recommendations will be much faster!")
        
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"⚠ Warmup failed after {elapsed:.1f}s: {e}")
        print(f"  (This is OK, system will warm up on first use)")
    
    print("="*60 + "\n")


# ============================================================================
# Core Recommendation Functions
# ============================================================================

def _get_precedents_for_reaction(
    reaction_smiles: str,
    family: str,
    k: int = 10,
    relax_settings: Optional[Dict[str, Any]] = None
) -> Tuple[str, str]:
    """Get precedents for a reaction and format them with reaction schemes."""
    from rdkit import Chem
    from rdkit.Chem import Draw
    from PIL import Image as PILImage, ImageDraw as PILImageDraw
    import base64
    from io import BytesIO
    
    # Helper functions for rendering
    def _grid(smiles_str: str) -> Optional[PILImage.Image]:
        if not smiles_str or not smiles_str.strip():
            return None
        parts = [s.strip() for s in smiles_str.split(".") if s.strip()]
        mols = []
        for smi in parts:
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    mols.append(mol)
            except Exception:
                pass
        if not mols:
            return None
        try:
            return Draw.MolsToGridImage(
                mols,
                molsPerRow=len(mols),
                subImgSize=(220, 220),
                returnPNG=False,
            )
        except Exception:
            return None
    
    def _img_data_uri(img: PILImage.Image) -> str:
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        return f"data:image/png;base64,{img_str}"
    
    def _render_img(reactants_smi: str, products_smi: str) -> str:
        try:
            left_img = _grid(reactants_smi)
            right_img = _grid(products_smi)
            
            if left_img is None and right_img is None:
                return ""
            if left_img is None:
                uri = _img_data_uri(right_img)
                return f"<img src='{uri}' width='340'/>" if uri else ""
            if right_img is None:
                uri = _img_data_uri(left_img)
                return f"<img src='{uri}' width='340'/>" if uri else ""
            
            # Combine images with arrow
            try:
                arrow_width = 60
                w = left_img.width + right_img.width + arrow_width
                h = max(left_img.height, right_img.height)
                canvas = PILImage.new("RGB", (w, h), (255, 255, 255))
                
                canvas.paste(left_img, (0, (h - left_img.height) // 2))
                canvas.paste(right_img, (left_img.width + arrow_width, (h - right_img.height) // 2))
                
                dr = PILImageDraw.Draw(canvas)
                y = h // 2
                arrow_start_x = left_img.width + 10
                arrow_end_x = left_img.width + arrow_width - 10
                
                dr.line((arrow_start_x, y, arrow_end_x, y), fill=(0, 0, 0), width=3)
                
                arrow_head_size = 6
                dr.polygon([
                    (arrow_end_x, y),
                    (arrow_end_x - arrow_head_size, y - arrow_head_size),
                    (arrow_end_x - arrow_head_size, y + arrow_head_size)
                ], fill=(0, 0, 0))
                
                uri = _img_data_uri(canvas)
                return f"<img src='{uri}' width='340'/>" if uri else ""
            except Exception:
                uri = _img_data_uri(left_img)
                return f"<img src='{uri}' width='340'/>" if uri else ""
        except Exception:
            return ""
    
    # Prepare search settings
    relax = relax_settings or {
        "use_drfp": True,
        "precompute_drfp": True,
        "drfp_weight": 0.7,
    }
    
    # Only featurize if DRFP is disabled (otherwise we only need reaction SMILES)
    feat = {}
    if not relax.get("use_drfp", False):
        try:
            rx_norm = smiles.normalize_reaction(reaction_smiles)
            
            # Extract reactants from normalized reaction
            reactants = [
                (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
                for r in (rx_norm.get("reactants") or [])
            ]
            
            if not reactants:
                return "No reactants found in reaction", ""
            
            # Pick electrophile and nucleophile using simple heuristic
            def is_electrophile(s: str) -> bool:
                t = (s or "").lower()
                return ("br" in t) or ("cl" in t) or (" i" in t) or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
            
            elec_smi, nuc_smi = "", ""
            if len(reactants) == 1:
                elec_smi, nuc_smi = reactants[0], ""
            else:
                r0, r1 = reactants[0], reactants[1]
                elec_smi, nuc_smi = (r0, r1) if is_electrophile(r0) else ((r1, r0) if is_electrophile(r1) else (r0, r1))
            
            feat = featurizers.molecular.featurize(elec_smi, nuc_smi)
        except Exception as e:
            return f"Could not featurize reaction: {e}", ""
    
    # Add reaction SMILES for DRFP similarity
    relax["reaction_smiles"] = reaction_smiles
    
    try:
        pack = precedent.knn(family, feat, k=k, relax=relax)
    except Exception as e:
        return f"Precedent search failed: {e}", ""
    
    precs = pack.get("precedents", [])
    
    if not precs:
        return "No precedents found", ""
    
    # Build summary
    summary = f"**Precedents Used by ML System**\n\n"
    summary += f"- **Total precedents:** {len(precs)}\n"
    summary += f"- **Reaction family:** {family}\n"
    summary += f"- **Similarity method:** DRFP (70% transformation + 30% substrate)\n"
    
    # Get similarity score range
    scores = [p.get("similarity", 0.0) for p in precs]
    if any(s > 0 for s in scores):
        summary += f"- **Similarity range:** {min(scores):.3f} - {max(scores):.3f}\n"
    
    # Build HTML table
    html_rows = []
    html_rows.append('<table style="border-collapse:collapse; width:100%; table-layout:fixed">')
    html_rows.append(
        "<colgroup>"
        "<col style='width:340px'/>"
        "<col style='width:70px'/>"
        "<col style='width:100px'/>"
        "<col style='width:65px'/>"
        "<col style='width:150px'/>"
        "<col style='width:115px'/>"
        "<col style='width:115px'/>"
        "<col style='width:65px'/>"
        "<col style='width:65px'/>"
        "</colgroup>"
    )
    html_rows.append(
        "<tr style='background-color:#f0f0f0; font-weight:bold'>"
        "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Reaction</th>"
        "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Similarity</th>"
        "<th style='padding:8px; text-align:left; border:1px solid #ddd'>ID</th>"
        "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Yield</th>"
        "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Catalyst</th>"
        "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Base</th>"
        "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Solvent</th>"
        "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Temp °C</th>"
        "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Time h</th>"
        "</tr>"
    )
    
    for i, prec in enumerate(precs, 1):
        # Extract SMILES
        reaction_smi = prec.get("reaction_smiles", "")
        reactants_smi = prec.get("reactants_smiles", "")
        products_smi = prec.get("products_smiles", "")
        
        if not reactants_smi or not products_smi:
            if ">>" in reaction_smi:
                parts = reaction_smi.split(">>")
                reactants_smi = parts[0] if len(parts) > 0 else ""
                products_smi = parts[1] if len(parts) > 1 else ""
        
        img_html = _render_img(reactants_smi, products_smi)
        
        sim_score = prec.get("similarity", 0.0)
        sim_display = f"{sim_score:.3f}" if sim_score > 0 else "-"
        sim_color = "#2e7d32" if sim_score >= 0.8 else "#f57c00" if sim_score >= 0.6 else "#666"
        
        reaction_id = prec.get("reaction_id", "")
        yield_val = prec.get("yield", "")
        core = prec.get("core", "")
        base_uid = prec.get("base_uid", "")
        solvent_uid = prec.get("solvent_uid", "")
        temp = prec.get("T_C", "")
        time_h = prec.get("time_h", "")
        
        try:
            base_label = reagent_lookup.enrich_reagent_info(base_uid, "base").get("name", base_uid) if base_uid else ""
        except Exception:
            base_label = base_uid
        
        try:
            solvent_label = reagent_lookup.enrich_reagent_info(solvent_uid, "solvent").get("name", solvent_uid) if solvent_uid else ""
        except Exception:
            solvent_label = solvent_uid
        
        row_style = "background-color:#fafafa" if i % 2 == 0 else ""
        html_rows.append(
            f"<tr style='{row_style}'>"
            f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top'>{img_html}</td>"
            f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:nowrap; text-align:center; font-weight:bold; color:{sim_color}'>{sim_display}</td>"
            f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{reaction_id[:35]}{'...' if len(reaction_id) > 35 else ''}</td>"
            f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:nowrap; text-align:center'>{yield_val}</td>"
            f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{core}</td>"
            f"<td title='{base_uid}' style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{base_label}</td>"
            f"<td title='{solvent_uid}' style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{solvent_label}</td>"
            f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:nowrap; text-align:center'>{temp}</td>"
            f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:nowrap; text-align:center'>{time_h}</td>"
            "</tr>"
        )
    
    html_rows.append("</table>")
    html = "\n".join(html_rows)
    
    return summary, html


def get_ml_recommendations(
    reaction_smiles: str,
    reaction_type: str,
    top_k: int = 3,  # Default to 3 recommendations
) -> Tuple[str, List[List[Any]], str, str]:
    """Get ML-based recommendations with enhanced formatting and precedents used."""
    
    if not reaction_smiles or not reaction_smiles.strip():
        return "**Error:** Please enter a reaction SMILES.", [], "", ""
    
    # Auto-detect reaction type if needed
    detection_result = detect_and_map_reaction_type(reaction_smiles.strip(), reaction_type)
    
    if not detection_result["success"]:
        return detection_result["message"], [], "", ""
    
    # Use detected ML family name
    family = detection_result["ml_family"]
    detection_message = detection_result["message"]
    detection_confidence = detection_result["confidence"]
    
    # If no ML family mapping available, show helpful error with alternatives
    if not family:
        msg = "**Cannot Proceed with ML Recommendations**\n\n"
        msg += "🔍 **Detection Result:**\n"
        msg += detection_message + "\n\n"
        msg += f"❌ **No ML model available** for: `{detection_result['detected_family']}`\n\n"
        msg += "**Available ML reaction types:**\n"
        for name, ml_name in ML_FAMILY_MAP.items():
            if ml_name:
                msg += f"  - {name} (`{ml_name}`)\n"
        msg += "\n**What to do:**\n"
        msg += "1. ✅ **Try rule-based recommendations** instead\n"
        msg += "2. 🔄 **Manually select** a supported reaction type from the dropdown\n"
        msg += "3. 📝 **Verify your SMILES** represents a supported reaction\n"
        msg += "4. 📖 **Check if your reaction** is similar to a supported type\n"
        return msg, [], "", ""
    
    # Warn if detection confidence is low
    if detection_confidence < 0.5 and reaction_type == "Auto-detect":
        print(f"⚠️ WARNING: Low detection confidence ({detection_confidence:.1%})")
        print(f"   Consider manually selecting the reaction type for better results.")
    
    try:
        print(f"\n{'='*60}")
        print(f"ML RECOMMENDATION TIMING BREAKDOWN")
        print(f"{'='*60}")
        print(f"Detection: {detection_result['method']}")
        print(f"Detected family: {detection_result['detected_family']}")
        print(f"ML family: {family}")
        print(f"Confidence: {detection_result['confidence']:.2%}")
        
        overall_start = time.time()
        
        # Prepare relax settings - use DRFP for universal similarity (same as precedent search)
        relax_settings = {
            "use_drfp": True,  # Enable DRFP for universal reaction similarity
            "precompute_drfp": True,  # Use precomputed fingerprints for speed
            "drfp_weight": 0.7,  # Balance reaction similarity (70%) with substrate features (30%)
            "use_rxn_insight": True,  # Keep fast detection enabled
        }
        
        print(f"Relax settings: {relax_settings}")
        print(f"Family: {family}")
        print(f"k: 10, limit: {top_k}")
        
        # Call recommend_conditions_structured
        print(f"Starting recommendation call...")
        rec_start = time.time()
        
        data = recommend.recommend_conditions_structured(
            reaction=reaction_smiles.strip(),
            reaction_type=family,
            k=10,
            limit=top_k * 2,  # Get more to filter to top 3
            relax=relax_settings,
            constraints=None,
        )
        
        rec_time = time.time() - rec_start
        print(f"Recommendation call completed in {rec_time:.1f}s")
        
        # Debug: Check what we got back
        print(f"DEBUG - Data keys: {list(data.keys())}")
        print(f"DEBUG - Detection: {data.get('detection', {})}")
        print(f"DEBUG - Recommendations count: {len(data.get('recommendations', []))}")
        if data.get('recommendations'):
            print(f"DEBUG - First recommendation keys: {list(data['recommendations'][0].keys())}")
        
        # Format results using new enhanced formatter
        format_start = time.time()
        
        # Build enhanced recommendations
        recommendations_data = []
        recommendations_list = data.get("recommendations", [])
        
        # Check if we have any recommendations
        if not recommendations_list:
            print(f"⚠ No recommendations returned from ML system")
            print(f"Data keys: {list(data.keys())}")
            print(f"Detection: {data.get('detection', {})}")
            
            # Build comprehensive error message with helpful guidance
            detection = data.get("detection", {})
            detected_type = detection.get("reaction_type") or detection_result.get("detected_family", "Unknown")
            
            error_msg = f"**No ML Recommendations Found**\n\n"
            
            # Show detection info
            if reaction_type == "Auto-detect":
                error_msg += f"🔍 **Auto-Detection Result:**\n"
                error_msg += detection_message + "\n\n"
            
            error_msg += f"**Reaction Type:** {detected_type}\n"
            error_msg += f"**ML Family:** {family}\n"
            if detection_result.get("confidence", 0) > 0:
                error_msg += f"**Detection Confidence:** {detection_result['confidence']:.1%}\n"
            error_msg += "\n"
            
            # Provide context-specific guidance
            error_msg += "**Why this happened:**\n\n"
            
            # Check if it's a supported family
            supported_families = [v for v in ML_FAMILY_MAP.values() if v]
            if family not in supported_families:
                error_msg += f"❌ **Unsupported reaction type**: `{family}` is not in the ML training dataset.\n\n"
                error_msg += "**Available ML families:**\n"
                for name, ml_name in ML_FAMILY_MAP.items():
                    if ml_name:
                        error_msg += f"  - {name} (`{ml_name}`)\n"
                error_msg += "\n**What to do:**\n"
                error_msg += "1. ✅ **Try rule-based recommendations** instead (click 'Get Rule Recommendations')\n"
                error_msg += "2. 🔄 **Manually select** a different reaction type from the dropdown\n"
                error_msg += "3. 📝 **Verify your SMILES** format is correct\n"
            else:
                # Family is supported but no precedents found
                error_msg += f"✅ Reaction type `{family}` is supported, but no precedents were found.\n\n"
                error_msg += "**Possible reasons:**\n"
                error_msg += "1. 📊 **No similar reactions** in the precedent database\n"
                error_msg += "2. 🔬 **Unusual substrates** or functional groups\n"
                error_msg += "3. 💾 **Dataset not loaded** (first-run issue)\n"
                error_msg += "4. 🎯 **Detection error** - Wrong family detected\n\n"
                
                error_msg += "**What to do:**\n"
                error_msg += "1. ✅ **Try rule-based recommendations** (click 'Get Rule Recommendations')\n"
                error_msg += "2. 🔄 **Manually select** the correct reaction type if auto-detect was wrong\n"
                error_msg += "3. 🧪 **Simplify substrates** - Try a simpler model reaction first\n"
                error_msg += "4. 📖 **Check literature** for similar transformations\n"
                
                # Check precedent info
                precedents = data.get("precedents", {})
                if isinstance(precedents, dict):
                    total_prec = precedents.get("total_considered", 0)
                    if total_prec > 0:
                        error_msg += f"\n💡 **Note:** {total_prec} precedents were found but none matched your substrates closely enough.\n"
                        error_msg += "   Try adjusting the reaction SMILES or selecting a different reaction type.\n"
            
            return error_msg, [], "", ""
        
        # Limit to top_k
        recommendations_list = recommendations_list[:top_k]
        
        for i, rec in enumerate(recommendations_list, 1):
            # Extract reagents
            reagents = []
            
            # Add starting materials from chemicals array
            for chem in rec.get("chemicals", []):
                if chem.get("role") == "starting_material":
                    reagents.append({
                        "id": f"SM{len([r for r in reagents if r.get('role') in ['electrophile', 'nucleophile']]) + 1}",
                        "name": chem.get("name"),
                        "abbreviation": None,
                        "cas": chem.get("cas"),
                        "smiles": chem.get("smiles"),
                        "inchi_key": None,
                        "role": "electrophile" if not reagents else "nucleophile",
                        "equivalents": {
                            "value": 1.0 if not reagents else 1.2,
                            "range": [1.0, 1.0] if not reagents else [1.1, 1.5],
                            "unit": "eq"
                        }
                    })
            
            # Add catalytic system components (metal_precursor, ligand) from chemicals array
            catalyst_count = 0
            for chem in rec.get("chemicals", []):
                role = chem.get("role")
                if role in ["metal_precursor", "ligand"]:
                    catalyst_count += 1
                    catalyst_name = chem.get("name") or chem.get("cas") or f"Catalyst_{catalyst_count}"
                    # Enrich catalyst using reagent lookup
                    enriched = output_formatter.enrich_reagent(
                        name=catalyst_name,
                        reagent_type="catalyst",
                        role=role,  # Preserve the specific catalytic role
                        equivalents=0.05 if role == "metal_precursor" else 0.1,
                        equiv_range=(0.01, 0.1) if role == "metal_precursor" else (0.05, 0.15),
                        reagent_id=f"CAT{catalyst_count}"
                    )
                    # Override with original data if enrichment didn't find it
                    if not enriched.get("cas") or enriched.get("cas") == catalyst_name:
                        enriched["name"] = chem.get("name") or catalyst_name
                        enriched["cas"] = chem.get("cas")
                        enriched["smiles"] = chem.get("smiles")
                    reagents.append(enriched)
            
            # Get summary data
            summary = rec.get("summary", {})
            
            # Enrich base from summary
            base_obj = summary.get("base", {})
            if base_obj and (base_obj.get("cas") or base_obj.get("name")):
                base_identifier = base_obj.get("cas") or base_obj.get("name")
                reagents.append(output_formatter.enrich_reagent(
                    name=base_identifier,
                    reagent_type="base",
                    role="base",
                    equivalents=2.0,
                    equiv_range=(1.5, 2.5),
                    reagent_id="BASE1"
                ))
            
            # Enrich solvent from summary
            solvent_obj = summary.get("solvent", {})
            if solvent_obj and (solvent_obj.get("cas") or solvent_obj.get("name")):
                solvent_identifier = solvent_obj.get("cas") or solvent_obj.get("name")
                reagents.append(output_formatter.enrich_reagent(
                    name=solvent_identifier,
                    reagent_type="solvent",
                    role="solvent",
                    reagent_id="SOLV1"
                ))
            
            # Extract DRFP similarity score if available (from precedent data)
            # When DRFP is enabled, the precedent should include a similarity score
            drfp_similarity = summary.get("drfp_similarity") or summary.get("similarity", 0.0)
            
            # Extract ML confidence from summary (fallback to DRFP similarity if not available)
            ml_confidence = summary.get("confidence", drfp_similarity or 0.5)
            
            # Use DRFP similarity as the primary similarity metric if available
            similarity_score = drfp_similarity if drfp_similarity > 0 else ml_confidence
            
            # Get support data
            support_data = summary.get("support", {})
            support_count = support_data.get("count", 1) if isinstance(support_data, dict) else support_data
            
            # Extract conditions from the conditions object
            conditions_obj = rec.get("conditions", {})
            temp = conditions_obj.get("temperature")
            time_val = conditions_obj.get("time")
            
            # Use defaults if not provided
            if temp is None:
                temp = 110
            if time_val is None:
                time_val = 12
            
            formatted_cond = output_formatter.format_conditions(
                temperature=temp,
                temp_range=(temp - 20, temp + 20),
                time_hours=time_val,
                time_range=(time_val / 2, time_val * 1.5),
                atmosphere="N₂",
            )
            
            # Build recommendation
            formatted_rec = output_formatter.format_recommendation(
                rank=i,
                confidence=ml_confidence,
                support=support_count,
                reaction_smiles=reaction_smiles.strip(),
                reagents=reagents,
                conditions=formatted_cond,
                similarity_score=similarity_score,  # Use DRFP similarity if available, else ML confidence
                expected_yield=75.0,
                yield_range=(60.0, 85.0),
            )
            
            recommendations_data.append(formatted_rec)
        
        # Build full output
        detection = data.get("detection", {})
        
        # Use our enhanced detection info if available
        detected_type_display = detection_result.get("detected_family", detection.get("family", "Unknown"))
        detection_confidence_display = detection_result.get("confidence", detection.get("confidence", 0.0))
        
        output = output_formatter.format_ml_output(
            reaction_smiles=reaction_smiles.strip(),
            requested_type=reaction_type if reaction_type != "Auto-detect" else None,
            detected_type=detected_type_display,
            detection_confidence=detection_confidence_display,
            recommendations_data=recommendations_data,
            processing_time_ms=(time.time() - overall_start) * 1000,
        )
        
        format_time = time.time() - format_start
        print(f"Formatting completed in {format_time:.1f}s")
        
        elapsed_time = time.time() - overall_start
        print(f"Total time: {elapsed_time:.1f}s")
        print(f"{'='*60}\n")
        
        # Format as JSON for display
        json_output = json.dumps(output, indent=2)
        
        # Create summary with detection info
        summary = f"**ML Recommendations Generated (DRFP-Based)**\n\n"
        
        # Show detection info if auto-detected
        if reaction_type == "Auto-detect":
            summary += detection_message + "\n\n"
        
        summary += f"**Detected Type:** {output['detection']['detected_reaction_type']}\n"
        summary += f"**Confidence:** {output['detection']['confidence']:.2%}\n"
        summary += f"**Recommendations:** {len(output['recommendations'])}\n"
        summary += f"**Similarity Method:** DRFP (Differential Reaction Fingerprints)\n"
        summary += f"**Processing time:** {elapsed_time:.1f} seconds\n\n"
        summary += "*Both ML recommendations and precedent search now use universal DRFP similarity for consistent, reaction-type-agnostic matching.*\n\n"
        summary += "**JSON Output:**\n```json\n" + json_output + "\n```"
        
        # Create table view
        table = []
        for rec in output['recommendations']:
            reagent_list = [f"{r.get('name', r.get('smiles', 'Unknown'))} ({r['role']})" 
                          for r in rec['reagents'] if r.get('role') not in ['electrophile', 'nucleophile']]
            temp = rec['conditions'].get('temperature', {}).get('value', 'N/A')
            table.append([
                rec['rank'],
                f"{rec['confidence']:.2%}",
                ", ".join(reagent_list),
                f"{temp}°C" if temp != 'N/A' else 'N/A',
            ])
        
        # Extract precedents used by the ML system
        # Call precedent search directly to get full precedent details with reaction schemes
        try:
            prec_summary, prec_html = _get_precedents_for_reaction(
                reaction_smiles.strip(),
                family,
                k=10,
                relax_settings=relax_settings
            )
        except Exception as e:
            print(f"Warning: Could not extract precedents: {e}")
            prec_summary = "Precedent information unavailable"
            prec_html = ""
        
        return summary, table, prec_summary, prec_html
        
    except Exception as e:
        import traceback
        
        # Build helpful error message
        error_msg = f"**ML Recommendation System Error**\n\n"
        error_msg += f"❌ **Error Type:** `{type(e).__name__}`\n"
        error_msg += f"**Message:** {str(e)}\n\n"
        
        # Add context
        if reaction_type == "Auto-detect":
            error_msg += "🔍 **Detection Info:**\n"
            error_msg += detection_message + "\n\n"
        
        error_msg += f"**Reaction Type:** {family}\n"
        error_msg += f"**SMILES:** {reaction_smiles[:100]}{'...' if len(reaction_smiles) > 100 else ''}\n\n"
        
        # Common issues and solutions
        error_msg += "**Common Issues & Solutions:**\n\n"
        
        if "dataset" in str(e).lower() or "file" in str(e).lower():
            error_msg += "1. 💾 **Dataset Loading Issue**\n"
            error_msg += "   - The ML precedent database may not be loaded\n"
            error_msg += "   - Try running a simple test reaction first to warm up the system\n"
            error_msg += "   - Check that dataset files exist in `data/` directory\n\n"
        
        if "featurize" in str(e).lower() or "smiles" in str(e).lower():
            error_msg += "2. 📝 **SMILES Parsing Error**\n"
            error_msg += "   - Verify your SMILES syntax is correct\n"
            error_msg += "   - Use standard reaction SMILES format: `reactants>>products`\n"
            error_msg += "   - Check for unusual characters or formatting\n\n"
        
        if "memory" in str(e).lower():
            error_msg += "3. 💻 **Memory Issue**\n"
            error_msg += "   - ML system may be running out of memory\n"
            error_msg += "   - Try reducing `top_k` parameter\n"
            error_msg += "   - Restart the application\n\n"
        
        error_msg += "**What to try:**\n"
        error_msg += "1. ✅ **Use rule-based recommendations** instead (click 'Get Rule Recommendations')\n"
        error_msg += "2. 🔄 **Manually select** the reaction type instead of auto-detect\n"
        error_msg += "3. 📝 **Simplify your SMILES** or try a simpler model reaction\n"
        error_msg += "4. 🔁 **Restart the application** if this persists\n\n"
        
        error_msg += "**Technical Details:**\n```\n" + traceback.format_exc() + "\n```"
        
        return error_msg, [], "", ""


def get_rule_recommendations(
    reaction_smiles: str,
    reaction_type: str,
) -> Tuple[str, List[List[Any]]]:
    """Get rule-based (SchemeConditionDB) recommendations with enhanced formatting."""
    
    if not reaction_smiles or not reaction_smiles.strip():
        return "**Error:** Please enter a reaction SMILES.", []
    
    if not SCDB_AVAILABLE:
        return "**Error:** SchemeConditionDB matcher not installed.", []
    
    overall_start = time.time()
    
    # Auto-detect reaction type if needed
    detection_result = detect_and_map_reaction_type(reaction_smiles.strip(), reaction_type)
    
    if not detection_result["success"]:
        return detection_result["message"], []
    
    # Use detected rule database name
    rule_db_name = detection_result["rule_db_name"]
    detection_message = detection_result["message"]
    
    # If no rule database available, show helpful error
    if not rule_db_name:
        msg = detection_result["message"]
        msg += f"\n\n**No rule database available** for detected family: {detection_result['detected_family']}"
        msg += "\n\nAvailable rule databases:\n"
        for name in RULE_DATABASES.keys():
            msg += f"  - {name}\n"
        return msg, []
    
    # Get database path
    db_path = RULE_DATABASES.get(rule_db_name)
    
    if not db_path:
        return f"**Error:** No rule database configured for '{rule_db_name}'.", []
    
    if not Path(db_path).exists():
        return f"**Error:** Database file not found: {db_path}", []
    
    try:
        # Load database
        db = scdb_load_db(db_path)
        
        # Match reaction
        result = scdb_match(db, reaction_smiles.strip())
        
        # Extract conditions from result
        conditions_dict = {}
        if hasattr(result, 'conditions'):
            conditions_dict = result.conditions
        elif hasattr(result, 'to_json_dict'):
            result_dict = result.to_json_dict()
            if 'conditions' in result_dict:
                conditions_dict = result_dict['conditions']
        
        if not conditions_dict:
            return "**Error:** No conditions found in match result.", []
        
        # Use formatter to expand conditions into 3 recommendations
        recommendations_data = output_formatter.expand_rule_conditions_to_recommendations(
            reaction_smiles=reaction_smiles.strip(),
            conditions_dict=conditions_dict,
            num_recommendations=3,
        )
        
        # Build full output
        output = output_formatter.format_rule_output(
            reaction_smiles=reaction_smiles.strip(),
            requested_type=reaction_type if reaction_type != "Auto-detect" else None,
            detected_type=detection_result["detected_family"],
            recommendations_data=recommendations_data,
            database_name=rule_db_name,
            processing_time_ms=(time.time() - overall_start) * 1000,
        )
        
        # Format as JSON for display
        json_output = json.dumps(output, indent=2)
        
        # Create summary with detection info
        summary = f"**Rule-Based Recommendations Generated**\n\n"
        
        # Show detection info if auto-detected
        if reaction_type == "Auto-detect":
            summary += detection_message + "\n\n"
        
        summary += f"Database: {output['meta']['model']}\n"
        summary += f"Detected Type: {output['detection']['detected_reaction_type']}\n"
        summary += f"Recommendations: {len(output['recommendations'])}\n\n"
        summary += f"*Processing time: {output['meta']['processing_time_ms']:.1f} ms*\n\n"
        summary += "**JSON Output:**\n```json\n" + json_output + "\n```"
        
        # Create table view
        table = []
        for rec in output['recommendations']:
            reagent_list = []
            for r in rec['reagents']:
                if r.get('role') not in ['electrophile', 'nucleophile']:
                    name = r.get('name') or r.get('abbreviation', 'Unknown')
                    reagent_list.append(f"{name} ({r['role']})")
            
            temp = rec['conditions'].get('temperature', {}).get('value', 'N/A')
            table.append([
                rec['rank'],
                f"{rec['confidence']:.2%}",
                ", ".join(reagent_list),
                f"{temp}°C" if temp != 'N/A' else 'N/A',
            ])
        
        return summary, table
        
    except Exception as e:
        import traceback
        error_msg = f"**Rule-Based Recommendation Error:**\n\n{str(e)}\n\n```\n{traceback.format_exc()}\n```"
        return error_msg, []


def get_both_recommendations(
    reaction_smiles: str,
    reaction_type: str,
    top_k: int = 3,  # Default to 3 recommendations
) -> Tuple[str, List[List[Any]], str, str, str, List[List[Any]]]:
    """Get both ML and rule-based recommendations."""
    
    ml_summary, ml_table, ml_prec_summary, ml_prec_html = get_ml_recommendations(reaction_smiles, reaction_type, top_k)
    rule_summary, rule_table = get_rule_recommendations(reaction_smiles, reaction_type)
    
    return ml_summary, ml_table, ml_prec_summary, ml_prec_html, rule_summary, rule_table


def search_precedents(
    reaction_smiles: str,
    k: int,
    catalyst_filter: str,
    similarity_threshold: float = 0.0,
) -> Tuple[str, str]:
    """Search for precedent reactions using reaction fingerprint similarity.
    
    Universal approach that works for all reaction types without type-specific filtering.
    Uses DRFP (Differential Reaction Fingerprints) to focus on the reaction transformation.
    
    Args:
        reaction_smiles: Query reaction SMILES
        k: Number of results to return
        catalyst_filter: Optional catalyst class filter
        similarity_threshold: Minimum similarity score (0.0-1.0)
    
    Returns:
        Tuple of (summary_markdown, html_table_with_images)
    """
    
    if not reaction_smiles or not reaction_smiles.strip():
        return "**Error:** Please enter a reaction SMILES.", ""
    
    try:
        print(f"\n{'='*60}")
        print(f"UNIVERSAL PRECEDENT SEARCH (REACTION FINGERPRINT)")
        print(f"{'='*60}")
        
        start_time = time.time()
        
        # Parse reaction to get reactants
        from chemtools.smiles import normalize_reaction
        from chemtools import featurizers
        
        print(f"Normalizing reaction...")
        norm_start = time.time()
        norm = normalize_reaction(reaction_smiles.strip())
        print(f"  Normalization: {time.time() - norm_start:.2f}s")
        
        reactants = [
            (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
            for r in (norm.get("reactants") or [])
        ]
        
        if not reactants:
            return "**Error:** No reactants found in reaction SMILES.", ""
        
        # Detect family
        print(f"Detecting reaction family...")
        det_start = time.time()
        fam_result = router.detect_family(reactants)
        family = fam_result.get("family", "Unknown")
        confidence = fam_result.get("confidence", 0.0)
        print(f"  Detection: {time.time() - det_start:.2f}s")
        print(f"  Family: {family} (confidence: {confidence:.2%})")
        
        # Prepare search parameters with DRFP enabled
        relax: Dict[str, Any] = {
            "reaction_smiles": reaction_smiles.strip(),
            "use_drfp": True,  # Enable DRFP for reaction similarity
            "precompute_drfp": True,  # Use precomputed fingerprints
            "drfp_weight": 0.7,  # Weight for DRFP similarity (0.7 = 70% DRFP, 30% substrate features)
        }
        
        # Only featurize if DRFP is disabled (otherwise empty features, DRFP uses reaction_smiles)
        features: Dict[str, Any] = {}
        if not relax.get("use_drfp", False):
            # Featurize substrates only when needed (DRFP disabled)
            print(f"Featurizing substrates...")
            feat_start = time.time()
            
            # Pick electrophile and nucleophile
            def is_electrophile(s: str) -> bool:
                t = (s or "").lower()
                return ("br" in t) or ("cl" in t) or (" i" in t) or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
            
            elec, nuc = "", ""
            if len(reactants) == 1:
                elec, nuc = reactants[0], ""
            else:
                r0, r1 = reactants[0], reactants[1]
                elec, nuc = (r0, r1) if is_electrophile(r0) else ((r1, r0) if is_electrophile(r1) else (r0, r1))
            
            features = featurizers.molecular.featurize(elec, nuc)
            
            # Drop nested role_aware to keep features hashable
            if isinstance(features, dict) and "role_aware" in features:
                features = {k: v for k, v in features.items() if k != "role_aware"}
            
            print(f"  Featurization: {time.time() - feat_start:.2f}s")
        else:
            print(f"Skipping featurization (DRFP enabled - using reaction SMILES directly)")
        
        # Apply catalyst filter
        catalyst_map = {
            "Pd (palladium)": "Pd",
            "Cu (copper)": "Cu",
            "Ni (nickel)": "Ni",
        }
        
        if catalyst_filter in catalyst_map:
            relax["catalyst_class"] = catalyst_map[catalyst_filter]
        
        # Search precedents with larger k for filtering by threshold
        search_k = max(int(k) * 3, 50)  # Get more to apply threshold
        print(f"Searching precedents with DRFP (k={search_k})...")
        search_start = time.time()
        
        pack = precedent.knn(family, features, k=search_k, relax=relax)
        
        search_time = time.time() - search_start
        print(f"  Search: {search_time:.2f}s")
        
        # Extract precedents with similarity scores
        precs = pack.get("precedents", [])
        
        if not precs:
            elapsed = time.time() - start_time
            summary = f"**No Precedents Found**\n\n"
            summary += f"- **Family:** {family}\n"
            summary += f"- **Detection confidence:** {confidence:.2%}\n"
            summary += f"- **Search time:** {elapsed:.1f}s\n\n"
            summary += "*Try adjusting the similarity threshold or catalyst filter.*"
            return summary, ""
        
        print(f"  Initial results: {len(precs)} precedents")
        
        # Apply similarity threshold filtering (universal for all reaction types)
        if similarity_threshold > 0.0:
            filtered_precs = []
            for prec in precs:
                # Get similarity score (added by precedent.knn when DRFP is enabled)
                sim_score = prec.get("similarity", 0.0)
                if sim_score >= similarity_threshold:
                    filtered_precs.append(prec)
            
            if filtered_precs:
                precs = filtered_precs
                print(f"  ✓ Filtered by similarity ≥ {similarity_threshold:.2f}: {len(precs)} precedents remain")
            else:
                print(f"  ⚠ No precedents above similarity threshold {similarity_threshold:.2f}, showing all")
        
        # Limit to requested k
        precs = precs[:int(k)]
        print(f"  ✓ Returning top {len(precs)} precedents")
        
        # Show similarity score distribution
        if precs:
            scores = [p.get("similarity", 0.0) for p in precs]
            if any(s > 0 for s in scores):
                print(f"  Similarity range: {min(scores):.3f} to {max(scores):.3f}")
        
        # Build summary
        elapsed = time.time() - start_time
        summary = f"**Universal Precedent Search Results**\n\n"
        summary += f"- **Detected family:** {family} (confidence: {confidence:.2%})\n"
        summary += f"- **Precedents found:** {len(precs)}\n"
        summary += f"- **Similarity method:** DRFP (Differential Reaction Fingerprints)\n"
        summary += f"- **Search strategy:** Reaction-level similarity (works for all reaction types)\n"
        
        if similarity_threshold > 0.0:
            summary += f"- **Similarity threshold:** ≥ {similarity_threshold:.2f}\n"
        
        if precs:
            scores = [p.get("similarity", 0.0) for p in precs]
            if any(s > 0 for s in scores):
                summary += f"- **Similarity range:** {min(scores):.3f} - {max(scores):.3f}\n"
        
        summary += f"- **Search time:** {elapsed:.1f}s\n"
        
        if catalyst_filter != "Any (no filter)":
            summary += f"- **Catalyst filter:** {catalyst_filter}\n"
        
        # Helper functions for rendering reaction images
        def _img_data_uri(img_obj: Any) -> str | None:
            """Convert PIL image to data URI."""
            try:
                import base64
                import io
                buf = io.BytesIO()
                img_obj.save(buf, format="PNG")
                b64 = base64.b64encode(buf.getvalue()).decode("ascii")
                return f"data:image/png;base64,{b64}"
            except Exception:
                return None
        
        def _render_img(reactants_smi: str, products_smi: str) -> str:
            """Render reaction scheme with arrow between reactants and products."""
            try:
                from chemtools.util.rdkit_helpers import rdkit_available as _rd_avail
                if not _rd_avail():
                    return ""
                from rdkit import Chem
                from rdkit.Chem import Draw
                
                def _grid(smi: str):
                    """Create grid image from SMILES."""
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
                    return Draw.MolsToGridImage(ms, molsPerRow=min(3, len(ms)), subImgSize=(220, 220))
                
                # Generate images for reactants and products
                left_img = _grid(reactants_smi)
                right_img = _grid(products_smi)
                
                if left_img is None and right_img is None:
                    return ""
                if left_img is None:
                    uri = _img_data_uri(right_img)
                    return f"<img src='{uri}' width='340'/>" if uri else ""
                if right_img is None:
                    uri = _img_data_uri(left_img)
                    return f"<img src='{uri}' width='340'/>" if uri else ""
                
                # Combine images with arrow
                try:
                    from PIL import Image as _Image, ImageDraw as _ImageDraw
                    
                    # Create canvas with arrow space
                    arrow_width = 60
                    w = left_img.width + right_img.width + arrow_width
                    h = max(left_img.height, right_img.height)
                    canvas = _Image.new("RGB", (w, h), (255, 255, 255))
                    
                    # Paste reactants and products
                    canvas.paste(left_img, (0, (h - left_img.height) // 2))
                    canvas.paste(right_img, (left_img.width + arrow_width, (h - right_img.height) // 2))
                    
                    # Draw arrow
                    dr = _ImageDraw.Draw(canvas)
                    y = h // 2
                    arrow_start_x = left_img.width + 10
                    arrow_end_x = left_img.width + arrow_width - 10
                    
                    # Arrow line
                    dr.line((arrow_start_x, y, arrow_end_x, y), fill=(0, 0, 0), width=3)
                    
                    # Arrow head
                    arrow_head_size = 6
                    dr.polygon([
                        (arrow_end_x, y),
                        (arrow_end_x - arrow_head_size, y - arrow_head_size),
                        (arrow_end_x - arrow_head_size, y + arrow_head_size)
                    ], fill=(0, 0, 0))
                    
                    uri = _img_data_uri(canvas)
                    return f"<img src='{uri}' width='340'/>" if uri else ""
                except Exception:
                    # Fallback to just showing reactants
                    uri = _img_data_uri(left_img)
                    return f"<img src='{uri}' width='340'/>" if uri else ""
            except Exception:
                return ""
        
        # Build HTML table with reaction images
        html_rows: List[str] = []
        html_rows.append('<table style="border-collapse:collapse; width:100%; table-layout:fixed">')
        html_rows.append(
            "<colgroup>"
            "<col style='width:340px'/>"    # Reaction image
            "<col style='width:70px'/>"     # Similarity
            "<col style='width:100px'/>"    # ID
            "<col style='width:65px'/>"     # Yield
            "<col style='width:150px'/>"    # Catalyst
            "<col style='width:115px'/>"    # Base
            "<col style='width:115px'/>"    # Solvent
            "<col style='width:65px'/>"     # Temp
            "<col style='width:65px'/>"     # Time
            "</colgroup>"
        )
        html_rows.append(
            "<tr style='background-color:#f0f0f0; font-weight:bold'>"
            "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Reaction</th>"
            "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Similarity</th>"
            "<th style='padding:8px; text-align:left; border:1px solid #ddd'>ID</th>"
            "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Yield</th>"
            "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Catalyst</th>"
            "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Base</th>"
            "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Solvent</th>"
            "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Temp °C</th>"
            "<th style='padding:8px; text-align:left; border:1px solid #ddd'>Time h</th>"
            "</tr>"
        )
        
        for i, prec in enumerate(precs, 1):
            # Extract SMILES for rendering
            reaction_smi = prec.get("reaction_smiles", "")
            reactants_smi = prec.get("reactants_smiles", "")
            products_smi = prec.get("products_smiles", "")
            
            # If reactants/products not separate, try to split reaction_smiles
            if not reactants_smi or not products_smi:
                if ">>" in reaction_smi:
                    parts = reaction_smi.split(">>")
                    reactants_smi = parts[0] if len(parts) > 0 else ""
                    products_smi = parts[1] if len(parts) > 1 else ""
            
            # Render reaction image
            img_html = _render_img(reactants_smi, products_smi)
            
            # Extract similarity score
            sim_score = prec.get("similarity", 0.0)
            sim_display = f"{sim_score:.3f}" if sim_score > 0 else "-"
            sim_color = "#2e7d32" if sim_score >= 0.8 else "#f57c00" if sim_score >= 0.6 else "#666"
            
            # Extract data
            reaction_id = prec.get("reaction_id", "")
            yield_val = prec.get("yield", "")
            core = prec.get("core", "")
            base_uid = prec.get("base_uid", "")
            solvent_uid = prec.get("solvent_uid", "")
            temp = prec.get("T_C", "")
            time_h = prec.get("time_h", "")
            
            # Try to get readable labels
            try:
                base_label = reagent_lookup.enrich_reagent_info(base_uid, "base").get("name", base_uid) if base_uid else ""
            except Exception:
                base_label = base_uid
            
            try:
                solvent_label = reagent_lookup.enrich_reagent_info(solvent_uid, "solvent").get("name", solvent_uid) if solvent_uid else ""
            except Exception:
                solvent_label = solvent_uid
            
            # Add table row
            row_style = "background-color:#fafafa" if i % 2 == 0 else ""
            html_rows.append(
                f"<tr style='{row_style}'>"
                f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top'>{img_html}</td>"
                f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:nowrap; text-align:center; font-weight:bold; color:{sim_color}'>{sim_display}</td>"
                f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{reaction_id[:35]}{'...' if len(reaction_id) > 35 else ''}</td>"
                f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:nowrap; text-align:center'>{yield_val}</td>"
                f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{core}</td>"
                f"<td title='{base_uid}' style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{base_label}</td>"
                f"<td title='{solvent_uid}' style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:normal; word-break:break-word; overflow-wrap:anywhere'>{solvent_label}</td>"
                f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:nowrap; text-align:center'>{temp}</td>"
                f"<td style='padding:8px; border:1px solid #ddd; vertical-align:top; white-space:nowrap; text-align:center'>{time_h}</td>"
                "</tr>"
            )
        
        html_rows.append("</table>")
        html = "\n".join(html_rows)
        
        print(f"Total time: {elapsed:.1f}s")
        print(f"{'='*60}\n")
        
        return summary, html
        
    except Exception as e:
        import traceback
        
        error_msg = f"**Precedent Search Error**\n\n"
        error_msg += f"❌ **Error Type:** `{type(e).__name__}`\n"
        error_msg += f"**Message:** {str(e)}\n\n"
        error_msg += "**Technical Details:**\n```\n" + traceback.format_exc() + "\n```"
        
        return error_msg, ""


# ============================================================================
# Gradio UI
# ============================================================================

def create_ui():
    """Create the Gradio interface."""
    
    with gr.Blocks(
        theme=THEME,
        title="Condition Agent - Reaction Recommendations",
        css="""
        .compact-table { 
            font-size: 0.9em; 
        }
        .summary-box { 
            padding: 1em; 
            background: #f8f9fa; 
            border-radius: 0.5em; 
            color: #1f2937 !important;
        }
        .summary-box p {
            color: #1f2937 !important;
        }
        .summary-box * {
            color: #1f2937 !important;
        }
        .markdown-text {
            color: #1f2937 !important;
        }
        /* Ensure all markdown content has good contrast */
        .gradio-container .prose {
            color: #1f2937 !important;
        }
        .gradio-container .prose * {
            color: inherit !important;
        }
        .gradio-container .prose h1,
        .gradio-container .prose h2,
        .gradio-container .prose h3,
        .gradio-container .prose h4 {
            color: #111827 !important;
        }
        .gradio-container .prose p,
        .gradio-container .prose li,
        .gradio-container .prose ul,
        .gradio-container .prose ol {
            color: #374151 !important;
        }
        .gradio-container .prose strong,
        .gradio-container .prose b {
            color: #1f2937 !important;
        }
        .gradio-container .prose code {
            background: #f3f4f6;
            color: #1f2937 !important;
            padding: 0.2em 0.4em;
            border-radius: 0.25em;
        }
        .gradio-container .prose pre {
            background: #f3f4f6 !important;
            color: #1f2937 !important;
        }
        .gradio-container .prose pre code {
            background: transparent !important;
        }
        """
    ) as demo:
        
        gr.Markdown(
            """
            # 🧪 Condition Agent - Reaction Recommendations
            
            Get condition recommendations using:
            - **ML-based** (DRFP similarity search) - Data-driven recommendations with precedents shown
            - **Rule-based** (SchemeConditionDB pattern matching) - Fast expert rules
            
            *ML recommendations automatically show the literature precedents used, with full reaction schemes.*
            """
        )
        
        # Input section
        with gr.Row():
            with gr.Column(scale=3):
                reaction_input = gr.Textbox(
                    label="Reaction SMILES",
                    placeholder="Enter reaction SMILES (e.g., Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1)",
                    lines=2,
                )
            
            with gr.Column(scale=1):
                reaction_type = gr.Dropdown(
                    label="Reaction Type",
                    choices=list(ML_FAMILY_MAP.keys()),
                    value="Auto-detect",
                )
                
                top_k = gr.Slider(
                    label="Number of ML Results",
                    minimum=1,
                    maximum=10,
                    value=3,
                    step=1,
                )
        
        # Example reactions
        with gr.Row():
            gr.Examples(
                examples=[
                    ["Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1", "C-N Coupling (Pd)"],
                    ["Clc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1", "C-N Coupling (Cu)"],
                    ["Brc1ccccc1.NC1CCCC1>>c1ccc(N2CCCC2)cc1", "C-N Coupling (Ni)"],
                    ["CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1", "Amide Formation"],
                    ["Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1", "Suzuki Coupling"],
                ],
                inputs=[reaction_input, reaction_type],
                label="Example Reactions",
            )
        
        # Action buttons
        with gr.Row():
            ml_btn = gr.Button("🤖 ML Recommendations", variant="primary", scale=1)
            rule_btn = gr.Button("📋 Rule-Based Recommendations", variant="primary", scale=1)
            both_btn = gr.Button("🔄 Both Methods", variant="secondary", scale=1)
            clear_btn = gr.ClearButton(scale=1)
        
        # ML Results
        gr.Markdown("## 🤖 ML-Based Recommendations")
        with gr.Row():
            with gr.Column():
                ml_summary = gr.Markdown(
                    label="ML Summary", 
                    elem_classes=["summary-box"],
                    value=""
                )
            
        with gr.Row():
            ml_table = gr.Dataframe(
                headers=["Rank", "Confidence", "Reagents", "Temperature"],
                label="ML Recommendations Summary",
                interactive=False,
                elem_classes=["compact-table"],
            )
        
        # Precedents used by ML
        gr.Markdown("### 📚 Precedents Used by ML System")
        with gr.Row():
            with gr.Column():
                ml_precedent_summary = gr.Markdown(
                    label="Precedent Summary",
                    elem_classes=["summary-box"],
                    value=""
                )
        
        with gr.Row():
            ml_precedent_table = gr.HTML(
                label="Similar Reactions from Literature (with structures)",
                value="<p style='color:#666; padding:20px; text-align:center'>Precedents will appear here after running ML recommendations</p>"
            )
        
        # Rule-Based Results
        gr.Markdown("---")
        gr.Markdown("## 📋 Rule-Based Recommendations")
        with gr.Row():
            with gr.Column():
                rule_summary = gr.Markdown(
                    label="Rule Summary", 
                    elem_classes=["summary-box"],
                    value=""
                )
        
        with gr.Row():
            rule_table = gr.Dataframe(
                headers=["Rank", "Confidence", "Reagents", "Temperature"],
                label="Rule-Based Recommendations Summary",
                interactive=False,
                elem_classes=["compact-table"],
            )
        
        # Wire up events
        ml_btn.click(
            fn=get_ml_recommendations,
            inputs=[reaction_input, reaction_type, top_k],
            outputs=[ml_summary, ml_table, ml_precedent_summary, ml_precedent_table],
        )
        
        rule_btn.click(
            fn=get_rule_recommendations,
            inputs=[reaction_input, reaction_type],
            outputs=[rule_summary, rule_table],
        )
        
        both_btn.click(
            fn=get_both_recommendations,
            inputs=[reaction_input, reaction_type, top_k],
            outputs=[ml_summary, ml_table, ml_precedent_summary, ml_precedent_table, rule_summary, rule_table],
        )
        
        clear_btn.add([reaction_input, ml_summary, ml_table, ml_precedent_summary, ml_precedent_table, rule_summary, rule_table])
        
        # Footer
        gr.Markdown(
            """
            ---
            
            ### 💡 How to Use
            
            #### Recommendation Tabs
            1. **Enter a reaction SMILES** in the format: `Reactants>>Product`
            2. **Select reaction type** or use Auto-detect
            3. **Click a button** to get recommendations:
               - **ML Recommendations**: DRFP-based precedent search (~5-10 seconds) - Shows precedents used
               - **Rule-Based**: Pattern-matched expert rules (very fast, <1 second)
               - **Both**: Get comprehensive recommendations from both methods
            
            **ML Recommendations automatically display:**
            - Recommended conditions (catalyst, base, solvent, temperature, time)
            - Literature precedents used by the ML system (with full reaction schemes)
            - DRFP similarity scores showing how similar each precedent is to your query
            
            ### 🔬 About DRFP (Differential Reaction Fingerprints)
            
            **What is DRFP?** A universal reaction similarity metric that focuses on the chemical transformation 
            rather than just substrate structures. It works across all reaction types without custom rules.
            
            **Why we use it:**
            - ✅ **Universal**: No need for reaction-specific filtering or rules
            - ✅ **Chemical**: Captures the actual transformation, not just reactant similarity
            - ✅ **Scalable**: Same approach works for Suzuki, Buchwald-Hartwig, Ullmann, and more
            - ✅ **Transparent**: Similarity scores (0.0-1.0) show match quality
            
            **Interpreting similarity scores:**
            - **0.8-1.0**: Very similar transformation (green in precedent table)
            - **0.6-0.8**: Moderately similar (orange in precedent table)
            - **<0.6**: Different transformation (gray in precedent table)
            
            ### 📚 Supported Reaction Types
            
            - **C-N Coupling (Cu)**: Ullmann-type reactions
            - **C-N Coupling (Pd)**: Buchwald-Hartwig reactions
            - **C-N Coupling (Ni)**: Nickel-catalyzed C-N couplings
            - **Amide Formation**: Peptide coupling, amidations
            - **Suzuki Coupling**: Suzuki-Miyaura cross-coupling
            
            ### ⚡ Performance Notes
            
            - **ML recommendations** use DRFP with precomputed fingerprints (~5-10 seconds total)
            - **Precedents** are retrieved automatically as part of ML recommendations
            - **Rule-based recommendations** are nearly instant (<1 second)
            - **First ML query** takes ~60s to load datasets (one-time startup cost)
            - **Subsequent queries** are very fast thanks to in-memory caching
            - **Tip**: Run a test query when UI first loads to warm up the system!
            """
        )
    
    return demo


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("STARTING SIMPLIFIED UI - CONDITION AGENT")
    print("="*70)
    
    # Pre-load datasets at startup to avoid 157s delay on first query
    print("\n[1/2] Pre-loading reaction datasets...")
    datasets_loaded = preload_datasets()
    
    # Create and launch UI
    print("[2/2] Starting Gradio interface...")
    demo = create_ui()
    
    print("\n" + "="*70)
    if datasets_loaded:
        print("✓ UI READY - All datasets pre-loaded in memory!")
        print("  First query will be fast (~0.5-2s)")
    else:
        print("⚠ UI READY - Datasets will load on first query")
        print("  First query may take ~60-120s")
    print("="*70 + "\n")
    
    demo.launch(
        server_name="127.0.0.1",
        server_port=7861,
        share=False,
    )
