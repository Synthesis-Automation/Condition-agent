#!/usr/bin/env python3

"""PyQt6 UI for the reagent registry generator workflow."""

from __future__ import annotations

import json
import sys
from copy import deepcopy
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from PyQt6.QtCore import QObject, QRunnable, QThreadPool, Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

MODULE_DIR = Path(__file__).resolve().parent
ROOT_DIR = MODULE_DIR.parent

# Add root directory to path so chemtools can be imported
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# Import from the new unified chemtools.reagent package
from chemtools.reagent import (
    dedupe_synonyms,
    load_families_registry_entries,
    normalize_cas,
    resolve_identity_from_cas,
    tokenize_all,
    ROLE_ALIASES,
)

DEFAULT_RESOLVER_TIMEOUT = 6.0  # Default timeout for identity resolution

# Default to the unified flattened registry directory used across the project
# Aligns with chemtools.reagent defaults (data/reagent_db)
DEFAULT_REGISTRY_DIR = (ROOT_DIR / "data" / "reagent_db").resolve()


def canonical_role(role: str) -> str:
    """Normalize legacy role labels to the canonical registry keys."""
    return ROLE_ALIASES.get(role, role)

DEFAULT_LLM_MODELS: Dict[str, List[str]] = {
    "aliyun": [
        "deepseek-r1-distill-qwen-7b",
        "deepseek-v3.2",
        "deepseek-v3.1",
        "deepseek-r1",
        "deepseek-r1-0528",
        "deepseek-v3",
        "deepseek-r1-distill-qwen-1.5b",
        "deepseek-r1-distill-qwen-14b",
        "deepseek-r1-distill-qwen-32b",
        "deepseek-r1-distill-llama-8b",
        "deepseek-r1-distill-llama-70b",
    ],
    "openai": [
        "gpt-5",
        "gpt-5-pro",
        "gpt-5-mini",
        "gpt-5-nano",
        "gpt-5-codex",
        "o3",
        "o3-pro",
        "o3-mini",
        "o4-mini",
        "o3-deep-research",
        "o4-mini-deep-research",
        "gpt-4o",
        "gpt-4o-mini",
        "gpt-4.1-mini",
        "gpt-4.1-nano",
    ],
}

DEFAULT_LLM_RECOMMENDED: Dict[str, Dict[str, str]] = {
    "aliyun": {
        "reasoning": "deepseek-r1",
        "fast": "deepseek-r1-distill-qwen-7b",
        "balanced": "deepseek-v3.2",
        "experimental": "deepseek-v3.2",
    },
    "openai": {
        "reasoning": "o3-mini",
        "fast": "gpt-4o",
        "balanced": "gpt-4o",
        "advanced": "gpt-5-mini",
    },
}

LLM_SUPPORT_ERROR: Optional[str] = None

try:
    from llmtools.clients import AVAILABLE_MODELS as LLM_AVAILABLE_MODELS
    from llmtools.clients import RECOMMENDED_MODELS as LLM_RECOMMENDED_MODELS
except Exception as exc:  # pragma: no cover - safeguard for missing optional deps
    LLM_AVAILABLE_MODELS = DEFAULT_LLM_MODELS
    LLM_RECOMMENDED_MODELS = DEFAULT_LLM_RECOMMENDED
    LLM_SUPPORT_ERROR = f"{exc.__class__.__name__}: {exc}"

try:
    from llmtools.reagent_review import LLMReviewOptions, review_taxonomy_proposal
    from llmtools.reagent_classifier import classify_role, assign_fields, verify_entry, VALID_ROLES
    from llmtools.clients import LLMClient
    LLM_SUPPORT_AVAILABLE = True
except Exception as exc:  # pragma: no cover - optional dependency
    LLMReviewOptions = None  # type: ignore[assignment]
    review_taxonomy_proposal = None  # type: ignore[assignment]
    classify_role = None  # type: ignore[assignment]
    assign_fields = None  # type: ignore[assignment]
    verify_entry = None  # type: ignore[assignment]
    VALID_ROLES = None  # type: ignore[assignment]
    LLMClient = None  # type: ignore[assignment]
    LLM_SUPPORT_AVAILABLE = False
    if not LLM_SUPPORT_ERROR:
        LLM_SUPPORT_ERROR = f"{exc.__class__.__name__}: {exc}"

ROLE_CONFIG: Dict[str, Dict[str, Any]] = {
    "ligand": {
        "filename": "ligand.json",
        "label": "Ligand",
        "hint": "Ligands including phosphines, NHCs, diimines, and ancillary donor sets.",
        "priority": 0,
        "default_family": "trialkyl_triaryl_phosphines",
    },
    "metal_catalyst": {
        "filename": "metal_catalyst.json",
        "label": "Metal catalyst",
        "hint": "Metal salts or pre-ligated complexes that generate the active catalyst.",
        "priority": 1,
        "default_family": "pd_ii_salts",
    },
    "base": {
        "filename": "base.json",
        "label": "Base",
        "hint": "Bronsted/Lewis bases spanning amides, alkoxides, carbonates, superbases.",
        "priority": 3,
        "default_family": "tertiary_amines_aliphatic",
    },
    "acid": {
        "filename": "acid.json",
        "label": "Acid (Bronsted/Lewis)",
        "hint": "Mineral, sulfonic, and superacids used as activators or promoters.",
        "priority": 4,
        "default_family": "mineral_acids",
    },
    "condensation_agent": {
        "filename": "condensation_agent.json",
        "label": "Condensation agent",
        "hint": "Carbodiimides, uronium/phosphonium activators, and similar condensers.",
        "priority": 5,
        "default_family": "carbodiimides",
    },
    "oxidant": {
        "filename": "oxidant.json",
        "label": "Oxidant",
        "hint": "Terminal and co-oxidants (peroxides, hypervalent iodine, oxone, O2, etc.).",
        "priority": 6,
        "default_family": "o2_gas",
    },
    "reductant": {
        "filename": "reductant.json",
        "label": "Reductant",
        "hint": "Hydrides, silanes, metal reductants, and organic electron donors.",
        "priority": 7,
        "default_family": "metal_powders",
    },
    "additive": {
        "filename": "additive.json",
        "label": "Additive / Modulator",
        "hint": "Phase-transfer agents, halide scavengers, fluoride sources, and related modifiers.",
        "priority": 8,
        "default_family": "quaternary_ammonium_ptc",
    },
    "solvent": {
        "filename": "solvent.json",
        "label": "Solvent",
        "hint": "Reaction media categorized by polarity, coordinating ability, and safety profile.",
        "priority": 9,
        "default_family": "hydrocarbons_aromatic",
    },
    "organo_catalyst": {
        "filename": "organo_catalyst.json",
        "label": "Organo-catalyst",
        "hint": "Small-molecule catalysts (cinchona, phosphoric acids, NHCs) providing organocatalysis.",
        "priority": 10,
        "default_family": "cinchona_alkaloids",
    },
    "enzyme": {
        "filename": "enzyme.json",
        "label": "Enzyme",
        "hint": "Biocatalysts supplied as isolated enzymes or whole-cell systems.",
        "priority": 11,
        "default_family": "oxidoreductase_general",
    },
    "other_reagent": {
        "filename": "other_reagent.json",
        "label": "Other Reagent",
        "hint": "Generic reagents or supports that do not fit existing roles.",
        "priority": 10,
        "default_family": "misc_general",
    },
}

ROLE_EMBED_FIELDS: Dict[str, Sequence[str]] = {
    "ligand": ("donors", "denticity"),
    "metal_catalyst": ("metal", "oxidation_states", "ligand_type"),
    "base": ("basicity", "nucleophilicity", "sterics"),
    "acid": ("acidity",),
    "condensation_agent": ("strength_band",),
    "oxidant": ("strength_band",),
    "reductant": ("strength_band",),
    "solvent": ("proticity", "polarity", "polarity_band", "coordination"),
    "organo_catalyst": ("activation_mode", "chirality"),
    "enzyme": ("source", "cofactor_requirement"),
    "other_reagent": (),
}

ROLE_FIELD_LABELS: Dict[str, str] = {
    "donors": "donor atoms",
    "denticity": "denticity",
    "basicity": "basicity",
    "nucleophilicity": "nucleophilicity",
    "sterics": "sterics",
    "strength_band": "strength",
    "proticity": "proticity",
    "polarity": "polarity",
    "polarity_band": "polarity band",
    "coordination": "coordination",
    "metal": "metal",
    "oxidation_states": "oxidation states",
    "ligand_type": "ligand type",
    "acidity": "acidity",
    "activation_mode": "activation mode",
    "chirality": "chirality",
    "source": "source",
    "cofactor_requirement": "cofactor requirement",
}

ROLE_PAYLOAD_FIELDS: Dict[str, Sequence[str]] = {
    "additive": (),
    "base": ("basicity", "nucleophilicity", "sterics"),
    "condensation_agent": ("strength_band",),
    "ligand": ("donors", "denticity"),
    "oxidant": ("strength_band",),
    "metal_catalyst": ("metal", "oxidation_states", "ligand_type"),
    "reductant": ("strength_band",),
    "solvent": ("proticity", "polarity", "coordination"),
    "organo_catalyst": ("activation_mode", "chirality"),
    "enzyme": ("source", "cofactor_requirement"),
    "other_reagent": (),
}



class ReagentRegistryStore:
    """Load and mutate the reagent registry using the new flattened schema."""

    def __init__(self, base_dir: Path) -> None:
        self.base_dir = Path(base_dir)
        if not self.base_dir.exists():
            raise FileNotFoundError(f"Registry directory {self.base_dir} does not exist")
        self.registry_files: Dict[str, Path] = {}
        for role, cfg in ROLE_CONFIG.items():
            filename = cfg.get("filename")
            if not filename:
                continue
            self.registry_files[role] = self.base_dir / filename
        self.role_entries: Dict[str, List[Dict[str, Any]]] = {}
        self.cas_index: Dict[str, List[Tuple[str, Dict[str, Any]]]] = {}
        self.family_index: Dict[Tuple[str, str], Dict[str, Any]] = {}
        self.family_tokens: Dict[Tuple[str, str], Set[str]] = {}
        self.family_examples: Dict[Tuple[str, str], Set[str]] = {}
        self._load_families()
        self._load_registry()

    def _load_families(self) -> None:
        entries: List[Dict[str, Any]] = []
        schema_path = self.base_dir / "reagent_schema" / "families_registry.json"

        if schema_path.exists():
            data = json.loads(schema_path.read_text(encoding="utf-8"))
            entries = data.get("entries", [])
        else:
            try:
                entries = load_families_registry_entries()
            except FileNotFoundError as exc:
                raise FileNotFoundError(
                    f"Families registry not found in {self.base_dir / 'reagent_schema'} "
                    f"and unable to load unified taxonomy: {exc}"
                ) from exc

        for entry in entries:
            role = entry.get("role")
            family = entry.get("family")
            if not role or not family:
                continue
            key = (role, family)
            self.family_index[key] = entry
            tokens: Set[str] = set()
            tokens.update(tokenize_all(entry.get("keywords", [])))
            tokens.update(tokenize_all([family, entry.get("definition", ""), entry.get("notes", "")]))
            self.family_tokens[key] = {tok for tok in tokens if tok}
            examples = {str(cas) for cas in entry.get("examples_pos", []) if cas}
            self.family_examples[key] = examples

    def _load_registry(self) -> None:
        for role, path in self.registry_files.items():
            if not path.exists():
                self.role_entries[role] = []
                continue
            content = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(content, list):
                raise ValueError(f"Registry file {path} must contain a list of entries.")
            self.role_entries[role] = content
            for entry in content:
                cas = entry.get("cas")
                if not cas:
                    continue
                self.cas_index.setdefault(str(cas), []).append((role, entry))

    def file_for_role(self, role: str) -> Path:
        cfg = ROLE_CONFIG.get(role)
        if not cfg or not cfg.get("filename"):
            raise KeyError(f"Unknown or unsupported role '{role}'")
        return self.registry_files.get(role, self.base_dir / cfg["filename"])

    def family_entry(self, role: str, family_id: str) -> Optional[Dict[str, Any]]:
        return self.family_index.get((role, family_id))

    def default_family(self, role: str) -> Optional[str]:
        cfg = ROLE_CONFIG.get(role)
        return cfg.get("default_family") if cfg else None

    def build_role_payload(self, role: str, family_id: str) -> Dict[str, Any]:
        payload: Dict[str, Any] = {"families": [family_id]}
        family = self.family_entry(role, family_id)
        if family:
            required = family.get("required_props") or {}
            for key, value in required.items():
                payload[key] = value
        expected_fields = ROLE_PAYLOAD_FIELDS.get(role, ())
        defaults = self._default_role_fields(role, family_id, expected_fields)
        for field in expected_fields:
            if field in payload:
                continue
            if field in defaults:
                payload[field] = defaults[field]
            else:
                payload[field] = None
        return payload

    def _default_role_fields(
        self, role: str, family_id: str, expected_fields: Sequence[str]
    ) -> Dict[str, Any]:
        entries = self.role_entries.get(role, [])
        for entry in entries:
            role_payload = entry.get("roles", {}).get(role, {})
            families = role_payload.get("families") or []
            if family_id not in families:
                continue
            defaults: Dict[str, Any] = {}
            for field in expected_fields:
                if field in role_payload:
                    defaults[field] = deepcopy(role_payload[field])
            if defaults:
                return defaults
        return {}

    def infer_family(self, role: str, cas: str, tokens: Set[str]) -> Optional[Tuple[str, List[str], Dict[str, Any]]]:
        cas_key = str(cas)
        best: Optional[Tuple[Tuple[int, int, int], str, List[str], Dict[str, Any]]] = None
        for (fam_role, fam_id), entry in self.family_index.items():
            if fam_role != role:
                continue
            candidate_tokens = self.family_tokens.get((fam_role, fam_id), set())
            matches = sorted(tokens & candidate_tokens)
            cas_match = 1 if cas_key in self.family_examples.get((fam_role, fam_id), set()) else 0
            if not matches and not cas_match:
                continue
            precedence = int(entry.get("precedence", 0) or 0)
            score = (cas_match, len(matches), precedence)
            reason = ["cas_match"] if cas_match else matches
            if not best or score > best[0]:
                best = (score, fam_id, reason, entry)
        if best:
            _, fam_id, reason, entry = best
            return fam_id, reason, entry
        return None

    def find_by_cas(self, cas: str, role: Optional[str] = None) -> Optional[Tuple[str, Optional[str], Dict[str, Any]]]:
        entries = self.cas_index.get(str(cas), [])
        for entry_role, entry in entries:
            if role and entry_role != role:
                continue
            role_data = entry.get("roles", {}).get(entry_role, {})
            families = role_data.get("families") or []
            family_id = families[0] if families else None
            return entry_role, family_id, entry
        return None

    def add_entry(self, role: str, entry: Dict[str, Any]) -> None:
        entries = self.role_entries.setdefault(role, [])
        entries.append(entry)
        entries.sort(key=lambda item: (item.get("name") or "").lower())
        cas = entry.get("cas")
        if cas:
            self.cas_index.setdefault(str(cas), []).append((role, entry))

    def save_role(self, role: str) -> Path:
        path = self.file_for_role(role)
        path.parent.mkdir(parents=True, exist_ok=True)
        entries = self.role_entries.get(role, [])
        path.write_text(json.dumps(entries, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return path


NUMERIC_RUN_PATTERN = re.compile(r"\d{3,}")

def _looks_like_abbreviation(value: str) -> bool:
    clean = value.strip()
    if not clean or len(clean) < 2:
        return False
    uppercase = sum(1 for ch in clean if ch.isalpha() and ch.isupper())
    if clean.isupper() and len(clean) <= 10:
        return True
    if len(clean) <= 6 and uppercase >= 2:
        return True
    if len(clean) <= 6 and any(ch.isdigit() for ch in clean):
        return True
    if len(clean) <= 7 and '-' in clean:
        return True
    return False


def _has_long_digit_run(value: str) -> bool:
    return bool(NUMERIC_RUN_PATTERN.search(value))


def infer_abbreviations(name: str, synonyms: Sequence[str]) -> List[str]:
    candidates: List[str] = []
    seen: Set[str] = set()
    name_key = name.lower()
    for syn in synonyms:
        clean = syn.strip()
        if not clean:
            continue
        if _has_long_digit_run(clean):
            continue
        key = clean.lower()
        if key == name_key or key in seen:
            continue
        if _looks_like_abbreviation(clean):
            candidates.append(clean)
            seen.add(key)
    if not candidates:
        return [name]
    if all(candidate.lower() != name_key for candidate in candidates):
        candidates.append(name)
    return candidates


def build_aliases(name: str, cas: str, abbreviations: Sequence[str], synonyms: Sequence[str]) -> List[str]:
    name_key = name.lower()
    cas_key = cas.lower()
    abbr_keys = {abbr.lower() for abbr in abbreviations}
    aliases: List[str] = []
    seen: Set[str] = set()
    for syn in synonyms:
        clean = syn.strip()
        if not clean:
            continue
        key = clean.lower()
        if _has_long_digit_run(clean):
            continue
        if key in seen or key == name_key or key == cas_key or key in abbr_keys:
            continue
        aliases.append(clean)
        seen.add(key)
    return aliases


def _format_embedding_value(value: Any) -> str:
    if isinstance(value, (list, tuple, set)):
        return ", ".join(str(item) for item in value if item or item == 0)
    return str(value)


def build_embedding_text(role: str, family_entry: Dict[str, Any], entry: Dict[str, Any], synonyms: Sequence[str]) -> str:
    pieces: List[str] = [f"type: {role.upper()}"]
    family_label = family_entry.get("definition") or family_entry.get("family") or ""
    family_id = family_entry.get("family", "")
    if family_label:
        pieces.append(f"family: {family_label} ({family_id})")
    else:
        pieces.append(f"family: {family_id}")
    role_payload = entry.get("roles", {}).get(role, {})
    for field in ROLE_EMBED_FIELDS.get(role, ()):  # primary fields first
        value = role_payload.get(field)
        if value is not None:
            pieces.append(f"{ROLE_FIELD_LABELS.get(field, field)}: {_format_embedding_value(value)}")
    for field, value in role_payload.items():
        if field == "families" or field in ROLE_EMBED_FIELDS.get(role, ()):  # already emitted
            continue
        if value is not None:
            pieces.append(f"{ROLE_FIELD_LABELS.get(field, field)}: {_format_embedding_value(value)}")
    pieces.append(f"name: {entry.get('name', '')}")
    abbreviations = entry.get("abbreviation") or []
    if abbreviations:
        pieces.append("abbr: " + "; ".join(str(item) for item in abbreviations))
    pieces.append(f"CAS: {entry.get('cas', '')}")
    if synonyms:
        pieces.append("synonyms: " + "; ".join(synonyms))
    notes = family_entry.get("notes")
    if notes:
        pieces.append("family_notes: " + notes)
    return " | ".join(pieces)


def build_registry_entry(
    *,
    entry_id: str,
    name: str,
    abbreviations: Sequence[str],
    aliases: Sequence[str],
    cas: str,
    smiles: Optional[str],
    inchi_key: Optional[str],
    role: str,
    role_payload: Dict[str, Any],
    family_entry: Dict[str, Any],
    synonyms: Sequence[str],
) -> Dict[str, Any]:
    entry: Dict[str, Any] = {
        "id": entry_id,
        "name": name,
        "abbreviation": list(abbreviations) if abbreviations else [],
        "aliases": list(aliases) if aliases else [],
        "cas": cas,
        "inchi_key": inchi_key,
        "smiles": smiles,
        "roles": {role: role_payload},
    }
    entry["embedding_text"] = build_embedding_text(role, family_entry, entry, synonyms)
    return entry



class RegistryGenerationError(RuntimeError):
    """Raised when registry generation fails for the requested input."""


def generate_taxonomy_entry_llm(
    *,
    cas: str,
    registry_dir: Path,
    llm_client: Any,  # LLMClient from llmtools.clients
    name_override: Optional[str] = None,
    smiles_override: Optional[str] = None,
    role_override: Optional[str] = None,
    resolver_timeout: float = DEFAULT_RESOLVER_TIMEOUT,
) -> Dict[str, Any]:
    """
    Pure LLM workflow for reagent taxonomy generation.
    
    Workflow:
    1. Resolve chemical identity from CAS (PubChem API)
    2. LLM classifies reagent role (or use role_override)
    3. LLM assigns family and role-specific fields
    4. LLM verifies entry for errors
    
    Args:
        cas: CAS registry number
        registry_dir: Path to reagent registry
        llm_client: Configured LLM client (from llmtools.clients)
        name_override: Optional name override
        role_override: Optional role override (skips LLM classification)
        smiles_override: Optional SMILES override
        resolver_timeout: Timeout for CAS resolution
        
    Returns:
        {
            "status": "ready_to_save" | "needs_review" | "error",
            "workflow": {
                "step1_identity": {...},
                "step2_role": {...},
                "step3_fields": {...},
                "step4_verification": {...}
            },
            "entry": {/* final entry to save */},
            "error": "..." (if status=error)
        }
    """
    if not LLM_SUPPORT_AVAILABLE:
        return {
            "status": "error",
            "error": f"LLM support not available: {LLM_SUPPORT_ERROR}",
        }
    
    if not cas:
        return {"status": "error", "error": "CAS number is required."}
    
    normalized_cas = normalize_cas(cas)
    
    if not registry_dir.exists():
        return {"status": "error", "error": f"Registry directory '{registry_dir}' does not exist."}
    
    workflow = {}
    
    # Step 1: Resolve identity
    try:
        resolved_identity = resolve_identity_from_cas(normalized_cas, timeout=resolver_timeout)
        if not resolved_identity:
            return {
                "status": "error",
                "error": f"Failed to resolve CAS {normalized_cas}. No data from PubChem.",
            }
        
        # Apply overrides
        if name_override:
            resolved_identity["name"] = name_override
        if smiles_override:
            resolved_identity["smiles"] = smiles_override
        
        workflow["step1_identity"] = {
            "status": "success",
            "identity": resolved_identity,
        }
        
    except Exception as exc:
        return {
            "status": "error",
            "error": f"Step 1 (Identity Resolution) failed: {exc}",
        }
    
    # Step 2: Classify role (or use override)
    try:
        llm_suggested_role = None
        
        if role_override:
            # Use the role provided by the user as primary
            role = role_override
            
            # Also get LLM suggestion for comparison/reference
            try:
                llm_result = classify_role(resolved_identity, llm_client)
                if llm_result.get("status") == "success":
                    llm_suggested_role = llm_result.get("role")
                    workflow["step2_role"] = {
                        "status": "success",
                        "role": role,
                        "source": "user_override",
                        "llm_suggestion": llm_suggested_role,
                        "llm_confidence": llm_result.get("confidence"),
                        "llm_reasoning": llm_result.get("reasoning"),
                    }
                else:
                    # LLM suggestion failed, but that's OK since user provided role
                    workflow["step2_role"] = {
                        "status": "success",
                        "role": role,
                        "source": "user_override",
                        "llm_suggestion_error": llm_result.get("error"),
                    }
            except Exception as llm_exc:
                # LLM suggestion failed, but that's OK since user provided role
                workflow["step2_role"] = {
                    "status": "success",
                    "role": role,
                    "source": "user_override",
                    "llm_suggestion_error": str(llm_exc),
                }
        else:
            # LLM classifies the role
            role_result = classify_role(resolved_identity, llm_client)
            workflow["step2_role"] = role_result
            
            if role_result.get("status") != "success":
                return {
                    "status": "error",
                    "error": f"Step 2 (Role Classification) failed: {role_result.get('error')}",
                    "workflow": workflow,
                }
            
            role = role_result["role"]
        
    except Exception as exc:
        return {
            "status": "error",
            "error": f"Step 2 (Role Classification) failed: {exc}",
            "workflow": workflow,
        }
    
    # Step 3: Assign fields
    try:
        fields_result = assign_fields(
            identity=resolved_identity,
            role=role,
            registry_dir=registry_dir,
            llm_client=llm_client,
        )
        workflow["step3_fields"] = fields_result
        
        if fields_result.get("status") != "success":
            return {
                "status": "error",
                "error": f"Step 3 (Field Assignment) failed: {fields_result.get('error')}",
                "workflow": workflow,
            }
        
    except Exception as exc:
        return {
            "status": "error",
            "error": f"Step 3 (Field Assignment) failed: {exc}",
            "workflow": workflow,
        }
    
    # Build entry following the reagent schema
    try:
        # Get family info from registry
        from chemtools.reagent import build_entry, build_embedding_text
        
        # Prepare abbreviations
        abbreviations = fields_result.get("abbreviations", [])
        primary_abbr = abbreviations[0] if abbreviations else ""
        
        # Prepare synonyms
        existing_synonyms = resolved_identity.get("synonyms", [])
        additional_synonyms = fields_result.get("additional_synonyms", [])
        all_synonyms = dedupe_synonyms(existing_synonyms + additional_synonyms)
        
        # Filter aliases to exclude: primary name, CAS, and abbreviations
        name_lower = resolved_identity.get("name", "").lower()
        cas_values = {normalized_cas, normalized_cas.replace("-", "")} if normalized_cas else set()
        abbr_lower = {abbr.lower() for abbr in abbreviations}
        
        aliases = [
            syn for syn in all_synonyms
            if syn.lower() != name_lower  # Not the primary name
            and syn not in cas_values  # Not CAS number
            and syn.lower() not in abbr_lower  # Not an abbreviation
        ]
        
        # Create family dict for embedding text
        family_dict = {
            "family_id": fields_result["family"],
            "label": fields_result["family"],
            **fields_result.get("fields", {})
        }
        
        # Build entry following reagent_schema.json structure
        entry = {
            "id": resolved_identity.get("inchi_key") or normalized_cas,  # prefer InChIKey; else CAS
            "name": resolved_identity.get("name"),
            "abbreviation": abbreviations if abbreviations else [],
            "aliases": aliases,  # Use filtered aliases (no name, CAS, or abbreviations)
            "cas": normalized_cas,
            "inchi_key": resolved_identity.get("inchi_key"),
            "smiles": resolved_identity.get("smiles"),
            "roles": {
                role: {
                    "families": [fields_result["family"]],
                    **fields_result.get("fields", {}),
                }
            },
        }
        
        # If user overrode the role and LLM suggested something different, add it as secondary role
        if llm_suggested_role and llm_suggested_role != role:
            try:
                # Get fields for the LLM-suggested role
                llm_fields_result = assign_fields(
                    identity=resolved_identity,
                    role=llm_suggested_role,
                    registry_dir=registry_dir,
                    llm_client=llm_client,
                )
                
                if llm_fields_result.get("status") == "success":
                    # Add LLM-suggested role as secondary
                    entry["roles"][llm_suggested_role] = {
                        "families": [llm_fields_result["family"]],
                        **llm_fields_result.get("fields", {}),
                        "_note": "LLM-suggested alternative role",
                    }
            except Exception:
                # If getting fields for suggested role fails, that's OK
                # Just use the user's role only
                pass
        
        # Build embedding text following schema
        embedding_entry = {
            "name": entry["name"],
            "abbr": primary_abbr,
            "cas": normalized_cas,
        }
        entry["embedding_text"] = build_embedding_text(role, family_dict, embedding_entry, all_synonyms)
        
        # Remove None values to keep JSON clean, but keep empty lists (required by schema)
        entry = {k: v for k, v in entry.items() if v is not None or isinstance(v, list)}
        
    except Exception as exc:
        return {
            "status": "error",
            "error": f"Failed to build entry: {exc}",
            "workflow": workflow,
        }
    
    # Step 4: Verify entry
    try:
        verification_result = verify_entry(entry, llm_client)
        workflow["step4_verification"] = verification_result
        
        if verification_result.get("status") != "success":
            return {
                "status": "error",
                "error": f"Step 4 (Verification) failed: {verification_result.get('error')}",
                "workflow": workflow,
                "entry": entry,
            }
        
        # Check if entry is approved
        approved = verification_result.get("approved", False)
        issues = verification_result.get("issues", [])
        
        # Count error-level issues
        errors = [iss for iss in issues if iss.get("severity") == "error"]
        
        if not approved or errors:
            return {
                "status": "needs_review",
                "workflow": workflow,
                "entry": entry,
                "message": f"Entry has {len(errors)} error(s) and {len(issues) - len(errors)} warning(s). Please review before saving.",
            }
        
    except Exception as exc:
        return {
            "status": "error",
            "error": f"Step 4 (Verification) failed: {exc}",
            "workflow": workflow,
            "entry": entry,
        }
    
    # Success!
    return {
        "status": "ready_to_save",
        "workflow": workflow,
        "entry": entry,
        "message": "Entry successfully generated and verified by LLM. Ready to save!",
    }


def generate_taxonomy_entry(
    *,
    cas: str,
    role: Optional[str],
    registry_dir: Path,
    allow_default_family: bool,
    dry_run: bool,
    resolver_timeout: float = DEFAULT_RESOLVER_TIMEOUT,
    name_override: Optional[str] = None,
    smiles_override: Optional[str] = None,
    llm_options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Execute the registry workflow and return a JSON-friendly result."""
    if not cas:
        raise RegistryGenerationError("CAS number is required.")
    normalized_cas = normalize_cas(cas)

    if not registry_dir.exists():
        raise RegistryGenerationError(f"Registry directory '{registry_dir}' does not exist.")

    if not role:
        raise RegistryGenerationError("A reagent role must be specified.")
    if role not in ROLE_CONFIG:
        raise RegistryGenerationError(f"Unsupported role '{role}'.")

    store = ReagentRegistryStore(registry_dir)

    resolved_identity = resolve_identity_from_cas(normalized_cas, timeout=resolver_timeout)
    name = name_override or (resolved_identity.get("name") if resolved_identity else None)
    if not name:
        raise RegistryGenerationError(
            "Unable to determine reagent name. Provide a name override or verify the CAS number."
        )

    resolved_synonyms: Sequence[str] = resolved_identity.get("synonyms", []) if resolved_identity else []
    synonyms = dedupe_synonyms([name, *resolved_synonyms])
    tokens = tokenize_all([name, *synonyms])

    family_id: Optional[str] = None
    family_reason: Optional[List[str]] = None
    family_entry: Optional[Dict[str, Any]] = None
    used_default = False
    debug_log: List[str] = []

    inference = store.infer_family(role, normalized_cas, tokens)
    if inference:
        inferred_family, reason_tokens, inferred_entry = inference
        family_id = inferred_family
        family_reason = reason_tokens
        family_entry = inferred_entry
    else:
        default_family = store.default_family(role) if allow_default_family else None
        if default_family:
            family_entry = store.family_entry(role, default_family)
            if family_entry:
                family_id = default_family
                used_default = True
                debug_log.append(
                    f"Used default family '{default_family}' for role '{role}'."
                )
        if not family_id:
            raise RegistryGenerationError(
                "Unable to determine a family for this reagent. Provide a family "
                "override or enable default fallback."
            )

    if not family_entry:
        raise RegistryGenerationError(f"Family metadata not found for '{family_id}' (role '{role}').")

    existing_same_role = store.find_by_cas(normalized_cas, role=role)
    if existing_same_role:
        existing_role, existing_family, data = existing_same_role
        result = {
            "cas": normalized_cas,
            "name": data.get("name"),
            "role": existing_role,
            "family_id": existing_family,
            "status": "exists",
            "registry_file": str(store.file_for_role(existing_role)),
        }
        if family_reason:
            result["family_tokens"] = family_reason
        if debug_log:
            result["debug_log"] = debug_log
        return result

    other_entries = store.cas_index.get(normalized_cas, [])
    existing_other_roles: List[str] = []
    for other_role, _entry in other_entries:
        if other_role != role:
            existing_other_roles.append(str(other_role))
    if existing_other_roles:
        roles_str = ", ".join(sorted(existing_other_roles))
        debug_log.append(
            f"CAS already exists under roles [{roles_str}]; continuing to add entry for '{role}'."
        )

    auto_resolve_source = resolved_identity.get("source") if resolved_identity else None
    resolved_smiles = resolved_identity.get("smiles") if resolved_identity else None
    smiles = smiles_override or resolved_smiles
    inchi_key = None
    if resolved_identity:
        inchi_key = resolved_identity.get("inchi_key") or resolved_identity.get("inchikey")

    abbreviations = list(infer_abbreviations(name, synonyms))
    aliases = build_aliases(name, normalized_cas, abbreviations, synonyms)

    role_payload = store.build_role_payload(role, family_id)
    entry_id = inchi_key or f"cas-{normalized_cas}"
    entry = build_registry_entry(
        entry_id=entry_id,
        name=name,
        abbreviations=abbreviations,
        aliases=aliases,
        cas=normalized_cas,
        smiles=smiles,
        inchi_key=inchi_key,
        role=role,
        role_payload=role_payload,
        family_entry=family_entry,
        synonyms=synonyms,
    )

    result: Dict[str, Any] = {
        "cas": normalized_cas,
        "name": name,
        "role": role,
        "family_id": family_id,
        "registry_file": str(store.file_for_role(role)),
        "dry_run": dry_run,
        "used_default_family": used_default,
    }
    if auto_resolve_source:
        result["auto_resolve_source"] = auto_resolve_source
    if resolved_smiles and not smiles_override:
        result["smiles_source"] = auto_resolve_source or "resolver"
    if family_reason:
        result["family_tokens"] = family_reason
    if debug_log:
        result["debug_log"] = debug_log

    if llm_options:
        if not LLM_SUPPORT_AVAILABLE or not LLMReviewOptions or not review_taxonomy_proposal:
            result["llm_review"] = {
                "enabled": True,
                "status": "error",
                "error": "LLM support is not available in this environment.",
            }
        else:
            # Filter unexpected keys to avoid dataclass errors.
            supported_keys = {"enabled", "provider", "model", "temperature", "max_tokens", "timeout"}
            filtered_opts = {k: v for k, v in llm_options.items() if k in supported_keys}
            options_obj = LLMReviewOptions(**filtered_opts)
            if options_obj.enabled:
                llm_context = {
                    "cas": normalized_cas,
                    "name": name,
                    "role": role,
                    "family_id": family_id,
                    "synonyms": list(synonyms),
                    "abbreviations": abbreviations,
                    "family_entry": family_entry,
                    "family_reason": family_reason or [],
                    "used_default_family": used_default,
                    "debug_log": debug_log,
                    "resolved_identity": resolved_identity or {},
                    "existing_other_roles": existing_other_roles,
                    "registry_file": str(store.file_for_role(role)),
                    "entry_preview": entry,
                }
                llm_review = review_taxonomy_proposal(llm_context, options_obj)
                result["llm_review"] = llm_review
                analysis = llm_review.get("analysis") if isinstance(llm_review, dict) else None
                if isinstance(analysis, dict):
                    adjustment_errors: List[str] = []
                    changes_summary: Dict[str, Any] = {}

                    proposed_role = analysis.get("proposed_role")
                    proposed_family = analysis.get("proposed_family")
                    suggested_synonyms = [
                        str(value).strip()
                        for value in analysis.get("suggested_synonyms", []) or []
                        if str(value).strip()
                    ]

                    combined_synonyms = dedupe_synonyms([*synonyms, *suggested_synonyms])
                    synonyms_added = [syn for syn in combined_synonyms if syn not in synonyms]

                    target_role = role
                    target_family = family_id

                    candidate_role = str(proposed_role).strip() if proposed_role else ""
                    candidate_family = str(proposed_family).strip() if proposed_family else ""

                    # ENHANCEMENT 1: Auto-upgrade "other_reagent" when LLM suggests better role
                    if role == "other_reagent" and candidate_role and candidate_role != role:
                        debug_log.append(
                            f"LLM auto-upgrade: 'other_reagent' 鈫?'{candidate_role}' "
                            f"(confidence: {analysis.get('confidence', 'N/A')})"
                        )
                        result["llm_auto_upgrade"] = {
                            "from": "other_reagent",
                            "to": candidate_role,
                            "reason": analysis.get("justification", "LLM recommendation"),
                            "confidence": analysis.get("confidence"),
                        }

                    if candidate_role and candidate_role != role:
                        if candidate_role in ROLE_CONFIG:
                            # Use suggested family if provided; otherwise keep current until validated.
                            family_candidate_for_role = candidate_family or target_family
                            if family_candidate_for_role and store.family_entry(candidate_role, family_candidate_for_role):
                                target_role = candidate_role
                                target_family = family_candidate_for_role
                                changes_summary["role"] = {"from": role, "to": target_role}
                                if target_family != family_id:
                                    changes_summary["family"] = {"from": family_id, "to": target_family}
                            else:
                                adjustment_errors.append(
                                    f"Suggested role '{candidate_role}' missing valid family '{family_candidate_for_role or 'unspecified'}'."
                                )
                        else:
                            adjustment_errors.append(f"Suggested role '{candidate_role}' is not recognized.")

                    if candidate_family and target_role == role and candidate_family != family_id:
                        if store.family_entry(target_role, candidate_family):
                            target_family = candidate_family
                            changes_summary["family"] = {"from": family_id, "to": target_family}
                        else:
                            adjustment_errors.append(
                                f"Suggested family '{candidate_family}' not found for role '{target_role}'."
                            )

                    should_apply = (
                        target_role != role
                        or target_family != family_id
                        or bool(synonyms_added)
                    )

                    # ENHANCEMENT 2: Apply LLM field suggestions
                    field_suggestions = analysis.get("field_suggestions") or {}
                    if isinstance(field_suggestions, dict) and field_suggestions:
                        changes_summary["field_suggestions_applied"] = field_suggestions
                        should_apply = True

                    if should_apply:
                        if target_role == role and target_family == family_id:
                            target_family_entry = family_entry
                        else:
                            target_family_entry = store.family_entry(target_role, target_family)
                        if not target_family_entry:
                            adjustment_errors.append(
                                f"Family '{target_family}' for role '{target_role}' missing; cannot build adjusted entry."
                            )
                        else:
                            updated_abbreviations = list(infer_abbreviations(name, combined_synonyms))
                            updated_aliases = build_aliases(
                                name,
                                normalized_cas,
                                updated_abbreviations,
                                combined_synonyms,
                            )
                            adjusted_role_payload = store.build_role_payload(target_role, target_family)
                            
                            # ENHANCEMENT 2: Merge LLM field suggestions into role payload
                            if isinstance(field_suggestions, dict):
                                for field_name, field_value in field_suggestions.items():
                                    if field_name and field_value is not None:
                                        adjusted_role_payload[field_name] = field_value
                                        debug_log.append(
                                            f"LLM field suggestion: {field_name} = {field_value}"
                                        )
                            
                            adjusted_entry = build_registry_entry(
                                entry_id=entry_id,
                                name=name,
                                abbreviations=updated_abbreviations,
                                aliases=updated_aliases,
                                cas=normalized_cas,
                                smiles=smiles,
                                inchi_key=inchi_key,
                                role=target_role,
                                role_payload=adjusted_role_payload,
                                family_entry=target_family_entry,
                                synonyms=combined_synonyms,
                            )
                            result["llm_adjusted_entry"] = adjusted_entry
                            if target_role != role:
                                changes_summary.setdefault("role", {"from": role, "to": target_role})
                            if target_family != family_id:
                                changes_summary.setdefault("family", {"from": family_id, "to": target_family})
                            if synonyms_added:
                                changes_summary["synonyms_added"] = synonyms_added
                            if changes_summary:
                                result["llm_applied_changes"] = changes_summary

                    if adjustment_errors:
                        result["llm_adjustment_errors"] = adjustment_errors
            else:
                result["llm_review"] = {"enabled": False, "status": "disabled"}

    if dry_run:
        result["status"] = "dry_run"
        result["entry_preview"] = entry
        return result

    store.add_entry(role, entry)
    path = store.save_role(role)
    result["status"] = "written"
    result["written_to"] = str(path)
    return result



class GenerationWorkerSignals(QObject):
    """Qt signals emitted by the registry generation worker."""

    finished = pyqtSignal(dict)
    error = pyqtSignal(str)


class GenerationWorker(QRunnable):
    """Run registry generation logic on a background thread."""

    def __init__(self, params: Dict[str, Any]) -> None:
        super().__init__()
        self.params = params
        self.signals = GenerationWorkerSignals()

    def run(self) -> None:  # pragma: no cover - UI worker thread
        try:
            workflow_mode = self.params.pop("workflow_mode", "legacy")
            
            if workflow_mode == "llm_workflow":
                # Pure LLM workflow
                cas = self.params["cas"]
                registry_dir = self.params["registry_dir"]
                name_override = self.params.get("name_override")
                role_override = self.params.get("role_override")
                provider = self.params["provider"]
                model = self.params["model"]
                
                # Create LLM client
                from llmtools.clients import LLMClient
                llm_client = LLMClient(provider=provider, model=model)
                
                result = generate_taxonomy_entry_llm(
                    cas=cas,
                    registry_dir=registry_dir,
                    llm_client=llm_client,
                    name_override=name_override,
                    role_override=role_override,
                )
            else:
                # Legacy workflow
                result = generate_taxonomy_entry(**self.params)
                
        except Exception as exc:  # noqa: BLE001 - we want to surface all failures
            self.signals.error.emit(str(exc))
        else:
            self.signals.finished.emit(result)


class RegistryGeneratorWindow(QMainWindow):
    """Main window that wraps the reagent registry workflow."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Reagent Registry Generator")
        self.thread_pool = QThreadPool.globalInstance()
        self._current_worker: Optional[GenerationWorker] = None
        self._last_result: Optional[Dict[str, Any]] = None
        self._llm_support_available = LLM_SUPPORT_AVAILABLE and bool(LLM_AVAILABLE_MODELS)

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        form_layout = QFormLayout()
        layout.addLayout(form_layout)

        self.cas_input = QLineEdit()
        self.cas_input.setPlaceholderText("e.g. 100-00-5")
        form_layout.addRow("CAS number", self.cas_input)

        self.name_input = QLineEdit()
        self.name_input.setPlaceholderText("Optional; overrides resolver result")
        form_layout.addRow("Name override", self.name_input)

        self.role_combo = QComboBox()
        self.role_combo.addItem("Select role", userData=None)
        
        # Add "Auto-detect (LLM)" option if LLM is available
        if self._llm_support_available:
            self.role_combo.addItem("馃 Auto-detect (LLM)", userData="__auto_detect__")
            self.role_combo.setItemData(
                1, 
                "Use LLM to automatically detect the reagent role", 
                Qt.ItemDataRole.ToolTipRole
            )

        role_keys = sorted(
            ROLE_CONFIG.keys(),
            key=lambda key: (ROLE_CONFIG[key].get("priority", 99), key),
        )
        for role_key in role_keys:
            cfg = ROLE_CONFIG[role_key]
            label = cfg.get("label") or role_key.replace("_", " ").title()
            index = self.role_combo.count()
            self.role_combo.addItem(label, userData=role_key)
            hint = cfg.get("hint")
            if hint:
                self.role_combo.setItemData(index, hint, Qt.ItemDataRole.ToolTipRole)
        default_index = self.role_combo.findData("other_reagent")
        if default_index != -1:
            self.role_combo.setCurrentIndex(default_index)
        form_layout.addRow("Reagent role", self.role_combo)
        
        # Workflow mode selection
        workflow_widget = QWidget()
        workflow_layout = QHBoxLayout(workflow_widget)
        workflow_layout.setContentsMargins(0, 0, 0, 0)
        
        self.workflow_mode_combo = QComboBox()
        self.workflow_mode_combo.addItem("Deterministic + LLM Review", userData="legacy")
        self.workflow_mode_combo.addItem("馃殌 Pure LLM Workflow (Recommended)", userData="llm_workflow")
        self.workflow_mode_combo.setCurrentIndex(1)  # Default to pure LLM
        self.workflow_mode_combo.currentIndexChanged.connect(self.on_workflow_mode_changed)
        workflow_layout.addWidget(self.workflow_mode_combo)
        workflow_layout.addStretch()
        form_layout.addRow("Workflow mode", workflow_widget)
        
        llm_widget = QWidget()
        llm_layout = QHBoxLayout(llm_widget)
        llm_layout.setContentsMargins(0, 0, 0, 0)

        self.llm_mode_combo = QComboBox()
        self.llm_mode_combo.addItem("No LLM", userData=False)
        self.llm_mode_combo.addItem("Use LLM", userData=True)
        self.llm_mode_combo.currentIndexChanged.connect(self.on_llm_mode_changed)
        llm_layout.addWidget(self.llm_mode_combo)

        self.llm_provider_combo = QComboBox()
        self.llm_provider_combo.currentIndexChanged.connect(self.on_llm_provider_changed)
        llm_layout.addWidget(self.llm_provider_combo)

        self.llm_model_combo = QComboBox()
        llm_layout.addWidget(self.llm_model_combo)

        llm_layout.addStretch()
        form_layout.addRow("LLM assistance", llm_widget)

        self._populate_llm_provider_options()
        self.on_workflow_mode_changed()  # Initialize workflow mode (will call on_llm_mode_changed)

        path_layout = QHBoxLayout()
        self.registry_path_input = QLineEdit(str(DEFAULT_REGISTRY_DIR))
        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self.on_browse_registry_dir)
        path_layout.addWidget(self.registry_path_input)
        path_layout.addWidget(browse_button)
        form_layout.addRow("Registry directory", path_layout)

        options_layout = QHBoxLayout()
        self.allow_default_checkbox = QCheckBox("Allow default family fallback")
        self.allow_default_checkbox.setChecked(True)
        options_layout.addWidget(self.allow_default_checkbox)
        options_layout.addStretch()
        layout.addLayout(options_layout)

        button_layout = QHBoxLayout()
        self.generate_button = QPushButton("Generate")
        self.generate_button.clicked.connect(self.on_generate_clicked)
        self.save_button = QPushButton("Save")
        self.save_button.setEnabled(False)
        self.save_button.clicked.connect(self.on_save_clicked)
        self.clear_button = QPushButton("Clear")
        self.clear_button.clicked.connect(self.clear_output)
        button_layout.addWidget(self.generate_button)
        button_layout.addWidget(self.save_button)
        button_layout.addWidget(self.clear_button)
        button_layout.addStretch()
        layout.addLayout(button_layout)

        self.status_label = QLabel()
        layout.addWidget(self.status_label)

        self.output_view = QPlainTextEdit()
        self.output_view.setReadOnly(False)
        self.output_view.setPlaceholderText("Generation results will appear here.")
        layout.addWidget(self.output_view, stretch=1)

        self.resize(820, 600)

    def _populate_llm_provider_options(self) -> None:
        """Populate provider and model combos based on available catalog."""
        self.llm_provider_combo.clear()
        self.llm_model_combo.clear()

        if not self._llm_support_available:
            self.llm_provider_combo.addItem("LLM unavailable", userData=None)
            self.llm_provider_combo.setEnabled(False)
            self.llm_model_combo.addItem("LLM unavailable", userData=None)
            self.llm_model_combo.setEnabled(False)
            return

        providers = sorted(LLM_AVAILABLE_MODELS.keys())
        if not providers:
            self.llm_provider_combo.addItem("No providers found", userData=None)
            self.llm_provider_combo.setEnabled(False)
            self.llm_model_combo.addItem("No models", userData=None)
            self.llm_model_combo.setEnabled(False)
            return

        for provider in providers:
            label = provider.replace("_", " ").title()
            self.llm_provider_combo.addItem(label, userData=provider)

        self.llm_provider_combo.setEnabled(True)
        
        # Set default to openai provider
        openai_index = next((i for i in range(self.llm_provider_combo.count()) 
                            if self.llm_provider_combo.itemData(i) == "openai"), 0)
        self.llm_provider_combo.setCurrentIndex(openai_index)
        self._populate_llm_model_options(self.llm_provider_combo.currentData())

    def _populate_llm_model_options(self, provider: Optional[str]) -> None:
        """Populate model combo for selected provider."""
        self.llm_model_combo.clear()

        if not provider:
            self.llm_model_combo.addItem("Select provider", userData=None)
            self.llm_model_combo.setEnabled(False)
            return

        models = LLM_AVAILABLE_MODELS.get(provider, [])
        if not models:
            self.llm_model_combo.addItem("No models", userData=None)
            self.llm_model_combo.setEnabled(False)
            return

        recommended = self._recommended_model_for_provider(provider)
        recommended_index = 0
        for idx, model in enumerate(models):
            display = model
            if recommended and model == recommended:
                display = f"{model} (recommended)"
                recommended_index = idx
            self.llm_model_combo.addItem(display, userData=model)

        self.llm_model_combo.setEnabled(True)
        self.llm_model_combo.setCurrentIndex(recommended_index)

    def _recommended_model_for_provider(self, provider: str) -> Optional[str]:
        """Return preferred model for provider using catalog heuristics."""
        preferences = LLM_RECOMMENDED_MODELS.get(provider, {})
        for key in ("balanced", "fast", "reasoning", "advanced"):
            candidate = preferences.get(key)
            if candidate:
                return candidate
        models = LLM_AVAILABLE_MODELS.get(provider, [])
        return models[0] if models else None

    def on_workflow_mode_changed(self, index: int = -1) -> None:  # noqa: ARG002
        """Handle workflow mode changes."""
        workflow_mode = self.workflow_mode_combo.currentData()
        
        # Pure LLM workflow requires LLM to be enabled
        if workflow_mode == "llm_workflow":
            # Auto-enable LLM and show provider/model options
            self.llm_mode_combo.blockSignals(True)
            self.llm_mode_combo.setCurrentIndex(1)  # "Use LLM"
            self.llm_mode_combo.blockSignals(False)
            self.llm_mode_combo.setEnabled(False)  # Can't disable LLM in pure mode
            self.on_llm_mode_changed()
            
            # Show role auto-detect option
            auto_detect_idx = self.role_combo.findData("__auto_detect__")
            if auto_detect_idx != -1:
                self.role_combo.setCurrentIndex(auto_detect_idx)
            
        else:  # legacy mode
            # Re-enable LLM toggle
            self.llm_mode_combo.setEnabled(True)
            
            # Set role back to other_reagent if it was auto-detect
            if self.role_combo.currentData() == "__auto_detect__":
                default_index = self.role_combo.findData("other_reagent")
                if default_index != -1:
                    self.role_combo.setCurrentIndex(default_index)

    def on_llm_mode_changed(self, index: int = -1) -> None:  # noqa: ARG002 - index unused
        """Toggle provider/model controls based on LLM mode selection."""
        # When the widget itself is disabled (e.g. worker running), skip updates.
        if not self.llm_mode_combo.isEnabled():
            return

        enabled = bool(self.llm_mode_combo.currentData())

        if enabled and not self._llm_support_available:
            QMessageBox.information(
                self,
                "LLM support unavailable",
                "LLM integration is not available in this environment.\n"
                "Install optional dependencies and configure API keys to enable it."
                f"\n\nDetected issue: {LLM_SUPPORT_ERROR or 'Unknown'}",
            )
            self.llm_mode_combo.blockSignals(True)
            self.llm_mode_combo.setCurrentIndex(0)
            self.llm_mode_combo.blockSignals(False)
            enabled = False

        provider_has_options = self.llm_provider_combo.count() > 0 and self.llm_provider_combo.itemData(0) is not None
        provider_ready = enabled and provider_has_options and self._llm_support_available
        self.llm_provider_combo.setEnabled(provider_ready)
        if provider_ready:
            self._populate_llm_model_options(self.llm_provider_combo.currentData())
        else:
            self.llm_provider_combo.setEnabled(False)

        model_has_options = self.llm_model_combo.count() > 0 and self.llm_model_combo.itemData(0) is not None
        model_ready = enabled and model_has_options and self._llm_support_available
        self.llm_model_combo.setEnabled(model_ready)

    def on_llm_provider_changed(self, index: int = -1) -> None:  # noqa: ARG002 - index unused
        """Refresh model options when the provider selection changes."""
        if not self._llm_support_available:
            return
        if not self.llm_provider_combo.isEnabled():
            return
        provider = self.llm_provider_combo.currentData()
        self._populate_llm_model_options(provider)

    def on_browse_registry_dir(self) -> None:
        """Let the user choose a registry directory."""
        directory = QFileDialog.getExistingDirectory(self, "Select registry directory", str(DEFAULT_REGISTRY_DIR))
        if directory:
            self.registry_path_input.setText(directory)

    def on_generate_clicked(self) -> None:
        """Kick off registry generation with the current form values."""
        cas = self.cas_input.text().strip()
        if not cas:
            self.show_error("CAS number is required.")
            return

        workflow_mode = self.workflow_mode_combo.currentData()
        role = self.role_combo.currentData()
        
        # For pure LLM workflow, role can be auto-detect
        if workflow_mode == "llm_workflow":
            if not role or role == "__auto_detect__":
                # Auto-detect is allowed in pure LLM mode
                pass
            # Role is optional in pure LLM mode
        else:
            # Legacy mode requires explicit role
            if not role or role == "__auto_detect__":
                self.show_error("Select a reagent role before generating (auto-detect only available in Pure LLM mode).")
                return

        registry_dir_text = self.registry_path_input.text().strip()
        if not registry_dir_text:
            self.show_error("Registry directory is required.")
            return

        registry_dir = Path(registry_dir_text).expanduser()
        name_override = self.name_input.text().strip() or None
        
        # Get LLM settings
        provider = self.llm_provider_combo.currentData()
        model = self.llm_model_combo.currentData()
        
        if workflow_mode == "llm_workflow":
            # Pure LLM workflow - validate LLM is configured
            if not provider or not model:
                self.show_error("Pure LLM workflow requires an LLM provider and model.")
                return
            
            self.status_label.setText("Running Pure LLM workflow...")
            self.output_view.clear()
            self.save_button.setEnabled(False)
            self._last_result = None
            self.set_inputs_enabled(False)
            
            # Pass role if user selected one (skip if auto-detect)
            role_to_pass = role if role and role != "__auto_detect__" else None
            
            params = {
                "workflow_mode": "llm_workflow",
                "cas": cas,
                "registry_dir": registry_dir,
                "name_override": name_override,
                "role_override": role_to_pass,
                "provider": provider,
                "model": model,
            }
            
        else:
            # Legacy workflow
            allow_default = self.allow_default_checkbox.isChecked()
            llm_options: Optional[Dict[str, Any]] = None
            
            if bool(self.llm_mode_combo.currentData()):
                if not provider or not model:
                    self.show_error("Select an LLM provider and model before generating.")
                    return
                llm_options = {
                    "enabled": True,
                    "provider": provider,
                    "model": model,
                    "temperature": 0.0,
                    "max_tokens": 800,
                    "timeout": 60,
                }

            self.status_label.setText("Running registry generation...")
            self.output_view.clear()
            self.save_button.setEnabled(False)
            self._last_result = None
            self.set_inputs_enabled(False)

            params = {
                "workflow_mode": "legacy",
                "cas": cas,
                "role": role,
                "registry_dir": registry_dir,
                "allow_default_family": allow_default,
                "dry_run": True,
                "resolver_timeout": DEFAULT_RESOLVER_TIMEOUT,
                "name_override": name_override,
            }
            if llm_options:
                params["llm_options"] = llm_options
        
        worker = GenerationWorker(params)
        worker.signals.finished.connect(self.on_generation_success)
        worker.signals.error.connect(self.on_generation_failure)
        self._current_worker = worker
        self.thread_pool.start(worker)

    def on_generation_success(self, result: Dict[str, Any]) -> None:
        """Handle successful completion of the registry worker."""
        self._last_result = result
        
        # Detect workflow mode from result structure
        is_pure_llm = "workflow" in result and "status" in result and "entry" in result
        
        if is_pure_llm:
            # Pure LLM workflow output
            status = result.get("status", "unknown")
            entry = result.get("entry", {})
            message = result.get("message", "")
            
            # Update status label with helpful message
            if status == "ready_to_save":
                self.status_label.setText("鉁?LLM Approved - Ready to Save")
            elif status == "needs_review":
                self.status_label.setText("鈿狅笍 Needs Review - Check entry before saving")
            else:
                self.status_label.setText(f"Status: {status}")
            
            # Show ONLY the entry that will be saved (simplified for user review)
            display_payload: Dict[str, Any] = {
                "status": status,
            }
            
            # Add message if present
            if message:
                display_payload["message"] = message
            
            # Show the entry - this is what will be saved
            display_payload["entry"] = entry
            
            # Add verification warnings if present (from needs_review status)
            workflow = result.get("workflow", {})
            verification = workflow.get("step4_verification", {})
            if isinstance(verification, dict):
                issues = verification.get("issues", [])
                if issues:
                    # Show issues for user review
                    display_payload["llm_verification_issues"] = issues
            
            self.output_view.setPlainText(json.dumps(display_payload, indent=2, ensure_ascii=False))
            self.set_inputs_enabled(True)
            
            # Enable save if ready or needs review (both have valid entries)
            self.save_button.setEnabled(status in ("ready_to_save", "needs_review"))
            
        else:
            # Legacy workflow output
            self.status_label.setText(f"Status: {result.get('status', 'ok')}")
            entry_preview = result.get('entry_preview')
            llm_review = result.get('llm_review')
            llm_adjusted_entry = result.get('llm_adjusted_entry')
            llm_auto_upgrade = result.get('llm_auto_upgrade')
            
            display_payload: Dict[str, Any] = {}
            
            # ENHANCEMENT 3: Enhanced output format for review
            if isinstance(llm_review, dict) and isinstance(llm_adjusted_entry, dict):
                # Show comprehensive LLM review output
                display_payload["review_summary"] = {
                    "original_role": result.get("role"),
                    "original_family": result.get("family_id"),
                    "llm_status": llm_review.get("analysis", {}).get("status") if isinstance(llm_review.get("analysis"), dict) else "unknown",
                    "confidence": llm_review.get("analysis", {}).get("confidence") if isinstance(llm_review.get("analysis"), dict) else None,
                    "justification": llm_review.get("analysis", {}).get("justification") if isinstance(llm_review.get("analysis"), dict) else None,
                }
                
                if llm_auto_upgrade:
                    display_payload["review_summary"]["auto_upgrade"] = llm_auto_upgrade
                
                if isinstance(entry_preview, dict):
                    display_payload["entry_original"] = entry_preview
                
                display_payload["entry_revised"] = llm_adjusted_entry
                display_payload["llm_review_details"] = {
                    "model": llm_review.get("model"),
                    "provider": llm_review.get("provider"),
                    "tokens_used": llm_review.get("total_tokens"),
                    "latency_ms": llm_review.get("latency_ms"),
                }
                
                if isinstance(llm_review.get("analysis"), dict):
                    alerts = llm_review["analysis"].get("alerts", [])
                    if alerts:
                        display_payload["llm_alerts"] = alerts
                
                if result.get("llm_applied_changes"):
                    display_payload["changes_applied"] = result["llm_applied_changes"]
                
                if result.get("llm_adjustment_errors"):
                    display_payload["adjustment_errors"] = result["llm_adjustment_errors"]
            else:
                # No LLM or simple output
                if isinstance(entry_preview, dict):
                    display_payload["entry"] = entry_preview
                elif isinstance(llm_adjusted_entry, dict):
                    display_payload["entry"] = llm_adjusted_entry
                else:
                    exclude_keys = {'dry_run', 'registry_file'}
                    display_payload = {k: v for k, v in result.items() if k not in exclude_keys}
                if llm_review:
                    display_payload["llm_review"] = llm_review
            
            self.output_view.setPlainText(json.dumps(display_payload, indent=2, ensure_ascii=False))
            self.set_inputs_enabled(True)
            has_entry = any(
                key in display_payload for key in ("entry_revised", "entry", "entry_original")
            )
            self.save_button.setEnabled(has_entry)
        
        self._current_worker = None

    def on_save_clicked(self) -> None:
        """Persist the edited generation result to disk."""
        payload_text = self.output_view.toPlainText().strip()
        if not payload_text:
            self.show_error("Generate an entry before saving.")
            return
        try:
            payload = json.loads(payload_text)
        except json.JSONDecodeError as exc:
            self.show_error(f"Output is not valid JSON: {exc}")
            return
        entry_to_save = self._extract_entry_from_payload(payload)
        if not entry_to_save:
            self.show_error(
                "Output must contain a registry entry. Provide a JSON object with fields like 'cas', 'name', "
                "and 'roles', or include it under 'entry_revised', 'llm_adjusted_entry', 'entry', or 'entry_preview'."
            )
            return
        registry_dir_text = self.registry_path_input.text().strip()
        if not registry_dir_text:
            self.show_error("Registry directory is required.")
            return
        registry_dir = Path(registry_dir_text).expanduser()
        try:
            store = ReagentRegistryStore(registry_dir)
        except Exception as exc:
            self.show_error(str(exc))
            return
        roles_payload = entry_to_save.get("roles")
        if not isinstance(roles_payload, dict) or not roles_payload:
            self.show_error("Entry must include a 'roles' mapping to determine destination file.")
            return
        role_for_save = canonical_role(next(iter(roles_payload.keys())))
        role_details = roles_payload.get(role_for_save) or {}
        families = role_details.get("families")
        family_for_save: Optional[str] = families[0] if isinstance(families, list) and families else None
        result_context = self._last_result or {}
        try:
            store.add_entry(role_for_save, entry_to_save)
            path = store.save_role(role_for_save)
        except Exception as exc:
            self.show_error(str(exc))
            return
        result_context['status'] = 'written'
        result_context['written_to'] = str(path)
        result_context['entry_preview'] = entry_to_save
        result_context['role'] = role_for_save
        if family_for_save:
            result_context['family_id'] = family_for_save
        if isinstance(result_context.get("llm_adjusted_entry"), dict):
            result_context['llm_adjusted_entry'] = entry_to_save
        self.output_view.setPlainText(json.dumps(entry_to_save, indent=2, ensure_ascii=False))
        self.status_label.setText("Status: written")
        self.save_button.setEnabled(False)
        self._last_result = result_context

    def on_generation_failure(self, message: str) -> None:
        """Display worker errors to the user."""
        self.set_inputs_enabled(True)
        self.status_label.setText("Generation failed.")
        self.save_button.setEnabled(False)
        QMessageBox.critical(self, "Reagent registry generation failed", message)
        self._current_worker = None

    def clear_output(self) -> None:
        """Clear the output preview."""
        self.output_view.clear()
        self.status_label.clear()
        self.save_button.setEnabled(False)
        self._last_result = None

    def set_inputs_enabled(self, enabled: bool) -> None:
        """Toggle input widgets while the worker is running."""
        for widget in (
            self.cas_input,
            self.name_input,
            self.role_combo,
            self.llm_mode_combo,
            self.llm_provider_combo,
            self.llm_model_combo,
            self.registry_path_input,
            self.allow_default_checkbox,
            self.generate_button,
            self.save_button,
            self.clear_button,
        ):
            widget.setEnabled(enabled)
        if enabled:
            self.on_llm_mode_changed()

    @staticmethod
    def _is_registry_entry(candidate: Any) -> bool:
        """Return True if candidate looks like a registry entry payload."""
        return (
            isinstance(candidate, dict)
            and "cas" in candidate
            and "name" in candidate
            and isinstance(candidate.get("roles"), dict)
            and candidate.get("roles")
        )

    def _extract_entry_from_payload(self, payload: Any) -> Optional[Dict[str, Any]]:
        """Extract the first plausible registry entry from a payload structure."""
        if self._is_registry_entry(payload):
            return payload
        if isinstance(payload, dict):
            for key in (
                "entry_revised",
                "llm_adjusted_entry",
                "candidate_entry",
                "entry",
                "entry_original",
                "entry_preview",
            ):
                value = payload.get(key)
                if self._is_registry_entry(value):
                    return value
        return None

    def show_error(self, message: str) -> None:
        """Show a modal error dialog."""
        QMessageBox.warning(self, "Input error", message)


def main() -> None:
    """Launch the PyQt6 application."""
    app = QApplication(sys.argv)
    window = RegistryGeneratorWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

