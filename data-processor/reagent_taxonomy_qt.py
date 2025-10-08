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
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

from reagent_taxonomy_generator import (  # type: ignore
    DEFAULT_RESOLVER_TIMEOUT,
    dedupe_synonyms,
    normalize_cas,
    resolve_identity_from_cas,
    tokenize_all,
)


DEFAULT_REGISTRY_DIR = (MODULE_DIR.parent / "data" / "reagents").resolve()

ROLE_CONFIG: Dict[str, Dict[str, Any]] = {
    "ligand": {
        "filename": "ligand.json",
        "label": "Ligand",
        "hint": "Ligands including phosphines, NHCs, diimines, and ancillary donor sets.",
        "priority": 0,
        "default_family": "trialkyl_triaryl_phosphines",
    },
    "metal_precursor": {
        "filename": "metal_precursor.json",
        "label": "Metal precursor",
        "hint": "Metal salts or complexes that generate the catalytically active species.",
        "priority": 1,
        "default_family": "pd_ii_salts",
    },
    "preformed_metal_catalyst": {
        "filename": "preformed_metal_catalyst.json",
        "label": "Preformed metal catalyst",
        "hint": "Precatalysts supplied with ligands; typically used as-is.",
        "priority": 2,
        "default_family": "pd_phosphine_complexes",
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
    "metal_precursor": ("metal", "oxidation_states"),
    "preformed_metal_catalyst": ("metal", "oxidation_states", "ligand_type"),
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
    "metal_precursor": ("metal", "oxidation_states"),
    "oxidant": ("strength_band",),
    "preformed_metal_catalyst": ("metal", "oxidation_states", "ligand_type"),
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
        schema_path = self.base_dir / "reagent_schema" / "families_registry.json"
        if not schema_path.exists():
            raise FileNotFoundError(f"Families registry not found: {schema_path}")
        data = json.loads(schema_path.read_text(encoding="utf-8"))
        for entry in data.get("entries", []):
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
        "abbreviation": list(abbreviations),
        "aliases": list(aliases),
        "cas": cas,
        "inchi_key": inchi_key,
        "smiles": smiles,
        "roles": {role: role_payload},
    }
    entry["embedding_text"] = build_embedding_text(role, family_entry, entry, synonyms)
    return entry



class RegistryGenerationError(RuntimeError):
    """Raised when registry generation fails for the requested input."""



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
    for other_role, _entry in other_entries:
        if other_role != role:
            debug_log.append(
                f"CAS already exists under role '{other_role}'; continuing to add entry for '{role}'."
            )
            break

    auto_resolve_source = resolved_identity.get("source") if resolved_identity else None
    resolved_smiles = resolved_identity.get("smiles") if resolved_identity else None
    smiles = smiles_override or resolved_smiles
    inchi_key = None
    if resolved_identity:
        inchi_key = resolved_identity.get("inchi_key") or resolved_identity.get("inchikey")

    abbreviations = infer_abbreviations(name, synonyms)
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

        role = self.role_combo.currentData()
        if not role:
            self.show_error("Select a reagent role before generating.")
            return

        registry_dir_text = self.registry_path_input.text().strip()
        if not registry_dir_text:
            self.show_error("Registry directory is required.")
            return

        registry_dir = Path(registry_dir_text).expanduser()
        allow_default = self.allow_default_checkbox.isChecked()
        name_override = self.name_input.text().strip() or None

        self.status_label.setText("Running registry generation...")
        self.output_view.clear()
        self.save_button.setEnabled(False)
        self._last_result = None
        self.set_inputs_enabled(False)

        params = {
            "cas": cas,
            "role": role,
            "registry_dir": registry_dir,
            "allow_default_family": allow_default,
            "dry_run": True,
            "resolver_timeout": DEFAULT_RESOLVER_TIMEOUT,
            "name_override": name_override,
        }
        worker = GenerationWorker(params)
        worker.signals.finished.connect(self.on_generation_success)
        worker.signals.error.connect(self.on_generation_failure)
        self._current_worker = worker
        self.thread_pool.start(worker)

    def on_generation_success(self, result: Dict[str, Any]) -> None:
        """Handle successful completion of the registry worker."""
        self._last_result = result
        self.status_label.setText(f"Status: {result.get('status', 'ok')}")
        entry_preview = result.get('entry_preview')
        if isinstance(entry_preview, dict):
            display_payload = entry_preview
        else:
            display_payload = {k: v for k, v in result.items() if k not in {'dry_run', 'registry_file'}}
        self.output_view.setPlainText(json.dumps(display_payload, indent=2, ensure_ascii=False))
        self.set_inputs_enabled(True)
        self.save_button.setEnabled(bool(entry_preview))
        self._current_worker = None

    def on_save_clicked(self) -> None:
        """Persist the edited generation result to disk."""
        payload_text = self.output_view.toPlainText().strip()
        if not payload_text:
            self.show_error("Generate an entry before saving.")
            return
        try:
            entry_preview = json.loads(payload_text)
        except json.JSONDecodeError as exc:
            self.show_error(f"Output is not valid JSON: {exc}")
            return
        if not isinstance(entry_preview, dict):
            self.show_error("Edited output must be a JSON object representing the registry entry.")
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
        result_context = self._last_result or {}
        role = result_context.get('role')
        if not role:
            self.show_error("Role missing from result; cannot determine destination file.")
            return
        try:
            store.add_entry(role, entry_preview)
            path = store.save_role(role)
        except Exception as exc:
            self.show_error(str(exc))
            return
        result_context['status'] = 'written'
        result_context['written_to'] = str(path)
        result_context['entry_preview'] = entry_preview
        self.output_view.setPlainText(json.dumps(entry_preview, indent=2, ensure_ascii=False))
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
            self.registry_path_input,
            self.allow_default_checkbox,
            self.generate_button,
            self.save_button,
            self.clear_button,
        ):
            widget.setEnabled(enabled)

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
