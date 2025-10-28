"""
Lightweight registry store for the flattened reagent database (data/reagent_db).

This mirrors the logic used by the reagent addition GUI but without any GUI
dependencies so it can be reused by backend services and LangChain tools.
"""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import re

from .taxonomy_store import ROLE_KEYWORDS_RAW, load_families_registry_entries
from .taxonomy_utils import tokenize_all
from .constants import DEFAULT_FAMILY_BY_ROLE

# -----------------------------------------------------------------------------
# Role configuration for flattened registry schema
# -----------------------------------------------------------------------------

ROLE_CONFIG: Dict[str, Dict[str, Any]] = {
    "ligand": {
        "filename": "ligand.json",
        "default_family": "generic_ligands",
    },
    "metal_precursor": {
        "filename": "metal_precursor.json",
        "default_family": "pd_ii_salts",
    },
    "preformed_metal_catalyst": {
        "filename": "preformed_metal_catalyst.json",
        "default_family": "pd_nhc_precatalysts",
    },
    "base": {
        "filename": "base.json",
        "default_family": "carbonates",
    },
    "acid": {
        "filename": "acid.json",
        "default_family": "triflic_acids",
    },
    "additive": {
        "filename": "additive.json",
        "default_family": "phase_transfer_catalysts",
    },
    "condensation_agent": {
        "filename": "condensation_agent.json",
        "default_family": "organophosphorus_anhydrides",
    },
    "oxidant": {
        "filename": "oxidant.json",
        "default_family": "peroxides",
    },
    "reductant": {
        "filename": "reductant.json",
        "default_family": "hydride_reductants",
    },
    "solvent": {
        "filename": "solvent.json",
        "default_family": "aprotic_polar",
    },
    "other_reagent": {
        "filename": "other_reagent.json",
        "default_family": "misc_reagents",
    },
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
    "acid": ("acidity",),
    "solvent": ("polarity_band", "proticity"),
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
    "solvent": ("polarity_band", "proticity"),
    "other_reagent": (),
}

# -----------------------------------------------------------------------------
# Helper utilities (ported from the GUI module without PyQt dependencies)
# -----------------------------------------------------------------------------

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
        if not clean or _has_long_digit_run(clean):
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
        if not clean or _has_long_digit_run(clean):
            continue
        key = clean.lower()
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
    family_label = family_entry.get("definition") or family_entry.get("label") or family_entry.get("family") or ""
    family_id = family_entry.get("family", "")
    if family_label:
        pieces.append(f"family: {family_label} ({family_id})")
    else:
        pieces.append(f"family: {family_id}")
    role_payload = entry.get("roles", {}).get(role, {})
    for field in ROLE_EMBED_FIELDS.get(role, ()):
        value = role_payload.get(field)
        if value is not None:
            pieces.append(f"{ROLE_FIELD_LABELS.get(field, field)}: {_format_embedding_value(value)}")
    for field, value in role_payload.items():
        if field == "families" or field in ROLE_EMBED_FIELDS.get(role, ()):
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


# -----------------------------------------------------------------------------
# Registry store
# -----------------------------------------------------------------------------

class RegistryStore:
    """Minimal store for the flattened reagent registry."""

    def __init__(self, base_dir: Optional[Path] = None) -> None:
        if base_dir is None:
            base_dir = Path("data/reagent_db")
        self.base_dir = Path(base_dir)
        if not self.base_dir.exists():
            raise FileNotFoundError(f"Registry directory '{self.base_dir}' does not exist.")

        self.role_entries: Dict[str, List[Dict[str, Any]]] = {}
        self.family_lookup: Dict[str, Tuple[str, Dict[str, Any]]] = {}
        self.family_tokens: Dict[Tuple[str, str], Set[str]] = {}
        self.family_examples: Dict[Tuple[str, str], Set[str]] = {}
        self.family_numeric_baseline: Dict[Tuple[str, str], Optional[Dict[str, Any]]] = {}
        self.cas_index: Dict[str, Tuple[str, str, Dict[str, Any]]] = {}

        self._load_families()
        self._load_registry()

    # ------------------------------------------------------------------ loading
    def _load_families(self) -> None:
        schema_path = self.base_dir / "reagent_schema" / "families_registry.json"
        entries: List[Dict[str, Any]] = []
        if schema_path.exists():
            data = json.loads(schema_path.read_text(encoding="utf-8"))
            entries = list(data.get("entries", []))
        else:
            try:
                entries = list(load_families_registry_entries())
            except FileNotFoundError:
                entries = []
        if not entries:
            def add_entry(role: str, family: str, definition: str, keywords: Sequence[str]) -> None:
                entries.append({
                    "role": role,
                    "family": family,
                    "label": definition,
                    "definition": definition,
                    "notes": "",
                    "keywords": list(keywords),
                    "required_props": {},
                    "examples_pos": [],
                })

            for role, cfg in ROLE_CONFIG.items():
                default_family = cfg.get("default_family") or f"{role}_general"
                add_entry(role, default_family, f"Auto-generated family for {role}", [])

            # Additional common ligand scaffolds for better suggestions.
            add_entry(
                "ligand",
                "trialkyl_triaryl_phosphines",
                "Auto family for phosphine ligands",
                ["phosphine", "phosphorus", "pr3", "ph3"],
            )
            add_entry(
                "ligand",
                "bipyridine_ligands",
                "Auto family for bipyridine ligands",
                ["bipyridine", "bipy", "bpy", "diimine"],
            )
            add_entry(
                "ligand",
                "phenanthroline_ligands",
                "Auto family for phenanthroline ligands",
                ["phenanthroline", "phen", "phenanthroline"],
            )
        for entry in entries:
            role = entry.get("role")
            family = entry.get("family")
            if not role or not family:
                continue
            self.family_lookup[family] = (role, entry)
            key = (role, family)
            tokens: Set[str] = set()
            tokens.update(tokenize_all(entry.get("keywords", [])))
            tokens.update(tokenize_all([family, entry.get("definition", ""), entry.get("notes", "")]))
            self.family_tokens[key] = {tok for tok in tokens if tok}
            examples = {str(cas) for cas in entry.get("examples_pos", []) if cas}
            self.family_examples[key] = examples
    def _load_registry(self) -> None:
        for role, cfg in ROLE_CONFIG.items():
            filename = cfg.get("filename")
            if not filename:
                continue
            path = self.base_dir / filename
            if not path.exists():
                self.role_entries[role] = []
                continue
            content = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(content, list):
                raise ValueError(f"Registry file '{path}' must contain a JSON list.")
            self.role_entries[role] = content
            for entry in content:
                cas = entry.get("cas")
                if not cas:
                    continue
                families = entry.get("roles", {}).get(role, {}).get("families", [])
                family = families[0] if families else None
                self.cas_index[str(cas)] = (role, family, entry)
                if family and self.family_numeric_baseline.get((role, family)) is None:
                    numeric = entry.get("numeric_features")
                    if numeric:
                        self.family_numeric_baseline[(role, family)] = deepcopy(numeric)

    def build_role_payload(self, role: str, family_id: str) -> Dict[str, Any]:
        payload: Dict[str, Any] = {"families": [family_id]}
        family = self.family_entry(role, family_id)
        if family:
            required_props = family.get("required_props") or {}
            for key, value in required_props.items():
                payload[key] = deepcopy(value)
        defaults = self._default_role_fields(role, family_id)
        for field in ROLE_PAYLOAD_FIELDS.get(role, ()):
            if field not in payload:
                payload[field] = deepcopy(defaults.get(field))
        return payload

    def _default_role_fields(self, role: str, family_id: str) -> Dict[str, Any]:
        entries = self.role_entries.get(role, [])
        for entry in entries:
            role_payload = entry.get("roles", {}).get(role, {})
            families = role_payload.get("families") or []
            if family_id not in families:
                continue
            defaults: Dict[str, Any] = {}
            for field in ROLE_PAYLOAD_FIELDS.get(role, ()):
                if field in role_payload:
                    defaults[field] = deepcopy(role_payload[field])
            if defaults:
                return defaults
        return {}

    # ------------------------------------------------------------------ queries
    def file_for_role(self, role: str) -> Path:
        cfg = ROLE_CONFIG.get(role)
        if not cfg or not cfg.get("filename"):
            raise KeyError(f"Unknown or unsupported role '{role}'.")
        return self.base_dir / cfg["filename"]

    def role_for_family(self, family_id: str) -> Optional[str]:
        data = self.family_lookup.get(family_id)
        return data[0] if data else None

    def family_data(self, role: str, family_id: str) -> Dict[str, Any]:
        data = self.family_lookup.get(family_id)
        if not data or data[0] != role:
            raise KeyError(f"Family '{family_id}' not found for role '{role}'.")
        return data[1]

    def numeric_baseline(self, role: str, family_id: str) -> Optional[Dict[str, Any]]:
        return self.family_numeric_baseline.get((role, family_id))

    def family_entry(self, role: str, family_id: str) -> Optional[Dict[str, Any]]:
        data = self.family_lookup.get(family_id)
        if not data or data[0] != role:
            return None
        return data[1]

    def default_family(self, role: str) -> Optional[str]:
        cfg = ROLE_CONFIG.get(role)
        return cfg.get("default_family") if cfg else DEFAULT_FAMILY_BY_ROLE.get(role)

    def list_families(self) -> List[Dict[str, Any]]:
        items: List[Dict[str, Any]] = []
        for family_id, (role, data) in sorted(self.family_lookup.items()):
            items.append(
                {
                    "role": role,
                    "family_id": family_id,
                    "label": data.get("label") or data.get("definition") or "",
                }
            )
        return items

    def family_token_overlap(self, role: str, family_id: str, tokens: Iterable[str]) -> bool:
        familial = self.family_tokens.get((role, family_id), set())
        return bool(familial & set(tokens))

    def find_by_cas(self, cas: str) -> Optional[Tuple[str, str, Dict[str, Any]]]:
        record = self.cas_index.get(cas)
        if record:
            role, family, entry = record
            return role, family, entry
        return None

    # ------------------------------------------------------------------ updates
    def add_entry(self, role: str, entry: Dict[str, Any]) -> None:
        entries = self.role_entries.setdefault(role, [])
        entries.append(entry)
        entries.sort(key=lambda item: (item.get("name") or "").lower())
        cas = entry.get("cas")
        families = entry.get("roles", {}).get(role, {}).get("families", [])
        family = families[0] if families else None
        if cas:
            self.cas_index[str(cas)] = (role, family, entry)
        tokens = tokenize_all([entry.get("name"), *entry.get("abbreviation", []), *entry.get("aliases", [])])
        if family:
            self.family_tokens.setdefault((role, family), set()).update(tokens)

    def save_role(self, role: str) -> Path:
        path = self.file_for_role(role)
        path.parent.mkdir(parents=True, exist_ok=True)
        entries = self.role_entries.get(role, [])
        path.write_text(json.dumps(entries, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return path

    # ------------------------------------------------------------------ heuristics
    def infer_family(self, role: str, tokens: Set[str]) -> Optional[Tuple[str, List[str]]]:
        best: Optional[Tuple[int, List[str], str]] = None
        for (fam_role, fam_id), fam_tokens in self.family_tokens.items():
            if fam_role != role:
                continue
            matches = sorted(tokens & fam_tokens)
            if not matches:
                continue
            score = len(matches)
            if not best or score > best[0]:
                best = (score, matches, fam_id)
        if best:
            _score, matches, fam_id = best
            return fam_id, matches
        return None

    def suggest_families(self, role: str, tokens: Iterable[str], limit: int = 5) -> List[Dict[str, Any]]:
        token_set = {tok for tok in tokens if tok}
        suggestions: List[Tuple[int, str, List[str]]] = []
        for (fam_role, fam_id), fam_tokens in self.family_tokens.items():
            if fam_role != role:
                continue
            matches = sorted(token_set & fam_tokens)
            if not matches:
                continue
            suggestions.append((len(matches), fam_id, matches))
        suggestions.sort(key=lambda item: (item[0], item[1]), reverse=True)

        results: List[Dict[str, Any]] = []
        for score, fam_id, matches in suggestions[:limit]:
            entry = self.family_lookup.get(fam_id)
            label = ""
            if entry:
                label = entry[1].get("label") or entry[1].get("definition") or fam_id
            results.append({
                "family_id": fam_id,
                "label": label,
                "matches": matches,
                "score": score,
            })
        return results

    def infer_role(self, name: str, synonyms: Sequence[str]) -> Optional[Tuple[str, str]]:
        texts = " ".join([name, *synonyms])
        for role, patterns in ROLE_KEYWORDS_RAW.items():
            for pattern in patterns:
                if re.search(pattern, texts, flags=re.IGNORECASE):
                    return role, pattern
        return None
