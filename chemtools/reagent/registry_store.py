"""
Lightweight registry store for the flattened reagent database (data/reagent_db/reagents.csv).

This mirrors the logic used by the reagent addition GUI but without any GUI
dependencies so it can be reused by backend services and LangChain tools.
"""

from __future__ import annotations

import csv
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import re

from .taxonomy_store import ROLE_KEYWORDS_RAW, load_families_registry_entries
from .taxonomy_utils import tokenize_all
from .constants import DEFAULT_FAMILY_BY_ROLE, ROLE_ALIASES, ROLE_PRIORITY



def _canonical_role(role: str) -> str:
    """Normalize legacy role names to the canonical registry labels."""
    return ROLE_ALIASES.get(role, role)


# -----------------------------------------------------------------------------
# Role configuration for flattened registry schema
# -----------------------------------------------------------------------------

ROLE_CONFIG: Dict[str, Dict[str, Any]] = {
    "ligand": {"default_family": "trialkyl_triaryl_phosphines"},
    "metal_catalyst": {"default_family": "pd_ii_salts"},
    "base": {"default_family": "tertiary_amines_aliphatic"},
    "acid": {"default_family": "mineral_acids"},
    "additive": {"default_family": "quaternary_ammonium_ptc"},
    "condensation_agent": {"default_family": "carbodiimides"},
    "oxidant": {"default_family": "o2_gas"},
    "reductant": {"default_family": "metal_powders"},
    "solvent": {"default_family": "hydrocarbons_aromatic"},
    "other_reagent": {"default_family": "misc_general"},
    "enzyme": {"default_family": None},
    "organo_catalyst": {"default_family": None},
}

DEFAULT_CSV_FIELDS = (
    "name",
    "abbreviation",
    "cas",
    "smiles",
    "formula",
    "type",
    "density",
    "mw",
    "bp",
    "mp",
    "volatile",
    "viscose",
    "role_1",
    "family_1",
    "tag_1",
    "role_2",
    "family_2",
    "tag_2",
)
_LEGACY_FIELDS = {"role", "family_id", "tag", "smile"}
_PROPERTY_FIELDS = {
    "formula",
    "type",
    "density",
    "mw",
    "bp",
    "mp",
    "volatile",
    "viscose",
}
_BASE_FIELDS = {
    "name",
    "abbreviation",
    "cas",
    "smiles",
    "role_1",
    "family_1",
    "tag_1",
    "role_2",
    "family_2",
    "tag_2",
}
_KNOWN_FIELDS = _BASE_FIELDS | _PROPERTY_FIELDS | _LEGACY_FIELDS

ROLE_PAYLOAD_FIELDS: Dict[str, Sequence[str]] = {
    "additive": (),
    "base": ("basicity", "nucleophilicity", "sterics"),
    "condensation_agent": ("strength_band",),
    "ligand": ("donors", "denticity"),
    "oxidant": ("strength_band",),
    "metal_catalyst": ("metal", "oxidation_states", "ligand_type"),
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
    "metal_catalyst": ("metal", "oxidation_states", "ligand_type"),
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


def _normalize_abbreviation(value: Any) -> str:
    if isinstance(value, (list, tuple, set)):
        for item in value:
            if item:
                return str(item).strip()
        return ""
    if value is None:
        return ""
    return str(value).strip()


def _extract_roles(row: Dict[str, str]) -> List[Tuple[str, str, str]]:
    roles: List[Tuple[str, str, str]] = []
    for idx in (1, 2):
        role = (row.get(f"role_{idx}") or "").strip()
        if not role:
            continue
        role = _canonical_role(role)
        family_id = (row.get(f"family_{idx}") or "").strip()
        tag = (row.get(f"tag_{idx}") or "").strip()
        roles.append((role, family_id, tag))
    return roles


def _row_to_entry(row: Dict[str, str]) -> Dict[str, Any]:
    name = (row.get("name") or "").strip()
    abbr = (row.get("abbreviation") or "").strip()
    cas = (row.get("cas") or "").strip()
    smile = (row.get("smiles") or row.get("smile") or "").strip()
    properties: Dict[str, Any] = {}
    for key in _PROPERTY_FIELDS:
        value = row.get(key)
        if value is None:
            continue
        text = str(value).strip()
        if text:
            properties[key] = text
    for key, value in row.items():
        if not key or key in _KNOWN_FIELDS:
            continue
        if value is None:
            continue
        text = str(value).strip()
        if text:
            properties[key] = text

    roles_payload: Dict[str, Dict[str, Any]] = {}
    primary_role = ""
    primary_family = ""
    primary_tag = ""
    for idx, (role, family_id, tag) in enumerate(_extract_roles(row)):
        payload: Dict[str, Any] = {"families": [family_id] if family_id else []}
        if tag:
            payload["tag"] = tag
        roles_payload[role] = payload
        if idx == 0:
            primary_role = role
            primary_family = family_id
            primary_tag = tag

    return {
        "name": name,
        "abbreviation": [abbr] if abbr else [],
        "aliases": [],
        "cas": cas,
        "inchi_key": None,
        "smiles": smile,
        "role": primary_role,
        "family_id": primary_family,
        "tag": primary_tag,
        "roles": roles_payload,
        "properties": properties,
    }


def _role_payload_from_entry(entry: Dict[str, Any], role: str) -> Tuple[str, str]:
    role_payload = (entry.get("roles") or {}).get(role, {})
    families = role_payload.get("families") or []
    family_id = str(families[0]).strip() if families else ""
    tag = (role_payload.get("tag") or "").strip()
    if role == entry.get("role"):
        if not family_id:
            family_id = (entry.get("family_id") or "").strip()
        if not tag:
            tag = (entry.get("tag") or "").strip()
    return family_id, tag


def _entry_to_csv_row(entry: Dict[str, Any], role: Optional[str] = None) -> Dict[str, str]:
    name = (entry.get("name") or "").strip()
    cas = (entry.get("cas") or "").strip()
    smile = (entry.get("smiles") or entry.get("smile") or "").strip()
    abbr = _normalize_abbreviation(entry.get("abbreviation"))
    roles_payload = entry.get("roles") or {}
    roles: List[str] = []
    primary = _canonical_role((entry.get("role") or "").strip())
    if primary:
        roles.append(primary)
    if isinstance(roles_payload, dict):
        for role_key in roles_payload.keys():
            if role_key and role_key not in roles:
                roles.append(role_key)

    if role:
        role = _canonical_role(role)
        if role in roles:
            roles.remove(role)
        roles.insert(0, role)

    if len(roles) > 2:
        roles.sort(key=lambda key: (0 if key == primary else 1, ROLE_PRIORITY.get(key, 99), key))
    role_1 = roles[0] if roles else ""
    role_2 = roles[1] if len(roles) > 1 else ""
    family_1, tag_1 = _role_payload_from_entry(entry, role_1) if role_1 else ("", "")
    family_2, tag_2 = _role_payload_from_entry(entry, role_2) if role_2 else ("", "")

    row = {
        "name": name,
        "abbreviation": abbr,
        "cas": cas,
        "smiles": smile,
        "role_1": role_1,
        "family_1": family_1,
        "tag_1": tag_1,
        "role_2": role_2,
        "family_2": family_2,
        "tag_2": tag_2,
    }
    properties = entry.get("properties") or {}
    if isinstance(properties, dict):
        for key, value in properties.items():
            if not key or key in _LEGACY_FIELDS:
                continue
            row[key] = str(value).strip() if value is not None else ""
    return row

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
    name: str,
    abbreviation: Optional[str],
    cas: str,
    smile: Optional[str],
    role: str,
    family_id: Optional[str],
    tag: Optional[str] = None,
    family_entry: Optional[Dict[str, Any]] = None,
    synonyms: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    role = _canonical_role(role)
    family_id = family_id or ""
    role_payload: Dict[str, Any] = {"families": [family_id] if family_id else []}
    if tag:
        role_payload["tag"] = tag
    entry: Dict[str, Any] = {
        "name": name,
        "abbreviation": [abbreviation] if abbreviation else [],
        "aliases": [],
        "cas": cas,
        "smiles": smile,
        "role": role,
        "family_id": family_id,
        "tag": tag or "",
        "roles": {role: role_payload},
    }
    if family_entry and synonyms is not None:
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

        self.registry_file = self.base_dir / "reagents.csv"
        self.entries: List[Dict[str, Any]] = []
        self.role_entries: Dict[str, List[Dict[str, Any]]] = {}
        self.family_lookup: Dict[str, Tuple[str, Dict[str, Any]]] = {}
        self.family_tokens: Dict[Tuple[str, str], Set[str]] = {}
        self.family_examples: Dict[Tuple[str, str], Set[str]] = {}
        self.family_numeric_baseline: Dict[Tuple[str, str], Optional[Dict[str, Any]]] = {}
        self.cas_index: Dict[str, List[Tuple[str, str, Dict[str, Any]]]] = {}
        self.csv_fields: List[str] = list(DEFAULT_CSV_FIELDS)

        self._load_families()
        self._load_registry()

    # ------------------------------------------------------------------ loading
    def _load_families(self) -> None:
        entries: List[Dict[str, Any]] = []
        try:
            entries = list(load_families_registry_entries())
        except FileNotFoundError:
            entries = []
        if not entries:
            def add_entry(role: str, family: str, definition: str, keywords: Sequence[str]) -> None:
                role = _canonical_role(role)
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
        self.entries = []
        self.role_entries = {}
        self.cas_index = {}
        if not self.registry_file.exists():
            return
        with self.registry_file.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            header = [field for field in (reader.fieldnames or []) if field and field not in _LEGACY_FIELDS]
            extras = [field for field in header if field not in DEFAULT_CSV_FIELDS]
            self.csv_fields = list(DEFAULT_CSV_FIELDS) + extras
            for row in reader:
                entry = _row_to_entry(row)
                self.entries.append(entry)
                roles_payload = entry.get("roles") or {}
                if isinstance(roles_payload, dict):
                    for role in roles_payload.keys():
                        if role:
                            self.role_entries.setdefault(role, []).append(entry)
                cas = entry.get("cas")
                if cas:
                    for role in roles_payload.keys():
                        family_id, _tag = _role_payload_from_entry(entry, role)
                        self.cas_index.setdefault(str(cas), []).append((role or "", family_id, entry))

    def build_role_payload(self, role: str, family_id: str) -> Dict[str, Any]:
        role = _canonical_role(role)
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
        role = _canonical_role(role)
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
        role = _canonical_role(role)
        if role and role not in ROLE_CONFIG:
            raise KeyError(f"Unknown or unsupported role '{role}'.")
        return self.registry_file

    def role_for_family(self, family_id: str) -> Optional[str]:
        data = self.family_lookup.get(family_id)
        return data[0] if data else None

    def family_data(self, role: str, family_id: str) -> Dict[str, Any]:
        role = _canonical_role(role)
        data = self.family_lookup.get(family_id)
        if not data or data[0] != role:
            raise KeyError(f"Family '{family_id}' not found for role '{role}'.")
        return data[1]

    def numeric_baseline(self, role: str, family_id: str) -> Optional[Dict[str, Any]]:
        role = _canonical_role(role)
        return self.family_numeric_baseline.get((role, family_id))

    def family_entry(self, role: str, family_id: str) -> Optional[Dict[str, Any]]:
        role = _canonical_role(role)
        data = self.family_lookup.get(family_id)
        if not data or data[0] != role:
            return None
        return data[1]

    def default_family(self, role: str) -> Optional[str]:
        role = _canonical_role(role)
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
        role = _canonical_role(role)
        familial = self.family_tokens.get((role, family_id), set())
        return bool(familial & set(tokens))

    def find_by_cas(self, cas: str, role: Optional[str] = None) -> Optional[Tuple[str, str, Dict[str, Any]]]:
        records = self.cas_index.get(cas, [])
        role = _canonical_role(role) if role else None
        for entry_role, family, entry in records:
            if role and entry_role != role:
                continue
            return entry_role, family, entry
        return None

    # ------------------------------------------------------------------ updates
    def add_entry(self, role: str, entry: Dict[str, Any]) -> None:
        role = _canonical_role(role)
        normalized = dict(entry)
        normalized["role"] = role
        if not normalized.get("family_id"):
            families = normalized.get("roles", {}).get(role, {}).get("families", [])
            if families:
                normalized["family_id"] = families[0]

        abbr = normalized.get("abbreviation")
        if abbr and not isinstance(abbr, list):
            normalized["abbreviation"] = [abbr]

        cas = normalized.get("cas")
        existing_entry = None
        if cas:
            existing_records = self.cas_index.get(str(cas), [])
            if existing_records:
                existing_entry = existing_records[0][2]

        target = existing_entry or normalized
        if existing_entry:
            if not target.get("name"):
                target["name"] = normalized.get("name")
            if normalized.get("smiles") and not target.get("smiles"):
                target["smiles"] = normalized.get("smiles")
            if normalized.get("abbreviation"):
                target_abbr = target.get("abbreviation") or []
                for item in normalized.get("abbreviation") or []:
                    if item and item not in target_abbr:
                        target_abbr.append(item)
                target["abbreviation"] = target_abbr
            target_roles = target.setdefault("roles", {})
            normalized_roles = normalized.get("roles") or {}
            role_payload = normalized_roles.get(role) or {"families": [normalized.get("family_id") or ""] if normalized.get("family_id") else []}
            if normalized.get("tag") and "tag" not in role_payload:
                role_payload["tag"] = normalized.get("tag")
            target_roles[role] = role_payload
            if not target.get("role"):
                target["role"] = role
            if target.get("role") == role:
                if not target.get("family_id"):
                    target["family_id"] = normalized.get("family_id") or ""
                if normalized.get("tag") and not target.get("tag"):
                    target["tag"] = normalized.get("tag")
        else:
            self.entries.append(target)

        entries = self.role_entries.setdefault(role, [])
        if target not in entries:
            entries.append(target)
            entries.sort(key=lambda item: (item.get("name") or "").lower())

        family = (normalized.get("family_id") or "").strip()
        properties = normalized.get("properties") or {}
        if isinstance(properties, dict):
            for key in properties.keys():
                if key and key not in self.csv_fields and key not in _LEGACY_FIELDS:
                    self.csv_fields.append(key)
        if cas:
            family_id, _tag = _role_payload_from_entry(target, role)
            existing_roles = {item[0] for item in self.cas_index.get(str(cas), [])}
            if role not in existing_roles:
                self.cas_index.setdefault(str(cas), []).append((role, family_id or family, target))
        tokens = tokenize_all([normalized.get("name"), *normalized.get("abbreviation", [])])
        if family:
            self.family_tokens.setdefault((role, family), set()).update(tokens)

    def save_registry(self) -> Path:
        path = self.registry_file
        path.parent.mkdir(parents=True, exist_ok=True)
        rows = [_entry_to_csv_row(entry) for entry in self.entries]
        rows = [row for row in rows if any(row.values())]
        rows.sort(key=lambda row: (row.get("role_1", ""), row.get("name", "").lower(), row.get("cas", "")))
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(self.csv_fields))
            writer.writeheader()
            for row in rows:
                normalized = {field: row.get(field, "") for field in self.csv_fields}
                writer.writerow(normalized)
        return path

    def save_role(self, role: str) -> Path:
        role = _canonical_role(role)
        if role and role not in ROLE_CONFIG:
            raise KeyError(f"Unknown or unsupported role '{role}'.")
        return self.save_registry()

    # ------------------------------------------------------------------ heuristics
    def infer_family(self, role: str, tokens: Set[str]) -> Optional[Tuple[str, List[str]]]:
        role = _canonical_role(role)
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
        role = _canonical_role(role)
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




