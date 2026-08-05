#utf-8
#!/usr/bin/env python3
"""
Simple Qt6 GUI wrapper for processing RDF files only.
Lets the user pick a folder containing RDF files and processes all RDF files in the folder.
Works with PySide6 (preferred) or PyQt6 if installed.

Long reactions are retained; no chemistry is discarded based on SMILES length.
"""
from __future__ import annotations

import os
import sys
import traceback
import re
import csv
import copy
from typing import List, Optional, Dict, Any, Set, Tuple
from pathlib import Path

# Add parent directory to path so we can import chemtools
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from PyQt6 import QtWidgets, QtCore
QtBinding = 'PyQt6'

# Bind Signal/Slot names across PySide6/PyQt6
if hasattr(QtCore, 'Signal') and hasattr(QtCore, 'Slot'):
    Signal = QtCore.Signal
    Slot = QtCore.Slot
elif hasattr(QtCore, 'pyqtSignal') and hasattr(QtCore, 'pyqtSlot'):
    Signal = QtCore.pyqtSignal  # type: ignore[attr-defined]
    Slot = QtCore.pyqtSlot      # type: ignore[attr-defined]
else:  # pragma: no cover
    Signal = None  # type: ignore
    Slot = None    # type: ignore

# Import processing functions from the existing modules
try:
    from process_reactions import parse_rdf, assemble_rows, _molblock_to_smiles
except Exception as e:
    print(f"Error: Cannot import processing helpers: {e}")
    sys.exit(1)


# ---------------------------- Taxonomy integration ----------------------------

import json as _json
from chemtools.reagent.constants import ROLE_ALIASES



FALLBACK_REGISTRY_ENTRIES = {
    "7664-93-9": {
        "Name": "Sulfuric acid",
        "Abbreviation": "H2SO4",
        "Role": "ACID",
        "CategoryHint": "acid",
        "Aliases": ["Sulfuric acid", "H2SO4"],
    },
}



REGISTRY_ROLE_CODE_MAP: Dict[str, str] = {
    "acid": "ACID",
    "additive": "ADDITIVE",
    "base": "BASE",
    "condensation_agent": "COUPLING_REAGENT",
    "ligand": "LIGAND",
    "oxidant": "OXIDANT",
    "metal_catalyst": "METAL_CATALYST",
    "reductant": "REDUCTANT",
    "solvent": "SOLVENT",
    "other_reagent": "OTHER",
    "organo_catalyst": "ORGANO_CATALYST",
    "enzyme": "ENZYME",
}

REGISTRY_CATEGORY_HINT: Dict[str, str] = {
    "acid": "acid",
    "additive": "additive",
    "base": "base",
    "condensation_agent": "coupling_reagent",
    "ligand": "ligand",
    "oxidant": "oxidant",
    "metal_catalyst": "catalyst",
    "reductant": "reductant",
    "solvent": "solvent",
    "other_reagent": "other",
    "organo_catalyst": "catalyst",
    "enzyme": "catalyst",
}

def canonical_role(role: str) -> str:
    """Normalize role labels to canonical registry roles."""
    key = (role or "").lower()
    return ROLE_ALIASES.get(key, key)

class _TaxonomyIndex:
    """Load the reagent registry and expose lookup utilities compatible with the
    previous taxonomy-based workflow."""

    def __init__(self, base_dir: str) -> None:
        raw_base = os.path.abspath(base_dir)
        if not os.path.isdir(raw_base) and os.path.basename(raw_base) == 'compound_taxonomy':
            candidate = os.path.join(os.path.dirname(raw_base), 'reagent_db')
            if os.path.isdir(candidate):
                raw_base = candidate
        if not os.path.isdir(raw_base):
            raise FileNotFoundError(f"Reagent registry directory {raw_base} does not exist")
        self.base_dir = raw_base
        self.cas_map: Dict[str, Dict[str, Any]] = {}
        self.name_to_cas: Dict[str, str] = {}
        self.name_to_role: Dict[str, str] = {}
        self.cas_to_family: Dict[str, str] = {}
        self.family_metal: Dict[str, str] = {}
        self._load_all()

    @staticmethod
    def _norm_name(s: str) -> str:
        import re
        if not s:
            return ""
        t = s.strip().lower()
        t = re.sub(r"\s+", " ", t)
        t = re.sub(r"[^a-z0-9\+\-\.\(\)\[\]/']", "", t)
        return t

    def _load_registry_file(self, filename: str) -> List[Dict[str, Any]]:
        path = os.path.join(self.base_dir, filename)
        if not os.path.exists(path):
            return []
        if filename.lower().endswith(".csv"):
            try:
                with open(path, "r", encoding="utf-8", newline="") as fh:
                    reader = csv.DictReader(fh)
                    return [row for row in reader]
            except Exception:
                return []
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = _json.load(fh)
        except Exception:
            return []
        return data if isinstance(data, list) else []

    def _iter_registry_files(self) -> List[str]:
        csv_path = os.path.join(self.base_dir, "substances.v1.csv")
        if os.path.exists(csv_path):
            return ["substances.v1.csv"]
        return []

    def _index_member(
        self,
        *,
        cas: str | None,
        name: str | None,
        abbr: str | None,
        synonyms: List[str] | None,
        role: str,
        category_hint: str,
        generic_core: str = "",
        family_id: str | None = None,
        role_payload: Optional[Dict[str, Any]] = None,
    ) -> None:
        cas = (cas or "").strip()
        name = (name or "").strip()
        abbr = (abbr or "").strip()
        norm_synonyms: List[str] = []
        seen_syn: Set[str] = set()
        for syn in synonyms or []:
            syn_str = str(syn).strip()
            if not syn_str:
                continue
            key = syn_str.lower()
            if key in seen_syn:
                continue
            seen_syn.add(key)
            norm_synonyms.append(syn_str)

        token_source = (abbr or name or cas).strip()
        token = token_source.replace(" ", "")
        if not token:
            token = cas.replace(" ", "") if cas else ""

        rec: Dict[str, Any] = {
            "Name": name or abbr or cas,
            "Abbreviation": abbr,
            "Token": token,
            "Role": role,
            "CategoryHint": category_hint,
            "GenericCore": generic_core,
        }
        if family_id:
            rec["FamilyID"] = family_id
        if role_payload:
            rec["RolePayload"] = role_payload
            metal_val = role_payload.get("metal")
            if isinstance(metal_val, list):
                metal_joined = ", ".join(str(m).strip() for m in metal_val if m)
            elif metal_val is not None:
                metal_joined = str(metal_val).strip()
            else:
                metal_joined = ""
            if metal_joined:
                rec["Metal"] = metal_joined
            if "oxidation_states" in role_payload:
                rec["OxidationStates"] = role_payload.get("oxidation_states")

        if cas:
            existing = self.cas_map.get(cas)
            should_replace = False
            if not existing:
                should_replace = True
            else:
                existing_core = (existing.get("GenericCore") or "").strip()
                new_core = generic_core.strip() if isinstance(generic_core, str) else str(generic_core).strip()
                existing_family = existing.get("FamilyID")
                if new_core and not existing_core:
                    should_replace = True
                elif family_id and not existing_family:
                    should_replace = True
                elif existing.get("Role") != role and role == "CAT":
                    should_replace = True
            if should_replace:
                self.cas_map[cas] = rec
            if family_id:
                self.cas_to_family[cas] = family_id
            else:
                self.cas_to_family.setdefault(cas, "")

        keys = [name, abbr, token] + norm_synonyms
        for k in keys:
            k_norm = self._norm_name(k or "")
            if not k_norm:
                continue
            if cas and k_norm not in self.name_to_cas:
                self.name_to_cas[k_norm] = cas
            self.name_to_role.setdefault(k_norm, role)

        if family_id and generic_core:
            self.family_metal.setdefault(family_id, generic_core)

    def _load_all(self) -> None:
        for filename in self._iter_registry_files():
            entries = self._load_registry_file(filename)
            if not entries:
                continue
            for entry in entries:
                cas = entry.get("cas")
                name = entry.get("name")
                abbr = (entry.get("abbreviation") or "").strip()
                canonical = canonical_role(entry.get("role_1") or "")
                family_id = entry.get("family_1") or None
                if not canonical:
                    canonical = canonical_role(entry.get("role_2") or "")
                    family_id = entry.get("family_2") or None
                if not canonical:
                    canonical = canonical_role(entry.get("role") or "")
                    family_id = entry.get("family_id") or None
                if not canonical:
                    continue
                role_code = REGISTRY_ROLE_CODE_MAP.get(canonical, canonical.upper())
                category_hint = REGISTRY_CATEGORY_HINT.get(canonical, canonical)
                generic_core = ""  # CSV registry does not store metal-specific payloads.
                self._index_member(
                    cas=cas,
                    name=name,
                    abbr=abbr,
                    synonyms=[],
                    role=role_code,
                    category_hint=category_hint,
                    generic_core=generic_core,
                    family_id=family_id,
                    role_payload=None,
                )

        self._apply_fallback_entries()

    def _apply_fallback_entries(self) -> None:
        for cas, data in FALLBACK_REGISTRY_ENTRIES.items():
            if cas in self.cas_map:
                continue
            name = data.get("Name")
            abbr = data.get("Abbreviation")
            aliases = data.get("Aliases") or []
            role = data.get("Role", "UNK")
            category_hint = data.get("CategoryHint") or role.lower()
            generic_core = data.get("GenericCore", "")
            family_id = data.get("FamilyID")
            token_override = data.get("Token") or None
            self._index_member(
                cas=cas,
                name=name,
                abbr=abbr,
                synonyms=aliases,
                role=role,
                category_hint=category_hint,
                generic_core=generic_core,
                family_id=family_id,
                role_payload=None,
            )
            if token_override and cas in self.cas_map:
                self.cas_map[cas]["Token"] = token_override

    # ---- Public helpers ----

    def cas_lookup(self, cas: str) -> Dict[str, Any]:
        return self.cas_map.get((cas or "").strip(), {})

    def role_for(self, name: str | None, cas: str | None) -> str:
        cas = (cas or "").strip()
        if cas:
            role = (self.cas_map.get(cas, {}).get("Role") or "").strip().upper()
            if role:
                return role
        n = self._norm_name(name or "")
        if n:
            role = (self.name_to_role.get(n) or "").strip().upper()
            if role:
                return role
        fallback_role = self._heuristic_role(name, cas)
        if fallback_role != "UNK":
            return fallback_role
        return "UNK"

    def _heuristic_role(self, name: str | None, cas: str | None) -> str:
        cas_key = (cas or "").strip()
        if cas_key in FALLBACK_REGISTRY_ENTRIES:
            return FALLBACK_REGISTRY_ENTRIES[cas_key].get("Role", "UNK").upper()
        norm_name = self._norm_name(name or "")
        if norm_name and "acid" in norm_name:
            return "ACID"
        return "UNK"

    def ligand_token_for(self, name: str | None, cas: str | None) -> str:
        cas = (cas or "").strip()
        if cas and cas in self.cas_map:
            tok = (
                self.cas_map[cas].get("Token")
                or self.cas_map[cas].get("Abbreviation")
                or self.cas_map[cas].get("Name")
                or ""
            ).strip()
            return tok.replace(" ", "")
        n = (name or "").strip()
        if n:
            nn = self._norm_name(n)
            cas2 = self.name_to_cas.get(nn)
            if cas2 and cas2 in self.cas_map:
                tok = (
                    self.cas_map[cas2].get("Token")
                    or self.cas_map[cas2].get("Abbreviation")
                    or self.cas_map[cas2].get("Name")
                    or ""
                ).strip()
                return tok.replace(" ", "")
            return n.replace(" ", "")
        return ""

    def metal_token_from_core_pairs(self, core_pairs: List[str], fallback_generic: List[str]) -> str:
        for p in core_pairs or []:
            nm, _, cs = p.partition("|")
            cs = cs.strip()
            if cs and cs in self.cas_map:
                gen = (self.cas_map[cs].get("GenericCore") or "").strip()
                if gen:
                    return gen
        for g in fallback_generic or []:
            if g:
                return g
        return ""


class ReactionMarkdownGenerator:  # taxonomy-aware local generator
    def __init__(self, taxonomy: _TaxonomyIndex | None = None) -> None:
        self.taxonomy = taxonomy or None
        # For backward compatibility with existing code path
        self.cas_map = (taxonomy.cas_map if taxonomy else {})

    def _safe_json_list(self, val: str):
        try:
            return _json.loads(val or "[]")
        except Exception:
            return []

    def _format_cas_list(self, items: Any) -> str:
        if not items:
            return ""
        if isinstance(items, str):
            return items.strip()
        if not isinstance(items, list):
            return str(items).strip()
        seen: Set[str] = set()
        out: List[str] = []
        for item in items:
            if item is None:
                continue
            val = str(item).strip()
            if not val or val in seen:
                continue
            seen.add(val)
            out.append(val)
        return ", ".join(out)

    def _pair_to_obj(self, item: str):
        if "|" in item:
            name, cas = item.split("|", 1)
            return {"name": name.strip(), "cas": cas.strip()}
        return {"name": item.strip(), "cas": ""}

    def _component_from_item(self, item: Any) -> Dict[str, Any]:
        if isinstance(item, dict):
            return {k: v for k, v in item.items()}
        if isinstance(item, str):
            return self._pair_to_obj(item)
        return {"name": str(item).strip(), "cas": ""}

    def _lookup_registry_record(self, name: str, cas: str) -> Tuple[Dict[str, Any], str]:
        cas_key = (cas or "").strip()
        record: Dict[str, Any] = {}
        if cas_key:
            record = self.cas_map.get(cas_key, {})
        if (not record) and name:
            taxonomy = getattr(self, "taxonomy", None)
            if taxonomy:
                norm_func = getattr(taxonomy, "_norm_name", None)
                norm_name = norm_func(name) if callable(norm_func) else name.lower().strip()
                alt_cas = taxonomy.name_to_cas.get(norm_name)
                if alt_cas:
                    maybe_record = self.cas_map.get(alt_cas, {})
                    if maybe_record:
                        record = maybe_record
                        cas_key = alt_cas
        return record or {}, cas_key

    def _normalize_component(self, component: Dict[str, Any]) -> Dict[str, Any]:
        base = {k: v for k, v in component.items() if v is not None}
        name = (base.get("name") or "").strip()
        cas = (base.get("cas") or "").strip()
        record, resolved_cas = self._lookup_registry_record(name, cas)

        abbreviation = (record.get("Abbreviation") or "").strip()
        token = (record.get("Token") or "").strip()
        registry_name = (record.get("Name") or "").strip()

        display_name = (
            abbreviation
            or token
            or name
            or registry_name
            or resolved_cas
            or "?"
        ).strip()
        if not display_name:
            display_name = "?"

        normalized = dict(base)
        if resolved_cas:
            normalized["cas"] = resolved_cas

        normalized["name"] = display_name

        if abbreviation:
            normalized["abbreviation"] = abbreviation
        else:
            normalized.pop("abbreviation", None)

        # NOTE: Removed "original_name" field to save ~8.5% space per record.
        # The full chemical name can be looked up via CAS number if needed.
        # Previous behavior: stored full name like "Palladium(II) acetate"
        # when it differed from display name "Pd(OAc)2".

        return normalized

    def _join_names(self, arr):
        if not arr:
            return ""
        out = []
        for it in arr:
            component = self._component_from_item(it)
            normalized = self._normalize_component(component)
            nm = (normalized.get("name") or "").strip()
            if nm:
                out.append(nm)
        return ", ".join(out)

    _METAL_KEYWORDS: Dict[str, str] = {
        "nickel": "Ni",
        "palladium": "Pd",
        "platinum": "Pt",
        "copper": "Cu",
        "ruthenium": "Ru",
        "rhodium": "Rh",
        "iridium": "Ir",
        "cobalt": "Co",
        "iron": "Fe",
        "manganese": "Mn",
        "chromium": "Cr",
        "molybdenum": "Mo",
        "tantalum": "Ta",
        "niobium": "Nb",
        "vanadium": "V",
        "zinc": "Zn",
        "silver": "Ag",
        "gold": "Au",
        "tin": "Sn",
        "antimony": "Sb",
        "bismuth": "Bi",
        "osmium": "Os",
        "rhenium": "Re"
    }


    _METAL_SYMBOLS_SET: Set[str] = ({sym.lower() for sym in _METAL_KEYWORDS.values()} | {"pd", "ni", "cu", "ir", "ru", "rh", "fe", "co", "mn", "mo", "pt", "au", "ag", "zn"})
    _PREFORMED_IGNORE_TERMS: Set[str] = {"allyl", "cl", "br", "chloride", "iodide", "bromide", "std", "standard", "precatalyst", "precursor", "complex", "hbpin", "b2pin2", "bpin", "cod", "pf6", "bf4", "otf", "counterion", "g1", "g2", "g3", "g4", "g5", "dimer"}
    _PREFORMED_IGNORE_SUBSTRINGS: Tuple[str, ...] = ("borylation", "precatalyst", "precursor", "complex")

    def _extract_names(self, items: Optional[List[Any]] | None) -> List[str]:
        out: List[str] = []
        if not items:
            return out
        for element in items:
            if isinstance(element, dict):
                name = (element.get("name") or element.get("abbr") or element.get("token") or "").strip()
            else:
                text_val = str(element)
                if "|" in text_val:
                    name = text_val.split("|", 1)[0].strip()
                else:
                    name = text_val.strip()
            if name:
                out.append(name)
        return out

    def _infer_metal_from_names(self, *item_lists: Optional[List[Any]]) -> str:
        names: List[str] = []
        for items in item_lists:
            names.extend(self._extract_names(items))
        for name in names:
            lower_name = name.lower()
            for keyword, symbol in self._METAL_KEYWORDS.items():
                if keyword in lower_name:
                    return symbol
        return ""


    @staticmethod
    def _sanitize_token_value(token: str) -> str:
        if not token:
            return ""
        cleaned = str(token)
        replacements = [
            ("\u2013", "-"),  # en dash
            ("\u2014", "-"),  # em dash
            ("\u2212", "-"),  # minus sign
            ("\u2215", "-"),  # division slash
            ("\u00B7", ""),   # middle dot
            ("\u00A0", " "),  # non-breaking space
        ]
        for src, dest in replacements:
            cleaned = cleaned.replace(src, dest)
        cleaned = cleaned.replace(" ", "")
        return cleaned.strip("-")

    @classmethod
    def _cleanup_preformed_segment(cls, segment: str) -> str:
        if not segment:
            return ""
        segment = (
            segment.replace("\u2013", "-")
            .replace("\u2014", "-")
            .replace("\u2212", "-")
        )
        segment = segment.replace("\u00B7", " ")
        segment = segment.strip()
        for ch in "[]{}":
            segment = segment.replace(ch, " ")
        segment = segment.replace(",", " ")
        segment = segment.strip()
        return cls._sanitize_token_value(segment)

    @classmethod
    def _is_ligand_candidate(cls, segment: str) -> bool:
        if not segment:
            return False
        lower = segment.lower()
        if not any(ch.isalpha() for ch in segment):
            return False
        if lower in cls._PREFORMED_IGNORE_TERMS:
            return False
        if any(sub in lower for sub in cls._PREFORMED_IGNORE_SUBSTRINGS):
            return False
        if lower in cls._METAL_SYMBOLS_SET:
            return False
        if lower.startswith('g') and lower[1:].isdigit():
            return False
        return True

    @classmethod
    def _parse_ligand_from_text(cls, text: str) -> str:
        if not text:
            return ""
        normalized = (
            text.replace("\u2013", "-")
            .replace("\u2014", "-")
            .replace("\u2212", "-")
        )
        for content in re.findall(r'\(([^()]+)\)', normalized):
            candidate = cls._cleanup_preformed_segment(content)
            if cls._is_ligand_candidate(candidate):
                return candidate
        plus_normalized = normalized.replace('\u00B7', '+').replace('\u00A0', '+').replace('/', '+')
        segments = [seg.strip() for seg in plus_normalized.split('+') if seg.strip()]
        if len(segments) > 1:
            for seg in segments[1:]:
                candidate = cls._cleanup_preformed_segment(seg)
                if cls._is_ligand_candidate(candidate):
                    return candidate
        hyphen_segments = [seg for seg in re.split(r'-', normalized) if seg]
        ligand_parts = [cls._cleanup_preformed_segment(seg) for seg in hyphen_segments]
        ligand_parts = [seg for seg in ligand_parts if cls._is_ligand_candidate(seg)]
        if ligand_parts:
            return '-'.join(ligand_parts)
        return ""

    def _extract_ligand_from_preformed(self, rec: Dict[str, Any], name_hint: str) -> str:
        candidates: List[str] = []
        if name_hint:
            candidates.append(name_hint)
        for key in ("Abbreviation", "Name", "Token"):
            value = rec.get(key)
            if isinstance(value, str) and value:
                candidates.append(value)
            elif isinstance(value, list):
                candidates.extend(str(v) for v in value if v)
        for text in candidates:
            ligand = self._parse_ligand_from_text(text)
            if ligand:
                return self._sanitize_token_value(ligand)
        return ""

    def _detect_preformed_components(self, items: Optional[List[Any]]) -> Optional[Dict[str, str]]:
        if not self.taxonomy or not items:
            return None
        for element in items:
            obj = element if isinstance(element, dict) else self._pair_to_obj(element)
            name = (obj.get('name') or obj.get('abbr') or obj.get('token') or '').strip()
            cas = (obj.get('cas') or '').strip()
            rec = None
            if cas and cas in self.taxonomy.cas_map:
                rec = self.taxonomy.cas_map.get(cas)
            elif name:
                norm_name = self.taxonomy._norm_name(name)
                alt_cas = self.taxonomy.name_to_cas.get(norm_name)
                if alt_cas:
                    cas = alt_cas
                    rec = self.taxonomy.cas_map.get(alt_cas)
            if not rec or rec.get('Role') != 'METAL_CATALYST':
                continue
            metal = (rec.get('GenericCore') or rec.get('Metal') or '').strip()
            ligand = self._extract_ligand_from_preformed(rec, name)
            token_value = rec.get('Token') or rec.get('Abbreviation') or rec.get('Name') or name
            if isinstance(token_value, list):
                token_value = token_value[0] if token_value else ''
            token_clean = self._sanitize_token_value(token_value or '')
            return {'metal': metal, 'ligand': ligand, 'token': token_clean}
        return None

    def _collect_full_catalytic_system(self, row: Dict[str, Any]) -> List[Any]:
        full_system = self._safe_json_list(row.get("FullCatalyticSystem", "[]"))
        if full_system:
            return list(full_system)
        core_pairs = self._safe_json_list(row.get("CatalystCoreDetail", "[]"))
        lig_list = self._safe_json_list(row.get("Ligand", "[]"))
        combined: List[Any] = []
        if core_pairs:
            combined.extend(core_pairs)
        if lig_list:
            combined.extend(lig_list)
        return combined


    def _compute_condition_core(
            self, row: Dict[str, Any], full_system_list: Optional[List[Any]] | None = None
        ) -> str:
            core_gen = self._safe_json_list(row.get("CatalystCoreGeneric", "[]"))
            core_pairs = self._safe_json_list(row.get("CatalystCoreDetail", "[]"))
            lig_list = self._safe_json_list(row.get("Ligand", "[]"))

            if full_system_list is None:
                full_system_list = self._safe_json_list(row.get("FullCatalyticSystem", "[]"))

            preformed_info = None
            if self.taxonomy:
                preformed_info = self._detect_preformed_components(full_system_list)
                if not preformed_info:
                    preformed_info = self._detect_preformed_components(core_pairs)

            metal = ""
            if self.taxonomy:
                metal = self.taxonomy.metal_token_from_core_pairs(core_pairs, core_gen)
            if not metal and preformed_info:
                metal = (preformed_info.get("metal") or "").strip()
            if not metal:
                metal = (core_gen[0] if core_gen else "").strip()

            alt_metal = self._infer_metal_from_names(core_pairs, full_system_list, lig_list)
            if alt_metal and (not metal or alt_metal.lower() != metal.lower()):
                metal = alt_metal

            lig_tok = ""
            if lig_list:
                first = lig_list[0]
                cas = ""
                nm = ""
                if isinstance(first, str) and "|" in first:
                    nm, cas = (first.split("|", 1) + [""])[:2]
                elif isinstance(first, dict):
                    nm = first.get("name") or first.get("abbr") or first.get("token") or ""
                    cas = first.get("cas") or ""
                else:
                    nm = first
                if self.taxonomy:
                    lig_tok = self.taxonomy.ligand_token_for(nm, cas)
                else:
                    rec = (getattr(self, "cas_map", {}) or {}).get((cas or "").strip()) or {}
                    lig_tok = (rec.get("Token") or rec.get("Abbreviation") or rec.get("Name") or nm).strip()
                if lig_tok:
                    lig_tok = self._sanitize_token_value(lig_tok)

            if not lig_tok and preformed_info:
                ligand_candidate = preformed_info.get("ligand") or ""
                if ligand_candidate:
                    lig_tok = self._sanitize_token_value(ligand_candidate)

            # If no ligand but a metal core is present, try suggestion via PairingHelper
            if not lig_tok and metal and self.taxonomy and not preformed_info:
                # Try to identify catalyst family from core CAS
                fam_id = None
                for p in core_pairs or []:
                    _, _, cs = p.partition("|")
                    cs = cs.strip()
                    if cs and cs in self.taxonomy.cas_to_family:
                        fam_id = self.taxonomy.cas_to_family.get(cs)
                        break
                try:
                    from conditioncore_pairing_helper_for_ref_only import PairingHelper  # type: ignore
                    from condition_registry.loader import TAXONOMY_PATH
                    cat_path = str(TAXONOMY_PATH)
                    ph = PairingHelper(cat_path)
                    if fam_id:
                        hint = (row.get("ReactionType") or "").strip() or None
                        pick = ph.suggest_for(fam_id, None, hint)
                        if pick and pick.get("ligand", {}).get("abbr"):
                            lig_tok = self._sanitize_token_value(pick["ligand"]["abbr"] or "")
                except Exception:
                    pass

            preformed_token = ""
            if preformed_info:
                preformed_token = self._sanitize_token_value(preformed_info.get("token") or "")

            metal_output = self._sanitize_token_value(metal)

            if metal_output and lig_tok:
                return f"{metal_output}/{lig_tok}"
            if lig_tok:
                return lig_tok
            if preformed_token:
                return preformed_token
            if metal_output:
                return metal_output

            # Fallback: highlight prominent reagents (acid/additive/base/etc.) when no metal/ligand core exists
            reag_list = self._safe_json_list(row.get("Reagent", "[]"))
            role_list = self._safe_json_list(row.get("ReagentRole", "[]"))

            def _token_from_item(item_obj: Dict[str, Any]) -> str:
                cas = (item_obj.get("cas") or "").strip()
                name = (item_obj.get("name") or "").strip()
                token = ""
                if self.taxonomy and cas:
                    rec = self.taxonomy.cas_map.get(cas, {})
                    token = (rec.get("Token") or rec.get("Abbreviation") or rec.get("Name") or "").strip()
                if not token:
                    token = name or cas
                return self._sanitize_token_value(token)

            fallback_roles = (
                ("ACID", "Acid"),
                ("BASE", "Base"),
                ("ADDITIVE", "Additive"),
                ("COUPLING_REAGENT", "Coupling reagent"),
                ("OXIDANT", "Oxidant"),
                ("REDUCTANT", "Reductant"),
                ("OTHER", "Other reagent"),
                ("ORGANO_CATALYST", "Organo-catalyst"),
                ("ENZYME", "Enzyme"),
            )

            # Collect all reagent tokens for uncatalyzed reactions
            reagent_tokens = []
            for role_key, label in fallback_roles:
                for idx, item in enumerate(reag_list or []):
                    role = (role_list[idx] if idx < len(role_list) else "").upper()
                    if role != role_key:
                        continue
                    obj = item if isinstance(item, dict) else self._pair_to_obj(item)
                    token = _token_from_item(obj)
                    if token and token not in reagent_tokens:
                        reagent_tokens.append(token)
            
            # Return all reagent names separated by "/"
            if reagent_tokens:
                return "/".join(reagent_tokens)

            return ""

    def generate_markdown_report(self, rows, output_path: str, source_folder: str):
        try:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(f"# Reactions Report ({source_folder})\n\n")
                f.write(f"Total reactions: {len(rows)}\n\n")
                for row in rows:
                    rid = row.get("ReactionID", "")
                    rtype = row.get("ReactionType", "")
                    reag_list = self._safe_json_list(row.get("Reagent", "[]"))
                    role_list = self._safe_json_list(row.get("ReagentRole", "[]"))
                    solv_list = self._safe_json_list(row.get("Solvent", "[]"))

                    full_system_list = self._collect_full_catalytic_system(row)
                    catalytic_system = self._join_names(full_system_list)

                    T = row.get("Temperature_C", "")
                    t = row.get("Time_h", "")
                    y = row.get("Yield_%", "")

                    reag_out = []
                    for i, item in enumerate(reag_list):
                        obj = self._normalize_component(self._component_from_item(item))
                        role = (role_list[i] if i < len(role_list) else "").upper() or "ADDITIVE"
                        seg = obj.get("name") or obj.get("cas") or "?"
                        if obj.get("cas"):
                            seg += f" ({obj['cas']})"
                        reag_out.append(f"{seg} [{role}]")

                    solv_out = []
                    for item in solv_list:
                        obj = self._normalize_component(self._component_from_item(item))
                        seg = obj.get("name") or obj.get("cas") or "?"
                        if obj.get("cas"):
                            seg += f" ({obj['cas']})"
                        solv_out.append(seg)

                    r_smi = row.get("ReactantSMILES", "")
                    p_smi = row.get("ProductSMILES", "")

                    f.write(f"## Reaction {rid}\n\n")
                    if rtype:
                        f.write(f"- Type: {rtype}\n")
                    if catalytic_system:
                        f.write(f"- Catalytic System: {catalytic_system}\n")
                    if y != "":
                        f.write(f"- Yield %: {y}\n")
                    if T != "":
                        f.write(f"- Temperature (C): {T}\n")
                    if t != "":
                        f.write(f"- Time (h): {t}\n")
                    if reag_out:
                        f.write(f"- Reagents: {', '.join(reag_out)}\n")
                    if solv_out:
                        f.write(f"- Solvents: {', '.join(solv_out)}\n")
                    if r_smi or p_smi:
                        f.write(f"- SMILES: {r_smi}>>{p_smi}\n")
                    f.write("\n")
        except Exception:
            pass

    def generate_csv_export(self, rows, output_path: str, source_folder: str, rdf_map: Optional[Dict[str, Dict[str, Any]]] = None):
        """Write every parsed reaction and its chemistry/provenance fields."""
        fieldnames = [
            "reaction_id", "duplicate_reaction_ids", "source_files", "reaction_type",
            "yield_pct", "temperature_c", "time_h", "reactant_smiles",
            "product_smiles", "reaction_smiles", "title", "authors", "citation",
            "reference", "reactant_cas", "product_cas", "reagent_cas",
            "catalyst_cas", "solvent_cas", "reactant_amd", "product_amd",
            "reagent_amd", "catalyst_amd", "solvent_amd",
            "experimental_procedure", "stages", "steps", "product_yield_1",
            "product_yield_2", "product_yield_3", "product_yield_4",
            "product_yield_5", "product_yield_6", "product_yield_7",
            "product_yields_json",
            "reactants_json", "products_json", "reagents_json", "catalysts_json",
            "solvents_json", "notes", "structure_warnings",
        ]

        def _format_notes(value: Any) -> str:
            if not value:
                return ""
            if isinstance(value, list):
                return " | ".join(str(item).strip() for item in value if str(item).strip())
            return str(value).strip()

        def _component_records(
            rdf_data: Dict[str, Any],
            role: str,
            mol_field: Optional[str] = None,
        ) -> List[Dict[str, Any]]:
            components: Dict[int, Dict[str, Any]] = {}
            pattern = re.compile(rf":{role}\((\d+)\):(CAS_RN|AMD|YIELD|CTAB|MOL)$")
            for key, values in (rdf_data.get("raw_fields") or {}).items():
                match = pattern.search(key)
                if not match:
                    continue
                index = int(match.group(1))
                field = match.group(2).lower()
                component = components.setdefault(index, {"index": index})
                normalized = list(values or [])
                if field in {"cas_rn", "amd"}:
                    component[field] = normalized[0] if len(normalized) == 1 else normalized
                elif field == "yield":
                    component["yield_pct"] = normalized[0] if len(normalized) == 1 else normalized

            if mol_field:
                for index, mol_block in enumerate(rdf_data.get(mol_field) or [], 1):
                    component = components.setdefault(index, {"index": index})
                    smiles = _molblock_to_smiles(mol_block)
                    if smiles:
                        component["smiles"] = smiles
                    else:
                        component["unparsed_mol_block"] = mol_block
            return [components[index] for index in sorted(components)]

        try:
            written = 0
            with open(output_path, "w", newline="", encoding="utf-8") as output_file:
                writer = csv.DictWriter(output_file, fieldnames=fieldnames)
                writer.writeheader()
                for idx, row in enumerate(rows, 1):
                    reaction_id = row.get("ReactionID", "")
                    rdf_data = (
                        rdf_map.get(reaction_id, {})
                        if rdf_map and reaction_id in rdf_map
                        else {}
                    )
                    reactant_smiles = row.get("ReactantSMILES", "")
                    product_smiles = row.get("ProductSMILES", "")
                    reaction_smiles = (
                        f"{reactant_smiles}>>{product_smiles}"
                        if reactant_smiles or product_smiles else ""
                    )
                    warnings: List[str] = []
                    if rdf_data.get("rct_mol") and not reactant_smiles:
                        warnings.append("reactant_mol_to_smiles_failed")
                    if rdf_data.get("pro_mol") and not product_smiles:
                        warnings.append("product_mol_to_smiles_failed")
                    if RDKIT_AVAILABLE and reactant_smiles and not _is_valid_smiles(reactant_smiles):
                        warnings.append("reactant_smiles_validation_failed")
                    if RDKIT_AVAILABLE and product_smiles and not _is_valid_smiles(product_smiles):
                        warnings.append("product_smiles_validation_failed")

                    product_yields = rdf_data.get("pro_yields", {})

                    def _yield_at(index: int) -> Any:
                        if isinstance(product_yields, dict):
                            value = product_yields.get(index)
                            if value is None:
                                value = product_yields.get(str(index))
                            return "" if value is None else value
                        if isinstance(product_yields, list) and len(product_yields) >= index:
                            return product_yields[index - 1]
                        return ""

                    entry = {
                        "reaction_id": reaction_id,
                        "duplicate_reaction_ids": self._format_cas_list(
                            rdf_data.get("duplicate_reaction_ids", [reaction_id])
                        ),
                        "source_files": self._format_cas_list(rdf_data.get("source_files", [])),
                        "reaction_type": row.get("ReactionType", ""),
                        "yield_pct": row.get("Yield_%", ""),
                        "temperature_c": row.get("Temperature_C", ""),
                        "time_h": row.get("Time_h", ""),
                        "reactant_smiles": reactant_smiles,
                        "product_smiles": product_smiles,
                        "reaction_smiles": reaction_smiles,
                        "title": rdf_data.get("title", "") or "",
                        "authors": rdf_data.get("authors", "") or "",
                        "citation": rdf_data.get("citation", "") or "",
                        "reference": row.get("Reference", ""),
                        "reactant_cas": self._format_cas_list(rdf_data.get("rct_cas", [])),
                        "product_cas": self._format_cas_list(rdf_data.get("pro_cas", [])),
                        "reagent_cas": self._format_cas_list(rdf_data.get("rgt_cas", [])),
                        "catalyst_cas": self._format_cas_list(rdf_data.get("cat_cas", [])),
                        "solvent_cas": self._format_cas_list(rdf_data.get("sol_cas", [])),
                        "reactant_amd": self._format_cas_list(rdf_data.get("rct_amd", [])),
                        "product_amd": self._format_cas_list(rdf_data.get("pro_amd", [])),
                        "reagent_amd": self._format_cas_list(rdf_data.get("rgt_amd", [])),
                        "catalyst_amd": self._format_cas_list(rdf_data.get("cat_amd", [])),
                        "solvent_amd": self._format_cas_list(rdf_data.get("sol_amd", [])),
                        "experimental_procedure": _format_notes(rdf_data.get("exp_proc", "")),
                        "stages": rdf_data.get("stages", "") or "",
                        "steps": rdf_data.get("steps", "") or "",
                        "product_yield_1": _yield_at(1),
                        "product_yield_2": _yield_at(2),
                        "product_yield_3": _yield_at(3),
                        "product_yield_4": _yield_at(4),
                        "product_yield_5": _yield_at(5),
                        "product_yield_6": _yield_at(6),
                        "product_yield_7": _yield_at(7),
                        "product_yields_json": _json.dumps(
                            product_yields, ensure_ascii=False, sort_keys=True
                        ),
                        "reactants_json": _json.dumps(_component_records(rdf_data, "RCT", "rct_mol"), ensure_ascii=False),
                        "products_json": _json.dumps(_component_records(rdf_data, "PRO", "pro_mol"), ensure_ascii=False),
                        "reagents_json": _json.dumps(_component_records(rdf_data, "RGT"), ensure_ascii=False),
                        "catalysts_json": _json.dumps(_component_records(rdf_data, "CAT"), ensure_ascii=False),
                        "solvents_json": _json.dumps(_component_records(rdf_data, "SOL"), ensure_ascii=False),
                        "notes": _format_notes(rdf_data.get("notes", "")),
                        "structure_warnings": " | ".join(warnings),
                    }
                    writer.writerow(entry)
                    written += 1
                    if idx % 100 == 0:
                        print(f"  Processed {idx}/{len(rows)} reactions for CSV...")
            print(f"  Successfully wrote {written} reactions to {output_path}")
        except Exception as error:
            print(f"Error writing CSV: {error}")

# Detect RDKit availability
try:
    from rdkit import Chem  # type: ignore
    RDKIT_AVAILABLE = True
except Exception:
    Chem = None  # type: ignore
    RDKIT_AVAILABLE = False


def _is_valid_smiles(smiles: str) -> bool:
    """Check if a SMILES string is chemically valid using RDKit."""
    if not smiles or not RDKIT_AVAILABLE:
        return False
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        return mol is not None
    except Exception:
        return False


RDF_LIST_FIELDS = (
    "rct_cas", "pro_cas", "rgt_cas", "cat_cas", "sol_cas",
    "rct_amd", "pro_amd", "rgt_amd", "cat_amd", "sol_amd",
    "exp_proc", "notes", "rct_mol", "pro_mol", "source_files",
    "reaction_types", "duplicate_reaction_ids",
)


def merge_rdf_records(primary: Dict[str, Any], duplicate: Dict[str, Any]) -> Dict[str, Any]:
    """Merge duplicate RDF records while preserving every distinct value.

    The first non-empty scalar remains the compatibility value. Conflicting
    scalar and per-product yield values are retained in explicit multi-value
    provenance fields instead of being silently overwritten.
    """
    merged = copy.deepcopy(primary)
    for field in RDF_LIST_FIELDS:
        values = list(merged.get(field) or [])
        for value in duplicate.get(field) or []:
            if value not in values:
                values.append(value)
        merged[field] = values

    merged_raw = merged.setdefault("raw_fields", {})
    for key, values in (duplicate.get("raw_fields") or {}).items():
        target = merged_raw.setdefault(key, [])
        for value in values or []:
            if value not in target:
                target.append(value)

    merged_yields = merged.setdefault("pro_yields", {})
    yield_values = merged.setdefault("product_yield_values", {})
    for index, value in (duplicate.get("pro_yields") or {}).items():
        if index not in merged_yields:
            merged_yields[index] = value
        elif merged_yields[index] != value:
            values = yield_values.setdefault(str(index), [merged_yields[index]])
            if value not in values:
                values.append(value)

    ignored = set(RDF_LIST_FIELDS) | {
        "raw_fields", "pro_yields", "product_yield_values", "field_values",
    }
    field_values = merged.setdefault("field_values", {})
    for field, value in duplicate.items():
        if field in ignored or value in (None, "", [], {}):
            continue
        current = merged.get(field)
        if current in (None, "", [], {}):
            merged[field] = copy.deepcopy(value)
        elif current != value:
            values = field_values.setdefault(field, [current])
            if value not in values:
                values.append(copy.deepcopy(value))

    if not field_values:
        merged.pop("field_values", None)
    if not yield_values:
        merged.pop("product_yield_values", None)
    return merged


class RDFWorker(QtCore.QObject):
    finished = Signal(bool, str) if Signal else None  # type: ignore[misc]
    progress = Signal(str) if Signal else None  # type: ignore[misc]

    def __init__(
        self,
        folder_path: str,
        output_md_path: str,
        output_csv_path: str,
        combine_all: bool = False,
    ):
        super().__init__()
        self.folder_path = folder_path
        self.output_md_path = output_md_path
        self.output_csv_path = output_csv_path
        self.combine_all = combine_all
        self.rdf_files = []

    def _emit(self, msg: str):
        """Emit progress message"""
        sig = getattr(self, 'progress', None)
        if sig:
            try:
                sig.emit(msg)
            except Exception:
                pass

    def _get_reaction_type_subfolders(self) -> List[Tuple[str, str]]:
        """Detect if we have multiple reaction type subfolders to process separately.
        
        Returns:
            List of (subfolder_path, reaction_type_name) tuples.
            Empty list if this is a single reaction folder.
        """
        if not os.path.isdir(self.folder_path):
            return []
        
        try:
            # Check if this folder has subfolders with RDF files
            immediate_subdirs = []
            has_rdf_in_root = False
            
            for entry in os.scandir(self.folder_path):
                if entry.is_file() and entry.name.lower().endswith('.rdf'):
                    has_rdf_in_root = True
                elif entry.is_dir():
                    # Check if this subfolder has RDF files
                    has_rdf = any(
                        f.lower().endswith('.rdf')
                        for root, _, files in os.walk(entry.path)
                        for f in files
                    )
                    if has_rdf:
                        immediate_subdirs.append((entry.path, entry.name))
            
            # If root has RDF files, treat as single folder
            # If we have 2+ subfolders with RDF files and no RDF in root, batch process
            if has_rdf_in_root or len(immediate_subdirs) <= 1:
                return []
            
            return sorted(immediate_subdirs, key=lambda x: x[1])
            
        except Exception as e:
            self._emit(f"Warning: Error detecting subfolders: {e}")
            return []
    
    def _batch_process_subfolders(self, subfolders: List[Tuple[str, str]]) -> None:
        """Process multiple reaction type subfolders separately."""
        import time
        
        total_folders = len(subfolders)
        successful = 0
        failed = 0
        results = []
        
        start_time = time.time()
        self._emit(f"\n{'='*70}")
        self._emit(f"BATCH PROCESSING STARTED: {total_folders} reaction type folders detected")
        self._emit(f"{'='*70}\n")
        
        for idx, (subfolder_path, reaction_type) in enumerate(subfolders, 1):
            folder_start_time = time.time()
            progress_pct = int((idx - 1) / total_folders * 100)
            
            self._emit(f"\n{'='*70}")
            self._emit(f"[{progress_pct}%] Folder {idx}/{total_folders}: {reaction_type}")
            self._emit(f"{'='*70}")
            self._emit(f"Path: {subfolder_path}")
            self._emit(f"Progress: {successful} successful, {failed} failed so far\n")
            
            try:
                # Process this subfolder
                success = self._process_single_reaction_type(subfolder_path, reaction_type)
                
                folder_elapsed = time.time() - folder_start_time
                
                if success:
                    successful += 1
                    results.append(f"✓ {reaction_type}")
                    self._emit(f"\n✓ {reaction_type} completed in {folder_elapsed:.1f}s")
                else:
                    failed += 1
                    results.append(f"✗ {reaction_type} (no reactions)")
                    self._emit(f"\n✗ {reaction_type} - no valid reactions found")
            except Exception as e:
                failed += 1
                folder_elapsed = time.time() - folder_start_time
                error_msg = str(e).split('\n')[0][:100]  # First line, truncated
                results.append(f"✗ {reaction_type} (error: {error_msg})")
                self._emit(f"\n✗ {reaction_type} failed after {folder_elapsed:.1f}s")
                self._emit(f"   Error: {error_msg}")
            
            # Estimated time remaining
            if idx < total_folders:
                elapsed_total = time.time() - start_time
                avg_time_per_folder = elapsed_total / idx
                remaining_folders = total_folders - idx
                eta_seconds = avg_time_per_folder * remaining_folders
                eta_minutes = int(eta_seconds / 60)
                eta_seconds_remainder = int(eta_seconds % 60)
                
                self._emit(f"\nEstimated time remaining: {eta_minutes}m {eta_seconds_remainder}s")
                self._emit(f"(Average: {avg_time_per_folder:.1f}s per folder)")
        
        # Summary
        self._emit(f"\n{'='*60}")
        self._emit("BATCH PROCESSING COMPLETE")
        self._emit(f"{'='*60}")
        self._emit(f"Total: {total_folders} | Success: {successful} | Failed: {failed}\n")
        self._emit("Results:")
        for result in results:
            self._emit(f"  {result}")
        
        if self.finished:
            if successful > 0:
                self.finished.emit(
                    True,
                    f"Batch processing complete.\n\n"
                    f"Processed {total_folders} reaction types\n"
                    f"Successful: {successful}\n"
                    f"Failed: {failed}\n\n"
                    f"Output files saved to:\n"
                    f"  - data-processor/reaction_dataset/<ReactionType>.csv"
                )
            else:
                self.finished.emit(
                    False,
                    f"Batch processing failed: No reactions found in any subfolder."
                )
    
    def _process_single_reaction_type(self, folder_path: str, reaction_type: str) -> bool:
        """Process a single reaction type folder.
        
        Returns:
            True if reactions were successfully processed, False otherwise.
        """
        # Temporarily set folder path
        original_folder = self.folder_path
        self.folder_path = folder_path
        
        try:
            # Find RDF files in this subfolder
            self._emit(f"[1/5] Scanning for RDF files...")
            self.rdf_files = self._find_rdf_files()
            
            if not self.rdf_files:
                self._emit(f"      No RDF files found")
                return False
            
            self._emit(f"      Found {len(self.rdf_files)} RDF files")
            
            # Process RDF files
            self._emit(f"[2/5] Parsing RDF files...")
            combined_rdf_map = self._process_rdf_files()
            
            if not combined_rdf_map:
                self._emit(f"      No valid reactions parsed")
                return False
            
            # Count MOL blocks
            rct_mol_count = sum(1 for v in combined_rdf_map.values() if v.get('rct_mol'))
            pro_mol_count = sum(1 for v in combined_rdf_map.values() if v.get('pro_mol'))
            self._emit(f"      Parsed {len(combined_rdf_map)} reactions")
            self._emit(f"      MOL blocks: {rct_mol_count} reactants, {pro_mol_count} products")
            
            # Skip reagent registry lookups; export will be CAS-only
            self._emit(f"[3/5] Preparing condition mappings (CAS-only)...")
            cas_map: Dict[str, Dict[str, str]] = {}
            
            # Create minimal TXT map
            txt_map = self._create_minimal_txt_map(combined_rdf_map)
            self._emit(f"      Mapped {len(txt_map)} reaction conditions")
            
            # Assemble rows
            self._emit(f"[4/5] Assembling reaction data structures...")
            rows = assemble_rows(txt_map, combined_rdf_map, cas_map, txt_preferred=False)
            self._emit(f"      Assembled {len(rows)} reaction records")
            
            if not rows:
                self._emit(f"      No valid reaction records assembled")
                return False
            
            # Set reaction type for all rows
            for row in rows:
                row['ReactionType'] = reaction_type

            rows, duplicate_count = self._deduplicate_rows(rows, combined_rdf_map)
            if duplicate_count:
                self._emit(
                    f"      Removed {duplicate_count} exact duplicate reaction records"
                )
            
            # Calculate output paths for this reaction type
            self._emit(f"[5/5] Generating output files...")
            repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            dataset_dir = os.path.join(repo_root, "data-processor", "reaction_dataset")
            os.makedirs(dataset_dir, exist_ok=True)
            
            # Use reaction_type as filename (sanitized)
            import re as _re
            safe_name = _re.sub(r'[^A-Za-z0-9_-]+', '', _re.sub(r'\s+', '_', reaction_type))
            
            output_csv = os.path.join(dataset_dir, f"{safe_name}.csv")
            
            # Generate outputs
            generator = ReactionMarkdownGenerator()
            
            source_name = f"RDF_{reaction_type}"
            
            self._emit(f"      Generating CSV export...")
            generator.generate_csv_export(rows, output_csv, source_name, rdf_map=combined_rdf_map)
            
            # Show file sizes for verification
            try:
                csv_size = os.path.getsize(output_csv) / 1024  # KB
                self._emit(f"\n      ✓ Successfully saved {len(rows)} reactions:")
                self._emit(f"        CSV:  {output_csv} ({csv_size:.1f} KB)")
            except Exception:
                self._emit(f"\n      ✓ Saved {len(rows)} reactions")
                self._emit(f"        CSV: {output_csv}")
            
            return True
            
        finally:
            # Restore original folder path
            self.folder_path = original_folder
    
    def _find_rdf_files(self) -> List[str]:
        """Find all RDF files in the specified folder and its subfolders (recursive)"""
        rdf_files = []
        try:
            # Walk through all subdirectories
            for root, dirs, files in os.walk(self.folder_path):
                for file in files:
                    if file.lower().endswith('.rdf'):
                        full_path = os.path.join(root, file)
                        rdf_files.append(full_path)
        except Exception as e:
            raise RuntimeError(f"Error scanning folder: {e}")
        
        return sorted(rdf_files)

    def _create_minimal_txt_map(self, rdf_map: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """Create a minimal TXT map from RDF data (since we only have RDF)"""
        txt_map: Dict[str, Dict[str, Any]] = {}
        # Lightweight regex patterns (mirrors logic in process_reactions but simplified)
        import re, math
        re_time = re.compile(r"(?P<val>\d+(?:\.\d+)?)\s*(?P<unit>h|hr|hrs|hour|hours|min|mins|minute|minutes|d|day|days)\b", re.I)
        re_temp = re.compile(r"(?P<val>-?\d+(?:\.\d+)?)\s*[^A-Za-z0-9]{0,3}C\b")
        re_rt = re.compile(r"\brt\b|room temperature", re.I)

        for rid, rdf_data in rdf_map.items():
            notes = rdf_data.get('notes') or []
            all_condition_lines: List[str] = []
            # Use notes lines as condition lines source (SciFinder often stores experimental snippets here)
            for ln in notes:
                if isinstance(ln, str) and ln.strip():
                    all_condition_lines.append(ln.strip())

            # Aggregate time and temperature heuristically from notes
            total_h = 0.0
            max_c = -math.inf
            had_rt = False
            for ln in all_condition_lines:
                # Skip DOI-like lines
                if re.search(r"\b10\.\d{4,9}/", ln):
                    continue
                # time
                for m in re_time.finditer(ln):
                    try:
                        val = float(m.group('val'))
                    except ValueError:
                        continue
                    unit = m.group('unit').lower()
                    if unit.startswith('min'):
                        total_h += val / 60.0
                    elif unit in ('d', 'day', 'days'):
                        total_h += val * 24.0
                    else:
                        total_h += val
                if re.search(r"\bovernight\b", ln, re.I):
                    total_h += 16.0
                # temperature
                for m in re_temp.finditer(ln):
                    try:
                        valc = float(m.group('val'))
                    except ValueError:
                        continue
                    if valc > max_c:
                        max_c = valc
                if re_rt.search(ln):
                    had_rt = True

            temperature_c = max_c if max_c != -math.inf else (25.0 if had_rt else None)
            time_h = round(total_h, 3) if total_h > 0 else None

            txt_map[rid] = {
                'original_text': [],
                'all_condition_lines': all_condition_lines,
                'time_h': time_h,
                'temperature_c': round(temperature_c, 1) if temperature_c is not None else None,
                'title': rdf_data.get('title', ''),
                'authors': rdf_data.get('authors', ''),
                'citation': rdf_data.get('citation', ''),
                'doi': '',
                'reagents': [],
                'catalysts': [],
                'solvents': [],
                'txt_yield': None,
            }
        
        return txt_map

    def _extract_temp_time_from_md(self, md_path: str) -> Dict[str, Dict[str, Optional[float]]]:
        """Parse a markdown file to map CAS Reaction Number -> {temperature_c, time_h}.

        Heuristics:
        - Within each block following a line like "CAS Reaction Number: <ID>",
          accumulate time across all occurrences and take the max temperature.
        - Recognize units h/hr/hrs/hour/min/mins/minute/day/days; minutes converted to hours; days*24.
        - Recognize temperatures like "80 C" or "80 °C"; recognize 'rt'/'room temperature' as 25 °C when no numeric temp.
        - Ignore 'reflux' for temperature.
        """
        result: Dict[str, Dict[str, Optional[float]]] = {}
        if not os.path.exists(md_path):
            return result

        import re, math
        re_time = re.compile(r"(?<![A-Za-z0-9])(\d+(?:\.\d+)?)\s*(h|hr|hrs|hour|hours|min|mins|minute|minutes|d|day|days)(?![A-Za-z0-9])", re.I)
        re_temp_c = re.compile(r"(-?\d+(?:\.\d+)?)\s*[^A-Za-z0-9]{0,3}C\b")
        re_rt = re.compile(r"\brt\b|room temperature", re.I)
        re_rid = re.compile(r"^\s*CAS Reaction Number:\s*(\S+)\s*$", re.I)

        current_id: Optional[str] = None
        agg_time: float = 0.0
        agg_max_c: float = -math.inf
        had_rt: bool = False

        def _flush():
            nonlocal current_id, agg_time, agg_max_c, had_rt
            if current_id:
                temp_c: Optional[float]
                if agg_max_c != -math.inf:
                    temp_c = round(agg_max_c, 1)
                elif had_rt:
                    temp_c = 25.0
                else:
                    temp_c = None
                time_h = round(agg_time, 3) if agg_time > 0 else None
                result[current_id] = {"temperature_c": temp_c, "time_h": time_h}
            current_id = None
            agg_time = 0.0
            agg_max_c = -math.inf
            had_rt = False

        try:
            with open(md_path, "r", encoding="utf-8", errors="ignore") as f:
                for raw in f:
                    line = raw.rstrip("\n")
                    m_id = re_rid.match(line)
                    if m_id:
                        # flush previous block
                        _flush()
                        current_id = m_id.group(1).strip()
                        continue
                    if not current_id:
                        continue
                    # Accumulate within current block
                    for m in re_time.finditer(line):
                        try:
                            val = float(m.group(1))
                        except Exception:
                            continue
                        unit = (m.group(2) or "").lower()
                        if unit.startswith("min"):
                            agg_time += val / 60.0
                        elif unit in ("d", "day", "days"):
                            agg_time += val * 24.0
                        else:
                            agg_time += val
                    mtemp = re_temp_c.findall(line)
                    for t in mtemp:
                        try:
                            v = float(t)
                        except Exception:
                            continue
                        if v > agg_max_c:
                            agg_max_c = v
                    if re_rt.search(line):
                        had_rt = True
            # flush last block
            _flush()
        except Exception:
            return result

        return result

    def _process_rdf_files(self) -> Dict[str, Dict[str, Any]]:
        """Process RDF files, losslessly merging records with the same ID."""
        combined_rdf_map: Dict[str, Dict[str, Any]] = {}
        total_files = len(self.rdf_files)
        
        for i, rdf_file in enumerate(self.rdf_files, 1):
            filename = os.path.basename(rdf_file)
            progress_pct = int(i / total_files * 100)
            
            # Show progress every file (for small batches) or every 10% (for large batches)
            if total_files <= 20 or i % max(1, total_files // 10) == 0 or i == total_files:
                self._emit(f"      [{progress_pct}%] File {i}/{total_files}: {filename}")
            
            try:
                # Parse individual RDF file
                rdf_map = parse_rdf(rdf_file)
                relative_source = os.path.relpath(rdf_file, self.folder_path)
                reaction_type = self._reaction_type_for_source(rdf_file)
                added = 0
                merged_count = 0
                for rid, data in rdf_map.items():
                    data['source_file'] = relative_source
                    data['source_files'] = [relative_source]
                    data['reaction_types'] = [reaction_type]
                    data['duplicate_reaction_ids'] = [rid]
                    if rid in combined_rdf_map:
                        combined_rdf_map[rid] = merge_rdf_records(
                            combined_rdf_map[rid], data,
                        )
                        merged_count += 1
                    else:
                        combined_rdf_map[rid] = data
                        added += 1
                
                # Only show details if reactions were found or errors occurred
                if len(rdf_map) > 0:
                    msg_tail = f" ({added} new"
                    if merged_count:
                        msg_tail += f", {merged_count} merged duplicates"
                    msg_tail += ")"
                    if total_files <= 20:  # Show details for small batches
                        self._emit(f"           {len(rdf_map)} reactions{msg_tail}")
            except Exception as e:
                self._emit(f"           Error: {str(e)[:80]}")
                continue
        
        self._emit(f"      Total: {len(combined_rdf_map)} unique reactions from {total_files} files")
        return combined_rdf_map

    def _reaction_type_for_source(self, rdf_file: str) -> str:
        """Return a stable category label for one source RDF file."""
        parts = Path(os.path.normpath(rdf_file)).parts
        lowered = [part.lower() for part in parts]
        if "original_dataset" in lowered:
            index = lowered.index("original_dataset")
            if index + 1 < len(parts):
                return parts[index + 1]
        relative = Path(os.path.relpath(rdf_file, self.folder_path))
        if self.combine_all and len(relative.parts) > 1:
            return relative.parts[0]
        return os.path.basename(os.path.normpath(self.folder_path)) or "Other"

    @staticmethod
    def _normalized_smiles(smiles: str) -> str:
        """Canonicalize a whole reaction side, including component ordering."""
        if not smiles:
            return ""
        if RDKIT_AVAILABLE:
            try:
                molecule = Chem.MolFromSmiles(smiles)
                if molecule is not None:
                    return Chem.MolToSmiles(molecule, isomericSmiles=True)
            except Exception:
                pass
        return "".join(smiles.split())

    def _deduplication_key(
        self,
        row: Dict[str, Any],
        rdf_data: Dict[str, Any],
    ) -> str:
        """Build an order-invariant key for an exact chemistry precedent."""
        reactant_smiles = self._normalized_smiles(row.get("ReactantSMILES", ""))
        product_smiles = self._normalized_smiles(row.get("ProductSMILES", ""))
        has_structure = bool(reactant_smiles or product_smiles)
        has_identified_sides = bool(
            rdf_data.get("rct_cas") and rdf_data.get("pro_cas")
        )
        if not has_structure and not has_identified_sides:
            # Do not collapse poorly specified records merely because their
            # condition fields are empty or similar.
            return f"reaction-id:{row.get('ReactionID', '')}"

        def normalized_list(field: str) -> List[str]:
            return sorted(
                {" ".join(str(value).split()) for value in rdf_data.get(field, []) if value},
                key=str.casefold,
            )

        payload = {
            "reactants": reactant_smiles,
            "products": product_smiles,
            "rct_cas": normalized_list("rct_cas"),
            "pro_cas": normalized_list("pro_cas"),
            "rgt_cas": normalized_list("rgt_cas"),
            "cat_cas": normalized_list("cat_cas"),
            "sol_cas": normalized_list("sol_cas"),
            "rct_amd": normalized_list("rct_amd"),
            "pro_amd": normalized_list("pro_amd"),
            "rgt_amd": normalized_list("rgt_amd"),
            "cat_amd": normalized_list("cat_amd"),
            "sol_amd": normalized_list("sol_amd"),
            "product_yields": sorted(
                (str(index), value)
                for index, value in (rdf_data.get("pro_yields") or {}).items()
            ),
            "temperature_c": row.get("Temperature_C", ""),
            "time_h": row.get("Time_h", ""),
            "experimental_procedure": normalized_list("exp_proc"),
            "stages": " ".join(str(rdf_data.get("stages") or "").split()),
            "steps": " ".join(str(rdf_data.get("steps") or "").split()),
            "notes": normalized_list("notes"),
            "title": " ".join(str(rdf_data.get("title") or "").split()),
            "authors": " ".join(str(rdf_data.get("authors") or "").split()),
            "citation": " ".join(str(rdf_data.get("citation") or "").split()),
        }
        return _json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))

    def _deduplicate_rows(
        self,
        rows: List[Dict[str, Any]],
        rdf_map: Dict[str, Dict[str, Any]],
    ) -> Tuple[List[Dict[str, Any]], int]:
        """Remove exact duplicate precedents and merge their provenance."""
        unique_rows: List[Dict[str, Any]] = []
        primary_by_key: Dict[str, Dict[str, Any]] = {}
        duplicate_ids: List[str] = []

        for row in rows:
            reaction_id = str(row.get("ReactionID") or "")
            rdf_data = rdf_map.get(reaction_id, {})
            key = self._deduplication_key(row, rdf_data)
            primary_row = primary_by_key.get(key)
            if primary_row is None:
                primary_by_key[key] = row
                unique_rows.append(row)
                continue

            primary_id = str(primary_row.get("ReactionID") or "")
            rdf_map[primary_id] = merge_rdf_records(
                rdf_map.get(primary_id, {}), rdf_data,
            )
            reaction_types = rdf_map[primary_id].get("reaction_types") or []
            primary_row["ReactionType"] = " | ".join(
                sorted(set(reaction_types), key=str.casefold)
            )
            duplicate_ids.append(reaction_id)

        for reaction_id in duplicate_ids:
            rdf_map.pop(reaction_id, None)
        return unique_rows, len(duplicate_ids)

    def _generate_outputs(self, rows: List[Dict[str, Any]], rdf_map: Optional[Dict[str, Dict[str, Any]]] = None) -> None:
        """Generate the CSV output without creating a Markdown report."""
        generator = ReactionMarkdownGenerator()
        source_name = f"RDF_Folder_{os.path.basename(self.folder_path)}"

        self._emit("Generating CSV export...")
        generator.generate_csv_export(rows, self.output_csv_path, source_name, rdf_map=rdf_map)
        
        # NOTE: JSONL export is no longer used; CSV is saved directly to data-processor/reaction_dataset/.

    @Slot() if Slot else (lambda f: f)
    def run(self):
        """Main processing function"""
        try:
            # Check if we should batch-process multiple subfolders
            subfolders = self._get_reaction_type_subfolders()
            
            if len(subfolders) > 1 and not self.combine_all:
                self._emit(f"Detected {len(subfolders)} reaction type subfolders. Processing each separately...")
                self._batch_process_subfolders(subfolders)
                return
            
            # Single folder processing (original behavior)
            # Find all RDF files
            self._emit("Scanning folder for RDF files...")
            self.rdf_files = self._find_rdf_files()
            
            if not self.rdf_files:
                if self.finished:
                    self.finished.emit(False, "No RDF files found in the selected folder.")
                return
            
            self._emit(f"Found {len(self.rdf_files)} RDF files.")
            
            # Process all RDF files and combine them
            self._emit("Processing RDF files...")
            combined_rdf_map = self._process_rdf_files()
            
            if not combined_rdf_map:
                if self.finished:
                    self.finished.emit(False, "No valid reactions found in RDF files.")
                return
            
            # Count MOL blocks for diagnostics
            rct_mol_count = sum(1 for v in combined_rdf_map.values() if v.get('rct_mol'))
            pro_mol_count = sum(1 for v in combined_rdf_map.values() if v.get('pro_mol'))
            self._emit(f"RDF parsed. Reactions with reactant MOL blocks: {rct_mol_count}; with product MOL blocks: {pro_mol_count}")
            self._emit(f"RDKit available: {RDKIT_AVAILABLE}")
            
            # Skip reagent registry lookups; export will be CAS-only
            self._emit("Skipping reagent registry lookups (CAS-only export).")
            cas_map: Dict[str, Dict[str, str]] = {}
            
            # Create minimal TXT map (since we only have RDF)
            self._emit("Creating minimal TXT mapping...")
            txt_map = self._create_minimal_txt_map(combined_rdf_map)
            
            # Assemble rows using the same pipeline as original
            self._emit("Assembling reaction rows...")
            rows = assemble_rows(txt_map, combined_rdf_map, cas_map, txt_preferred=False)
            self._emit(f"Assembled {len(rows)} rows")

            # Override Temperature_C and Time_h from dataset/temp_time.md when available
            here = os.path.dirname(os.path.abspath(__file__))
            md_path = os.path.join(here, 'dataset', 'temp_time.md')
            md_map = self._extract_temp_time_from_md(md_path)
            if md_map:
                overridden = 0
                for row in rows:
                    rid = row.get('ReactionID')
                    if not rid:
                        continue
                    mt = md_map.get(rid)
                    if not mt:
                        continue
                    t_c = mt.get('temperature_c')
                    t_h = mt.get('time_h')
                    if t_c is not None:
                        row['Temperature_C'] = t_c
                    if t_h is not None:
                        row['Time_h'] = t_h
                    if (t_c is not None) or (t_h is not None):
                        overridden += 1
                self._emit(f"Applied temp/time overrides from temp_time.md for {overridden} reactions.")

            # Folder-derived types are annotations; preserve all categories when
            # an exact duplicate appeared in more than one source folder.
            for row in rows:
                rdf_data = combined_rdf_map.get(str(row.get('ReactionID') or ''), {})
                reaction_types = rdf_data.get('reaction_types') or []
                if reaction_types:
                    row['ReactionType'] = " | ".join(
                        sorted(set(reaction_types), key=str.casefold)
                    )

            rows, duplicate_count = self._deduplicate_rows(rows, combined_rdf_map)
            if duplicate_count:
                self._emit(
                    f"Removed {duplicate_count} exact duplicate reaction records "
                    "after merging their IDs and source provenance."
                )
            
            # Count rows with SMILES for diagnostics
            smi_rows = sum(1 for r in rows if (r.get('ReactantSMILES') or r.get('ProductSMILES')))
            self._emit(f"Rows with SMILES: {smi_rows} / {len(rows)}")
            
            if smi_rows == 0:
                if not RDKIT_AVAILABLE:
                    self._emit("Note: RDKit is not available; SMILES generation from MOL blocks is disabled.")
                elif (rct_mol_count + pro_mol_count) == 0:
                    self._emit("Note: No MOL/CTAB blocks found in RDF content; SMILES cannot be generated.")
                else:
                    self._emit("Warning: MOL blocks found and RDKit available, but SMILES are empty. MOL data may be malformed.")
            
            # Generate outputs
            self._generate_outputs(rows, combined_rdf_map)
            
            if self.finished:
                self.finished.emit(
                    True,
                    (
                        f"Successfully processed {len(self.rdf_files)} RDF files with {len(rows)} reactions.\n\n"
                        f"CSV (reaction dataset): {self.output_csv_path}"
                    ),
                )
                
        except Exception as e:
            msg = f"Error: {e}\n\n{traceback.format_exc()}"
            if self.finished:
                self.finished.emit(False, msg)


class RDFProcessorWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SciFinder RDF Processor")
        self.resize(700, 500)
        
        # Input controls
        self.folder_edit = QtWidgets.QLineEdit()
        self.btn_folder = QtWidgets.QPushButton("Browse Folder...")
        self.combine_checkbox = QtWidgets.QCheckBox(
            "Combine all subfolders into one deduplicated CSV"
        )
        self.combine_checkbox.setChecked(False)
        
        # File list display
        self.file_list = QtWidgets.QListWidget()
        self.file_count_label = QtWidgets.QLabel("No folder selected")
        
        # Control buttons
        self.btn_run = QtWidgets.QPushButton("Process RDF Files")
        self.btn_quit = QtWidgets.QPushButton("Quit")
        
        # Log output
        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(150)
        
        # Setup layout
        self._setup_layout()
        
        # Connect signals
        self.btn_folder.clicked.connect(self.pick_folder)
        self.btn_run.clicked.connect(self.run_processing)
        self.btn_quit.clicked.connect(self.close)
        
        # Runtime state
        self.thread = None
        self.worker = None
        self.rdf_files = []
        
        # Initialize button states
        self.btn_run.setEnabled(False)

    def _setup_layout(self):
        """Setup the GUI layout"""
        # Main layout
        layout = QtWidgets.QVBoxLayout(self)
        
        # Title
        title = QtWidgets.QLabel("SciFinder RDF File Processor")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)
        
        # Form layout for inputs
        form = QtWidgets.QFormLayout()
        
        # Folder selection
        folder_box = QtWidgets.QHBoxLayout()
        folder_box.addWidget(self.folder_edit)
        folder_box.addWidget(self.btn_folder)
        folder_hint = QtWidgets.QLabel("(includes all subfolders)")
        folder_hint.setStyleSheet("font-style: italic; color: #666; font-size: 9px;")
        form.addRow("RDF Folder:", folder_box)
        form.addRow("", folder_hint)
        form.addRow("Output mode:", self.combine_checkbox)
        
        # Add note about file locations
        note_label = QtWidgets.QLabel(
            "Note: RDF files are scanned recursively. Enable combined mode to create\n"
            "      one cross-category file with exact duplicates merged.\n"
            "      CSV -> data-processor/reaction_dataset/{category}.csv"
        )
        note_label.setStyleSheet("font-style: italic; color: #666; font-size: 10px;")
        form.addRow("", note_label)
        
        
        layout.addLayout(form)
        
        # File list section
        file_group = QtWidgets.QGroupBox("RDF Files Found")
        file_layout = QtWidgets.QVBoxLayout(file_group)
        file_layout.addWidget(self.file_count_label)
        file_layout.addWidget(self.file_list)
        layout.addWidget(file_group)
        
        # Control buttons
        button_layout = QtWidgets.QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(self.btn_run)
        button_layout.addWidget(self.btn_quit)
        layout.addLayout(button_layout)
        
        # Log section
        log_group = QtWidgets.QGroupBox("Processing Log")
        log_layout = QtWidgets.QVBoxLayout(log_group)
        log_layout.addWidget(self.log)
        layout.addWidget(log_group)

    def log_msg(self, text: str):
        """Add a message to the log"""
        self.log.appendPlainText(text)

    def pick_folder(self):
        """Select folder containing RDF files"""
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self, 
            "Select folder with RDF files", 
            os.getcwd(), 
            options=QtWidgets.QFileDialog.Option.ShowDirsOnly
        )
        if path:
            self.folder_edit.setText(path)
            self._update_file_list()
            

    def _update_file_list(self):
        """Update the list of RDF files found in the selected folder and subfolders"""
        folder_path = self.folder_edit.text().strip()
        self.file_list.clear()
        self.rdf_files = []
        
        if not folder_path or not os.path.isdir(folder_path):
            self.file_count_label.setText("No valid folder selected")
            self.btn_run.setEnabled(False)
            return
        
        try:
            # Find RDF files recursively in all subfolders
            for root, dirs, files in os.walk(folder_path):
                for file in files:
                    if file.lower().endswith('.rdf'):
                        full_path = os.path.join(root, file)
                        self.rdf_files.append(full_path)
                        # Display relative path from selected folder for better clarity
                        rel_path = os.path.relpath(full_path, folder_path)
                        self.file_list.addItem(rel_path)
            
            # Update UI
            count = len(self.rdf_files)
            if count == 0:
                self.file_count_label.setText("No RDF files found in this folder or subfolders")
                self.btn_run.setEnabled(False)
            else:
                self.file_count_label.setText(f"Found {count} RDF file{'s' if count != 1 else ''} (including subfolders)")
                self.btn_run.setEnabled(True)
                
        except Exception as e:
            self.file_count_label.setText(f"Error reading folder: {e}")
            self.btn_run.setEnabled(False)

    def validate_inputs(self) -> Optional[str]:
        """Validate user inputs"""
        folder = self.folder_edit.text().strip()
        if not folder or not os.path.isdir(folder):
            return "Please select a valid folder containing RDF files."
        
        if not self.rdf_files:
            return "No RDF files found in the selected folder."
        
        return None

    def run_processing(self):
        """Start the RDF processing"""
        err = self.validate_inputs()
        if err:
            QtWidgets.QMessageBox.warning(self, "Invalid Input", err)
            return
        
        # Disable UI during processing
        combine_all = self.combine_checkbox.isChecked()
        self.setEnabled(False)
        self.log.clear()
        self.log_msg("Starting RDF processing...")
        
        # Calculate output paths
        # Name CSV using folder hierarchy, e.g., C_N_Coupling/2020-2022 -> C_N_Coupling.csv
        try:
            folder_path = self.folder_edit.text().strip()
            norm_folder = os.path.normpath(folder_path)
            import re as _re
            def _safe(s: str) -> str:
                s = _re.sub(r"\s+", "", s or "")  # remove spaces
                s = _re.sub(r"[^A-Za-z0-9_-]+", "", s)  # keep only safe chars
                return s
            
            # Determine the main category folder name
            # If the selected folder is under original_dataset/, use that as category
            folder_parts = norm_folder.split(os.sep)
            category_name = os.path.basename(norm_folder)
            
            # Check if this is a category folder (like C_N_Coupling, Suzuki, etc.)
            # by looking for "original_dataset" in path and taking the next folder
            if "original_dataset" in folder_parts:
                idx = folder_parts.index("original_dataset")
                if idx + 1 < len(folder_parts):
                    category_name = folder_parts[idx + 1]
            
            # CSV: Save directly to data-processor/reaction_dataset/
            # Use category name (all subfolders are automatically included)
            repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            dataset_dir = os.path.join(repo_root, "data-processor", "reaction_dataset")
            os.makedirs(dataset_dir, exist_ok=True)
            
            final_name = _safe(category_name) or "dataset"
            if combine_all:
                final_name += "_combined_deduplicated"
            output_csv = os.path.join(dataset_dir, final_name + ".csv")

        except Exception:
            # Fallback: use default paths
            repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            dataset_dir = os.path.join(repo_root, "data-processor", "reaction_dataset")
            os.makedirs(dataset_dir, exist_ok=True)
            output_csv = os.path.join(dataset_dir, "dataset.csv")
        
        # Create worker and thread
        self.worker = RDFWorker(
            folder_path=self.folder_edit.text().strip(),
            output_md_path="",
            output_csv_path=output_csv,
            combine_all=combine_all,
        )
        
        self.thread = QtCore.QThread(self)
        self.worker.moveToThread(self.thread)
        
        # Connect signals
        self.thread.started.connect(self.worker.run)
        
        sig = getattr(self.worker, 'finished', None)
        if sig:
            sig.connect(self.on_finished)
            sig.connect(self.thread.quit)
            sig.connect(self.worker.deleteLater)
        
        prog = getattr(self.worker, 'progress', None)
        if prog:
            prog.connect(self.log_msg)
        
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(lambda: self.setEnabled(True))
        self.thread.finished.connect(lambda: setattr(self, 'worker', None))
        self.thread.finished.connect(lambda: setattr(self, 'thread', None))
        
        # Start processing
        self.thread.start()

    def on_finished(self, success: bool, message: str):
        """Handle processing completion"""
        self.setEnabled(True)
        self.log_msg(message)
        
        if success:
            QtWidgets.QMessageBox.information(self, "Processing Complete", message)
        else:
            QtWidgets.QMessageBox.critical(self, "Processing Error", message)


def main():
    """Main application entry point"""
    if hasattr(QtWidgets, 'QApplication'):
        try:
            QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_EnableHighDpiScaling, True)
            QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_UseHighDpiPixmaps, True)
        except Exception:
            pass
    
    app = QtWidgets.QApplication(sys.argv)
    window = RDFProcessorWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()






