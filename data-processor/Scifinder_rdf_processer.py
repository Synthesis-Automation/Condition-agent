#utf-8
#!/usr/bin/env python3
"""
Simple Qt6 GUI wrapper for processing RDF files only.
Lets the user pick a folder containing RDF files and processes all RDF files in the folder.
Works with PySide6 (preferred) or PyQt6 if installed.
"""
from __future__ import annotations

import os
import sys
import traceback
import tempfile
import re
import urllib.request
import urllib.error
import urllib.parse
from typing import List, Optional, Dict, Any, Set, Tuple

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
    from process_reactions import parse_rdf, assemble_rows, load_cas_maps, build_condkey, build_condsig, build_famsig
except Exception as e:
    print(f"Error: Cannot import processing helpers: {e}")
    sys.exit(1)

try:
    # We import but will prefer our taxonomy-aware generator for this GUI tool.
    from reaction_markdown_generator import ReactionMarkdownGenerator as _ExternalReactionMarkdownGenerator  # type: ignore
except Exception:
    _ExternalReactionMarkdownGenerator = None  # type: ignore


# ---------------------------- Taxonomy integration ----------------------------

import json as _json

PUBCHEM_PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound"
PUBCHEM_USER_AGENT = "ConditionAgent/1.0 (SciFinder RDF Processor)"
PUBCHEM_TIMEOUT = 10
CAS_NUMBER_RE = re.compile(r'\d{2,7}-\d{2}-\d')


REGISTRY_FILE_EXCLUDES = {
    "not_determined_reagents.json",
}


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
    "metal_precursor": "METAL_PRECURSOR",
    "oxidant": "OXIDANT",
    "preformed_metal_catalyst": "PREFORMED_METAL_CATALYST",
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
    "metal_precursor": "catalyst",
    "oxidant": "oxidant",
    "preformed_metal_catalyst": "catalyst",
    "reductant": "reductant",
    "solvent": "solvent",
    "other_reagent": "other",
    "organo_catalyst": "catalyst",
    "enzyme": "catalyst",
}

CATALYST_ROLES: Set[str] = {
    "metal_precursor",
    "preformed_metal_catalyst",
}

class _TaxonomyIndex:
    """Load the reagent registry and expose lookup utilities compatible with the
    previous taxonomy-based workflow."""

    def __init__(self, base_dir: str) -> None:
        raw_base = os.path.abspath(base_dir)
        if not os.path.isdir(raw_base) and os.path.basename(raw_base) == 'compound_taxonomy':
            candidate = os.path.join(os.path.dirname(raw_base), 'reagents')
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
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = _json.load(fh)
        except Exception:
            return []
        return data if isinstance(data, list) else []

    def _iter_registry_files(self) -> List[str]:
        files: List[str] = []
        try:
            for entry in os.scandir(self.base_dir):
                if not entry.is_file():
                    continue
                name = entry.name
                if not name.lower().endswith('.json'):
                    continue
                if name in REGISTRY_FILE_EXCLUDES:
                    continue
                files.append(name)
        except (FileNotFoundError, NotADirectoryError, OSError):
            return []
        return sorted(files)

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
                abbr_field = entry.get("abbreviation") or []
                if isinstance(abbr_field, str):
                    abbr_list = [abbr_field.strip()] if abbr_field.strip() else []
                elif isinstance(abbr_field, list):
                    abbr_list = [str(a).strip() for a in abbr_field if a]
                else:
                    abbr_list = []
                abbr = abbr_list[0] if abbr_list else ""
                alias_field = entry.get("aliases") or []
                if isinstance(alias_field, str):
                    alias_list = [alias_field.strip()] if alias_field.strip() else []
                elif isinstance(alias_field, list):
                    alias_list = [str(a).strip() for a in alias_field if a]
                else:
                    alias_list = []
                extra_synonyms = entry.get("synonyms")
                if isinstance(extra_synonyms, list):
                    extra_list = [str(a).strip() for a in extra_synonyms if a]
                elif isinstance(extra_synonyms, str):
                    extra_list = [extra_synonyms.strip()] if extra_synonyms.strip() else []
                else:
                    extra_list = []
                role_map = entry.get("roles")
                if not isinstance(role_map, dict):
                    role_map = {}
                if not role_map:
                    continue
                synonyms = abbr_list[1:] + alias_list + extra_list
                for role_name, role_payload_raw in role_map.items():
                    role_key = str(role_name).strip()
                    if not role_key:
                        continue
                    role_payload = role_payload_raw if isinstance(role_payload_raw, dict) else {}
                    role_lookup = role_key.lower()
                    role_code = REGISTRY_ROLE_CODE_MAP.get(role_lookup, role_key.upper())
                    category_hint = REGISTRY_CATEGORY_HINT.get(role_lookup, role_lookup)
                    families = role_payload.get("families") or []
                    family_id = families[0] if families else None
                    metal_val = role_payload.get("metal")
                    if isinstance(metal_val, list):
                        metal = next((str(m).strip() for m in metal_val if m), "")
                    elif metal_val is None:
                        metal = ""
                    else:
                        metal = str(metal_val).strip()
                    generic_core = metal if role_lookup in CATALYST_ROLES else ""
                    self._index_member(
                        cas=cas,
                        name=name,
                        abbr=abbr,
                        synonyms=synonyms,
                        role=role_code,
                        category_hint=category_hint,
                        generic_core=generic_core,
                        family_id=family_id,
                        role_payload=role_payload,
                    )
                    if family_id and metal:
                        self.family_metal.setdefault(family_id, metal)

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

    def _pair_to_obj(self, item: str):
        if "|" in item:
            name, cas = item.split("|", 1)
            return {"name": name.strip(), "cas": cas.strip()}
        return {"name": item.strip(), "cas": ""}

    def _join_names(self, arr):
        if not arr:
            return ""
        out = []
        for it in arr:
            if isinstance(it, dict):
                nm = (it.get("name") or "").strip()
            elif isinstance(it, str):
                nm = it.split("|", 1)[0].strip()
            else:
                nm = str(it)
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
        for src, dest in (("‑", "-"), ("–", "-"), ("—", "-"), ("−", "-"), ("·", "")):
            cleaned = cleaned.replace(src, dest)
        cleaned = cleaned.replace(" ", "")
        return cleaned.strip('-')

    @classmethod
    def _cleanup_preformed_segment(cls, segment: str) -> str:
        if not segment:
            return ""
        segment = segment.replace("‑", "-").replace("–", "-").replace("—", "-")
        segment = segment.replace("·", " ")
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
        normalized = text.replace("‑", "-").replace("–", "-").replace("—", "-")
        for content in re.findall(r'\(([^()]+)\)', normalized):
            candidate = cls._cleanup_preformed_segment(content)
            if cls._is_ligand_candidate(candidate):
                return candidate
        plus_normalized = normalized.replace('·', '+').replace('/', '+')
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
            if not rec or rec.get('Role') != 'PREFORMED_METAL_CATALYST':
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
                    cat_path = os.path.join(self.taxonomy.base_dir, "taxonomy_catalysts_precursor.json")
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

            for role_key, label in fallback_roles:
                for idx, item in enumerate(reag_list or []):
                    role = (role_list[idx] if idx < len(role_list) else "").upper()
                    if role != role_key:
                        continue
                    obj = item if isinstance(item, dict) else self._pair_to_obj(item)
                    token = _token_from_item(obj)
                    if token:
                        return f"{label}: {token}"

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

                    disp_core = self._compute_condition_core(row, full_system_list)

                    T = row.get("Temperature_C", "")
                    t = row.get("Time_h", "")
                    y = row.get("Yield_%", "")

                    reag_out = []
                    for i, item in enumerate(reag_list):
                        obj = self._pair_to_obj(item)
                        role = (role_list[i] if i < len(role_list) else "").upper() or "ADDITIVE"
                        seg = obj.get("name") or obj.get("cas") or "?"
                        if obj.get("cas"):
                            seg += f" ({obj['cas']})"
                        reag_out.append(f"{seg} [{role}]")

                    solv_out = []
                    for item in solv_list:
                        obj = self._pair_to_obj(item)
                        seg = obj.get("name") or obj.get("cas") or "?"
                        if obj.get("cas"):
                            seg += f" ({obj['cas']})"
                        solv_out.append(seg)

                    r_smi = row.get("ReactantSMILES", "")
                    p_smi = row.get("ProductSMILES", "")

                    f.write(f"## Reaction {rid}\n\n")
                    if rtype:
                        f.write(f"- Type: {rtype}\n")
                    if disp_core:
                        f.write(f"- Condition Core: {disp_core}\n")
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

    def generate_jsonl_export(self, rows, output_path: str, source_folder: str):
        """Generate JSONL export with precomputed normalization and features."""
        
        # Import chemtools for preprocessing
        try:
            import sys
            import os
            repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if repo_root not in sys.path:
                sys.path.insert(0, repo_root)
            from chemtools.smiles import normalize_reaction
            from chemtools import router, reaction_similarity as rs
            from chemtools.featurizers import molecular as feat_molecular
            chemtools_available = True
            drfp_available = rs.drfp_available()
        except Exception as e:
            print(f"Warning: chemtools not available for preprocessing: {e}")
            chemtools_available = False
            drfp_available = False
        
        out_lines = []
        precompute_stats = {"success": 0, "failed": 0, "skipped": 0, "drfp_computed": 0}
        family_sources = {"scifinder": 0, "scifinder_partial": 0, "scifinder_unmapped": 0, "unknown": 0}
        
        # Collect DRFP fingerprints for binary storage
        drfp_fingerprints = []  # List of numpy arrays
        drfp_reaction_ids = []  # Corresponding reaction IDs
        
        for idx, row in enumerate(rows, 1):
            rid = row.get("ReactionID", "")
            rtype = row.get("ReactionType", "")
            full_system_list = self._collect_full_catalytic_system(row)
            condition_core = self._compute_condition_core(row, full_system_list)

            # Process catalytic system (same format as reagents)
            catalytic_system = [self._pair_to_obj(x) for x in full_system_list]

            reag_list = self._safe_json_list(row.get("Reagent", "[]"))
            role_list = self._safe_json_list(row.get("ReagentRole", "[]"))
            reagents = []
            for i, item in enumerate(reag_list):
                obj = self._pair_to_obj(item)
                role = (role_list[i] if i < len(role_list) else "").upper() or "ADDITIVE"
                obj["role"] = role
                reagents.append(obj)

            solv_list = self._safe_json_list(row.get("Solvent", "[]"))
            solvents = [self._pair_to_obj(x) for x in solv_list]

            def _num(x):
                try:
                    return float(x)
                except Exception:
                    return None
            conditions = {
                "temperature_c": _num(row.get("Temperature_C")),
                "time_h": _num(row.get("Time_h")),
                "yield_pct": _num(row.get("Yield_%")),
            }

            reactants_smi = row.get("ReactantSMILES", "")
            products_smi = row.get("ProductSMILES", "")
            smiles = {
                "reactants": reactants_smi,
                "products": products_smi,
            }
            
            # PRECOMPUTE: Normalize reaction, detect family, and featurize
            precomputed = {}
            if chemtools_available and reactants_smi and products_smi:
                try:
                    # Build reaction SMILES
                    reaction_smiles = f"{reactants_smi}>>{products_smi}"
                    
                    # 1. Normalize reaction
                    norm_result = normalize_reaction(reaction_smiles)
                    normalized_rxn = norm_result.get("normalized", "")
                    
                    # Extract normalized reactants
                    reactants_normalized = [
                        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
                        for r in (norm_result.get("reactants") or [])
                    ]
                    
                    if reactants_normalized:
                        # 2. Use existing reaction type from SciFinder (primary source)
                        # Map SciFinder naming to our family names
                        scifinder_type = (rtype or "").strip()
                        detected_family = "Unknown"
                        family_confidence = 0.0
                        family_source = "unknown"
                        
                        # Map common SciFinder reaction type names
                        scifinder_map = {
                            # C-N Coupling variations (unified dataset)
                            "buchwald": "C_N_Coupling",
                            "buchwald-hartwig": "C_N_Coupling",
                            "buchwald c-n": "C_N_Coupling",
                            "buchwald_cn": "C_N_Coupling",
                            "ullmann": "C_N_Coupling",
                            "ullmann c-n": "C_N_Coupling",
                            "ullmann_cn": "C_N_Coupling",
                            "goldberg": "C_N_Coupling",
                            "c-n coupling": "C_N_Coupling",
                            "c-n cross-coupling": "C_N_Coupling",
                            # C-O Coupling variations
                            "c-o coupling": "C_O_Coupling",
                            "c-o cross-coupling": "C_O_Coupling",
                            "ullmann c-o": "C_O_Coupling",
                            "ullmann_co": "C_O_Coupling",
                            "ullmann ether": "C_O_Coupling",
                            # C-S Coupling variations
                            "c-s coupling": "C_S_Coupling",
                            "c-s cross-coupling": "C_S_Coupling",
                            "ullmann c-s": "C_S_Coupling",
                            "ullmann_cs": "C_S_Coupling",
                            "ullmann thioether": "C_S_Coupling",
                            # Suzuki
                            "suzuki": "Suzuki_CC",
                            "suzuki-miyaura": "Suzuki_CC",
                            "suzuki coupling": "Suzuki_CC",
                            # Others
                            "sonogashira": "Sonogashira_CC",
                            "amide": "Amide_Coupling",
                            "amide coupling": "Amide_Coupling",
                            "amide formation": "Amide_Coupling",
                        }
                        
                        scifinder_lower = scifinder_type.lower()
                        if scifinder_lower in scifinder_map:
                            # Exact match - highest confidence
                            detected_family = scifinder_map[scifinder_lower]
                            family_source = "scifinder"
                            family_confidence = 1.0
                        elif scifinder_type:
                            # Try partial match
                            for key, family in scifinder_map.items():
                                if key in scifinder_lower or scifinder_lower in key:
                                    detected_family = family
                                    family_source = "scifinder_partial"
                                    family_confidence = 0.8
                                    break
                        
                        # Note: We do NOT run SMARTS detection here!
                        # SMARTS detection is expensive and only needed for user queries.
                        # For datasets with SciFinder metadata, we either:
                        #   1. Match SciFinder type to our families (above), or
                        #   2. Leave as "Unknown" - can be reviewed/corrected manually
                        # 
                        # If you need SMARTS detection, it should be done:
                        #   - On-demand for user input (in recommend.py, router.py)
                        #   - Not during bulk dataset preprocessing
                        
                        if detected_family == "Unknown" and scifinder_type:
                            # Keep SciFinder name even if we don't recognize it
                            # This helps manual review/mapping later
                            family_source = "scifinder_unmapped"
                            family_confidence = 0.0
                        
                        # 3. Featurize substrates (for categorical similarity fallback)
                        # Pick electrophile and nucleophile
                        def is_electrophile(s: str) -> bool:
                            t = (s or "").lower()
                            return ("br" in t) or ("cl" in t) or (" i" in t) or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
                        
                        elec_smi, nuc_smi = "", ""
                        if len(reactants_normalized) == 1:
                            elec_smi = reactants_normalized[0]
                        elif len(reactants_normalized) >= 2:
                            r0, r1 = reactants_normalized[0], reactants_normalized[1]
                            if is_electrophile(r0):
                                elec_smi, nuc_smi = r0, r1
                            elif is_electrophile(r1):
                                elec_smi, nuc_smi = r1, r0
                            else:
                                elec_smi, nuc_smi = r0, r1
                        
                        # Featurize (keep only hashable features)
                        features = feat_molecular.featurize(elec_smi, nuc_smi)
                        if isinstance(features, dict) and "role_aware" in features:
                            features = {k: v for k, v in features.items() if k != "role_aware"}
                        
                        # 4. Precompute DRFP fingerprint (optional but highly recommended)
                        # NOTE: DRFP is saved to a separate binary file, not embedded in JSONL
                        drfp_fp = None
                        if drfp_available and reaction_smiles:
                            try:
                                fp_array = rs.encode_drfp(reaction_smiles, n_bits=4096, radius=3)
                                if fp_array is not None:
                                    # Keep as numpy array for later binary storage
                                    drfp_fp = fp_array
                                    precompute_stats["drfp_computed"] += 1
                            except Exception:
                                pass  # DRFP computation failed, skip
                        
                        # Store precomputed data (without DRFP - saved separately)
                        precomputed = {
                            "reaction_smiles": reaction_smiles,
                            "normalized": normalized_rxn,
                            "reactants_normalized": reactants_normalized,
                            "detected_family": detected_family,
                            "family_confidence": round(family_confidence, 3),
                            "family_source": family_source,  # Track where family came from
                            "scifinder_type": scifinder_type,  # Keep original SciFinder name
                            "features": features,
                        }
                        
                        # DO NOT add DRFP to JSONL - it will be saved to binary .npz file
                        # This saves ~90% space (e.g., 670 MB → 70 MB + 12 MB .npz)
                        
                        # Collect DRFP for later binary storage
                        if drfp_fp is not None:
                            drfp_fingerprints.append(drfp_fp)
                            drfp_reaction_ids.append(rid)
                        
                        precompute_stats["success"] += 1
                        
                        # Track family source for statistics
                        family_sources[family_source] = family_sources.get(family_source, 0) + 1
                        
                        if idx % 100 == 0:
                            print(f"  Preprocessed {idx}/{len(rows)} reactions...")
                    else:
                        precompute_stats["skipped"] += 1
                        
                except Exception as e:
                    precompute_stats["failed"] += 1
                    if precompute_stats["failed"] <= 5:  # Only show first 5 errors
                        print(f"  Warning: Failed to preprocess reaction {rid}: {e}")
            elif not chemtools_available:
                precompute_stats["skipped"] += 1

            analysis_record = {
                "reaction_id": rid,
                "reaction_type": rtype,
                "condition_core": condition_core,
                "catalytic_system": catalytic_system,
                "reagents": reagents,
                "solvents": solvents,
                "conditions": conditions,
                "smiles": smiles,
                "reference": row.get("Reference") or {},
            }
            
            # Add precomputed data if available
            if precomputed:
                analysis_record["precomputed"] = precomputed
            
            out_lines.append(_json.dumps(analysis_record, ensure_ascii=False))

        # Write output file
        try:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write("\n".join(out_lines) + ("\n" if out_lines else ""))
            
            # Save DRFP fingerprints to binary NPZ file
            if drfp_fingerprints and drfp_reaction_ids:
                # Determine NPZ output path (same as JSONL but with .npz extension)
                npz_path = output_path.rsplit('.', 1)[0] + '_drfp.npz'
                
                try:
                    import numpy as np
                    
                    # Stack fingerprints into matrix
                    fps_matrix = np.vstack(drfp_fingerprints)
                    reaction_ids_array = np.array(drfp_reaction_ids, dtype=object)
                    
                    # Save as compressed NPZ
                    np.savez_compressed(
                        npz_path,
                        fps=fps_matrix,
                        reaction_ids=reaction_ids_array,
                        n_bits=np.array(4096, dtype='int32'),
                        radius=np.array(3, dtype='int32')
                    )
                    
                    npz_size_mb = os.path.getsize(npz_path) / (1024 * 1024)
                    print(f"\n✓ Saved {len(drfp_reaction_ids)} DRFP fingerprints to {npz_path}")
                    print(f"  Binary file size: {npz_size_mb:.2f} MB ({npz_size_mb/len(drfp_reaction_ids)*1000:.1f} KB per reaction)")
                except Exception as e:
                    print(f"\n⚠️  Warning: Failed to save DRFP binary file: {e}")
            
            # Print preprocessing statistics
            if chemtools_available:
                total = len(rows)
                print(f"\n{'='*60}")
                print(f"PREPROCESSING STATISTICS")
                print(f"{'='*60}")
                print(f"Total reactions:           {total}")
                print(f"Successfully preprocessed: {precompute_stats['success']} ({precompute_stats['success']/total*100:.1f}%)")
                print(f"Failed:                    {precompute_stats['failed']} ({precompute_stats['failed']/total*100:.1f}%)")
                print(f"Skipped:                   {precompute_stats['skipped']} ({precompute_stats['skipped']/total*100:.1f}%)")
                if drfp_available:
                    print(f"DRFP fingerprints:         {precompute_stats['drfp_computed']} ({precompute_stats['drfp_computed']/max(1,precompute_stats['success'])*100:.1f}%)")
                print(f"\nFamily Detection Sources:")
                print(f"  SciFinder exact match:   {family_sources['scifinder']} ({family_sources['scifinder']/max(1,precompute_stats['success'])*100:.1f}%)")
                print(f"  SciFinder partial match: {family_sources['scifinder_partial']} ({family_sources['scifinder_partial']/max(1,precompute_stats['success'])*100:.1f}%)")
                print(f"  SciFinder unmapped:      {family_sources['scifinder_unmapped']} ({family_sources['scifinder_unmapped']/max(1,precompute_stats['success'])*100:.1f}%)")
                print(f"  Unknown/Missing:         {family_sources['unknown']} ({family_sources['unknown']/max(1,precompute_stats['success'])*100:.1f}%)")
                print(f"{'='*60}")
                print(f"✓ Dataset saved with precomputed normalization and features!")
                print(f"✓ DRFP fingerprints saved to separate binary .npz file (saves ~90% space)")
                print(f"✓ Family names from SciFinder metadata (no SMARTS detection)")
                print(f"✓ Unmapped types can be added to scifinder_map for next import")
                print(f"✓ SMARTS detection will run on-demand for user queries only")
            
        except Exception as e:
            print(f"Error writing JSONL: {e}")

# Detect RDKit availability
try:
    from rdkit import Chem  # type: ignore
    RDKIT_AVAILABLE = True
except Exception:
    Chem = None  # type: ignore
    RDKIT_AVAILABLE = False


class RDFWorker(QtCore.QObject):
    finished = Signal(bool, str) if Signal else None  # type: ignore[misc]
    progress = Signal(str) if Signal else None  # type: ignore[misc]

    def __init__(self, folder_path: str, output_md_path: str, output_jsonl_path: str, process_unknowns: bool = False):
        super().__init__()
        self.folder_path = folder_path
        self.output_md_path = output_md_path
        self.output_jsonl_path = output_jsonl_path
        self.process_unknowns = process_unknowns
        self.rdf_files = []
        self._undetermined_loaded = False
        self._undetermined_existing_entries: List[Dict[str, Any]] = []
        self._undetermined_existing_map: Dict[str, Dict[str, Any]] = {}
        self._undetermined_new_map: Dict[str, Dict[str, Any]] = {}
        self._undetermined_file_path_cache: Optional[str] = None
        self._pubchem_cache_by_cas: Dict[str, Dict[str, Any]] = {}
        self._pubchem_cache_by_name: Dict[str, Dict[str, Any]] = {}
        self._pubchem_cache_by_cid: Dict[int, Dict[str, Any]] = {}
        self._pubchem_failed_keys: Set[str] = set()
        self._pubchem_error_count = 0
        self._pubchem_offline = False
        self._pubchem_last_error = ""

    def _emit(self, msg: str):
        """Emit progress message"""
        sig = getattr(self, 'progress', None)
        if sig:
            try:
                sig.emit(msg)
            except Exception:
                pass

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

    def _load_default_cas_maps(self) -> Dict[str, Dict[str, str]]:
        """Deprecated in this GUI workflow: use taxonomy-based mapping instead."""
        return {}

    def _load_taxonomy(self) -> _TaxonomyIndex:
        """Load the reagent registry and build a CAS/name map for roles and tokens."""
        # Default base dir relative to repo
        here = os.path.dirname(os.path.abspath(__file__))
        tax_dir = os.path.join(os.path.dirname(here), 'data', 'reagents')
        if not os.path.exists(tax_dir):
            # fallback to local data/reagents relative to this file
            tax_dir = os.path.join(here, 'data', 'reagents')
        try:
            idx = _TaxonomyIndex(tax_dir)
            return idx
        except Exception as e:
            raise RuntimeError(f"Failed to load taxonomy from {tax_dir}: {e}")

    @staticmethod
    def _normalize_name(value: str) -> str:
        """Normalize a compound name for lookups."""
        import re as _re
        if not value:
            return ""
        normalized = value.strip().lower()
        normalized = _re.sub(r"\s+", " ", normalized)
        normalized = _re.sub(r"[^a-z0-9\+\-\.\(\)\[\]/']", "", normalized)
        return normalized

    def _undetermined_file_path(self) -> str:
        if self._undetermined_file_path_cache:
            return self._undetermined_file_path_cache
        here = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(here)
        target_dir = os.path.join(repo_root, "data", "reagents")
        target_path = os.path.join(target_dir, "not_determined_reagents.json")
        legacy_path = os.path.join(here, "not_determined_reagents.json")
        if not os.path.exists(target_path) and os.path.exists(legacy_path):
            self._emit("Note: using legacy not_determined_reagents.json from data-processor directory; it will be rewritten to data/reagents on save.")
        self._undetermined_file_path_cache = target_path
        return target_path

    def _ensure_undetermined_cache_loaded(self) -> None:
        if self._undetermined_loaded:
            return
        self._undetermined_loaded = True
        self._undetermined_existing_entries = []
        self._undetermined_existing_map = {}
        primary_path = self._undetermined_file_path()
        load_paths = [primary_path]
        here = os.path.dirname(os.path.abspath(__file__))
        legacy_path = os.path.join(here, "not_determined_reagents.json")
        if primary_path != legacy_path and os.path.exists(legacy_path):
            load_paths.append(legacy_path)
        data = None
        for candidate in load_paths:
            if not os.path.exists(candidate):
                continue
            try:
                with open(candidate, "r", encoding="utf-8") as f:
                    data = _json.load(f)
                break
            except Exception as exc:
                self._emit(f"Note: could not load existing not_determined_reagents.json from {candidate}: {exc}")
                data = None
        if data is None:
            return

        taxonomy = getattr(self, "_taxonomy_index", None)
        if isinstance(data, list):
            for entry in data:
                if not isinstance(entry, dict):
                    continue
                cas = str(entry.get("cas") or "").strip()
                name = str(entry.get("name") or "").strip()
                base_key = self._make_undetermined_key(cas, name, taxonomy)
                upgraded = self._upgrade_undetermined_entry(entry, taxonomy)
                self._undetermined_existing_entries.append(upgraded)
                new_cas = str(upgraded.get("cas") or "").strip()
                new_name = str(upgraded.get("name") or "").strip()
                new_key = self._make_undetermined_key(new_cas, new_name, taxonomy)
                keys = {k for k in (base_key, new_key) if k}
                if not keys:
                    continue
                for k in keys:
                    self._undetermined_existing_map[k] = upgraded

    def _make_undetermined_key(self, cas: str, name: str, taxonomy: Optional[_TaxonomyIndex] = None) -> str:
        cas_norm = (cas or "").strip()
        if cas_norm:
            return f"CAS::{cas_norm}"
        basis = ""
        if taxonomy is not None:
            basis = taxonomy._norm_name(name or "")
        else:
            basis = self._normalize_name(name)
        if not basis:
            return ""
        return f"NAME::{basis}"

    def _record_pubchem_success(self) -> None:
        self._pubchem_last_error = ""

    def _record_pubchem_error(self, reason: str) -> None:
        self._pubchem_error_count += 1
        self._pubchem_last_error = reason
        if self._pubchem_error_count >= 3 and not self._pubchem_offline:
            self._pubchem_offline = True
            self._emit(
                f"Notice: PubChem service appears unavailable (last error: {reason}). Skipping further lookups for this run."
            )

    def _pubchem_json(self, url: str) -> Optional[Dict[str, Any]]:
        if self._pubchem_offline:
            return None
        try:
            req = urllib.request.Request(url, headers={"User-Agent": PUBCHEM_USER_AGENT})
            with urllib.request.urlopen(req, timeout=PUBCHEM_TIMEOUT) as resp:
                payload = resp.read()
            self._record_pubchem_success()
            return _json.loads(payload.decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                self._record_pubchem_success()
                return None
            self._emit(f"Note: PubChem request failed ({exc.code}): {url}")
            self._record_pubchem_error(f"HTTP {exc.code}")
        except Exception as exc:
            self._emit(f"Note: PubChem request error: {exc}: {url}")
            self._record_pubchem_error(str(exc))
        return None

    def _extract_cas_from_synonyms(self, synonyms: List[str]) -> str:
        for syn in synonyms or []:
            syn_stripped = (syn or "").strip()
            if CAS_NUMBER_RE.fullmatch(syn_stripped):
                return syn_stripped
        return ""

    def _choose_abbreviation(self, synonyms: List[str], canonical_name: str) -> str:
        candidates = []
        for syn in synonyms or []:
            s = (syn or "").strip()
            if not s:
                continue
            if CAS_NUMBER_RE.fullmatch(s):
                continue
            if len(s) <= 20:
                candidates.append(s)
        if candidates:
            return candidates[0]
        canonical = (canonical_name or "").strip()
        if canonical and len(canonical) <= 24:
            return canonical
        return ""

    def _extract_use_notes(self, record: Dict[str, Any]) -> str:
        def walk_sections(sections: List[Dict[str, Any]]) -> List[str]:
            snippets: List[str] = []
            for section in sections or []:
                heading = (section.get("TOCHeading") or "").lower()
                if "use" in heading or "application" in heading or "function" in heading:
                    for info in section.get("Information") or []:
                        val = info.get("Value") or {}
                        for item in val.get("StringWithMarkup") or []:
                            txt = (item.get("String") or "").strip()
                            if txt:
                                snippets.append(txt)
                child = walk_sections(section.get("Section") or [])
                if child:
                    snippets.extend(child)
            return snippets
        record_sections = record.get("Record", {}).get("Section") or []
        pieces = walk_sections(record_sections)
        unique: List[str] = []
        seen: Set[str] = set()
        for piece in pieces:
            if piece in seen:
                continue
            seen.add(piece)
            unique.append(piece)
            if len(unique) >= 3:
                break
        return " ".join(unique)

    def _fetch_pubchem_uses(self, cid: int) -> str:
        url = f"{PUBCHEM_VIEW_BASE}/{cid}/JSON"
        data = self._pubchem_json(url)
        if not data:
            return ""
        return self._extract_use_notes(data)

    def _fetch_pubchem_compound(self, *, cas: str = "", name: str = "") -> Dict[str, Any]:
        cas = (cas or "").strip()
        name = (name or "").strip()
        if not cas and not name:
            return {}
        cid = None
        if cas:
            data = self._pubchem_json(f"{PUBCHEM_PUG_BASE}/compound/xref/RN/{urllib.parse.quote(cas)}/JSON")
            if data and data.get("PC_Compounds"):
                cid = data["PC_Compounds"][0].get("id", {}).get("id", {}).get("cid")
        if cid is None and name:
            data = self._pubchem_json(f"{PUBCHEM_PUG_BASE}/compound/name/{urllib.parse.quote(name)}/cids/JSON")
            if data:
                cid_list = data.get("IdentifierList", {}).get("CID") or []
                if cid_list:
                    cid = cid_list[0]
        if cid is None:
            return {}
        if cid in self._pubchem_cache_by_cid:
            return self._pubchem_cache_by_cid[cid]
        prop_data = self._pubchem_json(f"{PUBCHEM_PUG_BASE}/compound/cid/{cid}/property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES/JSON")
        synonyms_data = self._pubchem_json(f"{PUBCHEM_PUG_BASE}/compound/cid/{cid}/synonyms/JSON")
        synonyms = []
        if synonyms_data:
            info_list = synonyms_data.get("InformationList", {}).get("Information") or []
            if info_list:
                synonyms = list(info_list[0].get("Synonym") or [])

        props = {}
        if prop_data:
            prop_list = prop_data.get("PropertyTable", {}).get("Properties") or []
            if prop_list:
                props = prop_list[0]
        canonical_name = str(props.get("IUPACName") or "").strip()
        synonyms = [s for s in synonyms if s]
        cas_candidate = cas or self._extract_cas_from_synonyms(synonyms)
        abbreviation = self._choose_abbreviation(synonyms, canonical_name)
        uses_text = self._fetch_pubchem_uses(cid)
        info = {
            "pubchem_cid": cid,
            "cas": cas_candidate,
            "name": canonical_name or (synonyms[1] if len(synonyms) > 1 else (synonyms[0] if synonyms else name or cas_candidate)),
            "canonical_smiles": props.get("CanonicalSMILES"),
            "isomeric_smiles": props.get("IsomericSMILES"),
            "molecular_formula": props.get("MolecularFormula"),
            "molecular_weight": props.get("MolecularWeight"),
            "iupac_name": props.get("IUPACName"),
            "synonyms": synonyms,
            "abbreviation": abbreviation,
            "typical_usage": uses_text,
            "data_source": "PubChem",
            "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
        }
        self._pubchem_cache_by_cid[cid] = info
        if cas_candidate:
            self._pubchem_cache_by_cas[cas_candidate] = info
        if name:
            norm_name = self._normalize_name(name)
            if norm_name:
                self._pubchem_cache_by_name[norm_name] = info
        for syn in synonyms:
            norm_syn = self._normalize_name(syn)
            if norm_syn:
                self._pubchem_cache_by_name.setdefault(norm_syn, info)
        return info

    def _lookup_compound_info(self, cas: str, name: str) -> Dict[str, Any]:
        cas_norm = (cas or "").strip()
        name_norm = (name or "").strip()
        if cas_norm and cas_norm in self._pubchem_cache_by_cas:
            return self._pubchem_cache_by_cas[cas_norm]
        name_key = self._normalize_name(name_norm)
        if name_key and name_key in self._pubchem_cache_by_name:
            return self._pubchem_cache_by_name[name_key]
        fail_keys = []
        if cas_norm:
            fail_key = f"cas:{cas_norm}"
            if fail_key in self._pubchem_failed_keys:
                cas_norm = ""
            else:
                fail_keys.append(fail_key)
        if name_key:
            fail_name = f"name:{name_key}"
            if fail_name in self._pubchem_failed_keys:
                name_key = ""
            else:
                fail_keys.append(fail_name)
        info: Dict[str, Any] = {}
        if cas_norm:
            info = self._fetch_pubchem_compound(cas=cas_norm)
        if not info and name_norm:
            info = self._fetch_pubchem_compound(name=name_norm)
        if info:
            return info
        for key in fail_keys:
            self._pubchem_failed_keys.add(key)
        return {}

    def _sanitize_compound_value(self, value: Any):
        if isinstance(value, str):
            val = value.strip()
            return val or None
        if isinstance(value, (int, float, bool)):
            return value
        if isinstance(value, list):
            cleaned = []
            for item in value:
                cleaned_item = self._sanitize_compound_value(item)
                if cleaned_item not in (None, "", [], {}):
                    cleaned.append(cleaned_item)
            return cleaned or None
        if isinstance(value, tuple):
            cleaned = []
            for item in value:
                cleaned_item = self._sanitize_compound_value(item)
                if cleaned_item not in (None, "", [], {}):
                    cleaned.append(cleaned_item)
            return cleaned or None
        if isinstance(value, set):
            cleaned = []
            for item in value:
                cleaned_item = self._sanitize_compound_value(item)
                if cleaned_item not in (None, "", [], {}):
                    cleaned.append(cleaned_item)
            return cleaned or None
        if isinstance(value, dict):
            cleaned_dict = {}
            for k, v in value.items():
                cleaned_val = self._sanitize_compound_value(v)
                if cleaned_val not in (None, "", [], {}):
                    cleaned_dict[str(k)] = cleaned_val
            return cleaned_dict or None
        return value

    def _build_undetermined_entry(self, *, cas: str, reported_name: str, role: str, info: Dict[str, Any], taxonomy: Optional[_TaxonomyIndex]) -> Dict[str, Any]:
        entry: Dict[str, Any] = {}
        cas_clean = (cas or "").strip()
        reported = (reported_name or "").strip()
        role_clean = (role or "UNK").strip().upper()
        canonical_name = str(info.get("name") or "").strip()
        display_name = canonical_name or reported or cas_clean
        if cas_clean:
            entry["cas"] = cas_clean
        if display_name:
            entry["name"] = display_name
        if reported and reported != display_name:
            entry["reported_name"] = reported
        if role_clean:
            entry["role"] = role_clean
        abbr = info.get("abbreviation") or info.get("abbr")
        if isinstance(abbr, str):
            abbr = abbr.strip()
        if abbr:
            entry["abbreviation"] = abbr
        synonyms = info.get("synonyms") or []
        if isinstance(synonyms, (list, tuple, set)):
            cleaned_syns = sorted({str(s).strip() for s in synonyms if str(s).strip()})[:5]
            if cleaned_syns:
                entry["synonyms"] = cleaned_syns
        copy_keys = [
            "compound_type",
            "category_hint",
            "generic_core",
            "token",
            "notes",
            "usage",
            "typical_usage",
            "data_sources",
            "sources",
            "canonical_smiles",
            "isomeric_smiles",
            "molecular_formula",
            "molecular_weight",
            "iupac_name",
            "pubchem_cid",
            "pubchem_url",
            "data_source",
        ]
        for key in copy_keys:
            if key not in info:
                continue
            sanitized = self._sanitize_compound_value(info.get(key))
            if sanitized not in (None, "", [], {}):
                entry[key] = sanitized
        smiles = info.get("smiles") or info.get("smile") or info.get("canonical_smiles") or info.get("isomeric_smiles")
        if isinstance(smiles, str):
            smiles = smiles.strip()
        if smiles:
            entry["smiles"] = smiles
        typical_usage = entry.get("typical_usage") or entry.get("usage")
        if isinstance(typical_usage, list):
            typical_usage = ", ".join(str(x).strip() for x in typical_usage if str(x).strip())
        elif isinstance(typical_usage, dict):
            typical_usage = ", ".join(f"{k}: {v}" for k, v in typical_usage.items() if v not in (None, "", [], {}))
        if not typical_usage:
            fallback = [entry.get("category_hint"), entry.get("compound_type"), role_clean if role_clean and role_clean != "UNK" else None]
            typical_usage = next((str(val).strip() for val in fallback if isinstance(val, str) and val.strip()), "")
        if typical_usage:
            entry["typical_usage"] = typical_usage
        if entry.get("usage") == entry.get("typical_usage"):
            entry.pop("usage", None)
        if info:
            entry["data_source"] = info.get("data_source") or "PubChem"
        cleaned_entry: Dict[str, Any] = {}
        for key, value in entry.items():
            if key in {"cas", "name", "role"}:
                cleaned_entry[key] = value
                continue
            if value in (None, "", [], {}):
                continue
            cleaned_entry[key] = value
        if "reported_name" in entry and entry["reported_name"]:
            cleaned_entry["reported_name"] = entry["reported_name"]
        return cleaned_entry

    def _upgrade_undetermined_entry(self, entry: Dict[str, Any], taxonomy: Optional[_TaxonomyIndex]) -> Dict[str, Any]:
        cas = str(entry.get("cas") or "").strip()
        reported = str(entry.get("reported_name") or entry.get("raw_name") or entry.get("name") or "").strip()
        lookup_name = str(entry.get("name") or reported).strip()
        role = (entry.get("role") or "UNK").strip().upper()
        if entry.get("data_source"):
            info = dict(entry)
        else:
            info = self._lookup_compound_info(cas, lookup_name)
        if not cas:
            cas_candidate = str((info or {}).get("cas") or "").strip()
            if cas_candidate:
                cas = cas_candidate
        built = self._build_undetermined_entry(
            cas=cas,
            reported_name=reported or lookup_name,
            role=role or "UNK",
            info=info or {},
            taxonomy=taxonomy
        )
        for key, value in entry.items():
            if key in {"cas", "name", "role"}:
                continue
            if value in (None, "", [], {}):
                continue
            built.setdefault(key, value)
        return built

    def _is_missing_from_taxonomy(self, name: str, cas: str, taxonomy: _TaxonomyIndex) -> bool:
        cas_norm = (cas or "").strip()
        if cas_norm and cas_norm in taxonomy.cas_map:
            return False
        norm_name = taxonomy._norm_name(name or "")
        if norm_name and norm_name in taxonomy.name_to_cas:
            return False
        return True

    def _record_undetermined_reagent(self, *, name: str, cas: str, taxonomy: _TaxonomyIndex, role: str) -> None:
        if not self.process_unknowns:
            return
        cas = (cas or "").strip()
        reported = (name or "").strip()
        if not cas and not reported:
            return
        info = self._lookup_compound_info(cas, reported)
        if not cas:
            cas_candidate = str((info or {}).get("cas") or "").strip()
            if cas_candidate:
                cas = cas_candidate
        canonical_name = str((info or {}).get("name") or reported).strip()
        candidate_keys = {
            self._make_undetermined_key(cas, reported, taxonomy),
            self._make_undetermined_key(cas, canonical_name, taxonomy),
        }
        candidate_keys = {k for k in candidate_keys if k}
        if not candidate_keys:
            return
        self._ensure_undetermined_cache_loaded()
        for key in candidate_keys:
            if key in self._undetermined_existing_map or key in self._undetermined_new_map:
                return
        chosen_key = sorted(candidate_keys)[0]
        entry = self._build_undetermined_entry(
            cas=cas,
            reported_name=reported,
            role=role or "UNK",
            info=info or {},
            taxonomy=taxonomy
        )
        self._undetermined_new_map[chosen_key] = entry
        for key in candidate_keys:
            if key != chosen_key:
                self._undetermined_new_map.setdefault(key, entry)

    def _persist_undetermined_reagents(self) -> None:
        if not self.process_unknowns or not self._undetermined_new_map:
            return
        path = self._undetermined_file_path()
        entries = list(self._undetermined_existing_entries)
        unique_new: Dict[Tuple[str, str], Dict[str, Any]] = {}
        for entry in self._undetermined_new_map.values():
            cas_key = str(entry.get("cas") or "").strip()
            name_key = str(entry.get("name") or "").strip()
            unique_new.setdefault((cas_key, name_key), entry)
        new_entries = sorted(unique_new.values(), key=lambda e: ((e.get("cas") or "").strip(), (e.get("name") or "").strip()))
        entries.extend(new_entries)
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w", encoding="utf-8") as f:
                _json.dump(entries, f, ensure_ascii=False, indent=2)
        except Exception as exc:
            self._emit(f"Note: failed to write undetermined reagents: {exc}")
            return
        self._undetermined_existing_entries = entries
        self._undetermined_existing_map = {}
        taxonomy = getattr(self, "_taxonomy_index", None)
        for entry in entries:
            cas_val = str(entry.get("cas") or "").strip()
            name_val = str(entry.get("name") or "").strip()
            reported_val = str(entry.get("reported_name") or "").strip()
            keys = {
                self._make_undetermined_key(cas_val, name_val, taxonomy),
                self._make_undetermined_key(cas_val, reported_val, taxonomy),
            }
            for key in keys:
                if key:
                    self._undetermined_existing_map[key] = entry
        added = len(new_entries)
        self._undetermined_new_map = {}
        self._emit(f"Captured {added} undetermined reagents in {path}")

    def _reassign_reagent_roles_via_taxonomy(self, rows: List[Dict[str, Any]], taxonomy: _TaxonomyIndex) -> None:
        """Ensure ReagentRole aligns with Reagent list using taxonomy role lookup when possible."""
        track_unknowns = self.process_unknowns
        if track_unknowns:
            self._ensure_undetermined_cache_loaded()
        for row in rows or []:
            try:
                reag_json = row.get('Reagent') or '[]'
                reag_list = []
                try:
                    reag_list = _json.loads(reag_json)
                except Exception:
                    reag_list = []
                roles: List[str] = []
                for item in reag_list:
                    name, cas = '', ''
                    if isinstance(item, str):
                        if '|' in item:
                            name, cas = (item.split('|', 1) + [''])[:2]
                            name = (name or '').strip()
                            cas = (cas or '').strip()
                        else:
                            name = item.strip()
                    elif isinstance(item, dict):
                        name = (item.get('name') or '').strip()
                        cas = (item.get('cas') or '').strip()
                    role = taxonomy.role_for(name, cas)
                    if not role or role == 'UNK':
                        if track_unknowns and self._is_missing_from_taxonomy(name, cas, taxonomy):
                            self._record_undetermined_reagent(name=name, cas=cas, taxonomy=taxonomy, role=role or 'UNK')
                        role = 'UNK'
                    roles.append(role)
                # Only write back if we can align 1:1
                if roles and (len(roles) == len(reag_list)):
                    row['ReagentRole'] = _json.dumps(roles, ensure_ascii=False)
            except Exception:
                continue

    def _inject_suggested_ligands_via_taxonomy(self, rows: List[Dict[str, Any]], taxonomy: _TaxonomyIndex) -> int:
        """Populate suggested ligand into row['Ligand'] (and FullCatalyticSystem) when missing.

        Returns the count of rows updated.
        """
        updated = 0
        # Build PairingHelper once
        try:
            from conditioncore_pairing_helper_for_ref_only import PairingHelper  # type: ignore
            cat_path = os.path.join(taxonomy.base_dir, 'taxonomy_catalysts_precursor.json')
            ph = PairingHelper(cat_path)
        except Exception:
            ph = None  # type: ignore

        for row in rows or []:
            try:
                # Skip if ligand already present
                lig_list = []
                try:
                    lig_list = _json.loads(row.get('Ligand') or '[]')
                except Exception:
                    lig_list = []
                if lig_list:
                    continue

                # Determine catalyst family from core CAS
                core_pairs = []
                try:
                    core_pairs = _json.loads(row.get('CatalystCoreDetail') or '[]')
                except Exception:
                    core_pairs = []
                fam_id = None
                for p in core_pairs or []:
                    if not isinstance(p, str):
                        continue
                    _, _, cs = p.partition('|')
                    cs = cs.strip()
                    if cs and cs in taxonomy.cas_to_family:
                        fam_id = taxonomy.cas_to_family.get(cs)
                        if fam_id:
                            break

                if not fam_id or not ph:
                    continue

                # Use reaction type as hint
                hint = (row.get('ReactionType') or '').strip() or None
                pick = ph.suggest_for(fam_id, None, hint)
                if not pick:
                    continue
                abbr = ((pick.get('ligand') or {}).get('abbr') or '').strip()
                if not abbr:
                    continue

                # Map abbr to CAS via taxonomy; fall back to name-only pair
                cas = taxonomy.name_to_cas.get(taxonomy._norm_name(abbr), '')
                pair = f"{abbr}|{cas}"

                # Update Ligand list
                new_lig = [pair]
                row['Ligand'] = _json.dumps(new_lig, ensure_ascii=False)

                # Update FullCatalyticSystem
                full = []
                try:
                    full = _json.loads(row.get('FullCatalyticSystem') or '[]')
                except Exception:
                    full = []
                if pair not in full:
                    full.append(pair)
                row['FullCatalyticSystem'] = _json.dumps(full, ensure_ascii=False)

                # Recompute keys
                try:
                    row['CondKey'] = build_condkey(row)
                    row['CondSig'] = build_condsig(row)
                    row['FamSig'] = build_famsig(row)
                except Exception:
                    pass

                updated += 1
            except Exception:
                continue
        return updated

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
        """Process all RDF files and combine them into a single RDF map"""
        combined_rdf_map: Dict[str, Dict[str, Any]] = {}
        seen_ids: set[str] = set()
        for i, rdf_file in enumerate(self.rdf_files, 1):
            filename = os.path.basename(rdf_file)
            self._emit(f"[{i}/{len(self.rdf_files)}] Processing {filename}...")
            try:
                # Parse individual RDF file
                rdf_map = parse_rdf(rdf_file)
                # Merge reactions without prefixing filename to the ID; keep first occurrence only
                added = 0
                skipped = 0
                for rid, data in rdf_map.items():
                    data['source_file'] = filename
                    if rid in seen_ids or rid in combined_rdf_map:
                        skipped += 1
                        continue
                    seen_ids.add(rid)
                    combined_rdf_map[rid] = data
                    added += 1
                msg_tail = f" (added {added}"
                if skipped:
                    msg_tail += f", skipped dups {skipped}"
                msg_tail += ")"
                self._emit(f"  Found {len(rdf_map)} reactions in {filename}{msg_tail}")
            except Exception as e:
                self._emit(f"  Error processing {filename}: {e}")
                continue
        return combined_rdf_map

    def _generate_outputs(self, rows: List[Dict[str, Any]], cas_map: Dict[str, Dict[str, str]]) -> None:
        """Generate Markdown and JSONL outputs using ReactionMarkdownGenerator"""
        self._emit("Generating Markdown report...")
        
        # Create taxonomy-aware generator instance
        taxonomy = getattr(self, '_taxonomy_index', None)
        generator = ReactionMarkdownGenerator(taxonomy=taxonomy)
        generator.cas_map = cas_map or (taxonomy.cas_map if taxonomy else {})
        
        # Generate markdown report
        source_name = f"RDF_Folder_{os.path.basename(self.folder_path)}"
        generator.generate_markdown_report(rows, self.output_md_path, source_name)
        
        # Generate JSONL export
        self._emit("Generating JSONL export...")
        generator.generate_jsonl_export(rows, self.output_jsonl_path, source_name)
        
        # NOTE: No need to call export_jsonl_for_chemtools anymore
        # We now save directly to data/reaction_dataset/ instead of copying
        # The DRFP binary file (.npz) is also saved alongside the JSONL automatically

    @Slot() if Slot else (lambda f: f)
    def run(self):
        """Main processing function"""
        try:
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
            
            # Load CAS mappings
            self._emit("Loading compound taxonomy (roles, ligands, catalysts)...")
            taxonomy = self._load_taxonomy()
            self._taxonomy_index = taxonomy
            cas_map = taxonomy.cas_map
            
            # Create minimal TXT map (since we only have RDF)
            self._emit("Creating minimal TXT mapping...")
            txt_map = self._create_minimal_txt_map(combined_rdf_map)
            
            # Assemble rows using the same pipeline as original
            self._emit("Assembling reaction rows...")
            rows = assemble_rows(txt_map, combined_rdf_map, cas_map, txt_preferred=False)
            self._emit(f"Assembled {len(rows)} rows")

            # Post-process reagent roles using taxonomy to improve role assignment
            self._emit("Assigning reagent roles via taxonomy...")
            self._reassign_reagent_roles_via_taxonomy(rows, taxonomy)
            if self.process_unknowns:
                self._persist_undetermined_reagents()

            # Inject suggested ligands into rows where Ligand list is empty
            self._emit("Suggesting and adding ligands for catalyst families lacking ligands...")
            try:
                n = self._inject_suggested_ligands_via_taxonomy(rows, taxonomy)
                self._emit(f"Added suggested ligands to {n} reactions.")
            except Exception as e:
                self._emit(f"Ligand suggestion phase skipped due to error: {e}")

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

            # Override ReactionType using the parent folder name (e.g., ...\Suzuki\2023-2025 -> 'Suzuki')
            try:
                norm_folder = os.path.normpath(self.folder_path)
                parent_dir = os.path.basename(os.path.dirname(norm_folder))
                if parent_dir:
                    for r in rows:
                        r['ReactionType'] = parent_dir
                self._emit(f"Reaction type set to folder category: {parent_dir}")
            except Exception:
                # Non-fatal; keep existing reaction types if path parsing fails
                pass
            
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
            self._generate_outputs(rows, cas_map)
            
            if self.finished:
                self.finished.emit(
                    True, 
                    f"Successfully processed {len(self.rdf_files)} RDF files with {len(rows)} reactions.\n\n"
                    f"📝 Markdown (records): {self.output_md_path}\n"
                    f"📊 JSONL (chemtools): {self.output_jsonl_path}\n"
                    f"🔬 DRFP binary: {self.output_jsonl_path.rsplit('.', 1)[0] + '_drfp.npz'}"
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
        self.output_md_edit = QtWidgets.QLineEdit()
        self.btn_output_md = QtWidgets.QPushButton("Save As...")
        self.unknowns_checkbox = QtWidgets.QCheckBox("Save unknown compounds to not_determined_reagents.json")
        
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
        self.btn_output_md.clicked.connect(self.pick_output)
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
        
        # Output file selection
        output_box = QtWidgets.QHBoxLayout()
        output_box.addWidget(self.output_md_edit)
        output_box.addWidget(self.btn_output_md)
        form.addRow("Output Markdown:", output_box)
        
        # Add note about file locations
        note_label = QtWidgets.QLabel(
            "Note: All RDF files in folder and subfolders will be combined\n"
            "      JSONL → data/reaction_dataset/{category}.jsonl\n"
            "      Markdown → selected folder/{category}.md"
        )
        note_label.setStyleSheet("font-style: italic; color: #666; font-size: 10px;")
        form.addRow("", note_label)
        
        form.addRow("Save unknowns:", self.unknowns_checkbox)
        
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
            
            # Suggest default output file
            if not self.output_md_edit.text().strip():
                default_output = os.path.join(path, "rdf_reactions_rich.md")
                self.output_md_edit.setText(default_output)

    def pick_output(self):
        """Select output markdown file location"""
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save Markdown Report As",
            os.getcwd(),
            "Markdown files (*.md);;All files (*.*)"
        )
        if path:
            if not path.lower().endswith('.md'):
                path += '.md'
            self.output_md_edit.setText(path)

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
        output_md = self.output_md_edit.text().strip()
        
        if not folder or not os.path.isdir(folder):
            return "Please select a valid folder containing RDF files."
        
        if not output_md:
            return "Please specify an output Markdown file location."
        
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
        self.setEnabled(False)
        self.log.clear()
        self.log_msg("Starting RDF processing...")
        
        # Calculate output paths
        # Name JSONL using folder hierarchy, e.g., C_N_Coupling/2020-2022 -> C_N_Coupling.jsonl
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
            
            # JSONL: Save directly to data/reaction_dataset/ for chemtools consumption
            # Use category name (all subfolders are automatically included)
            repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            dataset_dir = os.path.join(repo_root, "data", "reaction_dataset")
            os.makedirs(dataset_dir, exist_ok=True)
            
            final_name = _safe(category_name) or "dataset"
            # Use simple category name for compatibility with chemtools
            output_jsonl = os.path.join(dataset_dir, final_name + ".jsonl")
            
            # Markdown: Save to the selected folder (preserves subfolder structure)
            output_md = os.path.join(norm_folder, final_name + ".md")
            
            try:
                # Reflect the computed Markdown path back into the UI
                self.output_md_edit.setText(output_md)
            except Exception:
                pass
        except Exception:
            # Fallback: use default paths
            repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            dataset_dir = os.path.join(repo_root, "data", "reaction_dataset")
            os.makedirs(dataset_dir, exist_ok=True)
            output_jsonl = os.path.join(dataset_dir, "dataset.jsonl")
            
            # Markdown to original_dataset
            original_dataset_dir = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "original_dataset"
            )
            os.makedirs(original_dataset_dir, exist_ok=True)
            output_md = os.path.join(original_dataset_dir, "dataset.md")
        
        # Create worker and thread
        self.worker = RDFWorker(
            folder_path=self.folder_edit.text().strip(),
            output_md_path=output_md,
            output_jsonl_path=output_jsonl,
            process_unknowns=self.unknowns_checkbox.isChecked()
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
