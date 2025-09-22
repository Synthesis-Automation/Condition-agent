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
from typing import List, Optional, Dict, Any

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


class _TaxonomyIndex:
    """Load compound taxonomies and expose simple lookup utilities.

    Builds a lightweight CAS map compatible with assemble_rows(), plus name/synonym
    indices for role assignment and tokenization.
    """

    def __init__(self, base_dir: str) -> None:
        self.base_dir = base_dir
        self.cas_map: Dict[str, Dict[str, str]] = {}
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
        t = re.sub(r"[^a-z0-9\+\-\.\(\)\[\]/']", "", t)  # keep common chem punctuation
        return t

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
    ) -> None:
        cas = (cas or "").strip()
        name = (name or "").strip()
        abbr = (abbr or "").strip()
        synonyms = list(synonyms or [])

        token = (abbr or name).strip().replace(" ", "")
        rec = {
            "Name": name or abbr or cas,
            "Abbreviation": abbr,
            "Token": token,
            "Role": role,
            "CategoryHint": category_hint,
            "GenericCore": generic_core,
        }
        if family_id:
            rec["FamilyID"] = family_id
        if cas:
            self.cas_map[cas] = rec
        # Names -> CAS and role indices
        keys = [name, abbr] + synonyms
        for k in keys:
            k_norm = self._norm_name(k)
            if not k_norm:
                continue
            if cas and k_norm not in self.name_to_cas:
                self.name_to_cas[k_norm] = cas
            # prefer stable role labels; do not override existing
            self.name_to_role.setdefault(k_norm, role)
        if family_id and cas:
            self.cas_to_family[cas] = family_id
        # family metal map set elsewhere for catalysts

    def _load_json(self, rel: str) -> Dict[str, Any] | None:
        p = os.path.join(self.base_dir, rel)
        if not os.path.exists(p):
            return None
        try:
            with open(p, "r", encoding="utf-8") as f:
                return _json.load(f)
        except Exception:
            return None

    def _load_all(self) -> None:
        # Ligands
        lig = self._load_json("taxonomy_ligand.json")
        if lig:
            for fam in lig.get("families", []) or []:
                fam_id = fam.get("family_id") or None
                for em in fam.get("example_members", []) or []:
                    self._index_member(
                        cas=em.get("cas"),
                        name=em.get("name"),
                        abbr=em.get("abbr"),
                        synonyms=em.get("synonyms") or [],
                        role="LIG",
                        category_hint="ligand",
                    )

        # Bases
        bas = self._load_json("taxonomy_base.json")
        if bas:
            for fam in bas.get("families", []) or []:
                for em in fam.get("example_members", []) or []:
                    self._index_member(
                        cas=em.get("cas"),
                        name=em.get("name"),
                        abbr=em.get("abbr"),
                        synonyms=em.get("synonyms") or [],
                        role="BASE",
                        category_hint="base",
                    )

        # Coupling reagents / activators (amidation)
        cou = self._load_json("taxonomy_coupling_reagent.json")
        if cou:
            for fam in cou.get("families", []) or []:
                for em in fam.get("example_members", []) or []:
                    self._index_member(
                        cas=em.get("cas"),
                        name=em.get("name"),
                        abbr=em.get("abbr"),
                        synonyms=em.get("synonyms") or [],
                        role="COUPLING_REAGENT",
                        category_hint="coupling_reagent",
                    )

        # Solvents
        solv = self._load_json("taxonomy_solvent.json")
        if solv:
            for fam in solv.get("families", []) or []:
                for em in fam.get("example_members", []) or []:
                    self._index_member(
                        cas=em.get("cas"),
                        name=em.get("name"),
                        abbr=em.get("abbr"),
                        synonyms=em.get("synonyms") or [],
                        role="SOLVENT",
                        category_hint="solvent",
                    )

        # Catalyst precursors (metals)
        cat = self._load_json("taxonomy_catalysts_precursor.json")
        if cat:
            for fam in cat.get("families", []) or []:
                fam_id = fam.get("family_id")
                metal = (fam.get("metal") or "").strip()
                if fam_id and metal:
                    self.family_metal[fam_id] = metal
                for em in fam.get("example_members", []) or []:
                    self._index_member(
                        cas=em.get("cas"),
                        name=em.get("name"),
                        abbr=em.get("abbr"),
                        synonyms=em.get("synonyms") or [],
                        role="CAT",
                        category_hint="catalyst_precursor",
                        generic_core=metal,
                        family_id=fam_id,
                    )

    # ---- Public helpers ----

    def cas_lookup(self, cas: str) -> Dict[str, str]:
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
        return "UNK"

    def ligand_token_for(self, name: str | None, cas: str | None) -> str:
        cas = (cas or "").strip()
        if cas and cas in self.cas_map:
            tok = (self.cas_map[cas].get("Token") or self.cas_map[cas].get("Abbreviation") or self.cas_map[cas].get("Name") or "").strip()
            return tok.replace(" ", "")
        n = (name or "").strip()
        if n:
            nn = self._norm_name(n)
            cas2 = self.name_to_cas.get(nn)
            if cas2 and cas2 in self.cas_map:
                tok = (self.cas_map[cas2].get("Token") or self.cas_map[cas2].get("Abbreviation") or self.cas_map[cas2].get("Name") or "").strip()
                return tok.replace(" ", "")
            return n.replace(" ", "")
        return ""

    def metal_token_from_core_pairs(self, core_pairs: List[str], fallback_generic: List[str]) -> str:
        # Inspect core name|cas pairs to find a catalyst precursor, else use generic
        for p in core_pairs or []:
            nm, _, cs = p.partition("|")
            cs = cs.strip()
            if cs and cs in self.cas_map:
                gen = (self.cas_map[cs].get("GenericCore") or "").strip()
                if gen:
                    return gen
        # fallback to generic tag already computed
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

    def _compute_condition_core(self, row: Dict[str, Any]) -> str:
        core_gen = self._safe_json_list(row.get("CatalystCoreGeneric", "[]"))
        core_pairs = self._safe_json_list(row.get("CatalystCoreDetail", "[]"))
        lig_list = self._safe_json_list(row.get("Ligand", "[]"))

        metal = ""
        if self.taxonomy:
            metal = self.taxonomy.metal_token_from_core_pairs(core_pairs, core_gen)
        if not metal:
            metal = (core_gen[0] if core_gen else "").strip()

        lig_tok = ""
        if lig_list:
            first = lig_list[0]
            cas = ""
            nm = ""
            if "|" in first:
                nm, cas = (first.split("|", 1) + [""])[:2]
            else:
                nm = first
            if self.taxonomy:
                lig_tok = self.taxonomy.ligand_token_for(nm, cas)
            else:
                rec = (getattr(self, "cas_map", {}) or {}).get((cas or "").strip()) or {}
                lig_tok = (rec.get("Token") or rec.get("Abbreviation") or rec.get("Name") or nm).strip()
            lig_tok = lig_tok.replace(" ", "")

        # If no ligand but a metal core is present, try suggestion via PairingHelper
        if not lig_tok and metal and self.taxonomy:
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
                        lig_tok = (pick["ligand"]["abbr"] or "").strip().replace(" ", "")
            except Exception:
                pass

        return (f"{metal}/{lig_tok}" if metal and lig_tok else metal or lig_tok)

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

                    disp_core = self._compute_condition_core(row)

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
                        f.write(f"- SMILES: {r_smi} >> {p_smi}\n")
                    f.write("\n")
        except Exception:
            pass

    def generate_jsonl_export(self, rows, output_path: str, source_folder: str):
        out_lines = []
        for row in rows:
            rid = row.get("ReactionID", "")
            rtype = row.get("ReactionType", "")
            condition_core = self._compute_condition_core(row)

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

            smiles = {
                "reactants": row.get("ReactantSMILES", ""),
                "products": row.get("ProductSMILES", ""),
            }

            analysis_record = {
                "reaction_id": rid,
                "reaction_type": rtype,
                "condition_core": condition_core,
                "reagents": reagents,
                "solvents": solvents,
                "conditions": conditions,
                "smiles": smiles,
                "reference": row.get("Reference") or {},
            }
            out_lines.append(_json.dumps(analysis_record, ensure_ascii=False))

        try:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write("\n".join(out_lines) + ("\n" if out_lines else ""))
        except Exception:
            pass

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

    def __init__(self, folder_path: str, output_md_path: str, output_jsonl_path: str):
        super().__init__()
        self.folder_path = folder_path
        self.output_md_path = output_md_path
        self.output_jsonl_path = output_jsonl_path
        self.rdf_files = []

    def _emit(self, msg: str):
        """Emit progress message"""
        sig = getattr(self, 'progress', None)
        if sig:
            try:
                sig.emit(msg)
            except Exception:
                pass

    def _find_rdf_files(self) -> List[str]:
        """Find all RDF files in the specified folder"""
        rdf_files = []
        try:
            for file in os.listdir(self.folder_path):
                if file.lower().endswith('.rdf'):
                    full_path = os.path.join(self.folder_path, file)
                    if os.path.isfile(full_path):
                        rdf_files.append(full_path)
        except Exception as e:
            raise RuntimeError(f"Error scanning folder: {e}")
        
        return sorted(rdf_files)

    def _load_default_cas_maps(self) -> Dict[str, Dict[str, str]]:
        """Deprecated in this GUI workflow: use taxonomy-based mapping instead."""
        return {}

    def _load_taxonomy(self) -> _TaxonomyIndex:
        """Load compound taxonomies and build a CAS/name map for roles and tokens."""
        # Default base dir relative to repo
        here = os.path.dirname(os.path.abspath(__file__))
        tax_dir = os.path.join(os.path.dirname(here), 'data', 'compound_taxonomy')
        if not os.path.exists(tax_dir):
            # fallback to local data/compound_taxonomy relative to this file
            tax_dir = os.path.join(here, 'data', 'compound_taxonomy')
        try:
            idx = _TaxonomyIndex(tax_dir)
            return idx
        except Exception as e:
            raise RuntimeError(f"Failed to load taxonomy from {tax_dir}: {e}")

    def _reassign_reagent_roles_via_taxonomy(self, rows: List[Dict[str, Any]], taxonomy: _TaxonomyIndex) -> None:
        """Ensure ReagentRole aligns with Reagent list using taxonomy role lookup when possible."""
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
                        # retain prior heuristic role if any
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
        # Best-effort: copy JSONL to chemtools dataset dir for direct consumption
        export_jsonl_for_chemtools = None
        try:
            # Try plain import when running from this folder
            from chemtools_sink import export_jsonl_for_chemtools  # type: ignore
        except Exception:
            try:
                # Fallback: load via file path without requiring a package
                import importlib.util as _ilu, os as _os
                _p = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'chemtools_sink.py')
                _spec = _ilu.spec_from_file_location('chemtools_sink', _p)
                if _spec and _spec.loader:
                    _mod = _ilu.module_from_spec(_spec)
                    _spec.loader.exec_module(_mod)  # type: ignore
                    export_jsonl_for_chemtools = getattr(_mod, 'export_jsonl_for_chemtools', None)
            except Exception:
                export_jsonl_for_chemtools = None  # type: ignore
        if export_jsonl_for_chemtools:
            try:
                ok, msg = export_jsonl_for_chemtools(self.output_jsonl_path)
                if ok:
                    self._emit(f"Exported dataset JSONL for chemtools: {msg}")
                else:
                    self._emit(f"Note: could not export to chemtools dataset dir: {msg}")
            except Exception:
                # Silent failure to avoid blocking the primary outputs
                pass

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
                self.finished.emit(True, f"Successfully processed {len(self.rdf_files)} RDF files with {len(rows)} reactions.\nMarkdown: {self.output_md_path}\nJSONL: {self.output_jsonl_path}")
                
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
        form.addRow("RDF Folder:", folder_box)
        
        # Output file selection
        output_box = QtWidgets.QHBoxLayout()
        output_box.addWidget(self.output_md_edit)
        output_box.addWidget(self.btn_output_md)
        form.addRow("Output Markdown:", output_box)
        
        # Add note about JSONL
        note_label = QtWidgets.QLabel("Note: JSONL file will be automatically created alongside the Markdown file")
        note_label.setStyleSheet("font-style: italic; color: #666;")
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
        """Update the list of RDF files found in the selected folder"""
        folder_path = self.folder_edit.text().strip()
        self.file_list.clear()
        self.rdf_files = []
        
        if not folder_path or not os.path.isdir(folder_path):
            self.file_count_label.setText("No valid folder selected")
            self.btn_run.setEnabled(False)
            return
        
        try:
            # Find RDF files
            for file in os.listdir(folder_path):
                if file.lower().endswith('.rdf'):
                    full_path = os.path.join(folder_path, file)
                    if os.path.isfile(full_path):
                        self.rdf_files.append(full_path)
                        self.file_list.addItem(file)
            
            # Update UI
            count = len(self.rdf_files)
            if count == 0:
                self.file_count_label.setText("No RDF files found in this folder")
                self.btn_run.setEnabled(False)
            else:
                self.file_count_label.setText(f"Found {count} RDF file{'s' if count != 1 else ''}")
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
        # Name JSONL using upper folder names, e.g., Buchwald/2021-2025 -> Buchwald2021-2025.jsonl
        try:
            folder_path = self.folder_edit.text().strip()
            norm_folder = os.path.normpath(folder_path)
            parent_dir = os.path.basename(os.path.dirname(norm_folder)) or ""
            current_dir = os.path.basename(norm_folder) or ""
            import re as _re
            def _safe(s: str) -> str:
                s = _re.sub(r"\s+", "", s or "")  # remove spaces
                s = _re.sub(r"[^A-Za-z0-9_-]+", "", s)  # keep only safe chars
                return s
            base_name_no_ext = (_safe(parent_dir) + _safe(current_dir)) or "dataset"
            # Prefer to save next to the user-chosen Markdown directory; if none, use CWD
            md_hint = self.output_md_edit.text().strip()
            jsonl_dir = os.path.dirname(md_hint) if md_hint else os.getcwd()
            output_jsonl = os.path.join(jsonl_dir, base_name_no_ext + ".jsonl")
            # Ensure Markdown shares the same basename and location
            output_md = os.path.join(jsonl_dir, base_name_no_ext + ".md")
            try:
                # Reflect the computed Markdown path back into the UI
                self.output_md_edit.setText(output_md)
            except Exception:
                pass
        except Exception:
            # Fallback: derive both from any provided Markdown hint
            md_hint = self.output_md_edit.text().strip() or os.path.join(os.getcwd(), "dataset.md")
            output_md = os.path.splitext(md_hint)[0] + '.md'
            output_jsonl = os.path.splitext(output_md)[0] + '.jsonl'
        
        # Create worker and thread
        self.worker = RDFWorker(
            folder_path=self.folder_edit.text().strip(),
            output_md_path=output_md,
            output_jsonl_path=output_jsonl
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
