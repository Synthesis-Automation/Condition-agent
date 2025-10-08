#!/usr/bin/env python3

"""Command-line reagent taxonomy generator that updates compound taxonomy files."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
from pathlib import Path
from urllib.parse import quote
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    import requests
except Exception:  # pragma: no cover
    requests = None  # type: ignore[assignment]

DEFAULT_RESOLVER_TIMEOUT = 6.0

ROLE_FILES: Dict[str, str] = {
    "acid": "taxonomy_acid.json",
    "additive": "taxonomy_additive.json",
    "ligand": "taxonomy_ligand.json",
    "metal_precursor": "taxonomy_catalysts_precursor.json",
    "base": "taxonomy_base.json",
    "coupling_reagent": "taxonomy_coupling_reagent.json",
    "oxidant": "taxonomy_oxidant.json",
    "reductant": "taxonomy_reductant.json",
    "solvent": "taxonomy_solvent.json",
}

DEFAULT_FAMILY_BY_ROLE: Dict[str, str] = {
    "acid": "mineral_acids",
    "additive": "quaternary_ammonium_ptc",
    "ligand": "trialkyl_triaryl_phosphines",
    "metal_precursor": "pd_ii_salts",
    "base": "tertiary_amines_aliphatic",
    "coupling_reagent": "carbodiimides",
    "oxidant": "o2_gas",
    "reductant": "metal_powders",
    "solvent": "hydrocarbons_aromatic",
    "other_reagent": "misc_general",
}

ROLE_PRIORITY: Dict[str, int] = {
    "ligand": 0,
    "metal_precursor": 1,
    "base": 2,
    "acid": 3,
    "coupling_reagent": 4,
    "oxidant": 5,
    "reductant": 6,
    "additive": 7,
    "solvent": 8,
    "other_reagent": 9,
}

ROLE_KEYWORDS_RAW: Dict[str, Sequence[str]] = {
    "ligand": [
        r"phos",
        r"ligand",
        r"imidazol",
        r"carbene",
        r"\bbpy\b",
        r"phenanthroline",
        r"nacnac",
        r"\bpyox\b",
    ],
    "base": [
        r"\bbase\b",
        r"amine",
        r"amide",
        r"carbonate",
        r"hydroxide",
        r"hydride",
        r"\bdbu\b",
        r"\bdbn\b",
        r"organolithium",
    ],
    "acid": [
        r"\bacid\b",
        r"triflic",
        r"trifluoracetic",
        r"\btfa\b",
        r"\btfoh\b",
        r"sulfonic",
        r"camphorsulfonic",
        r"\bhbf4\b",
        r"boric\s+acid",
    ],
    "solvent": [
        r"solvent",
        r"benzene",
        r"toluene",
        r"hexane",
        r"acetonitrile",
        r"\bdmf\b",
        r"\bdmso\b",
        r"water",
        r"alcohol",
    ],
    "other_reagent": [
        r"\bother\b",
        r"\bmisc\b",
        r"support",
        r"carrier",
    ],
    "coupling_reagent": [
        r"coupling",
        r"carbodiimide",
        r"uronium",
        r"phosphonium",
        r"\bt3p\b",
        r"\bbtffh\b",
        r"activated ester",
    ],
    "metal_precursor": [
        r"palladium",
        r"nickel",
        r"copper",
        r"cobalt",
        r"rhodium",
        r"iridium",
        r"ruthenium",
        r"iron",
        r"precatalyst",
    ],
    "reductant": [
        r"reductant",
        r"reducing agent",
        r"borohydride",
        r"hydride",
        r"silane",
        r"\bheh\b",
        r"\bnabh4\b",
        r"\blialh4\b",
    ],
    "additive": [
        r"additive",
        r"phase[- ]?transfer",
        r"\bptc\b",
        r"crown\s+ether",
        r"\bcrown\b",
        r"\btbaf\b",
        r"\btbab\b",
        r"\bquaternary\s+ammonium\b",
        r"silver\s+triflate",
        r"fluoride\s+source",
    ],
    "oxidant": [
        r"oxidant",
        r"peroxide",
        r"hypervalent",
        r"oxone",
        r"persulfate",
        r"permanganate",
        r"selectfluor",
        r"\bozone\b",
        r"\bo2\b",
    ],
}

MANUAL_FAMILY_PATTERNS: Dict[str, Tuple[str, str]] = {
    r"\bpeppsi\b": ("metal_precursor", "pd_peppsi_nhc"),
    r"\bgrubbs\b": ("metal_precursor", "ru_metathesis_grubbs"),
    r"\bdppf\b": ("ligand", "bidentate_diphosphines"),
    r"\bbtffh\b": ("coupling_reagent", "acyl_halide_fluoride_generators"),
    r"\bt3p\b": ("coupling_reagent", "organophosphorus_anhydrides"),
    r"\bcdi\b": ("coupling_reagent", "imidazolide_formers"),
    r"\bdbu\b": ("base", "amidine_guanidine_bases"),
    r"\bdbn\b": ("base", "amidine_guanidine_bases"),
    r"\b(?:potassium|sodium|cesium|caesium|lithium|rubidium|ammonium)\s+acetate\b": ("base", "inorg_carboxylates"),
    r"\b(?:k|na|cs|li|rb|nh4)oa?c\b": ("base", "inorg_carboxylates"),
    r"\bnabh4\b": ("reductant", "complex_hydrides"),
    r"\blialh4\b": ("reductant", "complex_hydrides"),
    r"\bhantzsch\b": ("reductant", "organic_reductants"),
    r"\btfoh\b": ("acid", "superacids"),
    r"\btriflic\b": ("acid", "superacids"),
    r"\btfa\b": ("acid", "halogenated_carboxylic_acids"),
    r"\btrifluoroacetic\b": ("acid", "halogenated_carboxylic_acids"),
    r"\btsoh\b": ("acid", "sulfonic_acids"),
    r"\bp-toluenesulfonic\b": ("acid", "sulfonic_acids"),
    r"\btbaf\b": ("additive", "fluoride_sources_activators"),
    r"\btbab\b": ("additive", "quaternary_ammonium_ptc"),
    r"\btetrabutylammonium\s+bromide\b": ("additive", "quaternary_ammonium_ptc"),
    r"\bcrown-?6\b": ("additive", "crown_ethers"),
    r"\boxone\b": ("oxidant", "oxone"),
    r"\bselectfluor\b": ("oxidant", "hypervalent_iodine"),
}

for _stream in ("stdout", "stderr"):
    try:
        getattr(sys, _stream).reconfigure(encoding="utf-8")
    except Exception:
        pass

TOKEN_PATTERN = re.compile(r"[a-z0-9]{3,}")
CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")
CAS_INLINE_PATTERN = re.compile(r"\d{2,7}-\d{2}-\d")

def _sanitize_text(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (text or "").lower())

def _tokenize_text(text: str) -> Set[str]:
    tokens: Set[str] = set()
    if not text:
        return tokens
    lower = (text or "").lower()
    compact = re.sub(r"[^a-z0-9]+", "", lower)
    if len(compact) >= 3:
        tokens.add(compact)
    tokens.update(tok for tok in TOKEN_PATTERN.findall(lower) if len(tok) >= 3)
    return tokens

def tokenize_all(texts: Iterable[str]) -> Set[str]:
    tokens: Set[str] = set()
    for text in texts:
        tokens.update(_tokenize_text(text))
    return tokens

def _valid_cas_checksum(cas: str) -> bool:
    digits = cas.replace("-", "")
    total = 0
    multiplier = 1
    for digit in reversed(digits[:-1]):
        total += int(digit) * multiplier
        multiplier += 1
    return total % 10 == int(digits[-1])

def normalize_cas(cas: str) -> str:
    digits = re.sub(r"[^0-9]", "", (cas or ""))
    if len(digits) < 5:
        raise ValueError(f"CAS '{cas}' does not have enough digits")
    prefix = str(int(digits[:-3]))
    mid = digits[-3:-1]
    check = digits[-1]
    normalized = f"{prefix}-{mid}-{check}"
    if not CAS_PATTERN.fullmatch(normalized):
        raise ValueError(f"CAS '{cas}' is not in a valid format")
    if not _valid_cas_checksum(normalized):
        raise ValueError(f"CAS '{cas}' failed checksum validation")
    return normalized

def _http_get_json(session: Any, url: str, timeout: float) -> Optional[Dict[str, Any]]:
    try:
        response = session.get(url, timeout=timeout)
    except Exception:
        return None
    if response.status_code != 200:
        return None
    try:
        return response.json()
    except Exception:
        return None

def _http_get_text(session: Any, url: str, timeout: float) -> Optional[str]:
    try:
        response = session.get(url, timeout=timeout)
    except Exception:
        return None
    if response.status_code != 200:
        return None
    return response.text

def _normalized_cas_tokens(synonyms: Sequence[str]) -> Set[str]:
    """Return normalized CAS numbers detected within the provided synonyms."""
    tokens: Set[str] = set()
    for synonym in synonyms:
        if not synonym:
            continue
        inline_matches = CAS_INLINE_PATTERN.findall(synonym)
        if inline_matches:
            for match in inline_matches:
                try:
                    tokens.add(normalize_cas(match))
                except ValueError:
                    continue
            continue
        try:
            normalized = normalize_cas(synonym)
        except ValueError:
            continue
        tokens.add(normalized)
    return tokens

def _pubchem_cid_rn_map(session: Any, cids: Sequence[Any], timeout: float) -> Dict[Any, Set[str]]:
    if not cids:
        return {}
    cid_tokens: List[str] = []
    for cid in cids:
        try:
            cid_int = int(str(cid))
        except (TypeError, ValueError):
            continue
        if cid_int < 0:
            continue
        cid_tokens.append(str(cid_int))
    if not cid_tokens:
        return {}
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        + ",".join(cid_tokens)
        + "/xrefs/RN/JSON"
    )
    data = _http_get_json(session, url, timeout)
    mapping: Dict[Any, Set[str]] = {}
    info_list = data.get("InformationList", {}).get("Information", []) if data else []
    for block in info_list:
        cid_val = block.get("CID")
        rn_list = block.get("RN") or []
        normalized: Set[str] = set()
        for rn in rn_list:
            try:
                normalized.add(normalize_cas(rn))
            except ValueError:
                continue
        if normalized:
            mapping[cid_val] = normalized
    return mapping

def _resolve_via_pubchem(session: Any, cas: str, timeout: float) -> Optional[Dict[str, Any]]:
    base = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RN/{quote(cas)}"
    props_url = base + "/property/Title,IUPACName,IsomericSMILES,CanonicalSMILES/JSON"
    props = _http_get_json(session, props_url, timeout)
    entries = props.get("PropertyTable", {}).get("Properties", []) if props else []
    if not entries:
        return None

    syn_url = base + "/synonyms/JSON"
    syn_data = _http_get_json(session, syn_url, timeout)
    synonyms_by_cid: Dict[Any, List[str]] = {}
    if syn_data:
        info = syn_data.get("InformationList", {}).get("Information", [])
        for block in info:
            cid = block.get("CID")
            synonyms = block.get("Synonym", []) or []
            deduped_synonyms = dedupe_synonyms([syn for syn in synonyms if syn])
            if cid is not None:
                synonyms_by_cid[cid] = deduped_synonyms

    try:
        normalized_query = normalize_cas(cas)
    except ValueError:
        normalized_query = None

    selected_record = entries[0]
    selected_synonyms = synonyms_by_cid.get(selected_record.get("CID"), [])
    match_found = False
    if normalized_query:
        for record in entries:
            cid = record.get("CID")
            record_synonyms = synonyms_by_cid.get(cid, [])
            if not record_synonyms:
                continue
            if normalized_query in _normalized_cas_tokens(record_synonyms):
                selected_record = record
                selected_synonyms = record_synonyms
                match_found = True
                break
        if len(entries) > 1 and not match_found:
            rn_map = _pubchem_cid_rn_map(
                session, [record.get("CID") for record in entries], timeout
            )
            for record in entries:
                cid = record.get("CID")
                if cid is None:
                    continue
                rn_set = rn_map.get(cid, set())
                if normalized_query in rn_set:
                    selected_record = record
                    selected_synonyms = synonyms_by_cid.get(cid, [])
                    match_found = True
                    break
        if len(entries) > 1 and not match_found:
            return None

    smiles = selected_record.get("IsomericSMILES") or selected_record.get("CanonicalSMILES")
    primary_name = selected_record.get("Title") or selected_record.get("IUPACName")

    raw: List[str] = []
    if primary_name:
        raw.append(primary_name)
    raw.extend(selected_synonyms)
    deduped = dedupe_synonyms(raw)
    if not deduped and primary_name:
        deduped = [primary_name]
    if not deduped:
        return None
    if len(deduped) > 16:
        deduped = deduped[:16]

    name = primary_name or deduped[0]
    return {"name": name, "synonyms": deduped, "smiles": smiles}

def _resolve_via_cactus(session: Any, cas: str, timeout: float) -> Optional[Dict[str, Any]]:
    base = f"https://cactus.nci.nih.gov/chemical/structure/{quote(cas)}"
    names_text = _http_get_text(session, base + "/names", timeout)
    synonyms = []
    if names_text:
        synonyms = [line.strip() for line in names_text.splitlines() if line.strip()]
    deduped = dedupe_synonyms(synonyms)
    if not deduped:
        return None
    smiles_text = _http_get_text(session, base + "/smiles", timeout)
    smiles = (smiles_text.strip() if smiles_text else None) or None
    if len(deduped) > 16:
        deduped = deduped[:16]
    return {"name": deduped[0], "synonyms": deduped, "smiles": smiles}

def resolve_identity_from_cas(
    cas: str,
    *,
    timeout: float = DEFAULT_RESOLVER_TIMEOUT,
    session: Optional[Any] = None,
) -> Optional[Dict[str, Any]]:
    if requests is None:
        return None
    if not cas:
        return None
    own_session = session is None
    sess = session or requests.Session()
    if own_session:
        sess.headers.setdefault("User-Agent", "ConditionAgent/TaxonomyGenerator")
    try:
        identity = _resolve_via_pubchem(sess, cas, timeout)
        if identity:
            identity["source"] = "pubchem"
            return identity
        identity = _resolve_via_cactus(sess, cas, timeout)
        if identity:
            identity["source"] = "cactus"
            return identity
        return None
    finally:
        if own_session:
            sess.close()

class RoleHeuristics:
    def __init__(self, store: "TaxonomyStore") -> None:
        self.store = store
        self.manual_rules: List[Tuple[re.Pattern, str, str]] = [
            (re.compile(pattern, re.IGNORECASE), role, family)
            for pattern, (role, family) in MANUAL_FAMILY_PATTERNS.items()
        ]
        self.role_patterns: Dict[str, List[re.Pattern]] = {
            role: [re.compile(pat, re.IGNORECASE) for pat in patterns]
            for role, patterns in ROLE_KEYWORDS_RAW.items()
        }

    def infer_family(self, name: str, synonyms: Sequence[str]) -> Optional[Tuple[str, str, List[str]]]:
        texts = [name, *synonyms]
        manual = self._manual_match(texts)
        if manual:
            return manual
        candidate_tokens = tokenize_all(texts)
        sanitized_name = _sanitize_text(name)
        alias_tokens = {_sanitize_text(alias) for alias in synonyms if alias}
        best: Optional[Tuple[str, str, List[str], Tuple[int, int, int, int]]] = None
        for (role, family_id), family_tokens in self.store.family_tokens.items():
            matches = sorted(candidate_tokens & family_tokens)
            if not matches:
                continue
            full_match = 1 if sanitized_name and sanitized_name in family_tokens else 0
            alias_match = 1 if any(tok and tok in family_tokens for tok in alias_tokens) else 0
            priority = -ROLE_PRIORITY.get(role, 99)
            score = (len(matches), full_match, alias_match, priority)
            if not best or score > best[-1]:
                best = (role, family_id, matches, score)
        if best:
            role, family_id, matches, _score = best
            return role, family_id, matches
        return None

    def infer_role(self, name: str, synonyms: Sequence[str]) -> Optional[Tuple[str, str]]:
        texts = [name, *synonyms]
        joined = " ".join(t for t in texts if t)
        for role, patterns in self.role_patterns.items():
            for pattern in patterns:
                if pattern.search(joined):
                    return role, pattern.pattern
        return None

    def default_family_for_role(self, role: str) -> Optional[str]:
        return DEFAULT_FAMILY_BY_ROLE.get(role)

    def _manual_match(self, texts: Sequence[str]) -> Optional[Tuple[str, str, List[str]]]:
        joined = " ".join(t for t in texts if t)
        for pattern, role, family in self.manual_rules:
            if pattern.search(joined):
                return role, family, [pattern.pattern]
        return None

class TaxonomyStore:
    def __init__(self, base_dir: Path) -> None:
        self.base_dir = Path(base_dir)
        if not self.base_dir.exists():
            raise FileNotFoundError(f"Taxonomy directory {self.base_dir} does not exist")
        self.role_data: Dict[str, Dict[str, Any]] = {}
        self.family_lookup: Dict[str, Tuple[str, Dict[str, Any]]] = {}
        self.cas_index: Dict[str, Tuple[str, str, Dict[str, Any]]] = {}
        self.family_tokens: Dict[Tuple[str, str], Set[str]] = {}
        self.family_numeric_baseline: Dict[Tuple[str, str], Optional[Dict[str, Any]]] = {}
        self._load_all()

    def _load_all(self) -> None:
        for role, filename in ROLE_FILES.items():
            path = self.base_dir / filename
            if not path.exists():
                raise FileNotFoundError(f"Expected taxonomy file {path} for role '{role}'")
            data = json.loads(path.read_text(encoding="utf-8"))
            self.role_data[role] = data
            for family in data.get("families", []):
                family_id = family.get("family_id")
                if not family_id:
                    continue
                self.family_lookup[family_id] = (role, family)
                tokens = self._collect_family_tokens(family)
                self.family_tokens[(role, family_id)] = tokens
                baseline: Optional[Dict[str, Any]] = None
                for member in family.get("example_members", []):
                    cas = member.get("cas")
                    if cas:
                        self.cas_index[cas] = (role, family_id, member)
                    if baseline is None and member.get("numeric_features"):
                        baseline = dict(member["numeric_features"])
                self.family_numeric_baseline[(role, family_id)] = dict(baseline) if baseline else None

    def _collect_family_tokens(self, family: Dict[str, Any]) -> Set[str]:
        tokens: Set[str] = set()

        def add(text: Any) -> None:
            if isinstance(text, str):
                tokens.update(_tokenize_text(text))
            elif isinstance(text, (list, tuple, set)):
                for item in text:
                    add(item)

        add(family.get("family_id"))
        add(family.get("label"))
        add(family.get("aliases", []))
        for member in family.get("example_members", []):
            add(member.get("name"))
            add(member.get("abbr"))
            add(member.get("synonyms", []))
            add(member.get("aliases", []))
        return tokens

    def family_token_overlap(self, role: str, family_id: str, tokens: Set[str]) -> bool:
        if not tokens:
            return False
        family_tokens = self.family_tokens.get((role, family_id))
        if not family_tokens:
            return False
        return bool(family_tokens & tokens)

    def list_families(self, role: Optional[str] = None) -> List[Tuple[str, str, str]]:
        roles = [role] if role else ROLE_FILES.keys()
        result: List[Tuple[str, str, str]] = []
        for rl in roles:
            data = self.role_data.get(rl)
            if not data:
                continue
            for family in data.get("families", []):
                fid = family.get("family_id")
                label = family.get("label", "")
                if fid:
                    result.append((rl, fid, label))
        return sorted(result, key=lambda item: (item[0], item[1]))

    def role_for_family(self, family_id: str) -> Optional[str]:
        entry = self.family_lookup.get(family_id)
        return entry[0] if entry else None

    def family_data(self, role: str, family_id: str) -> Dict[str, Any]:
        entry = self.family_lookup.get(family_id)
        if not entry:
            raise KeyError(f"Unknown family '{family_id}'")
        entry_role, data = entry
        if entry_role != role:
            raise KeyError(f"Family '{family_id}' belongs to role '{entry_role}', not '{role}'")
        return data

    def numeric_baseline(self, role: str, family_id: str) -> Optional[Dict[str, Any]]:
        baseline = self.family_numeric_baseline.get((role, family_id))
        return dict(baseline) if baseline is not None else None

    def find_by_cas(self, cas: str) -> Optional[Tuple[str, str, Dict[str, Any]]]:
        entry = self.cas_index.get(cas)
        if entry:
            role, family_id, data = entry
            return role, family_id, data
        return None

    def add_entry(self, role: str, family_id: str, entry: Dict[str, Any]) -> None:
        family = self.family_data(role, family_id)
        members = family.setdefault("example_members", [])
        members.append(entry)
        members.sort(key=lambda m: (m.get("name") or m.get("abbr") or "").lower())
        self.cas_index[entry["cas"]] = (role, family_id, entry)
        tokens = tokenize_all([entry.get("name"), entry.get("abbr"), *entry.get("synonyms", [])])
        self.family_tokens.setdefault((role, family_id), set()).update(tokens)
        if self.family_numeric_baseline.get((role, family_id)) is None and entry.get("numeric_features"):
            self.family_numeric_baseline[(role, family_id)] = dict(entry["numeric_features"])

    def file_for_role(self, role: str) -> Path:
        filename = ROLE_FILES.get(role)
        if not filename:
            raise KeyError(f"Unknown role '{role}'")
        return self.base_dir / filename

    def save_role(self, role: str) -> Path:
        if role not in self.role_data:
            raise KeyError(f"No data cached for role '{role}'")
        data = self.role_data[role]
        data["updated"] = dt.date.today().isoformat()
        path = self.file_for_role(role)
        text = json.dumps(data, indent=2, ensure_ascii=False)
        path.write_text(text + "\n", encoding="utf-8")
        return path

EMBEDDING_FIELD_MAP: Dict[str, Sequence[Tuple[str, str]]] = {
    "ligand": (
        ("donor_set", "donor_set"),
        ("denticity", "denticity"),
        ("scaffold", "scaffold"),
        ("steric_bulk_band", "steric bulk band"),
        ("electron_donor_band", "electron donor band"),
        ("pi_acceptor_band", "pi acceptor band"),
        ("rigidity_band", "rigidity band"),
        ("use_cases", "use_cases"),
        ("best_for", "best_for"),
        ("avoid_when", "avoid_when"),
        ("ehs_green", "ehs"),
    ),
    "base": (
        ("acid_base_type", "acid_base_type"),
        ("strength_band", "strength_band"),
        ("nucleophilicity_band", "nucleophilicity"),
        ("steric_profile", "steric"),
        ("moisture_sensitivity", "moisture"),
        ("typical_solvent_classes", "compatible_solvent_classes"),
        ("use_cases", "use_cases"),
        ("best_for", "best_for"),
        ("avoid_when", "avoid_when"),
        ("ehs_green", "ehs"),
    ),
    "acid": (
        ("functional_roles", "functional_roles"),
        ("effect_bands", "effect_bands"),
        ("moisture_sensitivity", "moisture"),
        ("typical_solvent_classes", "solvent_classes"),
        ("use_cases", "use_cases"),
        ("best_for", "best_for"),
        ("avoid_when", "avoid_when"),
        ("ehs_green", "ehs"),
    ),
    "solvent": (
        ("proticity", "proticity"),
        ("polarity_band", "polarity"),
        ("coordinating_ability", "coordination"),
        ("use_cases", "use_cases"),
        ("best_for", "best_for"),
        ("avoid_when", "avoid_when"),
        ("ehs_green", "ehs"),
    ),
    "coupling_reagent": (
        ("activation_mode", "activation"),
        ("leaving_group", "leaving_group"),
        ("typical_bases", "bases"),
        ("typical_solvents", "solvents"),
        ("relative_reactivity_band", "reactivity"),
        ("racemization_risk_band", "racemization_risk"),
        ("water_compatibility", "water_compatibility"),
        ("use_cases", "use_cases"),
        ("best_for", "best_for"),
        ("avoid_when", "avoid_when"),
        ("ehs_green", "ehs"),
    ),
    "metal_precursor": (
        ("metal", "metal"),
        ("oxidation_states", "oxidation_states"),
        ("precursor_class", "class"),
        ("preferred_reactions", "preferred_reactions"),
        ("use_cases", "use_cases"),
        ("best_for", "best_for"),
        ("avoid_when", "avoid_when"),
        ("handling", "handling"),
        ("additional_ligand", "additional_ligand"),
        ("recommended_pairs", "recommended_pairs"),
        ("ehs_green", "ehs"),
    ),
    "reductant": (
        ("class", "class"),
        ("mechanism", "mechanism"),
        ("strength_index", "strength"),
        ("use_cases", "use_cases"),
        ("best_for", "best_for"),
        ("avoid_when", "avoid_when"),
        ("recommended_pairs", "recommended_pairs"),
        ("ehs_green", "ehs"),
    ),
    "additive": (
        ("functional_roles", "functional_roles"),
        ("effect_bands", "effect_bands"),
        ("moisture_sensitivity", "moisture"),
        ("typical_solvent_classes", "solvent_classes"),
        ("use_cases", "use_cases"),
        ("best_for", "best_for"),
        ("avoid_when", "avoid_when"),
        ("ehs_green", "ehs"),
    ),
    "oxidant": (
        ("class", "class"),
        ("mechanism", "mechanism"),
        ("strength_index", "strength"),
        ("use_cases", "use_cases"),
        ("best_for", "best_for"),
        ("avoid_when", "avoid_when"),
        ("recommended_pairs", "recommended_pairs"),
        ("ehs_green", "ehs"),
    ),
}

def _format_family_field(key: str, value: Any) -> str:
    if value is None:
        return ""
    if key == "ehs_green" and isinstance(value, dict):
        flags = [k for k, v in value.items() if v]
        return ", ".join(flags) if flags else "none"
    if isinstance(value, dict):
        if key == "additional_ligand":
            parts: List[str] = []
            usage = value.get("usage")
            if usage:
                parts.append(f"usage: {usage}")
            typical = value.get("typical")
            if typical:
                if isinstance(typical, (list, tuple, set)):
                    parts.append("typical: " + ", ".join(str(x) for x in typical))
                else:
                    parts.append(f"typical: {typical}")
            for subkey, subval in value.items():
                if subkey in {"usage", "typical"}:
                    continue
                if subval:
                    parts.append(f"{subkey}: {subval}")
            return "; ".join(parts)
        return "; ".join(f"{k}: {v}" for k, v in value.items() if v)
    if isinstance(value, (list, tuple, set)):
        return "; ".join(str(x) for x in value)
    return str(value)

def build_embedding_text(role: str, family: Dict[str, Any], entry: Dict[str, Any], synonyms: Sequence[str]) -> str:
    pieces: List[str] = [f"type: {role.upper()}"]
    label = family.get("label", family.get("family_id", ""))
    fid = family.get("family_id", "")
    if label:
        pieces.append(f"family: {label} ({fid})")
    else:
        pieces.append(f"family: {fid}")
    for key, alias in EMBEDDING_FIELD_MAP.get(role, []):
        if key not in family:
            continue
        val = _format_family_field(key, family[key])
        if val:
            pieces.append(f"{alias}: {val}")
    pieces.append(f"name: {entry['name']}")
    if entry.get("abbr"):
        pieces.append(f"abbr: {entry['abbr']}")
    pieces.append(f"CAS: {entry['cas']}")
    if synonyms:
        pieces.append("synonyms: " + "; ".join(synonyms))
    return " | ".join(pieces)

def build_entry(
    role: str,
    family: Dict[str, Any],
    cas: str,
    name: str,
    abbr: str,
    synonyms: List[str],
    smiles: Optional[str],
    numeric_features: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    entry: Dict[str, Any] = {
        "name": name,
        "abbr": abbr,
        "cas": cas,
        "synonyms": synonyms,
    }
    if numeric_features:
        entry["numeric_features"] = numeric_features
    if smiles:
        entry["smiles"] = smiles
    alias_candidates = [syn for syn in synonyms if syn.lower() not in {name.lower(), abbr.lower()}]
    if alias_candidates:
        entry["aliases"] = alias_candidates
    entry["embedding_text"] = build_embedding_text(role, family, entry, synonyms)
    return entry

def dedupe_synonyms(values: Iterable[str]) -> List[str]:
    result: List[str] = []
    seen: Set[str] = set()
    for value in values:
        if not value:
            continue
        clean = value.strip()
        if not clean:
            continue
        key = clean.lower()
        if key in seen:
            continue
        seen.add(key)
        result.append(clean)
    return result

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assign reagent role/family and update taxonomy.")
    parser.add_argument("--cas", help="CAS number for the reagent.")
    parser.add_argument("--name", help="Primary name for the reagent.")
    parser.add_argument("--abbr", help="Short abbreviation.")
    parser.add_argument("--synonym", action="append", default=[], help="Additional synonym (repeatable).")
    parser.add_argument("--role", choices=sorted(ROLE_FILES.keys()), help="Override detected role.")
    parser.add_argument("--family", help="Explicit family identifier to use.")
    parser.add_argument("--smiles", help="Optional SMILES string to store with the entry.")
    parser.add_argument("--taxonomy-dir", default="data/compound_taxonomy", help="Path to taxonomy directory.")
    parser.add_argument("--dry-run", action="store_true", help="Infer role/family but do not write files.")
    parser.add_argument("--list-families", action="store_true", help="Print known families and exit.")
    parser.add_argument("--allow-default-family", action="store_true", help="Allow fallback default family when inference fails.")
    parser.add_argument("--no-auto-resolve", action="store_true", help="Disable CAS-based identity lookup when name is missing.")
    parser.add_argument("--resolver-timeout", type=float, default=DEFAULT_RESOLVER_TIMEOUT, help="Timeout (seconds) for resolver HTTP requests.")
    parser.add_argument("--verbose", action="store_true", help="Emit extra debug details.")
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    store = TaxonomyStore(Path(args.taxonomy_dir))
    heuristics = RoleHeuristics(store)

    if args.list_families:
        for role in sorted(ROLE_FILES.keys()):
            print(f"[{role}]")
            for _role, fid, label in store.list_families(role):
                label_str = f" - {label}" if label else ""
                print(f"  {fid}{label_str}")
        return

    if not args.cas:
        raise SystemExit("--cas is required unless --list-families is used.")

    cas = normalize_cas(args.cas)
    resolved_identity: Optional[Dict[str, Any]] = None
    if not args.no_auto_resolve and not args.name:
        resolved_identity = resolve_identity_from_cas(cas, timeout=args.resolver_timeout)
        if args.verbose:
            if resolved_identity:
                print(f"# auto-resolved via {resolved_identity.get('source')}: {resolved_identity.get('name')}")
            else:
                print("# auto-resolve failed to supply a name")

    name = args.name or (resolved_identity.get("name") if resolved_identity else None)
    if not name:
        raise SystemExit("Unable to determine reagent name. Provide --name or use --no-auto-resolve to skip lookup.")

    abbr = args.abbr or name
    resolved_synonyms = resolved_identity.get("synonyms", []) if resolved_identity else []
    synonyms = dedupe_synonyms([name, abbr, *args.synonym, *resolved_synonyms])
    input_tokens = tokenize_all([name, *synonyms])
    role = args.role
    family_id = args.family
    used_default = False
    default_rejection_reason: Optional[str] = None
    family_reason: Optional[List[str]] = None
    role_reason: Optional[str] = None
    auto_resolve_source = resolved_identity.get("source") if resolved_identity else None
    resolved_smiles = resolved_identity.get("smiles") if resolved_identity else None

    if family_id:
        family_role = store.role_for_family(family_id)
        if not family_role:
            raise SystemExit(f"Unknown family '{family_id}'. Use --list-families to inspect available options.")
        if role and role != family_role:
            raise SystemExit(f"Provided role '{role}' conflicts with family '{family_id}' (expected role '{family_role}').")
        role = family_role

    inference = heuristics.infer_family(name, synonyms) if not family_id else None
    if inference:
        role, inferred_family, reason_tokens = inference
        family_id = family_id or inferred_family
        family_reason = reason_tokens

    if not role:
        role_inference = heuristics.infer_role(name, synonyms)
        if role_inference:
            role, pattern = role_inference
            role_reason = pattern

    if not family_id:
        if role:
            default_family = heuristics.default_family_for_role(role)
            if default_family and args.allow_default_family:
                if store.family_token_overlap(role, default_family, input_tokens):
                    family_id = default_family
                    used_default = True
                else:
                    tokens_sample = ', '.join(sorted(input_tokens)[:6]) or 'none'
                    family_tokens = store.family_tokens.get((role, default_family), set())
                    family_sample = ', '.join(sorted(family_tokens)[:6]) or 'none'
                    default_rejection_reason = (
                        f"default family '{default_family}' rejected: no token overlap "
                        f"(input tokens: {tokens_sample}; family tokens sample: {family_sample})"
                    )
        if not family_id:
            message = (
                "Unable to determine family. Provide --family explicitly or use --allow-default-family."
            )
            if default_rejection_reason:
                message += f" Automatic fallback was skipped because {default_rejection_reason}."
            raise SystemExit(message)

    if not role:
        role = store.role_for_family(family_id)
    if not role:
        raise SystemExit("Unable to determine role; please supply --role.")

    existing = store.find_by_cas(cas)
    if existing:
        existing_role, existing_family, data = existing
        result = {
            "cas": cas,
            "name": data.get("name"),
            "role": existing_role,
            "family_id": existing_family,
            "status": "exists",
        }
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return

    family_data = store.family_data(role, family_id)
    numeric = store.numeric_baseline(role, family_id)
    entry = build_entry(role, family_data, cas, name, abbr, synonyms, args.smiles or resolved_smiles, numeric)

    result = {
        "cas": cas,
        "name": name,
        "role": role,
        "family_id": family_id,
        "taxonomy_file": str(store.file_for_role(role)),
        "dry_run": args.dry_run,
        "used_default_family": used_default,
    }
    if auto_resolve_source:
        result["auto_resolve_source"] = auto_resolve_source
    if resolved_smiles and not args.smiles:
        result["smiles_source"] = auto_resolve_source or "resolver"
    if family_reason:
        result["family_tokens"] = family_reason
    if role_reason:
        result["role_pattern"] = role_reason

    if args.dry_run:
        result["status"] = "dry_run"
        result["entry_preview"] = entry
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return

    store.add_entry(role, family_id, entry)
    path = store.save_role(role)
    result["status"] = "written"
    result["written_to"] = str(path)
    print(json.dumps(result, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
