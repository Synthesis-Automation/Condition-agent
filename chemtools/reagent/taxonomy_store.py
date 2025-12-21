"""
Taxonomy store and role heuristics for reagent classification (v2).

Provides:
- TaxonomyStore: Manages reagent taxonomy data from v2 JSON files
- RoleHeuristics: Infers reagent roles and families from names/synonyms
- Pattern matching for automated classification
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from chemtools.taxonomy.reagent_v2 import (
    DEFAULT_REAGENT_V2_DIR,
    ReagentFamilyV2,
    ReagentRoleV2,
    ReagentTaxonomyV2,
    classify_reagent_v2,
)

from .constants import DEFAULT_FAMILY_BY_ROLE, ROLE_PRIORITY


# Role keyword patterns for automatic role detection
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
    "condensation_agent": [
        r"coupling",
        r"carbodiimide",
        r"uronium",
        r"phosphonium",
        r"\bt3p\b",
        r"\bbtffh\b",
        r"activated ester",
    ],
    "metal_catalyst": [
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


def _family_entry_from_v2(role: ReagentRoleV2, family: ReagentFamilyV2) -> Dict[str, Any]:
    """Convert a v2 family into the registry-style schema dict."""
    entry: Dict[str, Any] = {
        "role": role.id,
        "role_label": role.name,
        "family": family.id,
        "family_id": family.id,
        "label": family.name or family.id,
        "definition": family.description or "",
        "precedence": family.precedence,
        "keywords": list(family.allowlists.keywords),
        "notes": "",
        "examples_pos": [],
        "allowlists": {
            "cas": sorted(family.allowlists.cas),
            "names": sorted(family.allowlists.names),
            "keywords": list(family.allowlists.keywords),
        },
    }
    detect = family.detect
    if detect is not None:
        entry["detect"] = {
            "smarts": {
                "any": list(detect.smarts_any),
                "none": list(detect.smarts_none),
            }
        }
    return entry


def load_families_registry_entries(root: Optional[Path] = None) -> List[Dict[str, Any]]:
    """Return families registry entries derived from the v2 reagent taxonomy."""
    taxonomy = ReagentTaxonomyV2.from_path(root)
    entries: List[Dict[str, Any]] = []
    for family in taxonomy.families.values():
        role = taxonomy.roles.get(family.role_id)
        if role is None:
            continue
        entries.append(_family_entry_from_v2(role, family))

    def _sort_key(entry: Dict[str, Any]) -> Tuple[int, int, str]:
        role_id = entry.get("role", "")
        role_obj = taxonomy.roles.get(role_id)
        priority = role_obj.priority if role_obj else 100
        precedence = entry.get("precedence", 100)
        return (priority, precedence, entry.get("family", ""))

    entries.sort(key=_sort_key)
    return entries

# Manual patterns for specific reagents (higher priority)
MANUAL_FAMILY_PATTERNS: Dict[str, Tuple[str, str]] = {
    r"\bpeppsi\b": ("metal_catalyst", "pd_peppsi_nhc"),
    r"\bgrubbs\b": ("metal_catalyst", "ru_metathesis_grubbs"),
    r"\bdppf\b": ("ligand", "bidentate_diphosphines"),
    r"\bbtffh\b": ("condensation_agent", "acyl_halide_fluoride_generators"),
    r"\bt3p\b": ("condensation_agent", "organophosphorus_anhydrides"),
    r"\bcdi\b": ("condensation_agent", "imidazolide_formers"),
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

# Embedding field mapping for different reagent types
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
    "condensation_agent": (
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
    "metal_catalyst": (
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
    """Format a family field value for embedding text."""
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


def build_embedding_text(
    role: str,
    family: Dict[str, Any],
    entry: Dict[str, Any],
    synonyms: Sequence[str]
) -> str:
    """Build embedding text for a reagent entry."""
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


class RoleHeuristics:
    """Infer reagent roles and families from names and synonyms."""

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
        """
        Infer role and family from name and synonyms.
        
        Returns:
            (role, family_id, matching_tokens) or None
        """
        from .taxonomy_utils import tokenize_all, _sanitize_text
        
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
            priority = -self.store.role_priority(role)
            score = (len(matches), full_match, alias_match, priority)
            if not best or score > best[-1]:
                best = (role, family_id, matches, score)
        if best:
            role, family_id, matches, _score = best
            return role, family_id, matches
        return None

    def infer_role(self, name: str, synonyms: Sequence[str]) -> Optional[Tuple[str, str]]:
        """
        Infer role from keywords in name and synonyms.
        
        Returns:
            (role, matched_pattern) or None
        """
        texts = [name, *synonyms]
        joined = " ".join(t for t in texts if t)
        for role, patterns in self.role_patterns.items():
            for pattern in patterns:
                if pattern.search(joined):
                    return role, pattern.pattern
        return None

    def default_family_for_role(self, role: str) -> Optional[str]:
        """Get default family for a role."""
        data = self.store.role_data.get(role)
        if data and data.get("default_family_id"):
            return data.get("default_family_id")
        return DEFAULT_FAMILY_BY_ROLE.get(role)

    def _manual_match(self, texts: Sequence[str]) -> Optional[Tuple[str, str, List[str]]]:
        """Check manual patterns for exact matches."""
        joined = " ".join(t for t in texts if t)
        for pattern, role, family in self.manual_rules:
            if pattern.search(joined):
                return role, family, [pattern.pattern]
        return None


class TaxonomyStore:
    """Manages reagent taxonomy data from the v2 files."""

    def __init__(self, base_dir: Optional[Path] = None) -> None:
        """
        Initialize taxonomy store.

        Args:
            base_dir: Optional path to reagent taxonomy v2 directory.
        """
        self.base_dir = Path(base_dir) if base_dir is not None else DEFAULT_REAGENT_V2_DIR
        self.taxonomy: Optional[ReagentTaxonomyV2] = None
        self.role_data: Dict[str, Dict[str, Any]] = {}
        self.family_lookup: Dict[str, Tuple[str, Dict[str, Any]]] = {}
        self.cas_index: Dict[str, Tuple[str, str, Dict[str, Any]]] = {}
        self.family_tokens: Dict[Tuple[str, str], Set[str]] = {}
        self.family_numeric_baseline: Dict[Tuple[str, str], Optional[Dict[str, Any]]] = {}
        self._load_all()

    def _load_all(self) -> None:
        """Load taxonomy data from reagent taxonomy v2 files."""
        self.role_data.clear()
        self.family_lookup.clear()
        self.cas_index.clear()
        self.family_tokens.clear()
        self.family_numeric_baseline.clear()

        taxonomy = ReagentTaxonomyV2.from_path(self.base_dir)
        self.taxonomy = taxonomy
        self._load_from_v2(taxonomy)

    def _load_from_v2(self, taxonomy: ReagentTaxonomyV2) -> None:
        """Populate store data from reagent taxonomy v2."""
        for role_id, role in taxonomy.roles.items():
            families: List[Dict[str, Any]] = []
            for family in taxonomy.families.values():
                if family.role_id != role_id:
                    continue
                entry = _family_entry_from_v2(role, family)
                families.append(entry)
                family_id = entry["family_id"]
                self.family_lookup[family_id] = (role_id, entry)
                tokens = self._collect_family_tokens(entry)
                self.family_tokens[(role_id, family_id)] = tokens
                for cas in entry.get("allowlists", {}).get("cas", []):
                    if cas:
                        self.cas_index[str(cas)] = (role_id, family_id, {"cas": cas})
                self.family_numeric_baseline[(role_id, family_id)] = None

            families.sort(key=lambda item: (item.get("precedence", 100), item.get("family_id", "")))
            self.role_data[role_id] = {
                "role": role_id,
                "name": role.name,
                "description": role.description,
                "priority": role.priority,
                "default_family_id": role.default_family_id,
                "metadata": {},
                "families": families,
            }

    def _collect_family_tokens(self, family: Dict[str, Any]) -> Set[str]:
        """Collect all tokens from family for matching."""
        from .taxonomy_utils import _tokenize_text
        
        tokens: Set[str] = set()

        def add(text: Any) -> None:
            if isinstance(text, str):
                tokens.update(_tokenize_text(text))
            elif isinstance(text, (list, tuple, set)):
                for item in text:
                    add(item)

        add(family.get("family_id"))
        add(family.get("label"))
        add(family.get("definition"))
        add(family.get("keywords", []))
        allowlists = family.get("allowlists") or {}
        add(allowlists.get("names", []))
        add(allowlists.get("keywords", []))
        return tokens

    def family_token_overlap(self, role: str, family_id: str, tokens: Set[str]) -> bool:
        """Check if tokens overlap with family tokens."""
        if not tokens:
            return False
        family_tokens = self.family_tokens.get((role, family_id))
        if not family_tokens:
            return False
        return bool(family_tokens & tokens)

    def list_families(self, role: Optional[str] = None) -> List[Tuple[str, str, str]]:
        """
        List all families.
        
        Args:
            role: Optional role filter
            
        Returns:
            List of (role, family_id, label) tuples
        """
        roles = [role] if role else self.role_data.keys()
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

    def role_priority(self, role: str) -> int:
        """Return priority for a role."""
        data = self.role_data.get(role) or {}
        priority = data.get("priority")
        if priority is None:
            return ROLE_PRIORITY.get(role, 99)
        try:
            return int(priority)
        except (TypeError, ValueError):
            return ROLE_PRIORITY.get(role, 99)

    def suggest_families(self, role: str, tokens: Iterable[str], limit: int = 5) -> List[Dict[str, Any]]:
        """Suggest families matching the provided token set."""
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

    def role_for_family(self, family_id: str) -> Optional[str]:
        """Get role for a family."""
        entry = self.family_lookup.get(family_id)
        return entry[0] if entry else None

    def family_data(self, role: str, family_id: str) -> Dict[str, Any]:
        """Get family data dictionary."""
        entry = self.family_lookup.get(family_id)
        if not entry:
            raise KeyError(f"Unknown family '{family_id}'")
        entry_role, data = entry
        if entry_role != role:
            raise KeyError(f"Family '{family_id}' belongs to role '{entry_role}', not '{role}'")
        return data

    def numeric_baseline(self, role: str, family_id: str) -> Optional[Dict[str, Any]]:
        """Get numeric baseline for a family."""
        baseline = self.family_numeric_baseline.get((role, family_id))
        return dict(baseline) if baseline is not None else None

    def find_by_cas(self, cas: str) -> Optional[Tuple[str, str, Dict[str, Any]]]:
        """
        Find entry by CAS number.
        
        Returns:
            (role, family_id, entry_data) or None
        """
        entry = self.cas_index.get(cas)
        if entry:
            role, family_id, data = entry
            return role, family_id, data
        return None

    def classify_reagent(
        self,
        *,
        name: Optional[str] = None,
        cas: Optional[str] = None,
        smiles: Optional[str] = None,
    ) -> Optional[Dict[str, Any]]:
        """Classify reagent into a family/role using v2 rules."""
        if self.taxonomy is None:
            return None
        match = classify_reagent_v2(
            {"name": name, "cas": cas, "smiles": smiles},
            self.taxonomy.roles,
            self.taxonomy.families,
        )
        if match is None:
            return None
        return {
            "family_id": match.family_id,
            "role_id": match.role_id,
            "match_kind": match.match_kind,
            "precedence": match.precedence,
            "role_priority": match.role_priority,
        }

    def add_entry(self, role: str, family_id: str, entry: Dict[str, Any]) -> None:
        """Add a new entry to a family."""
        raise RuntimeError("Reagent taxonomy v2 is read-only and cannot be modified.")
        if self.family_numeric_baseline.get((role, family_id)) is None and entry.get("numeric_features"):
            self.family_numeric_baseline[(role, family_id)] = dict(entry["numeric_features"])

    def file_for_role(self, role: str) -> Path:
        """Get file path for a role."""
        if role not in self.role_data:
            raise KeyError(f"Unknown role '{role}'")
        return self.base_dir / "reagent_families.v2_cas.json"

    def save_role(self, role: str) -> Path:
        """Save role data to file."""
        raise RuntimeError("Reagent taxonomy v2 is read-only and cannot be modified.")
