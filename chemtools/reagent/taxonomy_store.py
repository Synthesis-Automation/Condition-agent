"""
Taxonomy store and role heuristics for reagent classification.

Provides:
- TaxonomyStore: Manages reagent taxonomy data from JSON files
- RoleHeuristics: Infers reagent roles and families from names/synonyms
- Pattern matching for automated classification
"""

from __future__ import annotations

import copy
import datetime as dt
import json
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from chemtools.taxonomy import load_registry
from chemtools.taxonomy.models import ReagentFamily, ReagentRole

from .constants import DEFAULT_FAMILY_BY_ROLE, ROLE_FILES, ROLE_PRIORITY

if TYPE_CHECKING:
    from chemtools.taxonomy.registry import TaxonomyRegistry

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


def _safe_load_registry(candidates: Iterable[Optional[Path]]) -> Optional[TaxonomyRegistry]:
    """Attempt to load the taxonomy registry from a list of candidate roots."""
    for root in candidates:
        if root is None:
            try:
                return load_registry()
            except Exception:
                continue
        try:
            return load_registry(root)
        except Exception:
            continue
    return None


def _family_entry_from_registry(role: ReagentRole, family: ReagentFamily) -> Dict[str, Any]:
    """Convert a registry ReagentFamily into the legacy family schema dict."""
    entry: Dict[str, Any] = {
        "role": role.id,
        "role_label": role.name,
        "family": family.id,
        "family_id": family.id,
        "label": family.name or family.id,
        "definition": family.description or "",
        "include_smarts": list(family.include_smarts or []),
        "exclude_smarts": list(family.exclude_smarts or []),
        "required_props": dict(family.required_props or {}),
        "precedence": family.precedence if family.precedence is not None else role.priority,
        "keywords": list(family.keywords or []),
        "examples_pos": list(family.examples_pos or []),
        "examples_neg": list(family.examples_neg or []),
        "notes": family.notes or "",
        "metadata": copy.deepcopy(family.metadata or {}),
    }

    example_members = entry["metadata"].get("example_members")
    if example_members is None:
        example_members = [
            {"cas": cas, "source": "taxonomy_examples"}
            for cas in entry["examples_pos"]
            if cas
        ]
    entry["example_members"] = copy.deepcopy(example_members)
    return entry


def load_families_registry_entries(root: Optional[Path] = None) -> List[Dict[str, Any]]:
    """Return families registry entries derived from the unified taxonomy."""
    registry = _safe_load_registry([root, None])
    if registry is None:
        raise FileNotFoundError("Unable to load taxonomy registry for reagent families.")

    entries: List[Dict[str, Any]] = []
    for family in registry.reagent_families.values():
        role = registry.reagent_roles.get(family.role_id)
        if role is None:
            continue
        entries.append(_family_entry_from_registry(role, family))

    # Maintain deterministic ordering by role priority and family precedence.
    def _sort_key(entry: Dict[str, Any]) -> Tuple[int, str]:
        role_id = entry.get("role", "")
        role_obj = registry.reagent_roles.get(role_id)
        priority = role_obj.priority if role_obj else 100
        precedence = entry.get("precedence")
        if precedence is None:
            precedence = 100
        return (priority, precedence, entry.get("family", ""))

    entries.sort(key=_sort_key)
    return entries

# Manual patterns for specific reagents (higher priority)
MANUAL_FAMILY_PATTERNS: Dict[str, Tuple[str, str]] = {
    r"\bpeppsi\b": ("metal_precursor", "pd_peppsi_nhc"),
    r"\bgrubbs\b": ("metal_precursor", "ru_metathesis_grubbs"),
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
            priority = -ROLE_PRIORITY.get(role, 99)
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
        return DEFAULT_FAMILY_BY_ROLE.get(role)

    def _manual_match(self, texts: Sequence[str]) -> Optional[Tuple[str, str, List[str]]]:
        """Check manual patterns for exact matches."""
        joined = " ".join(t for t in texts if t)
        for pattern, role, family in self.manual_rules:
            if pattern.search(joined):
                return role, family, [pattern.pattern]
        return None


class TaxonomyStore:
    """Manages reagent taxonomy data sourced from the unified registry (with legacy fallback)."""

    def __init__(self, base_dir: Optional[Path] = None) -> None:
        """
        Initialize taxonomy store.

        Args:
            base_dir: Optional path to legacy taxonomy directory (e.g., data/compound_taxonomy).
        """
        self.base_dir = Path(base_dir) if base_dir is not None else None
        self.role_data: Dict[str, Dict[str, Any]] = {}
        self.family_lookup: Dict[str, Tuple[str, Dict[str, Any]]] = {}
        self.cas_index: Dict[str, Tuple[str, str, Dict[str, Any]]] = {}
        self.family_tokens: Dict[Tuple[str, str], Set[str]] = {}
        self.family_numeric_baseline: Dict[Tuple[str, str], Optional[Dict[str, Any]]] = {}
        self._load_all()

    def _load_all(self) -> None:
        """Load taxonomy data from the registry or legacy JSON files."""
        self.role_data.clear()
        self.family_lookup.clear()
        self.cas_index.clear()
        self.family_tokens.clear()
        self.family_numeric_baseline.clear()

        legacy_paths = []
        if self.base_dir is not None:
            legacy_paths = [self.base_dir / filename for filename in ROLE_FILES.values()]

        if legacy_paths and any(path.exists() for path in legacy_paths):
            self._load_from_legacy(legacy_paths)
            return

        registry = _safe_load_registry([self.base_dir, None])
        if registry is None:
            raise FileNotFoundError(
                "Unable to load reagent taxonomy from registry or legacy JSON files."
            )
        self._load_from_registry(registry)

    def _load_from_registry(self, registry: "TaxonomyRegistry") -> None:
        """Populate store data from the unified taxonomy registry."""
        for role_id, role in registry.reagent_roles.items():
            families: List[Dict[str, Any]] = []
            for family in registry.reagent_families.values():
                if family.role_id != role_id:
                    continue
                entry = _family_entry_from_registry(role, family)
                families.append(entry)
                family_id = entry["family_id"]
                self.family_lookup[family_id] = (role_id, entry)
                tokens = self._collect_family_tokens(entry)
                self.family_tokens[(role_id, family_id)] = tokens
                baseline: Optional[Dict[str, Any]] = None
                for member in entry.get("example_members", []):
                    cas = member.get("cas")
                    if cas:
                        self.cas_index[str(cas)] = (role_id, family_id, member)
                    numeric = member.get("numeric_features")
                    if baseline is None and numeric:
                        baseline = dict(numeric)
                self.family_numeric_baseline[(role_id, family_id)] = dict(baseline) if baseline else None

            # Preserve ordering by precedence within role
            families.sort(key=lambda item: (item.get("precedence", 100), item.get("family_id", "")))
            self.role_data[role_id] = {
                "role": role_id,
                "name": role.name,
                "description": role.description,
                "priority": role.priority,
                "default_family_id": role.default_family_id,
                "metadata": copy.deepcopy(role.metadata or {}),
                "families": families,
            }

    def _load_from_legacy(self, legacy_paths: List[Path]) -> None:
        """Fallback loader for legacy taxonomy JSON files."""
        for role, filename in ROLE_FILES.items():
            path = None
            if self.base_dir is not None:
                candidate = self.base_dir / filename
                if candidate.exists():
                    path = candidate
            if path is None:
                raise FileNotFoundError(
                    f"Expected taxonomy file '{filename}' for role '{role}' in legacy directory."
                )

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
        add(family.get("aliases", []))
        for member in family.get("example_members", []):
            add(member.get("name"))
            add(member.get("abbr"))
            add(member.get("synonyms", []))
            add(member.get("aliases", []))
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

    def add_entry(self, role: str, family_id: str, entry: Dict[str, Any]) -> None:
        """Add a new entry to a family."""
        from .taxonomy_utils import tokenize_all
        
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
        """Get file path for a role."""
        if self.base_dir is None:
            raise RuntimeError("Legacy taxonomy directory not configured for write access.")
        filename = ROLE_FILES.get(role)
        if not filename:
            raise KeyError(f"Unknown role '{role}'")
        return self.base_dir / filename

    def save_role(self, role: str) -> Path:
        """Save role data to file."""
        if self.base_dir is None:
            raise RuntimeError("Cannot persist taxonomy data without a legacy directory.")
        if role not in self.role_data:
            raise KeyError(f"No data cached for role '{role}'")
        data = self.role_data[role]
        data["updated"] = dt.date.today().isoformat()
        path = self.file_for_role(role)
        text = json.dumps(data, indent=2, ensure_ascii=False)
        path.write_text(text + "\n", encoding="utf-8")
        return path
