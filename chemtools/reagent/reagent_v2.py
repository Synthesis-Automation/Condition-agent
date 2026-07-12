"""
Lightweight loader and classifier for reagent taxonomy v2.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import json
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from chemtools.core.rdkit import parse_smiles
from chemtools.core.smarts import compile_smarts

DEFAULT_REAGENT_V2_DIR = Path(__file__).resolve().parents[2] / "condition_registry" / "definitions"
_ROLES_FILE = "roles_families.v1.json"

_NAME_SANITIZE = re.compile(r"[^a-z0-9]+")


@dataclass(slots=True)
class ReagentRoleV2:
    id: str
    name: str
    priority: int
    default_family_id: Optional[str] = None
    description: Optional[str] = None


@dataclass(slots=True)
class ReagentAllowlists:
    cas: Set[str] = field(default_factory=set)
    names: Set[str] = field(default_factory=set)
    keywords: List[str] = field(default_factory=list)


@dataclass(slots=True)
class ReagentDetect:
    smarts_any: List[str] = field(default_factory=list)
    smarts_none: List[str] = field(default_factory=list)


@dataclass(slots=True)
class ReagentFamilyV2:
    id: str
    role_id: str
    name: str
    precedence: int
    description: Optional[str] = None
    notes: Optional[str] = None
    required_props: Dict[str, Any] = field(default_factory=dict)
    allowlists: ReagentAllowlists = field(default_factory=ReagentAllowlists)
    detect: Optional[ReagentDetect] = None


@dataclass(slots=True)
class ReagentMatchV2:
    family_id: str
    role_id: str
    match_kind: str
    precedence: int
    role_priority: int


class ReagentTaxonomyV2:
    """
    In-memory representation of the reagent taxonomy v2.
    """

    def __init__(
        self,
        roles: Dict[str, ReagentRoleV2],
        families: Dict[str, ReagentFamilyV2],
    ) -> None:
        self.roles = roles
        self.families = families

    @classmethod
    def from_path(cls, root: Optional[Path] = None) -> "ReagentTaxonomyV2":
        base = root or DEFAULT_REAGENT_V2_DIR
        roles, families = _load_taxonomy(base / _ROLES_FILE)
        return cls(roles=roles, families=families)

    def iter_roles(self) -> Iterable[ReagentRoleV2]:
        return self.roles.values()

    def iter_families(self, *, role_id: Optional[str] = None) -> Iterable[ReagentFamilyV2]:
        if role_id is None:
            return self.families.values()
        return (fam for fam in self.families.values() if fam.role_id == role_id)

    def get_role(self, role_id: str) -> Optional[ReagentRoleV2]:
        return self.roles.get(role_id)

    def get_family(self, family_id: str) -> Optional[ReagentFamilyV2]:
        return self.families.get(family_id)

    def classify(self, record: Dict[str, Optional[str]]) -> Optional[ReagentMatchV2]:
        return classify_reagent_v2(record, self.roles, self.families)


def _read_json(path: Path) -> Any:
    if not path.exists():
        raise FileNotFoundError(f"Reagent taxonomy file not found: {path}")
    for encoding in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            with path.open("r", encoding=encoding) as handle:
                return json.load(handle)
        except UnicodeDecodeError:
            continue
    text = path.read_text(encoding="utf-8", errors="replace")
    return json.loads(text)


def _parse_roles(payload: List[Dict[str, Any]]) -> Dict[str, ReagentRoleV2]:
    roles: Dict[str, ReagentRoleV2] = {}
    for entry in payload:
        role = ReagentRoleV2(
            id=entry["id"],
            name=entry.get("name", entry["id"]),
            priority=int(entry.get("priority", 100)),
            default_family_id=entry.get("default_family_id"),
            description=entry.get("description"),
        )
        roles[role.id] = role
    return roles


def _parse_families(payload: List[Dict[str, Any]]) -> Dict[str, ReagentFamilyV2]:
    families: Dict[str, ReagentFamilyV2] = {}
    for entry in payload:
        allowlists_raw = entry.get("allowlists") or {}
        allowlists = ReagentAllowlists(
            cas={str(cas) for cas in allowlists_raw.get("cas", []) if cas},
            names={_normalize_name(name) for name in allowlists_raw.get("names", []) if name},
            keywords=[str(tok).strip().lower() for tok in allowlists_raw.get("keywords", []) if tok],
        )
        detect = None
        detect_raw = entry.get("detect") or {}
        smarts_raw = detect_raw.get("smarts") or {}
        smarts_any = [str(pat) for pat in smarts_raw.get("any", []) if pat]
        smarts_none = [str(pat) for pat in smarts_raw.get("none", []) if pat]
        if smarts_any or smarts_none:
            detect = ReagentDetect(smarts_any=smarts_any, smarts_none=smarts_none)

        family = ReagentFamilyV2(
            id=entry["id"],
            role_id=entry["role_id"],
            name=entry.get("name", entry["id"]),
            precedence=int(entry.get("precedence", 100)),
            description=entry.get("description"),
            notes=entry.get("notes"),
            required_props=entry.get("required_props") or {},
            allowlists=allowlists,
            detect=detect,
        )
        families[family.id] = family
    return families


def _load_taxonomy(path: Path) -> Tuple[Dict[str, ReagentRoleV2], Dict[str, ReagentFamilyV2]]:
    payload = _read_json(path)
    if isinstance(payload, list):
        roles_payload = payload
        families_payload: List[Dict[str, Any]] = []
    elif isinstance(payload, dict):
        roles_payload = payload.get("roles") or []
        families_payload = payload.get("families") or []
    else:
        raise ValueError("Reagent taxonomy payload must be a list or object.")
    roles = _parse_roles(roles_payload)
    families = _parse_families(families_payload)
    return roles, families



def _normalize_name(name: str) -> str:
    normalized = _NAME_SANITIZE.sub(" ", (name or "").lower()).strip()
    return " ".join(normalized.split())


def _match_smarts(mol, smarts_any: Sequence[str], smarts_none: Sequence[str]) -> bool:
    if mol is None:
        return False
    if smarts_any:
        matched_any = False
        for pat in smarts_any:
            compiled = compile_smarts(pat, validate=False)
            if compiled and mol.HasSubstructMatch(compiled):
                matched_any = True
                break
        if not matched_any:
            return False
    for pat in smarts_none:
        compiled = compile_smarts(pat, validate=False)
        if compiled and mol.HasSubstructMatch(compiled):
            return False
    return True


def classify_reagent_v2(
    record: Dict[str, Optional[str]],
    roles: Dict[str, ReagentRoleV2],
    families: Dict[str, ReagentFamilyV2],
) -> Optional[ReagentMatchV2]:
    cas = str(record.get("cas") or "").strip()
    name = str(record.get("name") or "").strip()
    smiles = str(record.get("smiles") or "").strip()

    normalized_name = _normalize_name(name)
    mol = parse_smiles(smiles) if smiles else None

    candidates: List[ReagentMatchV2] = []

    for family in families.values():
        match_kind = None
        allowlists = family.allowlists

        if cas and cas in allowlists.cas:
            match_kind = "cas"
        elif mol is not None and family.detect is not None:
            if _match_smarts(mol, family.detect.smarts_any, family.detect.smarts_none):
                match_kind = "smarts"
        elif normalized_name:
            if allowlists.names and normalized_name in allowlists.names:
                match_kind = "name"
            elif allowlists.keywords:
                if all(tok in normalized_name for tok in allowlists.keywords):
                    match_kind = "keywords"

        if match_kind:
            role = roles.get(family.role_id)
            role_priority = role.priority if role else 100
            candidates.append(
                ReagentMatchV2(
                    family_id=family.id,
                    role_id=family.role_id,
                    match_kind=match_kind,
                    precedence=family.precedence,
                    role_priority=role_priority,
                )
            )

    if not candidates:
        return None

    strength_rank = {"cas": 0, "smarts": 1, "name": 2, "keywords": 3}

    def _sort_key(match: ReagentMatchV2) -> Tuple[int, int, int, str]:
        return (
            strength_rank.get(match.match_kind, 99),
            match.precedence,
            match.role_priority,
            match.family_id,
        )

    candidates.sort(key=_sort_key)
    return candidates[0]
