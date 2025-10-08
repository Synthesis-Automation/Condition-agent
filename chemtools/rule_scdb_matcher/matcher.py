from __future__ import annotations

"""Deterministic matcher supporting Scheme and Selector rule databases."""

from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from .ecn import essential_core_normalize
from .types import (
    MatchResult,
    RuleDB,
    SchemeConditionDB,
    SchemeEntry,
    SelectorRule,
    SelectorRuleDB,
)


@dataclass(slots=True)
class _Candidate:
    entry: SchemeEntry
    evaluation: Dict[str, Any]

    @property
    def sort_key(self) -> Tuple[int, int, str]:
        return (
            self.entry.priority,
            self.entry.specificity_score,
            self.entry.id,
        )


@dataclass(slots=True)
class _SelectorCandidate:
    rule: SelectorRule
    score: float
    evaluation: Dict[str, Any]

    @property
    def sort_key(self) -> Tuple[float, str]:
        return (self.score, self.rule.id)


def _parse_reactant_smiles(rxn_smiles: str) -> Tuple[str, ...]:
    if ">>" not in rxn_smiles:
        raise ValueError("Reaction SMILES must contain '>>'")
    reactant_part = rxn_smiles.split(">>", 1)[0]
    reactants = tuple(filter(None, (segment.strip() for segment in reactant_part.split("."))))
    if not reactants:
        raise ValueError("No reactants found in reaction SMILES")
    return reactants


_CONDITION_OPERATORS: Tuple[Tuple[str, str], ...] = (
    (" not in", "not_in"),
    (" in", "in"),
    (">=", ">="),
    ("<=", "<="),
    (">", ">"),
    ("<", "<"),
    ("!=", "!="),
    ("=", "=="),
)


def _ensure_sequence(value: Any) -> Tuple[Any, ...]:
    if value is None:
        return tuple()
    if isinstance(value, (list, tuple, set, frozenset)):
        return tuple(value)
    return (value,)


def _normalize_token(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return f"{float(value):g}"
    return str(value).strip()


def _value_to_float(value: Any) -> float | None:
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return None
    return None


def _parse_condition_key(key: str) -> Tuple[str, str]:
    cleaned = key.strip()
    for token, operator in _CONDITION_OPERATORS:
        if token in cleaned:
            path, _ = cleaned.split(token, 1)
            return path.strip(), operator
    return cleaned, "=="


def _resolve_feature_path(features: Mapping[str, Any], path: str) -> Any:
    current: Any = features
    for segment in path.split('.'):
        if isinstance(current, Mapping) and segment in current:
            current = current[segment]
        else:
            return None
    return current


def _value_equal(lhs: Any, rhs: Any) -> bool:
    left_num = _value_to_float(lhs)
    right_num = _value_to_float(rhs)
    if left_num is not None and right_num is not None:
        return abs(left_num - right_num) < 1e-9
    return _normalize_token(lhs) == _normalize_token(rhs)


def _evaluate_condition(actual: Any, operator: str, expected: Any) -> bool:
    if operator == "in":
        if actual is None:
            return False
        actual_tokens = {_normalize_token(item) for item in _ensure_sequence(actual)}
        expected_tokens = {_normalize_token(item) for item in _ensure_sequence(expected)}
        if not expected_tokens:
            return False
        return bool(actual_tokens & expected_tokens)

    if operator == "not_in":
        if actual is None:
            return True
        actual_tokens = {_normalize_token(item) for item in _ensure_sequence(actual)}
        expected_tokens = {_normalize_token(item) for item in _ensure_sequence(expected)}
        return not bool(actual_tokens & expected_tokens)

    if operator in {">=", ">", "<=", "<"}:
        actual_number = _value_to_float(actual)
        expected_number = _value_to_float(expected)
        if actual_number is None or expected_number is None:
            return False
        if operator == ">=":
            return actual_number >= expected_number
        if operator == ">":
            return actual_number > expected_number
        if operator == "<=":
            return actual_number <= expected_number
        return actual_number < expected_number

    if operator == "!=":
        if actual is None:
            return expected is not None
        actual_values = _ensure_sequence(actual)
        expected_values = _ensure_sequence(expected)
        return all(not _value_equal(item, candidate) for item in actual_values for candidate in expected_values)

    if actual is None:
        return False

    actual_values = _ensure_sequence(actual)
    expected_values = _ensure_sequence(expected)
    return any(_value_equal(item, candidate) for item in actual_values for candidate in expected_values)


def _evaluate_selector(rule: SelectorRule, features: Mapping[str, Any]) -> Tuple[bool, List[Dict[str, Any]]]:
    condition_results: List[Dict[str, Any]] = []
    matched = True
    for raw_key, expected in rule.conditions.items():
        path, operator = _parse_condition_key(raw_key)
        actual = _resolve_feature_path(features, path)
        passed = _evaluate_condition(actual, operator, expected)
        if not passed:
            matched = False
        condition_results.append(
            {
                "key": raw_key,
                "path": path,
                "operator": operator,
                "expected": expected,
                "actual": actual,
                "passed": passed,
            }
        )
    return matched, condition_results


def _selector_base_score(rule: SelectorRule) -> float:
    base = rule.rank_hint if rule.rank_hint is not None else 0.5
    return float(base + 0.001 * len(rule.conditions))


def _pick_best_selector(candidates: Sequence[_SelectorCandidate]) -> _SelectorCandidate | None:
    if not candidates:
        return None
    return sorted(candidates, key=lambda item: item.sort_key, reverse=True)[0]


def _to_plain(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _to_plain(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_to_plain(item) for item in value]
    return value


_LG_CLASS_PATTERNS = {
    "Cl": tuple(filter(None, (Chem.MolFromSmarts("[Cl]"),))),
    "Br": tuple(filter(None, (Chem.MolFromSmarts("[Br]"),))),
    "I": tuple(filter(None, (Chem.MolFromSmarts("[I]"),))),
    "OTf": tuple(filter(None, (Chem.MolFromSmarts("OS(=O)(=O)C(F)(F)F"),))),
    "OMs": tuple(filter(None, (Chem.MolFromSmarts("OS(=O)(=O)C"),))),
    "OFs": tuple(filter(None, (Chem.MolFromSmarts("OS(=O)(=O)F"),))),
}

_AMINE_CLASS_SMARTS: Dict[str, Tuple[str, ...]] = {
    "aniline_primary": ("[NH2]-[c]",),
    "aniline_secondary": ("[NH](-[c])-[c]",),
    "alkyl_primary": ("[NH2]-[CX4]",),
    "alkyl_secondary": ("[NH]([CX4])[CX4]",),
    "amide": ("[NX3;H1,H0]C(=O)[O,N]",),
    "sulfonamide": ("[NX3;H1]S(=O)(=O)[O,N,C,F]", "[NX3;H0]S(=O)(=O)[O,N,C,F]"),
    "carbamate": ("[NX3;H1]C(=O)O", "[NX3;H0]C(=O)O"),
}

_AMINE_CLASS_PATTERNS: Dict[str, Tuple[Chem.Mol, ...]] = {
    label: tuple(filter(None, (Chem.MolFromSmarts(smarts) for smarts in patterns)))
    for label, patterns in _AMINE_CLASS_SMARTS.items()
}

_ELECTRON_POOR_PATTERNS: Tuple[Chem.Mol, ...] = tuple(
    filter(
        None,
        (
            Chem.MolFromSmarts("[CX3](=O)[O,N,C]"),
            Chem.MolFromSmarts("[NX3+](=O)[O-]"),
            Chem.MolFromSmarts("C#N"),
            Chem.MolFromSmarts("C(F)(F)F"),
            Chem.MolFromSmarts("S(=O)(=O)[O,C,F]"),
        ),
    )
)


def _detect_amine_classes(mols: Sequence[Mol]) -> set[str]:
    detected: set[str] = set()
    for mol in mols:
        for label, patterns in _AMINE_CLASS_PATTERNS.items():
            if not patterns:
                continue
            if any(mol.HasSubstructMatch(pattern) for pattern in patterns):
                detected.add(label)
    return detected


def _compile_feature_sets(normalization) -> Tuple[Dict[str, set[str]], Dict[str, float]]:
    sanitized = normalization.sanitized_mols

    lg_classes: set[str] = set()
    electronics: set[str] = set()
    electrophile_ring_hetero_count = 0.0
    electrophile_ortho_sub_count = 0.0
    electrophile_meta_sub_count = 0.0
    nucleophile_ortho_sub_count = 0.0
    nucleophile_ring_hetero_count = 0.0
    electrophile_is_vinyl = False
    electrophile_is_primary_alkyl = False
    nucleophile_is_vinyl = False

    for idx, mol in enumerate(sanitized):
        matched_lg = False
        for label, patterns in _LG_CLASS_PATTERNS.items():
            if not patterns:
                continue
            if any(mol.HasSubstructMatch(pattern) for pattern in patterns):
                lg_classes.add(label)
                matched_lg = True

        # Detect vinyl groups (C=C with leaving group)
        from rdkit import Chem
        vinyl_pattern = Chem.MolFromSmarts("[C:1]=[C:2]-[Br,I,Cl,F]")
        if vinyl_pattern and mol.HasSubstructMatch(vinyl_pattern):
            if idx == 0:  # First reactant is typically the electrophile
                electrophile_is_vinyl = True
        
        # Detect vinyl boron (nucleophile)
        vinyl_boron_pattern = Chem.MolFromSmarts("[C:1]=[C:2]-[B]")
        if vinyl_boron_pattern and mol.HasSubstructMatch(vinyl_boron_pattern):
            nucleophile_is_vinyl = True

        # Detect primary alkyl halides (sp3 carbon with leaving group)
        primary_alkyl_pattern = Chem.MolFromSmarts("[CH2:1]-[Br,I]")
        if primary_alkyl_pattern and mol.HasSubstructMatch(primary_alkyl_pattern):
            if idx == 0:
                electrophile_is_primary_alkyl = True

        halogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in {"Cl", "Br", "I"}]
        boron_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "B"]
        
        if halogen_atoms:
            ring_info = mol.GetRingInfo().AtomRings()
            for halogen in halogen_atoms:
                for neighbor in halogen.GetNeighbors():
                    if neighbor.GetAtomicNum() != 6 or not neighbor.GetIsAromatic():
                        continue
                    ring = next((ring for ring in ring_info if neighbor.GetIdx() in ring), None)
                    if ring:
                        hetero = sum(1 for idx_r in ring if mol.GetAtomWithIdx(idx_r).GetAtomicNum() != 6)
                        electrophile_ring_hetero_count = max(electrophile_ring_hetero_count, float(hetero))
                        ortho = 0
                        meta = 0
                        ring_list = list(ring)
                        c_idx = neighbor.GetIdx()
                        c_position = ring_list.index(c_idx)
                        
                        for nbr in neighbor.GetNeighbors():
                            if nbr.GetIdx() == halogen.GetIdx() or not nbr.GetIsAromatic():
                                continue
                            if nbr.GetTotalNumHs() == 0:
                                ortho += 1
                        
                        # Count meta substituents (2 positions away in ring)
                        for i in [-2, 2]:
                            meta_pos = (c_position + i) % len(ring_list)
                            meta_atom = mol.GetAtomWithIdx(ring_list[meta_pos])
                            if meta_atom.GetTotalNumHs() == 0 and meta_atom.GetIsAromatic():
                                meta += 1
                        
                        electrophile_ortho_sub_count = max(electrophile_ortho_sub_count, float(ortho))
                        electrophile_meta_sub_count = max(electrophile_meta_sub_count, float(meta))

        # Check nucleophile (boron partner) for ring features
        if boron_atoms:
            ring_info = mol.GetRingInfo().AtomRings()
            for boron in boron_atoms:
                for neighbor in boron.GetNeighbors():
                    if neighbor.GetAtomicNum() not in {6, 7} or not neighbor.GetIsAromatic():
                        continue
                    ring = next((ring for ring in ring_info if neighbor.GetIdx() in ring), None)
                    if ring:
                        hetero = sum(1 for idx_r in ring if mol.GetAtomWithIdx(idx_r).GetAtomicNum() != 6)
                        nucleophile_ring_hetero_count = max(nucleophile_ring_hetero_count, float(hetero))
                        ortho = 0
                        for nbr in neighbor.GetNeighbors():
                            if nbr.GetIdx() == boron.GetIdx() or not nbr.GetIsAromatic():
                                continue
                            if nbr.GetTotalNumHs() == 0:
                                ortho += 1
                        nucleophile_ortho_sub_count = max(nucleophile_ortho_sub_count, float(ortho))

        if matched_lg and any(mol.HasSubstructMatch(pattern) for pattern in _ELECTRON_POOR_PATTERNS):
            electronics.add("electron-poor")

    set_features: Dict[str, set[str]] = {
        "electrophile.lg_class": lg_classes,
        "electrophile.electronics": electronics,
        "nucleophile.amine_class": _detect_amine_classes(sanitized),
    }

    # Add boolean features as sets (True/False as strings)
    if electrophile_is_vinyl:
        set_features["electrophile.is_vinyl"] = {True}
    else:
        set_features["electrophile.is_vinyl"] = {False}
    
    if electrophile_is_primary_alkyl:
        set_features["electrophile.is_primary_alkyl"] = {True}
    else:
        set_features["electrophile.is_primary_alkyl"] = {False}
    
    if nucleophile_is_vinyl:
        set_features["nucleophile.is_vinyl"] = {True}
    else:
        set_features["nucleophile.is_vinyl"] = {False}

    numeric_features: Dict[str, float] = {
        "electrophile.ring_hetero_count": electrophile_ring_hetero_count,
        "electrophile.ortho_sub_count": electrophile_ortho_sub_count,
        "electrophile.meta_sub_count": electrophile_meta_sub_count,
        "nucleophile.ortho_sub_count": nucleophile_ortho_sub_count,
        "nucleophile.ring_hetero_count": nucleophile_ring_hetero_count,
    }

    return set_features, numeric_features


def _check_numeric_requirement(value: float, requirement: Any) -> bool:
    if isinstance(requirement, (int, float)):
        return value == float(requirement)
    if isinstance(requirement, str):
        requirement = requirement.strip()
        if requirement.startswith(">="):
            try:
                threshold = float(requirement[2:].strip())
            except ValueError:
                return False
            return value >= threshold
        if requirement.startswith("<="):
            try:
                threshold = float(requirement[2:].strip())
            except ValueError:
                return False
            return value <= threshold
        try:
            expected_val = float(requirement)
        except ValueError:
            return False
        return value == expected_val
    return False


def _satisfies_set_requirement(observed: set, expected: Any) -> bool:
    """Check if observed set contains any of the expected values.
    
    Handles both string values and boolean values.
    """
    if isinstance(expected, list):
        # Convert expected values to the same type as observed
        expected_set = set()
        for value in expected:
            if isinstance(value, bool):
                expected_set.add(value)
            else:
                expected_set.add(str(value))
        return bool(expected_set & observed)
    
    # Single value expected
    if isinstance(expected, bool):
        return expected in observed
    return str(expected) in observed


def _requirements_satisfied(
    requirements: Mapping[str, Any] | None,
    set_features: Dict[str, set[str]],
    numeric_features: Dict[str, float],
) -> bool:
    if not requirements:
        return True
    for key, expected in requirements.items():
        if key in numeric_features:
            if not _check_numeric_requirement(numeric_features[key], expected):
                return False
            continue

        observed = set_features.get(key)
        if observed is None:
            return False
        if not _satisfies_set_requirement(observed, expected):
            return False
    return True


def _entry_applicable(
    entry: SchemeEntry,
    set_features: Dict[str, set[str]],
    numeric_features: Dict[str, float],
) -> bool:
    if not _requirements_satisfied(entry.feature_requirements, set_features, numeric_features):
        return False
    if entry.applies_to is None:
        return True
    for key, expected in entry.applies_to.items():
        observed = set_features.get(key, set())
        if not _satisfies_set_requirement(observed, expected):
            return False
    return True


def _match_entry(entry: SchemeEntry, masked_mols: Sequence[Mol], masked_sources: Sequence[int], kept_reactants: Sequence[str]) -> Tuple[bool, Dict[str, Any]]:
    evaluation: Dict[str, Any] = {
        "id": entry.id,
        "type": entry.type,
        "priority": entry.priority,
        "specificity": entry.specificity_score,
        "matches": [],
        "missing": [],
    }

    for compiled in entry.compiled_smarts:
        match_record = None
        for masked_idx, mol in enumerate(masked_mols):
            matches = mol.GetSubstructMatches(compiled.mol)
            if not matches:
                continue
            source_idx = masked_sources[masked_idx]
            match_record = {
                "smarts": compiled.pattern,
                "masked_index": masked_idx,
                "reactant_index": source_idx,
                "reactant_smiles": kept_reactants[source_idx],
                "match_count": len(matches),
                "atom_indices": list(matches[0]),
            }
            break
        if match_record is None:
            evaluation["missing"].append(compiled.pattern)
            return False, evaluation
        evaluation["matches"].append(match_record)

    return True, evaluation


def _pick_best(candidates: List[_Candidate]) -> _Candidate | None:
    if not candidates:
        return None
    return sorted(candidates, key=lambda cand: cand.sort_key, reverse=True)[0]


def _build_trace(normalization, evaluations: List[Dict[str, Any]], selected: Dict[str, Any] | None) -> Dict[str, Any]:
    anchor_hits = [
        {
            "anchor": hit.anchor_smarts,
            "reactant_index": hit.reactant_index,
            "atom_indices": list(hit.atom_indices),
        }
        for hit in normalization.anchor_hits
    ]
    return {
        "normalization": {
            "reactant_smiles": list(normalization.reactant_smiles),
            "kept_reactants": list(normalization.kept_reactants),
            "dropped_reactants": list(normalization.dropped_reactants),
            "smarts_bag": list(normalization.smarts_bag),
            "anchor_hits": anchor_hits,
        },
        "evaluations": evaluations,
        "selected": selected,
    }


def _match_scheme(db: SchemeConditionDB, rxn_smiles: str) -> MatchResult:
    """Return the best matching scheme/default for the provided reaction."""

    reactant_smiles = _parse_reactant_smiles(rxn_smiles)
    normalization = essential_core_normalize(reactant_smiles)
    set_features, numeric_features = _compile_feature_sets(normalization)

    evaluations: List[Dict[str, Any]] = []
    evaluation_by_id: Dict[str, Dict[str, Any]] = {}
    scheme_candidates: List[_Candidate] = []

    for entry in db.entries:
        if _entry_applicable(entry, set_features, numeric_features):
            matched, evaluation = _match_entry(
                entry,
                normalization.masked_mols,
                normalization.masked_source_indices,
                normalization.kept_reactants,
            )
        else:
            matched = False
            evaluation = {
                "id": entry.id,
                "type": entry.type,
                "priority": entry.priority,
                "specificity": entry.specificity_score,
                "matches": [],
                "missing": [],
                "matched": False,
                "filtered": True,
            }
        evaluation["matched"] = matched
        evaluations.append(evaluation)
        evaluation_by_id[entry.id] = evaluation
        if matched and entry.type == "scheme":
            scheme_candidates.append(_Candidate(entry=entry, evaluation=evaluation))

    selected_entry: SchemeEntry | None = None
    match_type: str | None = None

    best_scheme = _pick_best(scheme_candidates)
    if best_scheme:
        selected_entry = best_scheme.entry
        match_type = "scheme"
    else:
        default_candidates = [
            _Candidate(entry=entry, evaluation=evaluation_by_id[entry.id])
            for entry in db.entries
            if entry.type == "default_condition" and evaluation_by_id[entry.id]["matched"]
        ]
        best_default = _pick_best(default_candidates)
        if best_default:
            selected_entry = best_default.entry
            if selected_entry.reactant_smarts or selected_entry.applies_to or selected_entry.feature_requirements:
                match_type = "default"
            else:
                match_type = "global_default"

    selected_eval = evaluation_by_id[selected_entry.id] if selected_entry else None

    if selected_entry is None or match_type is None or selected_eval is None:
        trace = _build_trace(normalization, evaluations, None)
        raise RuntimeError("No matching scheme or default condition found")

    if match_type == "scheme":
        conditions = selected_entry.conditions or {}
    else:
        conditions = selected_entry.default_condition or {}

    trace = _build_trace(
        normalization,
        evaluations,
        {
            "id": selected_entry.id,
            "type": match_type,
            "priority": selected_entry.priority,
            "specificity": selected_entry.specificity_score,
        },
    )

    return MatchResult(
        reaction_type=db.reaction_type,
        match_type=match_type,
        entry_id=selected_entry.id,
        entry_name=selected_entry.name,
        priority=selected_entry.priority,
        conditions=conditions,
        trace=trace,
    )


def _match_selectors(db: SelectorRuleDB, features: Mapping[str, Any]) -> MatchResult:
    plain_features = _to_plain(features)
    evaluations: List[Dict[str, Any]] = []
    candidates: List[_SelectorCandidate] = []

    for rule in db.selectors:
        matched, condition_results = _evaluate_selector(rule, features)
        evaluation = {
            "id": rule.id,
            "rank_hint": rule.rank_hint,
            "matched": matched,
            "conditions": condition_results,
        }
        evaluations.append(evaluation)
        if matched:
            score = _selector_base_score(rule)
            candidates.append(_SelectorCandidate(rule=rule, score=score, evaluation=evaluation))

    selected_candidate = _pick_best_selector(candidates)

    if selected_candidate is not None:
        selected_rule = selected_candidate.rule
        conditions_payload = _to_plain(selected_rule.payload)
        if not isinstance(conditions_payload, dict):
            raise RuntimeError(f"Selector {selected_rule.id!r} payload must be a mapping")
        conditions = conditions_payload
        match_type = "selector"
        entry_id = selected_rule.id
        entry_name = selected_rule.raw.get("name") if isinstance(selected_rule.raw, Mapping) else None
        priority = int(round(selected_candidate.score * 100))
        selected_summary: Dict[str, Any] | None = {
            "id": selected_rule.id,
            "score": selected_candidate.score,
            "rank_hint": selected_rule.rank_hint,
            "evaluation": selected_candidate.evaluation,
        }
    else:
        conditions_payload = _to_plain(db.defaults)
        if not isinstance(conditions_payload, dict):
            raise RuntimeError("Selector defaults payload must be a mapping")
        conditions = conditions_payload
        match_type = "rule_default"
        entry_id = "DEFAULT"
        entry_name = None
        priority = 0
        selected_summary = None

    trace = {
        "features": plain_features,
        "evaluations": evaluations,
        "selected": selected_summary,
        "metadata": _to_plain(db.metadata),
        "defaults": _to_plain(db.defaults),
        "guardrails": [_to_plain(item) for item in db.guardrails],
        "priors": [_to_plain(item) for item in db.priors],
        "repairs": [_to_plain(item) for item in db.repairs],
    }

    return MatchResult(
        reaction_type=db.reaction_type,
        match_type=match_type,
        entry_id=entry_id,
        entry_name=entry_name,
        priority=priority,
        conditions=conditions,
        trace=trace,
    )


def match(
    db: RuleDB,
    rxn_smiles: str | None = None,
    *,
    features: Mapping[str, Any] | None = None,
) -> MatchResult:
    """Match a reaction against a rule database.

    Parameters
    ----------
    db:
        Loaded rule database (scheme or selector).
    rxn_smiles:
        Reaction SMILES required for SMARTS-driven scheme databases.
    features:
        Nested mapping of features required for selector databases.
    """

    if isinstance(db, SchemeConditionDB):
        if rxn_smiles is None:
            raise ValueError("rxn_smiles is required for SchemeConditionDB matching")
        return _match_scheme(db, rxn_smiles)

    if isinstance(db, SelectorRuleDB):
        if features is None:
            raise ValueError("features are required for SelectorRuleDB matching")
        return _match_selectors(db, features)

    raise TypeError(f"Unsupported rule database type: {type(db)!r}")
