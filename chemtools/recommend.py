from __future__ import annotations

from typing import Dict, Any, List, Tuple, Optional
from collections import Counter

from .smiles import normalize_reaction
from .router import detect_family
try:
    # Optional rxn-insight integration
    from .reaction_type_detector import detect_reaction_type as _rxn_detect  # type: ignore
    _HAS_RXN_INSIGHT = True
except Exception:
    _HAS_RXN_INSIGHT = False
from .featurizers import molecular as feat_molecular
from . import precedent, constraints, explain

try:  # Optional role-aware featurization
    from chem_feats import featurize_mol as _role_featurize_mol  # type: ignore
    _HAS_ROLE_FEATS = True
except Exception:  # pragma: no cover - optional dependency may be absent
    _role_featurize_mol = None  # type: ignore
    _HAS_ROLE_FEATS = False


_FAMILY_ROLE_EXPECTATIONS: Dict[str, List[Dict[str, Any]]] = {
    "Ullmann_CN": [
        {"label": "electrophile", "role": "aryl_halide", "required": True},
        {"label": "nucleophile", "role": "amine", "required": True},
    ],
    "Buchwald_CN": [
        {"label": "electrophile", "role": "aryl_halide", "required": True},
        {"label": "nucleophile", "role": "amine", "required": True},
    ],
    "Ullmann_O": [
        {"label": "electrophile", "role": "aryl_halide", "required": True},
        {"label": "nucleophile", "role": "alcohol", "required": True},
    ],
}


def _expected_roles_for_family(family: str | None) -> List[Dict[str, Any]]:
    fam = (family or "").strip()
    if not fam:
        return []
    return list(_FAMILY_ROLE_EXPECTATIONS.get(fam, []))


def _role_featurization_for_reactants(
    family: str,
    reactants: List[str],
) -> Dict[str, Any] | None:
    if not _HAS_ROLE_FEATS or _role_featurize_mol is None:
        return None

    entries: List[Dict[str, Any]] = []
    for smi in reactants:
        if not smi:
            continue
        try:
            data = _role_featurize_mol(smi)
        except Exception:
            continue
        entry = dict(data)
        entry["smiles"] = smi
        vec = entry.get("vector")
        try:
            if hasattr(vec, "tolist"):
                entry["vector"] = vec.tolist()  # type: ignore[attr-defined]
        except Exception:
            pass
        entries.append(entry)

    expected = _expected_roles_for_family(family)
    assignments: Dict[str, Dict[str, Any]] = {}
    unused = list(range(len(entries)))

    for spec in expected:
        label = str(spec.get("label") or "")
        role_token = spec.get("role")
        assigned_idx: Optional[int] = None
        matched_via_mask = False
        if role_token:
            for idx in list(unused):
                masks = entries[idx].get("masks") or {}
                try:
                    present = bool(masks.get(role_token))
                except Exception:
                    present = False
                if present:
                    assigned_idx = idx
                    matched_via_mask = True
                    break
        if assigned_idx is None and unused:
            assigned_idx = unused[0]
        if assigned_idx is not None:
            try:
                unused.remove(assigned_idx)
            except ValueError:
                pass
        assignments[label] = {
            "role": role_token,
            "index": assigned_idx,
            "matched_via_mask": matched_via_mask,
        }
        if assigned_idx is not None and 0 <= assigned_idx < len(entries):
            entries[assigned_idx] = dict(entries[assigned_idx])
            entries[assigned_idx]["assigned_label"] = label or None
            entries[assigned_idx]["assigned_role"] = role_token
            entries[assigned_idx]["assignment_reason"] = "mask" if matched_via_mask else "fallback"

    # Ensure every entry has explicit assignment metadata
    for idx, entry in enumerate(entries):
        entry.setdefault("assigned_label", None)
        entry.setdefault("assigned_role", None)
        entry.setdefault("assignment_reason", "unassigned")
        entry.setdefault("index", idx)

    return {
        "available": True,
        "family": family,
        "roles": expected,
        "reactants": entries,
        "assignments": assignments,
    }


def _electrophile_class_text(lg: str, elec_class: str) -> str:
    lg_norm = (lg or "").strip().upper()
    eclass = (elec_class or "").strip().lower()
    if eclass == "aryl":
        if lg_norm == "CL":
            return "aryl chloride"
        if lg_norm == "BR":
            return "aryl bromide"
        if lg_norm == "I":
            return "aryl iodide"
        if lg_norm in {"OTF", "OTS", "OMES", "OSF"}:
            return "aryl sulfonate"
        if lg_norm == "UNK":
            return "aryl electrophile"
    if eclass in {"alkyl", "aliphatic"}:
        if lg_norm == "CL":
            return "alkyl chloride"
        if lg_norm == "BR":
            return "alkyl bromide"
        if lg_norm == "I":
            return "alkyl iodide"
        return "alkyl electrophile"
    if eclass == "alkenyl" or eclass == "vinyl":
        if lg_norm == "OTF":
            return "vinyl triflate"
        if lg_norm:
            return f"vinyl {lg_norm.lower()}"
        return "vinyl electrophile"
    if lg_norm:
        return f"electrophile ({lg_norm})"
    return "unknown electrophile"


def _nucleophile_class_text(nuc_class: str) -> str:
    cls = (nuc_class or "").strip().lower()
    mapping = {
        "aniline": "primary aniline",
        "amine_primary": "primary amine",
        "amine_secondary": "secondary amine",
        "amine": "amine",
        "phenol": "phenol",
        "indole": "indole",
        "amide_deactivated": "amide",
    }
    if cls in mapping:
        return mapping[cls]
    if not cls:
        return "unknown nucleophile"
    return cls.replace("_", " ")


def _compose_rule_features(
    family: str,
    features: Dict[str, Any],
    role_pack: Dict[str, Any] | None,
) -> Dict[str, Any]:
    lg = features.get("LG")
    elec_class = features.get("elec_class")
    nuc_class = features.get("nuc_class")
    rule_feats: Dict[str, Any] = {
        "family": family,
        "bin": features.get("bin"),
        "LG": lg,
        "elec_class": elec_class,
        "nuc_class": nuc_class,
        "n_basicity": features.get("n_basicity"),
        "steric_alpha": features.get("steric_alpha"),
        "ortho_count": features.get("ortho_count"),
        "para_EWG": features.get("para_EWG"),
        "heteroaryl": features.get("heteroaryl"),
        "electrophile": {
            "class": _electrophile_class_text(str(lg or ""), str(elec_class or "")),
            "lg": lg,
            "type": elec_class,
        },
        "nucleophile": {
            "class": _nucleophile_class_text(str(nuc_class or "")),
            "basicity": features.get("n_basicity"),
            "steric_alpha": features.get("steric_alpha"),
        },
    }
    if isinstance(role_pack, dict):
        assignments = role_pack.get("assignments") or {}
        reactants = role_pack.get("reactants") or []
        role_summary: Dict[str, Any] = {}
        for label, info in assignments.items():
            idx = info.get("index")
            smi: Optional[str] = None
            if isinstance(idx, int) and 0 <= idx < len(reactants):
                smi = reactants[idx].get("smiles")
            role_summary[label] = {
                "role": info.get("role"),
                "smiles": smi,
                "matched_via_mask": info.get("matched_via_mask"),
            }
        if role_summary:
            rule_feats["role_assignments"] = role_summary
    return rule_feats


def _pick_electrophile_nucleophile(reactants: List[str]) -> Tuple[str, str]:
    def is_electrophile(s: str) -> bool:
        t = (s or "").lower()
        return (
            ("br" in t) or ("cl" in t) or (" i" in t)
            or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
        )
    if not reactants:
        return "", ""
    if len(reactants) == 1:
        return reactants[0], ""
    r0, r1 = reactants[0], reactants[1]
    if is_electrophile(r0):
        return r0, r1
    if is_electrophile(r1):
        return r1, r0
    return r0, r1


def _median(vals: List[float]) -> float | None:
    xs = [float(v) for v in vals if isinstance(v, (int, float))]
    if not xs:
        return None
    xs.sort()
    n = len(xs)
    mid = n // 2
    if n % 2 == 1:
        return xs[mid]
    return 0.5 * (xs[mid - 1] + xs[mid])


def _pick_with_constraints(cands: List[str], rules: Dict[str, Any]) -> Tuple[str | None, Dict[str, Any]]:
    if not cands:
        return None, {"allowed": [], "blocked": []}
    if not rules:
        return cands[0], {"allowed": cands, "blocked": []}
    out = constraints.apply_filter(cands, rules)
    allowed = out.get("allowed") or []
    return (allowed[0] if allowed else None), out


def recommend_from_reaction(
    reaction: str,
    k: int = 25,
    relax: Dict[str, Any] | None = None,
    constraint_rules: Dict[str, Any] | None = None,
    *,
    family_override: Optional[str] = None,
    max_variants: int = 3,
) -> Dict[str, Any]:
    """Recommend conditions from a reaction SMILES.

    Returns dict with keys: input, family, features, bin, recommendation, alternatives, pack, reasons.
    """
    relax = dict(relax or {})
    max_variants = max(1, int(max_variants or 1))
    fam_override_clean = (family_override.strip() if family_override else None)

    # 1) Normalize and extract reactants
    norm = normalize_reaction(reaction)
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]

    # 2) Detect family (try rxn-insight first when available; fallback to rules)
    detection_source = "user_supplied" if fam_override_clean else "auto"
    fam = fam_override_clean or "Unknown"
    rxn_auto: Dict[str, Any] | None = None
    rxn_smiles_norm = norm.get("normalized") or reaction
    # Allow callers to disable rxn-insight via relax['use_rxn_insight']=False
    # or globally via CHEMTOOLS_DISABLE_RXN_INSIGHT=1
    use_rxn_insight = relax.get("use_rxn_insight") if isinstance(relax, dict) else None
    if use_rxn_insight is None:
        import os as _os
        env_off = (_os.environ.get("CHEMTOOLS_DISABLE_RXN_INSIGHT", "").strip().lower() in {"1", "true", "yes", "on"})
        use_rxn_insight = not env_off
    if bool(use_rxn_insight) and _HAS_RXN_INSIGHT:
        try:
            rxn_auto = _rxn_detect(rxn_smiles_norm)  # type: ignore[misc]
        except Exception:
            rxn_auto = None
    auto_family = None
    if rxn_auto and rxn_auto.get("success") and rxn_auto.get("mapped_family"):
        auto_family = str(rxn_auto.get("mapped_family") or "Unknown")

    rule_info = detect_family(reactants)
    rule_family = rule_info.get("family") or "Unknown"

    if not fam_override_clean:
        fam = auto_family or rule_family or "Unknown"
    else:
        fam = fam_override_clean or auto_family or rule_family or "Unknown"

    # 3) Featurize substrates (Ullmann featurizer also used as fallback)
    elec, nuc = _pick_electrophile_nucleophile(reactants)
    features: Dict[str, Any] = {}
    if fam == "Ullmann_CN":
        features = feat_molecular.featurize(elec, nuc)
    else:
        features = feat_molecular.featurize(elec, nuc)

    role_pack = _role_featurization_for_reactants(fam, reactants)
    rule_features = _compose_rule_features(fam, features, role_pack)
    # Ensure features are hashable for caching (drop nested role-aware block)
    if isinstance(features, dict) and "role_aware" in features:
        try:
            features = {k: v for k, v in features.items() if k != "role_aware"}
        except Exception:
            features.pop("role_aware", None)  # type: ignore

    # 4) Retrieve precedents (enable DRFP unless explicitly disabled)
    relax.setdefault("reaction_smiles", norm.get("normalized") or reaction)
    relax.setdefault("use_drfp", True)
    relax.setdefault("precompute_drfp", True)
    relax.setdefault("precompute_scope", "candidates")
    pack = precedent.knn(family=fam, features=features, k=int(k), relax=relax)

    precs: List[Dict[str, Any]] = list(pack.get("precedents") or [])
    support = int(pack.get("support") or len(precs))

    # 5) Core vote (Laplace smoothing)
    core_counts = Counter([str(p.get("core") or "") for p in precs if p.get("core")])
    labels = list(core_counts.keys())
    alpha = 1.0
    denom = sum(core_counts.values()) + alpha * max(1, len(labels))
    scores = {L: (core_counts.get(L, 0) + alpha) / denom for L in labels}
    chosen_core = max(scores, key=scores.get) if scores else None
    core_vote_share = (core_counts.get(chosen_core, 0) / max(1, sum(core_counts.values()))) if chosen_core else 0.0

    # 6) Choose base and solvent (conditioned on chosen core when possible)
    if chosen_core:
        group = [p for p in precs if str(p.get("core") or "") == chosen_core]
    else:
        group = precs
    bases = [str(p.get("base_uid") or "") for p in group if p.get("base_uid")]
    solvents = [str(p.get("solvent_uid") or "") for p in group if p.get("solvent_uid")]
    base_counts = Counter(bases)
    solv_counts = Counter(solvents)
    base_list = [b for b, _ in base_counts.most_common()] or [str(p.get("base_uid") or "") for p in precs if p.get("base_uid")]
    solv_list = [s for s, _ in solv_counts.most_common()] or [str(p.get("solvent_uid") or "") for p in precs if p.get("solvent_uid")]

    base_pick, base_filter = _pick_with_constraints(base_list, constraint_rules or {})
    solv_pick, solv_filter = _pick_with_constraints(solv_list, constraint_rules or {})
    # Fallback to overall precedents if filtered out or empty
    if (not base_pick) and precs:
        all_bases = [str(p.get("base_uid") or "") for p in precs if p.get("base_uid")]
        base_pick, base_filter = _pick_with_constraints(list(dict.fromkeys(all_bases)), constraint_rules or {})
    if (not solv_pick) and precs:
        all_solv = [str(p.get("solvent_uid") or "") for p in precs if p.get("solvent_uid")]
        solv_pick, solv_filter = _pick_with_constraints(list(dict.fromkeys(all_solv)), constraint_rules or {})

    # 7) T/time median from the same-core group; fallback to all precedents
    def nums(key: str, items: List[Dict[str, Any]]):
        return [p.get(key) for p in items if isinstance(p.get(key), (int, float))]

    T_med = _median(nums("T_C", group) or nums("T_C", precs))
    t_med = _median(nums("time_h", group) or nums("time_h", precs))

    # 8) Confidence (simple): share of top core among all with any core label, capped to [0.3, 0.95]
    conf = 0.95 * core_vote_share if support >= 5 else 0.5 * core_vote_share
    conf = max(0.3, min(0.95, conf))

    # 9) Reasons from precedents
    reasons_pack = explain.for_pack(pack, features)

    recommendation = {
        "core": chosen_core,
        "base_uid": base_pick,
        "solvent_uid": solv_pick,
        "T_C": T_med,
        "time_h": t_med,
        "confidence": round(float(conf), 3),
    }

    alternatives = {
        "cores": core_counts.most_common(3),
        "bases": base_counts.most_common(3),
        "solvents": solv_counts.most_common(3),
    }

    # 10) Build formatted output similar to tests/format_of_recommended_reactions.json
    def _nice_family_text(f: str) -> str:
        t = (f or "").lower()
        if t.startswith("ullmann"):
            return "Ullmann"
        if t.startswith("suzuki"):
            return "Suzuki"
        if t.startswith("amide"):
            return "Amide"
        return f or "Unknown"

    def _lookup(uid: str) -> Dict[str, Any]:
        try:
            from .properties import lookup  # lazy import
            res = lookup(uid)
            if isinstance(res, dict) and res.get("found") and isinstance(res.get("record"), dict):
                return res["record"]
        except Exception:
            pass
        return {"uid": uid}

    # Reactant chemicals (starting materials)
    reactants_chems: List[Dict[str, Any]] = []
    for r in (norm.get("reactants") or []):
        smi = r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or ""
        reactants_chems.append({
            "name": None,
            "cas": None,
            "smiles": smi or None,
            "equivalents": None,
            "role": "starting_material",
        })

    # Catalyst system from precedents of chosen core (if available)
    def _cat_items() -> List[Dict[str, Any]]:
        items: List[Dict[str, Any]] = []
        src = next((p for p in group if p.get("catalyst") or p.get("full_system")), None)
        fs = None
        if src:
            cat = src.get("catalyst") or {}
            fs = src.get("full_system") or (cat.get("full_system") if isinstance(cat, dict) else None)
        # Fallback: empty
        if not isinstance(fs, list):
            return items
        def role_for(name: str) -> str:
            n = (name or "").lower()
            if any(tok in n for tok in ["pd", "palladium", "cu", "copper", "ni", "nickel", "ru", "ruthenium", "co", "cobalt"]):
                return "metal_precursor"
            return "ligand"
        for it in fs:
            nm = (it or {}).get("name")
            cs = (it or {}).get("cas")
            items.append({
                "name": nm,
                "cas": cs,
                "smiles": None,
                "equivalents": None,
                "role": role_for(str(nm or "")),
            })
        return items

    # Base and solvent chemicals using lookups
    def _chemical_payload(uid: Optional[str], role: str) -> Dict[str, Any] | None:
        if not uid:
            return None
        rec = _lookup(uid)
        return {
            "name": rec.get("name") or rec.get("token") or rec.get("uid") or uid,
            "abbreviation": rec.get("token") or None,
            "cas": rec.get("uid") or uid,
            "smiles": rec.get("smiles"),
            "equivalents": None,
            "role": role,
        }

    base_item = _chemical_payload(base_pick, "base")
    solv_item = _chemical_payload(solv_pick, "solvent")

    combo_counts = Counter((str(p.get("base_uid") or ""), str(p.get("solvent_uid") or "")) for p in group)
    overall_combo_counts = Counter((str(p.get("base_uid") or ""), str(p.get("solvent_uid") or "")) for p in precs)

    def _matching_precedents(b_key: str, s_key: str, src: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        return [
            p for p in src
            if str(p.get("base_uid") or "") == b_key and str(p.get("solvent_uid") or "") == s_key
        ]

    cat_items = _cat_items()
    cond_text = {
        "temperature": (f"{int(T_med)} °C" if isinstance(T_med, (int, float)) else None),
        "time": (f"{t_med} h" if isinstance(t_med, (int, float)) else None),
        "atmosphere": None,
    }

    core_group_size = len(group) if group else 0

    def _build_variant(b_uid: Optional[str], s_uid: Optional[str], rank: int) -> Dict[str, Any]:
        b_key = b_uid or ""
        s_key = s_uid or ""
        chems = list(reactants_chems) + list(cat_items)
        base_payload = _chemical_payload(b_uid, "base")
        solvent_payload = _chemical_payload(s_uid, "solvent")
        if base_payload:
            chems.append(base_payload)
        if solvent_payload:
            chems.append(solvent_payload)
        support_count = combo_counts.get((b_key, s_key), 0)
        if support_count == 0:
            support_count = overall_combo_counts.get((b_key, s_key), 0)
        denom = core_group_size if core_group_size else len(precs)
        support_fraction = (support_count / denom) if denom else 0.0
        matched = _matching_precedents(b_key, s_key, group or precs)
        precedent_examples = [
            {
                "reaction_id": p.get("reaction_id"),
                "reference": p.get("reference"),
                "yield_pct": p.get("yield_pct"),
                "core": p.get("core"),
            }
            for p in matched[:3]
            if p.get("reaction_id")
        ]
        summary = {
            "rank": rank,
            "core": chosen_core,
            "base": base_payload,
            "solvent": solvent_payload,
            "confidence": round(float(conf), 3),
            "support": {
                "count": support_count,
                "fraction_core": round(float(support_fraction), 3) if support_fraction else 0.0,
                "reference_population": core_group_size if core_group_size else len(precs),
            },
            "precedents": precedent_examples,
        }
        variant = {
            "rank": rank,
            "reaction": {"smiles": norm.get("normalized") or reaction},
            "chemicals": chems,
            "conditions": cond_text,
            "summary": summary,
            "combo": {"base_uid": b_uid, "solvent_uid": s_uid},
        }
        return variant

    combos: List[Tuple[str, str]] = []
    seen_combos: set[Tuple[str, str]] = set()

    def _add_combo(b: Optional[str], s: Optional[str]) -> None:
        key = (b or "", s or "")
        if key in seen_combos:
            return
        seen_combos.add(key)
        combos.append(key)

    _add_combo(base_pick, solv_pick)
    for combo, _ in combo_counts.most_common():
        _add_combo(combo[0], combo[1])
    for combo, _ in overall_combo_counts.most_common():
        _add_combo(combo[0], combo[1])
    for b, _ in base_counts.most_common():
        _add_combo(b, solv_pick)
    for s, _ in solv_counts.most_common():
        _add_combo(base_pick, s)
    if not combos:
        _add_combo(None, None)

    variants: List[Dict[str, Any]] = []
    for combo in combos:
        if len(variants) >= max_variants:
            break
        b_key, s_key = combo
        variants.append(_build_variant(b_key or None, s_key or None, len(variants) + 1))


    from datetime import datetime, timezone
    formatted = {
        "meta": {
            "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            "analysis_type": "recommendation",
            "status": "success",
        },
        "input": {
            "reaction_smiles": norm.get("normalized") or reaction,
            "provided_reaction_type": fam_override_clean,
            "selected_reaction_type": _nice_family_text(fam),
        },
        "detection": {
            "reaction_type": _nice_family_text(fam),
            "source": detection_source,
            "provided_reaction_type": fam_override_clean,
            "auto": (
                {
                    "rxn_insight_available": bool(rxn_auto.get("available")) if isinstance(rxn_auto, dict) else False,
                    "rxn_insight_class": (rxn_auto.get("rxn_class") if isinstance(rxn_auto, dict) else None),
                    "rxn_insight_name": (rxn_auto.get("rxn_name") if isinstance(rxn_auto, dict) else None),
                    "rxn_insight_confidence": (rxn_auto.get("confidence") if isinstance(rxn_auto, dict) else None),
                    "status": ("success" if (isinstance(rxn_auto, dict) and rxn_auto.get("success")) else "fallback"),
                    "reaction_type": _nice_family_text(auto_family) if auto_family else None,
                }
                if rxn_auto is not None
                else None
            ),
            "rule_based": {
                "reaction_type": _nice_family_text(rule_family),
                "raw": rule_info,
            },
        },
        "recommended_conditions": variants,
    }

    formatted["starting_materials"] = {
        "role_featurization": role_pack,
        "rule_features": rule_features,
    }

    return {
        "input_reaction": reaction,
        "family": fam,
        "features": features,
        "bin": features.get("bin"),
        "recommendation": recommendation,
        "alternatives": alternatives,
        "precedent_pack": pack,
        "reasons": reasons_pack.get("reasons"),
        "filters": {"base": base_filter, "solvent": solv_filter},
        "formatted": formatted,
        "role_featurization": role_pack,
        "rule_features": rule_features,
    }


def recommend_conditions_structured(
    reaction: str,
    reaction_type: Optional[str] = None,
    *,
    k: int = 50,
    limit: int = 5,
    relax: Optional[Dict[str, Any]] = None,
    constraints: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Structured condition recommendations for API/UI consumers."""
    limit = max(1, int(limit or 1))
    cfg_relax = dict(relax or {})
    result = recommend_from_reaction(
        reaction=reaction,
        k=int(k or 50),
        relax=cfg_relax,
        constraint_rules=constraints or {},
        family_override=reaction_type,
        max_variants=limit,
    )

    formatted = dict(result.get("formatted") or {})
    recommendations = list(formatted.get("recommended_conditions") or [])
    recommendations = recommendations[:limit]
    for idx, rec in enumerate(recommendations, start=1):
        rec.setdefault("rank", idx)
        summary = rec.get("summary")
        if isinstance(summary, dict):
            summary.setdefault("rank", idx)

    detection = dict(formatted.get("detection") or {})
    if reaction_type and not detection.get("source"):
        detection["source"] = "user_supplied"
    detection.setdefault("source", detection.get("source") or "auto")
    detection.setdefault("provided_reaction_type", reaction_type)

    meta = dict(formatted.get("meta") or {})
    meta.setdefault("strategy", "precedent_knn")
    meta["result_count"] = len(recommendations)

    pack = result.get("precedent_pack") or {}
    precedents = list(pack.get("precedents") or [])
    top_precedents = [
        {
            "reaction_id": p.get("reaction_id"),
            "core": p.get("core"),
            "yield_pct": p.get("yield_pct"),
        }
        for p in precedents[:10]
        if p.get("reaction_id")
    ]
    core_family = result.get("family")
    core_support = len([p for p in precedents if p.get("core") == core_family])
    precedent_summary = {
        "total_considered": len(precedents),
        "core_family": core_family,
        "core_support": core_support,
        "top_precedents": top_precedents,
    }

    return {
        "meta": meta,
        "input": formatted.get("input"),
        "detection": detection,
        "recommendations": recommendations,
        "alternatives": result.get("alternatives"),
        "precedents": precedent_summary,
        "filters": result.get("filters"),
        "role_featurization": result.get("role_featurization"),
        "rule_features": result.get("rule_features"),
        "starting_materials": formatted.get("starting_materials"),
    }



def _well_ids(n: int) -> List[str]:
    # Generate common plate IDs. For 24-well, 4 rows (A-D) x 6 columns (1-6).
    # General small helper that creates as close to square as possible grid.
    import math
    rows = int(math.sqrt(n))
    while rows > 1 and n % rows != 0:
        rows -= 1
    if rows <= 1:
        rows = min(8, n)  # cap rows to 8 for small plates
    cols = (n + rows - 1) // rows
    rows = max(1, min(8, rows))
    cols = max(1, cols)
    letters = [chr(ord('A') + i) for i in range(rows)]
    ids: List[str] = []
    for r in letters:
        for c in range(1, cols + 1):
            ids.append(f"{r}{c}")
            if len(ids) >= n:
                return ids
    return ids[:n]


def design_plate_from_reaction(
    reaction: str,
    plate_size: int = 24,
    relax: Dict[str, Any] | None = None,
    constraint_rules: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """Design a diversified plate across cores for a reaction.

    Returns dict with keys: csv (string), rows (list of dict), meta.
    """
    relax = dict(relax or {})

    # Build one precedent pack (larger k for variety)
    rec = recommend_from_reaction(reaction, k=max(plate_size, 50), relax=relax, constraint_rules=constraint_rules or {})
    pack = rec.get("precedent_pack") or {}
    precs: List[Dict[str, Any]] = list(pack.get("precedents") or [])

    # Core ranking by frequency
    core_counts = Counter([str(p.get("core") or "") for p in precs if p.get("core")])
    core_list = [c for c, _ in core_counts.most_common() if c]
    if not core_list:
        return {"csv": "well_id,core,base_uid,solvent_uid,additive_uids,T_C,time_h\n", "rows": [], "meta": {"error": "NO_PRECEDENTS"}}

    # Build per-core groups for base/solvent extraction
    by_core: Dict[str, List[Dict[str, Any]]] = {}
    for p in precs:
        c = str(p.get("core") or "")
        if not c:
            continue
        by_core.setdefault(c, []).append(p)

    # Ordered core sequence to fill plate, cycling if needed
    seq: List[str] = []
    while len(seq) < int(plate_size):
        for c in core_list:
            seq.append(c)
            if len(seq) >= int(plate_size):
                break

    # For each core, pick base/solvent and T/time from its group
    rows_out: List[Dict[str, Any]] = []
    for i, core in enumerate(seq):
        group = by_core.get(core, [])
        bases = [str(p.get("base_uid") or "") for p in group if p.get("base_uid")]
        solvents = [str(p.get("solvent_uid") or "") for p in group if p.get("solvent_uid")]
        base_counts = Counter(bases)
        solv_counts = Counter(solvents)
        base_list = [b for b, _ in base_counts.most_common()]
        solv_list = [s for s, _ in solv_counts.most_common()]
        b_pick, _bf = _pick_with_constraints(base_list, constraint_rules or {})
        s_pick, _sf = _pick_with_constraints(solv_list, constraint_rules or {})
        # Fallback across all precedents if needed
        if not b_pick:
            all_b = [str(p.get("base_uid") or "") for p in precs if p.get("base_uid")]
            b_pick, _ = _pick_with_constraints(list(dict.fromkeys(all_b)), constraint_rules or {})
        if not s_pick:
            all_s = [str(p.get("solvent_uid") or "") for p in precs if p.get("solvent_uid")]
            s_pick, _ = _pick_with_constraints(list(dict.fromkeys(all_s)), constraint_rules or {})
        # T/time medians per core
        def nums(key: str, items: List[Dict[str, Any]]):
            return [p.get(key) for p in items if isinstance(p.get(key), (int, float))]
        T_med = _median(nums("T_C", group) or nums("T_C", precs))
        t_med = _median(nums("time_h", group) or nums("time_h", precs))

        rows_out.append({
            "well_id": _well_ids(int(plate_size))[i],
            "core": core,
            "base_uid": b_pick or "",
            "solvent_uid": s_pick or "",
            "additive_uids": "",
            "T_C": T_med if T_med is not None else "",
            "time_h": t_med if t_med is not None else "",
        })

    # Render CSV
    header = ["well_id", "core", "base_uid", "solvent_uid", "additive_uids", "T_C", "time_h"]
    def _csv_escape(v: Any) -> str:
        s = "" if v is None else str(v)
        return ('"' + s.replace('"', '""') + '"') if ("," in s or '"' in s) else s
    lines = [",".join(header)]
    for row in rows_out:
        lines.append(
            ",".join(_csv_escape(row[k]) for k in header)
        )
    csv_text = "\n".join(lines) + "\n"

    meta = {
        "family": rec.get("family"),
        "bin": rec.get("bin"),
        "cores": core_counts.most_common(8),
        "precedent_support": int(pack.get("support") or len(precs)),
    }
    return {"csv": csv_text, "rows": rows_out, "meta": meta}
