import csv
import json
import html
import math
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List
from PyQt6 import QtCore, QtGui, QtWidgets


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

RUN_ALL_SOURCE_SENTINEL = "__run_all_recommendation__"
RUN_ALL_GROUPS: Tuple[str, ...] = ("literature", "motif", "rules")
_RECOMMENDER_CACHE: Dict[str, Any] = {}
_RECOMMEND_DM_CACHE: Dict[str, Any] = {}
PUBLIC_STRATEGIES: Tuple[str, ...] = ("motif", "rules", "literature", "similarity")


def _parse_reaction_smiles(reaction_smiles: str) -> Tuple[str, Optional[str], Optional[str]]:
    text = (reaction_smiles or "").strip()
    if not text:
        return "", None, None
    if ">>" in text:
        reactants_part, product = text.split(">>", 1)
        reactants = [r for r in reactants_part.split(".") if r]
        reactant_a = reactants[0] if reactants else ""
        reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
        product_smiles = product.strip() or None
        return reactant_a, reactant_b, product_smiles
    reactants = [r for r in text.split(".") if r]
    reactant_a = reactants[0] if reactants else ""
    reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
    return reactant_a, reactant_b, None


def _validate_reaction_smiles_input(reaction_smiles: str) -> Tuple[bool, str]:
    text = str(reaction_smiles or "").strip()
    if not text:
        return False, "Reaction SMILES is empty."

    try:
        from chemtools.smiles import normalize_reaction

        payload = normalize_reaction(text)
    except Exception as exc:
        return False, f"Failed to parse reaction SMILES: {exc}"

    reactants = list(payload.get("reactants") or [])
    if not reactants:
        return False, "Reaction SMILES must contain at least one reactant."

    invalid_inputs: List[str] = []
    for component in reactants + list(payload.get("agents") or []) + list(payload.get("products") or []):
        if not isinstance(component, dict):
            continue
        if str(component.get("error") or "").strip():
            bad = str(component.get("input") or "").strip()
            if bad:
                invalid_inputs.append(bad)

    if invalid_inputs:
        shown = ", ".join(invalid_inputs[:3])
        if len(invalid_inputs) > 3:
            shown += f" (+{len(invalid_inputs) - 3} more)"
        return False, f"Invalid SMILES token(s): {shown}"

    return True, ""


def _detect_csv_type(path: Path) -> str:
    parts = [part.lower() for part in path.parts]
    for label in ("protocols", "rules", "literature", "datasets", "motif", "experiments", "experiment", "experiements"):
        if label in parts:
            if label in ("literature", "datasets"):
                return "literature"
            if label in ("motif", "experiments", "experiment", "experiements"):
                return "motif"
            if label == "protocols":
                return "literature"
            return label
    if path.is_dir():
        return "directory"
    try:
        with path.open("r", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            header = next(reader, [])
    except Exception:
        return "unknown"
    header_lower = [h.strip().lower() for h in header]
    if "reaction_type_standardized" in header_lower:
        return "literature"
    if "reaction_id" in header_lower:
        return "literature"
    if "temperature_c" in header_lower:
        return "rules"
    if "reactant_1" in header_lower and "reactant_2" in header_lower:
        return "motif"
    return "unknown"


def _normalize_source_group_label(value: Any) -> str:
    text = str(value or "").strip().lower()
    if not text or text == "nan":
        return "unknown"
    if text in ("literature", "datasets", "dataset", "lit"):
        return "literature"
    if text in ("motif", "motifs", "experiments", "experiment", "experiements"):
        return "motif"
    if text == "rules":
        return "rules"
    if text in ("protocols", "protocol"):
        return "literature"
    if text == "precedent":
        return "precedent"
    return text


def _normalize_strategy_label(value: Any) -> str:
    text = str(value or "").strip().lower()
    if not text:
        return "motif"
    aliases = {
        "experiment": "motif",
        "experiments": "motif",
        "experimental": "motif",
        "experiment-based": "motif",
        "motif-based": "motif",
        "rule-based": "rules",
        "lit": "literature",
        "precedent": "similarity",
        "similarity-based": "similarity",
        "fingerprint": "similarity",
        "fingerprint-based": "similarity",
    }
    return aliases.get(text, text if text in PUBLIC_STRATEGIES else "motif")


def _recommendation_fingerprint(rec: Any) -> str:
    try:
        return json.dumps(asdict(rec), sort_keys=True, default=str)
    except Exception:
        return str(id(rec))


def _dedupe_recommendations(recs: List[Any]) -> List[Any]:
    out: List[Any] = []
    seen: set[str] = set()
    for rec in recs or []:
        key = _recommendation_fingerprint(rec)
        if key in seen:
            continue
        seen.add(key)
        out.append(rec)
    return out


def _normalize_recommendations_by_source(source_map: Dict[str, Any]) -> Dict[str, List[Any]]:
    normalized_map: Dict[str, List[Any]] = {}
    for key, items in (source_map or {}).items():
        normalized_key = _normalize_source_group_label(key)
        normalized_map.setdefault(normalized_key, []).extend(list(items or []))
    return normalized_map


def _resolve_db_path_for_source(db_path: str, source_group: str) -> str:
    normalized = _normalize_source_group_label(source_group)
    if normalized != "motif":
        return db_path
    path = Path(db_path)
    if not path.is_dir():
        return db_path
    for subdir in ("motif", "experiments"):
        candidate = path / subdir / "HTE_canonical.csv"
        if candidate.exists():
            return str(candidate)
    return db_path


def _recommender_cache_key(db_path: str) -> str:
    try:
        return str(Path(db_path).resolve())
    except Exception:
        return str(db_path)


def _get_cached_recommender(db_path: str) -> Any:
    from chemtools.recommend import HTERecommender

    key = _recommender_cache_key(db_path)
    recommender = _RECOMMENDER_CACHE.get(key)
    if recommender is None:
        recommender = HTERecommender(db_path)
        _RECOMMENDER_CACHE[key] = recommender
    return recommender


def _get_cached_data_manager(db_path: str) -> Any:
    from chemtools.recommend import RecommendationDataManager

    key = _recommender_cache_key(db_path)
    manager = _RECOMMEND_DM_CACHE.get(key)
    if manager is None:
        manager = RecommendationDataManager(base_db_path=db_path)
        _RECOMMEND_DM_CACHE[key] = manager
    return manager


def _format_float(value: Optional[float]) -> str:
    if value is None:
        return ""
    try:
        return f"{float(value):.2f}".rstrip("0").rstrip(".")
    except (TypeError, ValueError):
        return str(value)


def _safe_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def _format_items(value: Any, *, sep: str = ", ") -> str:
    if value is None:
        return "None"
    if isinstance(value, (list, tuple, set)):
        items = [str(v).strip() for v in value if str(v).strip()]
        return sep.join(items) if items else "None"
    text = str(value).strip()
    return text if text else "None"


def _extract_crk_section(crk_key: str, label: str) -> str:
    if not crk_key:
        return ""
    prefix = f"{label}: "
    parts = [p.strip() for p in crk_key.split(" | ") if p.strip()]
    for part in parts:
        if part.startswith(prefix):
            return part[len(prefix):].strip()
    return ""


def _collect_reaction_analysis(reaction_smiles: str) -> Dict[str, Any]:
    info: Dict[str, Any] = {
        "reaction_key": "",
        "product_broad_tags": [],
        "product_motifs_reactive": [],
        "reacted_motifs": [],
        "formed_motifs": [],
        "formed_motifs_center": [],
        "formed_motifs_context": [],
        "spectator_motifs": [],
        "spectator_groups": [],
        "bond_formed": "",
        "bond_broken": "",
    }
    if not reaction_smiles:
        return info
    try:
        from chemtools.featurizers.unified import featurize_reaction
    except Exception:
        return info

    try:
        payload = featurize_reaction(
            reaction_smiles,
            options={"confirm_coupling_products": True},
        )
    except Exception:
        return info

    if not isinstance(payload, dict):
        return info
    reaction = payload.get("reaction") if isinstance(payload.get("reaction"), dict) else payload
    if not isinstance(reaction, dict):
        return info
    aggregates = reaction.get("aggregates") or {}
    reaction_key = str(reaction.get("reaction_key") or "").strip()
    groups = aggregates.get("spectator_groups_ranked") or aggregates.get("spectator_groups_combined") or []
    cleaned = [str(group).strip() for group in groups if str(group).strip()]
    info.update(
        {
            "reaction_key": reaction_key,
            "product_broad_tags": list(reaction.get("product_broad_tags") or []),
            "product_motifs_reactive": list(reaction.get("product_motifs_reactive") or []),
            "reacted_motifs": list(aggregates.get("reacted_motifs") or []),
            "formed_motifs": list(aggregates.get("formed_motifs") or []),
            "formed_motifs_center": list(aggregates.get("formed_motifs_center") or []),
            "formed_motifs_context": list(aggregates.get("formed_motifs_context") or []),
            "spectator_motifs": list(aggregates.get("spectator_motifs") or []),
            "spectator_groups": cleaned,
            "bond_formed": _extract_crk_section(reaction_key, "bond_formed"),
            "bond_broken": _extract_crk_section(reaction_key, "bond_broken"),
        }
    )
    if os.environ.get("HTE_DEBUG_SPECTATORS") == "1":
        print(f"[HTE_DEBUG] reaction_smiles={reaction_smiles}")
        print(f"[HTE_DEBUG] spectator_groups_combined={cleaned}")
    return info


def _collect_reaction_spectator_groups(reaction_smiles: str) -> List[str]:
    info = _collect_reaction_analysis(reaction_smiles)
    return list(info.get("spectator_groups") or [])


def _format_score_stats(scores: List[float], *, max_items: int = 6) -> str:
    if not scores:
        return "None"
    preview = ", ".join(f"{score:.3f}" for score in scores[:max_items])
    if len(scores) > max_items:
        preview += f", ... (+{len(scores) - max_items} more)"
    return preview


def _reaction_image_path(prefix: str) -> Path:
    output_dir = PROJECT_ROOT / "results" / "visualization"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / f"{prefix}_{time.time_ns()}.png"


_REACTION_PREVIEW_RENDER_SIZE: Tuple[int, int] = (1800, 600)
_MOLECULE_PREVIEW_RENDER_SIZE: Tuple[int, int] = (900, 600)


def _fit_pixmap_to_screen(pixmap: QtGui.QPixmap) -> QtGui.QPixmap:
    if pixmap.isNull():
        return pixmap
    screen = QtGui.QGuiApplication.primaryScreen()
    if screen is None:
        return pixmap
    available = screen.availableGeometry()
    max_width = max(1, int(available.width() * 0.88))
    max_height = max(1, int(available.height() * 0.72))
    if pixmap.width() <= max_width and pixmap.height() <= max_height:
        return pixmap
    return pixmap.scaled(
        max_width,
        max_height,
        QtCore.Qt.AspectRatioMode.KeepAspectRatio,
        QtCore.Qt.TransformationMode.SmoothTransformation,
    )


def _format_nearby_groups(groups: List[str]) -> str:
    return ", ".join(groups) if groups else "None"


def _clean_reaction_key_text(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "none":
        return ""
    compact = text.replace(" ", "")
    if compact == "[]->[]||[]":
        return ""
    return text


def _format_matched_key(value: Any) -> str:
    if isinstance(value, tuple):
        parts = [str(item).strip() for item in value if str(item).strip()]
        if not parts:
            return ""
        if len(parts) == 1:
            return parts[0]
        return " + ".join(parts)
    text = str(value or "").strip()
    return text


def _table_columns_for_type(data_type: str) -> List[Tuple[str, str]]:
    if data_type in {"precedent", "similarity"}:
        return [
            ("Rank", "rank"),
            ("Similarity", "match_score"),
            ("Yield", "avg_yield"),
            ("Catalyst", "catalyst"),
            ("Ligand", "ligand"),
            ("Base", "base"),
            ("Solvent", "solvent"),
            ("Additive", "additive"),
            ("Reaction Type", "reaction_type"),
            ("Reaction ID", "reaction_id"),
            ("Reaction Key", "reaction_key"),
        ]
    base = [
        ("Rank", "rank"),
        ("Avg Z-Score", "avg_z_score"),
        ("Confidence", "confidence_score"),
        ("Success %", "success_rate"),
        ("Avg Yield", "avg_yield"),
        ("Experiments", "num_experiments"),
        ("Catalyst", "catalyst"),
        ("Ligand", "ligand"),
        ("Base", "base"),
        ("Solvent", "solvent"),
        ("Additive", "additive"),
        ("Condensation Agent", "coupling_reagent"),
        ("Spectator Groups", "spectator_groups"),
    ]
    columns = base + [
        ("Reaction Type", "reaction_type"),
        ("Reaction ID", "reaction_id"),
        ("Reaction Key", "reaction_key"),
        ("Reactant Types", "reactant_types"),
        ("Match Score", "match_score"),
    ]
    return columns


def _tab_label_for_source_group(source_group: str, strategy_label: str) -> str:
    normalized_source = _normalize_source_group_label(source_group)
    if normalized_source == "precedent":
        normalized_source = "similarity"
    return (normalized_source or "unknown").replace("_", " ").title()


def _reaction_info_from_stats(stats: Dict[str, Any]) -> Dict[str, Any]:
    analysis = stats.get("analysis") if isinstance(stats, dict) else None
    if not isinstance(analysis, dict):
        return {}
    return {
        "reaction_key": str(analysis.get("reaction_key") or "").strip(),
        "reacted_motifs": list(analysis.get("reacted_motifs") or []),
        "formed_motifs": list(analysis.get("formed_motifs") or []),
        "spectator_motifs": list(analysis.get("spectator_motifs") or []),
        "spectator_groups": list(analysis.get("spectator_groups") or []),
        "product_broad_tags": [],
        "product_motifs_reactive": [],
        "formed_motifs_center": [],
        "formed_motifs_context": [],
        "bond_formed": "",
        "bond_broken": "",
    }


class RecommendationWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal(bool, object, str, dict)
    progress = QtCore.pyqtSignal(str)

    def __init__(
        self,
        db_path: str,
        reaction_smiles: str,
        top_k: int,
        min_exp: int,
        reaction_filter: str,
        catalyst_filter: str,
        strategy: str,
        source_override: str,
        use_aryl_weighting: bool,
        prefer_mixfp_similarity: bool,
    ) -> None:
        super().__init__()
        self.db_path = db_path
        self.reaction_smiles = reaction_smiles
        self.top_k = top_k
        self.min_exp = min_exp
        self.reaction_filter = reaction_filter
        self.catalyst_filter = catalyst_filter
        self.strategy = _normalize_strategy_label(strategy)
        self.source_override = _normalize_source_group_label(source_override) if source_override else ""
        self.use_aryl_weighting = use_aryl_weighting
        self.prefer_mixfp_similarity = bool(prefer_mixfp_similarity)

    def _run_single(
        self,
        recommender: Any,
        reactant_a: str,
        reactant_b: Optional[str],
        product: Optional[str],
        source_group: Optional[str],
    ) -> Any:
        return recommender.recommend(
            reactant_a_smiles=reactant_a,
            reactant_b_smiles=reactant_b,
            product_smiles=product,
            top_k=self.top_k,
            min_experiments=self.min_exp,
            reaction_type_filter=self.reaction_filter or None,
            catalyst_filter=self.catalyst_filter or None,
            source_group=source_group,
            use_aryl_steric_electronic_weighting=self.use_aryl_weighting,
        )

    def _run_all_recommendations(
        self,
        recommender: Any,
        reactant_a: str,
        reactant_b: Optional[str],
        product: Optional[str],
    ) -> Any:
        self.progress.emit("Running all recommendation sources ...")
        baseline = self._run_single(
            recommender,
            reactant_a,
            reactant_b,
            product,
            None,
        )

        baseline_map = _normalize_recommendations_by_source(
            getattr(baseline, "recommendations_by_source", {}) or {}
        )
        merged_by_source: Dict[str, List[Any]] = {}

        for group_name in RUN_ALL_GROUPS:
            self.progress.emit(f"Running {group_name} recommendation ...")
            group_result = self._run_single(
                recommender,
                reactant_a,
                reactant_b,
                product,
                group_name,
            )
            group_map = _normalize_recommendations_by_source(
                getattr(group_result, "recommendations_by_source", {}) or {}
            )
            group_recs = list(group_map.get(group_name) or [])
            if not group_recs:
                group_recs = list(getattr(group_result, "recommendations", []) or [])
            merged_by_source[group_name] = _dedupe_recommendations(group_recs)

            if group_name == "literature":
                precedent_recs = list(group_map.get("precedent") or [])
                if precedent_recs:
                    merged_by_source["precedent"] = _dedupe_recommendations(precedent_recs)

        for key, items in baseline_map.items():
            if key not in merged_by_source and items:
                merged_by_source[key] = _dedupe_recommendations(list(items))

        source_order = ["literature", "rules", "motif", "precedent"]
        baseline.recommendations_by_source = merged_by_source
        return baseline

    def _run_similarity_only(
        self,
        recommender: Any,
        reactant_a: str,
        reactant_b: Optional[str],
        product: Optional[str],
    ) -> Any:
        self.progress.emit("Running similarity recommendation ...")
        result = self._run_single(
            recommender,
            reactant_a,
            reactant_b,
            product,
            "literature",
        )
        source_map = _normalize_recommendations_by_source(
            getattr(result, "recommendations_by_source", {}) or {}
        )
        precedent_recs = _dedupe_recommendations(list(source_map.get("precedent") or []))
        result.recommendations_by_source = {"similarity": precedent_recs}
        result.recommendations = precedent_recs
        return result

    def _run_motif_only(
        self,
        recommender: Any,
        reactant_a: str,
        reactant_b: Optional[str],
        product: Optional[str],
    ) -> Any:
        narrowed = self.source_override if self.source_override in {"motif", "rules"} else ""
        if narrowed:
            self.progress.emit(f"Running motif recommendation ({narrowed}) ...")
            result = self._run_single(
                recommender,
                reactant_a,
                reactant_b,
                product,
                narrowed,
            )
            source_map = _normalize_recommendations_by_source(getattr(result, "recommendations_by_source", {}) or {})
            group_recs = _dedupe_recommendations(list(source_map.get(narrowed) or list(getattr(result, "recommendations", []) or [])))
            result.recommendations_by_source = {narrowed: group_recs}
            result.recommendations = group_recs
            return result

        self.progress.emit("Running motif recommendation (motif + rules) ...")
        baseline = self._run_single(recommender, reactant_a, reactant_b, product, None)
        merged_by_source: Dict[str, List[Any]] = {}
        for group_name in ("motif", "rules"):
            group_result = self._run_single(
                recommender,
                reactant_a,
                reactant_b,
                product,
                group_name,
            )
            group_map = _normalize_recommendations_by_source(
                getattr(group_result, "recommendations_by_source", {}) or {}
            )
            group_recs = list(group_map.get(group_name) or [])
            if not group_recs:
                group_recs = list(getattr(group_result, "recommendations", []) or [])
            merged_by_source[group_name] = _dedupe_recommendations(group_recs)

        combined: List[Any] = []
        indices = {"motif": 0, "rules": 0}
        while self.top_k <= 0 or len(combined) < self.top_k:
            progressed = False
            for group_name in ("motif", "rules"):
                items = merged_by_source.get(group_name) or []
                idx = indices[group_name]
                if idx >= len(items):
                    continue
                combined.append(items[idx])
                indices[group_name] = idx + 1
                progressed = True
                if self.top_k > 0 and len(combined) >= self.top_k:
                    break
            if not progressed:
                break

        baseline.recommendations_by_source = merged_by_source
        baseline.recommendations = combined if combined else _dedupe_recommendations(list(getattr(baseline, "recommendations", []) or []))
        return baseline

    def _run_literature_only(
        self,
        recommender: Any,
        reactant_a: str,
        reactant_b: Optional[str],
        product: Optional[str],
    ) -> Any:
        self.progress.emit("Running literature recommendation ...")
        result = self._run_single(
            recommender,
            reactant_a,
            reactant_b,
            product,
            "literature",
        )
        source_map = _normalize_recommendations_by_source(
            getattr(result, "recommendations_by_source", {}) or {}
        )
        literature_recs = _dedupe_recommendations(list(source_map.get("literature") or list(getattr(result, "recommendations", []) or [])))
        result.recommendations_by_source = {"literature": literature_recs}
        result.recommendations = literature_recs
        return result

    def run(self) -> None:
        try:
            from chemtools.recommend import RecommendationRequest
            from chemtools.recommend.api import recommend as recommend_facade

            normalized_strategy = _normalize_strategy_label(self.strategy)
            source_override = self.source_override if self.source_override in {"literature", "motif", "rules"} else ""
            source_group = source_override or "any"

            self.progress.emit(
                f"Running {normalized_strategy} recommendation (datasets auto-load and reuse cache) ..."
            )
            dm = _get_cached_data_manager(self.db_path)
            request = RecommendationRequest(
                reaction_smiles=self.reaction_smiles,
                strategy=normalized_strategy,
                source_group=source_group,
                top_k=self.top_k,
                min_experiments=self.min_exp,
                reaction_type_filter=self.reaction_filter or None,
                catalyst_filter=self.catalyst_filter or None,
                hte_db_path=self.db_path,
                use_aryl_steric_electronic_weighting=self.use_aryl_weighting,
                prefer_mixfp_for_similarity=self.prefer_mixfp_similarity,
            )
            run_result = recommend_facade(request, data_manager=dm)
            result = getattr(run_result, "recommendation", None)
            if result is None:
                raise RuntimeError("No recommendation result returned by facade.")
            stats = {
                "loaded_resources": getattr(run_result, "loaded_resources", {}) or {},
                "plan": asdict(getattr(run_result, "plan")) if getattr(run_result, "plan", None) else {},
                "analysis": asdict(getattr(run_result, "analysis")) if getattr(run_result, "analysis", None) else {},
            }
            self.finished.emit(True, result, "OK", stats)
        except Exception as exc:
            self.finished.emit(False, None, str(exc), {})


class CacheWarmWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal(bool, object, str)
    progress = QtCore.pyqtSignal(str)

    def __init__(self, db_path: str, source_group: str) -> None:
        super().__init__()
        self.db_path = db_path
        self.source_group = source_group

    def run(self) -> None:
        try:
            from chemtools.recommend.api import warm_recommendation_cache

            self.progress.emit("Prebuilding HTE cache ...")
            dm = _get_cached_data_manager(self.db_path)
            summary = warm_recommendation_cache(
                source_group=self.source_group or "all",
                data_manager=dm,
                hte_db_path=self.db_path,
            )
            self.finished.emit(True, summary, "OK")
        except Exception as exc:
            self.finished.emit(False, None, str(exc))


class KMNRebuildWorker(QtCore.QObject):
    """Run ``scripts/A_build_kmn_index.py`` in a background thread.

    Uses ``--index-only`` by default (fast: re-serialises the FAISS index
    from existing weights without re-training).  Pass ``train=True`` to
    trigger a full re-train (takes much longer).
    """

    finished = QtCore.pyqtSignal(bool, str)  # (success, message)
    progress = QtCore.pyqtSignal(str)

    def __init__(self, *, train: bool = False) -> None:
        super().__init__()
        self._train = train

    def run(self) -> None:
        script = str(PROJECT_ROOT / "scripts" / "A_build_kmn_index.py")
        cmd = [sys.executable, script]
        if not self._train:
            cmd.append("--index-only")
        self.progress.emit("Running KMN index build ...")
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(PROJECT_ROOT),
            )
            if result.returncode == 0:
                self.finished.emit(True, (result.stdout or "").strip() or "KMN index rebuilt successfully.")
            else:
                err = (result.stderr or result.stdout or "Unknown error").strip()
                self.finished.emit(False, f"KMN rebuild failed (exit {result.returncode}):\n{err}")
        except Exception as exc:
            self.finished.emit(False, f"KMN rebuild error: {exc}")


class HTERecommenderWindow(QtWidgets.QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("HTE Recommender (Qt6)")
        self.resize(1300, 850)

        self.db_path_edit = QtWidgets.QLineEdit()
        self.db_path_edit.setPlaceholderText("Select a CSV/JSONL or a folder like data/HTE_db")
        self.data_type_label = QtWidgets.QLabel("Data type: unknown")
        self.data_type_label.setStyleSheet("color: #555;")

        self.reaction_smiles_edit = QtWidgets.QLineEdit()
        self.reaction_smiles_edit.setPlaceholderText("A.B>>P or A.B or A")

        self.top_k_spin = QtWidgets.QSpinBox()
        self.top_k_spin.setRange(0, 200)
        self.top_k_spin.setSpecialValueText("All")
        self.top_k_spin.setValue(10)
        self.top_k_spin.setToolTip("Max recommendations to display. Set to All (0) to show every matched condition.")

        self.min_exp_spin = QtWidgets.QSpinBox()
        self.min_exp_spin.setRange(1, 200)
        self.min_exp_spin.setValue(1)

        self.reaction_filter_edit = QtWidgets.QLineEdit()
        self.reaction_filter_edit.setPlaceholderText("Optional reaction type/category filter (e.g., Suzuki)")

        self.catalyst_filter_edit = QtWidgets.QLineEdit()
        self.catalyst_filter_edit.setPlaceholderText("Optional catalyst filter (e.g., Pd, Cu)")

        self.strategy_combo = QtWidgets.QComboBox()
        self.strategy_combo.addItems(
            [
                "motif",
                "rules",
                "literature",
                "similarity",
            ]
        )
        self.strategy_combo.setCurrentText("motif")
        self.strategy_combo.setToolTip(
            "Public recommendation strategy. motif = motif-source + rules by default. "
            "similarity uses literature precedents ranked by similarity."
        )

        self.source_group_combo = QtWidgets.QComboBox()
        self.source_group_combo.addItems(
            [
                "Auto",
                "literature",
                "motif",
                "rules",
            ]
        )
        self.source_group_combo.setCurrentText("Auto")
        self.source_group_combo.setToolTip(
            "Optional source override. Leave as Auto for strategy-driven source loading."
        )

        self.aryl_weighting_check = QtWidgets.QCheckBox("Aryl steric/electronic weighting")
        self.aryl_weighting_check.setToolTip("Reweight matches by aryl steric/electronic similarity (when available).")
        self.prefer_mixfp_check = QtWidgets.QCheckBox("Prefer MixFP for similarity")
        self.prefer_mixfp_check.setChecked(True)
        self.prefer_mixfp_check.setToolTip(
            "When similarity/literature precedents are used, prefer MixFP routing/blending when available; "
            "otherwise rely on DRFP/categorical similarity."
        )

        self.run_button = QtWidgets.QPushButton("Run Recommendation")
        self.prebuild_cache_button = QtWidgets.QPushButton("Prebuild HTE Cache")
        self.rebuild_kmn_button = QtWidgets.QPushButton("Rebuild KMN Index")
        self.clear_button = QtWidgets.QPushButton("Clear Results")
        self.export_json_button = QtWidgets.QPushButton("Export JSON")
        self.run_button.setMinimumSize(160, 34)
        self.prebuild_cache_button.setMinimumSize(150, 34)
        self.rebuild_kmn_button.setMinimumSize(150, 34)
        self.clear_button.setMinimumSize(140, 34)
        self.export_json_button.setMinimumSize(140, 34)
        self.prebuild_cache_button.setToolTip(
            "Prebuild the HTE reactant/transformation lookup index.\n"
            "Saved to: results/hte_cache/ (not data/kmn_index/).\n"
            "If the index already exists and data has not changed, this is instant."
        )
        self.rebuild_kmn_button.setToolTip(
            "Rebuild the KMN/FAISS similarity index (data/kmn_index/).\n"
            "Runs: scripts/A_build_kmn_index.py --index-only\n"
            "Required for similarity-based recommendations after data changes."
        )
        self.export_json_button.setEnabled(False)
        self.run_button.setStyleSheet(
            "QPushButton {"
            " background-color: #2b6cb0;"
            " color: white;"
            " font-weight: 600;"
            " border-radius: 4px;"
            " padding: 4px 10px;"
            "}"
            "QPushButton:hover { background-color: #2c5282; }"
            "QPushButton:pressed { background-color: #2a4365; }"
            "QPushButton:disabled { background-color: #9bb7d6; color: #f2f2f2; }"
        )

        self.summary = QtWidgets.QTextEdit()
        self.summary.setReadOnly(True)

        self.results_tabs = QtWidgets.QTabWidget()
        self._initialize_result_tabs()

        self.status = QtWidgets.QLabel("")
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setMinimumWidth(160)

        self._setup_layout()
        self._bind_signals()
        self._set_default_path()

        self.thread: Optional[QtCore.QThread] = None
        self.worker: Optional[RecommendationWorker] = None
        self.cache_thread: Optional[QtCore.QThread] = None
        self.cache_worker: Optional[CacheWarmWorker] = None
        self.kmn_thread: Optional[QtCore.QThread] = None
        self.kmn_worker: Optional["KMNRebuildWorker"] = None
        self._reaction_dialog: Optional[QtWidgets.QDialog] = None
        self._spectator_groups_summary: str = ""
        self._strategy_label: str = "motif"
        self._source_group_label: str = "Auto"
        self._aryl_weighting_enabled: bool = False
        self._prefer_mixfp_similarity_enabled: bool = True
        self._all_json_output: Optional[Dict[str, Any]] = None
        self._last_result_obj: Optional[object] = None
        self._last_export_context: Dict[str, Any] = {}

    def _setup_layout(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        title = QtWidgets.QLabel("HTE Recommender")
        title.setStyleSheet("font-size: 18px; font-weight: bold;")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        db_row = QtWidgets.QHBoxLayout()
        db_row.addWidget(self.db_path_edit)
        browse_file_btn = QtWidgets.QPushButton("Select CSV")
        browse_dir_btn = QtWidgets.QPushButton("Select Folder")
        db_row.addWidget(browse_file_btn)
        db_row.addWidget(browse_dir_btn)

        browse_file_btn.clicked.connect(self._choose_file)
        browse_dir_btn.clicked.connect(self._choose_dir)

        layout.addLayout(db_row)
        layout.addWidget(self.data_type_label)

        form = QtWidgets.QFormLayout()
        form.addRow("Reaction SMILES:", self.reaction_smiles_edit)

        options_row = QtWidgets.QHBoxLayout()
        options_row.addWidget(QtWidgets.QLabel("Top K:"))
        options_row.addWidget(self.top_k_spin)
        options_row.addSpacing(20)
        options_row.addWidget(QtWidgets.QLabel("Min Experiments:"))
        options_row.addWidget(self.min_exp_spin)
        options_row.addStretch()
        form.addRow("Options:", options_row)

        filters_row = QtWidgets.QHBoxLayout()
        filters_row.addWidget(QtWidgets.QLabel("Reaction:"))
        filters_row.addWidget(self.reaction_filter_edit)
        filters_row.addSpacing(12)
        filters_row.addWidget(QtWidgets.QLabel("Catalyst:"))
        filters_row.addWidget(self.catalyst_filter_edit)
        filters_row.addSpacing(12)
        filters_row.addWidget(QtWidgets.QLabel("Strategy:"))
        filters_row.addWidget(self.strategy_combo)
        filters_row.addSpacing(12)
        filters_row.addWidget(QtWidgets.QLabel("Source Override:"))
        filters_row.addWidget(self.source_group_combo)
        filters_row.addStretch()
        form.addRow("Filters:", filters_row)
        weighting_row = QtWidgets.QHBoxLayout()
        weighting_row.addWidget(self.aryl_weighting_check)
        weighting_row.addSpacing(12)
        weighting_row.addWidget(self.prefer_mixfp_check)
        weighting_row.addStretch()
        form.addRow("Weighting:", weighting_row)
        layout.addLayout(form)

        button_row = QtWidgets.QHBoxLayout()
        button_row.addWidget(self.run_button)
        button_row.addWidget(self.prebuild_cache_button)
        button_row.addWidget(self.rebuild_kmn_button)
        button_row.addWidget(self.clear_button)
        button_row.addWidget(self.export_json_button)
        button_row.addStretch()
        layout.addLayout(button_row)

        layout.addWidget(QtWidgets.QLabel("Summary:"))
        layout.addWidget(self.summary, stretch=1)
        layout.addWidget(QtWidgets.QLabel("Recommendations:"))
        layout.addWidget(self.results_tabs, stretch=3)
        status_row = QtWidgets.QHBoxLayout()
        status_row.addWidget(self.progress_bar)
        status_row.addWidget(self.status)
        status_row.addStretch()
        layout.addLayout(status_row)

    def _bind_signals(self) -> None:
        self.run_button.clicked.connect(self._run_recommendation)
        self.prebuild_cache_button.clicked.connect(self._run_prebuild_cache)
        self.rebuild_kmn_button.clicked.connect(self._run_rebuild_kmn_index)
        self.clear_button.clicked.connect(self._clear_results)
        self.export_json_button.clicked.connect(self._export_json_output)
        self.db_path_edit.textChanged.connect(self._update_data_type)

    def _set_default_path(self) -> None:
        default_path = PROJECT_ROOT / "data" / "HTE_db"
        if default_path.exists():
            self.db_path_edit.setText(str(default_path))

    def _choose_file(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select HTE CSV/JSONL",
            str(PROJECT_ROOT),
            "HTE Files (*.csv *.jsonl);;All Files (*)",
        )
        if path:
            self.db_path_edit.setText(path)

    def _choose_dir(self) -> None:
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select HTE Folder",
            str(PROJECT_ROOT),
        )
        if path:
            self.db_path_edit.setText(path)

    def _update_data_type(self) -> None:
        path_text = self.db_path_edit.text().strip()
        if not path_text:
            self.data_type_label.setText("Data type: unknown")
            return
        path = Path(path_text)
        if not path.exists():
            self.data_type_label.setText("Data type: missing path")
            return
        data_type = _detect_csv_type(path)
        self.data_type_label.setText(f"Data type: {data_type}")

    def _create_results_table(self) -> QtWidgets.QTableWidget:
        table = QtWidgets.QTableWidget(0, 0)
        table.setSortingEnabled(True)
        table.horizontalHeader().setStretchLastSection(True)
        table.setAlternatingRowColors(True)
        return table

    def _initialize_result_tabs(self) -> None:
        self.results_tabs.clear()
        strategy_label = self.strategy_combo.currentText().strip() if hasattr(self, "strategy_combo") else "motif"
        for source_group in ("literature", "rules", "motif", "similarity"):
            group_table = self._create_results_table()
            self.results_tabs.addTab(group_table, _tab_label_for_source_group(source_group, strategy_label))

    def _clear_results(self) -> None:
        self.summary.clear()
        self._initialize_result_tabs()
        self.status.setText("")
        self._spectator_groups_summary = ""
        self._all_json_output = None
        self._last_result_obj = None
        self._last_export_context = {}
        self.export_json_button.setEnabled(False)
        self.run_button.setEnabled(True)
        self.prebuild_cache_button.setEnabled(True)
        self.rebuild_kmn_button.setEnabled(True)
        self.run_button.setText("Run Recommendation")
        self.prebuild_cache_button.setText("Prebuild HTE Cache")
        self.rebuild_kmn_button.setText("Rebuild KMN Index")
        self.progress_bar.setVisible(False)
        if self._reaction_dialog:
            self._reaction_dialog.close()
            self._reaction_dialog = None

    def _run_recommendation(self) -> None:
        if self.cache_thread and self.cache_thread.isRunning():
            QtWidgets.QMessageBox.information(
                self,
                "HTE Recommender",
                "Cache prebuild is running. Please wait for it to finish.",
            )
            return
        db_path = self.db_path_edit.text().strip()
        reaction_smiles = self.reaction_smiles_edit.text().strip()
        if not db_path:
            QtWidgets.QMessageBox.warning(self, "Missing data", "Select a CSV/JSONL or folder.")
            return
        if not reaction_smiles:
            QtWidgets.QMessageBox.warning(self, "Missing input", "Provide a reaction SMILES.")
            return
        if not Path(db_path).exists():
            QtWidgets.QMessageBox.warning(self, "Invalid path", "The data path does not exist.")
            return
        is_valid, validation_error = _validate_reaction_smiles_input(reaction_smiles)
        if not is_valid:
            QtWidgets.QMessageBox.warning(self, "Invalid reaction SMILES", validation_error)
            return

        self._show_reaction_image(reaction_smiles)
        self.summary.clear()
        self.status.setText("Working...")
        self.run_button.setEnabled(False)
        self.prebuild_cache_button.setEnabled(False)
        self.rebuild_kmn_button.setEnabled(False)
        self.run_button.setText("Running...")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self._strategy_label = self.strategy_combo.currentText().strip()
        self._source_group_label = self.source_group_combo.currentText().strip()
        self._aryl_weighting_enabled = bool(self.aryl_weighting_check.isChecked())
        self._prefer_mixfp_similarity_enabled = bool(self.prefer_mixfp_check.isChecked())

        self.thread = QtCore.QThread()
        self.worker = RecommendationWorker(
            db_path=db_path,
            reaction_smiles=reaction_smiles,
            top_k=self.top_k_spin.value(),
            min_exp=self.min_exp_spin.value(),
            reaction_filter=self.reaction_filter_edit.text().strip(),
            catalyst_filter=self.catalyst_filter_edit.text().strip(),
            strategy=self._strategy_label,
            source_override=("" if self._source_group_label == "Auto" else self._source_group_label),
            use_aryl_weighting=self._aryl_weighting_enabled,
            prefer_mixfp_similarity=self._prefer_mixfp_similarity_enabled,
        )
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self._append_status)
        self.worker.finished.connect(self._on_finished)
        self.thread.start()

    def _append_status(self, message: str) -> None:
        self.status.setText(message)

    def _run_prebuild_cache(self) -> None:
        if self.thread and self.thread.isRunning():
            QtWidgets.QMessageBox.information(
                self,
                "HTE Recommender",
                "Recommendation is running. Please wait for it to finish.",
            )
            return
        if self.cache_thread and self.cache_thread.isRunning():
            QtWidgets.QMessageBox.information(
                self,
                "HTE Recommender",
                "Cache prebuild is already running.",
            )
            return

        db_path = self.db_path_edit.text().strip()
        if not db_path:
            QtWidgets.QMessageBox.warning(self, "Missing data", "Select a CSV/JSONL or folder.")
            return
        if not Path(db_path).exists():
            QtWidgets.QMessageBox.warning(self, "Invalid path", "The data path does not exist.")
            return

        selected_strategy = _normalize_strategy_label(self.strategy_combo.currentText().strip())
        selected_override = self.source_group_combo.currentText().strip()
        normalized_override = _normalize_source_group_label(selected_override)
        if selected_override != "Auto" and normalized_override in {"literature", "motif", "rules"}:
            source_group = normalized_override
        elif selected_strategy in {"literature", "similarity"}:
            source_group = "literature"
        elif selected_strategy == "rules":
            source_group = "rules"
        else:
            source_group = "all"

        # ── Quick disk-cache validity check (no data load) ────────────────
        # If the index already exists and data hasn't changed, skip the thread
        # and just report the status without showing the full rebuild spinner.
        try:
            from chemtools.recommend.recommender import check_hte_cache_status
            cache_status = check_hte_cache_status(db_path, source_group=source_group)
        except Exception:
            cache_status = {"valid": False, "targets": []}

        if cache_status.get("valid"):
            targets = cache_status.get("targets") or []
            sources = [t.get("status", "unknown") for t in targets]
            source_summary = ", ".join(f"{s}={sources.count(s)}" for s in dict.fromkeys(sources) if sources.count(s) > 0)
            reply = QtWidgets.QMessageBox.question(
                self,
                "HTE Cache Already Up To Date",
                (
                    f"The HTE index cache is already valid ({source_summary}).\n"
                    "Note: this cache is stored in results/hte_cache/, NOT data/kmn_index/.\n\n"
                    "Force rebuild anyway?"
                ),
                QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
                QtWidgets.QMessageBox.StandardButton.No,
            )
            if reply != QtWidgets.QMessageBox.StandardButton.Yes:
                self.status.setText(f"HTE cache already up to date ({source_summary}).")
                return

        self.status.setText("Prebuilding HTE cache...")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self.run_button.setEnabled(False)
        self.prebuild_cache_button.setEnabled(False)
        self.rebuild_kmn_button.setEnabled(False)
        self.prebuild_cache_button.setText("Prebuilding...")

        self.cache_thread = QtCore.QThread()
        self.cache_worker = CacheWarmWorker(
            db_path=db_path,
            source_group=source_group,
        )
        self.cache_worker.moveToThread(self.cache_thread)
        self.cache_thread.started.connect(self.cache_worker.run)
        self.cache_worker.progress.connect(self._append_status)
        self.cache_worker.finished.connect(self._on_prebuild_finished)
        self.cache_thread.start()

    def _on_prebuild_finished(self, success: bool, summary: object, message: str) -> None:
        if self.cache_thread:
            self.cache_thread.quit()
            self.cache_thread.wait()
        self.run_button.setEnabled(True)
        self.prebuild_cache_button.setEnabled(True)
        self.rebuild_kmn_button.setEnabled(True)
        self.prebuild_cache_button.setText("Prebuild HTE Cache")
        self.progress_bar.setVisible(False)

        if not success:
            self.status.setText("Error")
            QtWidgets.QMessageBox.critical(self, "HTE Recommender", message)
            return

        targets = []
        if isinstance(summary, dict):
            targets = list(summary.get("targets") or [])
        target_count = len(targets)
        total_rows = sum(int(item.get("num_rows", 0) or 0) for item in targets if isinstance(item, dict))
        total_elapsed = sum(float(item.get("elapsed_s", 0.0) or 0.0) for item in targets if isinstance(item, dict))
        cache_sources = [
            str(item.get("cache_source") or "").strip().lower()
            for item in targets
            if isinstance(item, dict)
        ]
        source_counts = {
            label: sum(1 for src in cache_sources if src == label)
            for label in ("memory", "disk", "rebuilt")
        }
        source_parts = [f"{label}={count}" for label, count in source_counts.items() if count > 0]
        source_summary = ", ".join(source_parts) if source_parts else "source=unknown"
        if source_counts.get("rebuilt", 0) == 0 and target_count > 0:
            self.status.setText(f"HTE cache already ready ({source_summary}).")
        else:
            self.status.setText(f"HTE cache ready ({source_summary}).")

        # ── KMN index status ──────────────────────────────────────────────
        try:
            from chemtools.util.faiss_router import is_index_built, INDEX_DIR
            kmn_built = is_index_built()
            kmn_dir = INDEX_DIR
        except Exception:
            kmn_built = False
            kmn_dir = str(PROJECT_ROOT / "data" / "kmn_index")
        kmn_line = (
            f"KMN index (data/kmn_index/): {'\u2713 exists' if kmn_built else '\u2717 MISSING or empty'}\n"
            f"  To rebuild: click 'Rebuild KMN Index' or run\n"
            f"  python scripts/A_build_kmn_index.py --index-only"
        )

        QtWidgets.QMessageBox.information(
            self,
            "HTE Cache Prebuild Complete",
            (
                "HTE reactant index rebuilt/verified.\n"
                f"  Saved to: results/hte_cache/\n"
                f"  Targets:  {target_count}\n"
                f"  Rows:     {total_rows}\n"
                f"  Source:   {source_summary}\n"
                f"  Elapsed:  {total_elapsed:.2f}s\n\n"
                + kmn_line
            ),
        )

    # ── KMN index rebuild ─────────────────────────────────────────────────

    def _run_rebuild_kmn_index(self) -> None:
        if self.kmn_thread and self.kmn_thread.isRunning():
            QtWidgets.QMessageBox.information(
                self, "HTE Recommender", "KMN rebuild is already running."
            )
            return
        if self.thread and self.thread.isRunning() or self.cache_thread and self.cache_thread.isRunning():
            QtWidgets.QMessageBox.information(
                self,
                "HTE Recommender",
                "Another operation is running. Please wait for it to finish.",
            )
            return

        reply = QtWidgets.QMessageBox.question(
            self,
            "Rebuild KMN Index",
            (
                "This will run:\n"
                "  python scripts/A_build_kmn_index.py --index-only\n\n"
                "This re-serialises the FAISS index from existing KMN weights "
                "(fast, no re-training).\n"
                "The index will be saved to:  data/kmn_index/\n\n"
                "Proceed?"
            ),
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.Yes,
        )
        if reply != QtWidgets.QMessageBox.StandardButton.Yes:
            return

        self.status.setText("Rebuilding KMN index...")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self.run_button.setEnabled(False)
        self.prebuild_cache_button.setEnabled(False)
        self.rebuild_kmn_button.setEnabled(False)
        self.rebuild_kmn_button.setText("Rebuilding KMN...")

        self.kmn_thread = QtCore.QThread()
        self.kmn_worker = KMNRebuildWorker(train=False)
        self.kmn_worker.moveToThread(self.kmn_thread)
        self.kmn_thread.started.connect(self.kmn_worker.run)
        self.kmn_worker.progress.connect(self._append_status)
        self.kmn_worker.finished.connect(self._on_kmn_rebuild_finished)
        self.kmn_thread.start()

    def _on_kmn_rebuild_finished(self, success: bool, message: str) -> None:
        if self.kmn_thread:
            self.kmn_thread.quit()
            self.kmn_thread.wait()
        self.run_button.setEnabled(True)
        self.prebuild_cache_button.setEnabled(True)
        self.rebuild_kmn_button.setEnabled(True)
        self.rebuild_kmn_button.setText("Rebuild KMN Index")
        self.progress_bar.setVisible(False)

        if not success:
            self.status.setText("KMN rebuild failed")
            QtWidgets.QMessageBox.critical(self, "KMN Rebuild Failed", message)
            return

        try:
            from chemtools.util.faiss_router import is_index_built
            kmn_built = is_index_built()
        except Exception:
            kmn_built = False

        self.status.setText("KMN index rebuilt." if kmn_built else "KMN rebuild done (verify manually).")
        QtWidgets.QMessageBox.information(
            self,
            "KMN Index Rebuild Complete",
            (
                f"{'KMN index exists and is non-empty.' if kmn_built else 'KMN rebuild ran, but index may be empty — check logs.'}\n"
                f"Saved to: data/kmn_index/\n\n"
                + (message[:500] if message else "")
            ),
        )

    def _on_finished(self, success: bool, result: object, message: str, stats: Dict[str, Any]) -> None:
        if self.thread:
            self.thread.quit()
            self.thread.wait()
        self.run_button.setEnabled(True)
        self.prebuild_cache_button.setEnabled(True)
        self.rebuild_kmn_button.setEnabled(True)
        self.run_button.setText("Run Recommendation")
        self.progress_bar.setVisible(False)

        if not success:
            self.status.setText("Error")
            QtWidgets.QMessageBox.critical(self, "HTE Recommender", message)
            return

        self.status.setText("Done")
        self._render_result(result, stats)

    def _populate_table(
        self,
        table: QtWidgets.QTableWidget,
        recs: List[Any],
        columns: List[Tuple[str, str]],
    ) -> None:
        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels([label for label, _ in columns])
        table.setRowCount(len(recs))

        for row_index, rec in enumerate(recs):
            rec_dict = asdict(rec)
            rec_dict["rank"] = row_index + 1
            reactant_types = []
            for item in list(getattr(rec, "reactant_types", []) or []):
                text = _safe_text(item).strip()
                if text:
                    reactant_types.append(text)
            rec_dict["reactant_types"] = " + ".join(reactant_types)
            for col_index, (_, key) in enumerate(columns):
                value = rec_dict.get(key)
                if isinstance(value, float):
                    cell_text = _format_float(value)
                else:
                    cell_text = str(value) if value is not None else ""
                item = QtWidgets.QTableWidgetItem(cell_text)
                table.setItem(row_index, col_index, item)

        table.resizeColumnsToContents()

    def _render_result(self, result: object, stats: Dict[str, Any]) -> None:
        data_type = _detect_csv_type(Path(self.db_path_edit.text().strip()))
        reaction_smiles = self.reaction_smiles_edit.text().strip()

        recs = list(getattr(result, "recommendations", []) or [])
        type_a = getattr(result, "reactant_a_type", "")
        type_b = getattr(result, "reactant_b_type", "")
        detected = getattr(result, "matched_motifs", None)
        reaction_info = _reaction_info_from_stats(stats) or _collect_reaction_analysis(reaction_smiles)
        spectator_summary = _format_nearby_groups(list(reaction_info.get("spectator_groups") or []))
        self._spectator_groups_summary = spectator_summary

        query_key = _clean_reaction_key_text(
            getattr(result, "query_reaction_key", None) or reaction_info.get("reaction_key", "")
        )
        reacted_motifs = getattr(result, "reacted_motifs", None) or reaction_info.get("reacted_motifs", [])
        formed_motifs = getattr(result, "formed_motifs", None) or reaction_info.get("formed_motifs", [])
        spectator_motifs = getattr(result, "spectator_motifs", None) or reaction_info.get("spectator_motifs", [])
        matched_key = _clean_reaction_key_text(_format_matched_key(detected))
        reaction_filter = self.reaction_filter_edit.text().strip() or None
        detected_reaction_type = getattr(result, "predicted_reaction_type", None)
        detected_confidence = getattr(result, "reaction_type_confidence", None)
        query_spectator_groups = (
            getattr(result, "query_spectator_groups", None)
            or tuple(reaction_info.get("spectator_groups") or [])
        )
        scoring_applied = bool(getattr(result, "spectator_scoring_applied", False))
        rows_with_groups = int(getattr(result, "spectator_rows_with_groups", 0) or 0)
        rows_total = int(getattr(result, "spectator_rows_total", 0) or 0)
        sim_avg = float(getattr(result, "spectator_similarity_avg", 0.0) or 0.0)
        sim_range = getattr(result, "spectator_similarity_range", (0.0, 0.0)) or (0.0, 0.0)
        try:
            sim_min = float(sim_range[0])
            sim_max = float(sim_range[1])
        except Exception:
            sim_min = 0.0
            sim_max = 0.0
        displayed_scores = [
            float(getattr(rec, "spectator_score", 0.0) or 0.0)
            for rec in recs
        ]

        if detected_reaction_type:
            if isinstance(detected_confidence, (int, float)):
                detected_display = f"{detected_reaction_type} ({detected_confidence:.2f})"
            else:
                detected_display = str(detected_reaction_type)
        else:
            detected_display = "None"
        summary_lines = [
            html.escape(f"Query Reaction Key: {query_key or 'None'}"),
            html.escape(f"Matched Reaction Key: {matched_key or 'None'}"),
            html.escape(f"Reaction Type Filter: {reaction_filter or 'None'}"),
            html.escape(f"Detected Reaction Type: {detected_display}"),
            html.escape(f"Reacted Motifs: {', '.join(reacted_motifs) if reacted_motifs else 'None'}"),
            html.escape(f"Formed Motifs: {', '.join(formed_motifs) if formed_motifs else 'None'}"),
            html.escape(f"Spectator Motifs: {', '.join(spectator_motifs) if spectator_motifs else 'None'}"),
        ]
        summary_lines.extend(
            [
                "<br><b>Reaction Key Generation</b>",
                html.escape(f"CRK: {_format_items(reaction_info.get('reaction_key'))}"),
                html.escape(f"Broad Product Tags: {_format_items(reaction_info.get('product_broad_tags'))}"),
                html.escape(f"Product Motifs (Reactive): {_format_items(reaction_info.get('product_motifs_reactive'))}"),
                html.escape(f"Reacted motifs: {_format_items(reaction_info.get('reacted_motifs'), sep='|')}"),
                html.escape(f"Formed motifs (all): {_format_items(reaction_info.get('formed_motifs'))}"),
                html.escape(f"Formed motifs (reaction-center): {_format_items(reaction_info.get('formed_motifs_center'))}"),
                html.escape(f"Formed motifs (context): {_format_items(reaction_info.get('formed_motifs_context'))}"),
                html.escape(f"Bond formed: {_format_items(reaction_info.get('bond_formed'))}"),
                html.escape(f"Bond broken: {_format_items(reaction_info.get('bond_broken'))}"),
                html.escape(f"Spectator groups: {_format_items(reaction_info.get('spectator_groups'))}"),
                html.escape(f"Spectator motifs: {_format_items(reaction_info.get('spectator_motifs'))}"),
                "<br><b>Spectator Ranking</b>",
                html.escape(f"Query spectator groups: {_format_items(query_spectator_groups)}"),
                html.escape(
                    "Scoring formula: match_score *= (1 - 0.70) + 0.70 * spectator_similarity"
                ),
                html.escape(f"Spectator scoring applied: {'Yes' if scoring_applied else 'No'}"),
                html.escape(f"Matched rows with spectator groups: {rows_with_groups}/{rows_total}"),
                html.escape(
                    f"Row spectator similarity: avg={sim_avg:.3f}, range=[{sim_min:.3f}, {sim_max:.3f}]"
                ),
                html.escape(
                    f"Displayed recommendation spectator scores: {_format_score_stats(displayed_scores)}"
                ),
            ]
        )
        if not query_key:
            summary_lines.append(html.escape(f"Reactant Types (fallback): A={type_a or 'None'} | B={type_b or 'None'}"))
        self.summary.setHtml("<br>".join(summary_lines))

        source_map = getattr(result, "recommendations_by_source", {}) or {}
        normalized_map: Dict[str, List[Any]] = {}
        for key, items in source_map.items():
            normalized_key = _normalize_source_group_label(key)
            normalized_map.setdefault(normalized_key, []).extend(items)
        if "precedent" in normalized_map and "similarity" not in normalized_map:
            normalized_map["similarity"] = list(normalized_map.pop("precedent") or [])
        source_map = normalized_map
        base_groups = ["literature", "rules", "motif", "similarity"]
        extra_groups = [g for g in sorted(source_map) if g not in base_groups]

        self._all_json_output = None
        self._last_result_obj = result
        self._last_export_context = {
            "reaction_smiles": reaction_smiles,
            "reaction_type_filter": self.reaction_filter_edit.text().strip() or None,
            "catalyst_filter": self.catalyst_filter_edit.text().strip() or None,
        }
        self.export_json_button.setEnabled(True)
        self.results_tabs.clear()

        added_any_tab = False
        for source_group in base_groups + extra_groups:
            group_recs = list(source_map.get(source_group) or [])
            if not group_recs:
                continue
            group_table = self._create_results_table()
            group_columns = _table_columns_for_type(source_group)
            self._populate_table(group_table, group_recs, group_columns)
            label = _tab_label_for_source_group(source_group, self._strategy_label)
            self.results_tabs.addTab(group_table, label)
            added_any_tab = True

        if not added_any_tab:
            fallback_table = self._create_results_table()
            self._populate_table(fallback_table, recs, _table_columns_for_type(data_type))
            self.results_tabs.addTab(fallback_table, "Results")

    def _export_json_output(self) -> None:
        if self._all_json_output is None:
            if self._last_result_obj is None:
                QtWidgets.QMessageBox.information(self, "Export JSON", "Run a recommendation first.")
                return
            try:
                from chemtools.formatters import format_hte_output

                result = self._last_result_obj
                source_map = getattr(result, "recommendations_by_source", {}) or {}
                normalized_map: Dict[str, List[Any]] = {}
                for key, items in source_map.items():
                    normalized_key = _normalize_source_group_label(key)
                    normalized_map.setdefault(normalized_key, []).extend(list(items or []))
                literature_recs = list(normalized_map.get("literature") or [])
                if not literature_recs:
                    literature_recs = list(getattr(result, "recommendations", []) or [])
                self._all_json_output = format_hte_output(
                    result,
                    recommendations=literature_recs,
                    reaction_smiles=str(self._last_export_context.get("reaction_smiles") or ""),
                    reaction_type_filter=self._last_export_context.get("reaction_type_filter"),
                    catalyst_filter=self._last_export_context.get("catalyst_filter"),
                    explanation=None,
                )
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self, "Export JSON", f"Failed to build export JSON: {exc}")
                return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export JSON",
            str(PROJECT_ROOT),
            "JSON Files (*.json);;All Files (*)",
        )
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(self._all_json_output, handle, indent=2, ensure_ascii=False)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Export JSON", f"Failed to save JSON: {exc}")

    def _show_reaction_image(self, reaction_smiles: str) -> None:
        if not reaction_smiles:
            return
        if self._reaction_dialog:
            try:
                self._reaction_dialog.close()
            except RuntimeError:
                pass
            self._reaction_dialog = None

        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("Reaction Preview")
        dialog.setWindowFlags(
            dialog.windowFlags()
            | QtCore.Qt.WindowType.WindowSystemMenuHint
            | QtCore.Qt.WindowType.WindowMinimizeButtonHint
        )
        dialog.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
        dialog.destroyed.connect(lambda: setattr(self, "_reaction_dialog", None))

        layout = QtWidgets.QVBoxLayout(dialog)
        info_label = QtWidgets.QLabel("")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)

        images_widget = QtWidgets.QWidget()
        images_layout = QtWidgets.QHBoxLayout(images_widget)
        images_layout.setContentsMargins(0, 0, 0, 0)
        images_layout.setSpacing(12)
        layout.addWidget(images_widget)

        try:
            from chemtools.visualization import render_molecule_image, render_reaction_image
        except Exception as exc:
            info_label.setText(f"Visualization unavailable: {exc}")
            dialog.show()
            self._reaction_dialog = dialog
            return

        reactant_a, reactant_b, product = _parse_reaction_smiles(reaction_smiles)
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            if product:
                info_label.setText("")
                output_path = temp_path / "reaction.png"
                try:
                    render_reaction_image(
                        reaction_smiles,
                        output_path,
                        size=_REACTION_PREVIEW_RENDER_SIZE,
                    )
                    pixmap = _fit_pixmap_to_screen(QtGui.QPixmap(str(output_path)))
                except Exception as exc:
                    info_label.setText(f"Unable to render reaction image: {exc}")
                else:
                    image_label = QtWidgets.QLabel()
                    image_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                    image_label.setPixmap(pixmap)
                    images_layout.addWidget(image_label)
            else:
                info_label.setText("Product SMILES missing; showing reactants only.")
                for label, smiles in (("A", reactant_a), ("B", reactant_b)):
                    if not smiles:
                        continue
                    output_path = temp_path / f"reactant_{label.lower()}.png"
                    try:
                        render_molecule_image(
                            smiles,
                            output_path,
                            size=_MOLECULE_PREVIEW_RENDER_SIZE,
                            legend=f"Reactant {label}",
                        )
                        pixmap = _fit_pixmap_to_screen(QtGui.QPixmap(str(output_path)))
                    except Exception as exc:
                        info_label.setText(f"Unable to render reactant images: {exc}")
                        break
                    image_label = QtWidgets.QLabel()
                    image_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                    image_label.setPixmap(pixmap)
                    images_layout.addWidget(image_label)

        if dialog.layout():
            dialog.layout().setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetFixedSize)
        dialog.adjustSize()
        dialog.show()
        self._reaction_dialog = dialog


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    window = HTERecommenderWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
