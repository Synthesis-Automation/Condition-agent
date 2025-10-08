from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, TYPE_CHECKING

from chemtools import smiles as chem_smiles

try:
    import numpy as np
except Exception:  # pragma: no cover - numpy should be present with molpipeline extras
    np = None  # type: ignore

try:  # Optional dependency: molpipeline and rdkit stack
    from molpipeline import Pipeline as MolPipeline
    from molpipeline.any2mol import AutoToMol
    from molpipeline.mol2any import MolToMorganFP, MolToRDKitPhysChem
    from molpipeline.mol2mol import ElementFilter, SaltRemover

    _MOLPIPELINE_IMPORT_ERROR: Optional[Exception] = None
except Exception as exc:  # pragma: no cover - executed only when molpipeline missing
    MolPipeline = None  # type: ignore[assignment]
    AutoToMol = MolToMorganFP = MolToRDKitPhysChem = ElementFilter = SaltRemover = None  # type: ignore[assignment]
    _MOLPIPELINE_IMPORT_ERROR = exc

if TYPE_CHECKING:  # pragma: no cover - typing only
    from numpy.typing import NDArray as _NDArray
else:  # pragma: no cover - runtime aliasing
    _NDArray = Any


@dataclass(frozen=True)
class MolPipelineEnvironment:
    """Snapshot of the molpipeline runtime environment."""

    available: bool
    version: Optional[str]
    rdkit_version: Optional[str]
    sklearn_version: Optional[str]
    shap_version: Optional[str]
    import_error: Optional[str]


def is_available() -> bool:
    """Return True when molpipeline was imported successfully."""

    return MolPipeline is not None


def environment_snapshot() -> MolPipelineEnvironment:
    """Collect basic version metadata for reporting to users/tests."""

    if not is_available():
        return MolPipelineEnvironment(
            available=False,
            version=None,
            rdkit_version=None,
            sklearn_version=None,
            shap_version=None,
            import_error=str(_MOLPIPELINE_IMPORT_ERROR) if _MOLPIPELINE_IMPORT_ERROR else None,
        )

    import importlib

    def _version(mod_name: str) -> Optional[str]:
        try:
            mod = importlib.import_module(mod_name)
        except Exception:
            return None
        return getattr(mod, "__version__", None)

    return MolPipelineEnvironment(
        available=True,
        version=_version("molpipeline"),
        rdkit_version=_version("rdkit"),
        sklearn_version=_version("sklearn"),
        shap_version=_version("shap"),
        import_error=None,
    )


def _require_molpipeline() -> None:
    if not is_available():
        raise RuntimeError(
            "molpipeline is not installed. Install chemtools with the optional extras or run "
            "`pip install molpipeline` first."
        )
    if np is None:
        raise RuntimeError("numpy is required when using MolPipeline integrations.")


def _coerce_to_smiles(items: Iterable[Any]) -> List[str]:
    """Best-effort conversion of ChemTools reagent items into SMILES strings."""

    out: List[str] = []
    for item in items:
        if item is None:
            continue
        if isinstance(item, str):
            s = item.strip()
            if s:
                out.append(s)
            continue
        if isinstance(item, Mapping):
            for key in ("smiles", "SMILES", "canonical_smiles", "isomeric_smiles"):
                val = item.get(key)  # type: ignore[arg-type]
                if isinstance(val, str) and val.strip():
                    out.append(val.strip())
                    break
            else:
                name = item.get("name")
                if isinstance(name, str) and name.strip():
                    norm = chem_smiles.normalize(name.strip())
                    candidate = (
                        norm.get("smiles_norm")
                        or norm.get("largest_smiles")
                        or name.strip()
                    )
                    if isinstance(candidate, str) and candidate.strip():
                        out.append(candidate.strip())
            continue
    return out


def build_standardization_steps() -> List[Tuple[str, Any]]:
    """Return molpipeline steps that standardize molecules consistently."""

    _require_molpipeline()
    assert AutoToMol is not None and ElementFilter is not None and SaltRemover is not None
    return [
        ("auto_to_mol", AutoToMol()),
        ("element_filter", ElementFilter()),
        ("salt_remover", SaltRemover()),
    ]


def build_physchem_pipeline(
    descriptor_list: Optional[Sequence[str]] = None,
    *,
    n_jobs: int = 1,
) -> Any:
    """Create a pipeline that emits selected RDKit physicochemical descriptors."""

    _require_molpipeline()
    assert MolPipeline is not None and MolToRDKitPhysChem is not None
    steps = build_standardization_steps()
    steps.append(
        (
            "physchem",
            MolToRDKitPhysChem(
                descriptor_list=list(descriptor_list) if descriptor_list else None,
                standardizer=None,
            ),
        )
    )
    return MolPipeline(steps, n_jobs=n_jobs)


def build_morgan_pipeline(
    *,
    n_bits: int = 2048,
    radius: int = 2,
    n_jobs: int = 1,
) -> Any:
    """Create a pipeline that returns Morgan fingerprints as numpy arrays."""

    _require_molpipeline()
    assert MolPipeline is not None and MolToMorganFP is not None
    steps = build_standardization_steps()
    steps.append(
        (
            "morgan_fp",
            MolToMorganFP(n_bits=n_bits, radius=radius),
        )
    )
    return MolPipeline(steps, n_jobs=n_jobs)


def transform_smiles(smiles_like: Iterable[Any], pipeline: Any) -> Any:
    """Normalize arbitrary inputs to SMILES and run the provided molpipeline pipeline."""

    _require_molpipeline()
    smiles_seq = _coerce_to_smiles(smiles_like)
    if not smiles_seq:
        raise ValueError("No SMILES strings could be extracted from the provided input.")
    return pipeline.transform(smiles_seq)


def fit_pipeline(smiles_like: Iterable[Any], y: Sequence[float], pipeline: Any) -> Any:
    """Fit a molpipeline model using ChemTools-style inputs."""

    _require_molpipeline()
    smiles_seq = _coerce_to_smiles(smiles_like)
    if not smiles_seq:
        raise ValueError("No SMILES strings could be extracted from the provided input.")
    pipeline.fit(smiles_seq, y)
    return pipeline


def _to_dense_array(matrix: Any) -> _NDArray:
    """Convert MolPipeline outputs (dense or sparse) into numpy arrays."""

    _require_molpipeline()
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=float)


def _extract_reagents(reaction: Mapping[str, Any]) -> List[Any]:
    if not isinstance(reaction, Mapping):
        return []
    if isinstance(reaction.get("reagents"), Sequence):
        return list(reaction["reagents"])  # type: ignore[index]
    if isinstance(reaction.get("Reagents"), Sequence):
        return list(reaction["Reagents"])  # type: ignore[index]
    raw = reaction.get("Reagent")
    if isinstance(raw, str):
        try:
            parsed = json.loads(raw)
            if isinstance(parsed, Sequence):
                return list(parsed)
        except Exception:
            return []
    if isinstance(raw, Sequence):
        return list(raw)
    return []


def collect_role_smiles(reaction: Mapping[str, Any]) -> Dict[str, List[str]]:
    """Collect SMILES grouped by reagent role from a reaction-like mapping."""

    roles: Dict[str, List[str]] = {}
    for reagent in _extract_reagents(reaction):
        if not isinstance(reagent, Mapping):
            continue
        role = (
            reagent.get("role")
            or reagent.get("Role")
            or reagent.get("category")
            or reagent.get("Category")
        )
        if not role:
            continue
        role_norm = str(role).strip().upper()
        if not role_norm:
            continue
        smi_list = _coerce_to_smiles([reagent])
        if smi_list:
            roles.setdefault(role_norm, []).extend(smi_list)
    return roles


class MolPipelineRoleAggregator:
    """Aggregate role-specific molecule features using pre-built MolPipeline objects."""

    def __init__(
        self,
        role_pipelines: Mapping[str, Any],
        *,
        aggregate: str = "mean",
        missing_strategy: str = "zeros",
        fallback_smiles: str = "C",
    ) -> None:
        _require_molpipeline()
        if np is None:
            raise RuntimeError("numpy is required when using MolPipeline integrations.")
        if not role_pipelines:
            raise ValueError("role_pipelines must not be empty.")
        self.role_pipelines = {role.upper(): pipe for role, pipe in role_pipelines.items()}
        self.aggregate = aggregate
        self.missing_strategy = missing_strategy
        self.fallback_smiles = fallback_smiles
        self.role_order = list(self.role_pipelines.keys())
        self._zero_vectors: Dict[str, _NDArray] = {}
        for role, pipeline in self.role_pipelines.items():
            arr = pipeline.transform([fallback_smiles])
            vec = self._flatten_feature(arr)
            self._zero_vectors[role] = np.zeros_like(vec, dtype=float)

    def _flatten_feature(self, matrix: Any) -> _NDArray:
        vec = _to_dense_array(matrix)
        if vec.ndim == 0:
            return vec.reshape(1)
        if vec.ndim == 1:
            return vec.astype(float, copy=True)
        vec = vec.reshape(vec.shape[0], -1)
        return vec[0].astype(float, copy=True)

    def _aggregate_matrix(self, matrix: Any) -> _NDArray:
        mat = _to_dense_array(matrix)
        if mat.ndim == 0:
            return mat.reshape(1)
        if mat.ndim == 1:
            return mat.astype(float, copy=True)
        mat = mat.reshape(mat.shape[0], -1)
        if mat.shape[0] == 0:
            return np.zeros(mat.shape[1], dtype=float)
        if self.aggregate == "mean":
            return mat.mean(axis=0)
        if self.aggregate == "sum":
            return mat.sum(axis=0)
        if self.aggregate == "first":
            return mat[0].astype(float, copy=True)
        raise ValueError(f"Unsupported aggregate strategy: {self.aggregate}")

    def featurize_roles(
        self,
        *,
        reaction: Optional[Mapping[str, Any]] = None,
        role_smiles: Optional[Mapping[str, Sequence[Any]]] = None,
    ) -> Dict[str, Optional[_NDArray]]:
        if np is None:
            raise RuntimeError("numpy is required when using MolPipeline integrations.")
        if role_smiles is None:
            if reaction is None:
                raise ValueError("Either reaction or role_smiles must be provided.")
            role_smiles = collect_role_smiles(reaction)
        result: Dict[str, Optional[_NDArray]] = {}
        for role in self.role_order:
            smiles_seq: Sequence[Any] = []
            if isinstance(role_smiles, Mapping):
                for key, value in role_smiles.items():
                    if str(key).strip().upper() == role:
                        smiles_seq = value
                        break
            smiles_list = _coerce_to_smiles(smiles_seq)
            if not smiles_list:
                if self.missing_strategy == "zeros":
                    result[role] = self._zero_vectors[role].copy()
                else:
                    result[role] = None
                continue
            matrix = self.role_pipelines[role].transform(smiles_list)
            result[role] = self._aggregate_matrix(matrix)
        return result

    def concatenate(
        self,
        *,
        reaction: Optional[Mapping[str, Any]] = None,
        role_smiles: Optional[Mapping[str, Sequence[Any]]] = None,
    ) -> _NDArray:
        features = self.featurize_roles(reaction=reaction, role_smiles=role_smiles)
        if np is None:
            raise RuntimeError("numpy is required when using MolPipeline integrations.")
        vectors: List[_NDArray] = []
        for role in self.role_order:
            vec = features.get(role)
            if vec is None:
                if self.missing_strategy == "zeros":
                    vec = self._zero_vectors[role]
                else:
                    raise ValueError(f"No features for role {role} and missing_strategy != zeros")
            vectors.append(np.asarray(vec, dtype=float))
        return np.concatenate(vectors, axis=0)


def build_default_role_aggregator(
    *,
    roles: Optional[Sequence[str]] = None,
    aggregate: str = "mean",
    missing_strategy: str = "zeros",
    n_jobs: int = 1,
    ligand_n_bits: int = 512,
    ligand_radius: int = 2,
) -> MolPipelineRoleAggregator:
    """Create a role aggregator with common defaults (ligand/base/solvent)."""

    _require_molpipeline()
    role_order = [r.upper() for r in (roles or ("LIGAND", "BASE", "SOLVENT"))]
    pipelines: Dict[str, Any] = {}
    if "LIGAND" in role_order:
        pipelines["LIGAND"] = build_morgan_pipeline(
            n_bits=ligand_n_bits,
            radius=ligand_radius,
            n_jobs=n_jobs,
        )
    if "BASE" in role_order:
        pipelines["BASE"] = build_physchem_pipeline(
            ["HeavyAtomMolWt", "TPSA", "MolLogP", "NumHAcceptors"],
            n_jobs=n_jobs,
        )
    if "SOLVENT" in role_order:
        pipelines["SOLVENT"] = build_physchem_pipeline(
            ["TPSA", "MolLogP", "MolMR"],
            n_jobs=n_jobs,
        )
    if not pipelines:
        raise ValueError("No pipelines created for the requested roles.")
    return MolPipelineRoleAggregator(
        pipelines,
        aggregate=aggregate,
        missing_strategy=missing_strategy,
    )


__all__ = [
    "MolPipelineEnvironment",
    "environment_snapshot",
    "is_available",
    "build_physchem_pipeline",
    "build_morgan_pipeline",
    "build_standardization_steps",
    "transform_smiles",
    "fit_pipeline",
    "collect_role_smiles",
    "MolPipelineRoleAggregator",
    "build_default_role_aggregator",
]
