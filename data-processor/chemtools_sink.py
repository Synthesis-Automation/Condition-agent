#!/usr/bin/env python3
from __future__ import annotations

"""
Utilities to export data-processor artifacts to chemtools' data locations.

Behavior:
- Dataset JSONL: Copy a produced JSONL to chemtools dataset dir so the API can read it.
  Destination dir resolution order:
    1) Explicit dataset_dir arg, if provided
    2) Env CHEMTOOLS_DATASET_DIR, if set
    3) Repo default: <repo_root>/data/reaction_dataset

- Registry JSONL: Copy a registry JSONL to chemtools registry path so resolver can read it.
  Destination path resolution order:
    1) Explicit registry_out arg, if provided
    2) Env CHEMTOOLS_REGISTRY_PATH, if set
    3) Repo default: <repo_root>/data/cas_registry_merged.jsonl

All operations are best-effort; errors are returned as strings for the caller to log.
"""

import os
import shutil
from typing import Optional, Tuple


def _repo_root_from_here(here: str) -> str:
    # data-processor lives directly under the repo root
    return os.path.dirname(os.path.abspath(here))


def export_jsonl_for_chemtools(jsonl_path: str, dataset_dir: Optional[str] = None) -> Tuple[bool, str]:
    try:
        src = os.path.abspath(jsonl_path)
        if not os.path.exists(src):
            return False, f"source JSONL does not exist: {src}"

        # Resolve destination directory
        if dataset_dir:
            dst_dir = dataset_dir
        else:
            env_dir = os.environ.get("CHEMTOOLS_DATASET_DIR", "").strip()
            if env_dir:
                dst_dir = env_dir
            else:
                repo_root = _repo_root_from_here(os.path.dirname(__file__))
                dst_dir = os.path.join(repo_root, "data", "reaction_dataset")

        dst_dir = os.path.abspath(dst_dir)
        os.makedirs(dst_dir, exist_ok=True)

        # Name by source filename; chemtools only requires .jsonl in the folder
        dst = os.path.join(dst_dir, os.path.basename(src))
        shutil.copy2(src, dst)
        return True, dst
    except Exception as e:
        return False, f"export_jsonl_for_chemtools failed: {e}"


def export_registry_for_chemtools(registry_path: str, registry_out: Optional[str] = None) -> Tuple[bool, str]:
    try:
        src = os.path.abspath(registry_path)
        if not os.path.exists(src):
            return False, f"source registry does not exist: {src}"

        # Resolve destination path
        if registry_out:
            dst = registry_out
        else:
            env_path = os.environ.get("CHEMTOOLS_REGISTRY_PATH", "").strip()
            if env_path:
                dst = env_path
            else:
                repo_root = _repo_root_from_here(os.path.dirname(__file__))
                dst = os.path.join(repo_root, "data", "cas_registry_merged.jsonl")

        dst = os.path.abspath(dst)
        os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
        shutil.copy2(src, dst)
        return True, dst
    except Exception as e:
        return False, f"export_registry_for_chemtools failed: {e}"

