"""Verified, investigation-local receipts for completed steps of an unfinished route check."""

from __future__ import annotations

import hashlib
from contextlib import contextmanager
from contextvars import ContextVar
from collections.abc import Iterator
from typing import Any

from ..core.store import InvestigationStore, canonical_bytes


_REPLAY = ContextVar("recompute_route_steps", default=False)


@contextmanager
def recompute_route_steps() -> Iterator[None]:
    """Keep explicit workspace replay independent of completed step receipts."""
    token = _REPLAY.set(True)
    try:
        yield
    finally:
        _REPLAY.reset(token)


class RouteCheckpoints:
    """Resume only the identical request, keeping completed route checks reproducible."""

    def __init__(self, store: InvestigationStore, request: dict[str, Any]) -> None:
        self.store = store
        self.key = hashlib.sha256(canonical_bytes(request)).hexdigest()
        self.references: list[str] = []
        self.saved: dict[int, tuple[str, dict[str, Any]]] = {}
        self.existing: set[str] = set()
        for event in store.events(verify_artifacts=False):
            if event.kind != "route_step_checkpoint":
                continue
            payload = store.read_artifact(event.artifact_ref)
            if payload["route_check_key"] == self.key:
                self.existing.add(event.artifact_ref)
                result = payload["result"]
                review = result["step_review"]
                retry = (review["precedent"].get("status") == "unavailable" or
                         review.get("exact_source_mapping", {}).get("status") == "unavailable")
                if retry:
                    self.saved.pop(payload["ordinal"], None)
                elif not _REPLAY.get():
                    self.saved[payload["ordinal"]] = (event.artifact_ref, result)

    def record(self, ordinal: int, result: dict[str, Any]) -> None:
        """Publish an independently verifiable step before evaluating the next one."""
        payload = {"schema_version": "route_step_checkpoint.v1", "route_check_key": self.key,
                   "ordinal": ordinal, "result": result}
        # Find an identical prior receipt even when a completed route is being replayed.
        reference = hashlib.sha256(canonical_bytes(payload)).hexdigest()
        if f"sha256:{reference}" not in self.existing:
            event = self.store.append("route_step_checkpoint", payload)
            self.references.append(event.artifact_ref)
        else:
            self.references.append(f"sha256:{reference}")
