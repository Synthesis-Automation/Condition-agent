"""Subprocess-based MCP rules client."""
from __future__ import annotations

import atexit
import json
import os
import queue
import shlex
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union

CommandType = Union[str, Iterable[str]]


class McpRulesClientError(RuntimeError):
    """Base error for MCP rule client failures."""


@dataclass
class McpRulesClient:
    server_cmd: CommandType
    rules_path: Union[str, Path]
    startup_timeout_s: float = 10.0

    def __post_init__(self) -> None:
        self._proc: Optional[subprocess.Popen[str]] = None
        self._id = 0
        self._lock = threading.Lock()
        self._pending = queue.Queue()

    # ------------------------------------------------------------------
    # lifecycle helpers
    # ------------------------------------------------------------------
    def start(self) -> None:
        if self._proc is not None and self._proc.poll() is None:
            return

        env = os.environ.copy()
        env.setdefault("RULES_PATH", str(self.rules_path))

        if isinstance(self.server_cmd, str):
            cmd_parts = shlex.split(self.server_cmd, posix=os.name != "nt")
        else:
            cmd_parts = list(self.server_cmd)
        if not cmd_parts:
            raise McpRulesClientError("Empty server command")
        if "--rules" not in " ".join(cmd_parts):
            cmd_parts.extend(["--rules", str(self.rules_path)])

        self._proc = subprocess.Popen(  # noqa: S603,S607 - local subprocess
            cmd_parts,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            env=env,
        )
        atexit.register(self.stop)

        # Wait for ping response to ensure readiness
        deadline = time.time() + float(self.startup_timeout_s)
        try:
            while time.time() < deadline:
                try:
                    result = self._rpc("ping", {})
                    if result.get("status") == "ok":
                        return
                except McpRulesClientError:
                    time.sleep(0.1)
            raise TimeoutError("Timed out waiting for MCP rules server to initialize")
        except Exception:
            self.stop()
            raise

    def stop(self) -> None:
        if self._proc is None:
            return
        if self._proc.poll() is None:
            try:
                self._proc.terminate()
                self._proc.wait(timeout=3)
            except Exception:
                self._proc.kill()
        self._proc = None

    # ------------------------------------------------------------------
    # RPC helpers
    # ------------------------------------------------------------------
    def _rpc(self, method: str, params: Dict[str, Any]) -> Dict[str, Any]:
        if self._proc is None or self._proc.stdin is None or self._proc.stdout is None:
            raise McpRulesClientError("MCP rules server is not running")
        if self._proc.poll() is not None:
            raise McpRulesClientError("MCP rules server exited unexpectedly")

        with self._lock:
            self._id += 1
            request = {"jsonrpc": "2.0", "id": self._id, "method": method, "params": params}
            try:
                self._proc.stdin.write(json.dumps(request) + "\n")
                self._proc.stdin.flush()
            except Exception as exc:  # pragma: no cover - I/O edge cases
                raise McpRulesClientError(f"Failed to write to MCP server: {exc}") from exc

            line = self._proc.stdout.readline()
            if not line:
                raise McpRulesClientError("MCP rules server returned no data")
            try:
                response = json.loads(line)
            except json.JSONDecodeError as exc:
                raise McpRulesClientError(f"Invalid JSON from MCP server: {exc}") from exc

            if "error" in response:
                raise McpRulesClientError(str(response["error"]))
            return response.get("result", {})

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def apply(self, reaction_smiles: str, features: Dict[str, Any], max_suggestions: int = 5) -> Dict[str, Any]:
        return self._rpc(
            "rules.apply",
            {
                "reaction_smiles": reaction_smiles,
                "features": features,
                "max_suggestions": max_suggestions,
            },
        )

    def merge(self, candidates: Dict[str, Any], strategy: str = "append_as_candidate") -> Dict[str, Any]:
        payload = {"candidates": candidates, "strategy": strategy}
        return self._rpc("rules.merge", payload)

    def audit(self, rule_id: str) -> Dict[str, Any]:
        return self._rpc("rules.audit", {"rule_id": rule_id})


__all__ = ["McpRulesClient", "McpRulesClientError"]
