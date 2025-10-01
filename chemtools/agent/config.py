"""Configuration helpers for MCP rules integration."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RULES_PATH = REPO_ROOT / "data" / "rules" / "buchwald_cn.json"
DEFAULT_SERVER_BIN = f"{os.sys.executable} \"{REPO_ROOT / 'chemtools' / 'integrations' / 'mcp' / 'server' / 'server.py'}\""


@dataclass
class RulesConfig:
    rules_path: Path
    mcp_server_bin: str
    mcp_startup_timeout_s: float
    log_level: str = "info"


def load_config(env: Optional[dict[str, str]] = None) -> RulesConfig:
    """Build a :class:`RulesConfig` from environment variables."""

    env = dict(env or os.environ)
    rules_path = Path(env.get("RULES_PATH", str(DEFAULT_RULES_PATH))).expanduser().resolve()
    server_bin = env.get("MCP_SERVER_BIN", DEFAULT_SERVER_BIN)
    timeout = float(env.get("MCP_STARTUP_TIMEOUT_S", "10"))
    log_level = env.get("MCP_LOG_LEVEL", "info")
    return RulesConfig(rules_path=rules_path, mcp_server_bin=server_bin, mcp_startup_timeout_s=timeout, log_level=log_level)
