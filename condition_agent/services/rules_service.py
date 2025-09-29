"""Service wrapper around the MCP rules client."""
from __future__ import annotations

import logging
from typing import Any, Dict, Iterable, Optional

from clients.mcp_rules_client import McpRulesClient, McpRulesClientError
from condition_agent.config import RulesConfig
from condition_agent.features.mapping import build_features

LOGGER = logging.getLogger(__name__)


class RulesService:
    """Thin orchestrator-facing wrapper."""

    def __init__(self, cfg: RulesConfig) -> None:
        self.cfg = cfg
        self._client = McpRulesClient(cfg.mcp_server_bin, cfg.rules_path, cfg.mcp_startup_timeout_s)
        self._client.start()

    # ------------------------------------------------------------------
    def suggest_conditions(self, reaction_smiles: str, reactants: Iterable[str], context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        features = build_features(list(reactants), context)
        LOGGER.debug("Calling rules.apply with features=%s", features)
        return self._client.apply(reaction_smiles, features)

    def merge_candidates(self, payload: Dict[str, Any], strategy: str = "append_as_candidate") -> Dict[str, Any]:
        LOGGER.info("Merging %s playbooks", len(payload.get("playbooks", [])))
        return self._client.merge(payload, strategy=strategy)

    def audit_rule(self, rule_id: str) -> Dict[str, Any]:
        return self._client.audit(rule_id)

    # ------------------------------------------------------------------
    def close(self) -> None:
        self._client.stop()


__all__ = ["RulesService", "McpRulesClientError"]
