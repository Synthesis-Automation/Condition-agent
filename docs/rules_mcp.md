# MCP Rules Integration

This document summarises how the vendored MCP rules server is exposed to the
Condition-Agent orchestrator.

1. The server lives in `mcp_rules_server/server.py` and loads rule families from
   JSON files (default: `data/rules/buchwald_cn.json`).
2. The client adapter in `clients/mcp_rules_client.py` launches the server as a
   subprocess and communicates via newline-delimited JSON-RPC.
3. The orchestrator consumes the client through `condition_agent.services.rules_service.RulesService`.
4. Helper utilities in `condition_agent.features.mapping` translate substrate
   analysis output to the features schema expected by the server.
5. Examples and tests are available in `examples/demo_rules_cli.py` and
   `tests/test_rules_integration.py`.
