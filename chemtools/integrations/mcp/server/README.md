# MCP Rules Server

This directory vendors a lightweight JSON-RPC server that exposes the
Condition Rule Library through a Model Context Protocol–friendly interface.
The server is designed to be launched as a subprocess and controlled via
stdin/stdout, keeping deployments hermetic and testable.

## Running Locally

```bash
python chemtools/integrations/mcp/server/server.py --rules data/rules/buchwald_cn.json
```

Once started, the server accepts newline-delimited JSON-RPC requests. Supported
methods are:

- `ping`
- `rules.apply`
- `rules.merge`
- `rules.audit`

See `plan/MCP_Server_Integration_Plan.md` for payload expectations.
