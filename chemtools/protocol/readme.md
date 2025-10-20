# Protocol Database Tools

## CLI Commands

### Build/Rebuild Index
```bash
# Build index with DRFP fingerprints
python -m chemtools.protocol.cli build --force

# Build with custom directory
python -m chemtools.protocol.cli build --force --protocol-dir "c:\Git-softwares\Condition-agent\data\protocol_db"

# View index statistics
python -m chemtools.protocol.cli stats
```

### Validate Protocols
```bash
# Validate all protocols (check SMARTS matching)
python -m chemtools.protocol.validate_protocols

# Validate specific file
python -m chemtools.protocol.validate_protocols --file "Suzuki_Protocol.json"

# Show only errors
python -m chemtools.protocol.validate_protocols --errors-only

# Export validation report
python -m chemtools.protocol.validate_protocols --output validation_report.json

# Verbose output with details
python -m chemtools.protocol.validate_protocols --verbose

# Fail with exit code 1 if any protocols are invalid (useful for CI/CD)
python -m chemtools.protocol.validate_protocols --fail-on-error
```

## Notes

replace [a] with [c,C,n,o,s]


