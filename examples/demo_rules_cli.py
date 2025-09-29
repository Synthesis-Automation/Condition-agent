from __future__ import annotations

import json
import sys

from pathlib import Path

import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from condition_agent.config import load_config
from condition_agent.services.rules_service import RulesService


def main() -> int:
    cfg = load_config()
    service = RulesService(cfg)
    try:
        reaction = "Clc1ccc(Cl)cc1.Nc1ccccc1>>Clc1ccc(Nc2ccccc2)cc1"
        reactants = ["Clc1ccc(Cl)cc1", "Nc1ccccc1"]
        result = service.suggest_conditions(reaction, reactants, context={"mode": "batch"})
        print(json.dumps(result, indent=2, ensure_ascii=False))
    finally:
        service.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
