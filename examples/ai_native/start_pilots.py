"""Record development baselines; the external agent chooses subsequent investigations.

Run from the repository: python -m examples.ai_native.start_pilots --output results/ai_native/run01
This script is a reproducibility aid, not an agent or an evaluation of agent quality.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from chem_coworker.scientific_workspace import ScientificWorkspace


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    """Prepare one workspace and save condition and route-search starting points."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--artifacts", default=str(Path(__file__).with_name("artifacts.local.example.json")))
    args = parser.parse_args()
    cases = json.loads(Path(__file__).with_name("development_cases.json").read_text("utf-8"))
    paths = json.loads(Path(args.artifacts).read_text("utf-8"))
    workspace = ScientificWorkspace.create(
        args.output, repository=ROOT, artifacts=paths,
        objective="Run the authored condition-transfer and route-revision development pilots",
        constraints=("Development only; no untouched evaluation claim", "Do not change chemistry rules during the run"),
        agent_metadata={"runtime": "external agent", "model": "not_recorded_by_script"},
    )
    workspace.store.attach_file(__file__, description="Pilot baseline invocation script")
    workspace.store.attach_file(Path(__file__).with_name("development_cases.json"), description="Authored development tasks")
    condition: dict[str, Any] = dict(cases["condition_transfer"])
    condition.pop("objective")
    for operation, arguments in (
        ("analyze_reaction", {"reaction_smiles": condition["reaction_smiles"]}),
        ("recommend_conditions", condition),
        ("plan_routes", {"settings": cases["route_revision"]["settings"]}),
    ):
        event = workspace.run(operation, arguments)
        print(json.dumps(workspace.call_summary(event)), flush=True)
    workspace.store.note("question", "Which retrieved substrate differences limit condition transfer?")
    workspace.store.note("question", "Which route issue warrants revising an earlier choice?")


if __name__ == "__main__":
    main()
