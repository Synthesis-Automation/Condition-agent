"""Start Condition Desk; use --workbench to enable the research tools."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

import uvicorn

from .main import DEFAULT_FRONTEND_DIST, create_app
from .runtime import LocalRecommendationRuntime


def main() -> None:
    """Launch one frontend with its matching API profile."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--index", default=None)
    parser.add_argument("--scientific-chat", action="store_true", help="Enable local agent conversations at /scientific")
    parser.add_argument("--chat-root", default="results/ai_native/conversations")
    parser.add_argument("--chat-artifacts", default="examples/ai_native/artifacts.local.example.json")
    parser.add_argument("--codex", default=None, help="Native Codex executable; otherwise auto-discover")
    parser.add_argument("--agent-model", default=None, help="Optional model override; otherwise use Codex configuration")
    parser.add_argument("--agent-timeout", type=float, default=900)
    parser.add_argument(
        "--workbench",
        action="store_true",
        help="Open the research workbench and enable its API",
    )
    parser.add_argument(
        "--build", action="store_true", help="Rebuild the frontend before starting"
    )
    args = parser.parse_args()
    scientific_service = None
    if args.scientific_chat:
        if args.host not in {"127.0.0.1", "localhost", "::1"}:
            parser.error("Scientific chat currently supports loopback hosts only")
        from chem_coworker.scientific_workspace.agent_runtime import CodexRuntime
        from chem_coworker.scientific_workspace.conversation import ConversationService

        try:
            scientific_service = ConversationService(
                args.chat_root,
                runtime=CodexRuntime(executable=args.codex, model=args.agent_model, timeout_seconds=args.agent_timeout),
                artifacts=json.loads(Path(args.chat_artifacts).read_text("utf-8")),
            )
        except (OSError, ValueError, RuntimeError) as exc:
            parser.error(str(exc))
        args.workbench = True
    if args.build:
        npm = shutil.which("npm.cmd") or shutil.which("npm")
        if npm is None:
            parser.error("npm is required for --build; install Node.js first")
        subprocess.run(
            [npm, "run", "build"], cwd=DEFAULT_FRONTEND_DIST.parent, check=True
        )
    entry = "workbench.html" if args.workbench else "index.html"
    if not (DEFAULT_FRONTEND_DIST / entry).is_file() and not args.scientific_chat:
        parser.error("Frontend build is missing. Run python -m app.web_api --build")
    try:
        uvicorn.run(
            create_app(
                runtime=LocalRecommendationRuntime(args.index),
                recommendation_only=not args.workbench,
                scientific_service=scientific_service,
            ),
            host=args.host,
            port=args.port,
        )
    finally:
        if scientific_service is not None:
            scientific_service.close()


if __name__ == "__main__":
    main()
