"""Start Condition Desk; use --workbench to enable the research tools."""

from __future__ import annotations

import argparse
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
    parser.add_argument(
        "--workbench",
        action="store_true",
        help="Open the research workbench and enable its API",
    )
    parser.add_argument(
        "--build", action="store_true", help="Rebuild the frontend before starting"
    )
    args = parser.parse_args()
    if args.build:
        npm = shutil.which("npm.cmd") or shutil.which("npm")
        if npm is None:
            parser.error("npm is required for --build; install Node.js first")
        subprocess.run(
            [npm, "run", "build"], cwd=DEFAULT_FRONTEND_DIST.parent, check=True
        )
    entry = "workbench.html" if args.workbench else "index.html"
    if not (DEFAULT_FRONTEND_DIST / entry).is_file():
        parser.error("Frontend build is missing. Run python -m app.web_api --build")
    uvicorn.run(
        create_app(
            runtime=LocalRecommendationRuntime(args.index),
            recommendation_only=not args.workbench,
        ),
        host=args.host,
        port=args.port,
    )


if __name__ == "__main__":
    main()
