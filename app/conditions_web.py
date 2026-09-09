"""Serve only the condition-recommendation workflow for a Linux deployment.

Run locally with ``python -m app.conditions_web`` or under a process manager
with ``uvicorn app.conditions_web:app``.
"""

from __future__ import annotations

import argparse

from app.web_api.main import create_app


app = create_app(recommendation_only=True)


def main() -> None:
    """Start the focused UI and API on a configurable interface."""
    import uvicorn

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
