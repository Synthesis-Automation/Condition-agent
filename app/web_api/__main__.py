"""Run the local reaction-recommender API and compiled frontend."""

from __future__ import annotations

import argparse

import uvicorn

from .main import create_app
from .runtime import LocalRecommendationRuntime


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--index", default=None)
    args = parser.parse_args()
    uvicorn.run(
        create_app(runtime=LocalRecommendationRuntime(args.index)),
        host="127.0.0.1",
        port=args.port,
    )


if __name__ == "__main__":
    main()
