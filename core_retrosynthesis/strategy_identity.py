"""Deterministic identity for one single-step retrosynthetic strategy."""

from __future__ import annotations

from .chemistry import digest


STRATEGY_ID_NAMESPACE = "STRAT1"


def build_strategy_id(
    operator_id: str,
    disconnection_site_key: str,
    synthon_signature: str,
) -> str:
    """Identify a strategy independently of its concrete handle realization.

    The identity deliberately excludes templates, named reaction families,
    precursor handles, conditions, and source labels.  It is defined only by
    the validated graph operator, product disconnection site, and normalized
    precursor skeleton.
    """

    components = {
        "operator_id": operator_id,
        "disconnection_site_key": disconnection_site_key,
        "synthon_signature": synthon_signature,
    }
    missing = tuple(name for name, value in components.items() if not value)
    if missing:
        raise ValueError(
            "strategy identity requires nonempty " + ", ".join(missing)
        )
    return digest(
        STRATEGY_ID_NAMESPACE,
        operator_id,
        disconnection_site_key,
        synthon_signature,
    )


__all__ = ["STRATEGY_ID_NAMESPACE", "build_strategy_id"]
