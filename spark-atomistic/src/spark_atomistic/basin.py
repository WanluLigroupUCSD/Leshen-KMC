# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Errata-2 common basin capability: uniformly disabled."""

from __future__ import annotations

from typing import Any

from .errors import Outcome


DISABLED_CHECKPOINT_RECORD = {"enabled": False, "reason": "v1-disabled"}


def basin_disabled_outcome(state_id: str | None = None) -> Outcome:
    return Outcome("BASIN_DISABLED", "", {"component": "basin", "details": {},
                   "requirement_id": "E2-BASIN-002", "retryable": False,
                   "search_or_event_id": None, "state_id": state_id})


def verify_basin_completeness(*_args: Any, **_kwargs: Any) -> tuple[bool, str]:
    return False, "v1-disabled"


def sample_finite_basin(*_args: Any, **kwargs: Any) -> tuple[Outcome, None]:
    return basin_disabled_outcome(kwargs.get("current_state_id")), None
