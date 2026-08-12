# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Log-safe rates and atomic serial residence-time selection."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

from .canonical import digest
from .catalog import Event
from .errors import DomainFailure
from .rng import PhiloxStream


@dataclass(frozen=True, slots=True)
class RateTable:
    origin_state_id: str
    event_ids: tuple[str, ...]
    destination_state_ids: tuple[str, ...]
    log_rates: tuple[float, ...]
    rates: tuple[float, ...]
    total_rate: float
    lost_rate_log_upper_bound: float | None

    def snapshot(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "schema": "spark-atomistic-rate-table-snapshot/1",
            "origin_state_id": self.origin_state_id,
            "event_ids": list(self.event_ids),
            "destination_state_ids": list(self.destination_state_ids),
            "log_rates": list(self.log_rates), "rates": list(self.rates),
            "total_rate_per_s": self.total_rate,
            "lost_rate_log_upper_bound": self.lost_rate_log_upper_bound,
        }
        return {"payload": payload, "payload_sha256": digest(payload)}


def _neumaier_sum(values: Sequence[float]) -> float:
    total = 0.0
    correction = 0.0
    for value in values:
        updated = total + value
        if abs(total) >= abs(value):
            correction += (total - updated) + value
        else:
            correction += (value - updated) + total
        total = updated
    return total + correction


def build_rate_table(events: Sequence[Event], origin_state_id: str) -> RateTable:
    directed = sorted((event for event in events if event.origin_state_id == origin_state_id),
                      key=lambda event: event.event_id)
    enabled = [event for event in directed if event.enabled]
    disabled = [event for event in directed if not event.enabled]
    if not enabled:
        raise DomainFailure("NO_ENABLED_EVENT", "state has no enabled event",
                            component="kinetics", requirement="KMC-002", state_id=origin_state_id)
    maximum_log = max(event.log_rate for event in enabled)
    if not math.isfinite(maximum_log):
        raise DomainFailure("RATE_INVALID", "nonfinite log rate",
                            component="rates", requirement="RATE-004", state_id=origin_state_id)
    try:
        rates = tuple(math.exp(event.log_rate) for event in enabled)
        total = _neumaier_sum(rates)
    except (ArithmeticError, ValueError, OverflowError) as exc:
        raise DomainFailure("RATE_INVALID", "rate arithmetic/domain failure",
                            component="rates", requirement="RATE-004", state_id=origin_state_id) from exc
    if not math.isfinite(total) or total <= 0.0 or any(not math.isfinite(rate) or rate <= 0.0 for rate in rates):
        raise DomainFailure("RATE_INVALID", "nonpositive or nonfinite selectable rate",
                            component="rates", requirement="RATE-004", state_id=origin_state_id)
    lost_log = None
    if disabled:
        try:
            disabled_max = max(event.log_rate for event in disabled)
            lost_log = disabled_max + math.log(sum(math.exp(event.log_rate - disabled_max)
                                                   for event in disabled))
        except (ArithmeticError, ValueError, OverflowError) as exc:
            raise DomainFailure("RATE_INVALID", "lost-rate bound domain failure",
                                component="rates", requirement="RATE-004",
                                state_id=origin_state_id) from exc
    return RateTable(origin_state_id, tuple(event.event_id for event in enabled),
                     tuple(event.destination_state_id for event in enabled),
                     tuple(event.log_rate for event in enabled), rates, total, lost_log)


@dataclass(frozen=True, slots=True)
class Selection:
    event_id: str
    destination_state_id: str
    selection_uniform: float
    time_uniform: float
    total_rate: float
    selected_rate: float
    delta_time: float
    rng_after: PhiloxStream
    rate_table_snapshot: dict[str, object]


def propose_serial_step(table: RateTable, rng: PhiloxStream) -> Selection:
    """Consume only a clone; caller commits clone after event verification succeeds."""
    clone = rng.clone()
    selection_uniform = clone.uniform()
    time_uniform = clone.uniform()
    threshold = selection_uniform * table.total_rate
    cumulative = 0.0
    selected_index = len(table.rates) - 1
    for index, rate in enumerate(table.rates):
        cumulative += rate
        if cumulative > threshold:
            selected_index = index
            break
    delta_time = -math.log(time_uniform) / table.total_rate
    if not math.isfinite(delta_time) or delta_time <= 0.0:
        raise DomainFailure("RATE_INVALID", "KMC time increment is not finite and positive",
                            component="kinetics", requirement="KMC-005",
                            object_id=table.event_ids[selected_index])
    return Selection(table.event_ids[selected_index], table.destination_state_ids[selected_index],
                     selection_uniform, time_uniform,
                     table.total_rate, table.rates[selected_index], delta_time, clone,
                     table.snapshot())
