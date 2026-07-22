"""Traceable parameters with provenance and uncertainty.

Every quantity that enters a prediction is a :class:`Param` carrying its value, an
uncertainty, a unit, an epistemic KIND, and a citation. KIND is kept explicit per
the project's scientific contract:

- ``measured`` : taken directly from a primary experimental source (cite DOI/PMID).
- ``derived``  : computed from measured quantities by a stated formula.
- ``fitted``   : calibrated to data inside this model (declare against what).
- ``assumed``  : not measured for this system; a documented placeholder/bound.

Distributions for Monte-Carlo propagation:
- measured / derived / fitted -> Normal(value, sd), truncated positive
  (sd taken from ``sd`` or, if absent, from (high-low)/2 read as ~1 SEM).
- assumed                     -> Uniform(low, high)  (honest ignorance within bounds).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np


class Kind(str, Enum):
    MEASURED = "measured"
    DERIVED = "derived"
    FITTED = "fitted"
    ASSUMED = "assumed"


@dataclass(frozen=True)
class Param:
    """A single sourced parameter with uncertainty and provenance."""

    name: str
    value: float
    unit: str
    kind: Kind
    low: float | None = None   # lower bound (assumed) or value-uncertainty (measured)
    high: float | None = None  # upper bound (assumed) or value+uncertainty (measured)
    sd: float | None = None    # explicit standard deviation (measured), overrides (high-low)/2
    source: str = ""           # DOI / PMID / primary citation
    note: str = ""

    def _sd(self) -> float:
        if self.sd is not None:
            return float(self.sd)
        if self.low is not None and self.high is not None:
            return (self.high - self.low) / 2.0
        return 0.0

    def interval(self) -> tuple[float, float]:
        """Reporting interval: explicit [low, high] or value +/- sd."""
        if self.low is not None and self.high is not None:
            return (float(self.low), float(self.high))
        sd = self._sd()
        return (self.value - sd, self.value + sd)

    def sample(self, rng: np.random.Generator, n: int) -> np.ndarray:
        """Draw ``n`` samples according to KIND (see module docstring)."""
        if self.kind is Kind.ASSUMED:
            lo = self.low if self.low is not None else self.value
            hi = self.high if self.high is not None else self.value
            return rng.uniform(lo, hi, n)
        sd = self._sd()
        if sd == 0.0:
            return np.full(n, float(self.value))
        return np.clip(rng.normal(self.value, sd, n), 1e-12, None)
