"""Acid-base speciation of a (poly)protic weak base — the basis of ion trapping.

Contract
--------
A weak base exists in protonation states differing by bound protons. Only the
*neutral* species permeates a membrane freely; charged states are trapped. The
neutral fraction as a function of pH sets how much drug accumulates in an acidic
compartment (Layer B). This is textbook Henderson-Hasselbalch — NOT novel; it is
the analytical backbone the engine is built on.

Convention
----------
``pKa`` is the list of protonation pKa values ordered from the FIRST proton added
to the neutral base (the most basic site, highest pKa) downward. For a base with
protonation constants p_1 >= p_2 >= ... the population of the state with ``k``
bound protons relative to the neutral form is::

    [state_k] / [neutral] = 10^( (p_1 + ... + p_k) - k*pH )

so the partition function (sum over all states, relative to neutral) is::

    Z(pH) = sum_{k=0..n} 10^( S_k - k*pH ),   S_0 = 0,  S_k = p_1 + ... + p_k

and the neutral fraction is ``f_neutral(pH) = 1 / Z(pH)``.

Chloroquine is diprotic: pKa ~ [10.1 (side-chain amine, protonates first),
8.1 (quinoline N)]. Provenance: experimental, literature (Warhurst et al.).
"""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


@dataclass(frozen=True)
class WeakBase:
    """A polyprotic weak base defined by its protonation pKa values.

    Parameters
    ----------
    name : str
    pKa : list[float]
        Protonation pKa values, ordered first-proton-first (most basic site,
        highest pKa) to last. Empty list = a neutral, non-ionizable molecule.
    """

    name: str
    pKa: tuple[float, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        pk = tuple(float(p) for p in self.pKa)
        object.__setattr__(self, "pKa", pk)
        # Warn (not fail) if not ordered high->low, since order sets the single-proton term.
        if any(pk[i] < pk[i + 1] for i in range(len(pk) - 1)):
            # Not an error, but the convention above assumes descending order.
            object.__setattr__(self, "_unordered", True)

    def _log_partition_terms(self, pH: np.ndarray | float) -> np.ndarray:
        """log10 of each term 10^(S_k - k*pH), returned as an array over k=0..n."""
        pH = np.asarray(pH, dtype=float)
        cumulative = np.concatenate([[0.0], np.cumsum(self.pKa)])  # S_0..S_n
        k = np.arange(len(cumulative))
        # broadcast: shape (..., n+1)
        return cumulative - np.multiply.outer(pH, k)

    def partition(self, pH: np.ndarray | float) -> np.ndarray:
        """Z(pH) = sum_k 10^(S_k - k*pH). Computed log-stably."""
        log_terms = self._log_partition_terms(pH)
        # log-sum-exp in base 10 for numerical stability at extreme pH
        m = np.max(log_terms, axis=-1, keepdims=True)
        z = np.power(10.0, m).squeeze(-1) * np.sum(np.power(10.0, log_terms - m), axis=-1)
        return z

    def neutral_fraction(self, pH: np.ndarray | float) -> np.ndarray:
        """Fraction of drug in the neutral (permeant) form at a given pH."""
        return 1.0 / self.partition(pH)

    def charged_fraction(self, pH: np.ndarray | float) -> np.ndarray:
        """Fraction of drug in any charged (trapped) form at a given pH."""
        return 1.0 - self.neutral_fraction(pH)
