"""Analytical oracle: steady-state ion-trapping accumulation ratio.

This is the *independent* check for the Layer-B ODE model. It is derived purely
from acid-base speciation (``speciation.py``), with no dynamics: at steady state
the neutral (permeant) species equilibrates across the membrane, so::

    [neutral]_in = [neutral]_out
    => C_total,in / C_total,out = f_neutral(pH_out) / f_neutral(pH_in)
                                 = Z(pH_in) / Z(pH_out)

For a diprotic base this is the classic result predicting ~10^3-10^4 chloroquine
accumulation in the digestive vacuole (Homewood/Warhurst, Nature 1972). The
dynamic flux-balance model in ``exposure.py`` MUST reproduce this number with the
transporter turned off — two independent derivations that have to agree.
"""

from __future__ import annotations

from .speciation import WeakBase


def accumulation_ratio(base: WeakBase, pH_inside: float, pH_outside: float) -> float:
    """Steady-state total-concentration ratio C_inside / C_outside (passive, transporter off).

    Parameters
    ----------
    base : WeakBase
    pH_inside : float
        pH of the accumulating compartment (e.g. digestive vacuole ~5.2).
    pH_outside : float
        pH of the source compartment (e.g. cytosol/plasma ~7.4).

    Returns
    -------
    float
        Accumulation ratio R = Z(pH_inside) / Z(pH_outside). R > 1 means the drug
        concentrates in the acidic compartment.
    """
    return float(base.partition(pH_inside) / base.partition(pH_outside))
