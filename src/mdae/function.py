"""Layer C -> D bridge (slice-1 stub): free concentration -> functional effect.

Kept deliberately minimal for the first slice. The point it exists to prove: the
functional readout depends on the FREE concentration at the target (``c_free``,
delivered by Layer B) and on a target constant (``Kd`` / threshold) that is a
property of the drug-target chemistry and does NOT change with resistance that
acts through exposure (e.g. PfCRT). Exposure and binding are separate knobs.

For chloroquine the relevant effect is interference with heme detoxification in
the DV; we model it two ways:

* fractional target occupancy ``theta = c_free / (c_free + Kd)`` (mass-action), and
* a coarse functional threshold: killing requires ``c_free`` above ``c_threshold``.

Both ``Kd`` and ``c_threshold`` are ILLUSTRATIVE for slice 1 (provenance: not fitted).
"""

from __future__ import annotations


def occupancy(c_free: float, Kd: float) -> float:
    """Fractional target occupancy theta = c_free / (c_free + Kd). c_free, Kd in uM."""
    return c_free / (c_free + Kd)


def above_threshold(c_free: float, c_threshold: float) -> bool:
    """Coarse Layer-C->D gate: is free concentration above the functional threshold?"""
    return c_free >= c_threshold
