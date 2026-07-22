"""Layer B — parasite exposure (the 'second ADME'), vertical slice 1.

Minimal two-compartment ion-trapping model: a source compartment (cytosol/plasma,
treated as a clamped reservoir supplied by Layer A) and the acidic digestive
vacuole (DV), coupled by (1) passive permeation of the neutral species and (2) an
optional saturable transporter (PfCRT) that effluxes the charged species out of
the DV.

Governing equation (DV total aqueous concentration C, in uM; time in s)::

    dC/dt = k_perm * (f_ext * C_ext  -  f_dv * C)          # passive, neutral-only
            - vmax * (cf * C) / (Km + cf * C)              # PfCRT efflux of charged form

    f_ext = neutral_fraction(pH_ext),  f_dv = neutral_fraction(pH_dv)
    cf    = charged_fraction(pH_dv)          (~1 in the acidic DV)
    k_perm [1/s]  lumps permeability * area / volume (estimated)

Passive steady state (vmax = 0):  C* = (f_ext / f_dv) * C_ext = R * C_ext, i.e. it
recovers the analytical oracle exactly. With PfCRT on, C* drops BELOW R*C_ext with
NO change to any binding constant — resistance as an *exposure* phenomenon.

Units & provenance are the contract; see ``parameters.py``. All rate parameters
here are ESTIMATED/illustrative for slice 1, not fitted to chloroquine data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .speciation import WeakBase


@dataclass(frozen=True)
class Compartment:
    """A well-mixed aqueous compartment."""

    name: str
    pH: float


@dataclass(frozen=True)
class IonTrapModel:
    """Two-compartment weak-base ion-trapping model with optional PfCRT efflux.

    Parameters
    ----------
    base : WeakBase
    source : Compartment
        Reservoir feeding the DV (cytosol/plasma). Held at constant ``c_ext``.
    vacuole : Compartment
        The accumulating acidic organelle (digestive vacuole).
    c_ext : float
        Total drug concentration in the source compartment [uM] (from Layer A).
    k_perm : float
        Lumped passive permeation rate constant [1/s] = P*A/V (ESTIMATED).
    vmax : float
        Transporter (PfCRT) maximal efflux rate [uM/s] (0 = sensitive strain).
    Km : float
        Transporter half-saturation for the charged DV species [uM].
    """

    base: WeakBase
    source: Compartment
    vacuole: Compartment
    c_ext: float
    k_perm: float
    vmax: float = 0.0
    Km: float = 1.0

    # ---- speciation shortcuts ------------------------------------------------
    @property
    def f_ext(self) -> float:
        return float(self.base.neutral_fraction(self.source.pH))

    @property
    def f_dv(self) -> float:
        return float(self.base.neutral_fraction(self.vacuole.pH))

    @property
    def cf_dv(self) -> float:
        return float(self.base.charged_fraction(self.vacuole.pH))

    # ---- dynamics ------------------------------------------------------------
    def _dCdt(self, C: float) -> float:
        passive = self.k_perm * (self.f_ext * self.c_ext - self.f_dv * C)
        charged = self.cf_dv * C
        efflux = self.vmax * charged / (self.Km + charged) if self.vmax > 0 else 0.0
        return passive - efflux

    def passive_steady_state(self) -> float:
        """Analytic C* for the transporter-off case: R * C_ext."""
        return (self.f_ext / self.f_dv) * self.c_ext

    def steady_state(self) -> float:
        """Steady-state DV concentration C* [uM], solving dC/dt = 0.

        Transporter off -> analytic. Transporter on -> bracketed root find in
        (0, passive_steady_state), where dC/dt is +ve at 0 and -ve at the passive SS.
        """
        C_passive = self.passive_steady_state()
        if self.vmax <= 0:
            return C_passive
        lo, hi = 0.0, C_passive
        # guard the bracket: dCdt(lo) > 0 and dCdt(hi) < 0 by construction
        return float(brentq(self._dCdt, lo, hi, xtol=1e-9, rtol=1e-10))

    def accumulation_ratio(self) -> float:
        """Effective DV/source accumulation ratio at steady state, C* / C_ext."""
        return self.steady_state() / self.c_ext

    def simulate(self, t_end: float, C0: float = 0.0, n: int = 400):
        """Integrate C(t) from C0 over [0, t_end] (s). Returns (t, C).

        Stiff-capable (BDF): the approach to a very large accumulation ratio is a
        slow mode, so an implicit solver is used.
        """
        sol = solve_ivp(
            lambda t, y: [self._dCdt(y[0])],
            (0.0, t_end),
            [C0],
            method="BDF",
            t_eval=np.linspace(0.0, t_end, n),
            rtol=1e-8,
            atol=1e-6,
        )
        if not sol.success:
            raise RuntimeError(f"integration failed: {sol.message}")
        return sol.t, sol.y[0]

    def equilibration_time(self) -> float:
        """Characteristic time constant of the passive slow mode [s]: 1/(k_perm*f_dv)."""
        return 1.0 / (self.k_perm * self.f_dv)
