"""Layer C — molecular function for an enzyme target (competitive inhibition).

For a competitive inhibitor of a Michaelis-Menten enzyme (e.g. pyrimethamine
competing with dihydrofolate at PfDHFR), the steady-state velocity is

    v(I) = Vmax * [S] / ( Km * (1 + [I]/Ki) + [S] )

so the fractional residual activity relative to no inhibitor is

    a(I) = v(I)/v(0) = (Km + [S]) / ( Km * (1 + [I]/Ki) + [S] )

and the half-inhibition concentration follows the Cheng-Prusoff relation

    IC50 = Ki * (1 + [S]/Km).

This is the Layer-C "binding -> function" conversion: given an affinity (Ki) and the
assay context ([S], Km), predict the functional output (residual activity, IC50). It
is target-class specific (competitive enzyme inhibition) and makes NO cellular-
exposure assumption -- that is Layer B. All concentrations share one unit; Ki, IC50,
[S], Km in the same unit (e.g. nM or uM).
"""

from __future__ import annotations


def residual_activity(inhibitor: float, ki: float, substrate: float, km: float) -> float:
    """v(I)/v(0) for a competitive inhibitor (dimensionless, in (0, 1])."""
    return (km + substrate) / (km * (1.0 + inhibitor / ki) + substrate)


def velocity(inhibitor: float, ki: float, substrate: float, km: float, vmax: float) -> float:
    """Absolute steady-state velocity v(I) for a competitive inhibitor."""
    return vmax * substrate / (km * (1.0 + inhibitor / ki) + substrate)


def ic50_cheng_prusoff(ki: float, substrate: float, km: float) -> float:
    """IC50 of a competitive inhibitor: Ki * (1 + [S]/Km)."""
    return ki * (1.0 + substrate / km)


def ki_from_ic50(ic50: float, substrate: float, km: float) -> float:
    """Invert Cheng-Prusoff to recover Ki from a measured IC50 and the assay [S], Km."""
    return ic50 / (1.0 + substrate / km)
