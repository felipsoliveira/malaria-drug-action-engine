"""Falsifiable benchmark: predict 7G8 chloroquine resistance from PfCRT transport,
calibrated ONLY on Dd2, with 7G8's phenotype (IC50) held out.

Model
-----
In the first-order (specificity-constant) regime of the ion-trapping ODE in
``exposure.py`` — valid when the digestive-vacuole free concentration is well
below the transporter Km — the steady state gives a fold-resistance

    fold(strain) = IC50(strain) / IC50(WT) = 1 + kappa * (Vmax / Km)_strain

Derivation: with first-order PfCRT efflux k_eff = Vmax_cell/Km of the charged
species, steady state of ``exposure.IonTrapModel`` yields
``fold = 1 + (k_eff * cf_dv) / (k_perm * f_dv)``. Writing Vmax_cell = alpha *
Vmax_oocyte collapses every unknown (alpha oocyte->cell scaling, passive
permeability k_perm, and the pH/pKa speciation cf_dv/f_dv) into the SINGLE lumped
constant ``kappa = alpha * cf_dv / (k_perm * f_dv)``. All of it cancels from the
*fold ratio* except kappa. (Km and Vmax are the measured per-strain inputs;
Vmax/Km is the transporter's specificity constant.)

Identifiability
---------------
``kappa`` is one scalar fixed by one datum (Dd2's fold-resistance). It is
**structurally identifiable under this assumed linear model** (one parameter, one
constraint; fold is strictly monotone in efficiency, so the inverse is unique).
This is NOT a general practical-identifiability claim: it holds only conditional
on the model form fold = 1 + kappa*(Vmax/Km) and on the linear-regime assumption.
We deliberately do NOT fit several free parameters to a single IC50.

Frozen held-out
---------------
``kappa`` is calibrated from Dd2 (and the WT reference) ONLY. 7G8's IC50 is never
used to fit it. 7G8's *transport kinetics* (Km, Vmax) enter as input — measured in
the SAME experiment as Dd2 (Summers 2014 Table 1), so Vmax/Km is comparable across
strains. Per the scientific contract, the held-out target is therefore 7G8's
PHENOTYPE (IC50 / fold). ``predict_ic50`` takes no observed 7G8 IC50 argument.

Not this benchmark
------------------
The ideal ion-trapping accumulation ratio (~10^4x, see :func:`ideal_accumulation_ratio`)
is the model's EQUILIBRIUM LIMIT from pH/pKa. It is not an experimental cellular
accumulation and is not what this benchmark tests.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .speciation import WeakBase
from . import pfcrt_data as D


# --- deterministic point model ----------------------------------------------
def transport_efficiency(vmax: float, km: float) -> float:
    """Specificity constant Vmax/Km (per-strain transport efficiency)."""
    return vmax / km


def calibrate_kappa(ic50_wt: float, ic50_calib: float, vmax_calib: float, km_calib: float) -> float:
    """The single fitted constant, from the calibration strain (Dd2) alone.

    kappa = (fold_obs_calib - 1) / (Vmax/Km)_calib ,  fold_obs_calib = IC50_calib/IC50_wt.
    """
    fold_obs = ic50_calib / ic50_wt
    return (fold_obs - 1.0) / transport_efficiency(vmax_calib, km_calib)


def predict_fold(kappa: float, vmax: float, km: float) -> float:
    """Predicted fold-resistance for a strain from its transport kinetics. No IC50 input."""
    return 1.0 + kappa * transport_efficiency(vmax, km)


def predict_ic50(kappa: float, ic50_wt: float, vmax: float, km: float) -> float:
    """Predicted absolute IC50 [nM]. The held-out observation is NOT an argument."""
    return predict_fold(kappa, vmax, km) * ic50_wt


def ideal_accumulation_ratio(pH_dv: float, pH_source: float) -> float:
    """Ideal equilibrium ion-trapping ratio Z(pH_dv)/Z(pH_source) for chloroquine.

    This is the model's thermodynamic LIMIT (transporter off), not experimental
    accumulation. With the parasite cytosol (~7.2) as source it is ~10^4x; with
    plasma (7.4) ~2x10^4x. Measured cellular accumulation ratios are far lower
    (hundreds), because of binding, efflux and non-equilibrium — do not conflate.
    """
    cq = WeakBase("chloroquine", pKa=(D.PKA1_CQ.value, D.PKA2_CQ.value))
    return float(cq.partition(pH_dv) / cq.partition(pH_source))


# --- uncertainty propagation & sensitivity ----------------------------------
_INPUTS = ("ic50_wt", "ic50_dd2", "km_dd2", "vmax_dd2", "km_7g8", "vmax_7g8")


@dataclass
class BenchmarkResult:
    kappa_median: float
    kappa_ci: tuple[float, float]
    fold_pred_median: float
    fold_pred_ci: tuple[float, float]
    ic50_pred_median: float
    ic50_pred_ci: tuple[float, float]         # 90% CI (5th, 95th pct)
    ic50_obs: float
    ic50_obs_ci: tuple[float, float]
    fold_error: float                         # predicted median / observed
    covered: bool                             # observed IC50 inside predicted 90% CI?
    sensitivity: list[tuple[str, float]]      # (input, |Spearman rho| vs predicted IC50), ranked


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    xr = np.argsort(np.argsort(x))
    yr = np.argsort(np.argsort(y))
    xr = xr - xr.mean()
    yr = yr - yr.mean()
    denom = np.sqrt((xr**2).sum() * (yr**2).sum())
    return float((xr * yr).sum() / denom) if denom else 0.0


def monte_carlo(n: int = 50_000, seed: int = 12345) -> BenchmarkResult:
    """Propagate parameter uncertainty; return the frozen 7G8 prediction + sensitivity.

    Reproducible: identical (n, seed) -> identical result.
    """
    rng = np.random.default_rng(seed)
    s = {
        "ic50_wt": D.IC50_WT.sample(rng, n),
        "ic50_dd2": D.IC50_DD2.sample(rng, n),
        "km_dd2": D.KM_DD2.sample(rng, n),
        "vmax_dd2": D.VMAX_DD2.sample(rng, n),
        "km_7g8": D.KM_7G8.sample(rng, n),
        "vmax_7g8": D.VMAX_7G8.sample(rng, n),
    }
    ic50_7g8_obs = D.IC50_7G8.sample(rng, n)

    kappa = (s["ic50_dd2"] / s["ic50_wt"] - 1.0) / (s["vmax_dd2"] / s["km_dd2"])
    fold_pred = 1.0 + kappa * (s["vmax_7g8"] / s["km_7g8"])
    ic50_pred = fold_pred * s["ic50_wt"]  # frozen: does not use ic50_7g8_obs

    def ci(a):
        return (float(np.percentile(a, 5)), float(np.percentile(a, 95)))

    ic50_pred_ci = ci(ic50_pred)
    obs_lo, obs_hi = D.IC50_7G8.interval()

    sens = sorted(
        ((k, abs(_spearman(v, ic50_pred))) for k, v in s.items()),
        key=lambda kv: kv[1], reverse=True,
    )

    return BenchmarkResult(
        kappa_median=float(np.median(kappa)),
        kappa_ci=ci(kappa),
        fold_pred_median=float(np.median(fold_pred)),
        fold_pred_ci=ci(fold_pred),
        ic50_pred_median=float(np.median(ic50_pred)),
        ic50_pred_ci=ic50_pred_ci,
        ic50_obs=float(D.IC50_7G8.value),
        ic50_obs_ci=(obs_lo, obs_hi),
        fold_error=float(np.median(ic50_pred) / D.IC50_7G8.value),
        covered=(ic50_pred_ci[0] <= D.IC50_7G8.value <= ic50_pred_ci[1]),
        sensitivity=sens,
    )
