"""Benchmark #1b — nonlinear PfCRT transport -> resistance index (RI). PRE-REGISTERED.

Implements exactly the frozen protocol in docs/benchmark_1b_preregistration.md
(git tag `prereg-1b`). Validated on SYNTHETIC data first; the real Table S2
phenotypes are plugged in unchanged and the benchmark is run once.

Frozen model:

    RI_s(C*) = 1 + (RI_Dd2 - 1) * [ Vmax_s /(Km_s + C*) ] / [ Vmax_Dd2 /(Km_Dd2 + C*) ]

- One calibration datum: RI_Dd2. C* -> 0 recovers the linear #1a model.
- C* is a SCENARIO variable (swept, not sampled).
- Monte Carlo runs over MEASURED parameters only: RI_Dd2 and the Km/Vmax of Dd2
  and the target strain (these are transport INPUTS). The held-out strain's own
  RI -- and its SEM -- NEVER enter the prediction or its interval. The held-out
  SEM characterizes the OBSERVATION (used only in the comparison), it does not
  widen the predicted interval.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .provenance import Param


@dataclass(frozen=True)
class StrainKinetics:
    """PfCRT CQ-transport kinetics for one variant (transport input, not phenotype)."""

    name: str
    km: Param    # uM
    vmax: Param  # pmol/oocyte/h


# --- frozen point model -----------------------------------------------------
def transport_term(vmax: float, km: float, cstar: float) -> float:
    """Saturable transport efficiency Vmax/(Km + C*)."""
    return vmax / (km + cstar)


def predict_ri(ri_dd2, km_dd2, vmax_dd2, km_s, vmax_s, cstar):
    """Predicted resistance index of a strain. Takes NO held-out RI (frozen)."""
    eff_s = transport_term(vmax_s, km_s, cstar)
    eff_dd2 = transport_term(vmax_dd2, km_dd2, cstar)
    return 1.0 + (ri_dd2 - 1.0) * (eff_s / eff_dd2)


# --- uncertainty (measured inputs only) -------------------------------------
@dataclass
class Prediction:
    strain: str
    cstar: float
    median: float
    ci: tuple[float, float]  # 90% (5th, 95th percentile)


def predict_mc(ri_dd2: Param, dd2: StrainKinetics, s: StrainKinetics, cstar: float,
               n: int = 50_000, seed: int = 12345) -> Prediction:
    """MC over the measured inputs (RI_Dd2, Dd2 & strain Km/Vmax). Held-out RI absent."""
    rng = np.random.default_rng(seed)
    ri = ri_dd2.sample(rng, n)
    kmd, vmd = dd2.km.sample(rng, n), dd2.vmax.sample(rng, n)
    kms, vms = s.km.sample(rng, n), s.vmax.sample(rng, n)
    pred = predict_ri(ri, kmd, vmd, kms, vms, cstar)
    return Prediction(s.name, cstar, float(np.median(pred)),
                      (float(np.percentile(pred, 5)), float(np.percentile(pred, 95))))


# --- frozen pass/fail evaluation --------------------------------------------
@dataclass
class StrainVerdict:
    strain: str
    cstar: float
    pred_median: float
    pred_ci: tuple[float, float]
    obs: float
    obs_ci: tuple[float, float]
    inside_ci: bool
    fold_error: float
    width_ratio: float
    passed: bool


def evaluate(pred: Prediction, obs_ri: Param) -> StrainVerdict:
    """Apply the frozen per-strain criteria. obs_ri is used ONLY here, not in prediction."""
    lo, hi = pred.ci
    obs = obs_ri.value
    inside = bool(lo <= obs <= hi)
    fold_error = max(pred.median / obs, obs / pred.median)
    width_ratio = hi / lo
    passed = inside and fold_error <= 1.5 and width_ratio <= 2.0
    return StrainVerdict(pred.strain, pred.cstar, pred.median, pred.ci, obs,
                         obs_ri.interval(), inside, fold_error, width_ratio, passed)


@dataclass
class BenchmarkVerdict:
    cstar: float
    strains: list[StrainVerdict]
    n_pass: int
    n_total: int
    max_fold_error: float
    adequate: bool
    rmse_log2: float
    geomean_fold_error: float


def run(ri_dd2: Param, dd2: StrainKinetics,
        heldouts: list[tuple[StrainKinetics, Param]], cstar: float,
        n: int = 50_000, seed: int = 12345) -> BenchmarkVerdict:
    """Predict & evaluate all held-out strains at one C* scenario (frozen criteria)."""
    verdicts = [evaluate(predict_mc(ri_dd2, dd2, kin, cstar, n, seed), obs)
                for kin, obs in heldouts]
    n_pass = sum(v.passed for v in verdicts)
    max_fe = max(v.fold_error for v in verdicts)
    adequate = (n_pass >= 3) and (max_fe <= 2.0)
    rmse = math.sqrt(sum((math.log2(v.pred_median) - math.log2(v.obs)) ** 2
                         for v in verdicts) / len(verdicts))
    geo = math.exp(sum(math.log(v.fold_error) for v in verdicts) / len(verdicts))
    return BenchmarkVerdict(cstar, verdicts, n_pass, len(verdicts), max_fe, adequate, rmse, geo)


def sweep(ri_dd2: Param, dd2: StrainKinetics,
          heldouts: list[tuple[StrainKinetics, Param]],
          cstars=(0.0, 200.0, 400.0, 600.0), n: int = 50_000, seed: int = 12345
          ) -> list[BenchmarkVerdict]:
    """Run the frozen benchmark across the C* scenario sweep (C*=0 is the #1a limit)."""
    return [run(ri_dd2, dd2, heldouts, c, n, seed) for c in cstars]
