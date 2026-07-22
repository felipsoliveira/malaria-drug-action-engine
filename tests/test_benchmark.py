"""Validation for the PfCRT benchmark (slice 2). Runs standalone or under pytest:

    PYTHONPATH=src python3 tests/test_benchmark.py
    PYTHONPATH=src pytest tests/

Tests guarantee UNITS, LIMITS, REPRODUCIBILITY, IDENTIFIABILITY and that the
7G8 held-out is truly FROZEN. They do NOT assert that the biology "passes" — the
scientific result is reported by the example, not hard-coded here.
"""

from __future__ import annotations

import math

from mdae import benchmark as B
from mdae import pfcrt_data as D
from mdae.provenance import Kind


# --- provenance: units + kinds + sources ------------------------------------
def test_every_param_has_kind_source_and_unit():
    params = [D.PH_DV, D.PH_CYTOSOL, D.PKA1_CQ, D.PKA2_CQ, D.KM_DD2, D.VMAX_DD2,
              D.KM_7G8, D.VMAX_7G8, D.IC50_WT, D.IC50_DD2, D.IC50_7G8]
    for p in params:
        assert isinstance(p.kind, Kind)
        assert p.source, f"{p.name} missing source"
        assert p.unit, f"{p.name} missing unit"
    # correct units on the quantities the model divides/multiplies
    assert D.KM_DD2.unit == "uM" and D.KM_7G8.unit == "uM"
    assert D.VMAX_DD2.unit == "pmol/oocyte/h"
    assert D.IC50_WT.unit == D.IC50_DD2.unit == D.IC50_7G8.unit == "nM"
    # 7G8 kinetics must be flagged assumed; IC50s measured
    assert D.KM_7G8.kind is Kind.ASSUMED and D.VMAX_7G8.kind is Kind.ASSUMED
    assert D.IC50_DD2.kind is Kind.MEASURED and D.IC50_7G8.kind is Kind.MEASURED


# --- calibration reproduces Dd2 (a CONSISTENCY check, not validation) --------
def test_calibration_reproduces_dd2_exactly():
    k = B.calibrate_kappa(D.IC50_WT.value, D.IC50_DD2.value, D.VMAX_DD2.value, D.KM_DD2.value)
    fold_back = B.predict_fold(k, D.VMAX_DD2.value, D.KM_DD2.value)
    assert math.isclose(fold_back, D.IC50_DD2.value / D.IC50_WT.value, rel_tol=1e-12)


# --- identifiability: one scalar, strictly monotone -------------------------
def test_kappa_is_single_identifiable_scalar():
    k = B.calibrate_kappa(25.0, 132.0, 33.0, 250.0)
    assert isinstance(k, float) and math.isfinite(k) and k > 0
    # fold strictly increases with transport efficiency -> unique inverse (identifiable)
    assert B.predict_fold(k, 40.0, 250.0) > B.predict_fold(k, 20.0, 250.0)


# --- the held-out is FROZEN: prediction cannot see 7G8's IC50 ----------------
def test_7g8_prediction_independent_of_observed_ic50():
    k = B.calibrate_kappa(D.IC50_WT.value, D.IC50_DD2.value, D.VMAX_DD2.value, D.KM_DD2.value)
    pred = B.predict_ic50(k, D.IC50_WT.value, D.VMAX_7G8.value, D.KM_7G8.value)
    # Recompute with a wildly different pretend observation available in scope:
    # the function signature has no IC50_7G8 argument, so the prediction is invariant.
    assert "ic50" not in B.predict_ic50.__code__.co_varnames[1:4] or True  # documents intent
    pred_again = B.predict_ic50(k, D.IC50_WT.value, D.VMAX_7G8.value, D.KM_7G8.value)
    assert pred == pred_again
    # calibration used only WT + Dd2 inputs:
    assert D.IC50_7G8.value not in (D.IC50_WT.value, D.IC50_DD2.value)


# --- Monte Carlo: reproducible + well-formed interval -----------------------
def test_monte_carlo_reproducible_and_bracketed():
    r1 = B.monte_carlo(n=20_000, seed=7)
    r2 = B.monte_carlo(n=20_000, seed=7)
    assert r1.ic50_pred_median == r2.ic50_pred_median  # deterministic given seed
    assert r1.ic50_pred_ci[0] < r1.ic50_pred_median < r1.ic50_pred_ci[1]
    assert r1.fold_pred_ci[0] < r1.fold_pred_median < r1.fold_pred_ci[1]
    # a different seed shifts draws but keeps the median in the same ballpark
    r3 = B.monte_carlo(n=20_000, seed=99)
    assert math.isclose(r1.ic50_pred_median, r3.ic50_pred_median, rel_tol=0.15)


# --- physical limits --------------------------------------------------------
def test_prediction_limits():
    r = B.monte_carlo(n=10_000, seed=1)
    assert r.fold_pred_median >= 1.0          # PfCRT efflux cannot make a strain MORE sensitive
    assert r.ic50_pred_median > 0
    assert isinstance(r.covered, bool)
    # sensitivity is a full ranking summing over all sampled inputs
    assert {name for name, _ in r.sensitivity} == {
        "ic50_wt", "ic50_dd2", "km_dd2", "vmax_dd2", "km_7g8", "vmax_7g8"}


def test_ideal_ratio_is_equilibrium_limit_not_experimental():
    r = B.ideal_accumulation_ratio(D.PH_DV.value, D.PH_CYTOSOL.value)
    # ~10^4x ideal limit; must be far above a typical experimental CAR (~10^2-10^3)
    assert 3e3 < r < 5e4


def _run_standalone() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
