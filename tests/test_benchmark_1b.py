"""Synthetic validation of Benchmark #1b BEFORE any real Table S2 data is used.

Proves (pre-registration steps 1-3): the C*=0 limit reproduces the #1a linear model;
no held-out RI (or its SEM) can enter the prediction or its interval; the pipeline
recovers a synthetic ground truth; and the frozen pass/fail logic is correct.

    PYTHONPATH=src python3 tests/test_benchmark_1b.py
    PYTHONPATH=src pytest tests/
"""

from __future__ import annotations

import inspect
import math

from mdae import benchmark_1b as B1b
from mdae import benchmark as B1a
from mdae.provenance import Param, Kind


def _p(name, value, sd=0.0, kind=Kind.MEASURED):
    return Param(name, value, "u", kind, sd=sd)


# --- step 2: C*=0 reproduces the linear #1a model ---------------------------
def test_cstar0_recovers_1a_linear_model():
    ri_dd2, kmd, vmd, kms, vms = 5.28, 232.0, 61.0, 117.0, 9.0
    got = B1b.predict_ri(ri_dd2, kmd, vmd, kms, vms, cstar=0.0)
    # #1a form: 1 + (RI-1) * (Vmax_s/Km_s) / (Vmax_Dd2/Km_Dd2)
    expected = 1 + (ri_dd2 - 1) * (vms / kms) / (vmd / kmd)
    assert math.isclose(got, expected, rel_tol=1e-12)
    # and it matches the #1a kappa formulation exactly
    kappa = B1a.calibrate_kappa(1.0, ri_dd2, vmd, kmd)  # ic50_wt=1 -> fold==RI
    assert math.isclose(got, B1a.predict_fold(kappa, vms, kms), rel_tol=1e-12)


def test_cstar_raises_denominator_and_changes_ratio():
    # sanity: at large C*, the term -> Vmax/C*, so the ratio -> Vmax_s/Vmax_Dd2
    ri = B1b.predict_ri(5.28, 232.0, 61.0, 117.0, 9.0, cstar=1e9)
    expected = 1 + (5.28 - 1) * (9.0 / 61.0)
    assert math.isclose(ri, expected, rel_tol=1e-6)


# --- step 3: the held-out RI (and its SEM) cannot leak into the prediction --
def test_predict_has_no_heldout_argument():
    params = list(inspect.signature(B1b.predict_ri).parameters)
    assert params == ["ri_dd2", "km_dd2", "vmax_dd2", "km_s", "vmax_s", "cstar"]
    mc_params = list(inspect.signature(B1b.predict_mc).parameters)
    assert "obs" not in mc_params and "ri_s" not in mc_params and "obs_ri" not in mc_params


def test_heldout_value_and_sem_do_not_affect_prediction():
    ri_dd2 = _p("ri_dd2", 5.28, sd=0.3)
    dd2 = B1b.StrainKinetics("Dd2", _p("km", 232, 11), _p("vmax", 61, 6))
    s = B1b.StrainKinetics("S", _p("km", 117, 6), _p("vmax", 9, 1))
    pred = B1b.predict_mc(ri_dd2, dd2, s, cstar=200.0, n=20_000, seed=3)
    # two DIFFERENT held-out observations, with very different SEMs:
    v_tight = B1b.evaluate(pred, _p("obs", 3.0, sd=0.01))
    v_wide = B1b.evaluate(pred, _p("obs", 9.0, sd=100.0))
    # prediction (median + CI) is byte-identical regardless of the observation:
    assert v_tight.pred_median == v_wide.pred_median
    assert v_tight.pred_ci == v_wide.pred_ci
    # the predicted CI width does NOT include the observation's SEM
    assert v_tight.pred_ci == pred.ci


# --- the machinery is correct: recover a synthetic ground truth --------------
def test_synthetic_ground_truth_recovery():
    cstar_true, ri_dd2_true = 350.0, 4.0
    dd2 = B1b.StrainKinetics("Dd2", _p("km", 232), _p("vmax", 61))  # zero-SEM -> deterministic
    strains = [B1b.StrainKinetics(nm, _p("km", km), _p("vmax", vm))
               for nm, km, vm in [("A", 117, 9), ("B", 293, 79), ("C", 191, 36)]]
    for s in strains:
        truth = B1b.predict_ri(ri_dd2_true, 232, 61, s.km.value, s.vmax.value, cstar_true)
        pred = B1b.predict_mc(_p("ri", ri_dd2_true), dd2, s, cstar_true, n=2000, seed=1)
        assert math.isclose(pred.median, truth, rel_tol=1e-9)


# --- frozen pass/fail logic --------------------------------------------------
def _pred(median, lo, hi):
    return B1b.Prediction("S", 0.0, median, (lo, hi))


def test_pass_requires_all_three_conditions():
    # inside CI, small error, narrow interval -> PASS
    assert B1b.evaluate(_pred(3.0, 2.2, 4.0), _p("o", 3.2)).passed
    # inside CI, small error, but GIANT interval (width>2x) -> FAIL
    assert not B1b.evaluate(_pred(3.0, 1.0, 3.0), _p("o", 3.0)).passed
    # inside CI, narrow, but fold error > 1.5x -> FAIL
    assert not B1b.evaluate(_pred(3.0, 2.5, 4.0), _p("o", 4.7)).passed
    # outside CI -> FAIL
    assert not B1b.evaluate(_pred(3.0, 2.5, 3.6), _p("o", 5.0)).passed


def test_adequacy_and_diagnostics():
    ri_dd2 = _p("ri", 4.0)
    dd2 = B1b.StrainKinetics("Dd2", _p("km", 232), _p("vmax", 61))
    # synthetic held-outs whose observed RI equals the model prediction -> all pass
    heldouts = []
    for nm, km, vm in [("A", 117, 9), ("B", 293, 79), ("C", 191, 36), ("D", 275, 57)]:
        s = B1b.StrainKinetics(nm, _p("km", km), _p("vmax", vm))
        truth = B1b.predict_ri(4.0, 232, 61, km, vm, 300.0)
        heldouts.append((s, _p(f"obs_{nm}", truth)))
    v = B1b.run(ri_dd2, dd2, heldouts, cstar=300.0, n=2000, seed=2)
    assert v.n_pass == 4 and v.adequate
    assert v.rmse_log2 < 1e-6 and math.isclose(v.geomean_fold_error, 1.0, rel_tol=1e-6)


def test_reproducible():
    ri_dd2 = _p("ri", 5.0, sd=0.2)
    dd2 = B1b.StrainKinetics("Dd2", _p("km", 232, 11), _p("vmax", 61, 6))
    s = B1b.StrainKinetics("S", _p("km", 117, 6), _p("vmax", 9, 1))
    a = B1b.predict_mc(ri_dd2, dd2, s, 200.0, n=10_000, seed=7)
    b = B1b.predict_mc(ri_dd2, dd2, s, 200.0, n=10_000, seed=7)
    assert a.median == b.median and a.ci == b.ci


def _run_standalone() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t(); print(f"PASS {t.__name__}")
        except AssertionError as e:
            failed += 1; print(f"FAIL {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1; print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
