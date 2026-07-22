"""Validation for slice 1. Runs under pytest OR standalone:

    PYTHONPATH=src python3 tests/test_slice1.py     # standalone, no pytest needed
    PYTHONPATH=src pytest tests/                      # if pytest installed
"""

from __future__ import annotations

import math

from mdae import WeakBase, IonTrapModel, Compartment, accumulation_ratio
from mdae import parameters as P


# --- speciation limits ------------------------------------------------------
def test_neutral_fraction_limits():
    cq = P.CHLOROQUINE
    # Very high pH -> essentially all neutral.
    assert cq.neutral_fraction(14.0) > 0.99
    # Very low pH -> essentially all charged (neutral fraction -> 0).
    assert cq.neutral_fraction(2.0) < 1e-6
    # Fractions are a partition of unity.
    assert math.isclose(cq.neutral_fraction(7.0) + cq.charged_fraction(7.0), 1.0, rel_tol=1e-12)


def test_monoprotic_matches_hasselbalch():
    # For a monoprotic base, R = (1+10^(pKa-pH_in))/(1+10^(pKa-pH_out)).
    b = WeakBase("mono", pKa=(9.0,))
    pH_in, pH_out = 5.0, 7.4
    expected = (1 + 10 ** (9.0 - pH_in)) / (1 + 10 ** (9.0 - pH_out))
    assert math.isclose(accumulation_ratio(b, pH_in, pH_out), expected, rel_tol=1e-10)


# --- the analytical oracle for chloroquine ----------------------------------
def test_chloroquine_accumulation_ratio():
    R = accumulation_ratio(P.CHLOROQUINE, P.DIGESTIVE_VACUOLE.pH, P.CYTOSOL_PLASMA.pH)
    # Diprotic weak base into DV -> ~10^4; classic result is 10^3-10^4.
    assert 1e4 < R < 5e4, f"R={R} outside expected 10^4 band"


# --- ODE model must reproduce the independent oracle (transporter off) ------
def test_ode_passive_matches_oracle():
    sens = IonTrapModel(
        base=P.CHLOROQUINE, source=P.CYTOSOL_PLASMA, vacuole=P.DIGESTIVE_VACUOLE,
        c_ext=P.C_EXT_UM, k_perm=P.K_PERM_PER_S, **P.PFCRT_SENSITIVE,
    )
    R_oracle = accumulation_ratio(P.CHLOROQUINE, P.DIGESTIVE_VACUOLE.pH, P.CYTOSOL_PLASMA.pH)
    R_model = sens.accumulation_ratio()
    # Two independent derivations (speciation vs flux balance) must agree.
    assert math.isclose(R_model, R_oracle, rel_tol=1e-9)


def test_simulation_converges_to_steady_state():
    sens = IonTrapModel(
        base=P.CHLOROQUINE, source=P.CYTOSOL_PLASMA, vacuole=P.DIGESTIVE_VACUOLE,
        c_ext=P.C_EXT_UM, k_perm=P.K_PERM_PER_S, **P.PFCRT_SENSITIVE,
    )
    t, C = sens.simulate(t_end=8 * sens.equilibration_time(), C0=0.0)
    assert math.isclose(C[-1], sens.steady_state(), rel_tol=1e-3)


# --- resistance is an EXPOSURE change, not a binding change ------------------
def test_pfcrt_lowers_exposure_without_touching_binding():
    common = dict(base=P.CHLOROQUINE, source=P.CYTOSOL_PLASMA, vacuole=P.DIGESTIVE_VACUOLE,
                  c_ext=P.C_EXT_UM, k_perm=P.K_PERM_PER_S)
    sens = IonTrapModel(**common, **P.PFCRT_SENSITIVE)
    res = IonTrapModel(**common, **P.PFCRT_RESISTANT)
    c_sens, c_res = sens.steady_state(), res.steady_state()
    # Exposure drops in the resistant strain...
    assert c_res < c_sens
    # ...while the drug (its pKa) and thus the passive accumulation ratio are unchanged:
    assert sens.f_dv == res.f_dv and sens.f_ext == res.f_ext
    # ...and the resistant strain falls below the functional threshold while sensitive stays above.
    assert c_sens > P.C_FUNCTIONAL_THRESHOLD_UM > c_res


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
