"""Synthetic validation of the Layer-C competitive-inhibition model (enzyme.py).

Machinery correctness BEFORE any real Ki/Km data (the #1b discipline): limits,
Cheng-Prusoff self-consistency, monotonicity, and Ki<->IC50 inversion.

    PYTHONPATH=src python3 tests/test_enzyme.py
    PYTHONPATH=src pytest tests/
"""

from __future__ import annotations

import math

from mdae import enzyme as E


def test_residual_activity_limits():
    ki, s, km = 10.0, 5.0, 2.0
    assert math.isclose(E.residual_activity(0.0, ki, s, km), 1.0, rel_tol=1e-12)  # no inhibitor
    assert E.residual_activity(1e9, ki, s, km) < 1e-6                             # saturating inhibitor
    # strictly decreasing in [I]
    assert E.residual_activity(1.0, ki, s, km) > E.residual_activity(50.0, ki, s, km)


def test_ic50_gives_half_activity():
    ki, s, km = 12.0, 8.0, 3.0
    ic50 = E.ic50_cheng_prusoff(ki, s, km)
    assert math.isclose(E.residual_activity(ic50, ki, s, km), 0.5, rel_tol=1e-12)


def test_cheng_prusoff_inversion():
    ki, s, km = 1.5, 100.0, 4.0
    ic50 = E.ic50_cheng_prusoff(ki, s, km)
    assert math.isclose(E.ki_from_ic50(ic50, s, km), ki, rel_tol=1e-12)
    # IC50 >= Ki always (competitive), equality only at [S]=0
    assert ic50 >= ki
    assert math.isclose(E.ic50_cheng_prusoff(ki, 0.0, km), ki, rel_tol=1e-12)


def test_velocity_consistency():
    ki, s, km, vmax = 20.0, 10.0, 5.0, 3.3
    # residual_activity is velocity(I)/velocity(0)
    a = E.residual_activity(7.0, ki, s, km)
    ratio = E.velocity(7.0, ki, s, km, vmax) / E.velocity(0.0, ki, s, km, vmax)
    assert math.isclose(a, ratio, rel_tol=1e-12)


def test_higher_ki_means_weaker_inhibition():
    # larger Ki (weaker binder) -> higher residual activity at the same [I]
    s, km, inh = 10.0, 4.0, 50.0
    assert E.residual_activity(inh, 5.0, s, km) < E.residual_activity(inh, 500.0, s, km)


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
