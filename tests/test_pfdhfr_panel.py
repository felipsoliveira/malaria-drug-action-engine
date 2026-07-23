"""PfDHFR panel: data sanity + ROUND-TRIP self-consistency.

These are UNITS / IMPLEMENTATION / PROVENANCE checks on the Sirawaraporn 1997 panel.
They are explicitly NOT a biological validation of binding->function: the Table-3 Ki
were derived from IC50 via Cheng-Prusoff, so the round-trip is an algebraic identity.

    PYTHONPATH=src python3 tests/test_pfdhfr_panel.py
"""

from __future__ import annotations

import math

from mdae import enzyme as E
from mdae import pfdhfr_data as D
from mdae.provenance import Kind


def test_panel_shape_units_and_provenance():
    assert len(D.PANEL) == 7
    for v in D.PANEL:
        assert v.kcat.unit == "1/s"
        assert v.km_dhf.unit == "uM"
        assert v.ki_pyr.unit == "nM"
        assert v.km_dhf.kind is Kind.MEASURED
        # Ki are DERIVED (from IC50), not directly measured -- must be tagged so.
        assert v.ki_pyr.kind is Kind.DERIVED
        assert "Sirawaraporn" in v.ki_pyr.source


def test_ki_increases_along_resistance_series():
    ki = [v.ki_pyr.value for v in D.PANEL]
    assert ki == sorted(ki)                    # 1.5 -> 859 nM, monotonic
    assert ki[0] == 1.5 and ki[-1] == 859.0


def test_round_trip_is_self_consistent_only():
    # ki -> IC50 -> ki reproduces the input exactly (pure algebra); IC50 > Ki since [S]>0.
    S = D.ASSAY_DHF_UM
    for v in D.PANEL:
        ki, km = v.ki_pyr.value, v.km_dhf.value
        ic50 = E.ic50_cheng_prusoff(ki, S, km)
        assert ic50 > ki
        assert math.isclose(E.ki_from_ic50(ic50, S, km), ki, rel_tol=1e-12)


def test_km_dependence_reorders_relative_to_ki():
    # The functional IC50 is NOT a pure function of Ki: low-Km variants get a larger
    # Cheng-Prusoff factor. Concretely, N51I+S108N (Ki=37, Km=6) has a HIGHER IC50 than
    # C59R+S108N (Ki=72, Km=24) despite a lower Ki -- the Km-dependence Layer C encodes.
    S = D.ASSAY_DHF_UM
    a = D.BY_NAME["N51I+S108N"]
    b = D.BY_NAME["C59R+S108N"]
    assert a.ki_pyr.value < b.ki_pyr.value
    assert E.ic50_cheng_prusoff(a.ki_pyr.value, S, a.km_dhf.value) > \
           E.ic50_cheng_prusoff(b.ki_pyr.value, S, b.km_dhf.value)


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
