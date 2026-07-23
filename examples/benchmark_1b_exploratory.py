"""Benchmark #1b — EXPLORATORY analysis ONLY (the confirmatory verdict is NOT EVALUABLE).

The pre-registered confirmatory benchmark (prereg-1b) is **not evaluable**: the
Table S2 resistance indices are a heterogeneous compilation (per-study CQS
reference, aggregated across labs/assays), which fails the comparability audit
(stop rule). See docs/benchmark_1b_result.md.

This script runs the frozen model over the heterogeneous RIs anyway, as a clearly
labelled EXPLORATORY analysis. It assigns NO pass/fail and NO adequacy verdict.
Read the numbers as hypothesis-generating only.

Run:  PYTHONPATH=src python3 examples/benchmark_1b_exploratory.py
"""

from __future__ import annotations

from mdae import benchmark_1b as B
from mdae import pfcrt_data as D


def _kin(name, km, vmax):
    return B.StrainKinetics(name, km, vmax)


def main() -> None:
    print("=" * 74)
    print("Benchmark #1b — EXPLORATORY ONLY (confirmatory outcome = NOT EVALUABLE)")
    print("Phenotypes (Table S2 RIs) are a heterogeneous compilation -> no PASS/FAIL.")
    print("=" * 74)

    dd2 = _kin("Dd2", D.KM_DD2, D.VMAX_DD2)
    heldouts = [
        (_kin("K1", D.KM_K1, D.VMAX_K1), D.RI_K1),
        (_kin("GB4", D.KM_GB4, D.VMAX_GB4), D.RI_GB4),
        (_kin("Ecu1110", D.KM_ECU, D.VMAX_ECU), D.RI_ECU),
        (_kin("7G8", D.KM_7G8, D.VMAX_7G8), D.RI_7G8),
    ]
    print(f"\nCalibration datum: RI_Dd2 = {D.RI_DD2.value} +/- {D.RI_DD2.sd}  (Table S2, n=8; heterogeneous)")

    for cstar in (0.0, 200.0, 400.0, 600.0):
        tag = " (= linear #1a limit)" if cstar == 0 else ""
        print(f"\n--- C* = {cstar:.0f} uM{tag} ---")
        print(f"  {'strain':8s} {'pred RI (90% CI)':>22s} {'obs RI':>10s} {'fold err':>9s}")
        for kin, obs in heldouts:
            p = B.predict_mc(D.RI_DD2, dd2, kin, cstar, n=40_000, seed=12345)
            fe = max(p.median / obs.value, obs.value / p.median)
            ci = f"{p.median:5.2f} ({p.ci[0]:4.2f}-{p.ci[1]:4.2f})"
            print(f"  {kin.name:8s} {ci:>22s} {obs.value:10.2f} {fe:8.2f}x")

    print("\n" + "-" * 74)
    print("EXPLORATORY reading (NOT a benchmark verdict):")
    print("  Even calibrated on the Table-S2 RI, the transport model tracks K1/GB4 but")
    print("  under-predicts 7G8 (and saturation, larger C*, worsens 7G8 further). But NONE")
    print("  of this counts as pass/fail: the RIs are not mutually comparable, so any")
    print("  agreement or disagreement could be an artifact of pooling heterogeneous assays.")
    print("  Confirmatory verdict stays: NOT EVALUABLE (phenotype heterogeneity).")


if __name__ == "__main__":
    main()
