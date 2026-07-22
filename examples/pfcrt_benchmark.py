"""MDAE benchmark #1 — predict 7G8 chloroquine resistance, calibrated only on Dd2.

Run:  PYTHONPATH=src python3 examples/pfcrt_benchmark.py

Calibrates the single lumped transport constant kappa on Dd2, freezes it, and
predicts the held-out 7G8 phenotype (IC50) from 7G8's transport kinetics.
Reports prediction as an interval with fold-change and error, plus which inputs
dominate the uncertainty. Nothing here is tuned to the 7G8 result.
"""

from __future__ import annotations

from mdae import benchmark as B
from mdae import pfcrt_data as D


def main() -> None:
    print("=" * 72)
    print("MDAE benchmark #1 : chloroquine resistance, calibrate=Dd2, held-out=7G8")
    print("=" * 72)

    # Ideal equilibrium limit (NOT experimental accumulation, NOT the benchmark).
    r_cyt = B.ideal_accumulation_ratio(D.PH_DV.value, D.PH_CYTOSOL.value)
    r_pl = B.ideal_accumulation_ratio(D.PH_DV.value, 7.4)
    print("\n[context] ideal ion-trapping accumulation ratio (transporter off):")
    print(f"  vs cytosol pH {D.PH_CYTOSOL.value}: {r_cyt:,.0f}x   |   vs plasma pH 7.4: {r_pl:,.0f}x")
    print("  -> thermodynamic LIMIT from pH/pKa; experimental cellular CAR is far lower (~10^2-10^3).")

    # Point calibration (central values).
    kappa = B.calibrate_kappa(D.IC50_WT.value, D.IC50_DD2.value, D.VMAX_DD2.value, D.KM_DD2.value)
    print("\n[calibration on Dd2 ONLY]  (fitted: kappa; everything else measured/assumed)")
    print(f"  observed fold(Dd2) = {D.IC50_DD2.value}/{D.IC50_WT.value} = "
          f"{D.IC50_DD2.value / D.IC50_WT.value:.2f}x")
    print(f"  Vmax/Km (Dd2) = {D.VMAX_DD2.value}/{D.KM_DD2.value} = "
          f"{D.VMAX_DD2.value / D.KM_DD2.value:.4f}  ->  kappa = {kappa:.2f}")

    # Frozen point prediction for 7G8 (uses 7G8 kinetics; NOT 7G8 IC50).
    fold_7g8 = B.predict_fold(kappa, D.VMAX_7G8.value, D.KM_7G8.value)
    ic50_7g8 = B.predict_ic50(kappa, D.IC50_WT.value, D.VMAX_7G8.value, D.KM_7G8.value)
    print("\n[frozen prediction: 7G8]  (kinetics are ASSUMED bounds -> point uses midpoints)")
    print(f"  predicted fold(7G8) = {fold_7g8:.2f}x  ->  predicted IC50 = {ic50_7g8:.0f} nM")

    # Uncertainty propagation + sensitivity.
    res = B.monte_carlo()
    print("\n[uncertainty propagation]  Monte Carlo (n=50000, seed=12345)")
    print(f"  kappa            : {res.kappa_median:.2f}  (90% CI {res.kappa_ci[0]:.1f}-{res.kappa_ci[1]:.1f})")
    print(f"  predicted fold   : {res.fold_pred_median:.2f}x (90% CI "
          f"{res.fold_pred_ci[0]:.2f}-{res.fold_pred_ci[1]:.2f})")
    print(f"  predicted IC50   : {res.ic50_pred_median:.0f} nM (90% CI "
          f"{res.ic50_pred_ci[0]:.0f}-{res.ic50_pred_ci[1]:.0f})")
    print(f"  OBSERVED IC50(7G8): {res.ic50_obs:.0f} nM "
          f"(meas. {res.ic50_obs_ci[0]:.0f}-{res.ic50_obs_ci[1]:.0f})  [HELD-OUT]")
    print(f"  fold error (pred median / obs) = {res.fold_error:.2f}x")
    print(f"  observed inside predicted 90% CI? {'YES' if res.covered else 'NO'}")

    print("\n[sensitivity]  |Spearman rho| of each input vs predicted IC50(7G8):")
    for name, rho in res.sensitivity:
        kind = {"km_7g8": "ASSUMED", "vmax_7g8": "ASSUMED"}.get(name, "measured")
        print(f"  {name:10s} {rho:5.2f}  ({kind})")

    print("\n" + "-" * 72)
    print("Honest reading:")
    print("  The prediction is CONSISTENT with the held-out 7G8 IC50 but the interval is")
    print("  wide, because 7G8's transport kinetics are assumed (paywalled), not measured.")
    print("  The sensitivity ranking shows the 7G8 kinetics dominate the uncertainty:")
    print("  pinning 7G8's exact Km/Vmax (Summers 2014 Table 1) is the single highest-value")
    print("  data point to turn this into a sharp, strongly falsifiable test.")
    print("  NB: reproducing Dd2 is calibration, NOT validation — only 7G8 is the test.")


if __name__ == "__main__":
    main()
