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
    print("\n[calibration on Dd2 ONLY]  (fitted: kappa; all other inputs measured)")
    print(f"  observed fold(Dd2) = {D.IC50_DD2.value}/{D.IC50_WT.value} = "
          f"{D.IC50_DD2.value / D.IC50_WT.value:.2f}x")
    print(f"  Vmax/Km (Dd2) = {D.VMAX_DD2.value}/{D.KM_DD2.value} = "
          f"{D.VMAX_DD2.value / D.KM_DD2.value:.4f}  ->  kappa = {kappa:.2f}")

    # Frozen point prediction for 7G8 (uses 7G8 kinetics; NOT 7G8 IC50).
    fold_7g8 = B.predict_fold(kappa, D.VMAX_7G8.value, D.KM_7G8.value)
    ic50_7g8 = B.predict_ic50(kappa, D.IC50_WT.value, D.VMAX_7G8.value, D.KM_7G8.value)
    print("\n[frozen prediction: 7G8]  (7G8 kinetics measured, same experiment as Dd2)")
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
    group = {"km_7g8": "7G8 kinetics", "vmax_7g8": "7G8 kinetics",
             "km_dd2": "Dd2 kinetics", "vmax_dd2": "Dd2 kinetics",
             "ic50_wt": "IC50 (WT ref)", "ic50_dd2": "IC50 (Dd2 calib)"}
    for name, rho in res.sensitivity:
        print(f"  {name:10s} {rho:5.2f}  ({group.get(name, '')})")

    print("\n" + "-" * 72)
    print("Honest reading:")
    print("  First FALSIFIABLE result. With Dd2 and 7G8 kinetics measured in the SAME")
    print(f"  experiment, the transport-only model (calibrated on Dd2) predicts "
          f"{res.ic50_pred_median:.0f} nM")
    print(f"  (90% CI {res.ic50_pred_ci[0]:.0f}-{res.ic50_pred_ci[1]:.0f}), but the held-out "
          f"observation is {res.ic50_obs:.0f} nM -- OUTSIDE the interval:")
    print(f"  transport efficiency alone UNDER-predicts 7G8 resistance by ~{1/res.fold_error:.1f}x.")
    print("  This is the benchmark working as intended -- it is falsifiable, and for 7G8 the")
    print("  transport-only hypothesis is REJECTED. It matches the literature: 7G8 is a known")
    print("  transport-vs-resistance OUTLIER (excluding it raises R^2 0.86->0.98). 7G8 carries")
    print("  resistance beyond what its (low) CQ transport efficiency explains -- flagging")
    print("  additional biology to model next, not a coincidence to sell.")
    print("  NB: reproducing Dd2 is calibration, NOT validation; only 7G8 is the test.")


if __name__ == "__main__":
    main()
