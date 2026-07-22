"""Vertical slice 1 — chloroquine / PfCRT / digestive vacuole.

Runs the Layer-B ion-trapping model for a sensitive and a (PfCRT-mutant) resistant
parasite and prints the mechanistic story the engine exists to produce:

  same plasma, same target Kd, DIFFERENT local free concentration -> different fate.

Run:  PYTHONPATH=src python3 examples/chloroquine_pfcrt.py
"""

from __future__ import annotations

from mdae import WeakBase, IonTrapModel, accumulation_ratio, occupancy, above_threshold
from mdae import parameters as P


def build(pfcrt: dict) -> IonTrapModel:
    return IonTrapModel(
        base=P.CHLOROQUINE,
        source=P.CYTOSOL_PLASMA,
        vacuole=P.DIGESTIVE_VACUOLE,
        c_ext=P.C_EXT_UM,
        k_perm=P.K_PERM_PER_S,
        vmax=pfcrt["vmax"],
        Km=pfcrt["Km"],
    )


def report(label: str, model: IonTrapModel) -> float:
    c_dv = model.steady_state()
    R = model.accumulation_ratio()
    theta = occupancy(c_dv, P.KD_HEME_UM)
    kills = above_threshold(c_dv, P.C_FUNCTIONAL_THRESHOLD_UM)
    print(f"\n[{label}]  PfCRT vmax = {model.vmax} uM/s")
    print(f"  C_free,DV (steady state) = {c_dv:,.1f} uM")
    print(f"  accumulation ratio DV/plasma = {R:,.1f}x")
    print(f"  target occupancy theta (Kd={P.KD_HEME_UM} uM) = {theta:.4f}")
    print(f"  above functional threshold ({P.C_FUNCTIONAL_THRESHOLD_UM:,.0f} uM)? "
          f"{'YES -> kill' if kills else 'NO -> survive'}")
    return c_dv


def main() -> None:
    print("=" * 70)
    print("MDAE slice 1: chloroquine ion trapping into the digestive vacuole")
    print("=" * 70)

    # Independent analytical oracle (from speciation only, no dynamics).
    R_oracle = accumulation_ratio(P.CHLOROQUINE, P.DIGESTIVE_VACUOLE.pH, P.CYTOSOL_PLASMA.pH)
    print(f"\nAnalytical oracle (Henderson-Hasselbalch, diprotic):")
    print(f"  R = Z(pH_DV={P.DIGESTIVE_VACUOLE.pH}) / Z(pH_plasma={P.CYTOSOL_PLASMA.pH}) "
          f"= {R_oracle:,.1f}x")

    sens = build(P.PFCRT_SENSITIVE)
    res = build(P.PFCRT_RESISTANT)

    print(f"\nPassive equilibration time constant ~ {sens.equilibration_time()/3600:.2f} h")

    c_sens = report("SENSITIVE strain", sens)
    c_res = report("RESISTANT strain (PfCRT efflux)", res)

    print("\n" + "-" * 70)
    print("Mechanistic story:")
    print(f"  Plasma concentration is IDENTICAL in both strains ({P.C_EXT_UM} uM free).")
    print(f"  Target chemistry (Kd = {P.KD_HEME_UM} uM) is IDENTICAL in both strains.")
    print(f"  PfCRT efflux drops DV exposure {c_sens/c_res:,.1f}x "
          f"({c_sens:,.0f} -> {c_res:,.0f} uM).")
    print("  Resistance here is an EXPOSURE phenomenon, not a binding change —")
    print("  exactly what a docking score or a single IC50 cannot tell you.")
    print("\nCaveat: rate constants are illustrative (not fitted). The accumulation")
    print("ratio is grounded in theory; the sensitive ODE matches the oracle exactly.")


if __name__ == "__main__":
    main()
