"""PfDHFR / pyrimethamine panel — round-trip self-consistency (NOT a benchmark).

Reconstructs the IC50 implied by each variant's published Ki and Km at the assay
[DHF]=100 uM, via Cheng-Prusoff. Because the Table-3 Ki were themselves derived
from IC50 by the same relation, this reproduces the input by algebra: it checks
UNITS / IMPLEMENTATION / PROVENANCE, not the biology. No pass/fail, no benchmark.

Run:  PYTHONPATH=src python3 examples/pfdhfr_panel.py
"""

from __future__ import annotations

from mdae import enzyme as E
from mdae import pfdhfr_data as D


def main() -> None:
    S = D.ASSAY_DHF_UM
    print("=" * 78)
    print("PfDHFR + pyrimethamine (Sirawaraporn 1997 Table 3) — ROUND-TRIP, NOT a benchmark")
    print(f"assay: [DHF]={S} uM, [NADPH]={D.ASSAY_NADPH_UM} uM, pH {D.ASSAY_PH}, {D.ASSAY_TEMP_C} C")
    print("=" * 78)
    print(f"\n  {'variant':26s} {'kcat':>6s} {'Km(uM)':>8s} {'Ki(nM)':>9s} "
          f"{'IC50=Ki(1+S/Km)':>17s} {'IC50/Ki':>8s}")
    for v in D.PANEL:
        ki, km, kcat = v.ki_pyr.value, v.km_dhf.value, v.kcat.value
        ic50 = E.ic50_cheng_prusoff(ki, S, km)
        print(f"  {v.name:26s} {kcat:6.1f} {km:8.1f} {ki:9.1f} {ic50:15.0f} nM {ic50/ki:7.1f}x")

    print("\n" + "-" * 78)
    print("Reading (self-consistency only):")
    print("  * The Ki span ~1.5 -> 859 nM across the resistance series (two routes to the quad).")
    print("  * IC50/Ki = 1 + [S]/Km varies 6-18x because Km varies (6-25 uM): the SAME Ki maps")
    print("    to different functional IC50 depending on Km -- that Km-dependence is the Layer-C")
    print("    conversion the model encodes. Useful as an illustration, NOT as validation:")
    print("  * these Ki were BACK-CALCULATED from IC50 with this very formula, so the round-trip")
    print("    is algebra. A falsifiable #2a must freeze Ki,Km and predict INDEPENDENT data")
    print("    (IC50 at a different [DHF], or a directly-measured Kd + a functional IC50).")


if __name__ == "__main__":
    main()
