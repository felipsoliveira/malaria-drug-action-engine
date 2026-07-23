"""PfDHFR + pyrimethamine kinetic panel (Sirawaraporn et al. 1997 PNAS, Table 3).

A SINGLE experimental system (same lab, same assay) covering the pyrimethamine-
resistance mutation series: kcat, Km(dihydrofolate), and Ki(pyrimethamine). Assay:
100 uM DHF, 100 uM NADPH, pH 7.0, 25 C (base method: Sirawaraporn et al. 1993 JBC).

Source: Sirawaraporn et al., PNAS 1997, Table 3 (PMC19755).

IMPORTANT provenance caveat (why a round-trip is not a benchmark): the Ki values in
Table 3 were DERIVED from the measured IC50 via the Cheng-Prusoff relation
IC50 = Ki*(1 + [S]/Km). Reconstructing IC50 from these Ki with `enzyme.py` therefore
reproduces the input by algebra -- it validates UNITS / IMPLEMENTATION / PROVENANCE,
NOT the biology of binding->function. A falsifiable #2a needs INDEPENDENT data
(e.g. an activity/IC50 curve at a different [DHF], or a directly-measured Kd/Ki plus a
functional IC50).

Do NOT merge with Sandefur et al. 2007 (PMC2020854): different expression construct;
WT Km shifts 13 -> 44.4 uM and kcat 88 -> 1.24 /s -- the two are not interchangeable
as a single panel.
"""

from __future__ import annotations

from dataclasses import dataclass

from .provenance import Param, Kind

_SIRA = "Sirawaraporn et al. 1997 PNAS Table 3 (PMC19755)"

# Standard assay context in which the Ki were defined (Cheng-Prusoff [S]).
ASSAY_DHF_UM = 100.0      # [S] dihydrofolate
ASSAY_NADPH_UM = 100.0
ASSAY_PH = 7.0
ASSAY_TEMP_C = 25.0


@dataclass(frozen=True)
class DHFRVariant:
    name: str
    kcat: Param     # 1/s
    km_dhf: Param   # uM
    ki_pyr: Param   # nM (DERIVED from IC50 via Cheng-Prusoff)


def _v(name, kcat, km, km_sd, ki, ki_sd) -> DHFRVariant:
    return DHFRVariant(
        name,
        Param(f"kcat_{name}", kcat, "1/s", Kind.MEASURED, source=_SIRA, note="no SEM reported in Table 3"),
        Param(f"Km_DHF_{name}", km, "uM", Kind.MEASURED, sd=km_sd, source=_SIRA),
        Param(f"Ki_pyr_{name}", ki, "nM", Kind.DERIVED, sd=ki_sd, source=_SIRA,
              note="derived from IC50 via Cheng-Prusoff; round-trip is NOT independent validation"),
    )


# Ordered along the resistance gradient (two evolutionary routes converge on the quadruple).
PANEL: list[DHFRVariant] = [
    _v("WT", 88.0, 13.0, 5.0, 1.5, 0.2),
    _v("S108N", 92.0, 25.0, 9.0, 13.0, 4.0),
    _v("N51I+S108N", 77.0, 6.0, 1.0, 37.0, 6.0),
    _v("C59R+S108N", 4.2, 24.0, 1.0, 72.0, 3.0),
    _v("N51I+C59R+S108N", 3.2, 6.2, 0.3, 120.0, 5.0),
    _v("C59R+S108N+I164L", 3.0, 18.0, 1.0, 383.0, 33.0),
    _v("N51I+C59R+S108N+I164L", 15.0, 14.0, 1.0, 859.0, 117.0),
]

BY_NAME = {v.name: v for v in PANEL}
