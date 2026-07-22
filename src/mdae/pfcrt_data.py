"""Sourced parameters for the chloroquine / PfCRT benchmark (slice 2).

Each value cites a primary source; the full traceable table lives in
``docs/parameters_pfcrt.md``. Transport kinetics are for CQ moved by PfCRT
expressed in Xenopus oocytes (acidic medium, ~pH 6.0). IC50 values are from
ISOGENIC lines differing only in the pfcrt allele (same assay) — so IC50
differences are attributable to PfCRT, controlling for genetic background.

CRITICAL — same-experiment matching: Vmax is an absolute rate whose scale varies
between assays/batches (e.g. Dd2 Vmax is 61 pmol/oocyte/h in Summers 2014 Table 1
but 33 in the 2019 phosphomimetic study). Because the benchmark compares the
transport efficiency Vmax/Km ACROSS strains, Dd2 and 7G8 MUST come from the SAME
experiment. We therefore take BOTH from Summers et al. 2014 PNAS Table 1
(verified in the ANU open-access PDF, p. E1762). 7G8 is thus MEASURED here, not
assumed. The paper also notes 7G8 is an outlier: excluding it raises the
resistance-index vs transport correlation R^2 from ~0.86 to ~0.98 (Vmax).
"""

from __future__ import annotations

from .provenance import Param, Kind

_SUMMERS = "Summers et al. 2014 PNAS 10.1073/pnas.1322965111, Table 1 (p.E1762; ANU open PDF)"

# --- Compartment pH (measured) ----------------------------------------------
PH_DV = Param(
    "pH_digestive_vacuole", 5.2, "pH units", Kind.MEASURED, low=5.0, high=5.6,
    source="Kuhn et al. 2007 Cell Microbiol 10.1111/j.1462-5822.2006.00847.x",
    note="pHluorin 5.18+/-0.05; 5.0-5.6 across methods. DV pH does NOT differ CQS vs CQR.",
)
PH_CYTOSOL = Param(
    "pH_parasite_cytosol", 7.15, "pH units", Kind.MEASURED, sd=0.07,
    source="Kuhn et al. 2007 Cell Microbiol",
    note="pHluorin 7.15+/-0.07. This (not plasma 7.4) is the DV's immediate source compartment.",
)

# --- Chloroquine acid-base (measured) ---------------------------------------
PKA1_CQ = Param(
    "chloroquine_pKa1", 10.1, "pKa", Kind.MEASURED, low=10.0, high=10.3,
    source="Warhurst et al.; CQ is a diprotic weak base",
    note="side-chain diethylamino N (first protonation of the neutral base)",
)
PKA2_CQ = Param(
    "chloroquine_pKa2", 8.1, "pKa", Kind.MEASURED, low=8.0, high=8.4,
    source="Warhurst et al.", note="quinoline ring N (second protonation)",
)

# --- PfCRT CQ transport kinetics, Xenopus oocytes, ~pH 6.0, SAME experiment --
# Both from Summers 2014 Table 1 (mean +/- SEM; n = number of experiments).
KM_DD2 = Param(
    "Km_CQ_PfCRT_Dd2", 232.0, "uM", Kind.MEASURED, sd=11.0,
    source=_SUMMERS, note="232+/-11 uM, n=15. NB other assays report ~250 (do not mix with 7G8).",
)
VMAX_DD2 = Param(
    "Vmax_CQ_PfCRT_Dd2", 61.0, "pmol/oocyte/h", Kind.MEASURED, sd=6.0,
    source=_SUMMERS, note="61+/-6 pmol/oocyte/h, n=15 (assay-dependent scale; matched to 7G8)",
)
KM_7G8 = Param(
    "Km_CQ_PfCRT_7G8", 117.0, "uM", Kind.MEASURED, sd=6.0,
    source=_SUMMERS, note="117+/-6 uM, n=6. High-affinity (lowest Km of the CQR variants).",
)
VMAX_7G8 = Param(
    "Vmax_CQ_PfCRT_7G8", 9.0, "pmol/oocyte/h", Kind.MEASURED, sd=1.0,
    source=_SUMMERS, note="9+/-1 pmol/oocyte/h, n=6. Low-capacity; 7G8 is a transport-vs-resistance outlier.",
)

# --- Isogenic-line CQ IC50 (measured, same assay, PfCRT-only difference) -----
# From PMC9067703 Table 1: transfectants in the GC03 background.
IC50_WT = Param(
    "IC50_CQ_isogenic_WT_GC03", 25.0, "nM", Kind.MEASURED, sd=1.1,
    source="PMC9067703 Table 1 (line C2GC03, wild-type PfCRT)", note="25.0+/-1.1 nM; sensitive reference",
)
IC50_DD2 = Param(
    "IC50_CQ_isogenic_Dd2", 132.0, "nM", Kind.MEASURED, sd=4.9,
    source="PMC9067703 Table 1 (line C4Dd2)", note="132+/-4.9 nM; CALIBRATION strain",
)
IC50_7G8 = Param(
    "IC50_CQ_isogenic_7G8", 84.0, "nM", Kind.MEASURED, sd=3.5,
    source="PMC9067703 Table 1 (line C67G8)",
    note="84+/-3.5 nM; HELD-OUT phenotype — never used to calibrate the model",
)
