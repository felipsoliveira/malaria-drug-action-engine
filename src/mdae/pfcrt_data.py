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

# Held-out strain kinetics (Summers 2014 Table 1, same experiment). These are
# transport INPUTS, not the held-out phenotype (which is the RI, from Table S2).
KM_K1 = Param("Km_CQ_PfCRT_K1", 293.0, "uM", Kind.MEASURED, sd=10.0, source=_SUMMERS, note="293+/-10, n=6")
VMAX_K1 = Param("Vmax_CQ_PfCRT_K1", 79.0, "pmol/oocyte/h", Kind.MEASURED, sd=17.0, source=_SUMMERS, note="79+/-17, n=6")
KM_GB4 = Param("Km_CQ_PfCRT_GB4", 275.0, "uM", Kind.MEASURED, sd=16.0, source=_SUMMERS, note="275+/-16, n=7")
VMAX_GB4 = Param("Vmax_CQ_PfCRT_GB4", 57.0, "pmol/oocyte/h", Kind.MEASURED, sd=5.0, source=_SUMMERS, note="57+/-5, n=7")
KM_ECU = Param("Km_CQ_PfCRT_Ecu1110", 191.0, "uM", Kind.MEASURED, sd=17.0, source=_SUMMERS, note="191+/-17, n=5")
VMAX_ECU = Param("Vmax_CQ_PfCRT_Ecu1110", 36.0, "pmol/oocyte/h", Kind.MEASURED, sd=4.0, source=_SUMMERS, note="36+/-4, n=5")

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

# --- Table S2 CQ resistance indices (Summers thesis Table 6.3 = PNAS Table S2) ---
# HETEROGENEOUS COMPILATION. Per the table's own note: each RI uses the CQS strain
# measured in that same study; when a study has several CQS, their mean is used; then
# results from DIFFERENT studies are aggregated. There is NO common CQS denominator,
# assay, or laboratory. This fails the pre-registered comparability audit -> the
# confirmatory Benchmark #1b is NOT EVALUABLE with these phenotypes (stop rule).
# Kept here only for the clearly-separated EXPLORATORY analysis. See
# docs/benchmark_1b_result.md. Source: Summers PhD thesis (ANU), Table 6.3, p.162.
_S2 = "Summers PhD thesis (ANU) Table 6.3 = PNAS Table S2; HETEROGENEOUS compilation (per-study CQS)"
RI_DD2 = Param("RI_CQ_Dd2", 11.13, "dimensionless", Kind.DERIVED, sd=1.40, source=_S2, note="n=8 studies; NOT mutually comparable")
RI_K1 = Param("RI_CQ_K1", 11.79, "dimensionless", Kind.DERIVED, sd=2.01, source=_S2, note="n=4 studies")
RI_GB4 = Param("RI_CQ_GB4", 8.56, "dimensionless", Kind.DERIVED, sd=1.70, source=_S2, note="n=3 studies")
RI_ECU = Param("RI_CQ_Ecu1110", 6.50, "dimensionless", Kind.DERIVED, sd=0.89, source=_S2, note="n=2 studies")
RI_7G8 = Param("RI_CQ_7G8", 6.63, "dimensionless", Kind.DERIVED, sd=0.60, source=_S2, note="n=9 studies")
RI_PH1 = Param("RI_CQ_Ph1", 5.04, "dimensionless", Kind.DERIVED, sd=0.0, source=_S2, note="n=1 study (no SEM)")
