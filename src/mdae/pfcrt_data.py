"""Sourced parameters for the chloroquine / PfCRT benchmark (slice 2).

Each value cites a primary source; the full traceable table lives in
``docs/parameters_slice1.md``. Transport kinetics are for CQ moved by PfCRT
expressed in Xenopus oocytes (assay pH 6.0). IC50 values are from ISOGENIC lines
differing only in the pfcrt allele (same assay) — so IC50 differences are
attributable to PfCRT, controlling for genetic background.

HONESTY on 7G8: its exact CQ Km/Vmax live in Summers et al. 2014 PNAS Table 1,
which is paywalled; no open source reproduces the numbers. 7G8 is documented as
belonging to the high-affinity (low Km) / low-capacity (low Vmax) group AND as an
OUTLIER in the transport-vs-resistance correlation. We therefore enter 7G8's
kinetics as ASSUMED bounds and propagate that ignorance — the predicted 7G8 IC50
is an interval, not a point. Swap in the measured Table-1 values to sharpen it.
"""

from __future__ import annotations

from .provenance import Param, Kind

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

# --- PfCRT CQ transport kinetics, Xenopus oocytes, pH 6.0 -------------------
# Dd2: MEASURED and robust (agree across Summers 2014, PMC6709616, PMC4276893).
KM_DD2 = Param(
    "Km_CQ_PfCRT_Dd2", 250.0, "uM", Kind.MEASURED, sd=40.0,
    source="Summers 2014 PNAS 10.1073/pnas.1322965111; confd PMC6709616, PMC4276893",
    note="250+/-40 uM (related work: 260+/-10, 250+/-30). CQ Km range across variants 117-293 uM.",
)
VMAX_DD2 = Param(
    "Vmax_CQ_PfCRT_Dd2", 33.0, "pmol/oocyte/h", Kind.MEASURED, sd=3.0,
    source="Summers 2014 PNAS; confd PMC6709616, PMC4276893", note="33+/-3 pmol/oocyte/h",
)

# 7G8: ASSUMED bounds (exact Table-1 values paywalled). High-affinity/low-capacity group.
KM_7G8 = Param(
    "Km_CQ_PfCRT_7G8", 183.0, "uM", Kind.ASSUMED, low=117.0, high=250.0,
    source="Summers 2014 PNAS group characterization (exact value paywalled)",
    note="high-affinity group -> Km below Dd2; bounded by observed variant range 117-250 uM",
)
VMAX_7G8 = Param(
    "Vmax_CQ_PfCRT_7G8", 20.0, "pmol/oocyte/h", Kind.ASSUMED, low=10.0, high=33.0,
    source="Summers 2014 PNAS group characterization (exact value paywalled)",
    note="low-capacity group -> Vmax below Dd2 (=33); 7G8 is a transport-vs-resistance OUTLIER",
)

# --- Isogenic-line CQ IC50 (measured, same assay, PfCRT-only difference) ----
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
