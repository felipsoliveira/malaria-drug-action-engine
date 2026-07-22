"""Parameters for slice 1 (chloroquine - PfCRT - digestive vacuole).

Doctrine: every parameter states provenance (measured / estimated / illustrative)
and, where it matters, uncertainty. Slice 1 is a MECHANISM demonstration and a
correctness check against the analytical oracle — it is NOT yet fitted to
chloroquine potency/resistance data. Do not read the absolute rate constants as
calibrated; the accumulation *ratio* is what is grounded in theory.
"""

from __future__ import annotations

from .speciation import WeakBase
from .exposure import Compartment

# --- Drug: chloroquine ------------------------------------------------------
# pKa ordered first-proton-first (most basic site protonates first, highest pKa).
# Provenance: experimental, literature (Warhurst et al.; diprotic weak base).
CHLOROQUINE = WeakBase(name="chloroquine", pKa=(10.1, 8.1))

# --- Compartments -----------------------------------------------------------
# Digestive vacuole pH ~5.0-5.4 (measured, literature). Cytosol/plasma ~7.4 (measured).
DIGESTIVE_VACUOLE = Compartment(name="digestive_vacuole", pH=5.2)  # measured, +/-0.2
CYTOSOL_PLASMA = Compartment(name="cytosol_plasma", pH=7.4)         # measured

# --- Concentrations / rates (ESTIMATED / ILLUSTRATIVE for slice 1) ----------
C_EXT_UM = 1.0        # free source-compartment chloroquine [uM] (illustrative)
K_PERM_PER_S = 1.0e4  # lumped P*A/V [1/s] (estimated; sets ~hour equilibration)

# PfCRT efflux of the charged DV species. vmax=0 => sensitive strain.
# Resistant values chosen to drop DV exposure below the functional threshold;
# illustrative, NOT fitted to a specific PfCRT haplotype.
PFCRT_SENSITIVE = dict(vmax=0.0, Km=500.0)
PFCRT_RESISTANT = dict(vmax=5.0, Km=500.0)  # [uM/s], [uM]

# --- Target (Layer C) -------------------------------------------------------
# Chloroquine-heme interaction / heme-detox interference. Illustrative.
KD_HEME_UM = 1.0           # illustrative target Kd [uM]
C_FUNCTIONAL_THRESHOLD_UM = 5.0e3  # free DV conc needed to interfere enough to kill [uM] (illustrative)
