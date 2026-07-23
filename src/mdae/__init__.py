"""Malaria Drug Action Engine (MDAE).

Mechanistic pipeline for *Plasmodium falciparum*:

    chemical structure
      -> C_free at the site of action (Layer B, this slice)
      -> molecular function            (Layer C)
      -> stage-specific killing         (Layer D)

Vertical slice 1 (this release): Layer B for a weak-base drug — passive ion
trapping into the acidic digestive vacuole plus a transporter (PfCRT) efflux
term, validated against the analytical Henderson-Hasselbalch accumulation ratio.

Doctrine: mechanistic where parameterized; every parameter carries provenance
and uncertainty; retrodict then predict.
"""

from .speciation import WeakBase
from .oracle import accumulation_ratio
from .exposure import IonTrapModel, Compartment
from .function import occupancy, above_threshold
from .provenance import Param, Kind
from . import benchmark, benchmark_1b, pfcrt_data, pfdhfr_data, enzyme

__version__ = "0.1.0"

__all__ = [
    "WeakBase",
    "accumulation_ratio",
    "IonTrapModel",
    "Compartment",
    "occupancy",
    "above_threshold",
    "Param",
    "Kind",
    "benchmark",
    "pfcrt_data",
    "__version__",
]
