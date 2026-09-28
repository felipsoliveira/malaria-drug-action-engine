"""Traceable physicochemical properties for antimalarial exposure models.

This module separates *compound exposure properties* from target-binding scores.
The fields are intentionally small and auditable: each numerical quantity is a
:class:`mdae.provenance.Param` carrying units, uncertainty, epistemic kind and
primary-source provenance.

The initial schema follows the kinds of experimentally measured inputs commonly
used in antimalarial PK/exposure modelling (pKa, logD, solubility, permeability,
fraction unbound, blood/plasma partitioning and intrinsic clearance).

No property in this module is interpreted as a binding affinity or efficacy
measurement.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .provenance import Kind, Param


def _validate_param(param: Param) -> None:
    if not param.name.strip():
        raise ValueError("parameter name must be non-empty")
    if not param.unit.strip():
        raise ValueError(f"{param.name}: unit must be explicit")
    if param.sd is not None and param.sd < 0:
        raise ValueError(f"{param.name}: sd must be >= 0")
    if (param.low is None) ^ (param.high is None):
        raise ValueError(f"{param.name}: low and high must be provided together")
    if param.low is not None and param.high is not None and param.low > param.high:
        raise ValueError(f"{param.name}: low must be <= high")
    if param.kind is not Kind.ASSUMED and not param.source.strip():
        raise ValueError(
            f"{param.name}: measured/derived/fitted values require provenance"
        )


def _validate_fraction(param: Param, field_name: str) -> None:
    _validate_param(param)
    if not 0.0 <= param.value <= 1.0:
        raise ValueError(f"{field_name}: value must be between 0 and 1")
    if param.low is not None and param.low < 0.0:
        raise ValueError(f"{field_name}: lower bound must be >= 0")
    if param.high is not None and param.high > 1.0:
        raise ValueError(f"{field_name}: upper bound must be <= 1")


def _validate_nonnegative(param: Param, field_name: str) -> None:
    _validate_param(param)
    if param.value < 0.0:
        raise ValueError(f"{field_name}: value must be >= 0")
    if param.low is not None and param.low < 0.0:
        raise ValueError(f"{field_name}: lower bound must be >= 0")


def _param_to_dict(param: Param) -> dict[str, object]:
    return {
        "name": param.name,
        "value": float(param.value),
        "unit": param.unit,
        "kind": param.kind.value,
        "low": None if param.low is None else float(param.low),
        "high": None if param.high is None else float(param.high),
        "sd": None if param.sd is None else float(param.sd),
        "source": param.source,
        "note": param.note,
    }


@dataclass(frozen=True)
class CompoundProperties:
    """Exposure-relevant, provenance-carrying properties for one compound.

    Fields are optional because early-stage compounds rarely have a complete
    property panel. Missingness stays explicit rather than being silently filled
    by default estimates.
    """

    name: str
    pka: tuple[Param, ...] = ()
    logd_7_4: Param | None = None
    solubility: Param | None = None
    permeability: Param | None = None
    fraction_unbound_plasma: Param | None = None
    fraction_unbound_medium: Param | None = None
    blood_to_plasma_ratio: Param | None = None
    intrinsic_clearance: Param | None = None

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("compound name must be non-empty")

        for p in self.pka:
            _validate_param(p)

        if self.logd_7_4 is not None:
            _validate_param(self.logd_7_4)
        if self.solubility is not None:
            _validate_nonnegative(self.solubility, "solubility")
        if self.permeability is not None:
            _validate_nonnegative(self.permeability, "permeability")
        if self.fraction_unbound_plasma is not None:
            _validate_fraction(self.fraction_unbound_plasma, "fraction_unbound_plasma")
        if self.fraction_unbound_medium is not None:
            _validate_fraction(self.fraction_unbound_medium, "fraction_unbound_medium")
        if self.blood_to_plasma_ratio is not None:
            _validate_nonnegative(self.blood_to_plasma_ratio, "blood_to_plasma_ratio")
            if self.blood_to_plasma_ratio.value == 0.0:
                raise ValueError("blood_to_plasma_ratio: value must be > 0")
        if self.intrinsic_clearance is not None:
            _validate_nonnegative(self.intrinsic_clearance, "intrinsic_clearance")

    def parameters(self) -> tuple[Param, ...]:
        """Return all present parameters in a stable order."""
        fields: Iterable[Param | None] = (
            self.logd_7_4,
            self.solubility,
            self.permeability,
            self.fraction_unbound_plasma,
            self.fraction_unbound_medium,
            self.blood_to_plasma_ratio,
            self.intrinsic_clearance,
        )
        return self.pka + tuple(p for p in fields if p is not None)

    def missing(self) -> tuple[str, ...]:
        """Return names of optional property groups that are not yet populated."""
        missing: list[str] = []
        if not self.pka:
            missing.append("pka")
        for name in (
            "logd_7_4",
            "solubility",
            "permeability",
            "fraction_unbound_plasma",
            "fraction_unbound_medium",
            "blood_to_plasma_ratio",
            "intrinsic_clearance",
        ):
            if getattr(self, name) is None:
                missing.append(name)
        return tuple(missing)

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable benchmark/input representation."""
        return {
            "name": self.name,
            "pka": [_param_to_dict(p) for p in self.pka],
            "logd_7_4": None if self.logd_7_4 is None else _param_to_dict(self.logd_7_4),
            "solubility": None if self.solubility is None else _param_to_dict(self.solubility),
            "permeability": None if self.permeability is None else _param_to_dict(self.permeability),
            "fraction_unbound_plasma": (
                None
                if self.fraction_unbound_plasma is None
                else _param_to_dict(self.fraction_unbound_plasma)
            ),
            "fraction_unbound_medium": (
                None
                if self.fraction_unbound_medium is None
                else _param_to_dict(self.fraction_unbound_medium)
            ),
            "blood_to_plasma_ratio": (
                None
                if self.blood_to_plasma_ratio is None
                else _param_to_dict(self.blood_to_plasma_ratio)
            ),
            "intrinsic_clearance": (
                None
                if self.intrinsic_clearance is None
                else _param_to_dict(self.intrinsic_clearance)
            ),
        }
