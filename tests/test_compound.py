import json

import pytest

from mdae.compound import CompoundProperties
from mdae.provenance import Kind, Param


SRC = "10.1186/s12936-019-3075-5"


def p(name, value, unit, *, kind=Kind.MEASURED, **kwargs):
    return Param(name, value, unit, kind, source=SRC if kind is not Kind.ASSUMED else "", **kwargs)


def test_minimal_record_keeps_missingness_explicit():
    c = CompoundProperties("example")
    assert c.parameters() == ()
    assert "pka" in c.missing()
    assert "intrinsic_clearance" in c.missing()


def test_complete_record_is_json_serializable():
    c = CompoundProperties(
        name="example",
        pka=(p("pKa1", 8.1, "pKa", sd=0.1),),
        logd_7_4=p("logD7.4", 2.2, "logD"),
        solubility=p("solubility", 50.0, "uM"),
        permeability=p("Papp", 12.0, "10^-6 cm/s"),
        fraction_unbound_plasma=p("fu_plasma", 0.12, "fraction"),
        fraction_unbound_medium=p("fu_medium", 0.25, "fraction"),
        blood_to_plasma_ratio=p("blood_plasma", 1.4, "ratio"),
        intrinsic_clearance=p("CLint", 3.0, "uL/min/mg"),
    )
    payload = c.to_dict()
    json.dumps(payload)
    assert payload["pka"][0]["kind"] == "measured"
    assert c.missing() == ()


def test_measured_property_requires_source():
    with pytest.raises(ValueError, match="require provenance"):
        CompoundProperties(
            "example",
            logd_7_4=Param("logD7.4", 2.0, "logD", Kind.MEASURED),
        )


def test_assumed_property_may_be_explicitly_unsourced():
    c = CompoundProperties(
        "example",
        solubility=Param(
            "solubility",
            100.0,
            "uM",
            Kind.ASSUMED,
            low=10.0,
            high=1000.0,
            note="placeholder bound for sensitivity analysis",
        ),
    )
    assert c.solubility.kind is Kind.ASSUMED


def test_unit_is_required():
    with pytest.raises(ValueError, match="unit must be explicit"):
        CompoundProperties("example", logd_7_4=p("logD", 2.0, ""))


@pytest.mark.parametrize("value", [-0.01, 1.01])
def test_fraction_unbound_must_be_physical(value):
    with pytest.raises(ValueError, match="between 0 and 1"):
        CompoundProperties(
            "example",
            fraction_unbound_plasma=p("fu", value, "fraction"),
        )


def test_nonnegative_properties_are_checked():
    with pytest.raises(ValueError, match="value must be >= 0"):
        CompoundProperties(
            "example",
            intrinsic_clearance=p("CLint", -1.0, "uL/min/mg"),
        )


def test_bounds_must_be_complete_and_ordered():
    with pytest.raises(ValueError, match="provided together"):
        CompoundProperties(
            "example",
            solubility=p("sol", 10.0, "uM", low=5.0),
        )

    with pytest.raises(ValueError, match="low must be <= high"):
        CompoundProperties(
            "example",
            solubility=p("sol", 10.0, "uM", low=20.0, high=5.0),
        )


def test_blood_to_plasma_ratio_must_be_positive():
    with pytest.raises(ValueError, match="value must be > 0"):
        CompoundProperties(
            "example",
            blood_to_plasma_ratio=p("bp", 0.0, "ratio"),
        )
