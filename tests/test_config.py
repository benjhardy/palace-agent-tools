import pytest

from palace_tools import schema
from palace_tools.config import build_config, count_frequencies, excitations, load_tags

TAGS = {"volumes": {"substrate": {"tag": 1}, "air": {"tag": 2}},
        "surfaces": {"ground": {"tag": 101}, "strip": {"tag": 102}, "port1": {"tag": 103}, "port2": {"tag": 104},
                     "cap": {"tag": 105}, "res": {"tag": 106}, "sym": {"tag": 107}, "abc": {"tag": 108}}}


def base(**kw):
    s = {"problem": "driven", "mesh": "m.msh", "units": "mm", "tags": TAGS, "order": 2,
         "materials": [{"tags": ["substrate"], "eps_r": 3.55, "tan_d": 0.0027}, {"tags": ["air"]}],
         "pec": ["ground", "strip"], "absorbing": {"tags": ["abc"]},
         "lumped_ports": [{"index": 1, "tags": ["port1"], "direction": "+Z", "R": 50},
                          {"index": 2, "tags": ["port2"], "direction": "+Z", "R": 50}],
         "frequency": {"min_ghz": 1, "max_ghz": 10, "step_ghz": 0.5, "adaptive_tol": 1e-3}}
    s.update(kw)
    return s


def test_driven_defaults():
    c = build_config(base())
    assert c["Model"]["L0"] == 1e-3
    lin = c["Solver"]["Linear"]
    assert lin["Type"] == "MUMPS" and lin["MGMaxLevels"] == 1 and lin["ComplexCoarseSolve"] is True
    assert lin["EstimatorMG"] is True and lin["EstimatorTol"] == 1e-4
    assert [p["Excitation"] for p in c["Boundaries"]["LumpedPort"]] == [1, 2]
    assert c["Boundaries"]["PEC"]["Attributes"] == [101, 102]
    assert c["Problem"]["OutputFormats"] == {"Paraview": False, "GridFunction": False}
    assert count_frequencies(c) == 19
    assert excitations(c) == {1: [("LumpedPort", 1)], 2: [("LumpedPort", 2)]}
    assert schema.validate(c)["valid"] is True


def test_auto_preset_switches_to_multigrid_for_large_models():
    c = build_config(base(solver={"unknowns": 2_000_000}))
    assert "MGMaxLevels" not in c["Solver"]["Linear"]


def test_lumped_element_cut_scaling():
    s = base(lumped_elements=[{"index": 3, "tags": ["cap"], "direction": "+Y", "C": 1.8e-12, "cut": "half-width"},
                              {"index": 4, "tags": ["res"], "direction": "+X", "R": 100, "L": 1e-9,
                               "cut": "half-length"}])
    lp = {p["Index"]: p for p in build_config(s)["Boundaries"]["LumpedPort"]}
    assert lp[3]["C"] == pytest.approx(0.9e-12)
    assert lp[4]["R"] == pytest.approx(50.0) and lp[4]["L"] == pytest.approx(0.5e-9)
    assert "Excitation" not in lp[3]


def test_impedance_sheet_from_r():
    c = build_config(base(impedance=[{"tags": ["res"], "R": 100.0, "length": 0.5, "width": 0.25}]))
    assert c["Boundaries"]["Impedance"][0]["Rs"] == pytest.approx(50.0)


def test_symmetry_planes_go_to_pec_or_pmc():
    c = build_config(base(symmetry=[{"tags": ["sym"], "type": "pmc"}]))
    assert c["Boundaries"]["PMC"]["Attributes"] == [107]
    c = build_config(base(symmetry=[{"tags": ["sym"], "type": "pec"}]))
    assert 107 in c["Boundaries"]["PEC"]["Attributes"]


def test_eigen_target_from_expected():
    s = {"problem": "eigenmode", "mesh": "c.msh", "units": "mm", "materials": [{"tags": [1]}], "pec": [101],
         "eigen": {"expected_ghz": 8.0, "n": 3}}
    c = build_config(s)
    assert c["Solver"]["Eigenmode"]["Target"] == 4.0
    assert c["Solver"]["Linear"]["MGMaxLevels"] == 1
    assert "ComplexCoarseSolve" not in c["Solver"]["Linear"]
    assert schema.validate(c)["valid"] is True


def test_wave_ports_and_points():
    s = base(lumped_ports=[], wave_ports=[{"index": 1, "tags": ["port1"], "offset": 2.0},
                                          {"index": 2, "tags": ["port2"]}],
             frequency={"points_ghz": [1, 2, 3]})
    c = build_config(s)
    assert c["Boundaries"]["WavePort"][0]["Offset"] == 2.0
    assert c["Solver"]["Driven"]["Samples"][0] == {"Type": "Point", "Freq": [1.0, 2.0, 3.0], "SaveStep": 0}
    assert c["Solver"]["Driven"]["AdaptiveTol"] == 0.0
    assert schema.validate(c)["valid"] is True


def test_errors():
    with pytest.raises(ValueError, match="units"):
        build_config({k: v for k, v in base().items() if k != "units"})
    with pytest.raises(KeyError, match="unknown physical group"):
        build_config(base(pec=["nope"]))
    with pytest.raises(ValueError, match="target_ghz"):
        build_config({"problem": "eigenmode", "mesh": "c.msh", "units": "mm", "materials": [{"tags": [1]}]})


def test_extra_deep_merge_and_amr():
    c = build_config(base(extra={"Solver": {"Linear": {"MaxIts": 999}}}, amr={"max_its": 3, "max_size": 400000}))
    assert c["Solver"]["Linear"]["MaxIts"] == 999 and c["Solver"]["Linear"]["Type"] == "MUMPS"
    assert c["Model"]["Refinement"]["Nonconformal"] is False
    assert "EstimatorTol" not in c["Solver"]["Linear"]  # keep the accurate estimator when AMR uses it


def test_load_tags_flat():
    assert load_tags({"a": 1, "b": {"tag": 2}}) == {"a": 1, "b": 2}
