import copy

import pytest

from palace_tools.config import build_config
from palace_tools.lint import RULES, lint

TAGS = {"volumes": {"substrate": {"tag": 1}, "air": {"tag": 2}},
        "surfaces": {"ground": {"tag": 101}, "strip": {"tag": 102}, "port1": {"tag": 103}, "port2": {"tag": 104},
                     "abc": {"tag": 105}}}


def good():
    return build_config({
        "problem": "driven", "mesh": "m.msh", "units": "mm", "tags": TAGS, "order": 2,
        "materials": [{"tags": ["substrate"], "eps_r": 3.55}, {"tags": ["air"]}],
        "pec": ["ground", "strip"], "absorbing": {"tags": ["abc"]},
        "lumped_ports": [{"index": 1, "tags": ["port1"], "direction": "+Z", "R": 50},
                         {"index": 2, "tags": ["port2"], "direction": "+Z", "R": 50}],
        "frequency": {"min_ghz": 1, "max_ghz": 10, "step_ghz": 0.5, "adaptive_tol": 1e-3}})


def ids(res, severity=None):
    return {f["id"] for f in res["findings"] if severity is None or f["severity"] == severity}


def test_good_config_is_clean():
    r = lint(good())
    assert r["ok"] and r["counts"]["error"] == 0 and r["counts"]["warning"] == 0


def test_every_rule_has_lesson_and_source():
    assert len(RULES) >= 35
    for rid, r in RULES.items():
        assert r["lesson"] and r["source"] and r["title"], rid


def test_superlu_default():
    c = good()
    c["Solver"]["Linear"]["Type"] = "Default"
    assert "SOL001" in ids(lint(c), "warning")
    c["Solver"]["Linear"].pop("Type")
    assert "SOL001" in ids(lint(c))


def test_mixed_ports_version_dependent():
    c = good()
    c["Boundaries"]["WavePort"] = [{"Index": 3, "Attributes": [105], "Excitation": 3}]
    c["Boundaries"].pop("Absorbing")
    assert "PRT001" in ids(lint(c), "error")
    assert "PRT001" in ids(lint(c, context={"palace_version": "0.17.0"}), "info")


def test_excitation_index_mismatch_and_grouped():
    c = good()
    c["Boundaries"]["LumpedPort"][1]["Excitation"] = 5
    assert "PRT002" in ids(lint(c), "error")
    c["Boundaries"]["LumpedPort"][1]["Excitation"] = 1
    r = lint(c)
    assert "PRT003" in ids(r) and "PRT002" not in ids(r)


def test_bool_excitation_single_port_ok():
    c = good()
    for p in c["Boundaries"]["LumpedPort"]:
        p.pop("Excitation")
    c["Boundaries"]["LumpedPort"][1]["Excitation"] = True
    assert "PRT002" not in ids(lint(c))


def test_duplicate_index_and_no_excitation():
    c = good()
    c["Boundaries"]["WavePort"] = [{"Index": 1, "Attributes": [105]}]
    assert "PRT005" in ids(lint(c), "error")
    c = good()
    for p in c["Boundaries"]["LumpedPort"]:
        p.pop("Excitation")
    assert "PRT006" in ids(lint(c), "error")


def test_lumped_values():
    c = good()
    c["Boundaries"]["LumpedPort"].append({"Index": 3, "Attributes": [105], "Direction": "+X", "C": 1.8})
    assert "PRT010" in ids(lint(c), "error")
    c["Boundaries"]["LumpedPort"][-1] = {"Index": 3, "Attributes": [105], "Direction": "+X", "C": 1e-12, "Cs": 1e-12}
    assert "PRT007" in ids(lint(c), "error")
    c["Boundaries"]["LumpedPort"][-1] = {"Index": 3, "Attributes": [105], "Direction": "+X", "R": 50, "C": 1e-12,
                                         "Excitation": 3}
    assert "PRT008" in ids(lint(c), "error")
    assert "PRT008" not in ids(lint(c, context={"palace_version": "0.18"}))


def test_unequal_r_in_016():
    c = good()
    c["Boundaries"]["LumpedPort"][1]["R"] = 400.0
    assert "PRT011" in ids(lint(c), "warning")
    assert "PRT011" not in ids(lint(c, context={"palace_version": "v0.17.0"}))


def test_units_and_output():
    c = good()
    c["Model"].pop("L0")
    c["Problem"].pop("Output")
    r = lint(c)
    assert "MOD001" in ids(r, "warning") and "MOD003" in ids(r, "error")


def test_frequency_units():
    c = good()
    c["Solver"]["Driven"]["Samples"][0]["MaxFreq"] = 1e10
    assert "DRV001" in ids(lint(c), "error")


def test_eigen_target():
    c = build_config({"problem": "eigenmode", "mesh": "c.msh", "units": "mm", "materials": [{"tags": [1]}],
                      "pec": [101], "eigen": {"target_ghz": 9.0, "n": 3}})
    assert "EIG001" in ids(lint(c, context={"spec": {"eigen": {"expected_ghz": 8.2}}}), "error")
    c["Solver"]["Eigenmode"].pop("Target")
    assert "EIG001" in ids(lint(c), "error")
    c["Solver"]["Eigenmode"]["Target"] = 4.0
    c["Solver"]["Linear"].pop("MGMaxLevels")
    c["Boundaries"]["LumpedPort"] = [{"Index": 1, "Attributes": [102], "Direction": "+Z", "C": 1e-12}]
    assert "SOL002" in ids(lint(c), "warning")


def test_amr_rules():
    c = good()
    c["Model"]["Refinement"] = {"MaxIts": 4}
    c["Solver"]["Driven"]["AdaptiveTol"] = 0.0
    c["Boundaries"]["WavePort"] = [{"Index": 3, "Attributes": [105]}]
    r = lint(c)
    assert {"AMR001", "AMR002", "PRT014"} <= ids(r, "warning")


def test_conductivity_thickness_and_double_bc():
    c = good()
    c["Boundaries"]["Conductivity"] = [{"Attributes": [102], "Conductivity": 5.8e7, "Thickness": 0.0005}]
    r = lint(c)
    assert "BND003" in ids(r, "warning")
    assert "BND007" in ids(r, "error")  # 102 is both PEC and Conductivity


def test_closed_box_adaptive():
    c = good()
    c["Boundaries"].pop("Absorbing")
    c["Solver"]["Linear"]["Type"] = "SuperLU"
    assert "BND002" in ids(lint(c), "warning")


def test_field_saving_with_many_excitations():
    c = good()
    c["Solver"]["Driven"]["Samples"][0]["SaveStep"] = 1
    r = lint(c)
    assert "PRT016" in ids(r, "warning") and "MOD004" in ids(r)


# ----------------------------------------------------------------------------------------- geometry rules

def geometry():
    return {"domain_bbox": [0, -5, 0, 20, 5, 5], "structure_bbox": [0, -0.6, 0, 20, 0.6, 0.508],
            "volumes": {"substrate": {"tag": 1}, "air": {"tag": 2}},
            "surfaces": {
                "ground": {"tag": 101, "bbox": [0, -5, 0, 20, 5, 0], "outer": True, "sides": ["-z"]},
                "strip": {"tag": 102, "bbox": [0, -0.6, 0.508, 20, 0.6, 0.508], "outer": False, "sides": [],
                          "planar": True, "axes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "lengths": [20, 1.2, 0]},
                "port1": {"tag": 103, "bbox": [0, -0.6, 0, 0, 0.6, 0.508], "outer": True, "sides": ["-x"],
                          "planar": True, "axes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "lengths": [0, 1.2, 0.508]},
                "port2": {"tag": 104, "bbox": [20, -0.6, 0, 20, 0.6, 0.508], "outer": True, "sides": ["+x"],
                          "planar": True, "axes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "lengths": [0, 1.2, 0.508]},
                "abc": {"tag": 105, "bbox": [0, -5, 0, 20, 5, 5], "outer": True, "sides": ["-y", "+y", "+z"]}}}


def test_geometry_clean_and_skipped_without():
    r = lint(good(), geometry())
    assert r["counts"]["error"] == 0 and not r["skipped_rules_need_geometry"]
    assert "PRT009" in lint(good())["skipped_rules_need_geometry"]


def test_direction_along_normal():
    c = good()
    c["Boundaries"]["LumpedPort"][0]["Direction"] = "+X"
    assert "PRT009" in ids(lint(c, geometry()), "error")
    c["Boundaries"]["LumpedPort"][0]["Direction"] = [0.0, 0.05, 1.0]  # 2.9 degrees off
    assert "PRT009" in ids(lint(c, geometry()), "error")


def test_absorbing_too_close():
    g = geometry()
    g["surfaces"]["abc"]["sides"] = ["+z"]
    g["surfaces"]["abc"]["bbox"] = [0, -5, 1.0, 20, 5, 1.0]  # lid 0.49 mm above the strip
    assert "BND001" in ids(lint(good(), g), "warning")


def test_lumped_port_on_absorbing_end_wall():
    g = geometry()
    g["surfaces"]["abc"]["sides"] = ["-x", "+x", "+z"]
    assert "PRT012" in ids(lint(good(), g), "warning")


def test_unknown_attribute_and_missing_material():
    c = good()
    c["Boundaries"]["PEC"]["Attributes"].append(999)
    c["Domains"]["Materials"] = c["Domains"]["Materials"][:1]
    r = lint(c, geometry())
    msgs = [f["message"] for f in r["findings"] if f["id"] == "BND005"]
    assert any("999" in m for m in msgs) and any("no material" in m for m in msgs)


def test_partial_wave_port_with_absorbing():
    c = good()
    c["Boundaries"]["LumpedPort"] = []
    c["Boundaries"]["WavePort"] = [{"Index": 1, "Attributes": [103], "Excitation": 1},
                                   {"Index": 2, "Attributes": [104], "Excitation": 2}]
    assert "PRT013" in ids(lint(c, geometry()), "warning")


def test_real_world_configs_from_projects_are_clean():
    """Configs from the source projects (copied shapes) give no errors."""
    wilk = {"Problem": {"Type": "Driven", "Verbose": 2, "Output": "output"},
            "Model": {"Mesh": "tuned.msh", "L0": 0.001},
            "Domains": {"Materials": [{"Attributes": [1], "Permittivity": 3.55, "LossTan": 0.0027},
                                      {"Attributes": [2], "Permittivity": 1.0}]},
            "Boundaries": {"PEC": {"Attributes": [11]}, "Absorbing": {"Attributes": [12], "Order": 1},
                           "Conductivity": [{"Attributes": [10], "Conductivity": 5.8e7}],
                           "WavePort": [{"Index": i, "Attributes": [20 + i], "Mode": 1, "Offset": 5.0, "Excitation": i}
                                        for i in (1, 2, 3)],
                           "WavePortPEC": {"Attributes": [10]}, "Impedance": [{"Attributes": [30], "Rs": 100.28}]},
            "Solver": {"Order": 2, "Device": "CPU",
                       "Driven": {"Samples": [{"Type": "Linear", "MinFreq": 4.0, "MaxFreq": 12.0, "FreqStep": 0.05,
                                               "SaveStep": 0}], "AdaptiveTol": 0.001},
                       "Linear": {"Type": "MUMPS", "KSPType": "GMRES", "Tol": 1e-08, "MaxIts": 200}}}
    r = lint(copy.deepcopy(wilk))
    assert r["counts"]["error"] == 0
    assert "PRT015" in ids(r)  # offsets with lossy conductors
