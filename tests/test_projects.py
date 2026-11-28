from configparser import ConfigParser
import pytest
from src import constants
from src.helpers import methods
from src.interfaces.mainwindow_interface import Widget


@pytest.fixture
def window(app):
    w = Widget(app)
    yield w
    w.close()


def save(window, path, monkeypatch):
    monkeypatch.setattr(methods, "save_file", lambda *a, **k: str(path))
    window.save_project()


def load(window, path, monkeypatch):
    monkeypatch.setattr(methods, "load_file", lambda *a, **k: str(path))
    window.load_project()


@pytest.mark.parametrize("name", ["100% complete", "%(name)s Türkçe boşluk"])
def test_project_raw_text_roundtrip(window, tmp_path, monkeypatch, name):
    p = tmp_path / "project"
    window.systemname_tipt.setText(name)
    save(window, p, monkeypatch)
    window.systemname_tipt.setText("changed")
    load(window, p, monkeypatch)
    assert window.systemname_tipt.text() == name


@pytest.mark.parametrize("change", ["number", "section", "index", "nan", "curve_count"])
def test_invalid_project_does_not_change_ui(window, tmp_path, monkeypatch, change):
    p = tmp_path / "project"
    window.systemname_tipt.setText("old")
    save(window, p, monkeypatch)
    parser = ConfigParser(interpolation=None)
    parser.read(p)
    parser.set(constants.CONFIG_SECTION_MAIN, "system_name", "bad")
    if change == "number":
        parser.set(constants.CONFIG_SECTION_SYSTEM, "period", "invalid")
    elif change == "section":
        parser.remove_section(constants.CONFIG_SECTION_SURFACE)
    elif change == "index":
        parser.set(constants.CONFIG_SECTION_MAIN, "operation_mode", "999")
    elif change == "nan":
        parser.set(constants.CONFIG_SECTION_SYSTEM, "period", "nan")
    else:
        parser.set(constants.CONFIG_SECTION_CURVE_COUNTS, "light_curves", "-1")
    with p.open("w") as f:
        parser.write(f)
    load(window, p, monkeypatch)
    assert window.systemname_tipt.text() == "old"


def test_project_save_is_atomic(window, tmp_path, monkeypatch):
    import os

    p = tmp_path / "project"
    p.write_text("original")
    monkeypatch.setattr(os, "replace", lambda *args: (_ for _ in ()).throw(OSError("denied")))
    save(window, p, monkeypatch)
    assert p.read_text() == "original"


def test_valid_keeps_survive_constraint_reapplication(window, tmp_path, monkeypatch):
    window.if3b_groupbox.setChecked(True)
    window.ifcgs_chk.setChecked(False)
    window.dc_widget.a3b_chk.setChecked(True)
    window.dc_widget.l1_chk.setChecked(True)
    p = tmp_path / "project"
    save(window, p, monkeypatch)
    window.dc_widget.a3b_chk.setChecked(False)
    window.dc_widget.l1_chk.setChecked(False)
    load(window, p, monkeypatch)
    assert window.dc_widget.a3b_chk.isChecked()
    assert window.dc_widget.l1_chk.isChecked()


def test_invalid_spot_does_not_change_project(window, tmp_path, monkeypatch):
    p = tmp_path / "project"
    window.systemname_tipt.setText("old")
    save(window, p, monkeypatch)
    parser = ConfigParser(interpolation=None)
    parser.read(p)
    parser.set(constants.CONFIG_SECTION_MAIN, "system_name", "bad")
    prefix = constants.CONFIG_SECTION_STAR1_SPOTS
    parser.set(constants.CONFIG_SECTION_SPOT_PARAMS, prefix + "count", "1")
    parser[prefix + "1"] = {
        "a": "False",
        "b": "False",
        "lat": "-1",
        "lon": "0",
        "radsp": "1",
        "temsp": "1",
        "tstart": "0",
        "tmax1": "0",
        "tmax2": "0",
        "tend": "0",
    }
    with p.open("w") as output:
        parser.write(output)
    load(window, p, monkeypatch)
    assert window.systemname_tipt.text() == "old"


def test_legacy_locale_encoded_project_loads(window, tmp_path, monkeypatch):
    from src.helpers import project_io

    monkeypatch.setattr(project_io.locale, "getpreferredencoding", lambda *a: "cp1254")
    p = tmp_path / "legacy"
    parser = window.project_parser()
    parser.set(constants.CONFIG_SECTION_MAIN, "system_name", "Türkçe")
    with p.open("w", encoding="cp1254") as output:
        parser.write(output)
    load(window, p, monkeypatch)
    assert window.systemname_tipt.text() == "Türkçe"


def add_curve(window, tmp_path):
    from src.interfaces.curveproperties_interface import CurveProperty, Widget

    path = tmp_path / "100% observation"
    path.write_text("0 1\n.1 .9\n")
    dialog = Widget("light", window)
    dialog.filepath = str(path)
    dialog.populate_data_from_filepath(str(path))
    dialog.time_combobox.setCurrentIndex(0)
    dialog.obs_combobox.setCurrentIndex(1)
    dialog.weight_combobox.setCurrentIndex(2)
    curve = CurveProperty("light", window)
    curve.populate_from_ui(dialog)
    window.loadobservations_widget.light_curves = [curve]
    return curve


def test_project_with_observations_roundtrips(window, tmp_path, monkeypatch):
    add_curve(window, tmp_path)
    p = tmp_path / "project"
    save(window, p, monkeypatch)
    window.loadobservations_widget.light_curves = []
    load(window, p, monkeypatch)
    assert len(window.loadobservations_widget.light_curves) == 1
    assert window.loadobservations_widget.light_curves[0].constant_weights


@pytest.mark.parametrize("field,value", [("ksd", "3"), ("wla", "1000000000"), ("col2_idx", "99")])
def test_invalid_curve_property_preserves_project(window, tmp_path, monkeypatch, field, value):
    add_curve(window, tmp_path)
    window.systemname_tipt.setText("old")
    p = tmp_path / "project"
    save(window, p, monkeypatch)
    parser = ConfigParser(interpolation=None)
    parser.read(p)
    parser.set(constants.CONFIG_SECTION_MAIN, "system_name", "bad")
    parser.set(constants.CONFIG_SECTION_LIGHT_CURVE_STUB + "1", field, value)
    with p.open("w") as output:
        parser.write(output)
    load(window, p, monkeypatch)
    assert window.systemname_tipt.text() == "old"
