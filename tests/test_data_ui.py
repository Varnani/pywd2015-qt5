from types import SimpleNamespace
from pathlib import Path
import pytest
import numpy as np
from PyQt5 import QtWidgets
from src.interfaces import eclipsetimings_interface,curveproperties_interface,starpositions_interface

@pytest.fixture
def parent(app):
    from src.interfaces.mainwindow_interface import Widget
    window=Widget(app)
    yield window
    window.close()

@pytest.mark.parametrize("text",["1 1\n", "1 1\n2 2\n", "1 1 0.5\n2 2 2\n"])
def test_eclipse_regular_shape_and_weights(parent,tmp_path,text):
    path=tmp_path/"eclipse.txt";path.write_text(text)
    w=parent.eclipsetimings_widget
    assert w.load_timings(str(path)) is not False
    data=w.get_data()
    assert data.shape==(3,len(text.splitlines()))
    if len(text.splitlines()[0].split())==2: assert np.all(data[2]==1)
    w.constant_weight_checkbox.setChecked(True)
    assert np.all(w.get_data()[2]==1)

@pytest.mark.parametrize("text",["", "1\n2\n", "1 nan\n", "1 inf\n", "x 1\n", "1 1 -1\n"])
def test_invalid_eclipse_preserves_old_state(parent,tmp_path,text):
    w=parent.eclipsetimings_widget
    old=tmp_path/"old";old.write_text("1 1\n")
    w.load_timings(str(old))
    path=tmp_path/"bad";path.write_text(text)
    assert w.load_timings(str(path)) is False
    assert w.filepath_label.text()==str(old)


def test_single_observation_curve(parent,tmp_path):
    path=tmp_path/"observation";path.write_text("1 2\n")
    curve=curveproperties_interface.CurveProperty("light",parent)
    curve.filepath=str(path);curve.col1_idx=0;curve.col2_idx=1;curve.col3_idx=2;curve.constant_weights=True
    x,y,w=curve.get_data()
    assert list(x)==[1] and list(y)==[2] and list(w)==[1]


def test_roche_uses_enabled_phase_control(parent,monkeypatch):
    w=parent.lc_starpositions_widget
    w.phase_spinbox.setValue(.25);parent.p0_ipt.setValue(3.7)
    w.roche_groupbox.setChecked(True)
    assert w.phase_spinbox.isEnabled()
    phases=[]
    monkeypatch.setattr(starpositions_interface.methods,"compute_roche_potentials",lambda w,e,q,phase,*a,**k:(phases.append(phase) or (1,2)))
    w.plot_roche()
    assert phases==[.25]

@pytest.mark.parametrize("start,end,step",[(0,1,0),(0,1,-1),(1,0,.1)])
def test_invalid_animation_does_not_run_motor(parent,monkeypatch,start,end,step):
    w=parent.lc_starpositions_widget
    parent.phase_radiobtn.setChecked(True)
    parent.lc_phs_start_ipt.setValue(start);parent.lc_phs_stop_ipt.setValue(end)
    fake=SimpleNamespace(value=lambda:step)
    monkeypatch.setattr(parent,"lc_phs_incr_ipt",fake)
    monkeypatch.setattr(parent,"get_lc_params",lambda:(_ for _ in ()).throw(AssertionError("motor input built")))
    assert w.render_stars()==1
