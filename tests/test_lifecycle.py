import time
from pathlib import Path
import pytest
from PyQt5 import QtCore, QtTest, QtWidgets
from src.interfaces.mainwindow_interface import Widget
from src.helpers.wd_utils import wd_io
from src import constants


def pump_until(predicate):
    end=time.monotonic()+3
    while not predicate() and time.monotonic()<end:
        QtTest.QTest.qWait(10)
    assert predicate()


@pytest.mark.parametrize("cancel",[True,False])
def test_dc_error_or_early_cancel_restores_ui(app,monkeypatch,cancel):
    class IO:
        def __init__(self,*a,**k):self.cancelled=False
        def fill_for_solution(self):return self
        def save(self):return self
        def run(self):
            for _ in range(15):
                if self.cancelled: raise RuntimeError("cancelled")
                time.sleep(.01)
            raise OSError("fixture failure")
        def cancel(self):self.cancelled=True
    monkeypatch.setattr(wd_io,"DCIO",IO)
    w=Widget(app);w.dc_path=".";w.dc_binary="fixture"
    dc=w.dc_widget
    for _ in range(2):
        dc.start_differential_correction()
        worker=dc.iterator
        assert worker is not None
        if cancel: dc.abort_differential_correction()
        assert dc.iterator is worker
        pump_until(lambda:dc.iterator is None)
        assert dc.rundc2015_btn.text()=="Run DC"
        assert dc.rundc2015_btn.isEnabled()
        assert dc.niter_spinbox.isEnabled()
        assert not w.dc_history_widget.solutions
    w.close()


def test_resources_use_module_location():
    assert Path(constants.MAIN_DIR)==Path(__file__).resolve().parents[1]
    assert Path(constants.MAIN_ICON_PATH).is_file()
    assert not Path(constants.PYWD_CONFIG_PATH).is_relative_to(Path(constants.MAIN_DIR))


def test_malformed_lc_table_rejected(tmp_path):
    from src.helpers.wd_utils.wd_containers import LCParameterContainer
    io=wd_io.LCIO(LCParameterContainer(),str(tmp_path),'unused')
    Path(io._get_output_path()).write_text('      JD               Phase     light 1       light 2\n1 2\n3\n\n')
    with pytest.raises(ValueError):io.read_synthetic_light_curve()


def test_progress_dialogs_are_disposed_after_repeated_jobs(app):
    from src.helpers.wd_utils.qt_job import run_with_progress,JobDialog
    class IO:
        timeout=60
        def _run_process(self):pass
    parent=QtWidgets.QMainWindow();parent.show();parent.activateWindow();app.processEvents()
    for _ in range(3):
        run_with_progress(IO())
        QtCore.QCoreApplication.sendPostedEvents(None,QtCore.QEvent.DeferredDelete)
    assert not parent.findChildren(JobDialog)
    parent.close()


def test_invalid_calendar_input_stays_inside_ui_boundary(app):
    window=Widget(app)
    # Direct slot call verifies the same boundary used by clicked.connect.
    window.t_time_year_ipt.setValue(2024)
    window.t_time_month_ipt.setValue(2)
    window.t_time_day_ipt.setValue(31)
    assert window.compute_jd()==1
    window.close()


def test_dc_time_limit_is_configurable(app,tmp_path):
    window=Widget(app);window.dc_path=str(tmp_path);window.dc_binary='unused'
    window.dc_widget.timeout_spinbox.setValue(120)
    window.dc_widget.prepare_iterator()
    assert window.dc_widget.iterator.dc_io.timeout==7200
    window.close()
