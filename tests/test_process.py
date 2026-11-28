import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pytest
from src.helpers.wd_utils.wd_io import LCIO
from src.helpers.wd_utils.wd_containers import LCParameterContainer


def job(tmp_path, code):
    io = LCIO(LCParameterContainer(), str(tmp_path), sys.executable)
    io.command_arguments = ["-c", code]
    io._input = "input"
    return io


def test_failure_never_returns_stale_output(tmp_path):
    (tmp_path / "lcout.active").write_text("old")
    io = job(tmp_path, "import sys; sys.exit(7)")
    with pytest.raises(RuntimeError, match="7"):
        io.save().run()
    assert (tmp_path / "lcout.active").read_text() == "old"


def test_missing_output_is_failure(tmp_path):
    io = job(tmp_path, "pass")
    with pytest.raises(RuntimeError, match="output"):
        io.save().run()


def test_each_job_owns_files_and_arguments(tmp_path):
    code = "from pathlib import Path; Path('lcout.active').write_text(Path('lcin.active').read_text())"
    a, b = job(tmp_path, code), job(tmp_path, code)
    a._input, b._input = "first", "second"
    a.save()
    b.save()
    assert a._get_input_path() != b._get_input_path()
    with ThreadPoolExecutor(2) as pool:
        assert all(f.result() for f in [pool.submit(a.run), pool.submit(b.run)])
    assert Path(a._get_output_path()).read_text() == "first"
    assert Path(b._get_output_path()).read_text() == "second"


def test_cancel_before_spawn(tmp_path):
    io = job(tmp_path, "raise AssertionError('must not run')")
    io.cancel()
    with pytest.raises(RuntimeError, match="ancel"):
        io.save().run()
    assert io.process is None


def test_timeout_reaps_child(tmp_path):
    io = job(tmp_path, "import time; time.sleep(10)")
    io.timeout = 0.1
    with pytest.raises(RuntimeError, match="ime"):
        io.save().run()
    assert io.process is None


def test_cancel_reaps_child(tmp_path):
    io = job(tmp_path, "import time; time.sleep(10)").save()
    timer = threading.Timer(0.2, io.cancel)
    timer.start()
    try:
        with pytest.raises(RuntimeError, match="ancel"):
            io.run()
    finally:
        timer.join()
    assert io.process is None


def test_lc_keeps_qt_event_loop_responsive(app, tmp_path):
    from PyQt5 import QtCore

    io = job(tmp_path, "import time; from pathlib import Path; time.sleep(.2); Path('lcout.active').write_text('ok')")
    ticks = []
    timer = QtCore.QTimer()
    timer.timeout.connect(lambda: ticks.append(1))
    timer.start(10)
    try:
        io.save().run()
    finally:
        timer.stop()
    assert len(ticks) >= 3


def test_worker_catches_failures_and_finishes(app):
    from src.interfaces.dc_interface import IterationWorker

    class Broken:
        def fill_for_solution(self):
            raise PermissionError("denied")

    worker = IterationWorker(Broken())
    errors = []
    worker.failed.connect(errors.append)
    worker.run()
    assert errors and "denied" in errors[0]
    assert worker.runtime is not None


def test_independent_application_processes_do_not_share_files(tmp_path):
    import subprocess

    script = """import sys
from pathlib import Path
from src.helpers.wd_utils.wd_io import LCIO
from src.helpers.wd_utils.wd_containers import LCParameterContainer
io=LCIO(LCParameterContainer(),sys.argv[1],sys.executable)
io.command_arguments=['-c',"from pathlib import Path; Path('lcout.active').write_text(Path('lcin.active').read_text())"]
io._input=sys.argv[2]
io.save().run()
assert Path(io._get_output_path()).read_text()==sys.argv[2]
"""
    processes = [
        subprocess.Popen(
            [sys.executable, "-c", script, str(tmp_path), token], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        for token in ("first", "second")
    ]
    for process in processes:
        stdout, stderr = process.communicate(timeout=10)
        assert process.returncode == 0, stderr.decode(errors="replace")


def test_auxiliary_data_is_copied_and_source_input_link_is_untouched(tmp_path):
    (tmp_path / "atmcof.dat").write_text("auxiliary")
    target = tmp_path / "target"
    target.write_text("original")
    link = tmp_path / "lcin.active"
    try:
        link.symlink_to(target)
    except OSError:
        pass
    io = job(tmp_path, "from pathlib import Path; Path('lcout.active').write_text(Path('atmcof.dat').read_text())")
    io.save().run()
    assert Path(io._get_output_path()).read_text() == "auxiliary"
    assert target.read_text() == "original"
