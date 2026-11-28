from types import SimpleNamespace
import pytest
from src.interfaces import dc_interface, history_interface


@pytest.mark.parametrize("mode", ["plaintext", "latex"])
def test_dc_export_uses_computed_result_and_fractional_kelvin(app, monkeypatch, tmp_path, mode):
    path = tmp_path / "result.txt"

    class Menu:
        def __init__(self, *args):
            pass

        def addAction(self, name):
            return SimpleNamespace(setObjectName=lambda name: None)

        def exec_(self, *args):
            return SimpleNamespace(objectName=lambda: mode)

    monkeypatch.setattr(dc_interface.QtWidgets, "QMenu", Menu)
    monkeypatch.setattr(dc_interface.methods, "save_file", lambda *a, **k: str(path))
    widget = SimpleNamespace(
        results=[[9, 0, 2.0, 0.5, 2.5, 0.01], [19, 0, 0.6, 0.00015, 0.60015, 0.000025]],
        dimensions=[[], []],
        solution_stats=[[0.1]],
        main_window=SimpleNamespace(),
    )
    dc_interface.Widget.export_data(widget)
    text = path.read_text()
    assert "2.5" in text
    assert "6001.5" in text
    assert "0.25" in text


def test_empty_history_does_not_open_export_dialog(app, monkeypatch):
    def unexpected(*a, **k):
        raise AssertionError("empty history opened a dialog")

    monkeypatch.setattr(history_interface.methods, "save_file", unexpected)
    history_interface.Widget.export_data(SimpleNamespace(solutions=[]))
